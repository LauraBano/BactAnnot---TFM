"""Parser for KofamScan detail-tsv results"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

from bactannot.schema import (
    ANNOTATION_COLUMNS,
    EVIDENCE_COLUMNS,
)

# Define el nombre con el que KofamScan aparecerá como fuente
SOURCE = "KofamScan"


# Define las columnas esperadas en la salida detail-tsv de KofamScan
KOFAMSCAN_COLUMNS = [
    "marker",
    "protein_id",
    "ko",
    "threshold",
    "score",
    "evalue",
    "description",
]


# Reconoce identificadores KEGG Orthology con el formato Kxxxxx
KO_RE = re.compile(r"^K\d{5}$")


# Localiza los bloques [EC:...] incluidos dentro de las descripciones de KOfam
EC_BLOCK_RE = re.compile(r"\[EC:([^\]]+)\]")


# Comprueba que un identificador tenga un formato EC válido
EC_ID_RE = re.compile(
    r"^\d+\."
    r"(?:\d+|-)\."
    r"(?:\d+|-)\."
    r"(?:\d+|-)$"
)


# Limpia los valores ausentes de KofamScan para representarlos siempre de la misma forma
def normalize_optional(
    value: str | None,
) -> str:
    """Normalize missing KofamScan values"""

    if value is None:
        return ""

    value = value.strip()

    if not value or value == "-":
        return ""

    return value


# Extrae los identificadores EC que aparecen explícitamente dentro de la descripción de un KO
def extract_ec_ids(
    description: str | None,
) -> list[str]:
    """Extract explicit EC identifiers from a KOfam KO description"""

    description = normalize_optional(description)

    if not description:
        return []

    ec_ids: set[str] = set()

    # Una misma descripción puede contener uno o varios bloques [EC:...]
    for block in EC_BLOCK_RE.findall(description):
        for token in re.split(
            r"[\s,]+",
            block.strip(),
        ):
            if not token:
                continue

            if EC_ID_RE.fullmatch(token) is None:
                continue

            ec_ids.add(f"EC:{token}")

    return sorted(ec_ids)


# Lee la salida detail-tsv de KofamScan y valida su estructura básica
def read_kofamscan_tsv(
    path: Path,
) -> list[dict[str, str]]:
    """Read a KofamScan detail-tsv result"""

    rows: list[dict[str, str]] = []

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.reader(
            handle,
            delimiter="\t",
            quotechar='"',
        )

        for line_number, values in enumerate(
            reader,
            start=1,
        ):
            if not values:
                continue

            if not any(value.strip() for value in values):
                continue

            if values[0].strip().startswith("#"):
                continue

            if len(values) != len(KOFAMSCAN_COLUMNS):
                raise ValueError(
                    "Unexpected KofamScan column count "
                    f"in {path} at line {line_number}: "
                    f"{len(values)}"
                )

            row = dict(
                zip(
                    KOFAMSCAN_COLUMNS,
                    values,
                    strict=True,
                )
            )

            row = {key: normalize_optional(value) for key, value in row.items()}

            marker = row["marker"]

            # KofamScan marca con * los hits que superan el umbral recomendado
            if marker not in {"", "*"}:
                raise ValueError(
                    f"Unexpected KofamScan marker in {path} at line {line_number}: {marker}"
                )

            protein_id = row["protein_id"]
            ko = row["ko"]

            if not protein_id:
                raise ValueError(
                    f"Missing KofamScan protein identifier in {path} at line {line_number}"
                )

            if not KO_RE.fullmatch(ko):
                raise ValueError(
                    f"Unexpected KofamScan KO identifier in {path} at line {line_number}: {ko}"
                )

            rows.append(row)

    return rows


# Carga la relación protein_id -> feature_uid definida previamente por Bakta
def load_protein_to_feature_index(
    features_path: Path,
    sample_id: str,
) -> dict[str, str]:
    """Map protein IDs to canonical BactAnnot feature UIDs"""

    index: dict[str, str] = {}

    with features_path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        required = {
            "sample_id",
            "feature_uid",
            "protein_id",
        }

        if reader.fieldnames is None:
            raise ValueError(f"Missing header in features table: {features_path}")

        missing = required.difference(reader.fieldnames)

        if missing:
            raise ValueError(
                f"Missing required columns in {features_path}: {', '.join(sorted(missing))}"
            )

        for row in reader:
            row_sample_id = row["sample_id"]
            protein_id = row["protein_id"]
            feature_uid = row["feature_uid"]

            if row_sample_id != sample_id:
                raise ValueError(
                    f"Unexpected sample_id in {features_path}: {row_sample_id} != {sample_id}"
                )

            if not protein_id:
                continue

            if protein_id in index:
                raise ValueError(f"Duplicate protein_id in features table: {protein_id}")

            index[protein_id] = feature_uid

    return index


# Recupera la feature_uid asociada a una proteína y falla si Bakta no la conoce
def get_feature_uid(
    protein_id: str,
    feature_index: dict[str, str],
) -> str:
    """Resolve a KofamScan protein to a feature UID"""

    try:
        return feature_index[protein_id]

    except KeyError as exc:
        raise ValueError(f"KofamScan protein is absent from features.tsv: {protein_id}") from exc


# Convierte los hits aceptados de KofamScan en anotaciones KEGG_KO y EC
def parse_kofamscan_annotations(
    input_path: Path,
    features_path: Path,
    sample_id: str,
) -> list[dict[str, Any]]:
    """Normalize accepted KofamScan KO assignments"""

    source_rows = read_kofamscan_tsv(input_path)

    feature_index = load_protein_to_feature_index(
        features_path,
        sample_id,
    )

    # Guarda cada anotación una sola vez usando feature, tipo e identificador
    annotations: dict[
        tuple[str, str, str],
        dict[str, Any],
    ] = {}

    for source_row in source_rows:
        # Solo usamos como anotaciones los hits que superan el umbral de KofamScan
        if source_row["marker"] != "*":
            continue

        protein_id = source_row["protein_id"]

        feature_uid = get_feature_uid(
            protein_id,
            feature_index,
        )

        ko = source_row["ko"]
        annotation_id = f"KEGG:{ko}"

        key = (
            feature_uid,
            "KEGG_KO",
            annotation_id,
        )

        # Conserva el KO aceptado como anotación funcional principal
        if key not in annotations:
            annotations[key] = {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "source": SOURCE,
                "annotation_type": "KEGG_KO",
                "annotation_id": annotation_id,
                "annotation_value": "",
                "description": source_row["description"],
            }

        # Extrae también los EC que KOfam incluye explícitamente en la descripción del KO
        for ec_id in extract_ec_ids(source_row["description"]):
            ec_key = (
                feature_uid,
                "EC",
                ec_id,
            )

            if ec_key not in annotations:
                annotations[ec_key] = {
                    "sample_id": sample_id,
                    "feature_uid": feature_uid,
                    "protein_id": protein_id,
                    "source": SOURCE,
                    "annotation_type": "EC",
                    "annotation_id": ec_id,
                    "annotation_value": annotation_id,
                    "description": "",
                }

                continue

            # Si varios KO llevan al mismo EC, conserva todos los KO de origen
            origins = {
                value for value in annotations[ec_key]["annotation_value"].split(",") if value
            }

            origins.add(annotation_id)

            annotations[ec_key]["annotation_value"] = ",".join(sorted(origins))

    return list(annotations.values())


# Conserva cada hit aceptado de KofamScan como evidencia del match contra el perfil HMM
def parse_kofamscan_evidence(
    input_path: Path,
    features_path: Path,
    sample_id: str,
) -> list[dict[str, Any]]:
    """Normalize accepted KofamScan HMM evidence"""

    source_rows = read_kofamscan_tsv(input_path)

    feature_index = load_protein_to_feature_index(
        features_path,
        sample_id,
    )

    output: list[dict[str, Any]] = []

    for source_row in source_rows:
        # Los hits sin * se conservan en la salida cruda pero no como evidencia canónica aceptada
        if source_row["marker"] != "*":
            continue

        protein_id = source_row["protein_id"]

        feature_uid = get_feature_uid(
            protein_id,
            feature_index,
        )

        output.append(
            {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "evidence_type": "hmm_match",
                "source": SOURCE,
                "accession": f"KEGG:{source_row['ko']}",
                "gene": "",
                "product": source_row["description"],
                "rank": "above_threshold",
                "score": source_row["score"],
                "native_threshold": source_row["threshold"],
                "evalue": source_row["evalue"],
                "identity": "",
                "query_coverage": "",
                "subject_coverage": "",
                "aa_coverage": "",
                "hmm_coverage": "",
                "start": "",
                "stop": "",
            }
        )

    return output


# Escribe las filas normalizadas usando las columnas de la tabla canónica indicada
def write_tsv(
    rows: list[dict[str, Any]],
    columns: list[str],
    output_path: Path,
) -> None:
    """Write canonical BactAnnot rows to a TSV file"""

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


# Ejecuta el parser completo y genera las tablas de annotations y evidence
def parse_kofamscan_bundle(
    input_path: Path,
    features_path: Path,
    sample_id: str,
    annotations_output: Path,
    evidence_output: Path,
) -> None:
    """Parse and write canonical KofamScan tables"""

    annotations = parse_kofamscan_annotations(
        input_path=input_path,
        features_path=features_path,
        sample_id=sample_id,
    )

    evidence = parse_kofamscan_evidence(
        input_path=input_path,
        features_path=features_path,
        sample_id=sample_id,
    )

    write_tsv(
        rows=annotations,
        columns=ANNOTATION_COLUMNS,
        output_path=annotations_output,
    )

    write_tsv(
        rows=evidence,
        columns=EVIDENCE_COLUMNS,
        output_path=evidence_output,
    )
