"""Parser for InterProScan protein annotation results"""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

from bactannot.schema import (
    ANNOTATION_COLUMNS,
    EVIDENCE_COLUMNS,
)

# Define el nombre con el que InterProScan aparecerá como fuente en annotations
SOURCE = "InterProScan"


# Define las columnas del TSV de InterProScan en el orden en que aparecen
INTERPROSCAN_COLUMNS = [
    "protein_id",
    "sequence_md5",
    "sequence_length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "stop",
    "score",
    "status",
    "date",
    "interpro_accession",
    "interpro_description",
    "go_terms",
    "pathways",
]


# Traduce cada base de datos analizada por InterProScan al tipo de anotación usado por BactAnnot
ANALYSIS_MAP = {
    "AntiFam": ("ANTIFAM", "ANTIFAM"),
    "CDD": ("CDD", "CDD"),
    "Coils": ("COILS", "COILS"),
    "FunFam": ("FUNFAM", "FUNFAM"),
    "Gene3D": ("GENE3D", "GENE3D"),
    "Hamap": ("HAMAP", "HAMAP"),
    "MobiDBLite": ("MOBIDBLITE", "MOBIDBLITE"),
    "NCBIfam": ("NCBIFAM", "NCBIFAM"),
    "PANTHER": ("PANTHER", "PANTHER"),
    "Pfam": ("PFAM", "PFAM"),
    "PIRSF": ("PIRSF", "PIRSF"),
    "PIRSR": ("PIRSR", "PIRSR"),
    "PRINTS": ("PRINTS", "PRINTS"),
    "ProSitePatterns": ("PROSITE", "PROSITE"),
    "ProSiteProfiles": ("PROSITE", "PROSITE"),
    "SFLD": ("SFLD", "SFLD"),
    "SMART": ("SMART", "SMART"),
    "SUPERFAMILY": ("SUPERFAMILY", "SUPERFAMILY"),
}


# Reconoce términos GO y, cuando aparece, la fuente indicada entre paréntesis
GO_TERM_RE = re.compile(r"^(GO:\d+)(?:\(([^)]*)\))?$")


# Limpia los valores ausentes de InterProScan para representarlos siempre de la misma forma
def normalize_optional(value: str | None) -> str:
    """Normalize missing InterProScan values"""

    if value is None:
        return ""

    value = value.strip()

    if not value or value == "-":
        return ""

    return value


# Lee el TSV de InterProScan y convierte cada línea en un diccionario con nombres de columna
def read_interproscan_tsv(
    path: Path,
) -> list[dict[str, str]]:
    """Read an InterProScan TSV result"""

    rows: list[dict[str, str]] = []

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line_number, line in enumerate(
            handle,
            start=1,
        ):
            line = line.rstrip("\n")

            if not line:
                continue

            values = line.split("\t")

            if not 11 <= len(values) <= 15:
                raise ValueError(
                    "Unexpected InterProScan column count "
                    f"in {path} at line {line_number}: "
                    f"{len(values)}"
                )

            # Completa con valores vacíos las columnas opcionales que InterProScan no haya escrito
            values.extend(["-"] * (len(INTERPROSCAN_COLUMNS) - len(values)))

            rows.append(
                dict(
                    zip(
                        INTERPROSCAN_COLUMNS,
                        values,
                        strict=True,
                    )
                )
            )

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
    """Resolve an InterProScan protein to a feature UID"""

    try:
        return feature_index[protein_id]

    except KeyError as exc:
        raise ValueError(f"InterProScan protein is absent from features.tsv: {protein_id}") from exc


# Convierte una firma de una base miembro de InterProScan al identificador canónico de BactAnnot
def canonical_signature(
    analysis: str,
    accession: str,
) -> tuple[str, str]:
    """Normalize a member-database signature accession"""

    try:
        annotation_type, prefix = ANALYSIS_MAP[analysis]

    except KeyError as exc:
        raise ValueError(f"Unsupported InterProScan analysis: {analysis}") from exc

    prefix_with_colon = f"{prefix}:"

    if accession.startswith(prefix_with_colon):
        annotation_id = accession
    else:
        annotation_id = f"{prefix}:{accession}"

    return annotation_type, annotation_id


# Separa los términos GO de una fila y conserva también las fuentes que InterProScan les atribuye
def parse_go_terms(
    value: str | None,
) -> list[tuple[str, set[str]]]:
    """Parse GO identifiers and their reported sources"""

    value = normalize_optional(value)

    if not value:
        return []

    parsed: list[tuple[str, set[str]]] = []

    for item in value.split("|"):
        item = item.strip()

        if not item:
            continue

        match = GO_TERM_RE.fullmatch(item)

        if match is None:
            raise ValueError(f"Unexpected InterProScan GO value: {item}")

        go_id = match.group(1)
        source_text = match.group(2)

        sources: set[str] = set()

        if source_text:
            sources.update(source.strip() for source in source_text.split(",") if source.strip())

        parsed.append(
            (
                go_id,
                sources,
            )
        )

    return parsed


# Guarda las relaciones en annotation_value mediante JSON
def serialize_go_origins(
    origins: set[
        tuple[
            str,
            str,
            str,
            str,
        ]
    ],
) -> str:
    """Serialize exact InterProScan GO provenance deterministically"""

    ordered_origins = sorted(
        origins,
        key=lambda origin: (
            tuple(value.casefold() for value in origin),
            origin,
        ),
    )

    payload = [
        {
            "analysis": analysis,
            "signature_accession": signature_accession,
            "interpro_accession": interpro_accession,
            "mapping_source": mapping_source,
        }
        for (
            analysis,
            signature_accession,
            interpro_accession,
            mapping_source,
        ) in ordered_origins
    ]

    return json.dumps(
        payload,
        ensure_ascii=False,
        separators=(",", ":"),
    )


# Convierte las firmas, entradas InterPro y términos GO al formato canónico de annotations
def parse_interproscan_annotations(
    input_path: Path,
    features_path: Path,
    sample_id: str,
) -> list[dict[str, Any]]:
    """Normalize InterProScan annotations"""

    source_rows = read_interproscan_tsv(input_path)

    feature_index = load_protein_to_feature_index(
        features_path,
        sample_id,
    )

    # Guarda las anotaciones por su identidad para evitar duplicados
    annotations: dict[
        tuple[str, str, str],
        dict[str, Any],
    ] = {}

    # Guarda por separado las fuentes asociadas a cada término GO
    go_origins: dict[
        tuple[str, str, str],
        set[
            tuple[
                str,
                str,
                str,
                str,
            ]
        ],
    ] = defaultdict(set)

    # Añade una anotación solo una vez y completa su descripción si aparece más tarde
    def add_annotation(
        *,
        feature_uid: str,
        protein_id: str,
        annotation_type: str,
        annotation_id: str,
        description: str = "",
    ) -> tuple[str, str, str]:

        key = (
            feature_uid,
            annotation_type,
            annotation_id,
        )

        if key not in annotations:
            annotations[key] = {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "source": SOURCE,
                "annotation_type": annotation_type,
                "annotation_id": annotation_id,
                "annotation_value": "",
                "description": description,
            }

        elif not annotations[key]["description"] and description:
            annotations[key]["description"] = description

        return key

    for source_row in source_rows:
        protein_id = normalize_optional(source_row["protein_id"])

        if not protein_id:
            continue

        feature_uid = get_feature_uid(
            protein_id,
            feature_index,
        )

        analysis = normalize_optional(source_row["analysis"])

        accession = normalize_optional(source_row["signature_accession"])

        description = normalize_optional(source_row["signature_description"])

        if not analysis or not accession:
            raise ValueError(f"InterProScan match is missing analysis/accession for {protein_id}")

        # Primero conserva la firma específica encontrada por la base miembro
        (
            annotation_type,
            annotation_id,
        ) = canonical_signature(
            analysis,
            accession,
        )

        add_annotation(
            feature_uid=feature_uid,
            protein_id=protein_id,
            annotation_type=annotation_type,
            annotation_id=annotation_id,
            description=description,
        )

        # Si esa firma está integrada en una entrada InterPro, conserva también esa anotación
        interpro_accession = normalize_optional(source_row["interpro_accession"])

        interpro_description = normalize_optional(source_row["interpro_description"])

        if interpro_accession:
            if interpro_accession.startswith("INTERPRO:"):
                interpro_id = interpro_accession
            else:
                interpro_id = f"INTERPRO:{interpro_accession}"

            add_annotation(
                feature_uid=feature_uid,
                protein_id=protein_id,
                annotation_type="INTERPRO",
                annotation_id=interpro_id,
                description=interpro_description,
            )

        # Añade los términos GO asociados y acumula sus fuentes si aparecen varias veces
        for go_id, mapping_sources in parse_go_terms(source_row["go_terms"]):
            key = add_annotation(
                feature_uid=feature_uid,
                protein_id=protein_id,
                annotation_type="GO",
                annotation_id=go_id,
            )

            if mapping_sources:
                for mapping_source in mapping_sources:
                    go_origins[key].add(
                        (
                            analysis,
                            accession,
                            interpro_accession,
                            mapping_source,
                        )
                    )

            else:
                go_origins[key].add(
                    (
                        analysis,
                        accession,
                        interpro_accession,
                        "",
                    )
                )

    # Guarda en annotation_value todas las fuentes que respaldan cada término GO
    for key, origins in go_origins.items():
        annotations[key]["annotation_value"] = serialize_go_origins(origins)

    return list(annotations.values())


# Conserva cada match de firma de InterProScan como una fila independiente de evidence
def parse_interproscan_evidence(
    input_path: Path,
    features_path: Path,
    sample_id: str,
) -> list[dict[str, Any]]:
    """Normalize InterProScan signature-match evidence"""

    source_rows = read_interproscan_tsv(input_path)

    feature_index = load_protein_to_feature_index(
        features_path,
        sample_id,
    )

    output: list[dict[str, Any]] = []

    for source_row in source_rows:
        protein_id = normalize_optional(source_row["protein_id"])

        if not protein_id:
            continue

        feature_uid = get_feature_uid(
            protein_id,
            feature_index,
        )

        analysis = normalize_optional(source_row["analysis"])

        accession = normalize_optional(source_row["signature_accession"])

        if not analysis or not accession:
            raise ValueError(f"InterProScan match is missing analysis/accession for {protein_id}")

        output.append(
            {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "evidence_type": "signature_match",
                "source": (f"InterProScan:{analysis}"),
                "accession": accession,
                "gene": "",
                "product": normalize_optional(source_row["signature_description"]),
                "rank": "",
                "score": normalize_optional(source_row["score"]),
                "evalue": "",
                "identity": "",
                "query_coverage": "",
                "subject_coverage": "",
                "aa_coverage": "",
                "hmm_coverage": "",
                "start": normalize_optional(source_row["start"]),
                "stop": normalize_optional(source_row["stop"]),
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
def parse_interproscan_bundle(
    input_path: Path,
    features_path: Path,
    sample_id: str,
    annotations_output: Path,
    evidence_output: Path,
) -> None:
    """Parse and write canonical InterProScan tables"""

    annotations = parse_interproscan_annotations(
        input_path=input_path,
        features_path=features_path,
        sample_id=sample_id,
    )

    evidence = parse_interproscan_evidence(
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
