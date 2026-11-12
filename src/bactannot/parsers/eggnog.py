"""Parser for eggNOG-mapper functional annotations"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

from bactannot.schema import (
    ANNOTATION_COLUMNS,
    EVIDENCE_COLUMNS,
)

# Define la fuente y cómo se traducen los campos de eggNOG al modelo de BactAnnot
SOURCE = "eggNOG-mapper"

# Reconoce identificadores COG que vienen incluidos dentro de los grupos de eggNOG
COG_OG_RE = re.compile(r"^(COG\d+)@[^|]+\|.+$")


# Indica qué tipo de anotación corresponde a cada columna de eggNOG
ANNOTATION_FIELD_MAP = {
    "eggNOG_OGs": "EGGNOG_OG",
    "COG_category": "COG_CATEGORY",
    "GOs": "GO",
    "EC": "EC",
    "KEGG_ko": "KEGG_KO",
    "KEGG_Pathway": "KEGG_PATHWAY",
    "KEGG_Module": "KEGG_MODULE",
    "KEGG_Reaction": "KEGG_REACTION",
    "KEGG_rclass": "KEGG_RCLASS",
    "BRITE": "KEGG_BRITE",
    "KEGG_TC": "KEGG_TC",
    "CAZy": "CAZY",
    "BiGG_Reaction": "BIGG_REACTION",
    "PFAMs": "PFAM_NAME",
}


# Limpia los valores vacíos de eggNOG para que todos se representen igual
def normalize_optional(value: str | None) -> str:
    """Normalize empty eggNOG values"""

    if value is None:
        return ""

    value = value.strip()

    if not value or value == "-":
        return ""

    return value


# Separa los campos de eggNOG que contienen varios valores separados por comas
def split_multivalue(value: str | None) -> list[str]:
    """Split a comma-separated eggNOG annotation field"""

    value = normalize_optional(value)

    if not value:
        return []

    return [item.strip() for item in value.split(",") if item.strip() and item.strip() != "-"]


# Separa cada campo según su formato y trata las categorías COG letra a letra
def split_annotation_values(
    field: str,
    value: str | None,
) -> list[str]:
    """Split an eggNOG annotation field according to its semantics"""

    value = normalize_optional(value)

    if not value:
        return []

    if field == "COG_category":
        return list(value)

    return [item.strip() for item in value.split(",") if item.strip() and item.strip() != "-"]


# Extrae los identificadores COG reales que vienen incluidos dentro de eggNOG_OGs
def extract_cog_ids(
    value: str | None,
) -> list[str]:
    """Extract COG accessions from eggNOG OG assignments"""

    cog_ids: set[str] = set()

    for raw_og in split_multivalue(value):
        match = COG_OG_RE.fullmatch(
            raw_og,
        )

        if match is None:
            continue

        cog_ids.add(match.group(1))

    return sorted(cog_ids)


# Convierte los identificadores de eggNOG al formato común usado por BactAnnot
def canonicalize_annotation_id(
    field: str,
    value: str,
) -> str:
    """Convert source-specific identifiers to BactAnnot canonical IDs"""

    value = value.strip()

    if field == "GOs":
        if value.startswith("GO:"):
            return value

        return f"GO:{value}"

    if field == "EC":
        if value.startswith("EC:"):
            return value

        return f"EC:{value}"

    if field == "KEGG_ko":
        if value.startswith("ko:"):
            value = value.split(":", 1)[1]

        if value.startswith("KEGG:"):
            return value

        return f"KEGG:{value}"

    if field in {
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "KEGG_rclass",
        "BRITE",
    }:
        if value.startswith("KEGG:"):
            return value

        return f"KEGG:{value}"

    return value


# Lee las tablas de eggNOG ignorando sus líneas de metadatos y recuperando las filas
def read_emapper_table(
    path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    """Read an eggNOG-mapper TSV while ignoring metadata comments"""

    header: list[str] | None = None
    rows: list[dict[str, str]] = []

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("##"):
                continue

            if line.startswith("#"):
                header = line[1:].split("\t")

                continue

            if header is None:
                raise ValueError(f"No eggNOG header found in: {path}")

            values = line.split("\t")

            if len(values) != len(header):
                raise ValueError(
                    f"Unexpected column count in {path}: {len(values)} != {len(header)}"
                )

            rows.append(
                dict(
                    zip(
                        header,
                        values,
                        strict=True,
                    )
                )
            )

    if header is None:
        raise ValueError(f"No eggNOG header found in: {path}")

    return header, rows


# Comprueba que una salida de eggNOG tenga las columnas que necesita el parser
def require_columns(
    header: list[str],
    required: set[str],
    path: Path,
) -> None:
    """Validate required columns in an eggNOG output table"""

    missing = sorted(required.difference(header))

    if missing:
        raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")


# Carga la relación protein_id -> feature_uid definida por Bakta para esta muestra
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


# Recupera la feature_uid de una proteína y falla si Bakta no la conoce
def get_feature_uid(
    protein_id: str,
    feature_index: dict[str, str],
) -> str:
    """Resolve an eggNOG query to a canonical feature UID"""

    try:
        return feature_index[protein_id]

    except KeyError as exc:
        raise ValueError(f"eggNOG query is absent from features.tsv: {protein_id}") from exc


# Convierte las anotaciones funcionales de eggNOG al formato canónico de BactAnnot
def parse_eggnog_annotations(
    annotations_path: Path,
    features_path: Path,
    sample_id: str,
) -> list[dict[str, Any]]:
    """Parse eggNOG functional annotations into canonical long format"""

    header, source_rows = read_emapper_table(annotations_path)

    require_columns(
        header,
        {
            "query",
            *ANNOTATION_FIELD_MAP.keys(),
        },
        annotations_path,
    )

    feature_index = load_protein_to_feature_index(
        features_path,
        sample_id,
    )

    output: list[dict[str, Any]] = []

    seen: set[tuple[str, str, str]] = set()

    for source_row in source_rows:
        protein_id = normalize_optional(source_row["query"])

        if not protein_id:
            continue

        feature_uid = get_feature_uid(
            protein_id,
            feature_index,
        )

        for field, annotation_type in ANNOTATION_FIELD_MAP.items():
            for raw_value in split_annotation_values(
                field,
                source_row.get(field),
            ):
                annotation_id = canonicalize_annotation_id(
                    field,
                    raw_value,
                )

                key = (
                    feature_uid,
                    annotation_type,
                    annotation_id,
                )

                if key in seen:
                    continue

                seen.add(key)

                output.append(
                    {
                        "sample_id": sample_id,
                        "feature_uid": feature_uid,
                        "protein_id": protein_id,
                        "source": SOURCE,
                        "annotation_type": annotation_type,
                        "annotation_id": annotation_id,
                        "annotation_value": "",
                        "description": "",
                    }
                )

        # Los COG se extraen también como identificadores propios a partir de eggNOG_OGs
        for cog_id in extract_cog_ids(source_row.get("eggNOG_OGs")):
            annotation_id = f"COG:{cog_id}"

            key = (
                feature_uid,
                "COG",
                annotation_id,
            )

            if key in seen:
                continue

            seen.add(key)

            output.append(
                {
                    "sample_id": sample_id,
                    "feature_uid": feature_uid,
                    "protein_id": protein_id,
                    "source": SOURCE,
                    "annotation_type": "COG",
                    "annotation_id": annotation_id,
                    "annotation_value": "",
                    "description": "",
                }
            )

    return output


# Combina los seed orthologs con la información funcional para construir evidence
def parse_eggnog_evidence(
    annotations_path: Path,
    seed_orthologs_path: Path,
    features_path: Path,
    sample_id: str,
) -> list[dict[str, Any]]:
    """Parse eggNOG seed-ortholog evidence into canonical long format"""

    annotation_header, annotation_rows = read_emapper_table(annotations_path)

    require_columns(
        annotation_header,
        {
            "query",
            "Preferred_name",
            "Description",
        },
        annotations_path,
    )

    seed_header, seed_rows = read_emapper_table(seed_orthologs_path)

    require_columns(
        seed_header,
        {
            "qseqid",
            "sseqid",
            "evalue",
            "bitscore",
            "qstart",
            "qend",
            "pident",
            "qcov",
            "scov",
        },
        seed_orthologs_path,
    )

    feature_index = load_protein_to_feature_index(
        features_path,
        sample_id,
    )

    # Guarda la fila funcional de cada proteína para recuperar su gen y descripción
    annotation_by_query: dict[
        str,
        dict[str, str],
    ] = {}

    for row in annotation_rows:
        query = normalize_optional(row["query"])

        if not query:
            continue

        if query in annotation_by_query:
            raise ValueError(f"Duplicate eggNOG annotation query: {query}")

        annotation_by_query[query] = row

    output: list[dict[str, Any]] = []

    seen_queries: set[str] = set()

    for seed_row in seed_rows:
        protein_id = normalize_optional(seed_row["qseqid"])

        if not protein_id:
            continue

        if protein_id in seen_queries:
            raise ValueError(f"Duplicate eggNOG seed query: {protein_id}")

        seen_queries.add(protein_id)

        feature_uid = get_feature_uid(
            protein_id,
            feature_index,
        )

        annotation_row = annotation_by_query.get(
            protein_id,
            {},
        )

        output.append(
            {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "evidence_type": "seed_ortholog",
                "source": SOURCE,
                "accession": normalize_optional(seed_row.get("sseqid")),
                "gene": normalize_optional(annotation_row.get("Preferred_name")),
                "product": normalize_optional(annotation_row.get("Description")),
                "rank": "",
                "score": normalize_optional(seed_row.get("bitscore")),
                "evalue": normalize_optional(seed_row.get("evalue")),
                "identity": normalize_optional(seed_row.get("pident")),
                "query_coverage": normalize_optional(seed_row.get("qcov")),
                "subject_coverage": normalize_optional(seed_row.get("scov")),
                "aa_coverage": "",
                "hmm_coverage": "",
                "start": normalize_optional(seed_row.get("qstart")),
                "stop": normalize_optional(seed_row.get("qend")),
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
def parse_eggnog_bundle(
    annotations_path: Path,
    seed_orthologs_path: Path,
    features_path: Path,
    sample_id: str,
    annotations_output: Path,
    evidence_output: Path,
) -> None:
    """Parse and write canonical eggNOG annotation and evidence tables"""

    annotations = parse_eggnog_annotations(
        annotations_path=annotations_path,
        features_path=features_path,
        sample_id=sample_id,
    )

    evidence = parse_eggnog_evidence(
        annotations_path=annotations_path,
        seed_orthologs_path=seed_orthologs_path,
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
