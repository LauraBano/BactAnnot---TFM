"""Parser for Bakta JSON annotation results"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from bactannot.schema import (
    ANNOTATION_COLUMNS,
    EVIDENCE_COLUMNS,
)

FEATURE_COLUMNS = (
    "sample_id",
    "feature_uid",
    "contig_id",
    "feature_type",
    "start",
    "stop",
    "strand",
    "bakta_id",
    "locus_tag",
    "protein_id",
    "gene",
    "product",
    "hypothetical",
    "pseudogene",
    "truncated",
    "aa_md5",
    "frame",
    "start_type",
    "rbs_motif",
)

FUNCTIONAL_XREF_TYPES = {
    "EC": "EC",
    "GO": "GO",
    "KEGG": "KEGG_KO",
    "PFAM": "PFAM",
    "RFAM": "RFAM",
    "NCBIFAM": "NCBIFAM",
}

PROTEIN_FEATURE_TYPES = frozenset(
    {
        "cds",
        "sorf",
    }
)


# Carga el JSON de Bakta y comprueba que contiene la lista de features esperada
def load_bakta_json(path: Path) -> dict[str, Any]:
    """Load a Bakta JSON result"""

    if not path.exists():
        raise FileNotFoundError(f"Bakta JSON file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if "features" not in data:
        raise TypeError("Invalid Bakta JSON: missing 'features' field")

    if not isinstance(data["features"], list):
        raise TypeError("Invalid Bakta JSON: 'features' must be a list")

    return data


# Crea un identificador único y reproducible para cada feature dentro de BactAnnot
def build_feature_uid(
    sample_id: str,
    contig_id: str,
    feature_type: str,
    start: int,
    stop: int,
    strand: str,
) -> str:
    """Build a deterministic BactAnnot feature identifier"""

    return f"{sample_id}|{contig_id}|{feature_type}|{start}|{stop}|{strand}"


# Sustituye los valores ausentes por cadenas vacías para poder escribirlos en las tablas
def normalize_optional(value: Any) -> Any:
    """Convert missing values to empty strings for tabular output"""

    if value is None:
        return ""

    return value


# Limpia los valores de texto y unifica las formas de representar datos ausentes
def normalize_text(value: Any) -> str:
    """Normalize optional textual Bakta values"""

    if value is None:
        return ""

    text = str(value).strip()

    if not text or text == "-":
        return ""

    return text


# Unifica la forma de escribir los nombres de las bases de datos
def normalize_xref_namespace(namespace: str) -> str:
    """Normalize Bakta database namespace spelling"""

    aliases = {
        "ec": "EC",
        "go": "GO",
        "kegg": "KEGG",
        "pfam": "PFAM",
        "rfam": "RFAM",
        "cog": "COG",
        "ncbifam": "NCBIFAM",
    }

    normalized = namespace.strip()

    return aliases.get(normalized.lower(), normalized)


# Elimina la versión de los identificadores PFAM y NCBIFAM para usar su forma canónica
def normalize_versioned_xref_identifier(
    namespace: str,
    identifier: str,
) -> str:
    """Remove model versions from supported canonical identifiers"""

    prefixes = {
        "PFAM": "PF",
        "NCBIFAM": "NF",
    }

    prefix = prefixes.get(namespace)

    if prefix is None:
        return identifier

    base, separator, version = identifier.rpartition(".")

    if not separator or not version.isdigit():
        return identifier

    if not base.startswith(prefix):
        return identifier

    accession_number = base[len(prefix) :]

    if not accession_number.isdigit():
        return identifier

    return base


# Normaliza una referencia de base de datos manteniendo una estructura namespace:identificador
def normalize_xref(xref: str) -> str:
    """Normalize a Bakta database cross-reference"""

    text = normalize_text(xref)

    if not text or ":" not in text:
        return text

    namespace, identifier = text.split(
        ":",
        maxsplit=1,
    )

    namespace = normalize_xref_namespace(namespace)

    identifier = identifier.strip()

    identifier = normalize_versioned_xref_identifier(
        namespace,
        identifier,
    )

    return f"{namespace}:{identifier}"


# Extrae de una referencia normalizada el nombre de la base de datos de origen
def get_xref_namespace(xref: str) -> str | None:
    """Return the normalized namespace of a database cross-reference"""

    xref = normalize_xref(xref)

    if ":" not in xref:
        return None

    namespace, _ = xref.split(":", maxsplit=1)

    return namespace or None


# Determina qué herramienta o fuente produce una evidencia de tipo expert
def get_expert_source(record: dict[str, Any]) -> str:
    """Resolve a stable source label for Bakta expert evidence"""

    source = normalize_text(record.get("source"))

    if source:
        return source

    expert_type = normalize_text(record.get("type"))

    known_sources = {
        "amrfinder": "AMRFinderPlus",
        "user_hmms": "UserHMMs",
    }

    if expert_type in known_sources:
        return known_sources[expert_type]

    if expert_type:
        return f"Bakta:{expert_type}"

    return "Bakta:expert"


# Recorre el JSON de forma recursiva y recupera todas las referencias a bases de datos
def collect_bakta_xrefs(value: Any) -> list[str]:
    """Collect every db_xrefs entry recursively from a Bakta structure"""

    xrefs: list[str] = []

    if isinstance(value, dict):
        for key, child in value.items():
            if key == "db_xrefs" and isinstance(child, list):
                xrefs.extend(xref for xref in child if isinstance(xref, str))

            xrefs.extend(collect_bakta_xrefs(child))

    elif isinstance(value, list):
        for item in value:
            xrefs.extend(collect_bakta_xrefs(item))

    return xrefs


# Busca dentro del JSON todos los registros asociados a una clave concreta
def find_nested_records(value: Any, key_name: str) -> list[dict[str, Any]]:
    """Find Bakta dictionaries stored under a key anywhere in a structure"""

    records: list[dict[str, Any]] = []

    if isinstance(value, dict):
        for key, child in value.items():
            if key == key_name:
                if isinstance(child, dict):
                    records.append(child)
                elif isinstance(child, list):
                    records.extend(item for item in child if isinstance(item, dict))

            records.extend(find_nested_records(child, key_name))

    elif isinstance(value, list):
        for item in value:
            records.extend(find_nested_records(item, key_name))

    return records


# Extrae el id del expert record y lo asigna al campo annotation_value
def get_expert_annotation_origins(
    feature: dict[str, Any],
) -> dict[tuple[str, str], set[str]]:
    """Map functional annotations to their originating Bakta expert records"""

    origins: dict[
        tuple[str, str],
        set[str],
    ] = {}

    for record in find_nested_records(
        feature,
        "expert",
    ):
        origin_accession = normalize_text(record.get("id"))

        if not origin_accession:
            continue

        db_xrefs = record.get(
            "db_xrefs",
            [],
        )

        if not isinstance(db_xrefs, list):
            continue

        for raw_xref in db_xrefs:
            if not isinstance(raw_xref, str):
                continue

            for (
                annotation_type,
                annotation_id,
            ) in canonicalize_functional_xref(raw_xref):
                key = (
                    annotation_type,
                    annotation_id,
                )

                origins.setdefault(
                    key,
                    set(),
                ).add(origin_accession)

    return origins


# Convierte una feature individual de Bakta al formato estructural usado por BactAnnot
def parse_feature(
    feature: dict[str, Any],
    sample_id: str,
) -> dict[str, Any]:
    """Convert one Bakta feature to the BactAnnot canonical schema"""

    required = (
        "id",
        "sequence",
        "type",
        "start",
        "stop",
        "strand",
        "product",
    )

    missing = [field for field in required if field not in feature]

    if missing:
        raise ValueError("Bakta feature is missing required fields: " + ", ".join(missing))

    feature_type = str(feature["type"])
    locus_tag = normalize_optional(feature.get("locus"))
    protein_id = locus_tag if feature_type in PROTEIN_FEATURE_TYPES else ""

    return {
        "sample_id": sample_id,
        "feature_uid": build_feature_uid(
            sample_id=sample_id,
            contig_id=str(feature["sequence"]),
            feature_type=feature_type,
            start=int(feature["start"]),
            stop=int(feature["stop"]),
            strand=str(feature["strand"]),
        ),
        "contig_id": feature["sequence"],
        "feature_type": feature_type,
        "start": feature["start"],
        "stop": feature["stop"],
        "strand": feature["strand"],
        "bakta_id": feature["id"],
        "locus_tag": locus_tag,
        "protein_id": protein_id,
        "gene": normalize_optional(feature.get("gene")),
        "product": feature["product"],
        "hypothetical": bool(feature.get("hypothetical", False)),
        "pseudogene": bool(feature.get("pseudogene", False)),
        "truncated": normalize_optional(feature.get("truncated")),
        "aa_md5": normalize_optional(feature.get("aa_hexdigest")),
        "frame": normalize_optional(feature.get("frame")),
        "start_type": normalize_optional(feature.get("start_type")),
        "rbs_motif": normalize_optional(feature.get("rbs_motif")),
    }


# Convierte todas las features de Bakta y comprueba que sus identificadores sean únicos
def parse_bakta_features(
    data: dict[str, Any],
    sample_id: str,
) -> list[dict[str, Any]]:
    """Parse all genomic features from a Bakta JSON result"""

    if not sample_id.strip():
        raise ValueError("sample_id cannot be empty")

    rows = [parse_feature(feature, sample_id) for feature in data["features"]]
    feature_uids = [row["feature_uid"] for row in rows]

    if len(feature_uids) != len(set(feature_uids)):
        raise ValueError("Duplicate feature_uid values detected in Bakta annotations")

    return rows


# Decide qué anotación canónica representa cada referencia funcional encontrada en Bakta
def canonicalize_functional_xref(xref: str) -> list[tuple[str, str]]:
    """Convert a Bakta db_xref into canonical BactAnnot annotations"""

    xref = normalize_xref(xref)

    if ":" not in xref:
        return []

    namespace, identifier = xref.split(":", maxsplit=1)
    identifier = identifier.strip()

    if not identifier:
        return []

    # Bakta usa COG tanto para identificadores COG como para categorías de una letra, por eso aquí distinguimos ambos casos
    if namespace == "COG":
        identifier = identifier.upper()

        if identifier.startswith("COG") and identifier[3:].isdigit():
            return [("COG", f"COG:{identifier}")]

        categories = []
        for category in identifier:
            if category.isalpha() and category not in categories:
                categories.append(category)

        return [("COG_CATEGORY", category) for category in categories]

    annotation_type = FUNCTIONAL_XREF_TYPES.get(namespace)

    if annotation_type is None:
        return []

    # Bakta guarda las anotaciones de KOfam dentro del campo KEGG, por eso aquí conservamos solo identificadores KO válidos
    if namespace == "KEGG":
        identifier = identifier.upper()

        if not (identifier.startswith("K") and identifier[1:].isdigit()):
            return []

        return [("KEGG_KO", f"KEGG:{identifier}")]

    if namespace == "PFAM":
        return [("PFAM", f"PFAM:{identifier.upper()}")]

    if namespace == "RFAM":
        return [("RFAM", f"RFAM:{identifier.upper()}")]

    if namespace == "NCBIFAM":
        return [("NCBIFAM", f"NCBIFAM:{identifier.upper()}")]

    return [(annotation_type, f"{namespace}:{identifier}")]


# Recupera los nombres descriptivos de los dominios Pfam asociados a una feature
def get_pfam_descriptions(feature: dict[str, Any]) -> dict[str, str]:
    """Map canonical Pfam identifiers to their human-readable names"""

    descriptions: dict[str, str] = {}

    for record in find_nested_records(feature, "pfams"):
        pfam_id = normalize_text(record.get("id"))
        name = normalize_text(record.get("name"))

        if not pfam_id or not name:
            continue

        raw_xref = pfam_id if ":" in pfam_id else f"PFAM:{pfam_id}"

        for annotation_type, annotation_id in canonicalize_functional_xref(raw_xref):
            if annotation_type == "PFAM":
                descriptions[annotation_id] = name

    return descriptions


# Recorre las features y construye las anotaciones funcionales canónicas de Bakta
def parse_bakta_annotations(
    data: dict[str, Any],
    sample_id: str,
) -> list[dict[str, Any]]:
    """Extract canonical functional annotations from Bakta features"""

    if not sample_id.strip():
        raise ValueError("sample_id cannot be empty")

    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()

    for feature in data["features"]:
        parsed_feature = parse_feature(feature, sample_id)
        feature_uid = parsed_feature["feature_uid"]
        protein_id = parsed_feature["protein_id"]
        pfam_descriptions = get_pfam_descriptions(feature)
        expert_origins = get_expert_annotation_origins(feature)

        raw_xrefs = collect_bakta_xrefs(feature)
        raw_xrefs.extend(pfam_descriptions)

        for raw_xref in raw_xrefs:
            for annotation_type, annotation_id in canonicalize_functional_xref(raw_xref):
                key = (
                    feature_uid,
                    annotation_type,
                    annotation_id,
                )

                if key in seen:
                    continue

                seen.add(key)

                origins = expert_origins.get(
                    (
                        annotation_type,
                        annotation_id,
                    ),
                    set(),
                )

                annotation_value = ",".join(
                    sorted(
                        origins,
                        key=str.casefold,
                    )
                )

                rows.append(
                    {
                        "sample_id": sample_id,
                        "feature_uid": feature_uid,
                        "protein_id": protein_id,
                        "source": "bakta",
                        "annotation_type": annotation_type,
                        "annotation_id": annotation_id,
                        "annotation_value": annotation_value,
                        "description": pfam_descriptions.get(annotation_id, ""),
                    }
                )

    return rows


# Recupera los identificadores de la feature que deben acompañar a cualquier evidencia
def build_evidence_context(
    feature: dict[str, Any],
    sample_id: str,
) -> dict[str, Any]:
    """Return identifiers shared by all evidence records"""

    parsed_feature = parse_feature(feature, sample_id)

    return {
        "sample_id": sample_id,
        "feature_uid": parsed_feature["feature_uid"],
        "protein_id": parsed_feature["protein_id"],
    }


# Elige el identificador de secuencia más informativo disponible dentro de una evidencia UPS
def get_ups_accession(record: dict[str, Any]) -> tuple[str, str]:
    """Select the most informative available UPS accession and its source"""

    candidates = (
        ("uniref100_id", "UniRef100"),
        ("uniparc_id", "UniParc"),
        ("ncbi_nrp_id", "RefSeq"),
    )

    for field, source in candidates:
        accession = normalize_text(record.get(field))

        if accession:
            return source, accession

    return "", ""


# Convierte las coincidencias UPS de secuencia exacta en filas de evidencia
def parse_ups_evidence(
    feature: dict[str, Any],
    sample_id: str,
) -> list[dict[str, Any]]:
    """Extract UPS exact-sequence lookup evidence from one Bakta feature"""

    context = build_evidence_context(feature, sample_id)
    rows: list[dict[str, Any]] = []

    for ups in find_nested_records(feature, "ups"):
        source, accession = get_ups_accession(ups)

        if not accession:
            continue

        rows.append(
            {
                **context,
                "evidence_type": "ups",
                "source": source,
                "accession": accession,
                "gene": "",
                "product": "",
                "rank": "",
                "score": "",
                "evalue": "",
                "identity": "",
                "query_coverage": "",
                "subject_coverage": "",
                "aa_coverage": "",
                "hmm_coverage": "",
                "start": "",
                "stop": "",
            }
        )

    return rows


# Convierte la información IPS basada en UniRef100 en filas de evidencia
def parse_ips_evidence(
    feature: dict[str, Any],
    sample_id: str,
) -> list[dict[str, Any]]:
    """Extract IPS/UniRef100 evidence from one Bakta feature"""

    context = build_evidence_context(feature, sample_id)
    rows: list[dict[str, Any]] = []

    for ips in find_nested_records(feature, "ips"):
        accession = normalize_text(ips.get("uniref100_id"))

        if not accession:
            continue

        rows.append(
            {
                **context,
                "evidence_type": "ips",
                "source": "UniRef100",
                "accession": accession,
                "gene": normalize_text(ips.get("gene")),
                "product": normalize_text(ips.get("product")),
                "rank": "",
                "score": "",
                "evalue": "",
                "identity": "",
                "query_coverage": "",
                "subject_coverage": "",
                "aa_coverage": "",
                "hmm_coverage": "",
                "start": "",
                "stop": "",
            }
        )

    return rows


# Convierte la información PSC basada en UniRef90 y sus métricas en filas de evidencia
def parse_psc_evidence(
    feature: dict[str, Any],
    sample_id: str,
) -> list[dict[str, Any]]:
    """Extract PSC/UniRef90 evidence from one Bakta feature"""

    context = build_evidence_context(feature, sample_id)
    rows: list[dict[str, Any]] = []

    for psc in find_nested_records(feature, "psc"):
        accession = normalize_text(psc.get("uniref90_id"))

        if not accession:
            continue

        rows.append(
            {
                **context,
                "evidence_type": "psc",
                "source": "UniRef90",
                "accession": accession,
                "gene": normalize_text(psc.get("gene")),
                "product": normalize_text(psc.get("product")),
                "rank": "",
                "score": normalize_optional(psc.get("score")),
                "evalue": normalize_optional(psc.get("evalue")),
                "identity": normalize_optional(psc.get("identity")),
                "query_coverage": normalize_optional(psc.get("query_cov")),
                "subject_coverage": normalize_optional(psc.get("subject_cov")),
                "aa_coverage": "",
                "hmm_coverage": "",
                "start": "",
                "stop": "",
            }
        )

    return rows


# Convierte la información PSCC basada en UniRef50 y sus métricas en filas de evidencia
def parse_pscc_evidence(
    feature: dict[str, Any],
    sample_id: str,
) -> list[dict[str, Any]]:
    """Extract PSCC/UniRef50 evidence from one Bakta feature"""

    context = build_evidence_context(feature, sample_id)
    rows: list[dict[str, Any]] = []

    for pscc in find_nested_records(feature, "pscc"):
        accession = normalize_text(pscc.get("uniref50_id"))

        if not accession:
            continue

        rows.append(
            {
                **context,
                "evidence_type": "pscc",
                "source": "UniRef50",
                "accession": accession,
                "gene": "",
                "product": normalize_text(pscc.get("product")),
                "rank": "",
                "score": normalize_optional(pscc.get("score")),
                "evalue": normalize_optional(pscc.get("evalue")),
                "identity": normalize_optional(pscc.get("identity")),
                "query_coverage": normalize_optional(pscc.get("query_cov")),
                "subject_coverage": normalize_optional(pscc.get("subject_cov")),
                "aa_coverage": "",
                "hmm_coverage": "",
                "start": "",
                "stop": "",
            }
        )

    return rows


# Conserva como evidencia los resultados especializados que Bakta guarda como expert
def parse_expert_evidence(
    feature: dict[str, Any],
    sample_id: str,
) -> list[dict[str, Any]]:
    """Extract specialized expert evidence from one Bakta feature"""

    context = build_evidence_context(feature, sample_id)
    rows: list[dict[str, Any]] = []

    for record in find_nested_records(feature, "expert"):
        rows.append(
            {
                **context,
                "evidence_type": "expert",
                "source": get_expert_source(record),
                "accession": normalize_text(record.get("id")),
                "gene": normalize_text(record.get("gene")),
                "product": normalize_text(record.get("product")),
                "rank": normalize_optional(record.get("rank")),
                "score": normalize_optional(record.get("score")),
                "evalue": normalize_optional(record.get("evalue")),
                "identity": normalize_optional(record.get("identity")),
                "query_coverage": normalize_optional(record.get("query_cov")),
                "subject_coverage": normalize_optional(record.get("subject_cov")),
                "aa_coverage": "",
                "hmm_coverage": "",
                "start": normalize_optional(record.get("start")),
                "stop": normalize_optional(record.get("stop")),
            }
        )

    return rows


# Convierte los dominios Pfam y sus métricas en filas de evidencia
def parse_pfam_evidence(
    feature: dict[str, Any],
    sample_id: str,
) -> list[dict[str, Any]]:
    """Extract Pfam domain evidence from one Bakta feature"""

    context = build_evidence_context(feature, sample_id)
    rows: list[dict[str, Any]] = []

    for record in find_nested_records(feature, "pfams"):
        rows.append(
            {
                **context,
                "evidence_type": "pfam",
                "source": "Pfam",
                "accession": normalize_text(record.get("id")),
                "gene": "",
                "product": normalize_text(record.get("name")),
                "rank": "",
                "score": normalize_optional(record.get("score")),
                "evalue": normalize_optional(record.get("evalue")),
                "identity": "",
                "query_coverage": "",
                "subject_coverage": "",
                "aa_coverage": normalize_optional(record.get("aa_cov")),
                "hmm_coverage": normalize_optional(record.get("hmm_cov")),
                "start": normalize_optional(record.get("start")),
                "stop": normalize_optional(record.get("stop")),
            }
        )

    return rows


# Reúne para cada feature todos los tipos de evidencia de Bakta
def parse_bakta_evidence(
    data: dict[str, Any],
    sample_id: str,
) -> list[dict[str, Any]]:
    """Extract all supported Bakta evidence records"""

    if not sample_id.strip():
        raise ValueError("sample_id cannot be empty")

    rows: list[dict[str, Any]] = []

    for feature in data["features"]:
        rows.extend(parse_ups_evidence(feature, sample_id))
        rows.extend(parse_ips_evidence(feature, sample_id))
        rows.extend(parse_psc_evidence(feature, sample_id))
        rows.extend(parse_pscc_evidence(feature, sample_id))
        rows.extend(parse_expert_evidence(feature, sample_id))
        rows.extend(parse_pfam_evidence(feature, sample_id))

    return rows


# Escritura de las tablas


# Escribe cualquier conjunto de filas usando las columnas del esquema que se le indique
def write_tsv(
    rows: Iterable[dict[str, Any]],
    output_path: Path,
    columns: tuple[str, ...],
) -> None:
    """Write dictionaries to a tab-separated file"""

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(row)


# Escribe la tabla estructural de features
def write_features_tsv(
    rows: Iterable[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write canonical Bakta features to a TSV file"""

    write_tsv(rows, output_path, FEATURE_COLUMNS)


# Escribe la tabla de anotaciones funcionales
def write_annotations_tsv(
    rows: Iterable[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write normalized functional annotations to a TSV file"""

    write_tsv(rows, output_path, ANNOTATION_COLUMNS)


# Escribe la tabla de evidencias
def write_evidence_tsv(
    rows: Iterable[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write Bakta evidence records to a TSV file"""

    write_tsv(rows, output_path, EVIDENCE_COLUMNS)


# Ejecuta únicamente la parte necesaria para generar features.tsv a partir del JSON
def parse_bakta_json(
    input_path: Path,
    sample_id: str,
    output_path: Path,
) -> None:
    """Parse a Bakta JSON file and write features.tsv"""

    data = load_bakta_json(input_path)
    rows = parse_bakta_features(data, sample_id)
    write_features_tsv(rows, output_path)


# Ejecuta el parser completo y genera features, annotations y evidence en una sola llamada
def parse_bakta_bundle(
    input_path: Path,
    sample_id: str,
    features_output: Path,
    annotations_output: Path,
    evidence_output: Path,
) -> None:
    """Parse all supported Bakta information"""

    data = load_bakta_json(input_path)

    features = parse_bakta_features(data, sample_id)
    annotations = parse_bakta_annotations(data, sample_id)
    evidence = parse_bakta_evidence(data, sample_id)

    write_features_tsv(features, features_output)
    write_annotations_tsv(annotations, annotations_output)
    write_evidence_tsv(evidence, evidence_output)
