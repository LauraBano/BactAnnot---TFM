"""Parser for PathoFact functional annotations"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

from bactannot.schema import (
    ANNOTATION_COLUMNS,
    EVIDENCE_COLUMNS,
)

# Define el nombre con el que PathoFact aparecerá como fuente en las tablas canónicas
SOURCE = "PathoFact2"


# Reconoce los distintos tipos de perfiles HMM que PathoFact puede devolver
PFAM_RE = re.compile(
    r"^pfam(\d+)$",
    re.IGNORECASE,
)

CDD_RE = re.compile(
    r"^cd(\d+)$",
    re.IGNORECASE,
)

SMART_RE = re.compile(
    r"^smart(\d+)$",
    re.IGNORECASE,
)

TIGR_RE = re.compile(
    r"^TIGR\d+$",
    re.IGNORECASE,
)

NF_RE = re.compile(
    r"^NF\d+$",
    re.IGNORECASE,
)

COG_RE = re.compile(
    r"^COG\d+$",
    re.IGNORECASE,
)

KOG_RE = re.compile(
    r"^KOG\d+$",
    re.IGNORECASE,
)


# Limpia valores ausentes o especiales para representarlos siempre como cadenas vacías
def normalize_optional(
    value: Any,
) -> str:
    if value is None:
        return ""

    text = str(value).strip()

    if not text or text == "-":
        return ""

    return text


# Convierte los identificadores de perfiles HMM de PathoFact al formato canónico de BactAnnot
def canonical_signature(
    accession: str,
) -> tuple[str, str]:
    """
    Convert a PathoFact2 HMM profile identifier into
    a BactAnnot canonical annotation identifier.

    Only equivalences that can be established from the
    source accession itself are normalized.
    """

    accession = accession.strip()

    match = PFAM_RE.fullmatch(accession)

    if match:
        pfam_id = f"PF{match.group(1)}"

        return (
            "PFAM",
            f"PFAM:{pfam_id}",
        )

    match = CDD_RE.fullmatch(accession)

    if match:
        cdd_id = f"cd{match.group(1)}"

        return (
            "CDD",
            f"CDD:{cdd_id}",
        )

    match = SMART_RE.fullmatch(accession)

    if match:
        smart_id = f"SM{match.group(1)}"

        return (
            "SMART",
            f"SMART:{smart_id}",
        )

    accession_upper = accession.upper()

    if TIGR_RE.fullmatch(accession) or NF_RE.fullmatch(accession):
        return (
            "NCBIFAM",
            (f"NCBIFAM:{accession_upper}"),
        )

    if COG_RE.fullmatch(accession):
        return (
            "COG",
            f"COG:{accession_upper}",
        )

    if KOG_RE.fullmatch(accession):
        return (
            "KOG",
            f"KOG:{accession_upper}",
        )

    return (
        "CDD_PROFILE",
        (f"CDD_PROFILE:{accession_upper}"),
    )


# Carga la relación protein_id -> feature_uid definida previamente por Bakta
def load_protein_to_feature_index(
    features_path: Path,
    sample_id: str,
) -> dict[str, str]:
    """
    Map canonical BactAnnot protein_id values to
    feature_uid values.
    """

    index: dict[str, str] = {}

    with features_path.open(
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

        if not reader.fieldnames:
            raise ValueError("Feature table has no header.")

        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(f"Feature table is missing required columns: {sorted(missing)}")

        for row in reader:
            if row["sample_id"] != sample_id:
                continue

            protein_id = row["protein_id"].strip()

            if not protein_id:
                continue

            if protein_id in index:
                raise ValueError(f"Duplicate protein_id in feature table: {protein_id}")

            index[protein_id] = row["feature_uid"]

    return index


# Reconstruye los clusters de MMseqs2 para saber qué proteínas originales representa cada secuencia
def load_cluster_map(
    cluster_path: Path,
) -> dict[str, list[str]]:
    """
    Map each representative PathoFact2 protein to all
    original proteins belonging to its MMseqs2 cluster.
    """

    clusters: dict[
        str,
        set[str],
    ] = defaultdict(set)

    member_to_rep: dict[
        str,
        str,
    ] = {}

    with cluster_path.open(
        encoding="utf-8",
    ) as handle:
        for line_number, line in enumerate(
            handle,
            start=1,
        ):
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) != 2:
                raise ValueError(
                    "Unexpected PathoFact2 "
                    "cluster row at line "
                    f"{line_number}: "
                    f"{len(fields)} columns"
                )

            representative, member = fields

            existing = member_to_rep.get(member)

            if existing is not None and existing != representative:
                raise ValueError(f"Protein belongs to more than one cluster: {member}")

            member_to_rep[member] = representative

            clusters[representative].add(member)

            clusters[representative].add(representative)

    return {representative: sorted(members) for representative, members in clusters.items()}


# Expande una proteína representante a todas las proteínas idénticas de su cluster
def expand_cluster(
    protein_id: str,
    clusters: dict[
        str,
        list[str],
    ],
) -> list[str]:
    """
    Expand a PathoFact2 representative protein to all
    identical proteins in its MMseqs2 cluster.
    """

    return clusters.get(
        protein_id,
        [protein_id],
    )


# Recupera la feature_uid de una proteína y falla si no existe en la tabla canónica de Bakta
def validate_protein(
    protein_id: str,
    feature_index: dict[
        str,
        str,
    ],
) -> str:
    try:
        return feature_index[protein_id]

    except KeyError as exc:
        raise ValueError(
            f"PathoFact2 protein_id is not present in the canonical feature table: {protein_id}"
        ) from exc


# Lee los hits HMM de virulencia o toxinas y conserva sus métricas principales
def read_hmm_hits(
    path: Path,
) -> list[dict[str, str]]:
    """
    Read the PyHMMER 'targets' table generated by
    PathoFact2.

    PathoFact2 uses:
        target protein -> column 1
        query profile  -> column 3
        full E-value   -> column 5
        bit score      -> column 6
    """

    rows: list[dict[str, str]] = []

    with path.open(
        encoding="utf-8",
    ) as handle:
        for line_number, line in enumerate(
            handle,
            start=1,
        ):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split(maxsplit=18)

            if len(fields) < 18:
                raise ValueError(
                    f"Unexpected PathoFact2 HMM row at line {line_number}: {len(fields)} fields"
                )

            rows.append(
                {
                    "protein_id": (fields[0]),
                    "accession": (fields[2]),
                    "evalue": (fields[4]),
                    "score": (fields[5]),
                    "description": normalize_optional(fields[18] if len(fields) > 18 else ""),
                }
            )

    return rows


# Lee las predicciones de machine learning de virulencia o toxinas
def read_ml_predictions(
    path: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    with path.open(
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        expected = {
            "Sequence",
            "Prediction",
            "Probability",
        }

        if not reader.fieldnames:
            raise ValueError("PathoFact2 ML table has no header.")

        if not expected.issubset(reader.fieldnames):
            raise ValueError("Unexpected PathoFact2 ML header.")

        for row in reader:
            rows.append(
                {
                    "protein_id": (row["Sequence"]),
                    "prediction": (row["Prediction"]),
                    "probability": (row["Probability"]),
                }
            )

    return rows


# Lee los resultados de DeepARG ya normalizados por PathoFact y recupera sus anotaciones AMR
def read_deeparg_amr(
    path: Path,
) -> list[dict[str, str]]:
    """
    Read a PathoFact2 DeepARG table normalized by argNorm.

    The table preserves the original DeepARG prediction
    fields and adds CARD ARO identifiers and resistance
    ontology mappings.
    """

    rows: list[dict[str, str]] = []

    with path.open(
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        expected = {
            "#ARG",
            "query-start",
            "query-end",
            "read_id",
            "predicted_ARG-class",
            "best-hit",
            "probability",
            "identity",
            "alignment-length",
            "alignment-bitscore",
            "alignment-evalue",
            "counts",
            "ARO",
            "confers_resistance_to",
            "resistance_to_drug_classes",
        }

        if not reader.fieldnames:
            raise ValueError("PathoFact2 DeepARG table has no header.")

        if not expected.issubset(reader.fieldnames):
            raise ValueError("Unexpected PathoFact2 DeepARG header.")

        for row in reader:
            rows.append(
                {
                    "protein_id": normalize_optional(row["read_id"]),
                    "arg": normalize_optional(row["#ARG"]),
                    "arg_class": normalize_optional(row["predicted_ARG-class"]),
                    "best_hit": normalize_optional(row["best-hit"]),
                    "probability": normalize_optional(row["probability"]),
                    "identity": normalize_optional(row["identity"]),
                    "alignment_length": normalize_optional(row["alignment-length"]),
                    "bitscore": normalize_optional(row["alignment-bitscore"]),
                    "evalue": normalize_optional(row["alignment-evalue"]),
                    "start": normalize_optional(row["query-start"]),
                    "stop": normalize_optional(row["query-end"]),
                    "aro": normalize_optional(row["ARO"]),
                    "confers_resistance_to": normalize_optional(row["confers_resistance_to"]),
                    "resistance_to_drug_classes": normalize_optional(
                        row["resistance_to_drug_classes"]
                    ),
                }
            )

    return rows


# Lee los resultados de AMRFinderPlus ya normalizados por PathoFact y recupera sus anotaciones AMR
def read_amrfinder_amr(
    path: Path,
) -> list[dict[str, str]]:
    """
    Read a PathoFact2 AMRFinderPlus table normalized by argNorm.

    The table preserves the original AMRFinderPlus fields
    and adds CARD ARO identifiers and resistance ontology
    mappings.
    """

    rows: list[dict[str, str]] = []

    with path.open(
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        expected = {
            "Protein identifier",
            "Gene symbol",
            "Sequence name",
            "Scope",
            "Element type",
            "Element subtype",
            "Class",
            "Subclass",
            "Method",
            "Target length",
            "Reference sequence length",
            "% Coverage of reference sequence",
            "% Identity to reference sequence",
            "Alignment length",
            "Accession of closest sequence",
            "Name of closest sequence",
            "HMM id",
            "HMM description",
            "ARO",
            "confers_resistance_to",
            "resistance_to_drug_classes",
        }

        if not reader.fieldnames:
            raise ValueError("PathoFact2 AMRFinderPlus table has no header.")

        if not expected.issubset(reader.fieldnames):
            raise ValueError("Unexpected PathoFact2 AMRFinderPlus header.")

        for row in reader:
            rows.append(
                {
                    "protein_id": normalize_optional(row["Protein identifier"]),
                    "gene_symbol": normalize_optional(row["Gene symbol"]),
                    "sequence_name": normalize_optional(row["Sequence name"]),
                    "scope": normalize_optional(row["Scope"]),
                    "element_type": normalize_optional(row["Element type"]),
                    "element_subtype": normalize_optional(row["Element subtype"]),
                    "amr_class": normalize_optional(row["Class"]),
                    "amr_subclass": normalize_optional(row["Subclass"]),
                    "method": normalize_optional(row["Method"]),
                    "target_length": normalize_optional(row["Target length"]),
                    "reference_length": normalize_optional(row["Reference sequence length"]),
                    "reference_coverage": normalize_optional(
                        row["% Coverage of reference sequence"]
                    ),
                    "identity": normalize_optional(row["% Identity to reference sequence"]),
                    "alignment_length": normalize_optional(row["Alignment length"]),
                    "closest_accession": normalize_optional(row["Accession of closest sequence"]),
                    "closest_sequence": normalize_optional(row["Name of closest sequence"]),
                    "hmm_id": normalize_optional(row["HMM id"]),
                    "hmm_description": normalize_optional(row["HMM description"]),
                    "aro": normalize_optional(row["ARO"]),
                    "confers_resistance_to": normalize_optional(row["confers_resistance_to"]),
                    "resistance_to_drug_classes": (
                        normalize_optional(row["resistance_to_drug_classes"])
                    ),
                }
            )

    return rows


# Normaliza los identificadores ARO de RGI al formato ARO:xxxx usado por BactAnnot
def normalize_rgi_aro(
    value: Any,
) -> str:
    """Convert an RGI CARD ARO number to its canonical identifier."""

    aro = normalize_optional(value)

    if not aro:
        return ""

    if aro.upper().startswith("ARO:"):
        return f"ARO:{aro.split(':', maxsplit=1)[1]}"

    return f"ARO:{aro}"


# Lee los resultados de RGI/CARD y adapta sus campos al formato que usará el parser
def read_rgi_amr(
    path: Path,
) -> list[dict[str, str]]:
    """
    Read a PathoFact2 RGI/CARD result table.

    RGI may preserve the complete FASTA header in ORF_ID.
    The canonical BactAnnot protein_id is therefore the
    first whitespace-delimited token of that field.
    """

    rows: list[dict[str, str]] = []

    with path.open(
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        expected = {
            "ORF_ID",
            "Contig",
            "Start",
            "Stop",
            "Orientation",
            "Cut_Off",
            "Pass_Bitscore",
            "Best_Hit_Bitscore",
            "Best_Hit_ARO",
            "Best_Identities",
            "ARO",
            "Model_type",
            "SNPs_in_Best_Hit_ARO",
            "Other_SNPs",
            "Drug Class",
            "Resistance Mechanism",
            "AMR Gene Family",
            "Predicted_DNA",
            "Predicted_Protein",
            "CARD_Protein_Sequence",
            "Percentage Length of Reference Sequence",
            "ID",
            "Model_ID",
            "Nudged",
            "Note",
            "Hit_Start",
            "Hit_End",
            "Antibiotic",
        }

        if not reader.fieldnames:
            raise ValueError("PathoFact2 RGI table has no header.")

        if not expected.issubset(reader.fieldnames):
            raise ValueError("Unexpected PathoFact2 RGI header.")

        for row in reader:
            orf_id = normalize_optional(row["ORF_ID"])

            protein_id = orf_id.split(maxsplit=1)[0] if orf_id else ""

            rows.append(
                {
                    "protein_id": protein_id,
                    "orf_id": orf_id,
                    "contig": normalize_optional(row["Contig"]),
                    "start": normalize_optional(row["Start"]),
                    "stop": normalize_optional(row["Stop"]),
                    "orientation": normalize_optional(row["Orientation"]),
                    "cutoff": normalize_optional(row["Cut_Off"]),
                    "pass_bitscore": normalize_optional(row["Pass_Bitscore"]),
                    "best_hit_bitscore": normalize_optional(row["Best_Hit_Bitscore"]),
                    "best_hit": normalize_optional(row["Best_Hit_ARO"]),
                    "identity": normalize_optional(row["Best_Identities"]),
                    "aro": normalize_rgi_aro(row["ARO"]),
                    "model_type": normalize_optional(row["Model_type"]),
                    "snps_in_best_hit": normalize_optional(row["SNPs_in_Best_Hit_ARO"]),
                    "other_snps": normalize_optional(row["Other_SNPs"]),
                    "drug_class": normalize_optional(row["Drug Class"]),
                    "resistance_mechanism": (normalize_optional(row["Resistance Mechanism"])),
                    "amr_gene_family": normalize_optional(row["AMR Gene Family"]),
                    "reference_coverage": (
                        normalize_optional(row["Percentage Length of Reference Sequence"])
                    ),
                    "model_id": normalize_optional(row["Model_ID"]),
                    "nudged": normalize_optional(row["Nudged"]),
                    "note": normalize_optional(row["Note"]),
                    "hit_start": normalize_optional(row["Hit_Start"]),
                    "hit_end": normalize_optional(row["Hit_End"]),
                    "antibiotic": normalize_optional(row["Antibiotic"]),
                }
            )

    return rows


# Lee los listados finales de PathoFact y devuelve las proteínas clasificadas en cada categoría
def read_classification_report(
    path: Path,
) -> set[str]:
    """
    Read a PathoFact2 high-confidence or ambiguous
    summary table and return representative protein IDs.
    """

    proteins: set[str] = set()

    with path.open(
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if not reader.fieldnames or "ProteinID" not in reader.fieldnames:
            raise ValueError("PathoFact2 summary table does not contain ProteinID.")

        for row in reader:
            protein_id = row["ProteinID"].strip()

            if protein_id:
                proteins.add(protein_id)

    return proteins


# Reúne todas las salidas de PathoFact y construye las tablas canónicas de annotations y evidence
def parse_pathofact(
    *,
    features_path: Path,
    cluster_path: Path,
    vf_ml_path: Path,
    vf_hmm_path: Path,
    vf_high_path: Path,
    vf_ambiguous_path: Path,
    tox_ml_path: Path,
    tox_hmm_path: Path,
    tox_high_path: Path,
    tox_ambiguous_path: Path,
    amr_deeparg_path: Path,
    amr_amrfinder_path: Path,
    amr_rgi_path: Path,
    amr_summary_path: Path,
    sample_id: str,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    feature_index = load_protein_to_feature_index(
        features_path,
        sample_id,
    )

    clusters = load_cluster_map(cluster_path)

    # PathoFact trabaja con representantes de clusters, por eso después expandimos cada resultado
    # a todas las proteínas originales que representa
    annotations_by_key: dict[
        tuple[str, str, str],
        dict[str, Any],
    ] = {}

    evidence: list[dict[str, Any]] = []

    # Añade anotaciones canónicas sin duplicarlas y conserva el estado especializado más fuerte
    def add_annotation(
        *,
        protein_id: str,
        annotation_type: str,
        annotation_id: str,
        annotation_value: str = "",
        description: str = "",
    ) -> None:
        feature_uid = validate_protein(
            protein_id,
            feature_index,
        )

        key = (
            feature_uid,
            annotation_type,
            annotation_id,
        )

        existing = annotations_by_key.get(key)

        if existing is None:
            annotations_by_key[key] = {
                "sample_id": (sample_id),
                "feature_uid": (feature_uid),
                "protein_id": (protein_id),
                "source": SOURCE,
                "annotation_type": (annotation_type),
                "annotation_id": (annotation_id),
                "annotation_value": (annotation_value),
                "description": (description),
            }

            return

        # Preserve the strongest specialized status.
        if annotation_value == "high_confidence":
            existing["annotation_value"] = "high_confidence"

        if not existing["description"] and description:
            existing["description"] = description

    # Procesa los perfiles HMM de virulencia y toxinas y los convierte en anotaciones y evidence
    hmm_inputs = (
        (
            vf_hmm_path,
            "PathoFact2:VF-HMM",
        ),
        (
            tox_hmm_path,
            "PathoFact2:TOX-HMM",
        ),
    )

    for hmm_path, source in hmm_inputs:
        for hit in read_hmm_hits(hmm_path):
            representative = hit["protein_id"]

            annotation_type, (annotation_id) = canonical_signature(hit["accession"])

            for protein_id in expand_cluster(
                representative,
                clusters,
            ):
                feature_uid = validate_protein(
                    protein_id,
                    feature_index,
                )

                add_annotation(
                    protein_id=protein_id,
                    annotation_type=(annotation_type),
                    annotation_id=(annotation_id),
                    description=(hit["description"]),
                )

                evidence.append(
                    {
                        "sample_id": (sample_id),
                        "feature_uid": (feature_uid),
                        "protein_id": (protein_id),
                        "evidence_type": ("hmm_match"),
                        "source": source,
                        "accession": (hit["accession"]),
                        "gene": "",
                        "product": (hit["description"]),
                        "rank": "",
                        "score": (hit["score"]),
                        "evalue": (hit["evalue"]),
                        "identity": "",
                        "query_coverage": "",
                        "subject_coverage": "",
                        "aa_coverage": "",
                        "hmm_coverage": "",
                        "start": "",
                        "stop": "",
                    }
                )

    # Conserva como evidence solo las predicciones positivas de los modelos de virulencia y toxinas
    ml_inputs = (
        (
            vf_ml_path,
            "PathoFact2:VF-ML",
            "PATHOFACT2:VF",
            ("virulence-associated protein"),
        ),
        (
            tox_ml_path,
            "PathoFact2:TOX-ML",
            "PATHOFACT2:TOX",
            ("toxin-associated protein"),
        ),
    )

    for (
        ml_path,
        source,
        accession,
        product,
    ) in ml_inputs:
        for prediction in read_ml_predictions(ml_path):
            if prediction["prediction"] != "1":
                continue

            representative = prediction["protein_id"]

            for protein_id in expand_cluster(
                representative,
                clusters,
            ):
                feature_uid = validate_protein(
                    protein_id,
                    feature_index,
                )

                evidence.append(
                    {
                        "sample_id": (sample_id),
                        "feature_uid": (feature_uid),
                        "protein_id": (protein_id),
                        "evidence_type": ("ml_prediction"),
                        "source": source,
                        "accession": (accession),
                        "gene": "",
                        "product": product,
                        "rank": "",
                        "score": (prediction["probability"]),
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

    # Añade las clasificaciones finales de virulencia y toxinas como high_confidence o ambiguous
    classification_inputs = (
        (
            vf_high_path,
            "VIRULENCE_ASSOCIATED",
            "PATHOFACT2:VF",
            "high_confidence",
            ("PathoFact2 virulence factor prediction"),
        ),
        (
            vf_ambiguous_path,
            "VIRULENCE_ASSOCIATED",
            "PATHOFACT2:VF",
            "ambiguous",
            ("PathoFact2 virulence factor prediction"),
        ),
        (
            tox_high_path,
            "TOXIN_ASSOCIATED",
            "PATHOFACT2:TOX",
            "high_confidence",
            ("PathoFact2 toxin prediction"),
        ),
        (
            tox_ambiguous_path,
            "TOXIN_ASSOCIATED",
            "PATHOFACT2:TOX",
            "ambiguous",
            ("PathoFact2 toxin prediction"),
        ),
    )

    for (
        report_path,
        annotation_type,
        annotation_id,
        status,
        description,
    ) in classification_inputs:
        representatives = read_classification_report(report_path)

        for representative in representatives:
            for protein_id in expand_cluster(
                representative,
                clusters,
            ):
                add_annotation(
                    protein_id=protein_id,
                    annotation_type=(annotation_type),
                    annotation_id=(annotation_id),
                    annotation_value=(status),
                    description=(description),
                )

    # Reúne DeepARG, AMRFinderPlus y RGI por separado y usa el resumen AMR solo como comprobación
    deeparg_hits = read_deeparg_amr(amr_deeparg_path)

    amrfinder_hits = read_amrfinder_amr(amr_amrfinder_path)

    rgi_hits = read_rgi_amr(amr_rgi_path)

    source_amr_proteins = {hit["protein_id"] for hit in (deeparg_hits + amrfinder_hits + rgi_hits)}

    # El resumen integrado debe contener exactamente las mismas proteínas detectadas por los métodos AMR
    summary_amr_proteins = read_classification_report(amr_summary_path)

    if source_amr_proteins != summary_amr_proteins:
        raise ValueError("PathoFact2 AMR source outputs do not match the integrated AMR summary.")

    # Convierte los resultados de DeepARG en anotaciones ARO y evidence específico del método
    for hit in deeparg_hits:
        representative = hit["protein_id"]

        for protein_id in expand_cluster(
            representative,
            clusters,
        ):
            feature_uid = validate_protein(
                protein_id,
                feature_index,
            )

            if hit["aro"]:
                add_annotation(
                    protein_id=protein_id,
                    annotation_type="ARO",
                    annotation_id=hit["aro"],
                    description=("PathoFact2 antimicrobial resistance prediction"),
                )

            evidence.append(
                {
                    "sample_id": sample_id,
                    "feature_uid": feature_uid,
                    "protein_id": protein_id,
                    "evidence_type": "amr_prediction",
                    "source": "PathoFact2:DeepARG",
                    "accession": hit["aro"],
                    "gene": hit["arg"],
                    "product": hit["best_hit"],
                    "rank": hit["arg_class"],
                    "score": hit["probability"],
                    "evalue": hit["evalue"],
                    "identity": hit["identity"],
                    "query_coverage": "",
                    "subject_coverage": "",
                    "aa_coverage": "",
                    "hmm_coverage": "",
                    "start": hit["start"],
                    "stop": hit["stop"],
                }
            )

    # Convierte los resultados de AMRFinderPlus en anotaciones ARO y evidence específico del método
    for hit in amrfinder_hits:
        representative = hit["protein_id"]

        for protein_id in expand_cluster(
            representative,
            clusters,
        ):
            feature_uid = validate_protein(
                protein_id,
                feature_index,
            )

            if hit["aro"]:
                add_annotation(
                    protein_id=protein_id,
                    annotation_type="ARO",
                    annotation_id=hit["aro"],
                    description=("PathoFact2 antimicrobial resistance prediction"),
                )

            evidence.append(
                {
                    "sample_id": sample_id,
                    "feature_uid": feature_uid,
                    "protein_id": protein_id,
                    "evidence_type": "amr_prediction",
                    "source": "PathoFact2:AMRFinderPlus",
                    "accession": hit["aro"],
                    "gene": hit["gene_symbol"],
                    "product": hit["sequence_name"],
                    "rank": hit["method"],
                    "score": "",
                    "evalue": "",
                    "identity": hit["identity"],
                    "query_coverage": "",
                    "subject_coverage": (hit["reference_coverage"]),
                    "aa_coverage": "",
                    "hmm_coverage": "",
                    "start": "",
                    "stop": "",
                }
            )

    # Convierte los resultados de RGI/CARD en anotaciones ARO y evidence específico del método
    for hit in rgi_hits:
        representative = hit["protein_id"]

        for protein_id in expand_cluster(
            representative,
            clusters,
        ):
            feature_uid = validate_protein(
                protein_id,
                feature_index,
            )

            if hit["aro"]:
                add_annotation(
                    protein_id=protein_id,
                    annotation_type="ARO",
                    annotation_id=hit["aro"],
                    description=("PathoFact2 antimicrobial resistance prediction"),
                )

            evidence.append(
                {
                    "sample_id": sample_id,
                    "feature_uid": feature_uid,
                    "protein_id": protein_id,
                    "evidence_type": "amr_prediction",
                    "source": "PathoFact2:RGI",
                    "accession": hit["aro"],
                    "gene": hit["best_hit"],
                    "product": hit["amr_gene_family"],
                    "rank": hit["cutoff"],
                    "score": hit["best_hit_bitscore"],
                    "evalue": "",
                    "identity": hit["identity"],
                    "query_coverage": "",
                    "subject_coverage": (hit["reference_coverage"]),
                    "aa_coverage": "",
                    "hmm_coverage": "",
                    "start": hit["hit_start"],
                    "stop": hit["hit_end"],
                }
            )

    annotations = list(annotations_by_key.values())

    # Ordena ambas tablas para obtener salidas deterministas y fáciles de comparar
    annotations.sort(
        key=lambda row: (
            row["protein_id"],
            row["annotation_type"],
            row["annotation_id"],
        )
    )

    evidence.sort(
        key=lambda row: (
            row["protein_id"],
            row["source"],
            row["accession"],
            row["score"],
        )
    )

    return (
        annotations,
        evidence,
    )


# Escribe las filas normalizadas usando las columnas de la tabla canónica indicada
def write_tsv(
    rows: list[dict[str, Any]],
    columns: list[str],
    output_path: Path,
) -> None:
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


# Ejecuta el parser completo de PathoFact y escribe annotations y evidence
def parse_pathofact_bundle(
    *,
    features_path: Path,
    cluster_path: Path,
    vf_ml_path: Path,
    vf_hmm_path: Path,
    vf_high_path: Path,
    vf_ambiguous_path: Path,
    tox_ml_path: Path,
    tox_hmm_path: Path,
    tox_high_path: Path,
    tox_ambiguous_path: Path,
    amr_deeparg_path: Path,
    amr_amrfinder_path: Path,
    amr_rgi_path: Path,
    amr_summary_path: Path,
    sample_id: str,
    annotations_output: Path,
    evidence_output: Path,
) -> None:
    annotations, evidence = parse_pathofact(
        features_path=features_path,
        cluster_path=cluster_path,
        vf_ml_path=vf_ml_path,
        vf_hmm_path=vf_hmm_path,
        vf_high_path=vf_high_path,
        vf_ambiguous_path=(vf_ambiguous_path),
        tox_ml_path=tox_ml_path,
        tox_hmm_path=tox_hmm_path,
        tox_high_path=tox_high_path,
        tox_ambiguous_path=tox_ambiguous_path,
        amr_deeparg_path=amr_deeparg_path,
        amr_amrfinder_path=amr_amrfinder_path,
        amr_rgi_path=amr_rgi_path,
        amr_summary_path=amr_summary_path,
        sample_id=sample_id,
    )

    write_tsv(
        annotations,
        ANNOTATION_COLUMNS,
        annotations_output,
    )

    write_tsv(
        evidence,
        EVIDENCE_COLUMNS,
        evidence_output,
    )
