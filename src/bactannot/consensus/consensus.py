"""Evidence-aware consensus orchestration and table export"""

from __future__ import annotations

import csv
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from bactannot.consensus.biomedical import (
    build_aro_claims_with_evidence,
    build_biomedical_classification_claims_with_evidence,
)
from bactannot.consensus.cog import build_cog_claims_with_evidence
from bactannot.consensus.domains import build_domain_claims_with_evidence
from bactannot.consensus.ec import (
    build_ec_claims_with_evidence,
    build_ec_source_comparisons,
)
from bactannot.consensus.go import (
    build_go_claims_with_evidence,
    build_go_comparison_view,
    build_go_source_comparisons,
)
from bactannot.consensus.go_ontology import (
    build_go_ancestor_map,
    get_bundled_go_ontology,
    load_go_ontology,
)
from bactannot.consensus.interpretation import (
    build_enzyme_activity_candidates,
    build_molecular_function_candidates,
    build_protein_function_interpretations,
    collect_claim_evidence_pairs,
    select_convergent_ec_activity_proposals,
    select_convergent_family_function_proposals,
    select_convergent_molecular_function_proposals,
    select_single_ko_function_proposals,
)
from bactannot.consensus.ko import (
    build_eggnog_kofam_comparisons,
    build_ko_claims_with_evidence,
)
from bactannot.consensus.policy import validate_annotation_type_coverage
from bactannot.consensus.relations import (
    build_ec_go_relations,
    build_ko_ec_relations,
    load_kofam_ko_ec_map,
)
from bactannot.integration.annotations import read_annotation_table
from bactannot.integration.evidence import read_evidence_table
from bactannot.schema import (
    ARO_EVIDENCE_CLAIM_COLUMNS,
    BIOMEDICAL_CLASSIFICATION_EVIDENCE_CLAIM_COLUMNS,
    COG_EVIDENCE_CLAIM_COLUMNS,
    DOMAIN_EVIDENCE_CLAIM_COLUMNS,
    EC_EVIDENCE_CLAIM_COLUMNS,
    EC_SOURCE_COMPARISON_COLUMNS,
    GO_EVIDENCE_CLAIM_COLUMNS,
    GO_SOURCE_COMPARISON_COLUMNS,
    KO_EC_RELATION_COLUMNS,
    KO_EGGNOG_KOFAM_COMPARISON_COLUMNS,
    KO_EVIDENCE_CLAIM_COLUMNS,
)

# Define las columnas mínimas de features.tsv que necesitamos para validar el consenso
REQUIRED_FEATURE_COLUMNS = {
    "sample_id",
    "feature_uid",
    "protein_id",
}

# Asocia cada salida con sus columnas existentes
EVIDENCE_CONSENSUS_TABLE_COLUMNS = {
    "ko_claims": KO_EVIDENCE_CLAIM_COLUMNS,
    "ec_claims": EC_EVIDENCE_CLAIM_COLUMNS,
    "go_claims": GO_EVIDENCE_CLAIM_COLUMNS,
    "domain_claims": DOMAIN_EVIDENCE_CLAIM_COLUMNS,
    "cog_claims": COG_EVIDENCE_CLAIM_COLUMNS,
    "aro_claims": ARO_EVIDENCE_CLAIM_COLUMNS,
    "biomedical_claims": BIOMEDICAL_CLASSIFICATION_EVIDENCE_CLAIM_COLUMNS,
    "ko_ec_relations": KO_EC_RELATION_COLUMNS,
    "ko_comparisons": KO_EGGNOG_KOFAM_COMPARISON_COLUMNS,
    "ec_comparisons": EC_SOURCE_COMPARISON_COLUMNS,
    "go_comparisons": [
        *GO_SOURCE_COMPARISON_COLUMNS,
        "reference_id",
    ],
    "ec_go_relations": [
        "sample_id",
        "feature_uid",
        "protein_id",
        "ec_id",
        "go_id",
        "resolved_go_id",
        "go_resolution_status",
        "relation",
        "reference_id",
    ],
    "protein_interpretations": [
        "sample_id",
        "feature_uid",
        "protein_id",
        "interpretation_status",
        "primary_claim_type",
        "primary_claim_id",
        "primary_function",
        "primary_scope",
        "primary_selection_rule",
        "primary_description_source",
        "primary_description_reference_id",
        "ko_proposal_ids",
        "ec_activity_proposal_ids",
        "mf_proposal_ids",
        "unlinked_ec_activity_ids",
        "family_proposal_ids",
    ],
}


# Carga la relación feature_uid -> protein_id definida por Bakta para esta muestra
def load_feature_to_protein_index(
    features_path: Path,
    sample_id: str,
) -> dict[str, str]:
    """Map canonical feature UIDs to their protein identifiers"""

    feature_index: dict[str, str] = {}

    with features_path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(f"Missing header in features table: {features_path}")

        missing = REQUIRED_FEATURE_COLUMNS.difference(reader.fieldnames)

        if missing:
            raise ValueError(
                f"Missing required columns in {features_path}: {', '.join(sorted(missing))}"
            )

        for row in reader:
            row_sample_id = row["sample_id"].strip()
            feature_uid = row["feature_uid"].strip()
            protein_id = row["protein_id"].strip()

            # Todas las features deben pertenecer a la muestra que estamos procesando
            if row_sample_id != sample_id:
                raise ValueError(
                    f"Unexpected sample_id in {features_path}: {row_sample_id} != {sample_id}"
                )

            if not feature_uid:
                raise ValueError(f"Missing feature_uid in {features_path}")

            # Cada feature_uid debe identificar una única feature
            if feature_uid in feature_index:
                raise ValueError(f"Duplicate feature_uid in features table: {feature_uid}")

            feature_index[feature_uid] = protein_id

    return feature_index


# Escribe la tabla final de consenso respetando el esquema canónico
def write_consensus_table(
    rows: list[dict[str, str]],
    output_path: Path,
    columns: list[str],
) -> None:
    """Write one consensus table using its explicit schema."""

    expected_columns = set(columns)

    for row in rows:
        if set(row) != expected_columns:
            raise ValueError(f"Unexpected consensus schema for {output_path}")

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


def build_evidence_claim_tables(
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    """Build the implemented namespace-specific claim tables."""

    validate_annotation_type_coverage(annotation_rows)

    return {
        "ko_claims": build_ko_claims_with_evidence(annotation_rows, evidence_rows),
        "ec_claims": build_ec_claims_with_evidence(annotation_rows, evidence_rows),
        "go_claims": build_go_claims_with_evidence(annotation_rows, evidence_rows),
        "domain_claims": build_domain_claims_with_evidence(annotation_rows, evidence_rows),
        "cog_claims": build_cog_claims_with_evidence(annotation_rows, evidence_rows),
        "aro_claims": build_aro_claims_with_evidence(annotation_rows, evidence_rows),
        "biomedical_claims": build_biomedical_classification_claims_with_evidence(
            annotation_rows,
            evidence_rows,
        ),
    }


def build_evidence_consensus_tables(
    annotations_path: Path,
    evidence_path: Path,
    features_path: Path,
    sample_id: str,
    ko_list_path: Path,
    output_dir: Path,
    go_ontology_path: Path | None = None,
) -> dict[str, Path]:
    """Build and write evidence-aware claims, relations and source comparisons"""

    feature_index = load_feature_to_protein_index(
        features_path,
        sample_id,
    )

    annotation_rows = read_annotation_table(
        annotations_path,
        sample_id,
        feature_index,
    )

    evidence_rows = read_evidence_table(
        evidence_path,
        sample_id,
        feature_index,
    )

    ko_ec_map = load_kofam_ko_ec_map(ko_list_path)

    reference_id = f"KOfam:ko_list:sha256:{sha256(ko_list_path.read_bytes()).hexdigest()}"

    tables = build_evidence_claim_tables(
        annotation_rows,
        evidence_rows,
    )

    tables["ko_ec_relations"] = build_ko_ec_relations(
        tables["ko_claims"],
        tables["ec_claims"],
        ko_ec_map,
        reference_id,
    )

    protein_annotation_rows = [row for row in annotation_rows if row["protein_id"].strip()]

    protein_claim_rows = [
        claim
        for name in (
            "ko_claims",
            "ec_claims",
            "go_claims",
            "domain_claims",
            "cog_claims",
            "aro_claims",
            "biomedical_claims",
        )
        for claim in tables[name]
        if claim["protein_id"].strip()
    ]

    # Comparaciones KO
    tables["ko_comparisons"] = build_eggnog_kofam_comparisons(
        protein_annotation_rows,
    )

    # Comparaciones EC
    ec_sources = sorted(
        {row["source"] for row in protein_annotation_rows if row["annotation_type"] == "EC"},
        key=lambda source: (
            source.casefold(),
            source,
        ),
    )

    tables["ec_comparisons"] = []

    for source_a, source_b in combinations(ec_sources, 2):
        tables["ec_comparisons"].extend(
            build_ec_source_comparisons(
                protein_annotation_rows,
                source_a,
                source_b,
            )
        )

    # Comparaciones GO
    go_annotation_rows = [row for row in protein_annotation_rows if row["annotation_type"] == "GO"]

    tables["go_comparisons"] = []
    tables["ec_go_relations"] = []

    mf_proposals: list[dict[str, str]] = []

    if go_annotation_rows:
        ontology_path = go_ontology_path

        if ontology_path is None:
            ontology_path = get_bundled_go_ontology()

        terms, alt_to_primary = load_go_ontology(ontology_path)
        ancestor_map = build_go_ancestor_map(terms)

        ontology_sha256 = sha256(ontology_path.read_bytes()).hexdigest()
        go_reference_id = f"GO:obo:sha256:{ontology_sha256}"

        tables["ec_go_relations"] = build_ec_go_relations(
            tables["ec_claims"],
            tables["go_claims"],
            terms,
            alt_to_primary,
            go_reference_id,
        )

        go_comparison_view = build_go_comparison_view(
            go_annotation_rows,
            terms,
            alt_to_primary,
            ancestor_map,
        )

        go_sources = sorted(
            {row["source"] for row in go_annotation_rows},
            key=lambda source: (
                source.casefold(),
                source,
            ),
        )

        for source_a, source_b in combinations(go_sources, 2):
            tables["go_comparisons"].extend(
                build_go_source_comparisons(
                    go_comparison_view,
                    ancestor_map,
                    source_a,
                    source_b,
                )
            )

        for row in tables["go_comparisons"]:
            row["reference_id"] = go_reference_id

        mf_candidates = build_molecular_function_candidates(
            [claim for claim in protein_claim_rows if claim["claim_type"] == "GO"],
            terms,
            alt_to_primary,
            go_reference_id,
        )

        mf_proposals = select_convergent_molecular_function_proposals(
            mf_candidates,
            terms,
        )

    # Interpretación funcional por proteína
    claim_evidence_pairs = collect_claim_evidence_pairs(
        protein_claim_rows,
        protein_annotation_rows,
        evidence_rows,
    )

    ko_proposals = select_single_ko_function_proposals(
        claim_evidence_pairs,
    )

    activity_candidates = build_enzyme_activity_candidates(
        [claim for claim in protein_claim_rows if claim["claim_type"] == "EC"],
        tables["ko_ec_relations"],
        tables["ec_go_relations"],
    )

    activity_proposals = select_convergent_ec_activity_proposals(
        activity_candidates,
        claim_evidence_pairs,
    )

    family_proposals = select_convergent_family_function_proposals(
        claim_evidence_pairs,
        evidence_rows,
    )

    tables["protein_interpretations"] = build_protein_function_interpretations(
        sample_id,
        feature_index,
        ko_proposals,
        activity_proposals,
        mf_proposals,
        family_proposals=family_proposals,
    )

    # Escritura de las tablas
    output_paths: dict[str, Path] = {}

    for name, columns in EVIDENCE_CONSENSUS_TABLE_COLUMNS.items():
        output_path = output_dir / f"{sample_id}.{name}.tsv"

        write_consensus_table(
            rows=tables[name],
            output_path=output_path,
            columns=columns,
        )

        output_paths[name] = output_path

    return output_paths
