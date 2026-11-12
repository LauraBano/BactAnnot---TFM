"""Canonical BactAnnot table schemas"""

ANNOTATION_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "source",
    "annotation_type",
    "annotation_id",
    "annotation_value",
    "description",
]


EVIDENCE_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "evidence_type",
    "source",
    "accession",
    "gene",
    "product",
    "rank",
    "score",
    "native_threshold",
    "evalue",
    "identity",
    "query_coverage",
    "subject_coverage",
    "aa_coverage",
    "hmm_coverage",
    "start",
    "stop",
]


KO_CLAIM_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "claim_type",
    "claim_id",
    "observing_sources",
    "observation_count",
]


KO_EVIDENCE_CLAIM_COLUMNS = [
    *KO_CLAIM_COLUMNS,
    "evidence_sources",
    "evidence_lines",
    "independent_evidence_count",
    "context_sources",
    "unresolved_sources",
    "evidence_status",
]


KO_EGGNOG_KOFAM_COMPARISON_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "eggnog_kos",
    "kofam_kos",
    "shared_kos",
    "eggnog_only_kos",
    "kofam_only_kos",
    "eggnog_kofam_relation",
]


EC_CLAIM_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "claim_type",
    "claim_id",
    "observing_sources",
    "observation_count",
]


EC_EVIDENCE_CLAIM_COLUMNS = [
    *EC_CLAIM_COLUMNS,
    "evidence_sources",
    "evidence_lines",
    "evidence_relations",
    "independent_evidence_count",
    "unresolved_sources",
    "evidence_status",
]


EC_SOURCE_COMPARISON_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "source_a",
    "source_b",
    "source_a_ecs",
    "source_b_ecs",
    "matched_ec_pairs",
    "source_a_unmatched_ecs",
    "source_b_unmatched_ecs",
    "relation",
]


GO_CLAIM_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "claim_type",
    "claim_id",
    "observing_sources",
    "observation_count",
]


GO_EVIDENCE_CLAIM_COLUMNS = [
    *GO_CLAIM_COLUMNS,
    "evidence_sources",
    "evidence_lines",
    "evidence_relations",
    "independent_evidence_count",
    "context_sources",
    "unresolved_sources",
    "evidence_status",
]


GO_SOURCE_COMPARISON_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "aspect",
    "source_a",
    "source_b",
    "source_a_go_ids",
    "source_b_go_ids",
    "source_a_unresolved_go_ids",
    "source_b_unresolved_go_ids",
    "matched_go_pairs",
    "source_a_unmatched_go_ids",
    "source_b_unmatched_go_ids",
    "relation",
]


DOMAIN_CLAIM_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "claim_type",
    "claim_id",
    "observing_sources",
    "observation_count",
]


DOMAIN_EVIDENCE_CLAIM_COLUMNS = [
    *DOMAIN_CLAIM_COLUMNS,
    "evidence_sources",
    "evidence_lines",
    "independent_evidence_count",
    "context_sources",
    "unresolved_sources",
    "evidence_status",
]


COG_CLAIM_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "claim_type",
    "claim_id",
    "observing_sources",
    "observation_count",
]


COG_EVIDENCE_CLAIM_COLUMNS = [
    *COG_CLAIM_COLUMNS,
    "evidence_sources",
    "evidence_lines",
    "independent_evidence_count",
    "context_sources",
    "unresolved_sources",
    "evidence_status",
]


ARO_CLAIM_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "claim_type",
    "claim_id",
    "observing_sources",
    "observation_count",
]


ARO_EVIDENCE_CLAIM_COLUMNS = [
    *ARO_CLAIM_COLUMNS,
    "evidence_sources",
    "evidence_lines",
    "independent_evidence_count",
    "context_sources",
    "unresolved_sources",
    "evidence_status",
]


BIOMEDICAL_CLASSIFICATION_CLAIM_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "claim_type",
    "claim_id",
    "classification_status",
    "observing_sources",
    "observation_count",
]


BIOMEDICAL_CLASSIFICATION_EVIDENCE_CLAIM_COLUMNS = [
    *BIOMEDICAL_CLASSIFICATION_CLAIM_COLUMNS,
    "evidence_sources",
    "evidence_lines",
    "independent_evidence_count",
    "context_sources",
    "unresolved_sources",
    "evidence_status",
]


KO_EC_RELATION_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "ko_id",
    "ec_id",
    "relation",
    "reference_id",
]


ARO_EVIDENCE_CLAIM_COLUMNS = [
    *ARO_CLAIM_COLUMNS,
    "evidence_sources",
    "evidence_lines",
    "evidence_relations",
    "independent_evidence_count",
    "context_sources",
    "unresolved_sources",
    "evidence_status",
]


CONCORDANCE_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "annotation_type",
    "source_1",
    "source_2",
    "source_1_ids",
    "source_2_ids",
    "shared_ids",
    "relation",
]


CONCORDANCE_SUMMARY_COLUMNS = [
    "sample_id",
    "annotation_type",
    "source_1",
    "source_2",
    "both",
    "exact",
    "partial",
    "disjoint",
    "only_source_1",
    "only_source_2",
    "total",
]
