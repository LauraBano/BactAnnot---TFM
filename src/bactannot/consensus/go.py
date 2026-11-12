"""Canonical GO claims for evidence-aware consensus"""

import json
from collections import defaultdict

from bactannot.consensus.go_ontology import (
    GoTerm,
    reduce_to_most_specific,
    resolve_go_id,
)
from bactannot.schema import (
    GO_CLAIM_COLUMNS,
    GO_EVIDENCE_CLAIM_COLUMNS,
    GO_SOURCE_COMPARISON_COLUMNS,
)

GO_COMPARISON_VIEW_COLUMNS = [
    "sample_id",
    "feature_uid",
    "protein_id",
    "source",
    "aspect",
    "specific_go_ids",
    "ontology_unresolved_go_ids",
]


def serialize_sources(
    sources: set[str],
) -> str:
    """Serialize distinct annotation sources deterministically"""

    return ";".join(
        sorted(
            sources,
            key=str.casefold,
        )
    )


GO_CONTEXT_SOURCES = {
    "bakta",
}

EGGNOG_GO_EVIDENCE_LINE = (
    "orthology_transfer",
    "eggNOG",
)

INTERPROSCAN_GO_EVIDENCE_LINE = (
    "signature_mapping",
    "InterProScan",
)


def serialize_evidence_lines(
    evidence_lines: set[tuple[str, str]],
) -> str:
    """Serialize distinct GO evidence lines deterministically"""

    return ";".join(
        f"{method}@{resource}"
        for method, resource in sorted(
            evidence_lines,
            key=lambda item: (
                item[0].casefold(),
                item[1].casefold(),
            ),
        )
    )


def serialize_relations(
    relations: set[str],
) -> str:
    """Serialize GO evidence relations deterministically"""

    return ";".join(
        sorted(
            relations,
            key=str.casefold,
        )
    )


def build_go_claims_from_annotations(
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build canonical GO claims from annotation observations"""

    observations_by_claim: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in annotation_rows:
        if row["annotation_type"] != "GO":
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            row["annotation_id"],
        )

        observations_by_claim[key].append(row)

    claims: list[dict[str, str]] = []

    for key, observations in observations_by_claim.items():
        (
            sample_id,
            feature_uid,
            claim_id,
        ) = key

        protein_ids = {row["protein_id"] for row in observations}

        if len(protein_ids) != 1:
            raise TypeError(f"Conflicting protein_id values for GO claim: {key}")

        protein_id = next(iter(protein_ids))

        observing_sources = {row["source"] for row in observations}

        claims.append(
            {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "claim_type": "GO",
                "claim_id": claim_id,
                "observing_sources": (serialize_sources(observing_sources)),
                "observation_count": str(len(observations)),
            }
        )

    claims.sort(
        key=lambda claim: (
            claim["sample_id"],
            claim["feature_uid"],
            claim["claim_id"],
        )
    )

    for claim in claims:
        if list(claim) != GO_CLAIM_COLUMNS:
            raise TypeError("GO claim does not match the canonical GO claim schema")

    return claims


def parse_interproscan_go_origins(
    annotation_value: str,
) -> list[dict[str, str]]:
    """Parse canonical InterProScan GO provenance"""

    try:
        payload = json.loads(annotation_value)

    except json.JSONDecodeError as exc:
        raise TypeError("Invalid InterProScan GO provenance JSON") from exc

    if not isinstance(
        payload,
        list,
    ):
        raise TypeError("InterProScan GO provenance must be a list")

    if not payload:
        raise TypeError("InterProScan GO provenance cannot be empty")

    required_fields = {
        "analysis",
        "signature_accession",
        "interpro_accession",
        "mapping_source",
    }

    origins: list[dict[str, str]] = []

    for origin in payload:
        if not isinstance(
            origin,
            dict,
        ):
            raise TypeError("InterProScan GO provenance origin must be an object")

        if set(origin) != required_fields:
            raise TypeError("Unexpected InterProScan GO provenance fields")

        origins.append({field: str(origin[field]) for field in required_fields})

    return origins


def collect_eggnog_go_evidence(
    claim: dict[str, str],
    observations: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Collect eggNOG seed-ortholog evidence for one GO claim"""

    if claim["claim_type"] != "GO":
        raise TypeError(f"Expected a GO claim, got {claim['claim_type']!r}")

    has_eggnog_observation = any(row["source"] == "eggNOG-mapper" for row in observations)

    if not has_eggnog_observation:
        return []

    matches = [
        row
        for row in evidence_rows
        if (
            row["sample_id"] == claim["sample_id"]
            and row["feature_uid"] == claim["feature_uid"]
            and row["protein_id"] == claim["protein_id"]
            and row["source"] == "eggNOG-mapper"
            and row["evidence_type"] == "seed_ortholog"
        )
    ]

    matches.sort(
        key=lambda row: (
            row["source"],
            row["evidence_type"],
            row["accession"],
        )
    )

    return matches


def collect_interproscan_go_evidence(
    claim: dict[str, str],
    observations: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Collect signature evidence explicitly linked to one GO claim"""

    if claim["claim_type"] != "GO":
        raise TypeError(f"Expected a GO claim, got {claim['claim_type']!r}")

    interpro_observations = [
        row
        for row in observations
        if (
            row["source"] == "InterProScan"
            and row["annotation_type"] == "GO"
            and row["annotation_id"] == claim["claim_id"]
        )
    ]

    if not interpro_observations:
        return []

    expected_evidence_keys: set[tuple[str, str, str]] = set()

    for observation in interpro_observations:
        origins = parse_interproscan_go_origins(observation["annotation_value"])

        for origin in origins:
            expected_evidence_keys.add(
                (
                    (f"InterProScan:{origin['analysis']}"),
                    "signature_match",
                    origin["signature_accession"],
                )
            )

    matches = [
        row
        for row in evidence_rows
        if (
            row["sample_id"] == claim["sample_id"]
            and row["feature_uid"] == claim["feature_uid"]
            and row["protein_id"] == claim["protein_id"]
            and (
                row["source"],
                row["evidence_type"],
                row["accession"],
            )
            in expected_evidence_keys
        )
    ]

    matched_keys = {
        (
            row["source"],
            row["evidence_type"],
            row["accession"],
        )
        for row in matches
    }

    missing_keys = expected_evidence_keys - matched_keys

    if missing_keys:
        raise TypeError(
            "InterProScan GO provenance "
            "cannot be resolved to signature evidence: "
            f"{sorted(missing_keys)}"
        )

    matches.sort(
        key=lambda row: (
            row["source"],
            row["evidence_type"],
            row["accession"],
        )
    )

    return matches


def summarize_go_claim_evidence(
    claim: dict[str, str],
    observations: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> dict[str, str]:
    """Summarize methodological evidence supporting one GO claim"""

    eggnog_evidence = collect_eggnog_go_evidence(
        claim,
        observations,
        evidence_rows,
    )

    interproscan_evidence = collect_interproscan_go_evidence(
        claim,
        observations,
        evidence_rows,
    )

    evidence_sources: set[str] = set()

    evidence_lines: set[tuple[str, str]] = set()

    evidence_relations: set[str] = set()

    if eggnog_evidence:
        evidence_sources.add("eggNOG-mapper")

        evidence_lines.add(EGGNOG_GO_EVIDENCE_LINE)

        evidence_relations.add("transferred_annotation")

    if interproscan_evidence:
        evidence_sources.add("InterProScan")

        evidence_lines.add(INTERPROSCAN_GO_EVIDENCE_LINE)

        evidence_relations.add("mapped_from_signature")

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    context_sources = (observing_sources & GO_CONTEXT_SOURCES) - evidence_sources

    unresolved_sources = observing_sources - evidence_sources - context_sources

    return {
        "evidence_sources": (serialize_sources(evidence_sources)),
        "evidence_lines": (serialize_evidence_lines(evidence_lines)),
        "evidence_relations": (serialize_relations(evidence_relations)),
        "independent_evidence_count": str(len(evidence_lines)),
        "context_sources": (serialize_sources(context_sources)),
        "unresolved_sources": (serialize_sources(unresolved_sources)),
    }


def classify_go_evidence_status(
    evidence: dict[str, str],
) -> str:
    """Classify the evidence status of one GO claim"""

    independent_count = int(evidence["independent_evidence_count"])

    if independent_count >= 2:
        return "multiple_independent"

    if independent_count == 1:
        return "single_independent"

    if evidence["unresolved_sources"]:
        return "unresolved"

    if evidence["context_sources"]:
        return "context_only"

    return "no_evidence"


def build_go_claims_with_evidence(
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build GO claims enriched with methodological evidence"""

    claims = build_go_claims_from_annotations(annotation_rows)

    observations_by_claim: dict[
        tuple[
            str,
            str,
            str,
            str,
        ],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in annotation_rows:
        if row["annotation_type"] != "GO":
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["annotation_id"],
        )

        observations_by_claim[key].append(row)

    evidence_by_feature: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in evidence_rows:
        key = (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
        )

        evidence_by_feature[key].append(row)

    evidence_claims: list[dict[str, str]] = []

    for claim in claims:
        claim_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
            claim["claim_id"],
        )

        feature_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
        )

        observations = observations_by_claim.get(
            claim_key,
            [],
        )

        candidate_evidence = evidence_by_feature.get(
            feature_key,
            [],
        )

        evidence = summarize_go_claim_evidence(
            claim,
            observations,
            candidate_evidence,
        )

        evidence_claim = {
            **claim,
            **evidence,
            "evidence_status": (classify_go_evidence_status(evidence)),
        }

        if list(evidence_claim) != GO_EVIDENCE_CLAIM_COLUMNS:
            raise TypeError("GO evidence claim does not match the canonical schema")

        evidence_claims.append(evidence_claim)

    return evidence_claims


def serialize_go_ids(
    go_ids: set[str],
) -> str:
    """Serialize distinct GO identifiers deterministically"""

    return ";".join(
        sorted(
            go_ids,
            key=str.casefold,
        )
    )


def build_go_comparison_view(
    annotation_rows: list[dict[str, str]],
    terms: dict[
        str,
        GoTerm,
    ],
    alt_to_primary: dict[
        str,
        str,
    ],
    ancestor_map: dict[
        str,
        frozenset[str],
    ],
) -> list[dict[str, str]]:
    """Build ontology-normalized GO sets for source comparison"""

    active_by_unit: dict[
        tuple[
            str,
            str,
            str,
            str,
            str,
        ],
        set[str],
    ] = defaultdict(set)

    unresolved_by_unit: dict[
        tuple[
            str,
            str,
            str,
            str,
            str,
        ],
        set[str],
    ] = defaultdict(set)

    all_units: set[
        tuple[
            str,
            str,
            str,
            str,
            str,
        ]
    ] = set()

    for row in annotation_rows:
        if row["annotation_type"] != "GO":
            continue

        resolution = resolve_go_id(
            row["annotation_id"],
            terms,
            alt_to_primary,
        )

        key = (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["source"],
            resolution.namespace,
        )

        all_units.add(key)

        if resolution.hierarchy_eligible:
            active_by_unit[key].add(resolution.resolved_id)

        else:
            unresolved_by_unit[key].add(resolution.resolved_id)

    output: list[dict[str, str]] = []

    for key in sorted(
        all_units,
        key=lambda item: tuple(value.casefold() for value in item),
    ):
        (
            sample_id,
            feature_uid,
            protein_id,
            source,
            aspect,
        ) = key

        active_go_ids = active_by_unit.get(
            key,
            set(),
        )

        specific_go_ids = reduce_to_most_specific(
            active_go_ids,
            ancestor_map,
        )

        unresolved_go_ids = unresolved_by_unit.get(
            key,
            set(),
        )

        row = {
            "sample_id": sample_id,
            "feature_uid": feature_uid,
            "protein_id": protein_id,
            "source": source,
            "aspect": aspect,
            "specific_go_ids": (serialize_go_ids(specific_go_ids)),
            "ontology_unresolved_go_ids": (serialize_go_ids(unresolved_go_ids)),
        }

        if list(row) != GO_COMPARISON_VIEW_COLUMNS:
            raise ValueError("GO comparison view row does not match the expected schema")

        output.append(row)

    return output


def classify_go_term_relation(
    go_a: str,
    go_b: str,
    ancestor_map: dict[
        str,
        frozenset[str],
    ],
    unresolved_a: set[str],
    unresolved_b: set[str],
) -> str:
    """Classify the ontology relation between two GO identifiers"""

    if go_a == go_b:
        return "exact"

    if go_a in unresolved_a or go_b in unresolved_b:
        return "disjoint"

    if go_a in ancestor_map.get(
        go_b,
        frozenset(),
    ):
        return "source_a_ancestor"

    if go_b in ancestor_map.get(
        go_a,
        frozenset(),
    ):
        return "source_a_descendant"

    return "disjoint"


def maximum_go_matching(
    source_a_specific: set[str],
    source_a_unresolved: set[str],
    source_b_specific: set[str],
    source_b_unresolved: set[str],
    ancestor_map: dict[
        str,
        frozenset[str],
    ],
) -> list[tuple[str, str, str]]:
    """Build a deterministic maximum matching between compatible GO terms"""

    source_a_ids = source_a_specific | source_a_unresolved

    source_b_ids = source_b_specific | source_b_unresolved

    adjacency: dict[
        str,
        list[tuple[str, str]],
    ] = {}

    relation_priority = {
        "exact": 0,
        "source_a_ancestor": 1,
        "source_a_descendant": 1,
    }

    for go_a in source_a_ids:
        candidates = []

        for go_b in source_b_ids:
            relation = classify_go_term_relation(
                go_a,
                go_b,
                ancestor_map,
                source_a_unresolved,
                source_b_unresolved,
            )

            if relation == "disjoint":
                continue

            candidates.append(
                (
                    go_b,
                    relation,
                )
            )

        candidates.sort(
            key=lambda item: (
                relation_priority[item[1]],
                item[0],
            )
        )

        adjacency[go_a] = candidates

    matched_b: dict[
        str,
        str,
    ] = {}

    def augment(
        go_a: str,
        visited_b: set[str],
    ) -> bool:
        for go_b, _ in adjacency[go_a]:
            if go_b in visited_b:
                continue

            visited_b.add(go_b)

            current_a = matched_b.get(go_b)

            if current_a is None or augment(
                current_a,
                visited_b,
            ):
                matched_b[go_b] = go_a

                return True

        return False

    for go_a in sorted(source_a_ids):
        augment(
            go_a,
            set(),
        )

    matches = []

    for go_b, go_a in matched_b.items():
        relation = classify_go_term_relation(
            go_a,
            go_b,
            ancestor_map,
            source_a_unresolved,
            source_b_unresolved,
        )

        matches.append(
            (
                go_a,
                go_b,
                relation,
            )
        )

    matches.sort(
        key=lambda item: (
            item[0],
            item[1],
            item[2],
        )
    )

    return matches


def classify_go_sets(
    source_a_specific: set[str],
    source_a_unresolved: set[str],
    source_b_specific: set[str],
    source_b_unresolved: set[str],
    ancestor_map: dict[
        str,
        frozenset[str],
    ],
) -> str:
    """Classify two non-empty GO comparison sets"""

    source_a_ids = source_a_specific | source_a_unresolved

    source_b_ids = source_b_specific | source_b_unresolved

    if source_a_specific == source_b_specific and source_a_unresolved == source_b_unresolved:
        return "exact"

    matches = maximum_go_matching(
        source_a_specific,
        source_a_unresolved,
        source_b_specific,
        source_b_unresolved,
        ancestor_map,
    )

    matched_count = len(matches)

    if matched_count == len(source_a_ids) == len(source_b_ids):
        return "hierarchical"

    if matched_count:
        return "partial"

    return "disjoint"


def parse_serialized_go_ids(
    value: str,
) -> set[str]:
    """Parse a semicolon-separated GO identifier set"""

    return {go_id for go_id in value.split(";") if go_id}


def serialize_go_matches(
    matches: list[tuple[str, str, str]],
) -> str:
    """Serialize matched GO identifier pairs deterministically"""

    return ";".join(f"{go_a}~{go_b}" for go_a, go_b, _ in matches)


def build_go_source_comparisons(
    comparison_view: list[dict[str, str]],
    ancestor_map: dict[
        str,
        frozenset[str],
    ],
    source_a: str,
    source_b: str,
) -> list[dict[str, str]]:
    """Compare ontology-normalized GO sets from two sources"""

    if source_a == source_b:
        raise ValueError("GO comparison requires two distinct sources")

    rows_by_unit: dict[
        tuple[
            str,
            str,
            str,
            str,
        ],
        dict[
            str,
            dict[str, str],
        ],
    ] = defaultdict(dict)

    for row in comparison_view:
        source = row["source"]

        if source not in {
            source_a,
            source_b,
        }:
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["aspect"],
        )

        if source in rows_by_unit[key]:
            raise ValueError(f"Duplicate GO comparison-view row for {key} / {source}")

        rows_by_unit[key][source] = row

    comparisons: list[dict[str, str]] = []

    for key in sorted(
        rows_by_unit,
        key=lambda item: tuple(value.casefold() for value in item),
    ):
        (
            sample_id,
            feature_uid,
            protein_id,
            aspect,
        ) = key

        source_rows = rows_by_unit[key]

        row_a = source_rows.get(source_a)

        row_b = source_rows.get(source_b)

        specific_a = parse_serialized_go_ids(row_a["specific_go_ids"]) if row_a else set()

        unresolved_a = (
            parse_serialized_go_ids(row_a["ontology_unresolved_go_ids"]) if row_a else set()
        )

        specific_b = parse_serialized_go_ids(row_b["specific_go_ids"]) if row_b else set()

        unresolved_b = (
            parse_serialized_go_ids(row_b["ontology_unresolved_go_ids"]) if row_b else set()
        )

        all_a = specific_a | unresolved_a

        all_b = specific_b | unresolved_b

        if all_a and not all_b:
            matches = []
            relation = "only_source_a"

        elif all_b and not all_a:
            matches = []
            relation = "only_source_b"

        elif all_a and all_b:
            matches = maximum_go_matching(
                specific_a,
                unresolved_a,
                specific_b,
                unresolved_b,
                ancestor_map,
            )

            relation = classify_go_sets(
                specific_a,
                unresolved_a,
                specific_b,
                unresolved_b,
                ancestor_map,
            )

        else:
            # A comparison-view unit should always
            # contain at least one GO identifier.
            raise ValueError(f"Empty GO comparison unit: {key}")

        matched_a = {go_a for go_a, _, _ in matches}

        matched_b = {go_b for _, go_b, _ in matches}

        comparison = {
            "sample_id": sample_id,
            "feature_uid": feature_uid,
            "protein_id": protein_id,
            "aspect": aspect,
            "source_a": source_a,
            "source_b": source_b,
            "source_a_go_ids": (serialize_go_ids(specific_a)),
            "source_b_go_ids": (serialize_go_ids(specific_b)),
            "source_a_unresolved_go_ids": (serialize_go_ids(unresolved_a)),
            "source_b_unresolved_go_ids": (serialize_go_ids(unresolved_b)),
            "matched_go_pairs": (serialize_go_matches(matches)),
            "source_a_unmatched_go_ids": (serialize_go_ids(all_a - matched_a)),
            "source_b_unmatched_go_ids": (serialize_go_ids(all_b - matched_b)),
            "relation": relation,
        }

        if list(comparison) != GO_SOURCE_COMPARISON_COLUMNS:
            raise ValueError("GO source comparison does not match the canonical schema")

        comparisons.append(comparison)

    return comparisons
