"""Canonical EC claims for evidence-aware consensus"""

from collections import defaultdict

from bactannot.consensus.evidence_registry import (
    get_evidence_channel,
)

EC_CLAIM_RELATIONS = {
    (
        "eggNOG-mapper",
        "seed_ortholog",
    ): "transferred_annotation",
    (
        "KofamScan",
        "hmm_match",
    ): "derived_from_ko",
    (
        "BlastRules",
        "expert",
    ): "transferred_annotation",
}


EC_ORIGIN_LINKED_CHANNELS = {
    (
        "KofamScan",
        "hmm_match",
    ),
    (
        "BlastRules",
        "expert",
    ),
}


def get_ec_claim_relation(
    source: str,
    evidence_type: str,
) -> str | None:
    """Return the EC-specific relation for an evidence channel"""

    return EC_CLAIM_RELATIONS.get(
        (
            source,
            evidence_type,
        )
    )


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


def serialize_evidence_lines(
    evidence_lines: set[tuple[str, str]],
) -> str:
    """Serialize distinct method-resource evidence lines deterministically."""

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
    """Serialize distinct EC evidence relations deterministically"""

    return ";".join(
        sorted(
            relations,
            key=str.casefold,
        )
    )


def serialize_ec_ids(
    ec_ids: set[str],
) -> str:
    """Serialize EC identifiers deterministically."""

    return ";".join(sorted(ec_ids))


def serialize_ec_matching(
    matching: list[tuple[str, str]],
) -> str:
    """Serialize one deterministic maximum EC matching."""

    return ";".join(f"{ec_a}~{ec_b}" for ec_a, ec_b in sorted(matching))


def build_ec_claims_from_annotations(
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build canonical EC claims from annotation observations"""

    observations_by_claim: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in annotation_rows:
        if row["annotation_type"] != "EC":
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
            raise ValueError(f"Conflicting protein_id values for EC claim: {key}")

        protein_id = next(iter(protein_ids))

        sources = {row["source"] for row in observations}

        claims.append(
            {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "claim_type": "EC",
                "claim_id": claim_id,
                "observing_sources": serialize_sources(sources),
                "observation_count": str(len(observations)),
            }
        )

    claims.sort(
        key=lambda row: (
            row["feature_uid"],
            row["claim_id"],
        )
    )

    return claims


def collect_ec_claim_observations(
    claim: dict[str, str],
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Collect annotation observations belonging to one EC claim"""

    if claim["claim_type"] != "EC":
        raise ValueError(f"Expected an EC claim, got {claim['claim_type']!r}")

    matches = [
        row
        for row in annotation_rows
        if (
            row["sample_id"] == claim["sample_id"]
            and row["feature_uid"] == claim["feature_uid"]
            and row["protein_id"] == claim["protein_id"]
            and row["annotation_type"] == "EC"
            and row["annotation_id"] == claim["claim_id"]
        )
    ]

    matches.sort(
        key=lambda row: (
            row["source"],
            row["annotation_value"],
        )
    )

    return matches


def parse_ec_annotation_origins(
    annotation_value: str,
) -> set[str]:
    """Parse explicit origin accessions stored for one EC observation"""

    return {value.strip() for value in annotation_value.split(",") if value.strip()}


def collect_ec_evidence(
    claim: dict[str, str],
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Collect methodological evidence supporting one EC claim"""

    observations = collect_ec_claim_observations(
        claim,
        annotation_rows,
    )

    matches: list[dict[str, str]] = []

    for row in evidence_rows:
        channel = get_evidence_channel(
            row["source"],
            row["evidence_type"],
        )

        if channel is None:
            continue

        relation = get_ec_claim_relation(
            row["source"],
            row["evidence_type"],
        )

        if relation is None:
            continue

        if not (
            row["sample_id"] == claim["sample_id"]
            and row["feature_uid"] == claim["feature_uid"]
            and row["protein_id"] == claim["protein_id"]
        ):
            continue

        producer_observations = [
            observation for observation in observations if observation["source"] == channel.producer
        ]

        if not producer_observations:
            continue

        channel_key = (
            row["source"],
            row["evidence_type"],
        )

        if channel_key in EC_ORIGIN_LINKED_CHANNELS:
            origins: set[str] = set()

            for observation in producer_observations:
                origins.update(parse_ec_annotation_origins(observation["annotation_value"]))

            if row["accession"] not in origins:
                continue

        matches.append(row)

    matches.sort(
        key=lambda row: (
            row["source"],
            row["evidence_type"],
            row["accession"],
        )
    )

    return matches


def summarize_ec_claim_evidence(
    claim: dict[str, str],
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> dict[str, str]:
    """Summarize methodological evidence supporting one EC claim"""

    matched_evidence = collect_ec_evidence(
        claim,
        annotation_rows,
        evidence_rows,
    )

    evidence_sources: set[str] = set()

    evidence_lines: set[tuple[str, str]] = set()

    relations: set[str] = set()

    resolved_producers: set[str] = set()

    for row in matched_evidence:
        channel = get_evidence_channel(
            row["source"],
            row["evidence_type"],
        )

        if channel is None:
            continue

        relation = get_ec_claim_relation(
            row["source"],
            row["evidence_type"],
        )

        if relation is None:
            continue

        evidence_sources.add(channel.source)

        evidence_lines.add(
            (
                channel.method_group,
                channel.resource_group,
            )
        )

        relations.add(relation)

        resolved_producers.add(channel.producer)

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    unresolved_sources = observing_sources - resolved_producers

    return {
        "evidence_sources": serialize_sources(evidence_sources),
        "evidence_lines": serialize_evidence_lines(evidence_lines),
        "evidence_relations": serialize_relations(relations),
        "independent_evidence_count": str(len(evidence_lines)),
        "unresolved_sources": serialize_sources(unresolved_sources),
    }


def classify_ec_evidence_status(
    evidence: dict[str, str],
) -> str:
    """Classify methodological evidence support for one EC claim"""

    independent_evidence_count = int(evidence["independent_evidence_count"])

    if independent_evidence_count >= 2:
        return "multiple_independent"

    if independent_evidence_count == 1:
        return "single_independent"

    if evidence["unresolved_sources"]:
        return "unresolved"

    return "no_evidence"


def build_ec_claims_with_evidence(
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build EC claims enriched with methodological evidence."""

    claims = build_ec_claims_from_annotations(annotation_rows)

    annotations_by_feature: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in annotation_rows:
        key = (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
        )

        annotations_by_feature[key].append(row)

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
        feature_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
        )

        candidate_annotations = annotations_by_feature.get(
            feature_key,
            [],
        )

        candidate_evidence = evidence_by_feature.get(
            feature_key,
            [],
        )

        evidence = summarize_ec_claim_evidence(
            claim,
            candidate_annotations,
            candidate_evidence,
        )

        evidence_claim = {
            **claim,
            **evidence,
            "evidence_status": (classify_ec_evidence_status(evidence)),
        }

        evidence_claims.append(evidence_claim)

    return evidence_claims


def parse_ec_components(
    ec_id: str,
) -> tuple[str, str, str, str]:
    """Parse one canonical EC identifier into its four components."""

    if not ec_id.startswith("EC:"):
        raise ValueError(f"Expected canonical EC identifier, got {ec_id!r}")

    components = ec_id.removeprefix("EC:").split(".")

    if len(components) != 4:
        raise ValueError(f"Expected four EC components, got {ec_id!r}")

    return (
        components[0],
        components[1],
        components[2],
        components[3],
    )


def classify_ec_relation(
    ec_a: str,
    ec_b: str,
) -> str:
    """Classify the hierarchical relation between two EC identifiers"""

    components_a = parse_ec_components(ec_a)

    components_b = parse_ec_components(ec_b)

    if components_a == components_b:
        return "exact"

    compatible = all(
        (component_a == component_b or component_a == "-" or component_b == "-")
        for component_a, component_b in zip(
            components_a,
            components_b,
            strict=True,
        )
    )

    if compatible:
        return "hierarchical_compatible"

    return "disjoint"


def maximum_ec_matching(
    ecs_a: set[str],
    ecs_b: set[str],
) -> list[tuple[str, str]]:
    """Find a maximum one-to-one matching between compatible EC claims."""

    adjacency: dict[
        str,
        list[str],
    ] = {}

    for ec_a in sorted(ecs_a):
        compatible: list[tuple[str, str]] = []

        for ec_b in sorted(ecs_b):
            relation = classify_ec_relation(
                ec_a,
                ec_b,
            )

            if relation == "disjoint":
                continue

            compatible.append(
                (
                    ec_b,
                    relation,
                )
            )

        compatible.sort(
            key=lambda item: (
                0 if item[1] == "exact" else 1,
                item[0],
            )
        )

        adjacency[ec_a] = [ec_b for ec_b, _ in compatible]

    matched_b: dict[
        str,
        str,
    ] = {}

    def augment(
        ec_a: str,
        visited_b: set[str],
    ) -> bool:
        for ec_b in adjacency[ec_a]:
            if ec_b in visited_b:
                continue

            visited_b.add(ec_b)

            if ec_b not in matched_b or augment(
                matched_b[ec_b],
                visited_b,
            ):
                matched_b[ec_b] = ec_a

                return True

        return False

    for ec_a in sorted(ecs_a):
        augment(
            ec_a,
            set(),
        )

    pairs = [
        (
            ec_a,
            ec_b,
        )
        for ec_b, ec_a in matched_b.items()
    ]

    pairs.sort()

    return pairs


def classify_ec_sets(
    ecs_a: set[str],
    ecs_b: set[str],
) -> tuple[
    str,
    list[tuple[str, str]],
]:
    """Classify the relation between two non-empty EC claim sets."""

    if not ecs_a or not ecs_b:
        raise ValueError("EC set comparison requires two non-empty sets")

    if ecs_a == ecs_b:
        return (
            "exact",
            [
                (
                    ec_id,
                    ec_id,
                )
                for ec_id in sorted(ecs_a)
            ],
        )

    matching = maximum_ec_matching(
        ecs_a,
        ecs_b,
    )

    matched_count = len(matching)

    if matched_count == 0:
        return (
            "disjoint",
            matching,
        )

    if matched_count == len(ecs_a) and matched_count == len(ecs_b):
        return (
            "hierarchical_compatible",
            matching,
        )

    return (
        "partial",
        matching,
    )


def build_ec_source_comparisons(
    annotation_rows: list[dict[str, str]],
    source_a: str,
    source_b: str,
) -> list[dict[str, str]]:
    """Build per-feature EC comparisons between two annotation sources."""

    if source_a == source_b:
        raise ValueError("EC source comparison requires two different sources")

    observations_by_feature: dict[
        tuple[str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in annotation_rows:
        if row["annotation_type"] != "EC":
            continue

        if row["source"] not in {
            source_a,
            source_b,
        }:
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
        )

        observations_by_feature[key].append(row)

    comparisons: list[dict[str, str]] = []

    for key, observations in observations_by_feature.items():
        (
            sample_id,
            feature_uid,
        ) = key

        protein_ids = {row["protein_id"] for row in observations}

        if len(protein_ids) != 1:
            raise ValueError(f"Conflicting protein_id values for EC feature: {key}")

        protein_id = next(iter(protein_ids))

        ecs_a = {row["annotation_id"] for row in observations if row["source"] == source_a}

        ecs_b = {row["annotation_id"] for row in observations if row["source"] == source_b}

        matching: list[tuple[str, str]] = []

        if ecs_a and not ecs_b:
            relation = "only_source_a"

        elif ecs_b and not ecs_a:
            relation = "only_source_b"

        else:
            (
                relation,
                matching,
            ) = classify_ec_sets(
                ecs_a,
                ecs_b,
            )

        matched_a = {ec_a for ec_a, _ in matching}

        matched_b = {ec_b for _, ec_b in matching}

        unmatched_a = ecs_a - matched_a

        unmatched_b = ecs_b - matched_b

        comparisons.append(
            {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "source_a": source_a,
                "source_b": source_b,
                "source_a_ecs": serialize_ec_ids(ecs_a),
                "source_b_ecs": serialize_ec_ids(ecs_b),
                "matched_ec_pairs": (serialize_ec_matching(matching)),
                "source_a_unmatched_ecs": (serialize_ec_ids(unmatched_a)),
                "source_b_unmatched_ecs": (serialize_ec_ids(unmatched_b)),
                "relation": relation,
            }
        )

    comparisons.sort(
        key=lambda row: (
            row["sample_id"],
            row["feature_uid"],
        )
    )

    return comparisons
