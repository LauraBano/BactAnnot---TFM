"""Evidence-aware consensus for biomedical annotations"""

from collections import defaultdict

from bactannot.schema import (
    ARO_CLAIM_COLUMNS,
    ARO_EVIDENCE_CLAIM_COLUMNS,
    BIOMEDICAL_CLASSIFICATION_CLAIM_COLUMNS,
    BIOMEDICAL_CLASSIFICATION_EVIDENCE_CLAIM_COLUMNS,
)

BIOMEDICAL_CLASSIFICATION_TYPES = frozenset(
    {
        "VIRULENCE_ASSOCIATED",
        "TOXIN_ASSOCIATED",
    }
)

BIOMEDICAL_CLASSIFICATION_IDS = {
    "VIRULENCE_ASSOCIATED": "PATHOFACT2:VF",
    "TOXIN_ASSOCIATED": "PATHOFACT2:TOX",
}

BIOMEDICAL_CLASSIFICATION_STATUSES = frozenset(
    {
        "high_confidence",
        "ambiguous",
    }
)


def serialize_sources(
    sources: set[str],
) -> str:
    """Serialize distinct sources deterministically"""

    return ";".join(
        sorted(
            sources,
            key=str.casefold,
        )
    )


def serialize_evidence_lines(
    evidence_lines: set[tuple[str, str]],
) -> str:
    """Serialize operational evidence lines deterministically"""

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


def build_aro_claims_from_annotations(
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build exact canonical ARO claims."""

    observations_by_claim: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in annotation_rows:
        if row["annotation_type"] != "ARO":
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            row["annotation_id"],
        )

        observations_by_claim[key].append(row)

    claims: list[dict[str, str]] = []

    for key in sorted(observations_by_claim):
        (
            sample_id,
            feature_uid,
            claim_id,
        ) = key

        observations = observations_by_claim[key]

        protein_ids = {row["protein_id"] for row in observations}

        if len(protein_ids) != 1:
            raise ValueError(f"Conflicting protein_id values for ARO claim: {key}")

        sources = {row["source"] for row in observations}

        claim = {
            "sample_id": sample_id,
            "feature_uid": feature_uid,
            "protein_id": next(iter(protein_ids)),
            "claim_type": "ARO",
            "claim_id": claim_id,
            "observing_sources": (serialize_sources(sources)),
            "observation_count": str(len(observations)),
        }

        if list(claim) != ARO_CLAIM_COLUMNS:
            raise ValueError("ARO claim does not match the expected schema")

        claims.append(claim)

    return claims


def get_aro_evidence_line(
    evidence_row: dict[str, str],
) -> tuple[str, str] | None:
    """Resolve the operational evidence line for an ARO claim"""

    if (
        evidence_row["source"] == "PathoFact2:DeepARG"
        and evidence_row["evidence_type"] == "amr_prediction"
    ):
        return (
            "deep_learning",
            "DeepARG",
        )

    if (
        evidence_row["source"] == "PathoFact2:RGI"
        and evidence_row["evidence_type"] == "amr_prediction"
    ):
        return (
            "card_detection_model",
            "CARD",
        )

    return None


def build_aro_evidence_index(
    evidence_rows: list[dict[str, str]],
) -> dict[
    tuple[str, str, str],
    list[dict[str, str]],
]:
    """Index native AMR evidence by exact canonical ARO claim"""

    index: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in evidence_rows:
        if (
            row["source"]
            not in {
                "PathoFact2:DeepARG",
                "PathoFact2:RGI",
            }
            or row["evidence_type"] != "amr_prediction"
        ):
            continue

        accession = row["accession"].strip()

        if not accession.startswith("ARO:"):
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            accession,
        )

        index[key].append(row)

    return dict(index)


def collect_aro_claim_evidence(
    claim: dict[str, str],
    evidence_index: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ],
) -> list[dict[str, str]]:
    """Collect native evidence supporting one exact ARO claim"""

    if claim["claim_type"] != "ARO":
        raise ValueError(f"Expected an ARO claim, got {claim['claim_type']!r}")

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    if "PathoFact2" not in observing_sources:
        return []

    key = (
        claim["sample_id"],
        claim["feature_uid"],
        claim["claim_id"],
    )

    matches = evidence_index.get(
        key,
        [],
    )

    return sorted(
        matches,
        key=lambda row: (
            row["source"].casefold(),
            row["accession"],
        ),
    )


def get_aro_claim_relation(
    source: str,
    evidence_type: str,
) -> str | None:
    """Describe how a supported AMR channel assigns an ARO identifier."""

    if evidence_type != "amr_prediction":
        return None

    if source == "PathoFact2:RGI":
        return "native_aro_assignment"

    if source == "PathoFact2:DeepARG":
        return "mapped_from_native_hit"

    return None


def summarize_aro_claim_evidence(
    claim: dict[str, str],
    evidence_index: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ],
) -> dict[str, str]:
    """Summarize operational evidence supporting one ARO claim."""

    matched_evidence = collect_aro_claim_evidence(
        claim,
        evidence_index,
    )

    evidence_lines: set[tuple[str, str]] = set()
    evidence_relations: set[str] = set()

    for row in matched_evidence:
        line = get_aro_evidence_line(row)

        if line is not None:
            evidence_lines.add(line)

            relation = get_aro_claim_relation(
                row["source"],
                row["evidence_type"],
            )

            if relation is not None:
                evidence_relations.add(relation)

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    resolved_sources: set[str] = set()

    if evidence_lines and "PathoFact2" in observing_sources:
        resolved_sources.add("PathoFact2")

    context_sources: set[str] = set()

    unresolved_sources = observing_sources - resolved_sources - context_sources

    return {
        "evidence_sources": serialize_sources(resolved_sources),
        "evidence_lines": serialize_evidence_lines(evidence_lines),
        "evidence_relations": ";".join(sorted(evidence_relations)),
        "independent_evidence_count": str(len(evidence_lines)),
        "context_sources": serialize_sources(context_sources),
        "unresolved_sources": serialize_sources(unresolved_sources),
    }


def classify_aro_evidence_status(
    evidence: dict[str, str],
) -> str:
    """Classify the evidence status of one ARO claim"""

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


def build_aro_claims_with_evidence(
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build ARO claims enriched with operational evidence"""

    claims = build_aro_claims_from_annotations(annotation_rows)

    evidence_index = build_aro_evidence_index(evidence_rows)

    output: list[dict[str, str]] = []

    for claim in claims:
        evidence = summarize_aro_claim_evidence(
            claim,
            evidence_index,
        )

        row = {
            **claim,
            **evidence,
            "evidence_status": (classify_aro_evidence_status(evidence)),
        }

        if list(row) != ARO_EVIDENCE_CLAIM_COLUMNS:
            raise ValueError("ARO evidence claim does not match the expected schema")

        output.append(row)

    return output


def build_biomedical_classification_claims_from_annotations(
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build canonical virulence and toxin classification claims."""

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
        annotation_type = row["annotation_type"]

        if annotation_type not in BIOMEDICAL_CLASSIFICATION_TYPES:
            continue

        expected_id = BIOMEDICAL_CLASSIFICATION_IDS[annotation_type]

        if row["annotation_id"] != expected_id:
            raise ValueError(
                "Unexpected biomedical "
                "classification identifier: "
                f"{row['annotation_id']!r} "
                f"for {annotation_type!r}"
            )

        key = (
            row["sample_id"],
            row["feature_uid"],
            annotation_type,
            row["annotation_id"],
        )

        observations_by_claim[key].append(row)

    claims: list[dict[str, str]] = []

    for key in sorted(observations_by_claim):
        (
            sample_id,
            feature_uid,
            claim_type,
            claim_id,
        ) = key

        observations = observations_by_claim[key]

        protein_ids = {row["protein_id"] for row in observations}

        if len(protein_ids) != 1:
            raise ValueError(
                f"Conflicting protein_id values for biomedical classification claim: {key}"
            )

        statuses = {row["annotation_value"] for row in observations if row["annotation_value"]}

        if len(statuses) != 1:
            raise ValueError(
                f"Expected exactly one classification status for biomedical claim: {key}"
            )

        classification_status = next(iter(statuses))

        if classification_status not in BIOMEDICAL_CLASSIFICATION_STATUSES:
            raise ValueError(
                f"Unexpected biomedical classification status: {classification_status!r}"
            )

        sources = {row["source"] for row in observations}

        claim = {
            "sample_id": sample_id,
            "feature_uid": feature_uid,
            "protein_id": next(iter(protein_ids)),
            "claim_type": claim_type,
            "claim_id": claim_id,
            "classification_status": (classification_status),
            "observing_sources": (serialize_sources(sources)),
            "observation_count": str(len(observations)),
        }

        if list(claim) != BIOMEDICAL_CLASSIFICATION_CLAIM_COLUMNS:
            raise ValueError("Biomedical classification claim does not match the expected schema")

        claims.append(claim)

    return claims


def get_biomedical_classification_evidence_line(
    evidence_row: dict[str, str],
) -> tuple[str, str, str] | None:
    """
    Resolve biomedical branch and operational evidence line.

    Returns:
        claim_type, method_group, resource_group
    """

    source = evidence_row["source"]
    evidence_type = evidence_row["evidence_type"]

    if source == "PathoFact2:VF-HMM" and evidence_type == "hmm_match":
        return (
            "VIRULENCE_ASSOCIATED",
            "profile_hmm",
            "PathoFact2_VF",
        )

    if source == "PathoFact2:VF-ML" and evidence_type == "ml_prediction":
        return (
            "VIRULENCE_ASSOCIATED",
            "machine_learning",
            "PathoFact2_VF",
        )

    if source == "PathoFact2:TOX-HMM" and evidence_type == "hmm_match":
        return (
            "TOXIN_ASSOCIATED",
            "profile_hmm",
            "PathoFact2_TOX",
        )

    if source == "PathoFact2:TOX-ML" and evidence_type == "ml_prediction":
        return (
            "TOXIN_ASSOCIATED",
            "machine_learning",
            "PathoFact2_TOX",
        )

    return None


def build_biomedical_classification_evidence_index(
    evidence_rows: list[dict[str, str]],
) -> dict[
    tuple[str, str, str],
    list[dict[str, str]],
]:
    """Index native evidence by biomedical branch and feature."""

    index: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in evidence_rows:
        resolved = get_biomedical_classification_evidence_line(row)

        if resolved is None:
            continue

        (
            claim_type,
            _method_group,
            _resource_group,
        ) = resolved

        key = (
            row["sample_id"],
            row["feature_uid"],
            claim_type,
        )

        index[key].append(row)

    return dict(index)


def collect_biomedical_classification_evidence(
    claim: dict[str, str],
    evidence_index: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ],
) -> list[dict[str, str]]:
    """Collect native evidence supporting one biomedical classification."""

    claim_type = claim["claim_type"]

    if claim_type not in BIOMEDICAL_CLASSIFICATION_TYPES:
        raise ValueError(f"Expected a biomedical classification claim, got {claim_type!r}")

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    if "PathoFact2" not in observing_sources:
        return []

    key = (
        claim["sample_id"],
        claim["feature_uid"],
        claim_type,
    )

    return sorted(
        evidence_index.get(
            key,
            [],
        ),
        key=lambda row: (
            row["source"].casefold(),
            row["evidence_type"],
            row["accession"],
        ),
    )


def summarize_biomedical_classification_evidence(
    claim: dict[str, str],
    evidence_index: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ],
) -> dict[str, str]:
    """Summarize evidence supporting one biomedical classification claim."""

    matched_evidence = collect_biomedical_classification_evidence(
        claim,
        evidence_index,
    )

    evidence_lines: set[tuple[str, str]] = set()

    for row in matched_evidence:
        resolved = get_biomedical_classification_evidence_line(row)

        if resolved is None:
            continue

        (
            evidence_claim_type,
            method_group,
            resource_group,
        ) = resolved

        if evidence_claim_type != claim["claim_type"]:
            continue

        evidence_lines.add(
            (
                method_group,
                resource_group,
            )
        )

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    resolved_sources: set[str] = set()

    if evidence_lines and "PathoFact2" in observing_sources:
        resolved_sources.add("PathoFact2")

    context_sources: set[str] = set()

    unresolved_sources = observing_sources - resolved_sources - context_sources

    return {
        "evidence_sources": (serialize_sources(resolved_sources)),
        "evidence_lines": (serialize_evidence_lines(evidence_lines)),
        "independent_evidence_count": str(len(evidence_lines)),
        "context_sources": (serialize_sources(context_sources)),
        "unresolved_sources": (serialize_sources(unresolved_sources)),
    }


def classify_biomedical_classification_evidence_status(
    evidence: dict[str, str],
) -> str:
    """Classify evidence status for one biomedical classification."""

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


def build_biomedical_classification_claims_with_evidence(
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build biomedical classification claims with operational evidence."""

    claims = build_biomedical_classification_claims_from_annotations(annotation_rows)

    evidence_index = build_biomedical_classification_evidence_index(evidence_rows)

    output: list[dict[str, str]] = []

    for claim in claims:
        evidence = summarize_biomedical_classification_evidence(
            claim,
            evidence_index,
        )

        row = {
            **claim,
            **evidence,
            "evidence_status": (classify_biomedical_classification_evidence_status(evidence)),
        }

        if list(row) != BIOMEDICAL_CLASSIFICATION_EVIDENCE_CLAIM_COLUMNS:
            raise ValueError(
                "Biomedical classification evidence claim does not match the expected schema"
            )

        output.append(row)

    return output
