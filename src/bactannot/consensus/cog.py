"""Evidence-aware consensus for COG annotations"""

from collections import defaultdict

from bactannot.schema import (
    COG_CLAIM_COLUMNS,
    COG_EVIDENCE_CLAIM_COLUMNS,
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


def build_cog_claims_from_annotations(
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build exact canonical COG claims"""

    observations_by_claim: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in annotation_rows:
        if row["annotation_type"] != "COG":
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            row["annotation_id"],
        )

        observations_by_claim[key].append(row)

    claims: list[dict[str, str]] = []

    for key in sorted(
        observations_by_claim,
        key=lambda item: tuple(value.casefold() for value in item),
    ):
        (
            sample_id,
            feature_uid,
            claim_id,
        ) = key

        observations = observations_by_claim[key]

        protein_ids = {row["protein_id"] for row in observations}

        if len(protein_ids) != 1:
            raise ValueError(f"Conflicting protein_id values for COG claim: {key}")

        sources = {row["source"] for row in observations}

        claim = {
            "sample_id": sample_id,
            "feature_uid": feature_uid,
            "protein_id": next(iter(protein_ids)),
            "claim_type": "COG",
            "claim_id": claim_id,
            "observing_sources": (serialize_sources(sources)),
            "observation_count": str(len(observations)),
        }

        if list(claim) != COG_CLAIM_COLUMNS:
            raise ValueError("COG claim does not match the expected schema")

        claims.append(claim)

    return claims


def build_cog_evidence_indexes(
    evidence_rows: list[dict[str, str]],
) -> tuple[
    dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ],
    dict[
        tuple[str, str, str, str],
        list[dict[str, str]],
    ],
]:
    """Index eggNOG and PathoFact2 evidence relevant to COG claims."""

    eggnog_by_feature: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    pathofact_by_claim: dict[
        tuple[str, str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in evidence_rows:
        if row["source"] == "eggNOG-mapper" and row["evidence_type"] == "seed_ortholog":
            key = (
                row["sample_id"],
                row["feature_uid"],
                row["protein_id"],
            )

            eggnog_by_feature[key].append(row)

            continue

        if (
            row["source"]
            not in {
                "PathoFact2:VF-HMM",
                "PathoFact2:TOX-HMM",
            }
            or row["evidence_type"] != "hmm_match"
        ):
            continue

        accession = row["accession"].strip().upper()

        if not (accession.startswith("COG") and accession[3:].isdigit()):
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            "COG",
            f"COG:{accession}",
        )

        pathofact_by_claim[key].append(row)

    return (
        dict(eggnog_by_feature),
        dict(pathofact_by_claim),
    )


def collect_cog_claim_evidence(
    claim: dict[str, str],
    eggnog_index: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ],
    pathofact_index: dict[
        tuple[str, str, str, str],
        list[dict[str, str]],
    ],
) -> list[dict[str, str]]:
    """Collect producer-specific evidence supporting one COG claim."""

    if claim["claim_type"] != "COG":
        raise ValueError(f"Expected a COG claim, got {claim['claim_type']!r}")

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    matches: list[dict[str, str]] = []

    if "eggNOG-mapper" in observing_sources:
        eggnog_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
        )

        matches.extend(
            eggnog_index.get(
                eggnog_key,
                [],
            )
        )

    if "PathoFact2" in observing_sources:
        pathofact_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["claim_type"],
            claim["claim_id"],
        )

        matches.extend(
            pathofact_index.get(
                pathofact_key,
                [],
            )
        )

    return sorted(
        matches,
        key=lambda row: (
            row["source"].casefold(),
            row["evidence_type"],
            row["accession"],
        ),
    )


def get_cog_evidence_line(
    evidence_row: dict[str, str],
) -> tuple[str, str] | None:
    """Resolve the operational evidence line for COG support."""

    if (
        evidence_row["source"] == "eggNOG-mapper"
        and evidence_row["evidence_type"] == "seed_ortholog"
    ):
        return (
            "orthology_transfer",
            "eggNOG",
        )

    if (
        evidence_row["source"]
        in {
            "PathoFact2:VF-HMM",
            "PathoFact2:TOX-HMM",
        }
        and evidence_row["evidence_type"] == "hmm_match"
    ):
        accession = evidence_row["accession"].strip().upper()

        if accession.startswith("COG") and accession[3:].isdigit():
            return (
                "profile_hmm",
                "COG",
            )

    return None


def get_cog_evidence_producer(
    evidence_row: dict[str, str],
) -> str | None:
    """Resolve the annotation producer represented by COG evidence."""

    if (
        evidence_row["source"] == "eggNOG-mapper"
        and evidence_row["evidence_type"] == "seed_ortholog"
    ):
        return "eggNOG-mapper"

    if (
        evidence_row["source"]
        in {
            "PathoFact2:VF-HMM",
            "PathoFact2:TOX-HMM",
        }
        and evidence_row["evidence_type"] == "hmm_match"
    ):
        return "PathoFact2"

    return None


def summarize_cog_claim_evidence(
    claim: dict[str, str],
    eggnog_index: dict[
        tuple[str, str, str],
        list[dict[str, str]],
    ],
    pathofact_index: dict[
        tuple[str, str, str, str],
        list[dict[str, str]],
    ],
) -> dict[str, str]:
    """Summarize operational evidence supporting one COG claim."""

    matched_evidence = collect_cog_claim_evidence(
        claim,
        eggnog_index,
        pathofact_index,
    )

    evidence_lines: set[tuple[str, str]] = set()

    evidence_sources: set[str] = set()

    for row in matched_evidence:
        line = get_cog_evidence_line(row)

        producer = get_cog_evidence_producer(row)

        if line is None or producer is None:
            continue

        evidence_lines.add(line)

        evidence_sources.add(producer)

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    context_sources: set[str] = set()

    unresolved_sources = observing_sources - evidence_sources - context_sources

    return {
        "evidence_sources": (serialize_sources(evidence_sources)),
        "evidence_lines": (serialize_evidence_lines(evidence_lines)),
        "independent_evidence_count": str(len(evidence_lines)),
        "context_sources": (serialize_sources(context_sources)),
        "unresolved_sources": (serialize_sources(unresolved_sources)),
    }


def classify_cog_evidence_status(
    evidence: dict[str, str],
) -> str:
    """Classify the evidence status of one COG claim."""

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


def build_cog_claims_with_evidence(
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build COG claims enriched with operational evidence."""

    claims = build_cog_claims_from_annotations(annotation_rows)

    (
        eggnog_index,
        pathofact_index,
    ) = build_cog_evidence_indexes(evidence_rows)

    output: list[dict[str, str]] = []

    for claim in claims:
        evidence = summarize_cog_claim_evidence(
            claim,
            eggnog_index,
            pathofact_index,
        )

        row = {
            **claim,
            **evidence,
            "evidence_status": (classify_cog_evidence_status(evidence)),
        }

        if list(row) != COG_EVIDENCE_CLAIM_COLUMNS:
            raise ValueError("COG evidence claim does not match the expected schema")

        output.append(row)

    return output
