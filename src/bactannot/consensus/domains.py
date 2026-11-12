"""Evidence-aware consensus for protein domain and family annotations"""

import re
from collections import defaultdict

from bactannot.schema import (
    DOMAIN_EVIDENCE_CLAIM_COLUMNS,
)

DOMAIN_CLAIM_TYPES = frozenset(
    {
        "CDD",
        "NCBIFAM",
        "PFAM",
        "SMART",
    }
)


PFAM_NATIVE_RE = re.compile(
    r"^(?:PFAM)?(PF\d+)(?:\.\d+)?$",
    re.IGNORECASE,
)

CDD_NATIVE_RE = re.compile(
    r"^(?:CDD:)?(cd\d+)$",
    re.IGNORECASE,
)

SMART_NATIVE_RE = re.compile(
    r"^(?:SMART)?(SM\d+)$",
    re.IGNORECASE,
)

NCBIFAM_NATIVE_RE = re.compile(
    r"^(?:NCBIFAM:)?((?:NF\d+|TIGR\d+))(?:\.\d+)?$",
    re.IGNORECASE,
)


def canonicalize_domain_evidence(
    evidence_row: dict[str, str],
) -> tuple[str, str] | None:
    """Map native domain evidence to a canonical domain claim"""

    source = evidence_row["source"]

    evidence_type = evidence_row["evidence_type"]

    accession = evidence_row["accession"].strip()

    if not accession:
        return None

    # InterProScan member-database signatures
    interpro_map = {
        (
            "InterProScan:CDD",
            "signature_match",
        ): (
            "CDD",
            CDD_NATIVE_RE,
            "CDD",
        ),
        (
            "InterProScan:NCBIfam",
            "signature_match",
        ): (
            "NCBIFAM",
            NCBIFAM_NATIVE_RE,
            "NCBIFAM",
        ),
        (
            "InterProScan:Pfam",
            "signature_match",
        ): (
            "PFAM",
            PFAM_NATIVE_RE,
            "PFAM",
        ),
        (
            "InterProScan:SMART",
            "signature_match",
        ): (
            "SMART",
            SMART_NATIVE_RE,
            "SMART",
        ),
    }

    mapping = interpro_map.get(
        (
            source,
            evidence_type,
        )
    )

    if mapping is not None:
        (
            claim_type,
            pattern,
            prefix,
        ) = mapping

        match = pattern.fullmatch(accession)

        if match is None:
            return None

        native_id = match.group(1)

        if claim_type == "CDD":
            native_id = "cd" + native_id[2:]

        else:
            native_id = native_id.upper()

        return (
            claim_type,
            f"{prefix}:{native_id}",
        )

    # PathoFact2 HMM profiles
    if (
        source
        in {
            "PathoFact2:VF-HMM",
            "PathoFact2:TOX-HMM",
        }
        and evidence_type == "hmm_match"
    ):
        match = re.fullmatch(
            r"pfam(\d+)",
            accession,
            re.IGNORECASE,
        )

        if match:
            return (
                "PFAM",
                "PFAM:PF" + match.group(1),
            )

        match = re.fullmatch(
            r"cd(\d+)",
            accession,
            re.IGNORECASE,
        )

        if match:
            return (
                "CDD",
                "CDD:cd" + match.group(1),
            )

        match = re.fullmatch(
            r"smart(\d+)",
            accession,
            re.IGNORECASE,
        )

        if match:
            return (
                "SMART",
                "SMART:SM" + match.group(1),
            )

        match = NCBIFAM_NATIVE_RE.fullmatch(accession)

        if match:
            return (
                "NCBIFAM",
                "NCBIFAM:" + match.group(1).upper(),
            )

        return None

    # Bakta native Pfam evidence
    if source == "Pfam" and evidence_type == "pfam":
        match = PFAM_NATIVE_RE.fullmatch(accession)

        if match is None:
            return None

        return (
            "PFAM",
            "PFAM:" + match.group(1).upper(),
        )

    # Bakta exposes NCBIfam models through
    # AMRFinderPlus expert evidence
    if source == "AMRFinderPlus" and evidence_type == "expert":
        match = NCBIFAM_NATIVE_RE.fullmatch(accession)

        if match is None:
            return None

        return (
            "NCBIFAM",
            "NCBIFAM:" + match.group(1).upper(),
        )

    return None


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


def build_domain_claims_from_annotations(
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build exact canonical domain and family claims"""

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
        claim_type = row["annotation_type"]

        if claim_type not in DOMAIN_CLAIM_TYPES:
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            claim_type,
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
            claim_type,
            claim_id,
        ) = key

        observations = observations_by_claim[key]

        protein_ids = {row["protein_id"] for row in observations if row["protein_id"]}

        if len(protein_ids) > 1:
            raise ValueError(f"Domain claim has inconsistent protein IDs: {key}")

        protein_id = next(iter(protein_ids)) if protein_ids else ""

        sources = {row["source"] for row in observations}

        claim = {
            "sample_id": sample_id,
            "feature_uid": feature_uid,
            "protein_id": protein_id,
            "claim_type": claim_type,
            "claim_id": claim_id,
            "observing_sources": (serialize_sources(sources)),
            "observation_count": str(len(observations)),
        }

        claims.append(claim)

    return claims


def build_domain_evidence_index(
    evidence_rows: list[dict[str, str]],
) -> dict[
    tuple[
        str,
        str,
        str,
        str,
    ],
    list[dict[str, str]],
]:
    """Index native evidence by canonical domain claim"""

    index: dict[
        tuple[
            str,
            str,
            str,
            str,
        ],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in evidence_rows:
        canonical = canonicalize_domain_evidence(row)

        if canonical is None:
            continue

        (
            claim_type,
            claim_id,
        ) = canonical

        key = (
            row["sample_id"],
            row["feature_uid"],
            claim_type,
            claim_id,
        )

        index[key].append(row)

    return dict(index)


def domain_evidence_belongs_to_source(
    annotation_source: str,
    evidence_row: dict[str, str],
) -> bool:
    """Check whether native evidence belongs to an annotation producer"""

    evidence_source = evidence_row["source"]

    if annotation_source == "InterProScan":
        return evidence_source.startswith("InterProScan:")

    if annotation_source == "PathoFact2":
        return evidence_source in {
            "PathoFact2:VF-HMM",
            "PathoFact2:TOX-HMM",
        }

    if annotation_source == "bakta":
        return evidence_source in {
            "Pfam",
            "AMRFinderPlus",
        }

    return False


def get_domain_evidence_line(
    evidence_row: dict[str, str],
) -> tuple[str, str] | None:
    """Resolve the operational evidence line for a domain claim."""

    canonical = canonicalize_domain_evidence(evidence_row)

    if canonical is None:
        return None

    claim_type, _ = canonical

    source = evidence_row["source"]
    evidence_type = evidence_row["evidence_type"]

    if source == "InterProScan:CDD" and evidence_type == "signature_match":
        return (
            "rpsblast_pssm",
            "CDD",
        )

    if (
        source
        in {
            "InterProScan:Pfam",
            "InterProScan:SMART",
            "InterProScan:NCBIfam",
        }
        and evidence_type == "signature_match"
    ):
        return (
            "profile_hmm",
            claim_type,
        )

    if (
        source
        in {
            "PathoFact2:VF-HMM",
            "PathoFact2:TOX-HMM",
        }
        and evidence_type == "hmm_match"
    ):
        return (
            "profile_hmm",
            claim_type,
        )

    if source == "Pfam" and evidence_type == "pfam":
        return (
            "profile_hmm",
            "PFAM",
        )

    if source == "AMRFinderPlus" and evidence_type == "expert" and claim_type == "NCBIFAM":
        return (
            "profile_hmm",
            "NCBIFAM",
        )

    return None


def serialize_evidence_lines(
    evidence_lines: set[tuple[str, str]],
) -> str:
    """Serialize operational evidence lines deterministically."""

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


def summarize_domain_claim_evidence(
    claim: dict[str, str],
    evidence_index: dict[
        tuple[
            str,
            str,
            str,
            str,
        ],
        list[dict[str, str]],
    ],
) -> dict[str, str]:
    """Summarize operational evidence supporting one domain claim"""

    candidates = collect_domain_claim_evidence(
        claim,
        evidence_index,
    )

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    evidence_sources: set[str] = set()

    evidence_lines: set[tuple[str, str]] = set()

    for source in observing_sources:
        source_lines = {
            line
            for row in candidates
            if domain_evidence_belongs_to_source(
                source,
                row,
            )
            and (line := get_domain_evidence_line(row)) is not None
        }

        if not source_lines:
            continue

        evidence_sources.add(source)

        evidence_lines.update(source_lines)

    context_sources: set[str] = set()

    unresolved_sources = observing_sources - evidence_sources - context_sources

    return {
        "evidence_sources": (serialize_sources(evidence_sources)),
        "evidence_lines": (serialize_evidence_lines(evidence_lines)),
        "independent_evidence_count": str(len(evidence_lines)),
        "context_sources": (serialize_sources(context_sources)),
        "unresolved_sources": (serialize_sources(unresolved_sources)),
    }


def classify_domain_evidence_status(
    evidence: dict[str, str],
) -> str:
    """Classify the evidence status of one domain claim."""

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


def build_domain_claims_with_evidence(
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build domain claims enriched with operational evidence."""

    claims = build_domain_claims_from_annotations(annotation_rows)

    evidence_index = build_domain_evidence_index(evidence_rows)

    output: list[dict[str, str]] = []

    for claim in claims:
        evidence = summarize_domain_claim_evidence(
            claim,
            evidence_index,
        )

        row = {
            **claim,
            **evidence,
            "evidence_status": (classify_domain_evidence_status(evidence)),
        }

        if list(row) != DOMAIN_EVIDENCE_CLAIM_COLUMNS:
            raise ValueError("Domain evidence claim does not match the expected schema")

        output.append(row)

    return output


def collect_domain_claim_evidence(
    claim: dict[str, str],
    evidence_index: dict[
        tuple[str, str, str, str],
        list[dict[str, str]],
    ],
) -> list[dict[str, str]]:
    """Collect producer-specific native evidence for one domain or family claim."""

    if claim["claim_type"] not in DOMAIN_CLAIM_TYPES:
        raise ValueError(f"Expected a domain or family claim, got {claim['claim_type']!r}")

    key = (
        claim["sample_id"],
        claim["feature_uid"],
        claim["claim_type"],
        claim["claim_id"],
    )

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    matches: list[dict[str, str]] = []

    for row in evidence_index.get(key, []):
        if any(domain_evidence_belongs_to_source(source, row) for source in observing_sources):
            matches.append(row)

    return matches
