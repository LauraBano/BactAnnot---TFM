"""Local reference relations for functional consensus"""

import csv
import re
import warnings
from collections import defaultdict
from pathlib import Path

from bactannot.consensus.ec import classify_ec_relation
from bactannot.consensus.go_ontology import GoTerm, resolve_go_id

KO_RE = re.compile(r"K\d{5}")
EC_BLOCK_RE = re.compile(r"\[EC:([^\]]*)\]")
EC_ID_RE = re.compile(r"\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)")


def load_kofam_ko_ec_map(
    path: Path,
) -> dict[str, set[str]]:
    """Load explicit KO-EC links, reporting unsupported EC tokens"""

    ko_ec_map: dict[str, set[str]] = {}

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"knum", "definition"}.difference(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"Missing KOfam KO-list columns in {path}: {', '.join(sorted(missing))}"
            )

        for row in reader:
            ko = (row["knum"] or "").strip()
            definition = row["definition"]

            if not KO_RE.fullmatch(ko) or definition is None or None in row:
                raise ValueError(f"Malformed KOfam KO-list row in {path} at line {reader.line_num}")

            ko_id = f"KEGG:{ko}"

            if ko_id in ko_ec_map:
                raise ValueError(f"Duplicate KOfam KO identifier in {path}: {ko}")

            ec_ids: set[str] = set()

            for block in EC_BLOCK_RE.findall(definition):
                for token in re.split(r"[\s,]+", block.strip()):
                    if not EC_ID_RE.fullmatch(token):
                        warnings.warn(
                            f"Skipping unsupported EC identifier {token!r} "
                            f"for {ko_id} in {path} "
                            f"at line {reader.line_num}",
                            UserWarning,
                            stacklevel=2,
                        )
                        continue

                    ec_ids.add(f"EC:{token}")

            ko_ec_map[ko_id] = ec_ids

    return ko_ec_map


def classify_ko_ec_relation(
    ko_id: str,
    ec_id: str,
    ko_ec_map: dict[str, set[str]],
) -> str:
    """Classify a KO-EC relation within the supplied reference map"""

    if not ko_id.startswith("KEGG:") or not KO_RE.fullmatch(ko_id.removeprefix("KEGG:")):
        raise ValueError(f"Expected a canonical KO identifier, got {ko_id!r}")

    if not ec_id.startswith("EC:") or not EC_ID_RE.fullmatch(ec_id.removeprefix("EC:")):
        raise ValueError(f"Expected a canonical EC identifier, got {ec_id!r}")

    reference_ecs = ko_ec_map.get(ko_id)

    if not reference_ecs:
        return "unmapped_ko"

    relations = {classify_ec_relation(reference_ec, ec_id) for reference_ec in reference_ecs}

    if "exact" in relations:
        return "documented"

    if "hierarchical_compatible" in relations:
        return "hierarchical_compatible"

    # Missing reference links do not prove biological incompatibility
    return "not_documented"


def build_ko_ec_relations(
    ko_claims: list[dict[str, str]],
    ec_claims: list[dict[str, str]],
    ko_ec_map: dict[str, set[str]],
    reference_id: str,
) -> list[dict[str, str]]:
    """Relate observed KO and EC claims within the same sample and feature"""

    if not reference_id.strip():
        raise ValueError("Missing KO-EC reference identifier")

    ec_claims_by_feature: dict[
        tuple[str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for claim in ec_claims:
        if claim["claim_type"] != "EC":
            raise ValueError(f"Expected an EC claim, got {claim['claim_type']!r}")

        key = (claim["sample_id"], claim["feature_uid"])
        ec_claims_by_feature[key].append(claim)

    relations: list[dict[str, str]] = []

    for ko_claim in ko_claims:
        if ko_claim["claim_type"] != "KEGG_KO":
            raise ValueError(f"Expected a KEGG_KO claim, got {ko_claim['claim_type']!r}")

        key = (ko_claim["sample_id"], ko_claim["feature_uid"])

        for ec_claim in ec_claims_by_feature.get(key, []):
            if ko_claim["protein_id"] != ec_claim["protein_id"]:
                raise ValueError(f"Conflicting protein_id values for KO-EC relations: {key}")

            relations.append(
                {
                    "sample_id": ko_claim["sample_id"],
                    "feature_uid": ko_claim["feature_uid"],
                    "protein_id": ko_claim["protein_id"],
                    "ko_id": ko_claim["claim_id"],
                    "ec_id": ec_claim["claim_id"],
                    "relation": classify_ko_ec_relation(
                        ko_claim["claim_id"],
                        ec_claim["claim_id"],
                        ko_ec_map,
                    ),
                    "reference_id": reference_id,
                }
            )

    relations.sort(
        key=lambda row: (
            row["sample_id"],
            row["feature_uid"],
            row["ko_id"],
            row["ec_id"],
        )
    )

    return relations


def classify_ec_go_relation(
    ec_id: str,
    go_id: str,
    terms: dict[str, GoTerm],
    alt_to_primary: dict[str, str],
) -> str:
    """Classify an EC against explicit xrefs of a resolved GO-MF term"""

    ec_pattern = r"EC:[0-9]+(?:\.(?:[0-9]+|-)){3}"

    if re.fullmatch(ec_pattern, ec_id) is None:
        raise ValueError(f"Expected canonical EC identifier, got {ec_id!r}")

    resolution = resolve_go_id(go_id, terms, alt_to_primary)

    if not resolution.hierarchy_eligible:
        return "obsolete_unresolved"

    if resolution.namespace != "molecular_function":
        return "not_molecular_function"

    term = terms[resolution.resolved_id]
    reference_relations: set[str] = set()

    for reference_ec in sorted(term.ec_xrefs):
        if re.fullmatch(ec_pattern, reference_ec) is None:
            warnings.warn(
                f"Skipping unsupported EC identifier {reference_ec!r} referenced by {term.go_id}",
                UserWarning,
                stacklevel=2,
            )
            continue

        reference_relations.add(classify_ec_relation(ec_id, reference_ec))

    if "exact" in reference_relations:
        return "documented"

    if "hierarchical_compatible" in reference_relations:
        return "hierarchical_compatible"

    if reference_relations:
        return "not_documented"

    return "unmapped_go"


def build_ec_go_relations(
    ec_claims: list[dict[str, str]],
    go_claims: list[dict[str, str]],
    terms: dict[str, GoTerm],
    alt_to_primary: dict[str, str],
    reference_id: str,
) -> list[dict[str, str]]:
    """Relate observed EC and GO-MF claims on the same protein."""

    if not reference_id.strip():
        raise ValueError("Missing EC-GO reference_id")

    ecs_by_protein: dict[tuple[str, str, str], set[str]] = defaultdict(set)

    for claim in ec_claims:
        if claim["claim_type"] != "EC":
            raise ValueError(f"Expected an EC claim, got {claim['claim_type']!r}")

        if not claim["protein_id"].strip():
            continue

        key = (claim["sample_id"], claim["feature_uid"], claim["protein_id"])
        ecs_by_protein[key].add(claim["claim_id"])

    relations: list[dict[str, str]] = []

    for claim in go_claims:
        if claim["claim_type"] != "GO":
            raise ValueError(f"Expected a GO claim, got {claim['claim_type']!r}")

        if not claim["protein_id"].strip():
            continue

        key = (claim["sample_id"], claim["feature_uid"], claim["protein_id"])

        if key not in ecs_by_protein:
            continue

        resolution = resolve_go_id(claim["claim_id"], terms, alt_to_primary)

        if resolution.namespace != "molecular_function":
            continue

        for ec_id in sorted(ecs_by_protein[key]):
            relations.append(
                {
                    "sample_id": claim["sample_id"],
                    "feature_uid": claim["feature_uid"],
                    "protein_id": claim["protein_id"],
                    "ec_id": ec_id,
                    "go_id": claim["claim_id"],
                    "resolved_go_id": resolution.resolved_id,
                    "go_resolution_status": resolution.status,
                    "relation": classify_ec_go_relation(
                        ec_id, claim["claim_id"], terms, alt_to_primary
                    ),
                    "reference_id": reference_id,
                }
            )

    relations.sort(
        key=lambda row: (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["ec_id"],
            row["go_id"],
        )
    )

    return relations
