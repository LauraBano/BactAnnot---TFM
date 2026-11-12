"""Functional interpretation of evidence-aware claims"""

import json
import re
from collections import defaultdict

from bactannot.consensus.biomedical import (
    BIOMEDICAL_CLASSIFICATION_TYPES,
    build_aro_evidence_index,
    build_biomedical_classification_evidence_index,
    collect_aro_claim_evidence,
    collect_biomedical_classification_evidence,
)
from bactannot.consensus.cog import (
    build_cog_evidence_indexes,
    collect_cog_claim_evidence,
)
from bactannot.consensus.domains import (
    DOMAIN_CLAIM_TYPES,
    build_domain_evidence_index,
    canonicalize_domain_evidence,
    collect_domain_claim_evidence,
)
from bactannot.consensus.ec import (
    collect_ec_evidence,
    get_ec_claim_relation,
    parse_ec_components,
)
from bactannot.consensus.go import (
    collect_eggnog_go_evidence,
    collect_interproscan_go_evidence,
)
from bactannot.consensus.go_ontology import (
    GoTerm,
    build_go_ancestor_map,
    reduce_to_most_specific,
    resolve_go_id,
)
from bactannot.consensus.ko import (
    collect_direct_ko_evidence,
    collect_transferred_ko_evidence,
)
from bactannot.schema import EVIDENCE_COLUMNS

UNSPECIFIED_FUNCTION_RE = re.compile(
    r"\b(?:hypothetical|uncharacteri[sz]ed|unknown|DUF\d+|UPF\d+)\b",
    re.IGNORECASE,
)


def group_claim_evidence(
    claim_evidence_pairs: list[tuple[dict[str, str], list[dict[str, str]]]],
) -> dict[tuple[str, ...], set[tuple[str, str]]]:
    """Group claims by their shared canonical evidence records"""

    claims_by_evidence: dict[
        tuple[str, ...],
        set[tuple[str, str]],
    ] = {}

    for claim, evidence_rows in claim_evidence_pairs:
        claim_feature = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
        )

        claim_key = (
            claim["claim_type"],
            claim["claim_id"],
        )

        for row in evidence_rows:
            evidence_feature = (
                row["sample_id"],
                row["feature_uid"],
                row["protein_id"],
            )

            if evidence_feature != claim_feature:
                raise ValueError("Claim and evidence refer to different biological features")

            evidence_key = tuple(row[column] for column in EVIDENCE_COLUMNS)

            claims_by_evidence.setdefault(
                evidence_key,
                set(),
            ).add(claim_key)

    return dict(sorted(claims_by_evidence.items()))


def collect_claim_evidence_pairs(
    claim_rows: list[dict[str, str]],
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[tuple[dict[str, str], list[dict[str, str]]]]:
    """Collect native evidence for claims using their existing resolvers."""

    observations_by_claim: dict[
        tuple[str, ...],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in annotation_rows:
        if row["annotation_type"] not in {"EC", "GO"}:
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["annotation_type"],
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

    domain_index = build_domain_evidence_index(evidence_rows)

    cog_eggnog_index, cog_pathofact_index = build_cog_evidence_indexes(
        evidence_rows,
    )

    aro_index = build_aro_evidence_index(evidence_rows)

    biomedical_index = build_biomedical_classification_evidence_index(
        evidence_rows,
    )

    pairs: list[tuple[dict[str, str], list[dict[str, str]]]] = []

    for claim in claim_rows:
        feature_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
        )

        claim_type = claim["claim_type"]

        claim_key = feature_key + (
            claim_type,
            claim["claim_id"],
        )

        observations = observations_by_claim.get(claim_key, [])
        candidate_evidence = evidence_by_feature.get(feature_key, [])

        if claim_type == "KEGG_KO":
            matches = collect_direct_ko_evidence(
                claim,
                candidate_evidence,
            )

            matches.extend(
                collect_transferred_ko_evidence(
                    claim,
                    candidate_evidence,
                )
            )

        elif claim_type == "EC":
            matches = collect_ec_evidence(
                claim,
                observations,
                candidate_evidence,
            )

        elif claim_type == "GO":
            matches = collect_eggnog_go_evidence(
                claim,
                observations,
                candidate_evidence,
            )

            matches.extend(
                collect_interproscan_go_evidence(
                    claim,
                    observations,
                    candidate_evidence,
                )
            )

        elif claim_type in DOMAIN_CLAIM_TYPES:
            matches = collect_domain_claim_evidence(
                claim,
                domain_index,
            )

        elif claim_type == "COG":
            matches = collect_cog_claim_evidence(
                claim,
                cog_eggnog_index,
                cog_pathofact_index,
            )

        elif claim_type == "ARO":
            matches = collect_aro_claim_evidence(
                claim,
                aro_index,
            )

        elif claim_type in BIOMEDICAL_CLASSIFICATION_TYPES:
            matches = collect_biomedical_classification_evidence(
                claim,
                biomedical_index,
            )

        else:
            raise ValueError(f"Unsupported claim type for evidence collection: {claim_type!r}")

        pairs.append((claim, matches))

    return pairs


def build_enzyme_activity_candidates(
    ec_claims: list[dict[str, str]],
    ko_ec_relations: list[dict[str, str]],
    ec_go_relations: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build observed EC activity candidates with documented KO and MF GO links."""

    ko_relations_by_ec: dict[
        tuple[str, str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in ko_ec_relations:
        if row["relation"] != "documented":
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["ec_id"],
        )

        ko_relations_by_ec[key].append(row)

    go_relations_by_ec: dict[
        tuple[str, str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in ec_go_relations:
        if row["relation"] != "documented":
            continue

        key = (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["ec_id"],
        )

        go_relations_by_ec[key].append(row)

    candidates: list[dict[str, str]] = []

    for claim in ec_claims:
        if claim["claim_type"] != "EC":
            raise ValueError(f"Expected an EC claim, got {claim['claim_type']!r}")

        key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
            claim["claim_id"],
        )

        ko_rows = ko_relations_by_ec.get(key, [])
        go_rows = go_relations_by_ec.get(key, [])

        ko_ids = {row["ko_id"] for row in ko_rows}
        mf_go_ids = {row["go_id"] for row in go_rows}

        reference_ids = {row["reference_id"] for row in ko_rows + go_rows}

        candidates.append(
            {
                **claim,
                "documented_ko_ids": ";".join(sorted(ko_ids)),
                "documented_mf_go_ids": ";".join(sorted(mf_go_ids)),
                "relation_reference_ids": ";".join(sorted(reference_ids)),
            }
        )

    candidates.sort(
        key=lambda row: (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["claim_id"],
        )
    )

    return candidates


def build_molecular_function_candidates(
    go_claims: list[dict[str, str]],
    terms: dict[str, GoTerm],
    alt_to_primary: dict[str, str],
    reference_id: str,
) -> list[dict[str, str]]:
    """Prepare resolved molecular-function proposals from original GO claims."""

    candidates: list[dict[str, str]] = []

    for claim in go_claims:
        if claim["claim_type"] != "GO":
            raise ValueError(f"Expected a GO claim, got {claim['claim_type']!r}")

        resolution = resolve_go_id(
            claim["claim_id"],
            terms,
            alt_to_primary,
        )

        if not resolution.hierarchy_eligible:
            continue

        if resolution.namespace != "molecular_function":
            continue

        if resolution.resolved_id == "GO:0003674":
            continue

        candidates.append(
            {
                **claim,
                "resolved_go_id": resolution.resolved_id,
                "go_resolution_status": resolution.status,
                "reference_id": reference_id,
            }
        )

    candidates.sort(
        key=lambda row: (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["resolved_go_id"],
            row["claim_id"],
        )
    )

    return candidates


def build_family_function_candidates(
    claim_evidence_pairs: list[tuple[dict[str, str], list[dict[str, str]]]],
) -> list[dict[str, str]]:
    """Prepare ortholog-group and signature proposals with native descriptions."""

    candidates: list[dict[str, str]] = []

    for claim, evidence_rows in claim_evidence_pairs:
        claim_type = claim["claim_type"]

        if claim_type not in {"KEGG_KO", "COG"} and claim_type not in DOMAIN_CLAIM_TYPES:
            continue

        if claim_type in DOMAIN_CLAIM_TYPES:
            candidate_scope = "domain_or_family"
        else:
            candidate_scope = "ortholog_group"

        descriptions: set[tuple[str, str, str]] = set()

        for row in evidence_rows:
            ko_definition = (
                claim_type == "KEGG_KO"
                and row["source"] == "KofamScan"
                and row["evidence_type"] == "hmm_match"
                and row["accession"] == claim["claim_id"]
            )

            signature_definition = (
                claim_type in DOMAIN_CLAIM_TYPES
                and row["source"].startswith("InterProScan:")
                and row["evidence_type"] == "signature_match"
            )

            pfam_definition = (
                claim_type == "PFAM" and row["source"] == "Pfam" and row["evidence_type"] == "pfam"
            )

            if not (ko_definition or signature_definition or pfam_definition):
                continue

            description = row["product"].strip()

            if description:
                descriptions.add(
                    (
                        row["source"],
                        row["accession"],
                        description,
                    )
                )

        description_records = [
            {
                "source": source,
                "accession": accession,
                "description": description,
            }
            for source, accession, description in sorted(descriptions)
        ]

        candidates.append(
            {
                **claim,
                "candidate_scope": candidate_scope,
                "native_descriptions": json.dumps(
                    description_records,
                    ensure_ascii=False,
                ),
            }
        )

    candidates.sort(
        key=lambda row: (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["claim_type"],
            row["claim_id"],
        )
    )

    return candidates


def select_convergent_family_function_proposals(
    claim_evidence_pairs: list[tuple[dict[str, str], list[dict[str, str]]]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Select signature descriptions that agree with transferred products."""

    blocked_features: set[tuple[str, str, str]] = set()
    family_pairs: list[tuple[dict[str, str], list[dict[str, str]]]] = []

    for claim, matches in claim_evidence_pairs:
        feature_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
        )

        if not claim["protein_id"].strip():
            continue

        if claim["claim_type"] in {"KEGG_KO", "EC"}:
            blocked_features.add(feature_key)

        elif claim["claim_type"] in DOMAIN_CLAIM_TYPES:
            matching_rows = [
                row
                for row in matches
                if (
                    row["sample_id"],
                    row["feature_uid"],
                    row["protein_id"],
                )
                == feature_key
                and canonicalize_domain_evidence(row) == (claim["claim_type"], claim["claim_id"])
            ]

            family_pairs.append((claim, matching_rows))

    transferred_by_feature: dict[
        tuple[str, str, str],
        set[tuple[str, str, str]],
    ] = defaultdict(set)

    for row in evidence_rows:
        if row["source"] != "eggNOG-mapper" or row["evidence_type"] != "seed_ortholog":
            continue

        accession = row["accession"].strip()
        description = row["product"].strip()

        if not row["protein_id"].strip() or not accession or not description:
            continue

        feature_key = (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
        )

        transferred_by_feature[feature_key].add((row["source"], accession, description))

    proposals: list[dict[str, str]] = []

    for candidate in build_family_function_candidates(family_pairs):
        feature_key = (
            candidate["sample_id"],
            candidate["feature_uid"],
            candidate["protein_id"],
        )

        if feature_key in blocked_features:
            continue

        native_records = json.loads(candidate["native_descriptions"])

        native_products = {row["description"].strip().casefold() for row in native_records}

        transferred_records = transferred_by_feature.get(feature_key, set())

        transferred_products = {description.casefold() for _, _, description in transferred_records}

        if len(native_products) != 1 or native_products != transferred_products:
            continue

        description = native_records[0]["description"].strip()

        if description.casefold() in {
            "protein",
            "conserved protein",
            "predicted protein",
            "putative protein",
        }:
            continue

        if UNSPECIFIED_FUNCTION_RE.search(description):
            continue

        proposal = dict(candidate)
        proposal["proposed_function"] = description
        proposal["description_source"] = ";".join(sorted({row["source"] for row in native_records}))
        proposal["selection_rule"] = "signature_transferred_product_agreement"
        proposal["transferred_descriptions"] = json.dumps(
            [
                {
                    "source": source,
                    "accession": accession,
                    "description": product,
                }
                for source, accession, product in sorted(transferred_records)
            ],
            ensure_ascii=False,
        )

        proposals.append(proposal)

    return proposals


def select_single_ko_function_proposals(
    claim_evidence_pairs: list[tuple[dict[str, str], list[dict[str, str]]]],
) -> list[dict[str, str]]:
    """Select functions from the only KO with direct and transferred support"""

    ko_pairs_by_feature: dict[
        tuple[str, str, str],
        list[tuple[dict[str, str], list[dict[str, str]]]],
    ] = defaultdict(list)

    for claim, matches in claim_evidence_pairs:
        if claim["claim_type"] != "KEGG_KO":
            continue

        feature_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
        )

        ko_pairs_by_feature[feature_key].append((claim, matches))

    proposals: list[dict[str, str]] = []

    for feature_key in sorted(ko_pairs_by_feature):
        pairs = ko_pairs_by_feature[feature_key]

        convergent_pairs: list[tuple[dict[str, str], list[dict[str, str]]]] = []

        for claim, matches in pairs:
            if not collect_direct_ko_evidence(claim, matches):
                continue

            if not collect_transferred_ko_evidence(claim, matches):
                continue

            convergent_pairs.append((claim, matches))

        if len(convergent_pairs) != 1:
            continue

        claim, matches = convergent_pairs[0]

        candidate = build_family_function_candidates(convergent_pairs)[0]

        descriptions = {
            row["description"].strip()
            for row in json.loads(candidate["native_descriptions"])
            if row["source"] == "KofamScan"
            and row["accession"] == claim["claim_id"]
            and row["description"].strip()
        }

        if len(descriptions) != 1:
            continue

        description = next(iter(descriptions))

        if description.casefold() in {
            "protein",
            "conserved protein",
            "predicted protein",
            "putative protein",
        }:
            continue

        if UNSPECIFIED_FUNCTION_RE.search(description):
            continue

        proposal = dict(candidate)
        proposal["proposed_function"] = description
        proposal["description_source"] = "KofamScan"
        if len(pairs) == 1:
            proposal["selection_rule"] = "single_ko_direct_and_transferred"
        else:
            proposal["selection_rule"] = "single_convergent_ko_with_alternatives"

        proposals.append(proposal)

    return proposals


def select_convergent_ec_activity_proposals(
    activity_candidates: list[dict[str, str]],
    claim_evidence_pairs: list[tuple[dict[str, str], list[dict[str, str]]]],
) -> list[dict[str, str]]:
    """Select complete EC activities with KO-derived and transferred support."""

    evidence_by_claim: dict[
        tuple[str, str, str, str],
        list[dict[str, str]],
    ] = {}

    for claim, matches in claim_evidence_pairs:
        if claim["claim_type"] != "EC":
            continue

        claim_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
            claim["claim_id"],
        )

        evidence_by_claim[claim_key] = matches

    proposals: list[dict[str, str]] = []

    for candidate in activity_candidates:
        if candidate["claim_type"] != "EC":
            raise ValueError("Expected an EC activity candidate")

        if "-" in parse_ec_components(candidate["claim_id"]):
            continue

        claim_key = (
            candidate["sample_id"],
            candidate["feature_uid"],
            candidate["protein_id"],
            candidate["claim_id"],
        )

        matches = evidence_by_claim.get(claim_key, [])

        relations = {get_ec_claim_relation(row["source"], row["evidence_type"]) for row in matches}

        if not {
            "derived_from_ko",
            "transferred_annotation",
        }.issubset(relations):
            continue

        proposal = dict(candidate)
        proposal["selection_rule"] = "ec_ko_derived_and_transferred"
        proposals.append(proposal)

    proposals.sort(
        key=lambda row: (
            row["sample_id"],
            row["feature_uid"],
            row["protein_id"],
            row["claim_id"],
        )
    )

    return proposals


def select_convergent_molecular_function_proposals(
    mf_candidates: list[dict[str, str]],
    terms: dict[str, GoTerm],
) -> list[dict[str, str]]:
    """Select supported MF concepts and reduce eligible is_a ancestors."""

    if not mf_candidates:
        return []

    reference_ids = {row["reference_id"] for row in mf_candidates}

    if len(reference_ids) != 1:
        raise ValueError("MF candidates must use one GO reference")

    concepts: dict[
        tuple[str, str, str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for candidate in mf_candidates:
        if candidate["claim_type"] != "GO":
            raise ValueError("Expected a molecular function candidate")

        concept_key = (
            candidate["sample_id"],
            candidate["feature_uid"],
            candidate["protein_id"],
            candidate["resolved_go_id"],
        )

        concepts[concept_key].append(candidate)

    eligible_by_feature: dict[
        tuple[str, str, str],
        set[str],
    ] = defaultdict(set)

    for concept_key, observations in concepts.items():
        relations: set[str] = set()

        for candidate in observations:
            relations.update(candidate["evidence_relations"].split(";"))

        if not {
            "transferred_annotation",
            "mapped_from_signature",
        }.issubset(relations):
            continue

        sample_id, feature_uid, protein_id, go_id = concept_key
        feature_key = (sample_id, feature_uid, protein_id)

        eligible_by_feature[feature_key].add(go_id)

    is_a_map = build_go_ancestor_map(terms, is_a_only=True)

    proposals: list[dict[str, str]] = []

    for feature_key in sorted(eligible_by_feature):
        selected_go_ids = reduce_to_most_specific(
            eligible_by_feature[feature_key],
            is_a_map,
        )

        sample_id, feature_uid, protein_id = feature_key

        for go_id in sorted(selected_go_ids):
            concept_key = (sample_id, feature_uid, protein_id, go_id)
            observations = concepts[concept_key]

            original_go_ids = {candidate["claim_id"] for candidate in observations}

            proposals.append(
                {
                    "sample_id": sample_id,
                    "feature_uid": feature_uid,
                    "protein_id": protein_id,
                    "resolved_go_id": go_id,
                    "original_go_ids": ";".join(sorted(original_go_ids)),
                    "proposed_function": terms[go_id].name or go_id,
                    "selection_rule": "mf_transfer_and_signature_is_a_reduced",
                    "reference_id": observations[0]["reference_id"],
                }
            )

    return proposals


def build_protein_function_interpretations(
    sample_id: str,
    feature_index: dict[str, str],
    ko_proposals: list[dict[str, str]],
    activity_proposals: list[dict[str, str]],
    mf_proposals: list[dict[str, str]],
    family_proposals: list[dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    """Assemble protein summaries from previously selected proposals."""

    if family_proposals is None:
        family_proposals = []

    for proposal in family_proposals:
        if (
            proposal["claim_type"] not in DOMAIN_CLAIM_TYPES
            or proposal["candidate_scope"] != "domain_or_family"
            or not proposal["proposed_function"].strip()
        ):
            raise ValueError("Invalid domain or family description proposal")

    feature_keys = {
        (sample_id, feature_uid, protein_id)
        for feature_uid, protein_id in feature_index.items()
        if protein_id
    }

    ko_index: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    ec_index: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    mf_index: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    family_index: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)

    for proposals, index in (
        (ko_proposals, ko_index),
        (activity_proposals, ec_index),
        (mf_proposals, mf_index),
        (family_proposals, family_index),
    ):
        for proposal in proposals:
            feature_key = (
                proposal["sample_id"],
                proposal["feature_uid"],
                proposal["protein_id"],
            )

            if feature_key not in feature_keys:
                raise ValueError("Proposal does not match the sample features")

            index[feature_key].append(proposal)

    interpretations: list[dict[str, str]] = []

    for feature_key in sorted(feature_keys):
        row_sample_id, feature_uid, protein_id = feature_key
        kos = ko_index.get(feature_key, [])
        activities = ec_index.get(feature_key, [])
        functions = mf_index.get(feature_key, [])
        families = family_index.get(feature_key, [])

        if len(kos) > 1:
            raise ValueError("Expected at most one KO function proposal per protein")

        row = {
            "sample_id": row_sample_id,
            "feature_uid": feature_uid,
            "protein_id": protein_id,
            "interpretation_status": "no_selected_proposal",
            "primary_claim_type": "",
            "primary_claim_id": "",
            "primary_function": "",
            "primary_scope": "",
            "primary_selection_rule": "",
            "primary_description_source": "",
            "primary_description_reference_id": "",
            "ko_proposal_ids": ";".join(sorted({item["claim_id"] for item in kos})),
            "ec_activity_proposal_ids": ";".join(sorted({item["claim_id"] for item in activities})),
            "mf_proposal_ids": ";".join(sorted({item["resolved_go_id"] for item in functions})),
            "unlinked_ec_activity_ids": "",
            "family_proposal_ids": ";".join(sorted({item["claim_id"] for item in families})),
        }

        if kos:
            proposal = kos[0]
            unlinked_ec_ids = {
                activity["claim_id"]
                for activity in activities
                if proposal["claim_id"] not in activity["documented_ko_ids"].split(";")
            }

            if unlinked_ec_ids:
                row["interpretation_status"] = "ko_ec_relation_unresolved"
                row["unlinked_ec_activity_ids"] = ";".join(sorted(unlinked_ec_ids))
            else:
                row["interpretation_status"] = "ko_function_proposed"
                row["primary_claim_type"] = "KEGG_KO"
                row["primary_claim_id"] = proposal["claim_id"]
                row["primary_function"] = proposal["proposed_function"]
                row["primary_scope"] = "protein_function"
                row["primary_selection_rule"] = proposal["selection_rule"]
                row["primary_description_source"] = proposal["description_source"]

        elif len(activities) == 1:
            proposal = activities[0]
            row["interpretation_status"] = "ec_activity_proposed"
            row["primary_claim_type"] = "EC"
            row["primary_claim_id"] = proposal["claim_id"]
            row["primary_function"] = proposal["claim_id"]
            row["primary_scope"] = "enzyme_activity"
            row["primary_selection_rule"] = proposal["selection_rule"]

            matching_functions = [
                function
                for function in functions
                if function["resolved_go_id"] in proposal["documented_mf_go_ids"].split(";")
                and function["reference_id"] in proposal["relation_reference_ids"].split(";")
            ]

            if len(matching_functions) == 1:
                function = matching_functions[0]
                row["primary_function"] = function["proposed_function"]
                row["primary_description_source"] = function["resolved_go_id"]
                row["primary_description_reference_id"] = function["reference_id"]

        elif activities:
            row["interpretation_status"] = "multiple_ec_activities_proposed"

        elif families:
            descriptions = {item["proposed_function"].strip() for item in families}

            if len({description.casefold() for description in descriptions}) > 1:
                row["interpretation_status"] = "multiple_family_descriptions_proposed"

            else:
                row["interpretation_status"] = "family_description_proposed"
                row["primary_function"] = min(descriptions)
                row["primary_scope"] = "domain_or_family"

                row["primary_selection_rule"] = ";".join(
                    sorted({item["selection_rule"] for item in families})
                )

                row["primary_description_source"] = ";".join(
                    sorted(
                        {
                            source
                            for item in families
                            for source in item["description_source"].split(";")
                            if source
                        }
                    )
                )

                claim_keys = {(item["claim_type"], item["claim_id"]) for item in families}

                if len(claim_keys) == 1:
                    row["primary_claim_type"], row["primary_claim_id"] = next(iter(claim_keys))

        elif functions:
            row["interpretation_status"] = "molecular_function_properties_only"

        interpretations.append(row)

    return interpretations
