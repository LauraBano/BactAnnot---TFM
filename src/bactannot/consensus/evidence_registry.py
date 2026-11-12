"""Evidence provenance registry for functional consensus"""

from dataclasses import dataclass


@dataclass(frozen=True)
class EvidenceChannel:
    """Describe the methodological provenance of one evidence channel"""

    source: str
    producer: str
    evidence_type: str
    method_group: str
    resource_group: str
    native_acceptance_rule: str


EGGNOG_KO_CHANNEL = EvidenceChannel(
    source="eggNOG-mapper",
    producer="eggNOG-mapper",
    evidence_type="seed_ortholog",
    method_group="orthology_transfer",
    resource_group="eggNOG",
    native_acceptance_rule="native_orthology_mapping",
)


KOFAM_KO_CHANNEL = EvidenceChannel(
    source="KofamScan",
    producer="KofamScan",
    evidence_type="hmm_match",
    method_group="profile_hmm",
    resource_group="KOfam",
    native_acceptance_rule="family_specific_threshold",
)


BAKTA_BLASTRULES_CHANNEL = EvidenceChannel(
    source="BlastRules",
    producer="bakta",
    evidence_type="expert",
    method_group="sequence_homology",
    resource_group="BlastRules",
    native_acceptance_rule="native_expert_rule_match",
)


EVIDENCE_REGISTRY = {
    (
        EGGNOG_KO_CHANNEL.source,
        EGGNOG_KO_CHANNEL.evidence_type,
    ): EGGNOG_KO_CHANNEL,
    (
        KOFAM_KO_CHANNEL.source,
        KOFAM_KO_CHANNEL.evidence_type,
    ): KOFAM_KO_CHANNEL,
    (
        BAKTA_BLASTRULES_CHANNEL.source,
        BAKTA_BLASTRULES_CHANNEL.evidence_type,
    ): BAKTA_BLASTRULES_CHANNEL,
}


def get_evidence_channel(
    source: str,
    evidence_type: str,
) -> EvidenceChannel | None:
    return EVIDENCE_REGISTRY.get(
        (
            source,
            evidence_type,
        )
    )
