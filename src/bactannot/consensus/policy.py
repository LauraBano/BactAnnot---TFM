"""Explicit semantic policy for annotation namespaces."""

CONSENSUS = "consensus"
SPECIALIST_CLAIM = "specialist_claim"
CONTEXT = "context"
NORMALIZATION_PENDING = "normalization_pending"


ANNOTATION_TYPE_POLICY = {
    # --------------------------------------------------------
    # Evidence-aware consensus already implemented
    # --------------------------------------------------------
    "ARO": CONSENSUS,
    "CDD": CONSENSUS,
    "COG": CONSENSUS,
    "EC": CONSENSUS,
    "GO": CONSENSUS,
    "KEGG_KO": CONSENSUS,
    "NCBIFAM": CONSENSUS,
    "PFAM": CONSENSUS,
    "SMART": CONSENSUS,
    "TOXIN_ASSOCIATED": CONSENSUS,
    "VIRULENCE_ASSOCIATED": CONSENSUS,
    # --------------------------------------------------------
    # Valid biological claims without a comparable
    # multi-source consensus in the current model
    # --------------------------------------------------------
    "ANTIFAM": SPECIALIST_CLAIM,
    "CAZY": SPECIALIST_CLAIM,
    "CDD_PROFILE": SPECIALIST_CLAIM,
    "COILS": SPECIALIST_CLAIM,
    "FUNFAM": SPECIALIST_CLAIM,
    "GENE3D": SPECIALIST_CLAIM,
    "HAMAP": SPECIALIST_CLAIM,
    "INTERPRO": SPECIALIST_CLAIM,
    "KEGG_REACTION": SPECIALIST_CLAIM,
    "KEGG_TC": SPECIALIST_CLAIM,
    "KOG": SPECIALIST_CLAIM,
    "MOBIDBLITE": SPECIALIST_CLAIM,
    "PANTHER": SPECIALIST_CLAIM,
    "PFAM_NAME": SPECIALIST_CLAIM,
    "PIRSF": SPECIALIST_CLAIM,
    "PIRSR": SPECIALIST_CLAIM,
    "PRINTS": SPECIALIST_CLAIM,
    "PROSITE": SPECIALIST_CLAIM,
    "RFAM": SPECIALIST_CLAIM,
    "SFLD": SPECIALIST_CLAIM,
    "SUPERFAMILY": SPECIALIST_CLAIM,
    # --------------------------------------------------------
    # Hierarchical or interpretative functional context.
    # These annotations are preserved but must not increase
    # the number of independent evidence lines.
    # --------------------------------------------------------
    "COG_CATEGORY": CONTEXT,
    "EGGNOG_OG": CONTEXT,
    "KEGG_BRITE": CONTEXT,
    "KEGG_MODULE": CONTEXT,
    "KEGG_PATHWAY": CONTEXT,
    "KEGG_RCLASS": CONTEXT,
    # --------------------------------------------------------
    # Preserved annotation requiring semantic normalization
    # before it can safely participate in downstream support
    # interpretation.
    # --------------------------------------------------------
    "BIGG_REACTION": NORMALIZATION_PENDING,
}


def get_annotation_type_policy(
    annotation_type: str,
) -> str:
    """Return the explicit semantic policy for an annotation type."""

    try:
        return ANNOTATION_TYPE_POLICY[annotation_type]

    except KeyError as exc:
        raise ValueError(f"No annotation policy defined for {annotation_type!r}") from exc


def validate_annotation_type_coverage(
    annotation_rows: list[dict[str, str]],
) -> None:
    """Ensure every observed annotation type has an explicit policy."""

    observed_types = {row["annotation_type"] for row in annotation_rows if row["annotation_type"]}

    unknown_types = sorted(observed_types - ANNOTATION_TYPE_POLICY.keys())

    if unknown_types:
        raise ValueError(
            "Annotation types without an explicit semantic policy: " + ", ".join(unknown_types)
        )
