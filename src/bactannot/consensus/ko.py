"""Canonical KEGG KO claims for evidence-aware consensus"""

from collections import defaultdict

from bactannot.consensus.evidence_registry import get_evidence_channel

KO_CONTEXT_SOURCES = {
    "bakta",
}

KO_CLAIM_RELATIONS = {
    (
        "eggNOG-mapper",
        "seed_ortholog",
    ): "transferred_annotation",
    (
        "KofamScan",
        "hmm_match",
    ): "direct_match",
}


def get_ko_claim_relation(
    source: str,
    evidence_type: str,
) -> str | None:
    return KO_CLAIM_RELATIONS.get(
        (
            source,
            evidence_type,
        )
    )


# Convierte un conjunto de herramientas en una única cadena de texto ordenada
def serialize_sources(
    # genera un conjunto desordenado de elementos únicos
    sources: set[str],
) -> str:
    """Serialize distinct annotation sources deterministically"""

    return ";".join(
        # ordena el conjunto y devuelve una lista
        sorted(
            sources,
            # permite comparación entre mayúsculas y minúsculas
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


# Genera una lista de diccionarios con una fila por afirmación funcional única
def build_ko_claims_from_annotations(
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build canonical KEGG KO claims from annotation observations"""

    # crea un diccionario para agrupar anotaciones
    observations_by_claim: dict[
        # genera una tupla con los campos sample_id, feature_id, annotation_id
        tuple[str, str, str],
        # los valores de la clave (tupla) son listas de diccionarios, cada lista corresponde a una observación
        list[dict[str, str]],
    ] = defaultdict(list)

    # bucle recorre individualmente las filas de annotations.tsv, en cada iteración row contiene un diccionario
    for row in annotation_rows:
        if row["annotation_type"] != "KEGG_KO":
            continue

        # construye la clave tupla con tres elementos que en conjunto identifica una afirmación funcional concreta
        key = (
            row["sample_id"],
            row["feature_uid"],
            row["annotation_id"],
        )

        # crea una entrada en el diccionario para cada nueva claim, observaciones (fuentes) diferentes de una misma claim se añaden como valor (lista) de la misma claim
        observations_by_claim[key].append(row)

    # crea una lista vacía para llenarla con los valores del diccionario
    claims: list[dict[str, str]] = []

    # recorre simultáneamente las claves y valores del diccionario
    for key, observations in observations_by_claim.items():
        # desempaqueta la tupla y asigna cada elemento a una variable
        (
            sample_id,
            feature_uid,
            claim_id,
        ) = key

        # comprensión de conjuntos, construye un conjunto a partir de una colección
        # obtiene los protein_id que aparecen dentro de un mismo grupo de observaciones
        protein_ids = {row["protein_id"] for row in observations}

        # comprueba que se han agrupado las observaciones por protein_id y no varía
        if len(protein_ids) != 1:
            raise ValueError(f"Conflicting protein_id values for KEGG KO claim: {key}")

        # next() recupera el siguiente elemento, iter() itera sobre el conjunto de protein_ids
        protein_id = next(iter(protein_ids))

        # recorre las observaciones y construye un conjunto con los nombres de las herramientas
        sources = {row["source"] for row in observations}

        # construye el diccionario del claim
        claims.append(
            {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "claim_type": "KEGG_KO",
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


def collect_direct_ko_evidence(
    claim: dict[str, str],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Collect evidence records that directly support one KEGG KO claim"""

    if claim["claim_type"] != "KEGG_KO":
        raise ValueError(f"Expected a KEGG_KO claim, got {claim['claim_type']!r}")

    matches: list[dict[str, str]] = []

    for row in evidence_rows:
        channel = get_evidence_channel(
            row["source"],
            row["evidence_type"],
        )

        if channel is None:
            continue

        claim_relation = get_ko_claim_relation(
            row["source"],
            row["evidence_type"],
        )

        if claim_relation != "direct_match":
            continue

        if (
            row["sample_id"] == claim["sample_id"]
            and row["feature_uid"] == claim["feature_uid"]
            and row["protein_id"] == claim["protein_id"]
            and row["accession"] == claim["claim_id"]
        ):
            matches.append(row)

    matches.sort(
        key=lambda row: (
            row["source"],
            row["evidence_type"],
            row["accession"],
        )
    )

    return matches


def collect_transferred_ko_evidence(
    claim: dict[str, str],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Collect transferred evidence supporting one KEGG KO claim"""

    if claim["claim_type"] != "KEGG_KO":
        raise ValueError(f"Expected a KEGG_KO claim, got {claim['claim_type']!r}")

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    matches: list[dict[str, str]] = []

    for row in evidence_rows:
        channel = get_evidence_channel(
            row["source"],
            row["evidence_type"],
        )

        if channel is None:
            continue

        claim_relation = get_ko_claim_relation(
            row["source"],
            row["evidence_type"],
        )

        if claim_relation != "transferred_annotation":
            continue

        if channel.producer not in observing_sources:
            continue

        if (
            row["sample_id"] == claim["sample_id"]
            and row["feature_uid"] == claim["feature_uid"]
            and row["protein_id"] == claim["protein_id"]
        ):
            matches.append(row)

    matches.sort(
        key=lambda row: (
            row["source"],
            row["evidence_type"],
            row["accession"],
        )
    )

    return matches


def summarize_ko_claim_evidence(
    claim: dict[str, str],
    evidence_rows: list[dict[str, str]],
) -> dict[str, str]:
    """Summarize linked methodological evidence for one KEGG KO claim"""

    direct_evidence = collect_direct_ko_evidence(
        claim,
        evidence_rows,
    )

    transferred_evidence = collect_transferred_ko_evidence(
        claim,
        evidence_rows,
    )

    matched_evidence = direct_evidence + transferred_evidence

    evidence_lines: set[tuple[str, str]] = set()
    evidence_sources: set[str] = set()
    resolved_producers: set[str] = set()

    for row in matched_evidence:
        channel = get_evidence_channel(
            row["source"],
            row["evidence_type"],
        )

        if channel is None:
            continue

        evidence_lines.add(
            (
                channel.method_group,
                channel.resource_group,
            )
        )

        evidence_sources.add(channel.source)
        resolved_producers.add(channel.producer)

    observing_sources = {source for source in claim["observing_sources"].split(";") if source}

    context_sources = (observing_sources & KO_CONTEXT_SOURCES) - resolved_producers

    unresolved_sources = observing_sources - resolved_producers - context_sources

    return {
        "evidence_sources": serialize_sources(evidence_sources),
        "evidence_lines": serialize_evidence_lines(evidence_lines),
        "independent_evidence_count": str(len(evidence_lines)),
        "context_sources": serialize_sources(context_sources),
        "unresolved_sources": serialize_sources(unresolved_sources),
    }


def classify_ko_evidence_status(
    evidence: dict[str, str],
) -> str:
    """Classify the evidence status of one KEGG KO claim"""

    independent_evidence_count = int(evidence["independent_evidence_count"])

    if independent_evidence_count >= 2:
        return "multiple_independent"

    if independent_evidence_count == 1:
        return "single_independent"

    if evidence["unresolved_sources"]:
        return "unresolved"

    if evidence["context_sources"]:
        return "context_only"

    return "no_evidence"


def build_ko_claims_with_evidence(
    annotation_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build KEGG KO claims enriched with methodological evidence."""

    claims = build_ko_claims_from_annotations(annotation_rows)

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
        evidence_key = (
            claim["sample_id"],
            claim["feature_uid"],
            claim["protein_id"],
        )

        candidate_evidence = evidence_by_feature.get(
            evidence_key,
            [],
        )

        evidence = summarize_ko_claim_evidence(
            claim,
            candidate_evidence,
        )

        evidence_claim = {
            **claim,
            **evidence,
            "evidence_status": (classify_ko_evidence_status(evidence)),
        }

        evidence_claims.append(evidence_claim)

    return evidence_claims


def classify_eggnog_kofam_relation(
    eggnog_kos: set[str],
    kofam_kos: set[str],
) -> str:
    """Classify the relation between eggNOG and Kofam KO sets"""

    if not eggnog_kos and not kofam_kos:
        return "no_eggnog_or_kofam"

    if eggnog_kos and not kofam_kos:
        return "only_eggnog"

    if kofam_kos and not eggnog_kos:
        return "only_kofam"

    if eggnog_kos == kofam_kos:
        return "exact"

    if eggnog_kos & kofam_kos:
        return "partial"

    return "disjoint"


def serialize_ko_ids(
    ko_ids: set[str],
) -> str:
    """Serialize KEGG KO identifiers deterministically"""

    return ";".join(sorted(ko_ids))


def build_eggnog_kofam_comparisons(
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Build per-feature eggNOG versus Kofam KO comparisons"""

    observations_by_feature: dict[
        tuple[str, str],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in annotation_rows:
        if row["annotation_type"] != "KEGG_KO":
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
            raise ValueError(f"Conflicting protein_id values for KEGG KO feature: {key}")

        protein_id = next(iter(protein_ids))

        eggnog_kos = {
            row["annotation_id"] for row in observations if row["source"] == "eggNOG-mapper"
        }

        kofam_kos = {row["annotation_id"] for row in observations if row["source"] == "KofamScan"}

        shared_kos = eggnog_kos & kofam_kos

        eggnog_only_kos = eggnog_kos - kofam_kos

        kofam_only_kos = kofam_kos - eggnog_kos

        comparisons.append(
            {
                "sample_id": sample_id,
                "feature_uid": feature_uid,
                "protein_id": protein_id,
                "eggnog_kos": serialize_ko_ids(eggnog_kos),
                "kofam_kos": serialize_ko_ids(kofam_kos),
                "shared_kos": serialize_ko_ids(shared_kos),
                "eggnog_only_kos": serialize_ko_ids(eggnog_only_kos),
                "kofam_only_kos": serialize_ko_ids(kofam_only_kos),
                "eggnog_kofam_relation": (
                    classify_eggnog_kofam_relation(
                        eggnog_kos,
                        kofam_kos,
                    )
                ),
            }
        )

    comparisons.sort(
        key=lambda row: (
            row["sample_id"],
            row["feature_uid"],
        )
    )

    return comparisons
