"""Pairwise annotation concordance diagnostics"""

from __future__ import annotations

import csv
from collections import defaultdict
from collections.abc import Iterable
from itertools import combinations
from pathlib import Path

from bactannot.integration.annotations import read_annotation_table
from bactannot.schema import CONCORDANCE_COLUMNS, CONCORDANCE_SUMMARY_COLUMNS

# Representa una comparación entre source_1, source_2 y un tipo de anotación
Comparison = tuple[str, str, str]


# Limita la concordancia a tipos de anotación con identificadores comparables entre herramientas
COMPARABLE_ANNOTATION_TYPES = frozenset(
    {
        "CDD",
        "COG",
        "EC",
        "GO",
        "KEGG_KO",
        "NCBIFAM",
        "PFAM",
        "SMART",
    }
)


# Define las columnas mínimas de features.tsv necesarias para validar las anotaciones
REQUIRED_FEATURE_COLUMNS = {
    "sample_id",
    "feature_uid",
    "protein_id",
}


# Carga la relación feature_uid -> protein_id definida por Bakta para esta muestra
def load_feature_to_protein_index(
    features_path: Path,
    sample_id: str,
) -> dict[str, str]:
    """Map canonical feature UIDs to their protein identifiers"""

    feature_index: dict[str, str] = {}

    with features_path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(f"Missing header in features table: {features_path}")

        missing = REQUIRED_FEATURE_COLUMNS.difference(reader.fieldnames)

        if missing:
            raise ValueError(
                f"Missing required columns in {features_path}: {', '.join(sorted(missing))}"
            )

        for row in reader:
            row_sample_id = row["sample_id"].strip()
            feature_uid = row["feature_uid"].strip()
            protein_id = row["protein_id"].strip()

            # Todas las features deben pertenecer a la muestra que estamos analizando
            if row_sample_id != sample_id:
                raise ValueError(
                    f"Unexpected sample_id in {features_path}: {row_sample_id} != {sample_id}"
                )

            if not feature_uid:
                raise ValueError(f"Missing feature_uid in {features_path}")

            # Cada feature_uid debe identificar una única feature
            if feature_uid in feature_index:
                raise ValueError(f"Duplicate feature_uid in features table: {feature_uid}")

            feature_index[feature_uid] = protein_id

    return feature_index


# Detecta automáticamente qué pares de herramientas pueden compararse en esta ejecución
def derive_evaluable_comparisons(
    annotation_rows: list[dict[str, str]],
) -> tuple[Comparison, ...]:
    """Derive protein-level comparisons evaluable in the current run"""

    # Guarda qué herramientas han producido cada tipo comparable de anotación
    sources_by_type: dict[
        str,
        set[str],
    ] = defaultdict(set)

    for row in annotation_rows:
        # Concordance trabaja solo con features que producen proteínas
        if not row["protein_id"].strip():
            continue

        annotation_type = row["annotation_type"]

        # Ignora tipos que no hemos definido como directamente comparables
        if annotation_type not in COMPARABLE_ANNOTATION_TYPES:
            continue

        sources_by_type[annotation_type].add(row["source"])

    output: list[Comparison] = []

    # Para cada tipo genera todos los pares posibles entre las fuentes que realmente lo anotaron
    for annotation_type in sorted(
        sources_by_type,
        key=lambda value: (
            value.casefold(),
            value,
        ),
    ):
        sources = sorted(
            sources_by_type[annotation_type],
            key=lambda value: (
                value.casefold(),
                value,
            ),
        )

        for source_1, source_2 in combinations(
            sources,
            2,
        ):
            output.append(
                (
                    source_1,
                    source_2,
                    annotation_type,
                )
            )

    return tuple(output)


# Convierte un conjunto de identificadores en una cadena ordenada y reproducible
def serialize_ids(
    values: set[str],
) -> str:
    """Serialize annotation identifiers deterministically"""

    return ";".join(
        sorted(
            values,
            key=lambda value: (
                value.casefold(),
                value,
            ),
        )
    )


# Clasifica la relación entre los conjuntos de anotaciones aportados por dos fuentes
def resolve_relation(
    source_1_ids: set[str],
    source_2_ids: set[str],
) -> str | None:
    """Classify the relationship between two annotation sets"""

    if source_1_ids and source_2_ids:
        # Las dos fuentes aportan exactamente el mismo conjunto de identificadores
        if source_1_ids == source_2_ids:
            return "exact"

        # Comparten al menos un identificador pero cada conjunto no es idéntico
        if source_1_ids & source_2_ids:
            return "partial"

        # Ambas fuentes anotan este tipo pero no comparten ningún identificador
        return "disjoint"

    # Solo la primera fuente aporta identificadores para esta feature
    if source_1_ids:
        return "only_source_1"

    # Solo la segunda fuente aporta identificadores para esta feature
    if source_2_ids:
        return "only_source_2"

    # Si ninguna aporta anotaciones no hay una comparación que registrar
    return None


# Compara por pares los conjuntos de identificadores funcionales de cada feature
def build_concordance(
    annotation_rows: list[dict[str, str]],
    *,
    comparisons: Iterable[Comparison] | None = None,
) -> list[dict[str, str]]:
    """Build pairwise annotation concordance diagnostics"""

    # Si no se indica una matriz concreta, deriva las comparaciones posibles de los propios datos
    if comparisons is None:
        selected_comparisons = derive_evaluable_comparisons(annotation_rows)
    else:
        selected_comparisons = tuple(dict.fromkeys(comparisons))

    # Conserva solo las combinaciones source + tipo que participan en alguna comparación
    relevant_source_types = {
        (source, annotation_type)
        for source_1, source_2, annotation_type in selected_comparisons
        for source in (
            source_1,
            source_2,
        )
    }

    # Permite comprobar que cada feature mantiene un único protein_id
    protein_ids: dict[
        tuple[str, str],
        set[str],
    ] = defaultdict(set)

    # Reúne todos los annotation_id de cada feature, tipo y fuente como un conjunto
    annotation_sets: dict[
        tuple[str, str, str, str],
        set[str],
    ] = defaultdict(set)

    for row in annotation_rows:
        sample_id = row["sample_id"]
        feature_uid = row["feature_uid"]
        protein_id = row["protein_id"]
        source = row["source"]
        annotation_type = row["annotation_type"]
        annotation_id = row["annotation_id"]

        feature_key = (
            sample_id,
            feature_uid,
        )

        protein_ids[feature_key].add(protein_id)

        # Las anotaciones que no participan en ninguna comparación no necesitan indexarse
        if (
            source,
            annotation_type,
        ) not in relevant_source_types:
            continue

        annotation_sets[
            (
                sample_id,
                feature_uid,
                annotation_type,
                source,
            )
        ].add(annotation_id)

    # Convierte el conjunto de protein_id de cada feature en un único valor ya validado
    resolved_protein_ids: dict[
        tuple[str, str],
        str,
    ] = {}

    for feature_key, values in protein_ids.items():
        if len(values) != 1:
            raise ValueError(
                f"Conflicting protein_id values for feature {feature_key}: {sorted(values)}"
            )

        resolved_protein_ids[feature_key] = next(iter(values))

    output_rows: list[dict[str, str]] = []

    # Evalúa cada pareja de fuentes por separado para cada tipo de anotación
    for (
        source_1,
        source_2,
        annotation_type,
    ) in selected_comparisons:
        # Selecciona las features donde al menos una de las dos fuentes aporta ese tipo
        feature_keys = {
            (
                sample_id,
                feature_uid,
            )
            for (
                sample_id,
                feature_uid,
                current_type,
                source,
            ) in annotation_sets
            if current_type == annotation_type
            and source
            in {
                source_1,
                source_2,
            }
        }

        for (
            sample_id,
            feature_uid,
        ) in feature_keys:
            protein_id = resolved_protein_ids[
                (
                    sample_id,
                    feature_uid,
                )
            ]

            # Concordance se limita a las features proteicas
            if not protein_id:
                continue

            # Recupera el conjunto completo de IDs aportado por la primera fuente
            source_1_ids = annotation_sets.get(
                (
                    sample_id,
                    feature_uid,
                    annotation_type,
                    source_1,
                ),
                set(),
            )

            # Recupera el conjunto completo de IDs aportado por la segunda fuente
            source_2_ids = annotation_sets.get(
                (
                    sample_id,
                    feature_uid,
                    annotation_type,
                    source_2,
                ),
                set(),
            )

            relation = resolve_relation(
                source_1_ids,
                source_2_ids,
            )

            if relation is None:
                continue

            # La intersección muestra exactamente qué identificadores comparten ambas fuentes
            shared_ids = source_1_ids & source_2_ids

            output_rows.append(
                {
                    "sample_id": sample_id,
                    "feature_uid": feature_uid,
                    "protein_id": protein_id,
                    "annotation_type": annotation_type,
                    "source_1": source_1,
                    "source_2": source_2,
                    "source_1_ids": serialize_ids(source_1_ids),
                    "source_2_ids": serialize_ids(source_2_ids),
                    "shared_ids": serialize_ids(shared_ids),
                    "relation": relation,
                }
            )

    # Ordena las comparaciones para obtener una salida estable y reproducible
    return sorted(
        output_rows,
        key=lambda row: (
            row["sample_id"],
            row["feature_uid"],
            row["annotation_type"],
            row["source_1"],
            row["source_2"],
        ),
    )


# Define todas las relaciones posibles que puede producir la comparación por pares
CONCORDANCE_RELATIONS = (
    "exact",
    "partial",
    "disjoint",
    "only_source_1",
    "only_source_2",
)


# Resume cuántas features caen en cada relación para cada par de fuentes y tipo de anotación
def build_concordance_summary(
    concordance_rows: list[dict[str, str]],
) -> list[dict[str, str | int]]:
    """Summarize pairwise concordance diagnostics"""

    # Inicializa a cero todos los tipos de relación para cada grupo
    grouped_counts: dict[
        tuple[str, str, str, str],
        dict[str, int],
    ] = defaultdict(lambda: {relation: 0 for relation in CONCORDANCE_RELATIONS})

    for row in concordance_rows:
        relation = row["relation"]

        if relation not in CONCORDANCE_RELATIONS:
            raise ValueError(f"Unknown concordance relation: {relation}")

        # Cada resumen corresponde a una muestra, tipo y pareja concreta de fuentes
        group_key = (
            row["sample_id"],
            row["annotation_type"],
            row["source_1"],
            row["source_2"],
        )

        grouped_counts[group_key][relation] += 1

    output_rows: list[dict[str, str | int]] = []

    for (
        sample_id,
        annotation_type,
        source_1,
        source_2,
    ), counts in grouped_counts.items():
        exact = counts["exact"]
        partial = counts["partial"]
        disjoint = counts["disjoint"]
        only_source_1 = counts["only_source_1"]
        only_source_2 = counts["only_source_2"]

        # both cuenta las features donde las dos fuentes aportaron alguna anotación
        both = exact + partial + disjoint

        # total incluye también las features anotadas únicamente por una de las dos fuentes
        total = both + only_source_1 + only_source_2

        output_rows.append(
            {
                "sample_id": sample_id,
                "annotation_type": annotation_type,
                "source_1": source_1,
                "source_2": source_2,
                "both": both,
                "exact": exact,
                "partial": partial,
                "disjoint": disjoint,
                "only_source_1": only_source_1,
                "only_source_2": only_source_2,
                "total": total,
            }
        )

    # Mantiene siempre el mismo orden en la tabla resumen
    return sorted(
        output_rows,
        key=lambda row: (
            row["sample_id"],
            row["annotation_type"],
            row["source_1"],
            row["source_2"],
        ),
    )


# Escribe la tabla detallada con una comparación por feature, tipo y pareja de fuentes
def write_concordance_table(
    rows: list[dict[str, str]],
    output_path: Path,
) -> None:
    """Write the canonical concordance table"""

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CONCORDANCE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


# Escribe la tabla resumen con los recuentos de cada tipo de relación
def write_concordance_summary_table(
    rows: list[dict[str, str | int]],
    output_path: Path,
) -> None:
    """Write the canonical concordance summary table"""

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CONCORDANCE_SUMMARY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


# Lee las anotaciones integradas, construye la concordancia y escribe ambas tablas
def build_concordance_tables(
    annotations_path: Path,
    features_path: Path,
    sample_id: str,
    concordance_output_path: Path,
    summary_output_path: Path,
) -> None:
    """Build and write canonical concordance tables"""

    # features.tsv sigue siendo la referencia para validar feature_uid y protein_id
    feature_index = load_feature_to_protein_index(
        features_path,
        sample_id,
    )

    # Reutiliza las mismas validaciones de annotations que usa la integración
    annotation_rows = read_annotation_table(
        annotations_path,
        sample_id,
        feature_index,
    )

    # Construye primero la comparación detallada por feature
    concordance_rows = build_concordance(
        annotation_rows,
    )

    # Resume después esas relaciones para cada pareja de fuentes
    summary_rows = build_concordance_summary(
        concordance_rows,
    )

    write_concordance_table(
        concordance_rows,
        concordance_output_path,
    )

    write_concordance_summary_table(
        summary_rows,
        summary_output_path,
    )
