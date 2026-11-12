"""Build report summaries from canonical BactAnnot tables"""

from pathlib import Path

import pandas as pd

COG_CATEGORY_DESCRIPTIONS = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}


EC_CLASS_DESCRIPTIONS = {
    "1": "Oxidoreductases",
    "2": "Transferases",
    "3": "Hydrolases",
    "4": "Lyases",
    "5": "Isomerases",
    "6": "Ligases",
    "7": "Translocases",
}


KEGG_ANNOTATION_TYPES = {
    "ko": "KEGG_KO",
    "pathway": "KEGG_PATHWAY",
    "module": "KEGG_MODULE",
    "reaction": "KEGG_REACTION",
    "brite": "KEGG_BRITE",
}


CAZY_CLASS_DESCRIPTIONS = {
    "GH": "Glycoside Hydrolases",
    "GT": "GlycosylTransferases",
    "PL": "Polysaccharide Lyases",
    "CE": "Carbohydrate Esterases",
    "AA": "Auxiliary Activities",
    "CBM": "Carbohydrate-Binding Modules",
}


TOOL_CONTRIBUTION_TYPES = (
    "COG",
    "GO",
    "EC",
    "KEGG_KO",
    "KEGG_PATHWAY",
    "PFAM",
    "INTERPRO",
    "SMART",
    "CAZY",
    "ARO",
    "VIRULENCE_ASSOCIATED",
    "TOXIN_ASSOCIATED",
)


# Carga las tablas canónicas de BactAnnot como DataFrames
def load_report_tables(
    features_path: Path,
    annotations_path: Path,
    evidence_path: Path,
    interpretations_path: Path,
    concordance_path: Path,
) -> dict[str, pd.DataFrame]:
    """Load canonical BactAnnot tables for reporting"""

    paths = {
        "features": features_path,
        "annotations": annotations_path,
        "evidence": evidence_path,
        "interpretations": interpretations_path,
        "concordance": concordance_path,
    }

    return {
        name: pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
        for name, path in paths.items()
    }


# Calcula las métricas principales directamente desde las tablas canónicas
def build_overview(
    features: pd.DataFrame,
    annotations: pd.DataFrame,
    interpretations: pd.DataFrame,
    concordance: pd.DataFrame,
) -> dict[str, object]:
    """Build overview metrics from canonical BactAnnot tables"""

    # Conservamos solo las features que representan proteínas
    protein_features = features[features["protein_id"].fillna("").ne("")]

    protein_feature_ids = set(protein_features["feature_uid"])

    # Las métricas de cobertura se calculan solo sobre proteínas
    protein_annotations = annotations[annotations["feature_uid"].isin(protein_feature_ids)]

    feature_count = len(features)

    protein_count = protein_features["feature_uid"].nunique()

    annotated_protein_count = protein_annotations["feature_uid"].nunique()

    # Cuenta proteínas distintas anotadas por cada herramienta
    coverage_by_source = (
        protein_annotations[["source", "feature_uid"]]
        .drop_duplicates()
        .groupby("source")
        .size()
        .sort_index()
    )

    protein_coverage_by_source = {
        str(source): int(count) for source, count in coverage_by_source.items()
    }

    # Cuenta proteínas distintas cubiertas por cada namespace
    coverage_by_annotation_type = (
        protein_annotations[["annotation_type", "feature_uid"]]
        .drop_duplicates()
        .groupby("annotation_type")
        .size()
        .sort_index()
    )

    protein_coverage_by_annotation_type = {
        str(annotation_type): int(count)
        for annotation_type, count in coverage_by_annotation_type.items()
    }

    # Resume las proteínas hipotéticas y su anotación externa
    hypothetical_features = protein_features[
        protein_features["hypothetical"].astype(str).eq("True")
    ]

    hypothetical_total = hypothetical_features["feature_uid"].nunique()

    external_feature_ids = set(
        protein_annotations.loc[
            protein_annotations["source"].ne("bakta"),
            "feature_uid",
        ]
    )

    externally_annotated = hypothetical_features[
        hypothetical_features["feature_uid"].isin(external_feature_ids)
    ]["feature_uid"].nunique()

    hypothetical_proteins = {
        "total": int(hypothetical_total),
        "externally_annotated": int(externally_annotated),
        "without_external_annotation": int(hypothetical_total - externally_annotated),
    }

    # Resume los estados de interpretación y el alcance de las propuestas principales
    interpretation_status_counts = {
        str(status): int(count)
        for status, count in (
            interpretations["interpretation_status"].value_counts().sort_index().items()
        )
    }

    primary_scope_counts = {
        str(scope): int(count)
        for scope, count in (interpretations["primary_scope"].value_counts().sort_index().items())
    }

    # Resume las relaciones de concordancia
    concordance_relation_counts = {
        str(relation): int(count)
        for relation, count in (concordance["relation"].value_counts().sort_index().items())
    }

    # Filtra una sola vez las anotaciones procedentes de PathoFact
    pathofact_annotations = protein_annotations[protein_annotations["source"].eq("PathoFact2")]

    amr_proteins = pathofact_annotations.loc[
        pathofact_annotations["annotation_type"].eq("ARO"),
        "feature_uid",
    ].nunique()

    virulence_high_confidence = pathofact_annotations.loc[
        pathofact_annotations["annotation_type"].eq("VIRULENCE_ASSOCIATED")
        & pathofact_annotations["annotation_value"].eq("high_confidence"),
        "feature_uid",
    ].nunique()

    virulence_ambiguous = pathofact_annotations.loc[
        pathofact_annotations["annotation_type"].eq("VIRULENCE_ASSOCIATED")
        & pathofact_annotations["annotation_value"].eq("ambiguous"),
        "feature_uid",
    ].nunique()

    toxin_high_confidence = pathofact_annotations.loc[
        pathofact_annotations["annotation_type"].eq("TOXIN_ASSOCIATED")
        & pathofact_annotations["annotation_value"].eq("high_confidence"),
        "feature_uid",
    ].nunique()

    toxin_ambiguous = pathofact_annotations.loc[
        pathofact_annotations["annotation_type"].eq("TOXIN_ASSOCIATED")
        & pathofact_annotations["annotation_value"].eq("ambiguous"),
        "feature_uid",
    ].nunique()

    pathofact = {
        "amr_proteins": int(amr_proteins),
        "virulence": {
            "high_confidence": int(virulence_high_confidence),
            "ambiguous": int(virulence_ambiguous),
        },
        "toxin": {
            "high_confidence": int(toxin_high_confidence),
            "ambiguous": int(toxin_ambiguous),
        },
    }

    return {
        "feature_count": int(feature_count),
        "protein_count": int(protein_count),
        "annotated_protein_count": int(annotated_protein_count),
        "protein_coverage_by_source": (protein_coverage_by_source),
        "protein_coverage_by_annotation_type": (protein_coverage_by_annotation_type),
        "hypothetical_proteins": hypothetical_proteins,
        "interpretation_status_counts": interpretation_status_counts,
        "primary_scope_counts": primary_scope_counts,
        "concordance_relation_counts": (concordance_relation_counts),
        "pathofact": pathofact,
    }


# Construye una tabla de una fila por proteína para el Feature Explorer
def build_feature_table(
    features: pd.DataFrame,
    annotations: pd.DataFrame,
) -> pd.DataFrame:
    """Build the protein-level Feature Explorer table"""

    # Conservamos solo las features codificantes
    protein_features = features[features["protein_id"].fillna("").ne("")].copy()

    # Agrupa IDs únicos por feature y tipo de anotación
    annotation_summary = (
        annotations[
            [
                "feature_uid",
                "annotation_type",
                "annotation_id",
            ]
        ]
        .drop_duplicates()
        .groupby(
            [
                "feature_uid",
                "annotation_type",
            ]
        )["annotation_id"]
        .agg(lambda values: "; ".join(sorted(values)))
        .unstack(fill_value="")
        .reset_index()
    )

    # Añade las anotaciones resumidas a la tabla estructural
    result = protein_features.merge(
        annotation_summary,
        on="feature_uid",
        how="left",
    )

    # Extrae estados interpretables de las anotaciones de PathoFact
    pathofact_annotations = annotations[annotations["source"].eq("PathoFact2")]

    amr_features = pathofact_annotations.loc[
        pathofact_annotations["annotation_type"].eq("ARO"),
        "feature_uid",
    ]

    result["amr"] = result["feature_uid"].isin(amr_features)

    virulence = (
        pathofact_annotations.loc[
            pathofact_annotations["annotation_type"].eq("VIRULENCE_ASSOCIATED"),
            [
                "feature_uid",
                "annotation_value",
            ],
        ]
        .drop_duplicates()
        .groupby("feature_uid")["annotation_value"]
        .agg(lambda values: "; ".join(sorted(value for value in values if value)))
    )

    result["virulence"] = result["feature_uid"].map(virulence).fillna("")

    toxin = (
        pathofact_annotations.loc[
            pathofact_annotations["annotation_type"].eq("TOXIN_ASSOCIATED"),
            [
                "feature_uid",
                "annotation_value",
            ],
        ]
        .drop_duplicates()
        .groupby("feature_uid")["annotation_value"]
        .agg(lambda values: "; ".join(sorted(value for value in values if value)))
    )

    result["toxin"] = result["feature_uid"].map(toxin).fillna("")

    # Las proteínas sin un namespace concreto deben mostrar cadena vacía
    annotation_columns = [
        column for column in annotation_summary.columns if column != "feature_uid"
    ]

    result[annotation_columns] = result[annotation_columns].fillna("")

    return result.reset_index(drop=True)


def build_protein_interpretation_table(
    features: pd.DataFrame,
    annotations: pd.DataFrame,
    interpretations: pd.DataFrame,
) -> pd.DataFrame:
    """Attach existing functional interpretations to their protein records."""

    table = build_feature_table(features, annotations)

    result = table.merge(
        interpretations,
        on=["sample_id", "feature_uid", "protein_id"],
        how="left",
        validate="one_to_one",
    )

    if len(interpretations) != len(table) or result["interpretation_status"].isna().any():
        raise ValueError("Protein interpretations must match protein features")

    return result


# Resume la composición funcional global mediante categorías COG
def build_functional_landscape(
    annotations: pd.DataFrame,
) -> dict[str, object]:
    """Build global functional summaries for the report"""

    # Conserva únicamente las categorías funcionales COG
    cog_categories = annotations.loc[
        annotations["annotation_type"].eq("COG_CATEGORY"),
        [
            "feature_uid",
            "annotation_id",
        ],
    ].drop_duplicates()

    # Cuenta proteínas únicas asociadas a cada categoría COG
    cog_category_counts = (
        cog_categories.groupby("annotation_id")["feature_uid"]
        .nunique()
        .sort_values(
            ascending=False,
        )
    )

    cog_category_summary = [
        {
            "category": str(category),
            "description": COG_CATEGORY_DESCRIPTIONS.get(
                str(category),
                "Unknown category",
            ),
            "protein_count": int(count),
        }
        for category, count in cog_category_counts.items()
    ]

    # Resume las anotaciones Gene Ontology sin duplicar soporte entre fuentes
    go_annotations = annotations.loc[
        annotations["annotation_type"].eq("GO"),
        [
            "feature_uid",
            "annotation_id",
        ],
    ].drop_duplicates()

    go_protein_count = go_annotations["feature_uid"].nunique()
    go_term_count = go_annotations["annotation_id"].nunique()
    go_association_count = len(go_annotations)

    # Cuenta cuántas proteínas distintas están asociadas a cada término GO
    go_term_counts = (
        go_annotations.groupby("annotation_id")["feature_uid"]
        .nunique()
        .reset_index(name="protein_count")
        .sort_values(
            ["protein_count", "annotation_id"],
            ascending=[False, True],
        )
        .reset_index(drop=True)
    )

    go_top_terms = [
        {
            "go_id": str(row.annotation_id),
            "protein_count": int(row.protein_count),
        }
        for row in go_term_counts.head(20).itertuples(index=False)
    ]

    # Resume las anotaciones enzimáticas EC sin duplicar soporte entre fuentes
    ec_annotations = annotations.loc[
        annotations["annotation_type"].eq("EC"),
        [
            "feature_uid",
            "annotation_id",
        ],
    ].drop_duplicates()

    ec_protein_count = ec_annotations["feature_uid"].nunique()
    ec_number_count = ec_annotations["annotation_id"].nunique()
    ec_association_count = len(ec_annotations)

    # Extrae la clase principal a partir del primer nivel del número EC
    ec_with_class = ec_annotations.copy()

    ec_with_class["ec_class"] = (
        ec_with_class["annotation_id"].str.removeprefix("EC:").str.split(".").str[0]
    )

    ec_class_counts = (
        ec_with_class[ec_with_class["ec_class"].isin(EC_CLASS_DESCRIPTIONS)]
        .drop_duplicates(
            [
                "feature_uid",
                "ec_class",
            ]
        )
        .groupby("ec_class")["feature_uid"]
        .nunique()
        .sort_index()
    )

    ec_class_summary = [
        {
            "ec_class": str(ec_class),
            "description": EC_CLASS_DESCRIPTIONS.get(
                str(ec_class),
                "Unknown class",
            ),
            "protein_count": int(count),
        }
        for ec_class, count in ec_class_counts.items()
    ]

    # Resume los distintos niveles funcionales asociados a KEGG
    kegg_summary = {}

    for key, annotation_type in KEGG_ANNOTATION_TYPES.items():
        kegg_annotations = annotations.loc[
            annotations["annotation_type"].eq(annotation_type),
            [
                "feature_uid",
                "annotation_id",
            ],
        ].drop_duplicates()

        kegg_entity_counts = (
            kegg_annotations.groupby("annotation_id")["feature_uid"]
            .nunique()
            .reset_index(name="protein_count")
            .sort_values(
                [
                    "protein_count",
                    "annotation_id",
                ],
                ascending=[
                    False,
                    True,
                ],
            )
            .reset_index(drop=True)
        )

        kegg_summary[key] = {
            "protein_count": int(kegg_annotations["feature_uid"].nunique()),
            "entity_count": int(kegg_annotations["annotation_id"].nunique()),
            "association_count": len(kegg_annotations),
            "top_entities": [
                {
                    "annotation_id": str(row.annotation_id),
                    "protein_count": int(row.protein_count),
                }
                for row in kegg_entity_counts.head(20).itertuples(index=False)
            ],
        }

    # Colapsa las representaciones map/ko de una misma ruta KEGG
    # únicamente para la capa de reporting
    pathway_annotations = (
        annotations.loc[
            annotations["annotation_type"].eq("KEGG_PATHWAY"),
            [
                "feature_uid",
                "annotation_id",
            ],
        ]
        .drop_duplicates()
        .copy()
    )

    pathway_codes = pathway_annotations["annotation_id"].str.extract(
        r"^KEGG:(?:map|ko)(\d{5})$",
        expand=False,
    )

    # Si algún identificador no sigue map/ko, se conserva sin modificar
    pathway_annotations["pathway_code"] = pathway_codes.fillna(pathway_annotations["annotation_id"])

    normalized_pathways = pathway_annotations[
        [
            "feature_uid",
            "pathway_code",
        ]
    ].drop_duplicates()

    pathway_counts = (
        normalized_pathways.groupby("pathway_code")["feature_uid"]
        .nunique()
        .reset_index(name="protein_count")
        .sort_values(
            [
                "protein_count",
                "pathway_code",
            ],
            ascending=[
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )

    kegg_summary["pathway"] = {
        "protein_count": int(normalized_pathways["feature_uid"].nunique()),
        "entity_count": int(normalized_pathways["pathway_code"].nunique()),
        "association_count": len(normalized_pathways),
        "top_entities": [
            {
                "pathway_code": str(row.pathway_code),
                "protein_count": int(row.protein_count),
            }
            for row in pathway_counts.head(20).itertuples(index=False)
        ],
    }

    # Resume las entidades InterPro a nivel global
    interpro_annotations = annotations.loc[
        annotations["annotation_type"].eq("INTERPRO"),
        [
            "feature_uid",
            "annotation_id",
            "description",
        ],
    ].drop_duplicates()

    interpro_pairs = interpro_annotations[
        [
            "feature_uid",
            "annotation_id",
        ]
    ].drop_duplicates()

    interpro_counts = (
        interpro_pairs.groupby("annotation_id")["feature_uid"]
        .nunique()
        .reset_index(name="protein_count")
        .sort_values(
            [
                "protein_count",
                "annotation_id",
            ],
            ascending=[
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )

    # InterProScan proporciona descripciones interpretables de las entradas InterPro
    interpro_descriptions = (
        interpro_annotations.loc[
            interpro_annotations["description"].ne(""),
            [
                "annotation_id",
                "description",
            ],
        ]
        .drop_duplicates()
        .groupby("annotation_id")["description"]
        .agg(lambda values: "; ".join(sorted(set(values))))
    )

    interpro_top_entities = [
        {
            "annotation_id": str(row.annotation_id),
            "description": str(
                interpro_descriptions.get(
                    row.annotation_id,
                    "",
                )
            ),
            "protein_count": int(row.protein_count),
        }
        for row in interpro_counts.head(20).itertuples(index=False)
    ]

    # Resume Pfam combinando el soporte de todas las fuentes
    pfam_annotations = annotations.loc[
        annotations["annotation_type"].eq("PFAM"),
        [
            "feature_uid",
            "source",
            "annotation_id",
            "description",
        ],
    ].drop_duplicates()

    pfam_pairs = pfam_annotations[
        [
            "feature_uid",
            "annotation_id",
        ]
    ].drop_duplicates()

    pfam_counts = (
        pfam_pairs.groupby("annotation_id")["feature_uid"]
        .nunique()
        .reset_index(name="protein_count")
        .sort_values(
            [
                "protein_count",
                "annotation_id",
            ],
            ascending=[
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )

    # Para el nombre de la familia Pfam se prioriza la descripción
    # procedente de InterProScan frente a productos específicos de otras fuentes
    pfam_descriptions = (
        pfam_annotations.loc[
            pfam_annotations["source"].eq("InterProScan") & pfam_annotations["description"].ne(""),
            [
                "annotation_id",
                "description",
            ],
        ]
        .drop_duplicates()
        .groupby("annotation_id")["description"]
        .agg(lambda values: "; ".join(sorted(set(values))))
    )

    pfam_top_entities = [
        {
            "annotation_id": str(row.annotation_id),
            "description": str(
                pfam_descriptions.get(
                    row.annotation_id,
                    "",
                )
            ),
            "protein_count": int(row.protein_count),
        }
        for row in pfam_counts.head(20).itertuples(index=False)
    ]

    domains = {
        "interpro": {
            "protein_count": int(interpro_pairs["feature_uid"].nunique()),
            "entity_count": int(interpro_pairs["annotation_id"].nunique()),
            "association_count": len(interpro_pairs),
            "top_entities": interpro_top_entities,
        },
        "pfam": {
            "protein_count": int(pfam_pairs["feature_uid"].nunique()),
            "entity_count": int(pfam_pairs["annotation_id"].nunique()),
            "association_count": len(pfam_pairs),
            "top_entities": pfam_top_entities,
        },
    }

    # Resume las anotaciones CAZy y sus clases funcionales principales
    cazy_annotations = (
        annotations.loc[
            annotations["annotation_type"].eq("CAZY"),
            [
                "feature_uid",
                "annotation_id",
            ],
        ]
        .drop_duplicates()
        .copy()
    )

    cazy_protein_count = cazy_annotations["feature_uid"].nunique()
    cazy_family_count = cazy_annotations["annotation_id"].nunique()
    cazy_association_count = len(cazy_annotations)

    # Extrae la clase CAZy a partir del prefijo de cada familia
    cazy_annotations["cazy_class"] = (
        cazy_annotations["annotation_id"]
        .str.extract(
            r"^([A-Za-z]+)",
            expand=False,
        )
        .str.upper()
    )

    cazy_class_counts = (
        cazy_annotations.loc[
            cazy_annotations["cazy_class"].isin(CAZY_CLASS_DESCRIPTIONS),
            [
                "feature_uid",
                "cazy_class",
            ],
        ]
        .drop_duplicates()
        .groupby("cazy_class")["feature_uid"]
        .nunique()
        .sort_values(
            ascending=False,
        )
    )

    cazy_class_summary = [
        {
            "cazy_class": str(cazy_class),
            "description": CAZY_CLASS_DESCRIPTIONS.get(
                str(cazy_class),
                "Unknown class",
            ),
            "protein_count": int(count),
        }
        for cazy_class, count in cazy_class_counts.items()
    ]

    cazy_family_counts = (
        cazy_annotations.groupby("annotation_id")["feature_uid"]
        .nunique()
        .reset_index(name="protein_count")
        .sort_values(
            [
                "protein_count",
                "annotation_id",
            ],
            ascending=[
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )

    cazy_top_families = [
        {
            "family": str(row.annotation_id),
            "protein_count": int(row.protein_count),
        }
        for row in cazy_family_counts.head(20).itertuples(index=False)
    ]

    cazy = {
        "protein_count": int(cazy_protein_count),
        "family_count": int(cazy_family_count),
        "association_count": int(cazy_association_count),
        "class_summary": cazy_class_summary,
        "top_families": cazy_top_families,
    }

    return {
        "cog_category_counts": {
            str(category): int(count) for category, count in cog_category_counts.items()
        },
        "cog_category_summary": cog_category_summary,
        "go": {
            "protein_count": int(go_protein_count),
            "term_count": int(go_term_count),
            "association_count": int(go_association_count),
            "top_terms": go_top_terms,
        },
        "ec": {
            "protein_count": int(ec_protein_count),
            "number_count": int(ec_number_count),
            "association_count": int(ec_association_count),
            "class_summary": ec_class_summary,
        },
        "kegg": kegg_summary,
        "domains": domains,
        "cazy": cazy,
    }


# Resume la contribución funcional de cada herramienta
def build_tool_contribution(
    annotations: pd.DataFrame,
) -> dict[str, object]:
    """Build source-by-annotation-type protein coverage summaries"""

    # Una proteína cuenta una sola vez por fuente y tipo de anotación
    source_type_features = annotations[
        [
            "source",
            "annotation_type",
            "feature_uid",
        ]
    ].drop_duplicates()

    source_type_counts = (
        source_type_features.groupby(
            [
                "source",
                "annotation_type",
            ]
        )["feature_uid"]
        .nunique()
        .reset_index(name="protein_count")
        .sort_values(
            [
                "source",
                "annotation_type",
            ]
        )
        .reset_index(drop=True)
    )

    contribution_records = [
        {
            "source": str(row.source),
            "annotation_type": str(row.annotation_type),
            "protein_count": int(row.protein_count),
        }
        for row in source_type_counts.itertuples(index=False)
    ]

    sources = sorted(str(source) for source in source_type_features["source"].unique())

    annotation_types = sorted(
        str(annotation_type) for annotation_type in source_type_features["annotation_type"].unique()
    )

    # Construye una matriz completa para los namespaces seleccionados
    # del heatmap, incluyendo explícitamente las combinaciones ausentes
    heatmap_index = pd.MultiIndex.from_product(
        [
            sources,
            TOOL_CONTRIBUTION_TYPES,
        ],
        names=[
            "source",
            "annotation_type",
        ],
    )

    heatmap_counts = (
        source_type_counts.set_index(
            [
                "source",
                "annotation_type",
            ]
        )["protein_count"]
        .reindex(
            heatmap_index,
            fill_value=0,
        )
        .reset_index()
    )

    heatmap_records = [
        {
            "source": str(row.source),
            "annotation_type": str(row.annotation_type),
            "protein_count": int(row.protein_count),
        }
        for row in heatmap_counts.itertuples(index=False)
    ]

    return {
        "sources": sources,
        "annotation_types": annotation_types,
        "records": contribution_records,
        "heatmap_annotation_types": list(TOOL_CONTRIBUTION_TYPES),
        "heatmap_records": heatmap_records,
    }


# Resume el soporte de las entidades funcionales del consenso
def build_annotation_support(
    consensus: pd.DataFrame,
) -> dict[str, object]:
    """Build entity-level annotation support summaries."""

    # Los TSV se cargan como texto, por lo que convertimos explícitamente
    # los campos numéricos del consenso
    source_counts = pd.to_numeric(
        consensus["source_count"],
        errors="raise",
    )

    method_counts = pd.to_numeric(
        consensus["method_count"],
        errors="raise",
    )

    support_status_counts = {
        str(status): int(count)
        for status, count in (consensus["support_status"].value_counts().sort_index().items())
    }

    source_count_distribution = {
        str(int(source_count)): int(count)
        for source_count, count in (source_counts.value_counts().sort_index().items())
    }

    method_count_distribution = {
        str(int(method_count)): int(count)
        for method_count, count in (method_counts.value_counts().sort_index().items())
    }

    return {
        "entity_count": (len(consensus)),
        "support_status_counts": support_status_counts,
        "source_count_distribution": source_count_distribution,
        "method_count_distribution": method_count_distribution,
        "single_source_entity_count": int(source_counts.eq(1).sum()),
        "multi_source_entity_count": int(source_counts.ge(2).sum()),
        "multimethod_entity_count": int(method_counts.ge(2).sum()),
    }


# Resume las comparaciones de concordancia entre herramientas
def build_concordance_report(
    concordance: pd.DataFrame,
) -> dict[str, object]:
    """Build pairwise concordance summaries for reporting"""

    relation_counts = (
        concordance.groupby(
            [
                "annotation_type",
                "source_1",
                "source_2",
                "relation",
            ]
        )
        .size()
        .unstack(
            fill_value=0,
        )
        .reset_index()
    )

    relation_columns = [
        "exact",
        "partial",
        "disjoint",
        "only_source_1",
        "only_source_2",
    ]

    for column in relation_columns:
        if column not in relation_counts.columns:
            relation_counts[column] = 0

    relation_counts["both"] = (
        relation_counts["exact"] + relation_counts["partial"] + relation_counts["disjoint"]
    )

    relation_counts["total"] = (
        relation_counts["both"]
        + relation_counts["only_source_1"]
        + relation_counts["only_source_2"]
    )

    # Los porcentajes exact/partial/disjoint se calculan solo
    # sobre features anotadas por ambas fuentes
    relation_counts["exact_percent"] = (
        relation_counts["exact"]
        .div(relation_counts["both"].where(relation_counts["both"].ne(0)))
        .mul(100)
        .fillna(0)
    )

    relation_counts["partial_percent"] = (
        relation_counts["partial"]
        .div(relation_counts["both"].where(relation_counts["both"].ne(0)))
        .mul(100)
        .fillna(0)
    )

    relation_counts["disjoint_percent"] = (
        relation_counts["disjoint"]
        .div(relation_counts["both"].where(relation_counts["both"].ne(0)))
        .mul(100)
        .fillna(0)
    )

    relation_counts = relation_counts.sort_values(
        [
            "annotation_type",
            "source_1",
            "source_2",
        ]
    ).reset_index(drop=True)

    comparisons = [
        {
            "annotation_type": str(row.annotation_type),
            "source_1": str(row.source_1),
            "source_2": str(row.source_2),
            "both": int(row.both),
            "exact": int(row.exact),
            "partial": int(row.partial),
            "disjoint": int(row.disjoint),
            "only_source_1": int(row.only_source_1),
            "only_source_2": int(row.only_source_2),
            "total": int(row.total),
            "exact_percent": float(row.exact_percent),
            "partial_percent": float(row.partial_percent),
            "disjoint_percent": float(row.disjoint_percent),
        }
        for row in relation_counts.itertuples(index=False)
    ]

    return {
        "comparison_count": len(comparisons),
        "comparisons": comparisons,
    }


# Agrupa las fuentes de evidencia por herramienta principal para el reporting
def _evidence_tool_name(source: str) -> str:
    """Return the main tool name for a canonical evidence source."""

    source = str(source)

    if source.startswith("InterProScan"):
        return "InterProScan"

    if source.startswith("PathoFact2"):
        return "PathoFact2"

    return source


# Resume la evidencia conservada por herramienta y tipo de evidencia
def build_evidence_context(
    evidence: pd.DataFrame,
) -> dict[str, object]:
    """Build source-level evidence summaries for the report."""

    data = evidence.copy()

    if data.empty:
        return {
            "evidence_count": 0,
            "protein_count": 0,
            "described_evidence_count": 0,
            "described_protein_count": 0,
            "tool_summary": [],
            "evidence_type_summary": [],
        }

    data["tool"] = data["source"].map(_evidence_tool_name)

    data["has_description"] = data["product"].fillna("").str.strip().ne("")

    tool_summary = (
        data.groupby("tool")
        .agg(
            protein_count=(
                "feature_uid",
                "nunique",
            ),
            evidence_count=(
                "feature_uid",
                "size",
            ),
            evidence_type_count=(
                "evidence_type",
                "nunique",
            ),
            described_evidence_count=(
                "has_description",
                "sum",
            ),
        )
        .reset_index()
        .sort_values(
            [
                "protein_count",
                "tool",
            ],
            ascending=[
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )

    tool_records = [
        {
            "tool": str(row.tool),
            "protein_count": int(row.protein_count),
            "evidence_count": int(row.evidence_count),
            "evidence_type_count": int(row.evidence_type_count),
            "described_evidence_count": int(row.described_evidence_count),
        }
        for row in tool_summary.itertuples(index=False)
    ]

    evidence_type_summary = (
        data.groupby(
            [
                "tool",
                "evidence_type",
            ]
        )
        .agg(
            protein_count=(
                "feature_uid",
                "nunique",
            ),
            evidence_count=(
                "feature_uid",
                "size",
            ),
        )
        .reset_index()
        .sort_values(
            [
                "tool",
                "protein_count",
                "evidence_type",
            ],
            ascending=[
                True,
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )

    evidence_type_records = [
        {
            "tool": str(row.tool),
            "evidence_type": str(row.evidence_type),
            "protein_count": int(row.protein_count),
            "evidence_count": int(row.evidence_count),
        }
        for row in evidence_type_summary.itertuples(index=False)
    ]

    described = data.loc[data["has_description"]]

    return {
        "evidence_count": len(data),
        "protein_count": int(data["feature_uid"].nunique()),
        "described_evidence_count": int(data["has_description"].sum()),
        "described_protein_count": int(described["feature_uid"].nunique()),
        "tool_summary": tool_records,
        "evidence_type_summary": (evidence_type_records),
    }


# Reúne evidencia contextual para proteínas hipotéticas de Bakta
def build_hypothetical_evidence(
    features: pd.DataFrame,
    evidence: pd.DataFrame,
) -> pd.DataFrame:
    """Build source-level evidence context for hypothetical proteins."""

    proteins = features.loc[features["protein_id"].fillna("").ne("")].copy()

    hypothetical = proteins.loc[proteins["hypothetical"].astype(str).eq("True")].copy()

    hypothetical_ids = set(hypothetical["feature_uid"])

    evidence_subset = evidence.loc[evidence["feature_uid"].isin(hypothetical_ids)].copy()

    columns = [
        "feature_uid",
        "protein_id",
        "gene",
        "bakta_product",
        "tool",
        "evidence_types",
        "accessions",
        "functional_descriptions",
        "evidence_count",
    ]

    if evidence_subset.empty:
        return pd.DataFrame(columns=columns)

    evidence_subset["tool"] = evidence_subset["source"].map(_evidence_tool_name)

    evidence_summary = (
        evidence_subset.groupby(
            [
                "feature_uid",
                "tool",
            ]
        )
        .agg(
            evidence_types=(
                "evidence_type",
                lambda values: "; ".join(sorted({str(value) for value in values if str(value)})),
            ),
            accessions=(
                "accession",
                lambda values: "; ".join(sorted({str(value) for value in values if str(value)})),
            ),
            functional_descriptions=(
                "product",
                lambda values: " | ".join(
                    sorted({str(value).strip() for value in values if str(value).strip()})
                ),
            ),
            evidence_count=(
                "feature_uid",
                "size",
            ),
        )
        .reset_index()
    )

    metadata = hypothetical[
        [
            "feature_uid",
            "protein_id",
            "gene",
            "product",
        ]
    ].rename(
        columns={
            "product": "bakta_product",
        }
    )

    result = metadata.merge(
        evidence_summary,
        on="feature_uid",
        how="inner",
    )

    result["evidence_count"] = result["evidence_count"].astype(int)

    return (
        result[columns]
        .sort_values(
            [
                "protein_id",
                "tool",
            ]
        )
        .reset_index(drop=True)
    )


def build_biological_highlights(
    features: pd.DataFrame,
    annotations: pd.DataFrame,
) -> dict:
    """Build biologically relevant summaries for the report"""

    proteins = features.loc[features["protein_id"].ne("")].copy()

    # Usa el campo canónico de Bakta en lugar de inferir el estado
    # a partir del texto de product.
    hypothetical = proteins.loc[proteins["hypothetical"].astype(str).eq("True")].copy()

    hypothetical_ids = set(hypothetical["feature_uid"])

    external_annotations = annotations.loc[
        annotations["feature_uid"].isin(hypothetical_ids) & annotations["source"].ne("bakta")
    ].copy()

    externally_annotated_ids = set(external_annotations["feature_uid"])

    source_summary = (
        external_annotations[["feature_uid", "source"]]
        .drop_duplicates()
        .groupby("source")
        .size()
        .sort_values(ascending=False)
        .rename("protein_count")
        .reset_index()
        .to_dict("records")
    )

    if external_annotations.empty:
        representative_hypothetical = []
    else:
        per_feature = (
            external_annotations.groupby("feature_uid")
            .agg(
                external_source_count=(
                    "source",
                    "nunique",
                ),
                annotation_type_count=(
                    "annotation_type",
                    "nunique",
                ),
                external_annotation_count=(
                    "annotation_id",
                    "nunique",
                ),
            )
            .reset_index()
        )

        source_labels = (
            external_annotations[["feature_uid", "source"]]
            .drop_duplicates()
            .groupby("feature_uid")["source"]
            .agg(lambda values: ";".join(sorted(values)))
            .rename("external_sources")
            .reset_index()
        )

        metadata_columns = [
            column
            for column in (
                "feature_uid",
                "protein_id",
                "gene",
                "product",
            )
            if column in hypothetical.columns
        ]

        representative_hypothetical = (
            hypothetical[metadata_columns]
            .merge(
                per_feature,
                on="feature_uid",
                how="inner",
            )
            .merge(
                source_labels,
                on="feature_uid",
                how="left",
            )
            .sort_values(
                [
                    "external_source_count",
                    "annotation_type_count",
                    "external_annotation_count",
                    "feature_uid",
                ],
                ascending=[
                    False,
                    False,
                    False,
                    True,
                ],
            )
            .head(10)
            .to_dict("records")
        )

    aro = annotations.loc[annotations["annotation_type"].eq("ARO")].copy()

    aro_unique = aro[
        [
            "feature_uid",
            "annotation_id",
            "description",
        ]
    ].drop_duplicates()

    aro_top = (
        aro_unique.groupby(
            [
                "annotation_id",
                "description",
            ],
            dropna=False,
        )["feature_uid"]
        .nunique()
        .sort_values(ascending=False)
        .rename("protein_count")
        .reset_index()
        .head(20)
        .to_dict("records")
    )

    def classification_summary(
        annotation_type: str,
    ) -> dict:
        subset = annotations.loc[annotations["annotation_type"].eq(annotation_type)].copy()

        if subset.empty:
            return {
                "protein_count": 0,
                "status_counts": [],
            }

        unique = subset[
            [
                "feature_uid",
                "annotation_value",
            ]
        ].drop_duplicates()

        status_counts = (
            unique.groupby("annotation_value")["feature_uid"]
            .nunique()
            .sort_values(ascending=False)
            .rename("protein_count")
            .reset_index()
            .rename(
                columns={
                    "annotation_value": "status",
                }
            )
            .to_dict("records")
        )

        return {
            "protein_count": int(unique["feature_uid"].nunique()),
            "status_counts": status_counts,
        }

    return {
        "hypothetical_proteins": {
            "total": len(hypothetical_ids),
            "with_external_annotation": len(externally_annotated_ids),
            "without_external_annotation": (len(hypothetical_ids) - len(externally_annotated_ids)),
            "source_summary": source_summary,
            "representative_features": (representative_hypothetical),
        },
        "amr": {
            "protein_count": int(aro["feature_uid"].nunique()),
            "aro_count": int(aro["annotation_id"].nunique()),
            "top_aro": aro_top,
        },
        "virulence": classification_summary("VIRULENCE_ASSOCIATED"),
        "toxin": classification_summary("TOXIN_ASSOCIATED"),
    }


def build_feature_explorer(
    features: pd.DataFrame,
    annotations: pd.DataFrame,
    interpretations: pd.DataFrame,
    concordance: pd.DataFrame,
) -> pd.DataFrame:
    """Build the protein-level table used by the Feature Explorer"""

    table = build_protein_interpretation_table(
        features,
        annotations,
        interpretations,
    ).copy()

    # Conserva el estado hipotético definido explícitamente por Bakta.
    table["hypothetical"] = table["hypothetical"].astype(str).eq("True")

    # Summarize annotation sources per protein.
    sources = (
        annotations[["feature_uid", "source"]]
        .drop_duplicates()
        .groupby("feature_uid")["source"]
        .agg(lambda values: ";".join(sorted(values)))
        .rename("annotation_sources")
        .reset_index()
    )

    table = table.merge(
        sources,
        on="feature_uid",
        how="left",
    )

    table["annotation_sources"] = table["annotation_sources"].fillna("")

    # Summarize concordance relations per protein.
    if concordance.empty:
        relation_counts = pd.DataFrame(columns=["feature_uid"])
    else:
        relation_counts = (
            concordance.groupby(
                [
                    "feature_uid",
                    "relation",
                ]
            )
            .size()
            .unstack(fill_value=0)
            .reset_index()
        )

    expected_relations = (
        "exact",
        "partial",
        "disjoint",
        "only_source_1",
        "only_source_2",
    )

    for relation in expected_relations:
        if relation not in relation_counts.columns:
            relation_counts[relation] = 0

    relation_counts = relation_counts.rename(
        columns={relation: (f"concordance_{relation}") for relation in expected_relations}
    )

    table = table.merge(
        relation_counts,
        on="feature_uid",
        how="left",
    )

    for relation in expected_relations:
        column = f"concordance_{relation}"

        table[column] = table[column].fillna(0).astype(int)

    table["has_disjoint_concordance"] = table["concordance_disjoint"].gt(0)

    return table


# Resume las métricas biológicas principales a nivel de proteína
def build_biological_overview(
    features: pd.DataFrame,
    annotations: pd.DataFrame,
) -> dict[str, int]:
    """Build biologically meaningful protein-level overview metrics."""

    proteins = features.loc[features["protein_id"].fillna("").ne("")].copy()
    protein_ids = set(proteins["feature_uid"])

    protein_annotations = annotations.loc[annotations["feature_uid"].isin(protein_ids)].copy()

    external_annotations = protein_annotations.loc[protein_annotations["source"].ne("bakta")]

    hypothetical = proteins.loc[proteins["hypothetical"].astype(str).eq("True")]
    hypothetical_ids = set(hypothetical["feature_uid"])

    external_hypothetical = external_annotations.loc[
        external_annotations["feature_uid"].isin(hypothetical_ids)
    ]

    domain_types = {
        "CDD",
        "CDD_PROFILE",
        "COILS",
        "FUNFAM",
        "GENE3D",
        "HAMAP",
        "INTERPRO",
        "KOG",
        "MOBIDBLITE",
        "NCBIFAM",
        "PANTHER",
        "PFAM",
        "PFAM_NAME",
        "PIRSF",
        "PRINTS",
        "PROSITE",
        "SFLD",
        "SMART",
        "SUPERFAMILY",
    }

    def protein_count(annotation_types: set[str]) -> int:
        return int(
            protein_annotations.loc[
                protein_annotations["annotation_type"].isin(annotation_types),
                "feature_uid",
            ].nunique()
        )

    return {
        "protein_count": int(proteins["feature_uid"].nunique()),
        "annotated_protein_count": int(protein_annotations["feature_uid"].nunique()),
        "externally_annotated_protein_count": int(external_annotations["feature_uid"].nunique()),
        "hypothetical_protein_count": int(hypothetical["feature_uid"].nunique()),
        "hypothetical_with_external_annotation_count": int(
            external_hypothetical["feature_uid"].nunique()
        ),
        "go_protein_count": protein_count({"GO"}),
        "enzyme_protein_count": protein_count({"EC"}),
        "kegg_ko_protein_count": protein_count({"KEGG_KO"}),
        "pathway_protein_count": protein_count({"KEGG_PATHWAY"}),
        "reaction_protein_count": protein_count({"KEGG_REACTION"}),
        "domain_family_protein_count": protein_count(domain_types),
        "amr_protein_count": protein_count({"ARO"}),
    }


# Reúne AMR, virulencia y toxinas en tablas orientadas a interpretación biológica
def build_biomedical_context(
    features: pd.DataFrame,
    annotations: pd.DataFrame,
    evidence: pd.DataFrame,
) -> dict[str, object]:
    """Build biologically interpretable AMR, virulence and toxin tables."""

    # ---------------------------------------------------------
    # 1. Prepara una tabla con una fila por proteína
    # ---------------------------------------------------------

    protein_table = build_feature_table(
        features,
        annotations,
    )

    protein_ids = set(protein_table["feature_uid"])

    # ---------------------------------------------------------
    # 2. Selecciona las anotaciones de resistencia antimicrobiana
    # ---------------------------------------------------------

    amr_annotations = annotations.loc[
        annotations["feature_uid"].isin(protein_ids) & annotations["annotation_type"].eq("ARO"),
        [
            "feature_uid",
            "annotation_id",
        ],
    ].copy()

    amr_annotations = amr_annotations.drop_duplicates().rename(
        columns={
            "annotation_id": "aro_id",
        }
    )

    # Añade el gen y el producto anotados por Bakta.
    protein_metadata = protein_table[
        [
            "feature_uid",
            "protein_id",
            "gene",
            "product",
        ]
    ].copy()

    protein_metadata = protein_metadata.rename(
        columns={
            "gene": "bakta_gene",
            "product": "bakta_product",
        }
    )

    amr = amr_annotations.merge(
        protein_metadata,
        on="feature_uid",
        how="left",
    )

    # ---------------------------------------------------------
    # 3. Selecciona las evidencias AMR de PathoFact2
    # ---------------------------------------------------------

    amr_evidence = evidence.loc[
        evidence["feature_uid"].isin(protein_ids)
        & evidence["evidence_type"].eq("amr_prediction")
        & evidence["accession"].fillna("").ne(""),
        [
            "feature_uid",
            "accession",
            "source",
            "gene",
            "product",
            "rank",
        ],
    ].copy()

    # PathoFact2:DeepARG -> DeepARG
    # PathoFact2:RGI     -> RGI
    amr_evidence["method"] = amr_evidence["source"].str.removeprefix("PathoFact2:")

    # Conserva el valor de "rank" junto al método que lo produjo.
    # No se interpreta como una categoría común entre herramientas.
    amr_evidence["source_rank"] = ""

    has_rank = amr_evidence["rank"].fillna("").ne("")

    amr_evidence.loc[
        has_rank,
        "source_rank",
    ] = (
        amr_evidence.loc[
            has_rank,
            "method",
        ]
        + ": "
        + amr_evidence.loc[
            has_rank,
            "rank",
        ]
    )

    # ---------------------------------------------------------
    # 4. Resume las evidencias de cada proteína y ARO
    # ---------------------------------------------------------
    #
    # Puede haber varias filas de evidencia para la misma combinación:
    #
    # feature_uid + ARO
    #
    # Aquí se reúnen en una sola fila para que el report sea legible.

    evidence_rows = []

    grouped_evidence = amr_evidence.groupby(
        [
            "feature_uid",
            "accession",
        ]
    )

    for (
        feature_uid,
        accession,
    ), group in grouped_evidence:
        # Métodos que detectaron este ARO.
        method_values = group["method"]
        method_values = method_values[method_values.ne("")]
        method_values = method_values.drop_duplicates().sort_values()
        methods = "; ".join(method_values)

        # Símbolos génicos propuestos por las herramientas.
        gene_values = group["gene"]
        gene_values = gene_values[gene_values.ne("")]
        gene_values = gene_values.drop_duplicates().sort_values()
        evidence_genes = "; ".join(gene_values)

        # Descripciones funcionales propuestas por las herramientas.
        product_values = group["product"]
        product_values = product_values[product_values.ne("")]
        product_values = product_values.drop_duplicates().sort_values()
        evidence_products = " | ".join(product_values)

        # Valores originales de rank, indicando qué método los produjo.
        rank_values = group["source_rank"]
        rank_values = rank_values[rank_values.ne("")]
        rank_values = rank_values.drop_duplicates().sort_values()
        source_ranks = "; ".join(rank_values)

        evidence_rows.append(
            {
                "feature_uid": feature_uid,
                "aro_id": accession,
                "methods": methods,
                "evidence_genes": (evidence_genes),
                "evidence_products": (evidence_products),
                "source_ranks": (source_ranks),
            }
        )

    evidence_summary = pd.DataFrame(
        evidence_rows,
        columns=[
            "feature_uid",
            "aro_id",
            "methods",
            "evidence_genes",
            "evidence_products",
            "source_ranks",
        ],
    )

    # ---------------------------------------------------------
    # 5. Une anotaciones ARO y evidencias AMR
    # ---------------------------------------------------------

    amr = amr.merge(
        evidence_summary,
        on=[
            "feature_uid",
            "aro_id",
        ],
        how="left",
    )

    amr = amr.fillna("")

    # Cuenta cuántos métodos distintos apoyan cada ARO.
    amr["method_count"] = 0

    has_methods = amr["methods"].ne("")

    amr.loc[
        has_methods,
        "method_count",
    ] = (
        amr.loc[
            has_methods,
            "methods",
        ].str.count("; ")
        + 1
    )

    amr["method_count"] = amr["method_count"].astype(int)

    # Coloca primero las anotaciones observadas por más métodos.
    amr = amr.sort_values(
        [
            "method_count",
            "protein_id",
            "aro_id",
        ],
        ascending=[
            False,
            True,
            True,
        ],
    ).reset_index(drop=True)

    # ---------------------------------------------------------
    # 6. Construye la tabla de virulencia
    # ---------------------------------------------------------

    virulence_columns = [
        "feature_uid",
        "protein_id",
        "gene",
        "product",
        "virulence",
    ]

    for column in (
        "KEGG_KO",
        "KEGG_PATHWAY",
        "INTERPRO",
        "PFAM",
    ):
        if column in protein_table.columns:
            virulence_columns.append(column)

    virulence = protein_table.loc[
        protein_table["virulence"].fillna("").ne(""),
        virulence_columns,
    ].copy()

    virulence_order = pd.Categorical(
        virulence["virulence"],
        categories=[
            "high_confidence",
            "ambiguous",
        ],
        ordered=True,
    )

    virulence["_status_order"] = virulence_order

    virulence = (
        virulence.sort_values(
            [
                "_status_order",
                "protein_id",
            ]
        )
        .drop(columns="_status_order")
        .reset_index(drop=True)
    )

    # ---------------------------------------------------------
    # 7. Construye la tabla de toxinas
    # ---------------------------------------------------------

    toxin_columns = [
        "feature_uid",
        "protein_id",
        "gene",
        "product",
        "toxin",
    ]

    for column in (
        "KEGG_KO",
        "KEGG_PATHWAY",
        "INTERPRO",
        "PFAM",
    ):
        if column in protein_table.columns:
            toxin_columns.append(column)

    toxin = protein_table.loc[
        protein_table["toxin"].fillna("").ne(""),
        toxin_columns,
    ].copy()

    toxin_order = pd.Categorical(
        toxin["toxin"],
        categories=[
            "high_confidence",
            "ambiguous",
        ],
        ordered=True,
    )

    toxin["_status_order"] = toxin_order

    toxin = (
        toxin.sort_values(
            [
                "_status_order",
                "protein_id",
            ]
        )
        .drop(columns="_status_order")
        .reset_index(drop=True)
    )

    # ---------------------------------------------------------
    # 8. Resume las principales métricas biomédicas
    # ---------------------------------------------------------

    amr_protein_count = int(amr["feature_uid"].nunique())

    aro_count = int(amr["aro_id"].nunique())

    multimethod_amr = amr.loc[amr["method_count"].ge(2)]

    amr_multimethod_protein_count = int(multimethod_amr["feature_uid"].nunique())

    virulence_protein_count = int(virulence["feature_uid"].nunique())

    virulence_high_confidence = virulence.loc[
        virulence["virulence"].eq("high_confidence"),
        "feature_uid",
    ]

    virulence_high_confidence_count = int(virulence_high_confidence.nunique())

    toxin_protein_count = int(toxin["feature_uid"].nunique())

    toxin_high_confidence = toxin.loc[
        toxin["toxin"].eq("high_confidence"),
        "feature_uid",
    ]

    toxin_high_confidence_count = int(toxin_high_confidence.nunique())

    # ---------------------------------------------------------
    # 9. Devuelve las tablas y el resumen
    # ---------------------------------------------------------

    return {
        "summary": {
            "amr_protein_count": (amr_protein_count),
            "aro_count": (aro_count),
            "amr_multimethod_protein_count": (amr_multimethod_protein_count),
            "virulence_protein_count": (virulence_protein_count),
            "virulence_high_confidence_count": (virulence_high_confidence_count),
            "toxin_protein_count": (toxin_protein_count),
            "toxin_high_confidence_count": (toxin_high_confidence_count),
        },
        "amr": amr,
        "virulence": virulence,
        "toxin": toxin,
    }


# Combina información funcional complementaria para proteínas hipotéticas de Bakta
def build_hypothetical_characterization(
    features: pd.DataFrame,
    annotations: pd.DataFrame,
    evidence: pd.DataFrame,
    consensus: pd.DataFrame,
) -> pd.DataFrame:
    """Build an integrated biological view of Bakta hypothetical proteins."""

    table = build_feature_table(features, annotations)
    table = table.loc[table["hypothetical"].astype(str).eq("True")].copy()

    hypothetical_ids = set(table["feature_uid"])

    external = annotations.loc[
        annotations["feature_uid"].isin(hypothetical_ids) & annotations["source"].ne("bakta")
    ].copy()

    sources = (
        external[["feature_uid", "source"]]
        .drop_duplicates()
        .groupby("feature_uid")["source"]
        .agg(lambda x: "; ".join(sorted(set(x))))
    )

    table["external_sources"] = table["feature_uid"].map(sources).fillna("")
    table["has_external_annotation"] = table["external_sources"].ne("")

    kofam = (
        annotations.loc[
            annotations["feature_uid"].isin(hypothetical_ids)
            & annotations["source"].eq("KofamScan")
            & annotations["annotation_type"].eq("KEGG_KO")
            & annotations["description"].fillna("").ne(""),
            ["feature_uid", "description"],
        ]
        .drop_duplicates()
        .groupby("feature_uid")["description"]
        .agg(lambda x: " | ".join(sorted(set(x))))
    )

    table["kofam_descriptions"] = table["feature_uid"].map(kofam).fillna("")

    for annotation_type, column in (
        ("INTERPRO", "interpro_descriptions"),
        ("PFAM", "pfam_descriptions"),
    ):
        descriptions = (
            annotations.loc[
                annotations["feature_uid"].isin(hypothetical_ids)
                & annotations["source"].eq("InterProScan")
                & annotations["annotation_type"].eq(annotation_type)
                & annotations["description"].fillna("").ne(""),
                ["feature_uid", "description"],
            ]
            .drop_duplicates()
            .groupby("feature_uid")["description"]
            .agg(lambda x: " | ".join(sorted(set(x))))
        )

        table[column] = table["feature_uid"].map(descriptions).fillna("")

    eggnog = (
        evidence.loc[
            evidence["feature_uid"].isin(hypothetical_ids)
            & evidence["source"].eq("eggNOG-mapper")
            & evidence["product"].fillna("").ne(""),
            ["feature_uid", "product"],
        ]
        .drop_duplicates()
        .groupby("feature_uid")["product"]
        .agg(lambda x: " | ".join(sorted(set(x))))
    )

    table["eggnog_descriptions"] = table["feature_uid"].map(eggnog).fillna("")

    support = consensus.loc[consensus["feature_uid"].isin(hypothetical_ids)].copy()

    support["source_count_numeric"] = pd.to_numeric(
        support["source_count"],
        errors="coerce",
    ).fillna(0)

    multi_source = support.loc[support["source_count_numeric"].ge(2)].groupby("feature_uid").size()

    table["multi_source_entity_count"] = (
        table["feature_uid"].map(multi_source).fillna(0).astype(int)
    )

    biological_columns = [
        column
        for column in (
            "KEGG_KO",
            "EC",
            "KEGG_PATHWAY",
            "KEGG_REACTION",
            "GO",
            "INTERPRO",
            "PFAM",
            "CAZY",
        )
        if column in table.columns
    ]

    table["functional_dimension_count"] = (
        table[biological_columns].fillna("").ne("").sum(axis=1).astype(int)
    )

    return table.sort_values(
        ["functional_dimension_count", "multi_source_entity_count", "protein_id"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
