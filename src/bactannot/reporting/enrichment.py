"""Optional biological enrichment for BactAnnot reports.

This module adds current biological names and ontology context from
external resources. It does not modify BactAnnot canonical output tables.
"""

import urllib.request

import pandas as pd


# Enriches KEGG identifiers with current names from the KEGG REST API
def build_kegg_enrichment(
    annotations: pd.DataFrame,
) -> pd.DataFrame:
    """Resolve KEGG identifiers to current KEGG names for reporting."""

    supported_types = {
        "KEGG_KO",
        "KEGG_PATHWAY",
        "KEGG_MODULE",
        "KEGG_REACTION",
    }

    kegg_annotations = annotations.loc[
        annotations["annotation_type"].isin(supported_types)
        & annotations["annotation_id"].fillna("").ne(""),
        [
            "annotation_type",
            "annotation_id",
        ],
    ].drop_duplicates()

    output_columns = [
        "annotation_type",
        "annotation_id",
        "kegg_database",
        "name",
        "resolved",
    ]

    if kegg_annotations.empty:
        return pd.DataFrame(columns=output_columns)

    database_names = {
        "KEGG_KO": "ko",
        "KEGG_PATHWAY": "pathway",
        "KEGG_MODULE": "module",
        "KEGG_REACTION": "reaction",
    }

    names_by_database = {}

    for database in (
        "ko",
        "pathway",
        "module",
        "reaction",
    ):
        url = "https://rest.kegg.jp/list/" + database

        request = urllib.request.Request(
            url,
            headers={
                "Accept": "text/plain",
                "User-Agent": "BactAnnot/0.1",
            },
        )

        with urllib.request.urlopen(
            request,
            timeout=60,
        ) as response:
            response_text = response.read().decode("utf-8")

        database_names_lookup = {}

        for line in response_text.splitlines():
            if not line.strip():
                continue

            fields = line.split(
                "\t",
                maxsplit=1,
            )

            if len(fields) != 2:
                continue

            external_id = fields[0]
            name = fields[1]

            if ":" in external_id:
                external_id = external_id.split(":", maxsplit=1)[1]

            database_names_lookup[external_id] = name

        names_by_database[database] = database_names_lookup

    rows = []

    for row in kegg_annotations.itertuples(index=False):
        annotation_type = str(row.annotation_type)
        annotation_id = str(row.annotation_id)

        kegg_database = database_names[annotation_type]

        lookup_id = annotation_id.removeprefix("KEGG:")

        # KEGG pathway lists use reference map IDs.
        # BactAnnot may contain either mapXXXXX or koXXXXX.
        if annotation_type == "KEGG_PATHWAY" and lookup_id.startswith("ko"):
            lookup_id = "map" + lookup_id[2:]

        name = names_by_database[kegg_database].get(
            lookup_id,
            "",
        )

        rows.append(
            {
                "annotation_type": annotation_type,
                "annotation_id": annotation_id,
                "kegg_database": kegg_database,
                "name": name,
                "resolved": bool(name),
            }
        )

    result = pd.DataFrame(
        rows,
        columns=output_columns,
    )

    return result.sort_values(
        [
            "annotation_type",
            "annotation_id",
        ]
    ).reset_index(drop=True)


# Enriches GO identifiers with the current Gene Ontology
# and maps them to the prokaryote GO slim for broad reporting
def build_go_enrichment(
    annotations: pd.DataFrame,
) -> pd.DataFrame:
    """Resolve GO identifiers and broad prokaryote GO slim categories."""

    go_ids = (
        annotations.loc[
            annotations["annotation_type"].eq("GO")
            & annotations["annotation_id"].fillna("").ne(""),
            "annotation_id",
        ]
        .drop_duplicates()
        .sort_values()
        .tolist()
    )

    output_columns = [
        "go_id",
        "canonical_go_id",
        "name",
        "aspect",
        "definition",
        "status",
        "original_is_obsolete",
        "replacement_go_id",
        "consider_go_ids",
        "prokaryote_slim_go_ids",
        "prokaryote_slim_names",
        "resolved",
    ]

    if not go_ids:
        return pd.DataFrame(columns=output_columns)

    # The current GO basic ontology is acyclic and is recommended
    # by the Gene Ontology Consortium for most annotation tools.
    url = "https://current.geneontology.org/ontology/go-basic.obo"

    request = urllib.request.Request(
        url,
        headers={
            "Accept": "text/plain",
            "User-Agent": "BactAnnot/0.1",
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=60,
    ) as response:
        ontology_text = response.read().decode("utf-8")

    # ---------------------------------------------------------
    # 1. Parse the GO OBO file
    # ---------------------------------------------------------
    #
    # Only the fields needed by the report are retained.

    terms = {}
    current_term = None

    for raw_line in ontology_text.splitlines():
        line = raw_line.strip()

        if line == "[Term]":
            if current_term is not None and current_term["id"]:
                terms[current_term["id"]] = current_term

            current_term = {
                "id": "",
                "name": "",
                "namespace": "",
                "definition": "",
                "alt_ids": [],
                "parents": set(),
                "subsets": set(),
                "is_obsolete": False,
                "replaced_by": [],
                "consider": [],
            }
            continue

        # A new stanza means the previous [Term] block is complete.
        if line.startswith("["):
            if current_term is not None and current_term["id"]:
                terms[current_term["id"]] = current_term

            current_term = None
            continue

        if current_term is None:
            continue

        if line.startswith("id: "):
            current_term["id"] = line.removeprefix("id: ")

        elif line.startswith("name: "):
            current_term["name"] = line.removeprefix("name: ")

        elif line.startswith("namespace: "):
            current_term["namespace"] = line.removeprefix("namespace: ")

        elif line.startswith("def: "):
            definition_text = line.removeprefix("def: ")

            definition_parts = definition_text.split(
                '"',
                maxsplit=2,
            )

            if len(definition_parts) >= 2:
                current_term["definition"] = definition_parts[1]

        elif line.startswith("alt_id: "):
            alternate_id = line.removeprefix("alt_id: ")
            current_term["alt_ids"].append(alternate_id)

        elif line.startswith("is_a: "):
            parent_id = line.removeprefix("is_a: ").split()[0]
            current_term["parents"].add(parent_id)

        elif line.startswith("relationship: part_of "):
            parent_id = line.removeprefix("relationship: part_of ").split()[0]
            current_term["parents"].add(parent_id)

        elif line.startswith("subset: "):
            subset_name = line.removeprefix("subset: ")
            current_term["subsets"].add(subset_name)

        elif line == "is_obsolete: true":
            current_term["is_obsolete"] = True

        elif line.startswith("replaced_by: "):
            replacement_id = line.removeprefix("replaced_by: ")
            current_term["replaced_by"].append(replacement_id)

        elif line.startswith("consider: "):
            consider_id = line.removeprefix("consider: ")
            current_term["consider"].append(consider_id)

    # Save the final term in the file.
    if current_term is not None and current_term["id"]:
        terms[current_term["id"]] = current_term

    # ---------------------------------------------------------
    # 2. Build a lookup for alternate GO identifiers
    # ---------------------------------------------------------

    alternate_ids = {}

    for primary_id, term in terms.items():
        for alternate_id in term["alt_ids"]:
            alternate_ids[alternate_id] = primary_id

    aspect_names = {
        "biological_process": "Biological Process",
        "molecular_function": "Molecular Function",
        "cellular_component": "Cellular Component",
    }

    rows = []

    # ---------------------------------------------------------
    # 3. Resolve every GO identifier present in BactAnnot
    # ---------------------------------------------------------

    for go_id in go_ids:
        lookup_id = go_id
        status = "unresolved"

        if go_id in alternate_ids:
            lookup_id = alternate_ids[go_id]
            status = "alternate_id"

        original_term = terms.get(lookup_id)

        canonical_go_id = ""
        name = ""
        aspect = ""
        definition = ""
        original_is_obsolete = False
        replacement_go_id = ""
        consider_go_ids = ""
        prokaryote_slim_go_ids = ""
        prokaryote_slim_names = ""
        resolved = False

        if original_term is not None:
            original_is_obsolete = bool(original_term["is_obsolete"])

            consider_go_ids = "; ".join(sorted(original_term["consider"]))

            current_term = original_term
            current_go_id = lookup_id

            # If an obsolete term has one explicit current replacement,
            # use that replacement for report-level interpretation.
            if original_is_obsolete:
                valid_replacements = []

                for replacement_id in original_term["replaced_by"]:
                    replacement_term = terms.get(replacement_id)

                    if replacement_term is not None and not replacement_term["is_obsolete"]:
                        valid_replacements.append(replacement_id)

                if len(valid_replacements) == 1:
                    replacement_go_id = valid_replacements[0]
                    current_go_id = replacement_go_id
                    current_term = terms[replacement_go_id]
                    status = "obsolete_replaced"
                    resolved = True
                else:
                    status = "obsolete_unresolved"

            else:
                resolved = True

                if status != "alternate_id":
                    status = "current"

            if resolved:
                canonical_go_id = current_go_id
                name = str(current_term["name"])
                aspect = aspect_names.get(
                    str(current_term["namespace"]),
                    str(current_term["namespace"]),
                )
                definition = str(current_term["definition"])

                # -------------------------------------------------
                # 4. Map the current term to prokaryote GO slim terms
                # -------------------------------------------------
                #
                # GO slim categories are broad parent terms. The
                # traversal uses is_a and part_of relationships from
                # go-basic. A protein may map to several slim terms.

                terms_to_visit = [current_go_id]
                visited_terms = set()
                slim_ids = set()

                while terms_to_visit:
                    term_id = terms_to_visit.pop()

                    if term_id in visited_terms:
                        continue

                    visited_terms.add(term_id)

                    term = terms.get(term_id)

                    if term is None:
                        continue

                    if "goslim_prokaryote" in term["subsets"]:
                        slim_ids.add(term_id)

                    for parent_id in term["parents"]:
                        if parent_id not in visited_terms:
                            terms_to_visit.append(parent_id)

                sorted_slim_ids = sorted(slim_ids)

                slim_names = []

                for slim_id in sorted_slim_ids:
                    slim_term = terms.get(slim_id)

                    if slim_term is not None:
                        slim_names.append(str(slim_term["name"]))

                prokaryote_slim_go_ids = "; ".join(sorted_slim_ids)
                prokaryote_slim_names = " | ".join(slim_names)

        rows.append(
            {
                "go_id": go_id,
                "canonical_go_id": canonical_go_id,
                "name": name,
                "aspect": aspect,
                "definition": definition,
                "status": status,
                "original_is_obsolete": (original_is_obsolete),
                "replacement_go_id": (replacement_go_id),
                "consider_go_ids": (consider_go_ids),
                "prokaryote_slim_go_ids": (prokaryote_slim_go_ids),
                "prokaryote_slim_names": (prokaryote_slim_names),
                "resolved": resolved,
            }
        )

    return pd.DataFrame(
        rows,
        columns=output_columns,
    )


# Summarizes protein coverage using broad prokaryote GO slim categories
def build_go_slim_profile(
    annotations: pd.DataFrame,
    go_enrichment: pd.DataFrame,
) -> pd.DataFrame:
    """Build a protein-level summary of prokaryote GO slim categories."""

    output_columns = [
        "aspect",
        "slim_go_id",
        "slim_name",
        "protein_count",
    ]

    go_annotations = annotations.loc[
        annotations["annotation_type"].eq("GO") & annotations["annotation_id"].fillna("").ne(""),
        [
            "feature_uid",
            "annotation_id",
        ],
    ].drop_duplicates()

    if go_annotations.empty or go_enrichment.empty:
        return pd.DataFrame(columns=output_columns)

    # ---------------------------------------------------------
    # 1. Convert the semicolon-separated GO slim mappings
    #    into one row per original GO term and slim category
    # ---------------------------------------------------------

    slim_mapping_rows = []

    for row in go_enrichment.itertuples(index=False):
        if not bool(row.resolved):
            continue

        slim_ids_text = str(row.prokaryote_slim_go_ids)
        slim_names_text = str(row.prokaryote_slim_names)

        if not slim_ids_text:
            continue

        slim_ids = slim_ids_text.split("; ")
        slim_names = slim_names_text.split(" | ")

        for index, slim_id in enumerate(slim_ids):
            slim_name = ""

            if index < len(slim_names):
                slim_name = slim_names[index]

            slim_mapping_rows.append(
                {
                    "go_id": str(row.go_id),
                    "aspect": str(row.aspect),
                    "slim_go_id": slim_id,
                    "slim_name": slim_name,
                }
            )

    slim_mapping = pd.DataFrame(
        slim_mapping_rows,
        columns=[
            "go_id",
            "aspect",
            "slim_go_id",
            "slim_name",
        ],
    )

    if slim_mapping.empty:
        return pd.DataFrame(columns=output_columns)

    # ---------------------------------------------------------
    # 2. Attach slim categories to the proteins annotated
    #    with each original GO identifier
    # ---------------------------------------------------------

    mapped_annotations = go_annotations.merge(
        slim_mapping,
        left_on="annotation_id",
        right_on="go_id",
        how="inner",
    )

    mapped_annotations = mapped_annotations[
        [
            "feature_uid",
            "aspect",
            "slim_go_id",
            "slim_name",
        ]
    ].drop_duplicates()

    # ---------------------------------------------------------
    # 3. Count distinct proteins in every slim category
    # ---------------------------------------------------------

    profile = (
        mapped_annotations.groupby(
            [
                "aspect",
                "slim_go_id",
                "slim_name",
            ]
        )["feature_uid"]
        .nunique()
        .reset_index(name="protein_count")
        .sort_values(
            [
                "aspect",
                "protein_count",
                "slim_name",
            ],
            ascending=[
                True,
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )

    return profile[output_columns]
