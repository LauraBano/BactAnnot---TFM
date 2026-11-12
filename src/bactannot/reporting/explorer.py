"""Interactive protein explorer for BactAnnot reports."""

import re
from html import escape
from pathlib import Path

import pandas as pd

from bactannot.reporting.report import build_feature_explorer

DETAIL_TOOLS = (
    "Bakta",
    "eggNOG-mapper",
    "KofamScan",
    "InterProScan",
    "PathoFact2",
)

INTERPRETATION_STATUS_LABELS = {
    "ko_function_proposed": "KO function proposed",
    "ec_activity_proposed": "Enzyme activity proposed",
    "family_description_proposed": "Domain/family description proposed",
    "multiple_family_descriptions_proposed": "Multiple domain/family descriptions proposed",
    "molecular_function_properties_only": "Molecular function properties only",
    "multiple_ec_activities_proposed": "Multiple enzyme activities proposed",
    "ko_ec_relation_unresolved": "KO-EC relation unresolved",
    "no_selected_proposal": "No proposal selected",
}

PRIMARY_SCOPE_LABELS = {
    "protein_function": "Protein function",
    "enzyme_activity": "Enzyme activity",
    "domain_or_family": "Domain/family description",
    "": "No primary proposal",
}

CONSENSUS_DETAIL_TABLES = (
    "ko_claims",
    "ec_claims",
    "domain_claims",
    "ko_ec_relations",
    "ec_go_relations",
)


def _tool_name(source: str) -> str:
    """Return the main BactAnnot tool for an annotation/evidence source."""

    source = str(source)

    if source == "bakta":
        return "Bakta"

    if source.startswith("eggNOG-mapper"):
        return "eggNOG-mapper"

    if source.startswith("KofamScan"):
        return "KofamScan"

    if source.startswith("InterProScan"):
        return "InterProScan"

    if source.startswith("PathoFact2"):
        return "PathoFact2"

    return source


def _safe_filename(value: str) -> str:
    """Convert a protein identifier into a safe HTML filename."""

    cleaned = re.sub(
        r"[^A-Za-z0-9_.-]",
        "_",
        str(value),
    )

    if not cleaned:
        cleaned = "protein"

    return f"{cleaned}.html"


def _table_html(dataframe: pd.DataFrame) -> str:
    """Render a compact HTML table and hide columns that are fully empty."""

    if dataframe.empty:
        return '<p class="empty-note">No records for this protein.</p>'

    table = dataframe.copy()

    visible_columns = []

    for column in table.columns:
        values = table[column].fillna("").astype(str).str.strip()

        if values.ne("").any():
            visible_columns.append(column)

    if not visible_columns:
        return '<p class="empty-note">No records for this protein.</p>'

    table = table[visible_columns]

    return (
        '<div class="table-wrap">'
        + table.to_html(
            index=False,
            escape=True,
            border=0,
        )
        + "</div>"
    )


def _vertical_summary_html(values: dict[str, object]) -> str:
    """Render one field per row and collapse very long multi-value fields."""

    rows = []

    for field, value in values.items():
        text = str(value).strip()

        if not text:
            continue

        parts = [part.strip() for part in text.split(";") if part.strip()]

        if len(parts) > 6:
            items = "".join(f"<li>{escape(part)}</li>" for part in parts)

            value_html = (
                '<details class="value-details">'
                f"<summary>{len(parts)} values</summary>"
                f'<ul class="value-list">{items}</ul>'
                "</details>"
            )
        else:
            value_html = escape(text)

        rows.append(f"<tr><th>{escape(field)}</th><td>{value_html}</td></tr>")

    if not rows:
        return '<p class="empty-note">No information available.</p>'

    return (
        '<div class="table-wrap vertical-summary-wrap">'
        '<table class="vertical-summary">'
        "<tbody>" + "".join(rows) + "</tbody></table></div>"
    )


def _consensus_detail_html(
    explorer_row: pd.Series,
    consensus_rows: dict[str, pd.DataFrame],
) -> str:
    """Display stored claims and relationships without inferring new support."""

    if not consensus_rows:
        return (
            '<p class="empty-note">Technical claim and relationship tables '
            "were not loaded. Native tool records remain available below.</p>"
        )

    primary_type = str(explorer_row.get("primary_claim_type", ""))
    primary_id = str(explorer_row.get("primary_claim_id", ""))
    selected_ids = {
        value
        for field in (
            "ko_proposal_ids",
            "ec_activity_proposal_ids",
            "unlinked_ec_activity_ids",
            "family_proposal_ids",
        )
        for value in str(explorer_row.get(field, "")).split(";")
        if value
    }
    if primary_id:
        selected_ids.add(primary_id)

    claims = []

    for name in ("ko_claims", "ec_claims", "domain_claims"):
        if name not in consensus_rows:
            continue

        table = consensus_rows[name].copy()
        if name == "domain_claims":
            table = table.loc[table["claim_id"].isin(selected_ids)].copy()
        if table.empty:
            continue

        table["Interpretation role"] = "Other observed claim"
        table.loc[table["claim_id"].isin(selected_ids), "Interpretation role"] = (
            "Additional proposal"
        )
        primary = table["claim_type"].eq(primary_type) & table["claim_id"].eq(primary_id)
        table.loc[primary, "Interpretation role"] = "Primary proposal"
        claims.append(table)

    sections = [
        "<h3>KO, EC and proposed domain/family claims</h3>",
        (
            "<p>Claims group observations about the same identifier. Observing sources "
            "and supporting evidence are distinct. Other observed claims have not been "
            "promoted to proposals by this report. Domain/family rows are limited to "
            "the proposal identifiers recorded in the interpretation.</p>"
        ),
    ]

    if claims:
        claim_table = pd.concat(claims, ignore_index=True).fillna("")
        columns = [
            "claim_type",
            "claim_id",
            "Interpretation role",
            "observing_sources",
            "evidence_sources",
            "evidence_lines",
            "evidence_relations",
            "independent_evidence_count",
            "context_sources",
            "unresolved_sources",
            "evidence_status",
        ]
        sections.append(
            _table_html(
                claim_table[[column for column in columns if column in claim_table.columns]]
            )
        )
        sections.append(
            "<p>Evidence lines and their counts are copied from the consensus tables. "
            "They describe provenance under its evidence policy, not a probability "
            "of correctness or a count of independent experiments.</p>"
        )
    elif any(name in consensus_rows for name in ("ko_claims", "ec_claims", "domain_claims")):
        sections.append('<p class="empty-note">No matching claims in the loaded tables.</p>')
    else:
        sections.append('<p class="empty-note">Claim tables were not loaded.</p>')

    sections.extend(
        [
            "<h3>Relationships between observed identifiers</h3>",
            (
                "<p>These relationships can include identifiers outside the primary proposal. "
                "A documented link records a mapping in the cited reference; hierarchical "
                "compatibility allows different levels of specificity. "
                "<code>not_documented</code> does not demonstrate biological incompatibility. "
                "<code>unmapped_ko</code> and <code>unmapped_go</code> indicate missing mappings; "
                "<code>obsolete_unresolved</code> indicates that the GO term could not be "
                "resolved for comparison. None is converted into a conflict by the report.</p>"
            ),
        ]
    )

    for name, title in (("ko_ec_relations", "KO–EC"), ("ec_go_relations", "EC–GO")):
        sections.append(f"<h4>{title}</h4>")
        if name not in consensus_rows:
            sections.append('<p class="empty-note">Relationship table was not loaded.</p>')
            continue
        table = consensus_rows[name]
        columns = [
            column
            for column in table.columns
            if column
            not in (
                "sample_id",
                "feature_uid",
                "protein_id",
            )
        ]
        sections.append(_table_html(table[columns]))

    return "\n".join(sections)


def _write_detail_page(
    feature: pd.Series,
    explorer_row: pd.Series,
    annotation_rows: pd.DataFrame,
    evidence_rows: pd.DataFrame,
    concordance_rows: pd.DataFrame,
    output_path: Path,
    consensus_rows: dict[str, pd.DataFrame] | None = None,
) -> None:
    """Write the functional interpretation and source records for one protein."""

    protein_id = str(feature["protein_id"])

    metadata_columns = [
        "sample_id",
        "feature_uid",
        "protein_id",
        "gene",
        "product",
        "contig_id",
        "feature_type",
        "start",
        "stop",
        "strand",
        "hypothetical",
        "pseudogene",
        "truncated",
        "bakta_id",
        "locus_tag",
        "aa_md5",
        "frame",
        "start_type",
        "rbs_motif",
    ]

    metadata = {
        column: feature.get(column, "") for column in metadata_columns if column in feature.index
    }

    status = str(explorer_row["interpretation_status"])
    scope = str(explorer_row["primary_scope"])
    interpretation = {
        "Primary proposal": explorer_row["primary_function"] or "No primary proposal",
        "Scope": PRIMARY_SCOPE_LABELS.get(scope, scope),
        "Interpretation": INTERPRETATION_STATUS_LABELS.get(status, status),
    }
    decision_columns = (
        "primary_claim_type",
        "primary_claim_id",
        "primary_selection_rule",
        "primary_description_source",
        "primary_description_reference_id",
        "ko_proposal_ids",
        "ec_activity_proposal_ids",
        "mf_proposal_ids",
        "unlinked_ec_activity_ids",
        "family_proposal_ids",
    )
    decision = {column: explorer_row.get(column, "") for column in decision_columns}

    summary_columns = [
        "annotation_sources",
        "COG_CATEGORY",
        "GO",
        "EC",
        "KEGG_KO",
        "KEGG_PATHWAY",
        "KEGG_MODULE",
        "KEGG_REACTION",
        "INTERPRO",
        "PFAM",
        "SMART",
        "CAZY",
        "amr",
        "virulence",
        "toxin",
        "concordance_exact",
        "concordance_partial",
        "concordance_disjoint",
    ]

    summary = {
        column: explorer_row.get(column, "")
        for column in summary_columns
        if column in explorer_row.index
    }

    concordance_columns = [
        "annotation_type",
        "source_1",
        "source_2",
        "source_1_ids",
        "source_2_ids",
        "shared_ids",
        "relation",
    ]

    concordance_columns = [
        column for column in concordance_columns if column in concordance_rows.columns
    ]

    page = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>BactAnnot — {escape(protein_id)}</title>",
        '<link rel="stylesheet" href="detail.css">',
        "</head>",
        "<body>",
        f"<h1>{escape(protein_id)}</h1>",
        '<p class="subtitle">Functional interpretation and source records for one protein.</p>',
        '<div class="note">',
        (
            "Proposals are copied from the consensus output. Tool annotations "
            "and evidence below provide context; they are not all supporting "
            "evidence for the primary proposal."
        ),
        "</div>",
        "<h2>Functional interpretation</h2>",
        _vertical_summary_html(interpretation),
        "<details>",
        "<summary>Decision details and additional proposals</summary>",
        '<div class="detail-content">',
        _vertical_summary_html(decision),
        "</div>",
        "</details>",
        "<details>",
        "<summary>Claim support and functional relationships</summary>",
        '<div class="detail-content">',
        _consensus_detail_html(explorer_row, consensus_rows or {}),
        "</div>",
        "</details>",
        "<h2>Structural identity</h2>",
        _vertical_summary_html(metadata),
        "<h2>Integrated functional summary</h2>",
        _vertical_summary_html(summary),
        "<details>",
        "<summary>Pairwise concordance</summary>",
        '<div class="detail-content">',
        _table_html(concordance_rows[concordance_columns]),
        "</div>",
        "</details>",
    ]

    for tool in DETAIL_TOOLS:
        tool_annotation_rows = annotation_rows.loc[
            annotation_rows["source"].map(_tool_name).eq(tool)
        ].copy()

        tool_evidence_rows = evidence_rows.loc[
            evidence_rows["source"].map(_tool_name).eq(tool)
        ].copy()

        annotation_columns = [
            "source",
            "annotation_type",
            "annotation_id",
            "annotation_value",
            "description",
        ]

        annotation_columns = [
            column for column in annotation_columns if column in tool_annotation_rows.columns
        ]

        evidence_columns = [
            "source",
            "evidence_type",
            "accession",
            "gene",
            "product",
            "rank",
            "score",
            "native_threshold",
            "evalue",
            "identity",
            "query_coverage",
            "subject_coverage",
            "aa_coverage",
            "hmm_coverage",
            "start",
            "stop",
        ]

        evidence_columns = [
            column for column in evidence_columns if column in tool_evidence_rows.columns
        ]

        page.extend(
            [
                "<details>",
                (
                    "<summary>"
                    f"{escape(tool)} — "
                    f"{len(tool_annotation_rows)} annotation rows, "
                    f"{len(tool_evidence_rows)} evidence rows"
                    "</summary>"
                ),
                '<div class="detail-content">',
                "<h3>Annotations</h3>",
                _table_html(tool_annotation_rows[annotation_columns]),
                "<h3>Evidence</h3>",
                _table_html(tool_evidence_rows[evidence_columns]),
                "</div>",
                "</details>",
            ]
        )

    page.extend(
        [
            "</body>",
            "</html>",
        ]
    )

    output_path.write_text(
        "\n".join(page),
        encoding="utf-8",
    )


def build_protein_explorer_html(
    features: pd.DataFrame,
    annotations: pd.DataFrame,
    evidence: pd.DataFrame,
    interpretations: pd.DataFrame,
    concordance: pd.DataFrame,
    output_dir: Path,
    detail_link_prefix: str,
    consensus_tables: dict[str, pd.DataFrame] | None = None,
) -> str:
    """Build a searchable protein table and write one detail page per protein."""

    protein_explorer = build_feature_explorer(
        features,
        annotations,
        interpretations,
        concordance,
    )

    consensus_tables = consensus_tables or {}
    consensus_groups = {
        name: table.groupby(["sample_id", "feature_uid", "protein_id"], sort=False)
        for name, table in consensus_tables.items()
        if name in CONSENSUS_DETAIL_TABLES
    }

    output_dir = Path(output_dir)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for old_page in output_dir.glob("*.html"):
        old_page.unlink()

    detail_css = """
body {
    max-width: 1500px;
    margin: 0 auto;
    padding: 1.5rem;
    font-family: system-ui, sans-serif;
    color: #212529;
}
h1, h2, h3 { line-height: 1.2; }
.subtitle { color: #6c757d; }
.note {
    padding: 0.8rem 1rem;
    background: #f8f9fa;
    border-left: 4px solid #6c757d;
    margin: 1rem 0;
}
details {
    border: 1px solid #dee2e6;
    border-radius: 0.45rem;
    margin: 0.8rem 0;
}
summary {
    cursor: pointer;
    padding: 0.75rem 1rem;
    font-weight: 700;
    background: #f8f9fa;
}
.detail-content { padding: 0.8rem 1rem 1rem 1rem; }
.table-wrap {
    overflow: auto;
    max-height: 34rem;
    border: 1px solid #e9ecef;
    margin: 0.5rem 0 1rem 0;
}
table {
    border-collapse: collapse;
    width: 100%;
    font-size: 0.82rem;
}
th, td {
    padding: 0.38rem 0.45rem;
    border-bottom: 1px solid #e9ecef;
    text-align: left;
    vertical-align: top;
}
th {
    position: sticky;
    top: 0;
    background: white;
}
.empty-note {
    color: #6c757d;
    font-size: 0.9rem;
}
.vertical-summary-wrap {
    max-height: none;
    overflow-x: hidden;
}
.vertical-summary {
    table-layout: fixed;
}
.vertical-summary th {
    position: static;
    width: 220px;
    background: #f8f9fa;
}
.vertical-summary td {
    overflow-wrap: anywhere;
}
.value-details {
    border: 0;
    margin: 0;
}
.value-details summary {
    display: inline-block;
    padding: 0.15rem 0.45rem;
    background: #f1f3f5;
    border-radius: 0.3rem;
    font-weight: 600;
}
.value-list {
    margin: 0.5rem 0 0 1.2rem;
    columns: 2;
    column-gap: 2rem;
}
@media (max-width: 800px) {
    .vertical-summary th { width: 150px; }
    .value-list { columns: 1; }
}
"""

    (output_dir / "detail.css").write_text(
        detail_css,
        encoding="utf-8",
    )

    protein_features = features.loc[features["protein_id"].fillna("").ne("")].copy()

    feature_rows = protein_features.set_index("feature_uid")
    explorer_rows = protein_explorer.set_index("feature_uid")

    annotations_by_feature = annotations.groupby(
        "feature_uid",
        sort=False,
    )

    evidence_by_feature = evidence.groupby(
        "feature_uid",
        sort=False,
    )

    if concordance.empty:
        concordance_by_feature = None
    else:
        concordance_by_feature = concordance.groupby(
            "feature_uid",
            sort=False,
        )

    detail_file_by_feature = {}

    for feature_uid in protein_explorer["feature_uid"]:
        if feature_uid not in feature_rows.index:
            continue

        feature = feature_rows.loc[feature_uid].copy()
        feature["feature_uid"] = feature_uid

        explorer_row = explorer_rows.loc[feature_uid]

        protein_id = str(feature["protein_id"])
        detail_filename = _safe_filename(protein_id)

        detail_file_by_feature[feature_uid] = detail_filename

        if feature_uid in annotations_by_feature.groups:
            annotation_rows = annotations_by_feature.get_group(feature_uid)
        else:
            annotation_rows = annotations.iloc[0:0]

        if feature_uid in evidence_by_feature.groups:
            evidence_rows = evidence_by_feature.get_group(feature_uid)
        else:
            evidence_rows = evidence.iloc[0:0]

        if concordance_by_feature is not None and feature_uid in concordance_by_feature.groups:
            concordance_rows = concordance_by_feature.get_group(feature_uid)
        else:
            concordance_rows = concordance.iloc[0:0]

        _write_detail_page(
            feature=feature,
            explorer_row=explorer_row,
            annotation_rows=annotation_rows,
            evidence_rows=evidence_rows,
            concordance_rows=concordance_rows,
            output_path=(output_dir / detail_filename),
            consensus_rows={
                name: (
                    groups.get_group((str(feature["sample_id"]), feature_uid, protein_id))
                    if (str(feature["sample_id"]), feature_uid, protein_id) in groups.indices
                    else consensus_tables[name].iloc[0:0]
                )
                for name, groups in consensus_groups.items()
            },
        )

    table_rows = []

    for row in protein_explorer.itertuples(index=False):
        protein_id = str(row.protein_id)

        if not protein_id:
            continue

        feature_uid = str(row.feature_uid)
        detail_filename = detail_file_by_feature[feature_uid]
        detail_url = detail_link_prefix.rstrip("/") + "/" + detail_filename

        gene = str(getattr(row, "gene", ""))
        product = str(getattr(row, "product", ""))
        hypothetical = bool(getattr(row, "hypothetical", False))
        amr = bool(getattr(row, "amr", False))
        virulence = str(getattr(row, "virulence", ""))
        toxin = str(getattr(row, "toxin", ""))
        primary_function = str(row.primary_function)
        primary_scope = str(row.primary_scope)
        status = str(row.interpretation_status)
        scope_label = PRIMARY_SCOPE_LABELS.get(primary_scope, primary_scope)
        status_label = INTERPRETATION_STATUS_LABELS.get(status, status)

        ko = str(getattr(row, "KEGG_KO", ""))
        interpro = str(getattr(row, "INTERPRO", ""))
        pfam = str(getattr(row, "PFAM", ""))
        go = str(getattr(row, "GO", ""))
        ec = str(getattr(row, "EC", ""))
        cazy = str(getattr(row, "CAZY", ""))

        search_text = (
            f"{protein_id} {feature_uid} {gene} {product} "
            f"{primary_function} {scope_label} {status_label} "
            f"{ko} {interpro} {pfam} {go} {ec} {cazy} "
            f"{virulence} {toxin}"
        ).lower()

        table_rows.append(
            '<tr class="protein-row" '
            f'data-search="{escape(search_text, quote=True)}" '
            f'data-hypothetical="{str(hypothetical).lower()}" '
            f'data-amr="{str(amr).lower()}" '
            f'data-virulence="{escape(virulence, quote=True)}" '
            f'data-toxin="{escape(toxin, quote=True)}" '
            f'data-interpretation="{escape(status, quote=True)}">'
            "<td>"
            f'<a href="{escape(detail_url, quote=True)}" target="_blank">'
            f"{escape(protein_id)}"
            "</a>"
            "</td>"
            f"<td>{escape(gene) if gene else '—'}</td>"
            f"<td>{escape(primary_function or 'No primary proposal')}</td>"
            f"<td>{escape(scope_label)}</td>"
            f"<td>{escape(status_label)}</td>"
            f"<td>{escape(product) if product else '—'}</td>"
            f"<td>{'Yes' if hypothetical else 'No'}</td>"
            f"<td>{'Yes' if amr else 'No'}</td>"
            f"<td>{escape(virulence) if virulence else '—'}</td>"
            f"<td>{escape(toxin) if toxin else '—'}</td>"
            "</tr>"
        )

    total_proteins = len(table_rows)
    interpretation_options = "".join(
        f'<option value="{escape(status, quote=True)}">{escape(label)}</option>'
        for status, label in INTERPRETATION_STATUS_LABELS.items()
    )

    explorer_html = f"""
<div class="bactannot-explorer">
<style>
.bactannot-explorer {{ margin: 1rem 0; }}
.explorer-controls {{
    display: grid;
    grid-template-columns: minmax(260px, 2fr) repeat(5, minmax(120px, 1fr));
    gap: 0.55rem;
    align-items: end;
    margin-bottom: 0.75rem;
}}
.explorer-controls label {{
    display: block;
    font-size: 0.78rem;
    font-weight: 600;
    margin-bottom: 0.15rem;
}}
.explorer-controls input,
.explorer-controls select,
.explorer-controls button {{
    width: 100%;
    padding: 0.42rem 0.5rem;
    border: 1px solid #ced4da;
    border-radius: 0.35rem;
    background: white;
}}
.explorer-summary {{
    margin: 0.25rem 0 0.65rem 0;
    color: #6c757d;
    font-size: 0.88rem;
}}
.explorer-table-wrap {{
    max-height: 42rem;
    overflow: auto;
    border: 1px solid #dee2e6;
    border-radius: 0.4rem;
}}
.explorer-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.80rem;
}}
.explorer-table th,
.explorer-table td {{
    padding: 0.38rem 0.45rem;
    border-bottom: 1px solid #e9ecef;
    vertical-align: top;
    text-align: left;
}}
.explorer-table th {{
    position: sticky;
    top: 0;
    z-index: 2;
    background: white;
}}
.explorer-table td:nth-child(3) {{ min-width: 260px; }}
.explorer-table td:nth-child(5),
.explorer-table td:nth-child(6) {{
    max-width: 230px;
    overflow-wrap: anywhere;
}}
.explorer-pagination {{
    display: flex;
    justify-content: center;
    align-items: center;
    gap: 0.75rem;
    margin-top: 0.75rem;
}}
.explorer-pagination button {{
    padding: 0.35rem 0.8rem;
    border: 1px solid #ced4da;
    border-radius: 0.35rem;
    background: white;
}}
.explorer-pagination button:disabled {{
    opacity: 0.45;
    cursor: not-allowed;
}}
@media (max-width: 1100px) {{
    .explorer-controls {{
        grid-template-columns: 1fr 1fr;
    }}
}}
</style>

<div class="explorer-controls">
    <div>
        <label for="protein-search">Search</label>
        <input id="protein-search" type="search"
               placeholder="protein, proposed function, gene, Bakta product, KO, GO, EC...">
    </div>
    <div>
        <label for="filter-hypothetical">Hypothetical</label>
        <select id="filter-hypothetical">
            <option value="all">All</option>
            <option value="true">Yes</option>
            <option value="false">No</option>
        </select>
    </div>
    <div>
        <label for="filter-amr">AMR</label>
        <select id="filter-amr">
            <option value="all">All</option>
            <option value="true">Yes</option>
            <option value="false">No</option>
        </select>
    </div>
    <div>
        <label for="filter-virulence">Virulence</label>
        <select id="filter-virulence">
            <option value="all">All</option>
            <option value="high_confidence">High confidence</option>
            <option value="ambiguous">Ambiguous</option>
            <option value="none">None</option>
        </select>
    </div>
    <div>
        <label for="filter-toxin">Toxin</label>
        <select id="filter-toxin">
            <option value="all">All</option>
            <option value="high_confidence">High confidence</option>
            <option value="ambiguous">Ambiguous</option>
            <option value="none">None</option>
        </select>
    </div>
    <div>
        <label for="filter-interpretation">Interpretation</label>
        <select id="filter-interpretation">
            <option value="all">All</option>
            {interpretation_options}
        </select>
    </div>
    <div>
        <label for="protein-page-size">Rows</label>
        <select id="protein-page-size">
            <option value="10">10</option>
            <option value="15" selected>15</option>
            <option value="25">25</option>
            <option value="50">50</option>
        </select>
    </div>
</div>

<div class="explorer-summary">
    Showing <strong id="protein-range">1–15</strong> of
    <strong id="visible-protein-count">{total_proteins}</strong> matching proteins
    (<strong>{total_proteins}</strong> total).
    Click a protein ID to open its interpretation and source records in a new tab.
</div>

<div class="explorer-table-wrap">
<table class="explorer-table">
<thead>
<tr>
    <th>Protein</th>
    <th>Gene</th>
    <th>Primary proposal</th>
    <th>Scope</th>
    <th>Interpretation</th>
    <th>Bakta product</th>
    <th>Hypothetical</th>
    <th>AMR</th>
    <th>Virulence</th>
    <th>Toxin</th>
</tr>
</thead>
<tbody>
{"".join(table_rows)}
</tbody>
</table>
</div>

<div class="explorer-pagination">
    <button id="protein-prev" type="button">Previous</button>
    <span id="protein-page-label"></span>
    <button id="protein-next" type="button">Next</button>
</div>

<script>
(function () {{
    const search = document.getElementById("protein-search");
    const hypothetical = document.getElementById("filter-hypothetical");
    const amr = document.getElementById("filter-amr");
    const virulence = document.getElementById("filter-virulence");
    const toxin = document.getElementById("filter-toxin");
    const interpretation = document.getElementById("filter-interpretation");
    const pageSize = document.getElementById("protein-page-size");
    const rows = Array.from(document.querySelectorAll(".protein-row"));
    const visibleCount = document.getElementById("visible-protein-count");
    const rangeLabel = document.getElementById("protein-range");
    const pageLabel = document.getElementById("protein-page-label");
    const previousButton = document.getElementById("protein-prev");
    const nextButton = document.getElementById("protein-next");

    let matchingRows = rows;
    let currentPage = 1;

    function valueMatches(value, selected) {{
        if (selected === "all") return true;
        if (selected === "none") return value === "";
        return value === selected;
    }}

    function renderPage() {{
        const size = Number(pageSize.value);
        const totalPages = Math.max(1, Math.ceil(matchingRows.length / size));

        if (currentPage > totalPages) currentPage = totalPages;

        rows.forEach((row) => {{
            row.style.display = "none";
        }});

        const start = (currentPage - 1) * size;
        const end = Math.min(start + size, matchingRows.length);

        matchingRows.slice(start, end).forEach((row) => {{
            row.style.display = "";
        }});

        visibleCount.textContent = matchingRows.length;

        if (matchingRows.length === 0) {{
            rangeLabel.textContent = "0";
        }} else {{
            rangeLabel.textContent = `${{start + 1}}–${{end}}`;
        }}

        pageLabel.textContent = `Page ${{currentPage}} of ${{totalPages}}`;
        previousButton.disabled = currentPage <= 1;
        nextButton.disabled = currentPage >= totalPages;
    }}

    function applyFilters() {{
        const query = search.value.trim().toLowerCase();

        matchingRows = rows.filter((row) => {{
            return (
                (!query || row.dataset.search.includes(query))
                && valueMatches(row.dataset.hypothetical, hypothetical.value)
                && valueMatches(row.dataset.amr, amr.value)
                && valueMatches(row.dataset.virulence, virulence.value)
                && valueMatches(row.dataset.toxin, toxin.value)
                && valueMatches(row.dataset.interpretation, interpretation.value)
            );
        }});

        currentPage = 1;
        renderPage();
    }}

    search.addEventListener("input", applyFilters);
    hypothetical.addEventListener("change", applyFilters);
    amr.addEventListener("change", applyFilters);
    virulence.addEventListener("change", applyFilters);
    toxin.addEventListener("change", applyFilters);
    interpretation.addEventListener("change", applyFilters);
    pageSize.addEventListener("change", () => {{
        currentPage = 1;
        renderPage();
    }});
    previousButton.addEventListener("click", () => {{
        currentPage -= 1;
        renderPage();
    }});
    nextButton.addEventListener("click", () => {{
        currentPage += 1;
        renderPage();
    }});

    applyFilters();
}})();
</script>
</div>
"""

    return explorer_html
