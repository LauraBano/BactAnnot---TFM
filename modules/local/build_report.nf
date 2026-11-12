process BUILD_REPORT {

    tag "${sample_id}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(sample_id),
          path(features),
          path(annotations),
          path(evidence),
          path(consensus_dir),
          path(concordance)

    path report_template
    path source_dir

    output:
    tuple val(sample_id),
          path("${sample_id}.report", type: 'dir'),
          emit: reports

    script:
    """
    export QUARTO_PYTHON="\$(command -v python)"
    export MPLBACKEND=Agg
    export PYTHONPATH="\$PWD/${source_dir}"

    mkdir -p "${sample_id}.report"
    cp "${report_template}" "${sample_id}.report/index.qmd"

    quarto render "${sample_id}.report/index.qmd" --to html \\
        -P "sample_id:'${sample_id}'" \\
        -P "features_file:\$PWD/${features}" \\
        -P "annotations_file:\$PWD/${annotations}" \\
        -P "evidence_file:\$PWD/${evidence}" \\
        -P "consensus_dir:\$PWD/${consensus_dir}" \\
        -P "concordance_file:\$PWD/${concordance}" \\
        -P "report_output_dir:\$PWD/${sample_id}.report"

    test -s "${sample_id}.report/index.html"
    rm "${sample_id}.report/index.qmd"
    """
}