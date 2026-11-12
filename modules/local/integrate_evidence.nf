process INTEGRATE_EVIDENCE {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(evidence_tables),
          path(features)

    output:
    tuple val(sample_id),
          path("${sample_id}.integrated.evidence.tsv"),
          emit: evidence

    script:
    def evidence_args = evidence_tables
        .collect { evidence_file -> "\"${evidence_file}\"" }
        .join(" ")

    """
    bactannot integrate-evidence \
        --features "${features}" \
        --evidence ${evidence_args} \
        --sample-id "${sample_id}" \
        --outdir .
    """
}