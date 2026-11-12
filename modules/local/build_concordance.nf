process BUILD_CONCORDANCE {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(annotations),
          path(features)

    output:
    tuple val(sample_id),
          path("${sample_id}.concordance.tsv"),
          path("${sample_id}.concordance_summary.tsv"),
          emit: concordance

    script:
    """
    bactannot build-concordance \
        --features "${features}" \
        --annotations "${annotations}" \
        --sample-id "${sample_id}" \
        --outdir .
    """
}