process PARSE_INTERPROSCAN {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(interproscan_tsv),
          path(features)

    output:
    tuple val(sample_id),
          path("${sample_id}.interproscan.annotations.tsv"),
          emit: annotations

    tuple val(sample_id),
          path("${sample_id}.interproscan.evidence.tsv"),
          emit: evidence

    script:
    """
    bactannot parse-interproscan \
        --input "${interproscan_tsv}" \
        --features "${features}" \
        --sample-id "${sample_id}" \
        --outdir .
    """
}