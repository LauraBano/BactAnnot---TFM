process PARSE_KOFAMSCAN {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(kofamscan_tsv),
          path(features)

    output:
    tuple val(sample_id),
          path("${sample_id}.kofamscan.annotations.tsv"),
          emit: annotations

    tuple val(sample_id),
          path("${sample_id}.kofamscan.evidence.tsv"),
          emit: evidence

    script:
    """
    bactannot parse-kofamscan \
        --input "${kofamscan_tsv}" \
        --features "${features}" \
        --sample-id "${sample_id}" \
        --outdir .
    """
}