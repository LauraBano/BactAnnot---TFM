process PARSE_BAKTA {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bakta_json)

    output:
    tuple val(sample_id),
          path("${sample_id}.bakta.features.tsv"),
          emit: features

    tuple val(sample_id),
          path("${sample_id}.bakta.annotations.tsv"),
          emit: annotations

    tuple val(sample_id),
          path("${sample_id}.bakta.evidence.tsv"),
          emit: evidence

    script:
    """
    bactannot parse-bakta \
        --input "${bakta_json}" \
        --sample-id "${sample_id}" \
        --outdir .
    """
}