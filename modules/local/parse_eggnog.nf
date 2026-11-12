process PARSE_EGGNOG {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(annotations),
          path(seed_orthologs),
          path(features)

    output:
    tuple val(sample_id),
          path("${sample_id}.eggnog.annotations.tsv"),
          emit: annotations

    tuple val(sample_id),
          path("${sample_id}.eggnog.evidence.tsv"),
          emit: evidence

    script:
    """
    bactannot parse-eggnog \
        --annotations "${annotations}" \
        --seed-orthologs "${seed_orthologs}" \
        --features "${features}" \
        --sample-id "${sample_id}" \
        --outdir .
    """
}