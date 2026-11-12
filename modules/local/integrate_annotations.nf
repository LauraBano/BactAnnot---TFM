process INTEGRATE_ANNOTATIONS {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(annotation_tables),
          path(features)

    output:
    tuple val(sample_id),
          path("${sample_id}.integrated.annotations.tsv"),
          emit: annotations

    script:
    def annotation_args = annotation_tables
        .collect { annotation_file -> "\"${annotation_file}\"" }
        .join(" ")

    """
    bactannot integrate-annotations \
        --features "${features}" \
        --annotations ${annotation_args} \
        --sample-id "${sample_id}" \
        --outdir .
    """
}