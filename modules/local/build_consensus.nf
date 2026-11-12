process BUILD_CONSENSUS {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(annotations),
          path(evidence),
          path(features)

    path ko_list

    output:
    tuple val(sample_id),
          path("${sample_id}.consensus", type: 'dir'),
          emit: consensus

    script:
    """
    bactannot build-consensus \
        --features "${features}" \
        --annotations "${annotations}" \
        --evidence "${evidence}" \
        --sample-id "${sample_id}" \
        --ko-list "${ko_list}" \
        --outdir "${sample_id}.consensus"
    """
}