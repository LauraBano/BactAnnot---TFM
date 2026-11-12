process KOFAMSCAN {

    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(faa)
    val kofam_db

    output:
    tuple val(sample_id),
          path("${sample_id}.kofamscan.detail.tsv"),
          emit: results

    script:
    """
    mkdir -p kofamscan_tmp

    exec_annotation \
        --profile "${kofam_db}/profiles/prokaryote.hal" \
        --ko-list "${kofam_db}/ko_list" \
        --cpu ${task.cpus} \
        --tmp-dir kofamscan_tmp \
        --format detail-tsv \
        -o "${sample_id}.kofamscan.detail.tsv" \
        "${faa}"
    """
}