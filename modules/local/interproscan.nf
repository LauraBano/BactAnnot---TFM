process INTERPROSCAN {

    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(faa)
    val interproscan_path

    output:
    tuple val(sample_id),
          path("${sample_id}.interproscan.tsv"),
          emit: results

    script:
    """
    mkdir -p interproscan_tmp

    "${interproscan_path}" \
        -i "${faa}" \
        -f TSV \
        --goterms \
        -cpu ${task.cpus} \
        -T interproscan_tmp \
        -o "${sample_id}.interproscan.tsv"
    """
}