process EGGNOG {

    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(faa)
    val eggnog_db

    output:
    tuple val(sample_id),
          path("${sample_id}.emapper.annotations"),
          path("${sample_id}.emapper.seed_orthologs"),
          path("${sample_id}.emapper.hits"),
          emit: results

    script:
    """
    mkdir -p eggnog_output

    emapper.py \
        -i "${faa}" \
        --itype proteins \
        -m diamond \
        --data_dir "${eggnog_db}" \
        --output "${sample_id}" \
        --output_dir eggnog_output \
        --cpu ${task.cpus} \
        --override

    cp \
        "eggnog_output/${sample_id}.emapper.annotations" \
        "${sample_id}.emapper.annotations"

    cp \
        "eggnog_output/${sample_id}.emapper.seed_orthologs" \
        "${sample_id}.emapper.seed_orthologs"

    cp \
        "eggnog_output/${sample_id}.emapper.hits" \
        "${sample_id}.emapper.hits"
    """
}