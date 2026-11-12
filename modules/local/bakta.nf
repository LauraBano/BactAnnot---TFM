process BAKTA {

    tag "${sample_id}"

    cpus 4
    memory '8 GB'

    input:
    tuple val(sample_id), path(fasta)
    val bakta_db

    output:
    tuple val(sample_id),
          path("${sample_id}.json"),
          path("${sample_id}.faa"),
          emit: results

    script:
    """
    bakta \
        --db "${bakta_db}" \
        --output bakta_output \
        --force \
        --prefix "${sample_id}" \
        --threads ${task.cpus} \
        "${fasta}"

    cp "bakta_output/${sample_id}.json" "${sample_id}.json"
    cp "bakta_output/${sample_id}.faa" "${sample_id}.faa"
    """
}