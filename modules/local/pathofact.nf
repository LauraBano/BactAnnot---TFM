process PATHOFACT {

    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(faa)
        path pathofact_config
        path pathofact_prepare_script
        path pathofact_assets
        val pathofact_path
        val pathofact_db
        val pathofact_conda_prefix

    output:
    tuple val(sample_id),
        path("${sample_id}.pathofact.clusters.tsv"),
        path("${sample_id}.pathofact.vf.ml.tsv"),
        path("${sample_id}.pathofact.vf.hmm"),
        path("${sample_id}.pathofact.vf.high.tsv"),
        path("${sample_id}.pathofact.vf.ambiguous.tsv"),
        path("${sample_id}.pathofact.tox.ml.tsv"),
        path("${sample_id}.pathofact.tox.hmm"),
        path("${sample_id}.pathofact.tox.high.tsv"),
        path("${sample_id}.pathofact.tox.ambiguous.tsv"),
        path("${sample_id}.pathofact.amr.deeparg.tsv"),
        path("${sample_id}.pathofact.amr.amrfinder.tsv"),
        path("${sample_id}.pathofact.amr.rgi.tsv"),
        path("${sample_id}.pathofact.amr.summary.tsv"),
        emit: results

    /*
    * Build a task-local PathoFact2 runtime so BactAnnot can
    * apply version-controlled compatibility adaptations
    * without modifying the external PathoFact2 installation.
    */

    script:
    """
    mkdir -p \
        pathofact_input \
        pathofact_runtime

    gzip -c "${faa}" \
        > "pathofact_input/${sample_id}.faa.gz"

    ln -s \
        "${pathofact_path}/Main.smk" \
        pathofact_runtime/Main.smk

    cp -a \
        "${pathofact_path}/src" \
        pathofact_runtime/src

    cp -a \
        "${pathofact_path}/rules" \
        pathofact_runtime/rules

    cp -a \
        "${pathofact_path}/conda_envs" \
        pathofact_runtime/conda_envs

    python "${pathofact_prepare_script}" \
        --runtime pathofact_runtime \
        --assets "${pathofact_assets}"

    printf '%s,%s,protein\\n' \
        "${sample_id}" \
        "\$PWD/pathofact_input/${sample_id}.faa.gz" \
        > pathofact_runtime/list_of_samples.csv

    cp \
        "${pathofact_config}" \
        pathofact_runtime/Config.yaml

    sed -i \
        -e "s|__WORK_DIR__|\$PWD/pathofact_runtime|g" \
        -e "s|__SAMPLES_INFO__|\$PWD/pathofact_runtime/list_of_samples.csv|g" \
        -e "s|__PATHOFACT_PATH__|${pathofact_path}|g" \
        -e "s|__PATHOFACT_DB__|${pathofact_db}|g" \
        pathofact_runtime/Config.yaml

    cd pathofact_runtime

    snakemake \
        -s Main.smk \
        --use-conda \
        --conda-prefix "${pathofact_conda_prefix}" \
        --conda-frontend conda \
        --cores ${task.cpus} \
        --rerun-incomplete \
        --printshellcmds

    cd ..

    cp \
        "pathofact_runtime/results/Cluster/${sample_id}/clean_cluster.tsv" \
        "${sample_id}.pathofact.clusters.tsv"

    cp \
        "pathofact_runtime/results/ML_predictions/VF/${sample_id}_classifier_prediction.tsv" \
        "${sample_id}.pathofact.vf.ml.tsv"

    cp \
        "pathofact_runtime/results/HMM_virulence/${sample_id}.hmmscan" \
        "${sample_id}.pathofact.vf.hmm"

    cp \
        "pathofact_runtime/results/summary/${sample_id}/VF_hits_rep_prot.tsv" \
        "${sample_id}.pathofact.vf.high.tsv"

    cp \
        "pathofact_runtime/results/summary/${sample_id}/ambiguous_VF_hits_rep_prot.tsv" \
        "${sample_id}.pathofact.vf.ambiguous.tsv"

    cp \
        "pathofact_runtime/results/ML_predictions/TOX/${sample_id}_classifier_prediction.tsv" \
        "${sample_id}.pathofact.tox.ml.tsv"

    cp \
        "pathofact_runtime/results/HMM_toxin/${sample_id}.hmmscan" \
        "${sample_id}.pathofact.tox.hmm"

    cp \
        "pathofact_runtime/results/summary/${sample_id}/TOX_hits_rep_prot.tsv" \
        "${sample_id}.pathofact.tox.high.tsv"

    cp \
        "pathofact_runtime/results/summary/${sample_id}/ambiguous_TOX_hits_rep_prot.tsv" \
        "${sample_id}.pathofact.tox.ambiguous.tsv"

    cp \
        "pathofact_runtime/results/DeepARG/${sample_id}.argNorm_normalised.tsv" \
        "${sample_id}.pathofact.amr.deeparg.tsv"

    cp \
        "pathofact_runtime/results/amrfinder/${sample_id}.argNorm_normalised.tsv" \
        "${sample_id}.pathofact.amr.amrfinder.tsv"

    cp \
        "pathofact_runtime/results/RGI/${sample_id}.tsv" \
        "${sample_id}.pathofact.amr.rgi.tsv"

    cp \
        "pathofact_runtime/results/summary/${sample_id}/amr_hits_rep_prot.tsv" \
        "${sample_id}.pathofact.amr.summary.tsv"
    """
}