process PARSE_PATHOFACT {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(clusters),
          path(vf_ml),
          path(vf_hmm),
          path(vf_high),
          path(vf_ambiguous),
          path(tox_ml),
          path(tox_hmm),
          path(tox_high),
          path(tox_ambiguous),
          path(amr_deeparg),
          path(amr_amrfinder),
          path(amr_rgi),
          path(amr_summary),
          path(features)

    output:
    tuple val(sample_id),
          path("${sample_id}.pathofact.annotations.tsv"),
          emit: annotations

    tuple val(sample_id),
          path("${sample_id}.pathofact.evidence.tsv"),
          emit: evidence

    script:
    """
    bactannot parse-pathofact \
        --features "${features}" \
        --clusters "${clusters}" \
        --vf-ml "${vf_ml}" \
        --vf-hmm "${vf_hmm}" \
        --vf-high "${vf_high}" \
        --vf-ambiguous "${vf_ambiguous}" \
        --tox-ml "${tox_ml}" \
        --tox-hmm "${tox_hmm}" \
        --tox-high "${tox_high}" \
        --tox-ambiguous "${tox_ambiguous}" \
        --amr-deeparg "${amr_deeparg}" \
        --amr-amrfinder "${amr_amrfinder}" \
        --amr-rgi "${amr_rgi}" \
        --amr-summary "${amr_summary}" \
        --sample-id "${sample_id}" \
        --outdir .
    """
}