"""Prepare a task-local PathoFact2 runtime for BactAnnot."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

AMR_ENV_FILES = (
    "amrfinder_plus_env.yaml",
    "rgi_env.yaml",
    "argnorm_env.yaml",
    "DeepARG_env.yaml",
)


def replace_exact_once(
    text: str,
    old: str,
    new: str,
    label: str,
) -> str:
    """Replace one expected PathoFact2 snippet."""

    count = text.count(old)

    if count != 1:
        raise RuntimeError(
            f"Unexpected PathoFact2 {label}: expected 1 compatibility-patch target, found {count}."
        )

    return text.replace(
        old,
        new,
        1,
    )


def patch_vf_classifier(runtime: Path) -> None:
    """Handle non-standard amino acids in the VF classifier."""

    path = runtime / "src" / "vf_prediction2.py"

    text = path.read_text(encoding="utf-8")

    old = "sequence = sequence #.replace('-', '')"
    new = "sequence = sequence.replace('-', '')"

    count = text.count(old)

    if count != 4:
        raise RuntimeError(
            "Unexpected PathoFact2 vf_prediction2.py: "
            "expected 4 compatibility-patch targets, "
            f"found {count}."
        )

    path.write_text(
        text.replace(old, new),
        encoding="utf-8",
    )


def patch_amr_rules(runtime: Path) -> None:
    """Adapt PathoFact2 AMR rules to the validated runtimes."""

    path = runtime / "rules" / "Annotations.smk"

    text = path.read_text(encoding="utf-8")

    text = replace_exact_once(
        text,
        (
            'input: s=dirout+"/Cluster/{s}/samples_rep_seq.fasta", '
            'd=config["DeepARG"]["database"]'
            '+"/database/v2/deeparg.fasta"'
        ),
        ('input: s=dirout+"/Cluster/{s}/samples_rep_seq.fasta"'),
        "DeepARG input",
    )

    text = replace_exact_once(
        text,
        'conda: "deeparg_env"',
        ('conda: config["work_dir"] +"/conda_envs/DeepARG_env.yaml"'),
        "DeepARG Conda environment",
    )

    text = replace_exact_once(
        text,
        (
            "deeparg predict --model LS --type prot "
            "--input {input.s} {params.p} "
            "-d {params.d} --out {params.o}"
        ),
        (
            "deeparg predict --model LS --type prot "
            "--input-file {input.s} {params.p} "
            "--hf-model-path {params.d} "
            "--hf-model-subfolder LS "
            "--diamond-path "
            "{params.d}/bin/linux-x86_64/diamond "
            "--output-file {params.o}"
        ),
        "DeepARG command",
    )

    text = replace_exact_once(
        text,
        'conda: "amrfinder_plus_env"',
        ('conda: config["work_dir"] +"/conda_envs/amrfinder_plus_env.yaml"'),
        "AMRFinderPlus Conda environment",
    )

    text = replace_exact_once(
        text,
        ('params: s="{s}", p=config["rgi"]["params"], o=dirout+"/RGI/{s}"'),
        (
            'params: s="{s}", '
            'p=config["rgi"]["params"], '
            'o=dirout+"/RGI/{s}", '
            'c=config["rgi"]["card_json"]'
        ),
        "RGI parameters",
    )

    text = replace_exact_once(
        text,
        (
            "rgi main -i {output.tmp} "
            "--output_file {params.o} "
            "-t protein -n {threads} "
            "{params.p}  &> {log}"
        ),
        (
            "rgi load -i {params.c} "
            "--local &> {log}\n"
            "              "
            "rgi main -i {output.tmp} "
            "--output_file {params.o} "
            "-t protein -n {threads} "
            "{params.p} --local "
            ">> {log} 2>&1"
        ),
        "RGI command",
    )

    path.write_text(
        text,
        encoding="utf-8",
    )


def patch_amr_reports(runtime: Path) -> None:
    """
    Preserve empty trailing DeepARG fields in the
    PathoFact2 AMR report parser.
    """

    path = runtime / "src" / "reports.py"

    text = path.read_text(encoding="utf-8")

    old = (
        '            if not line.startswith("#"):\n'
        "                line=line.rstrip()\n"
        '                line=line.split("\\t")\n'
        "                prob=str(round(float(line[6]),2))\n"
    )

    new = (
        '            if not line.startswith("#"):\n'
        '                line=line.rstrip("\\n")\n'
        '                line=line.split("\\t")\n'
        "                prob=str(round(float(line[6]),2))\n"
    )

    text = replace_exact_once(
        text,
        old,
        new,
        "DeepARG report parser",
    )

    path.write_text(
        text,
        encoding="utf-8",
    )


def install_conda_envs(
    runtime: Path,
    assets: Path,
) -> None:
    """Install BactAnnot-compatible PathoFact2 environments."""

    source_dir = assets / "conda_envs"
    target_dir = runtime / "conda_envs"

    for filename in AMR_ENV_FILES:
        source = source_dir / filename
        target = target_dir / filename

        if not source.is_file():
            raise FileNotFoundError(f"Missing PathoFact asset: {source}")

        shutil.copy2(
            source,
            target,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Prepare a task-local PathoFact2 runtime for BactAnnot.")
    )

    parser.add_argument(
        "--runtime",
        type=Path,
        required=True,
        help="Task-local PathoFact2 runtime directory.",
    )

    parser.add_argument(
        "--assets",
        type=Path,
        required=True,
        help="BactAnnot PathoFact asset directory.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    runtime = args.runtime.resolve()
    assets = args.assets.resolve()

    if not runtime.is_dir():
        raise NotADirectoryError(f"PathoFact runtime not found: {runtime}")

    if not assets.is_dir():
        raise NotADirectoryError(f"PathoFact assets not found: {assets}")

    patch_vf_classifier(runtime)
    patch_amr_rules(runtime)
    patch_amr_reports(runtime)

    install_conda_envs(
        runtime,
        assets,
    )

    print("Prepared BactAnnot-compatible PathoFact2 task-local runtime.")


if __name__ == "__main__":
    main()
