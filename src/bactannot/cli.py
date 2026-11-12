"""Command-line interface for BactAnnot."""

from __future__ import annotations

import argparse
from pathlib import Path

from bactannot import __version__
from bactannot.concordance.concordance import build_concordance_tables
from bactannot.consensus.consensus import build_evidence_consensus_tables
from bactannot.integration.annotations import integrate_annotation_tables
from bactannot.integration.evidence import integrate_evidence_tables
from bactannot.parsers.bakta import parse_bakta_bundle
from bactannot.parsers.eggnog import parse_eggnog_bundle
from bactannot.parsers.interproscan import parse_interproscan_bundle
from bactannot.parsers.kofamscan import parse_kofamscan_bundle
from bactannot.parsers.pathofact import parse_pathofact_bundle


# Construye la interfaz principal y define todos los subcomandos disponibles
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bactannot",
        description=("Integrate and analyze bacterial functional annotations."),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    # Crea un subcomando independiente para cada etapa que puede ejecutar BactAnnot
    subparsers = parser.add_subparsers(dest="command")

    # Subcomando para normalizar la salida de Bakta
    parse_bakta = subparsers.add_parser(
        "parse-bakta",
        help="Normalize a Bakta JSON annotation result.",
    )

    parse_bakta.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Bakta JSON input file.",
    )

    parse_bakta.add_argument(
        "--sample-id",
        required=True,
        help="BactAnnot sample identifier.",
    )

    parse_bakta.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    # Subcomando para normalizar las anotaciones y la evidencia de eggNOG
    parse_eggnog = subparsers.add_parser(
        "parse-eggnog",
        help=("Normalize eggNOG-mapper functional annotations."),
    )

    parse_eggnog.add_argument(
        "--annotations",
        type=Path,
        required=True,
        help="eggNOG-mapper annotations file.",
    )

    parse_eggnog.add_argument(
        "--seed-orthologs",
        type=Path,
        required=True,
        help="eggNOG-mapper seed orthologs file.",
    )

    parse_eggnog.add_argument(
        "--features",
        type=Path,
        required=True,
        help="Canonical BactAnnot features TSV.",
    )

    parse_eggnog.add_argument(
        "--sample-id",
        required=True,
        help="BactAnnot sample identifier.",
    )

    parse_eggnog.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    # Subcomando para normalizar las anotaciones y la evidencia de InterProScan
    parse_interproscan = subparsers.add_parser(
        "parse-interproscan",
        help=("Normalize InterProScan protein annotation results."),
    )

    parse_interproscan.add_argument(
        "--input",
        type=Path,
        required=True,
        help="InterProScan TSV input file.",
    )

    parse_interproscan.add_argument(
        "--features",
        type=Path,
        required=True,
        help="Canonical BactAnnot features TSV.",
    )

    parse_interproscan.add_argument(
        "--sample-id",
        required=True,
        help="BactAnnot sample identifier.",
    )

    parse_interproscan.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    # Subcomando para normalizar los KO y la evidencia HMM de KofamScan
    parse_kofamscan = subparsers.add_parser(
        "parse-kofamscan",
        help="Normalize KofamScan KEGG Orthology assignments.",
    )

    parse_kofamscan.add_argument(
        "--input",
        type=Path,
        required=True,
        help="KofamScan detail-tsv input file.",
    )

    parse_kofamscan.add_argument(
        "--features",
        type=Path,
        required=True,
        help="Canonical BactAnnot features TSV.",
    )

    parse_kofamscan.add_argument(
        "--sample-id",
        required=True,
        help="BactAnnot sample identifier.",
    )

    parse_kofamscan.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    # Subcomando para normalizar las salidas de virulencia, toxinas y AMR de PathoFact
    parse_pathofact = subparsers.add_parser(
        "parse-pathofact",
        help=("Normalize PathoFact2 virulence, toxin and AMR annotations."),
    )

    parse_pathofact.add_argument(
        "--features",
        type=Path,
        required=True,
        help=("Canonical BactAnnot features TSV generated from Bakta."),
    )

    parse_pathofact.add_argument(
        "--clusters",
        type=Path,
        required=True,
        help="PathoFact2 clean_cluster.tsv file.",
    )

    # Archivos de predicción de factores de virulencia
    parse_pathofact.add_argument(
        "--vf-ml",
        type=Path,
        required=True,
        help=("PathoFact2 virulence-factor ML prediction TSV."),
    )

    parse_pathofact.add_argument(
        "--vf-hmm",
        type=Path,
        required=True,
        help=("PathoFact2 virulence-factor HMM results."),
    )

    parse_pathofact.add_argument(
        "--vf-high",
        type=Path,
        required=True,
        help=("PathoFact2 high-confidence VF summary TSV."),
    )

    parse_pathofact.add_argument(
        "--vf-ambiguous",
        type=Path,
        required=True,
        help=("PathoFact2 ambiguous VF summary TSV."),
    )

    # Archivos de predicción de toxinas
    parse_pathofact.add_argument(
        "--tox-ml",
        type=Path,
        required=True,
        help=("PathoFact2 toxin ML prediction TSV."),
    )

    parse_pathofact.add_argument(
        "--tox-hmm",
        type=Path,
        required=True,
        help="PathoFact2 toxin HMM results.",
    )

    parse_pathofact.add_argument(
        "--tox-high",
        type=Path,
        required=True,
        help=("PathoFact2 high-confidence toxin summary TSV."),
    )

    parse_pathofact.add_argument(
        "--tox-ambiguous",
        type=Path,
        required=True,
        help=("PathoFact2 ambiguous toxin summary TSV."),
    )

    # Resultados de los métodos AMR y resumen integrado de PathoFact
    parse_pathofact.add_argument(
        "--amr-deeparg",
        type=Path,
        required=True,
        help=("PathoFact2 DeepARG argNorm-normalized TSV."),
    )

    parse_pathofact.add_argument(
        "--amr-amrfinder",
        type=Path,
        required=True,
        help=("PathoFact2 AMRFinderPlus argNorm-normalized TSV."),
    )

    parse_pathofact.add_argument(
        "--amr-rgi",
        type=Path,
        required=True,
        help=("PathoFact2 RGI/CARD TSV."),
    )

    parse_pathofact.add_argument(
        "--amr-summary",
        type=Path,
        required=True,
        help=("PathoFact2 integrated AMR summary TSV."),
    )

    parse_pathofact.add_argument(
        "--sample-id",
        required=True,
        help="BactAnnot sample identifier.",
    )

    parse_pathofact.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    # Subcomando para reunir todas las tablas canónicas de annotations
    integrate_annotations = subparsers.add_parser(
        "integrate-annotations",
        help="Integrate canonical annotation tables.",
    )

    integrate_annotations.add_argument(
        "--features",
        type=Path,
        required=True,
        help="Canonical BactAnnot features TSV.",
    )

    integrate_annotations.add_argument(
        "--annotations",
        type=Path,
        nargs="+",
        required=True,
        help="Canonical annotation TSV files to integrate.",
    )

    integrate_annotations.add_argument(
        "--sample-id",
        required=True,
        help="BactAnnot sample identifier.",
    )

    integrate_annotations.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    # Subcomando para reunir todas las tablas canónicas de evidence
    integrate_evidence = subparsers.add_parser(
        "integrate-evidence",
        help="Integrate canonical evidence tables.",
    )

    integrate_evidence.add_argument(
        "--features",
        type=Path,
        required=True,
        help="Canonical BactAnnot features TSV.",
    )

    integrate_evidence.add_argument(
        "--evidence",
        type=Path,
        nargs="+",
        required=True,
        help="Canonical evidence TSV files to integrate.",
    )

    integrate_evidence.add_argument(
        "--sample-id",
        required=True,
        help="BactAnnot sample identifier.",
    )

    integrate_evidence.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    # Subcomando para construir el consenso funcional a partir de las tablas integradas
    build_consensus = subparsers.add_parser(
        "build-consensus",
        help="Build functional claims and relation tables.",
    )

    build_consensus.add_argument(
        "--features",
        type=Path,
        required=True,
        help="Canonical BactAnnot features TSV.",
    )

    build_consensus.add_argument(
        "--annotations",
        type=Path,
        required=True,
        help="Integrated canonical annotation TSV.",
    )

    build_consensus.add_argument(
        "--evidence",
        type=Path,
        required=True,
        help="Integrated canonical evidence TSV.",
    )

    build_consensus.add_argument(
        "--sample-id",
        required=True,
        help="BactAnnot sample identifier.",
    )

    build_consensus.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    build_consensus.add_argument(
        "--ko-list",
        type=Path,
        required=True,
        help="Local KOfam ko_list reference catalog.",
    )

    # Subcomando para comparar por pares las anotaciones de las distintas fuentes
    build_concordance = subparsers.add_parser(
        "build-concordance",
        help="Build pairwise functional concordance diagnostics.",
    )

    build_concordance.add_argument(
        "--features",
        type=Path,
        required=True,
        help="Canonical BactAnnot features TSV.",
    )

    build_concordance.add_argument(
        "--annotations",
        type=Path,
        required=True,
        help="Integrated canonical annotation TSV.",
    )

    build_concordance.add_argument(
        "--sample-id",
        required=True,
        help="BactAnnot sample identifier.",
    )

    build_concordance.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory.",
    )

    return parser


# Prepara las salidas y ejecuta el parser completo de Bakta
def run_parse_bakta(
    args: argparse.Namespace,
) -> None:
    # Crea el directorio de salida si todavía no existe
    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    features_output = args.outdir / f"{args.sample_id}.bakta.features.tsv"

    annotations_output = args.outdir / f"{args.sample_id}.bakta.annotations.tsv"

    evidence_output = args.outdir / f"{args.sample_id}.bakta.evidence.tsv"

    # Delega la normalización real al parser de Bakta
    parse_bakta_bundle(
        input_path=args.input,
        sample_id=args.sample_id,
        features_output=features_output,
        annotations_output=annotations_output,
        evidence_output=evidence_output,
    )

    print(f"Features:    {features_output}")
    print(f"Annotations: {annotations_output}")
    print(f"Evidence:    {evidence_output}")


# Prepara las salidas y ejecuta el parser completo de eggNOG
def run_parse_eggnog(
    args: argparse.Namespace,
) -> None:
    # Crea el directorio de salida si todavía no existe
    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    annotations_output = args.outdir / f"{args.sample_id}.eggnog.annotations.tsv"

    evidence_output = args.outdir / f"{args.sample_id}.eggnog.evidence.tsv"

    # Delega la normalización real al parser de eggNOG
    parse_eggnog_bundle(
        annotations_path=args.annotations,
        seed_orthologs_path=args.seed_orthologs,
        features_path=args.features,
        sample_id=args.sample_id,
        annotations_output=annotations_output,
        evidence_output=evidence_output,
    )

    print(f"Annotations: {annotations_output}")

    print(f"Evidence: {evidence_output}")


# Prepara las salidas y ejecuta el parser completo de InterProScan
def run_parse_interproscan(
    args: argparse.Namespace,
) -> None:
    # Crea el directorio de salida si todavía no existe
    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    annotations_output = args.outdir / (f"{args.sample_id}.interproscan.annotations.tsv")

    evidence_output = args.outdir / (f"{args.sample_id}.interproscan.evidence.tsv")

    # Delega la normalización real al parser de InterProScan
    parse_interproscan_bundle(
        input_path=args.input,
        features_path=args.features,
        sample_id=args.sample_id,
        annotations_output=annotations_output,
        evidence_output=evidence_output,
    )

    print(f"Annotations: {annotations_output}")

    print(f"Evidence: {evidence_output}")


# Prepara las salidas y ejecuta el parser completo de KofamScan
def run_parse_kofamscan(
    args: argparse.Namespace,
) -> None:
    # Crea el directorio de salida si todavía no existe
    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    annotations_output = args.outdir / f"{args.sample_id}.kofamscan.annotations.tsv"

    evidence_output = args.outdir / f"{args.sample_id}.kofamscan.evidence.tsv"

    # Delega la normalización real al parser de KofamScan
    parse_kofamscan_bundle(
        input_path=args.input,
        features_path=args.features,
        sample_id=args.sample_id,
        annotations_output=annotations_output,
        evidence_output=evidence_output,
    )

    print(f"Annotations: {annotations_output}")

    print(f"Evidence: {evidence_output}")


# Prepara las salidas y ejecuta el parser completo de PathoFact
def run_parse_pathofact(
    args: argparse.Namespace,
) -> None:
    # Crea el directorio de salida si todavía no existe
    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    annotations_output = args.outdir / (f"{args.sample_id}.pathofact.annotations.tsv")

    evidence_output = args.outdir / (f"{args.sample_id}.pathofact.evidence.tsv")

    # Delega la normalización real al parser de PathoFact
    parse_pathofact_bundle(
        features_path=args.features,
        cluster_path=args.clusters,
        vf_ml_path=args.vf_ml,
        vf_hmm_path=args.vf_hmm,
        vf_high_path=args.vf_high,
        vf_ambiguous_path=args.vf_ambiguous,
        tox_ml_path=args.tox_ml,
        tox_hmm_path=args.tox_hmm,
        tox_high_path=args.tox_high,
        tox_ambiguous_path=args.tox_ambiguous,
        amr_deeparg_path=args.amr_deeparg,
        amr_amrfinder_path=args.amr_amrfinder,
        amr_rgi_path=args.amr_rgi,
        amr_summary_path=args.amr_summary,
        sample_id=args.sample_id,
        annotations_output=annotations_output,
        evidence_output=evidence_output,
    )

    print(f"Annotations: {annotations_output}")

    print(f"Evidence: {evidence_output}")


# Integra las tablas de annotations y escribe una única salida canónica
def run_integrate_annotations(
    args: argparse.Namespace,
) -> None:
    # Crea el directorio de salida si todavía no existe
    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    annotations_output = args.outdir / f"{args.sample_id}.integrated.annotations.tsv"

    # Delega la integración a la capa específica de annotations
    integrate_annotation_tables(
        annotation_paths=args.annotations,
        features_path=args.features,
        sample_id=args.sample_id,
        output_path=annotations_output,
    )

    print(f"Annotations: {annotations_output}")


# Integra las tablas de evidence y escribe una única salida canónica
def run_integrate_evidence(
    args: argparse.Namespace,
) -> None:
    # Crea el directorio de salida si todavía no existe
    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    evidence_output = args.outdir / f"{args.sample_id}.integrated.evidence.tsv"

    # Delega la integración a la capa específica de evidence
    integrate_evidence_tables(
        evidence_paths=args.evidence,
        features_path=args.features,
        sample_id=args.sample_id,
        output_path=evidence_output,
    )

    print(f"Evidence: {evidence_output}")


# Construye y escribe la tabla final de consensus
def run_build_consensus(
    args: argparse.Namespace,
) -> None:
    """Build functional claim and relation tables."""

    outputs = build_evidence_consensus_tables(
        annotations_path=args.annotations,
        evidence_path=args.evidence,
        features_path=args.features,
        sample_id=args.sample_id,
        ko_list_path=args.ko_list,
        output_dir=args.outdir,
    )

    for name, path in outputs.items():
        print(f"{name}: {path}")


# Construye y escribe las tablas detallada y resumen de concordance
def run_build_concordance(
    args: argparse.Namespace,
) -> None:
    """Build pairwise functional concordance diagnostics."""

    # Crea el directorio de salida si todavía no existe
    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    concordance_output = args.outdir / f"{args.sample_id}.concordance.tsv"

    summary_output = args.outdir / f"{args.sample_id}.concordance_summary.tsv"

    # Delega la comparación entre fuentes a concordance.py
    build_concordance_tables(
        annotations_path=args.annotations,
        features_path=args.features,
        sample_id=args.sample_id,
        concordance_output_path=concordance_output,
        summary_output_path=summary_output,
    )

    print(f"Concordance: {concordance_output}")
    print(f"Summary:     {summary_output}")


# Lee el comando solicitado y dirige la ejecución a la función correspondiente
def main() -> None:
    # Construye la CLI y convierte los argumentos escritos por el usuario en args
    parser = build_parser()
    args = parser.parse_args()

    # Cada subcomando se dirige a su función run_* correspondiente

    if args.command == "parse-bakta":
        run_parse_bakta(args)
        return

    if args.command == "parse-eggnog":
        run_parse_eggnog(args)
        return

    if args.command == "parse-interproscan":
        run_parse_interproscan(args)
        return

    if args.command == "parse-kofamscan":
        run_parse_kofamscan(args)
        return

    if args.command == "parse-pathofact":
        run_parse_pathofact(args)
        return

    if args.command == "integrate-annotations":
        run_integrate_annotations(args)
        return

    if args.command == "integrate-evidence":
        run_integrate_evidence(args)
        return

    if args.command == "build-consensus":
        run_build_consensus(args)
        return

    if args.command == "build-concordance":
        run_build_concordance(args)
        return

    parser.print_help()


# Permite ejecutar esta CLI directamente como módulo principal
if __name__ == "__main__":
    main()
