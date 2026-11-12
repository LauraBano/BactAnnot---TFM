"""Integration of canonical BactAnnot evidence tables"""

from __future__ import annotations

import csv
from pathlib import Path

from bactannot.schema import EVIDENCE_COLUMNS

# Define las columnas mínimas de features.tsv que necesitamos para validar las evidencias
REQUIRED_FEATURE_COLUMNS = {
    "sample_id",
    "feature_uid",
    "protein_id",
}


# Carga la relación feature_uid -> protein_id definida por Bakta para esta muestra
def load_feature_to_protein_index(
    features_path: Path,
    sample_id: str,
) -> dict[str, str]:
    """Map canonical feature UIDs to their protein identifiers"""

    feature_index: dict[str, str] = {}

    with features_path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(f"Missing header in features table: {features_path}")

        missing = REQUIRED_FEATURE_COLUMNS.difference(reader.fieldnames)

        if missing:
            raise ValueError(
                f"Missing required columns in {features_path}: {', '.join(sorted(missing))}"
            )

        for row in reader:
            row_sample_id = row["sample_id"].strip()
            feature_uid = row["feature_uid"].strip()
            protein_id = row["protein_id"].strip()

            # Todas las features deben pertenecer a la muestra que estamos integrando
            if row_sample_id != sample_id:
                raise ValueError(
                    f"Unexpected sample_id in {features_path}: {row_sample_id} != {sample_id}"
                )

            if not feature_uid:
                raise ValueError(f"Missing feature_uid in {features_path}")

            # Cada feature_uid debe identificar una única feature
            if feature_uid in feature_index:
                raise ValueError(f"Duplicate feature_uid in features table: {feature_uid}")

            feature_index[feature_uid] = protein_id

    return feature_index


# Lee una tabla de evidence y comprueba que siga correctamente el modelo canónico
def read_evidence_table(
    path: Path,
    sample_id: str,
    feature_index: dict[str, str],
) -> list[dict[str, str]]:
    """Read and validate one canonical evidence table"""

    rows: list[dict[str, str]] = []

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # Exigimos exactamente el mismo esquema para poder integrar evidencias de herramientas distintas
        if reader.fieldnames != EVIDENCE_COLUMNS:
            raise ValueError(f"Unexpected evidence schema in {path}: {reader.fieldnames}")

        for line_number, row in enumerate(
            reader,
            start=2,
        ):
            # Convierte cualquier valor ausente en cadena vacía y mantiene solo las columnas canónicas
            canonical_row = {column: row.get(column) or "" for column in EVIDENCE_COLUMNS}

            if canonical_row["sample_id"] != sample_id:
                raise ValueError(
                    f"Unexpected sample_id in {path} "
                    f"at line {line_number}: "
                    f"{canonical_row['sample_id']} != "
                    f"{sample_id}"
                )

            feature_uid = canonical_row["feature_uid"]

            if not feature_uid.strip():
                raise ValueError(f"Missing feature_uid in {path} at line {line_number}")

            # La evidencia solo es válida si apunta a una feature conocida por Bakta
            if feature_uid not in feature_index:
                raise ValueError(
                    f"Unknown feature_uid in {path} at line {line_number}: {feature_uid}"
                )

            expected_protein_id = feature_index[feature_uid]

            # El protein_id debe coincidir con el asignado a esa feature en features.tsv
            if canonical_row["protein_id"] != expected_protein_id:
                raise ValueError(
                    f"protein_id mismatch in {path} "
                    f"at line {line_number}: "
                    f"{canonical_row['protein_id']} != "
                    f"{expected_protein_id}"
                )

            # Toda evidencia necesita indicar qué tipo es, qué fuente la produjo y qué referencia encontró
            for column in (
                "evidence_type",
                "source",
                "accession",
            ):
                if not canonical_row[column].strip():
                    raise ValueError(f"Missing {column} in {path} at line {line_number}")

            rows.append(canonical_row)

    return rows


# Reúne las tablas de evidence de todas las herramientas en una única tabla canónica
def integrate_evidence(
    evidence_paths: list[Path],
    features_path: Path,
    sample_id: str,
) -> list[dict[str, str]]:
    """Integrate canonical evidence tables without losing provenance"""

    if not evidence_paths:
        raise ValueError("At least one evidence table is required")

    # Evita integrar accidentalmente dos veces el mismo archivo
    if len(evidence_paths) != len(set(evidence_paths)):
        raise ValueError("The same evidence table was provided more than once")

    # features.tsv actúa como referencia para validar feature_uid y protein_id
    feature_index = load_feature_to_protein_index(
        features_path,
        sample_id,
    )

    integrated: list[dict[str, str]] = []

    # En evidence consideramos duplicada únicamente una fila canónica completamente idéntica
    seen: set[tuple[str, ...]] = set()

    for path in evidence_paths:
        # Cada tabla se valida por separado antes de incorporarla al conjunto integrado
        rows = read_evidence_table(
            path,
            sample_id,
            feature_index,
        )

        for row in rows:
            # Todos los campos forman parte de la identidad de una fila de evidencia
            key = tuple(row[column] for column in EVIDENCE_COLUMNS)

            if key in seen:
                raise ValueError(f"Duplicate canonical evidence row: {key}")

            seen.add(key)
            integrated.append(row)

    # Ordena el resultado completo para obtener una salida estable entre ejecuciones
    integrated.sort(key=lambda row: tuple(row[column] for column in EVIDENCE_COLUMNS))

    return integrated


# Escribe la tabla integrada respetando el esquema canónico de evidence
def write_evidence_table(
    rows: list[dict[str, str]],
    output_path: Path,
) -> None:
    """Write an integrated canonical evidence table"""

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=EVIDENCE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


# Ejecuta la integración completa y escribe el resultado final
def integrate_evidence_tables(
    evidence_paths: list[Path],
    features_path: Path,
    sample_id: str,
    output_path: Path,
) -> None:
    """Integrate and write canonical evidence tables"""

    rows = integrate_evidence(
        evidence_paths=evidence_paths,
        features_path=features_path,
        sample_id=sample_id,
    )

    write_evidence_table(
        rows=rows,
        output_path=output_path,
    )
