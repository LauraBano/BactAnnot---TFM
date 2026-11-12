"""Minimal Gene Ontology support for GO consensus comparisons"""

from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path


@dataclass
class GoTerm:
    """Minimal GO term representation"""

    go_id: str
    namespace: str
    alt_ids: set[str] = field(default_factory=set)
    parents: set[str] = field(default_factory=set)
    is_obsolete: bool = False
    replaced_by: list[str] = field(default_factory=list)
    ec_xrefs: set[str] = field(default_factory=set)
    is_a_parents: set[str] = field(default_factory=set)
    name: str = ""


@dataclass(frozen=True)
class GoResolution:
    """Resolved GO identifier for ontology-aware comparison"""

    original_id: str
    resolved_id: str
    namespace: str
    status: str
    hierarchy_eligible: bool


def get_bundled_go_ontology() -> Path:
    """Return the bundled go-basic.obo path"""

    resource = files("bactannot").joinpath("resources").joinpath("go").joinpath("go-basic.obo")

    return Path(str(resource))


def load_go_ontology(
    path: Path | None = None,
) -> tuple[
    dict[str, GoTerm],
    dict[str, str],
]:
    """Load the GO fields required for BactAnnot comparisons."""

    if path is None:
        path = get_bundled_go_ontology()

    terms: dict[
        str,
        GoTerm,
    ] = {}

    alt_to_primary: dict[
        str,
        str,
    ] = {}

    current: GoTerm | None = None

    def store_current() -> None:
        if current is None:
            return

        if not current.go_id:
            return

        terms[current.go_id] = current

        for alt_id in current.alt_ids:
            alt_to_primary[alt_id] = current.go_id

    with path.open(
        encoding="utf-8",
    ) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if line == "[Term]":
                store_current()

                current = GoTerm(
                    go_id="",
                    namespace="",
                )

                continue

            if line.startswith("[") and line != "[Term]":
                store_current()
                current = None
                continue

            if current is None:
                continue

            if line.startswith("id: "):
                current.go_id = line[len("id: ") :].strip()

            elif line.startswith("namespace: "):
                current.namespace = line[len("namespace: ") :].strip()

            elif line.startswith("name: "):
                current.name = line[len("name: ") :].strip()

            elif line.startswith("alt_id: "):
                current.alt_ids.add(line[len("alt_id: ") :].strip())

            elif line.startswith("xref: EC:"):
                ec_id = line[len("xref: ") :].split(maxsplit=1)[0]

                current.ec_xrefs.add(ec_id)

            elif line.startswith("is_a: "):
                parent = (
                    line[len("is_a: ") :]
                    .split(
                        " ! ",
                        1,
                    )[0]
                    .strip()
                )

                current.parents.add(parent)
                current.is_a_parents.add(parent)

            elif line.startswith("relationship: part_of "):
                parent = (
                    line[len("relationship: part_of ") :]
                    .split(
                        " ! ",
                        1,
                    )[0]
                    .strip()
                )

                current.parents.add(parent)

            elif line == "is_obsolete: true":
                current.is_obsolete = True

            elif line.startswith("replaced_by: "):
                current.replaced_by.append(line[len("replaced_by: ") :].strip())

    store_current()

    return (
        terms,
        alt_to_primary,
    )


def resolve_go_id(
    go_id: str,
    terms: dict[str, GoTerm],
    alt_to_primary: dict[str, str],
) -> GoResolution:
    """Resolve a GO identifier conservatively."""

    was_alt = False

    if go_id in terms:
        primary_id = go_id

    elif go_id in alt_to_primary:
        primary_id = alt_to_primary[go_id]
        was_alt = True

    else:
        raise ValueError(f"GO identifier is absent from the ontology: {go_id}")

    term = terms[primary_id]

    if not term.is_obsolete:
        return GoResolution(
            original_id=go_id,
            resolved_id=primary_id,
            namespace=term.namespace,
            status=("alt_id" if was_alt else "active"),
            hierarchy_eligible=True,
        )

    if len(term.replaced_by) == 1:
        replacement_id = term.replaced_by[0]

        replacement = terms.get(replacement_id)

        if replacement is not None and not replacement.is_obsolete:
            return GoResolution(
                original_id=go_id,
                resolved_id=replacement_id,
                namespace=(replacement.namespace),
                status=("alt_replaced_by" if was_alt else "replaced_by"),
                hierarchy_eligible=True,
            )

    return GoResolution(
        original_id=go_id,
        resolved_id=primary_id,
        namespace=term.namespace,
        status="obsolete_unresolved",
        hierarchy_eligible=False,
    )


def build_go_ancestor_map(
    terms: dict[str, GoTerm],
    *,
    is_a_only: bool = False,
) -> dict[str, frozenset[str]]:
    """Build GO ancestry, optionally restricted to is_a relations"""

    cache: dict[
        str,
        frozenset[str],
    ] = {}

    visiting: set[str] = set()

    def get_ancestors(
        go_id: str,
    ) -> frozenset[str]:
        if go_id in cache:
            return cache[go_id]

        if go_id in visiting:
            raise ValueError(f"Cycle detected in GO hierarchy at {go_id}")

        visiting.add(go_id)

        result: set[str] = set()

        term = terms[go_id]

        if is_a_only:
            parent_ids = term.is_a_parents
        else:
            parent_ids = term.parents

        for parent_id in parent_ids:
            parent = terms.get(parent_id)

            if parent is None:
                continue

            if parent.is_obsolete:
                continue

            result.add(parent_id)

            result.update(get_ancestors(parent_id))

        visiting.remove(go_id)

        resolved = frozenset(result)

        cache[go_id] = resolved

        return resolved

    for go_id, term in terms.items():
        if term.is_obsolete:
            continue

        get_ancestors(go_id)

    return cache


def reduce_to_most_specific(
    go_ids: set[str],
    ancestor_map: dict[
        str,
        frozenset[str],
    ],
) -> set[str]:
    """Remove GO terms implied by more specific observed terms."""

    redundant: set[str] = set()

    for descendant_id in go_ids:
        redundant.update(
            go_ids.intersection(
                ancestor_map.get(
                    descendant_id,
                    frozenset(),
                )
            )
        )

    return go_ids - redundant
