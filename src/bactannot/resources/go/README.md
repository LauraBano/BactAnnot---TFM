# Gene Ontology resource

BactAnnot includes a local copy of the Gene Ontology `go-basic.obo` resource for reproducible and offline GO consensus analysis.

- Resource: `go-basic.obo`
- GO release: `2026-10-13`
- Data version: `releases/2026-10-13`
- Source: Gene Ontology Consortium
- Original resource: `go-basic.obo`
- License: CC BY 4.0
- SHA256: `<SHA256>`

## Purpose

The bundled ontology is used only by the BactAnnot GO comparison layer.

It provides:

- GO identifier resolution;
- alternative-ID normalization;
- safe `replaced_by` normalization for obsolete terms;
- `is_a` relationships;
- `part_of` relationships;
- ancestor traversal;
- reduction to the most specific observed GO terms.

The original annotations produced by Bakta, eggNOG-mapper and InterProScan are never rewritten using this resource.

## Reproducibility

The bundled ontology is intentionally versioned rather than downloaded automatically during normal BactAnnot execution.

This ensures that identical canonical inputs are interpreted against the same GO ontology release even when a newer Gene Ontology version becomes available.

A different ontology may be supplied explicitly in future BactAnnot interfaces when required.

SHA256: 893b32f8705da492df48eb8637e900145a2582e325745012eec1e4c483dd6c72