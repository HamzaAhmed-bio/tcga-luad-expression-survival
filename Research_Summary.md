# Lung adenocarcinoma: reproducible expression differences and survival follow-up

A public-data reanalysis of TCGA-LUAD compared tumor and non-tumor tissue in 50 matched patients. Of 19,026 expression-filtered genes, 13,821 passed FDR <0.05 and convergence checks; 4,528 also had an absolute log2 fold change of at least 1. Ten expression-selected genes showed the same direction and passed ten-test FDR <0.05 in the independent GSE31210 microarray cohort. Age- and stage-adjusted overall-survival analyses included 465 TCGA patients (169 deaths) and 204 GEO patients (30 deaths). No candidate passed survival FDR <0.05 in either cohort. These findings support reproducible tissue-expression differences, but do not establish prognostic biomarkers.

## Main findings

- 13,821 genes passed FDR <0.05 and convergence checks; 4,528 also met |log2FC| >=1.
- All ten frozen expression candidates agreed in direction in GSE31210 and passed exploratory external expression FDR <0.05.
- All ten passed paired sensitivity tests in both 50 total pairs and 39 same-plate pairs.
- No corrected-significant survival association was established in TCGA or GEO.
- 24 of 50 Hallmark sets had enrichment FDR <0.05.

## Interpretation

The project supports reproducible tissue-expression differences. It does not establish prognostic biomarkers, causality, clinical utility, or novel discoveries. External expression matching is unresolved, survival diagnostics identify limitations, and the GEO survival cohort has only 30 deaths. See Research_Report.html for complete methods, plots, diagnostics, and limitations.

## Artifact map

- Results: results/ contains the full differential-expression table, paired raw matrix, frozen candidates, survival results, and pathway results.
- QC and visuals: results/figures/ contains PNG and SVG charts; sample_qc.csv and validation_receipt.json record checks.
- Provenance: metadata/, requirements-lock.txt, docs/analysis_plan_v2.md, and docs/decision_log.md preserve selection and implementation details.
- Interpretation: this summary, Research_Report.html, and docs/learning_guide.md explain supported findings and claim boundaries.

AI-assisted code and analysis prepared for Hamza Ahmed; review and understanding are needed before independent presentation.
