# Analysis plan v1

## Objective and starting evidence
Estimate tumor-associated expression differences in TCGA-LUAD and subsequently investigate exploratory prognostic associations. Live GDC metadata is observed evidence; eligibility, count quality, and clinical completeness are not established. Do not interpret the inventory as an analysis-ready cohort.

## Scientific model
Human bulk RNA-seq; patient is the biological replicate. Primary differential-expression analysis uses patients with one eligible primary-tumor and one eligible solid-tissue-normal sample. Exclude recurrent tumors from this comparison. Audit aliquots, sample duplicates, file versions, and GDC annotations before selecting files; never treat multiple files from one sample as independent replicates. Record every exclusion and its reason. Resolve selection by documented quality and provenance, without looking at differential-expression or survival results.

Use harmonized STAR gene-level raw integer counts, not TPM/FPKM, for count-based testing. Verify the appropriate column against GDC pipeline documentation and actual file headers. Preserve Ensembl IDs and their versions; audit mappings before gene-symbol aggregation. Verify genome and annotation release from downloaded provenance instead of assuming them from the project name.

## Endpoints and methods
1. Paired tumor-minus-normal differential expression: DESeq2 design `~ patient + condition`, with normal as baseline. Check full-rank design and pairing. Patient-level fixed covariates are absorbed by patient effects; do not add collinear age/sex terms. Assess batch identifiability before adjustment. Proposed expression filter: at least 10 counts in at least as many samples as the smaller condition group, finalized before testing. Report log2 fold change, uncertainty, adjusted p-value, and filtering status for all genes. Use BH FDR < 0.05; effect-size thresholds are descriptive and must be labeled. Save package versions, normalization and dispersion diagnostics, library sizes, PCA, and documented exclusions.
2. Pathways: use a ranked statistic and versioned gene sets, documenting ID mapping and the tested gene universe; correct across tested pathways. Do not interpret enrichment as causality.
3. Overall survival: one primary tumor per patient. Reconcile death, last-follow-up, time origin, and days; preserve missingness and exclude unusable endpoints transparently. Evaluate events and stage/age completeness before specifying model complexity. Prefer continuous standardized expression in Cox models with justified clinical covariates; check proportional hazards and nonlinear effects. Report hazard ratios, confidence intervals, events, and multiplicity. Kaplan-Meier grouping is descriptive and must not use an optimized cutoff selected on the same outcome data.
4. Validation: choose an independent LUAD dataset only after checking assay, tissue, sample independence, and clinical endpoint compatibility. Freeze candidate selection, transformations, coefficients if applicable, and cutoffs before evaluation. Any internal split must be by patient, with outcome-informed selection confined to training data. Reusing TCGA for selection and testing supports exploratory associations only.

## Validity gates
- Query response must be complete; preserve query, retrieval time, file UUIDs, and MD5 checksums.
- Unique selected samples, exact matrix-to-metadata alignment, no unexplained duplicates.
- Confirm count type, gene annotation compatibility, pairing, and independent-patient denominators.
- Resolve confounding or leave affected comparisons untested. Document outliers; do not drop samples simply to improve significance.
- No survival claims until endpoints, missingness, event counts, and model assumptions pass review.
- No independent-validation claim without a genuinely independent, compatible cohort.

## Limitations and claim boundaries
Adjacent non-tumor tissue is not healthy-donor tissue. Bulk differences can reflect cell composition, tumor purity, and technical effects. Observational survival associations do not demonstrate causal mechanisms or clinical utility. Negative findings remain valid project outputs. Do not claim novel, validated biomarkers or improved admissions chances from these analyses.

## Implementation requirements and unresolved items
Use R/Bioconductor for count modeling and survival; Python is optional for data preparation. Neither Rscript nor Python was found on the current shell PATH during setup; this does not establish that they are absent from the machine. No software was installed and no analysis workflow was executed. Next: inspect local analysis software and perform the duplicate/aliquot audit before expression downloads. Clinical endpoint coverage, batch metadata, reference release, and external validation cohort remain unknown.
