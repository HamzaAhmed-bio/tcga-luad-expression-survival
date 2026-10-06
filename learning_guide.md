# Understand and present this project

## The question

We asked which genes differ between tumor and non-tumor lung tissue, and whether those genes are associated with survival among patients who already have lung cancer. The first question produced reproducible differences. The second did not produce corrected-significant findings under our models.

## Explain the project in one minute

"This AI-assisted project reanalyzed public TCGA lung adenocarcinoma data. The pipeline audited patient and sample identities, compared gene expression in 50 matched tumor/non-tumor pairs, and selected ten genes using expression results alone. All ten showed the same expression direction in an independent GEO dataset. Survival analyses adjusting for age and stage did not establish corrected-significant prognostic associations. The main lesson is that a strong tumor marker is not automatically a prognosis marker."

Use this description after reviewing and understanding the code and results. Describe your actual role accurately and acknowledge assistance when asked or required.

## Concepts to explain

1. **Patient versus sample versus file:** One patient can supply several tissue samples. A sample may have multiple sequencing files. Counting all files as independent patients exaggerates evidence.
2. **Matched analysis:** Each tumor is compared with non-tumor tissue from the same person. The patient term accounts for stable patient-specific differences.
3. **Counts versus TPM:** Raw read counts support the count-based expression model. Log2(TPM+1) was used as a continuous variable for tumor-only survival analysis.
4. **Log2 fold change:** +1 means twice as high; -1 means half as high, on the model's normalized scale. A large fold change alone does not establish clinical importance.
5. **False discovery rate:** Testing many genes creates opportunities for chance findings. Adjusted p-values address multiplicity. A nominal p-value below 0.05 is insufficient for a ten-gene survival screen.
6. **Hazard ratio:** Here it compares instantaneous death rates per one standard deviation of expression, adjusting for age and stage. It is not an individual survival probability or proof of causality.
7. **External replication:** Test in another cohort without choosing genes or cutoffs to make that cohort look favorable. Different assays and populations limit numerical comparison.
8. **Negative result:** None of the ten genes passed the corrected survival threshold. That does not prove an association is absent; limited events and model assumptions affect detection.

## Questions to practice

- Why exclude ambiguous patient-condition groups?
- Why use raw counts for PyDESeq2 and not GEO array intensities?
- Why can patient age not be added straightforwardly to a fixed-patient paired expression model?
- Why freeze genes before survival testing?
- Why is AGER's TCGA p=0.026 insufficient? Its ten-test q=0.240.
- Why is FABP4's GEO p=0.0059 not validated prognosis? Its q=0.0586, and TCGA did not replicate it.
- What changed in 39 same-plate pairs? All ten retained the same direction and passed exploratory paired sensitivity tests.
- What remains uncertain? Cell composition, purity, treatment, missingness, residual batch effects, external matching, model shape, and low external event counts.

## Application wording

After reviewing the work and being able to explain it:

"Completed an AI-assisted, reproducible analysis of TCGA-LUAD and GSE31210, including sample auditing, paired differential expression, pathway enrichment, and age/stage-adjusted survival testing; documented external expression concordance and negative prognostic results."

Do not claim novel cancer biomarker discovery, a validated clinical predictor, laboratory experiments, or independent authorship of every script. The project is strongest when you can defend its choices and limitations.

## Your next step

Open the report, explain Figures 1 and 5 in your own words, then trace one gene from frozen_candidates.csv to the two survival files. Understanding why the project retained negative results is more useful than memorizing gene names.
