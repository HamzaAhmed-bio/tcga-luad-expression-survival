# Data and methods provenance

Data retrieved October 5, 2026 (America/Chicago); UTC timestamps may read October 6. The data are public and deidentified.

| Source | Use | Snapshot |
|---|---|---|
| [GDC TCGA-LUAD](https://portal.gdc.cancer.gov/projects/TCGA-LUAD) | Harmonized gene counts and sample/aliquot identity | Query and complete JSON response in metadata; 601 files audited; 546 selected and MD5-verified |
| [GDC mRNA pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/) | Count interpretation | GRCh38, unstranded counts; every downloaded file header states GENCODE v36 |
| [Liu et al., TCGA-CDR, Cell 2018](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018) | Curated overall survival, age, stage | Official Supplemental Table S1, GDC file 1b5f413e-a8d1-4d10-92eb-7c4ae739ed81; historical clinical freeze |
| [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210) | Independent expression and survival follow-up | 226 primary lung adenocarcinomas plus 20 normal tissue samples; source-normalized MAS5 microarrays |
| [GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) | Probe-to-gene annotation | GEO annotation dated August 9, 2016; ambiguous multi-gene mappings excluded |
| [PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/generated/pydeseq2.dds.DeseqDataSet.html) | Paired negative-binomial differential expression | 0.5.4; Python implementation, not the R package |
| [lifelines](https://lifelines.readthedocs.io/en/latest/fitters/regression/CoxPHFitter.html) | Cox regression and Kaplan-Meier plots | Exact version in requirements-lock.txt |
| [GSEApy](https://gseapy.readthedocs.io/en/latest/gseapy_example.html) | Preranked gene-set enrichment | Hallmark 2020 GMT from Enrichr; saved content hash; 1,000 permutations and fixed seed |

The original GSE31210 publications are identified in its series metadata by PubMed IDs 22080568 and 23028479. Refer to the source study for recruitment and clinical endpoint definitions. The independent cohort concerns pathological stage I-II disease and differs from TCGA in stage distribution, assay, geography, and treatment context; replication is not proof of a transportable clinical predictor.

Current GDC clinical metadata were inspected but not used to reconstruct survival because follow-up availability was incomplete. An attempted Xena survival URL returned AccessDenied and was not used. The official TCGA-CDR source above was used instead. No controlled-access data, FASTQ processing, clinical deployment, or Workbench registered run was involved.
