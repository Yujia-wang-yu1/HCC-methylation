# HCC MISP scoring and analysis

- [Supplementary tables](data/Results.xlsx) contains Tables S1–S33. The Index sheet locates each table. Table S2 includes the HM450 probe counts, and Table S3 contains the fixed 30-gene MISP weights.
- [Figure source data](data/Figure_data.xlsx) supports the [plotting scripts](code/).
- [MISP scoring](code/score_misp.py) applies the fixed weights to a new cohort. [Study analysis](analysis/run.py) runs the analyses listed below from prepared inputs.

## Software

Use Python 3.13. The analysis was checked with NumPy 2.5.1, pandas 3.0.5 and SciPy 1.18.0. Scoring and plotting use openpyxl 3.1.5 to read the workbooks. Plotting uses Matplotlib 3.11.1.

The R stages use R 4.6.1. TCGA composition uses MASS 7.3.65 and EpiDISH 2.28.0 or its corresponding source helper. ICGC survival uses survival 3.8.6.

## Plot a figure or calculate MISP

From the repository root, run a figure script or one of the scoring commands below.

```sh
python code/figure2.py --output figures
python code/score_misp.py --assay methylation --input promoter_beta.tsv --output methylation_scores.tsv
python code/score_misp.py --assay rna --input log_expression.tsv --output rna_suppression_scores.tsv
```

For scoring, use one tab-separated matrix per assay and cohort, with unique gene symbols in the first column and samples in the remaining columns. Supply gene-level TSS200 promoter beta values or normalized, transformed gene-level RNA expression as specified in Methods 2.2 and Supplementary Methods S3. Pass `--standardized` for inputs already expressed as gene-wise z scores from the intended calibration cohort.

The scorer applies gene-wise cohort standardization and reads the weights from Table S3. The RNA suppression score reverses the sign of weighted expression. The full weight denominator is used for partially observed samples. A sample with no available contributions receives a missing score. To score an analysis subset, score the full calibration cohort first and then select the subset.

## Reproduce study analyses

The analysis starts from prepared gene profiles and covariates. The cohort sources are TCGA-LIHC (GDC), [GEO GSE56588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56588) and [GSE63898](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63898), and ICGC LIRI-JP (legacy release 28). The selection uses the frozen GO Biological Process, MSigDB Hallmark, Reactome and KEGG snapshots described in the manuscript. Methods 2.2–2.6 and Supplementary Methods S1–S3 and S8–S12 document processing and selection.

Request the prepared files from the corresponding author, subject to the applicable source terms. Keep the paths in the input list below relative to a common input directory. The analysis uses cell-composition and purity covariates for the specified TCGA models and unadjusted comparisons for GEO external evaluation.

```sh
python analysis/run.py fixed --inputs ../misp-inputs --output ../misp-run
```

Replace `fixed` with the stage you need and supply your input and output paths.

| Stage | Analysis |
| --- | --- |
| `fixed` | Fixed-score reconstruction, paired correlations and intervals (S13, S21a), ICGC immune associations (S15), methylation–immune and incremental R² models (S21b, S23), and exploratory TCGA clinical and TMB tests. |
| `matched` | Matched random gene-set registry, comparisons and empirical tail probabilities. |
| `criteria` | TCGA composition-criterion sensitivity and alternative set weights with unadjusted GEO comparisons (S28–S29). |
| `selection-check` | Numerical spot checks for the TCGA selection models. |
| `selection` | Historical and revised three-window TCGA gene models, pathway tests, selection and robustness (S30). Use a fresh output directory after changing inputs. Full selection can take several hours. |
| `transfer` | Fixed and reselected program transfer and unadjusted GEO intervals (S33). Pass `--selection-dir ../misp-run/results_reselection` to use newly calculated weights. |
| `clinical-covariates` | TCGA scores with composition, purity and clinical covariates (S32a and S32c). |
| `tcga-composition` | TCGA DNA-based composition estimates. The source-helper option uses the EpiDISH 2.28.0 helper with its original GPL-2 notice. |
| `icgc-survival` | Recorded overall survival models for the ICGC cohorts (S6a). |

Normal-cell and single-cell reference summaries are provided in S26–S27. Published immune classification results, including the adjusted TCGA comparisons, are provided in S31 and S32b.

<details>
<summary>Prepared input files and SHA-256 checksums</summary>

Use trusted copies of the two Python pickle inputs. Hashes identify the exact inputs used to check the scripts.

| Relative path | Stage | Bytes | SHA-256 |
| --- | --- | ---: | --- |
| `MISP_Code/data/coverage_matched_random_sets.tsv` | matched | 190377 | `b1e3f5913e30f77d7c6efcc08531faecef87ee3fc8d1dd17afbf28ee7b2dc2b3` |
| `MISP_Code/data/EpiDISH_2.28.0_epidish.R` | tcga-composition (alternative to installed EpiDISH 2.28.0) | 8261 | `2017666e789a4fe47a0c0f01dadfa9e54c2ff17aa78be9568b9e86ae7782ad8d` |
| `MISP_Code/data/fixed_MISP_weights.tsv` | criteria, fixed, matched | 798 | `fb82c844191a08ed08a807936b797274c99c9f8ec6a22d3689a7894355a34da8` |
| `MISP_Code/data/gene_selection_audit.tsv` | criteria | 2104274 | `8e81fe4e3179930bc492fca227907bb9043d278670393be1932691483c11f43f` |
| `MISP_Code/data/GEO_scores.tsv.gz` | criteria, fixed, matched | 6149 | `e8792538fe40b3ee3e22ddab4ca3b747271f01f1c0cc07f93cd50ae08c3cca66` |
| `MISP_Code/data/GEO_standardized.npz` | criteria, fixed, matched | 18802708 | `4bcf6b4811019dda4d57549ff2ad62fc42e0dad8922d17874d3ab9a8fe3a92d5` |
| `MISP_Code/data/ICGC_scores.tsv.gz` | criteria, fixed, matched | 13647 | `f6ca9b33e6daa744c66707d2e8a3d9df85ad99f295b43467ecf7b296f5e88a90` |
| `MISP_Code/data/ICGC_standardized.npz` | criteria, fixed, matched | 10053603 | `49e506dbb29885fa2cfa30ddfc0fe32a66982c8456879a881fb5c62d1e60accd` |
| `MISP_Code/data/ICGC_survival_inputs.tsv` | icgc-survival | 13946 | `4883d271c3c4289eeb8cf76353025fce5f4c14daacb9ebb484d625e12a6edad1` |
| `MISP_Code/data/immune_and_Inflamed_expression.npz` | criteria, fixed, matched | 1338266 | `ca0f8525f0ca9175ffd6d864c1618f637931fc189366df68c3a9dcdb817a9d6a` |
| `MISP_Code/data/immune_definitions.json` | criteria, fixed, matched | 3078 | `13f7fbeaea2d64ee5d4fe85bfdf62e3029d5e4429b14c5a133b229fc37aecc29` |
| `MISP_Code/data/Inflamed_gene_list.tsv` | fixed | 131 | `8ec8dd6efaa6f5a784d1a72594cd13f82e50195d8df909f2d71fa21a9784bf1d` |
| `MISP_Code/data/random_background_matching_features.tsv` | matched | 926689 | `a085cad68e7cc19c98037542a3028e7b3825eb259d1a00dfbf1125cc41f5380a` |
| `MISP_Code/data/random_background_technical_coverage.tsv` | matched | 416320 | `6432cf0130b6d0f4e823abf23a3498d43af39a250c76ca3ebad48ab965835f9f` |
| `MISP_Code/data/TCGA_background_weight_coefficients.tsv` | matched | 332954 | `193b5110bbdd188575aeafa5cbd8ebe27f8c727212e9fdca4f1e28710436e410` |
| `MISP_Code/data/TCGA_centBloodSub.m_QC_noMISP.tsv` | tcga-composition | 557776 | `cfadd919e853e6160d737c9d621984b46e1cc5a6daff0ba6525ca053b67b4dc0` |
| `MISP_Code/data/TCGA_centBloodSub.m_reference.tsv` | tcga-composition | 13829 | `b00d471b0e3721f11c83eb9f39c15f2c05a5ab1eb27de35e9a35d59ea095c7bc` |
| `MISP_Code/data/TCGA_centEpiFibIC.m_QC_noMISP.tsv` | tcga-composition | 2024774 | `87bb9aa3755ad485c2a2fbc048c28ad6d0e7e03e0a52c184b1c6d8c647ba982a` |
| `MISP_Code/data/TCGA_centEpiFibIC.m_reference.tsv` | tcga-composition | 34748 | `8499d48b647268c0d00e0a8648b487865840e7d67044e99da2dbcf33526b7727` |
| `MISP_Code/data/TCGA_DNA.npz` | criteria, fixed, matched | 18464486 | `173653a9655b1d83651cf2fb6efc141db16520e23c946876d7ff2b174c369ac4` |
| `MISP_Code/data/TCGA_DNA_composition_proxies.tsv` | criteria, fixed, matched, tcga-composition | 104910 | `358a0d9aaaa4175a7d0145dd8b9d497f33ed9359dce74560052ccfec169bd693` |
| `MISP_Code/data/TCGA_RNA.npz` | criteria, fixed, matched | 15091406 | `46448056148626a5e65e4a7a75bab5d62e8b925f891dbcd83e607a5d90c6f151` |
| `MISP_Code/data/TCGA_scores_covariates.tsv.gz` | clinical-covariates, criteria, fixed, matched | 118443 | `755a75f44b4539c69797059d216b3ef1c7406c71cef1aefd87e2501091119b81` |
| `MISP_Code/data/TCGA_TMB_variant_filter_counts.tsv` | criteria, fixed, matched | 14334 | `f79a768f6e0cf4eaa64d168eeb36cc796768f32dfe7244903b92a21345e91c74` |
| `MISP_Extended/Inputs/frozen_gene_sets.pkl` | selection, selection-check | 9139948 | `c1286facf1c92b875df2b3c9b3a5c8293914463099154355399746c619344b52` |
| `MISP_Extended/Inputs/full_gene_inputs.pkl` | selection, selection-check | 202026243 | `217af9ae9266d3103e4449da6b8c725474eaa98dec5e3cbdce8a7e0141d9dab1` |
| `MISP_Extended/Inputs/Selected_Transfer/evaluation_covariates.npz` | transfer | 15470 | `31b550bff8ef8017d18c36b3a9ede8e32098b089195199d80b7bf2e33f6b7ce6` |
| `MISP_Extended/Inputs/Selected_Transfer/GEO_metadata_and_fixed30_QA.tsv` | transfer | 12562 | `2937d2c50fd341a9a1b38f1d1655f20c920ebc3376523a9f95cd6b09aada9dad` |
| `MISP_Extended/Inputs/Selected_Transfer/ICGC_metadata_and_fixed30_QA.tsv` | transfer | 28157 | `e292075dd69d642c159b7c502deaf065be462863b66685e7e8abf9858d28685b` |
| `MISP_Extended/Inputs/Selected_Transfer/immune_definitions.json` | transfer | 3078 | `13f7fbeaea2d64ee5d4fe85bfdf62e3029d5e4429b14c5a133b229fc37aecc29` |
| `MISP_Extended/Inputs/Selected_Transfer/immune_endpoint_expression.npz` | transfer | 1258148 | `b81cb48f1dc05363b0654fc4aa4ea70952029223c6940679713b4b8f7badebff` |
| `MISP_Extended/Inputs/Selected_Transfer/standardized_program_profiles.npz` | transfer | 1468011 | `fd3fb3a319da1aaf330828821600411773c9ce0e392ad4b1daecc58dc008113d` |
| `MISP_Extended/Inputs/Selected_Transfer/supplied_selection/historical_baseline/three_track_selected_genes.tsv` | transfer (with prepared weights) | 6677 | `7e811c26149dbab890709bea6e53e15ced2a2a9637619ae5b32b8f653d3d2985` |
| `MISP_Extended/Inputs/Selected_Transfer/supplied_selection/joint_DNA_ABSOLUTE/three_track_selected_genes.tsv` | transfer (with prepared weights) | 10840 | `e25b59cc3e74f7244ec442f939c9ddd06c0793684a269368b143c480229d965d` |
| `MISP_Extended/Inputs/Selected_Transfer/supplied_selection/joint_DNA_ABSOLUTE_age_sex/three_track_selected_genes.tsv` | transfer (with prepared weights) | 10449 | `0c78db2e2c1bf5ee2e2532419842f713fe4d6ddc2052f06e1d711dfc26e609ce` |
| `MISP_Extended/Inputs/Selected_Transfer/TCGA_metadata_and_fixed30_QA.tsv` | transfer | 53518 | `9dfde53afc0ad2828d7df4142682d79462be0f4972adb88d31436ee62c8cc9d2` |
| `MISP_Extended/Vendor/reference_project/expected_historical_summary.json` | selection | 679 | `240a99d8bc12c6fcb3019e1c916be59774fa5db8683b97fd9d9109dd3796c29c` |
| `MISP_Extended/Vendor/reference_project/fixed_primary_TCGA_scores.tsv` | selection, selection-check | 21094 | `a63110b86741d341295a7ac9a4a031db445fb2b9c2c32edc62c84036174ec0e9` |
| `MISP_Extended/Vendor/reference_project/GPL13667_exact_symbol_probe_counts.json` | selection | 317890 | `eef97b4a57249b326e4d002886628bc6c8263882521e256b8f8d0ff8826899a7` |
| `MISP_Extended/Vendor/reference_project/ICGC_exact_gene_symbols.json` | selection | 298189 | `8b751965468a2c32c5c45d488e71d8dc96efd0d2c72b2c3979197e71ac7d3d8a` |
| `MISP_Extended/Vendor/reference_project/phase1C_output/genomewide_methylation_candidate_tiers.tsv` | selection | 248620 | `6b95f5c55767626a306723311c5dce756f8e5eccdbb9753f073c7e93cb58d97a` |
| `MISP_Extended/Vendor/reference_project/phase1C_output/genomewide_ranked_pathway_enrichment.tsv` | selection | 21177 | `ab0dd34d7d3bdc7634bc8bdc47ed1db51fa33bcef3d03d83d35ff5156ee337a9` |
| `MISP_Extended/Vendor/reference_project/phase1C_output/intermediate/association_promoter_union.tsv.gz` | selection | 13084 | `fa8e084a43a6b15fdfaae7c0f9d4f2caa217bafede0f83250ec0ff8e623cb212` |
| `MISP_Extended/Vendor/reference_project/phase1C_output/intermediate/association_TSS1500.tsv.gz` | selection | 12384 | `4b297bdc0a797776f10ae32597ffe7ffc48cf5d2c61f531e8ee0ca1105788403` |
| `MISP_Extended/Vendor/reference_project/phase1C_output/intermediate/association_TSS200.tsv.gz` | selection | 13045 | `b17824986256cf936149d8cc1260e69e73ab8f1225ff753daf9e5a2249064f58` |
| `MISP_Extended/Vendor/reference_project/phase2A_output/genomewide_purity_adjusted_models.tsv.gz` | selection | 51180 | `60a8af7d6023bce43eafedda8ccfd46bf60249f949066e7111b4564f0ffdff9b` |
| `MISP_Extended/Vendor/reference_project/phase2A_output/locked_immune_methylation_program_genes.tsv` | selection, selection-check | 798 | `fb82c844191a08ed08a807936b797274c99c9f8ec6a22d3689a7894355a34da8` |
| `MISP_Extended/Vendor/reference_project/phase2A_output/MISP_score_robustness.tsv` | selection | 17925 | `cc0f9a25e6ba1d723de5a1bb25ef032b626cf4e5d6529ad91d2cbe1885057f7c` |
| `MISP_Extended/Vendor/reference_project/smoke_model_targets.tsv` | selection, selection-check | 47541 | `9374c32538af87aedd096872b115e4ba1cce3ce140dbc5885e21df361b14b1b9` |

</details>

<details>
<summary>Original numerical modules</summary>

- [phase1C_genomewide.py](analysis/MISP_Extended/Vendor/phase1C_genomewide.py) was copied from `scripts/python/phase1C/phase1C_genomewide.py`. SHA-256 `ec15d3bb8cc29c8a3c4d753d0e0f6db4f9e8d5ff76a129db9b72b49a3142dceb`.
- [phase1C_integrate.py](analysis/MISP_Extended/Vendor/phase1C_integrate.py) was copied from `scripts/python/phase1C/phase1C_integrate.py`. SHA-256 `8d065881220dd081ef7078f8076e960130a216a5925634227e4355e3b6006ddf`.
- [phase2A_tcga_models.py](analysis/MISP_Extended/Vendor/phase2A_tcga_models.py) was copied from `scripts/phase2A/phase2A_tcga_models.py`. SHA-256 `534e5a1981d35700dbcf6cf4d4e352c569e24409af16f56058187e8596a7b8eb`.
- [phase2A_lock_program.py](analysis/MISP_Extended/Vendor/phase2A_lock_program.py) was copied from `scripts/phase2A/phase2A_lock_program.py`. SHA-256 `c0df489d3839b728a1cb470f956e5bc6911f14435f0abd4fb0456b204e1c6bad`.

</details>
