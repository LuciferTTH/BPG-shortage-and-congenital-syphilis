R codes for  project of association between global benzathine benzylpenicillin shortages and resurgence of congenital syphilis: a multi-country cross-sectional study

## Script contents

| File | Main work |
|---|---|
| 01_descriptive_analysis.R | Input and panel checks, annual rates by shortage-report status, and included versus unclassified country summaries. |
| 02_primary_and_subgroup_analysis.R | Primary weighted model, counterfactual contrasts, main and subgroup figures, country-selection figure, and DAG adjustment-set check. |
| 03_sensitivity_analyses.R | Alternative specifications and exclusions, WHO ANC screening and service adjustment, QBA, alternative break years, permutation analysis, and sensitivity figures. |
| 04_additional_analyses.R | WHO-GBD comparison, IQVIA data reconstruction and mixed models, HIV service-data audit and adjusted models, and publication table and figure inventory. |

## Required and optional inputs

All four files currently require both prepared panels in data/processed_data/: data_combined.Rdata, containing an object named data_combined, and data_combined_hiv.Rdata, containing an object named data_combined_hiv. The panels are not included in this public archive. The data README lists their expected contents and the optional raw and curated inputs.

Optional sections run only when their required files are present. The WHO service sensitivity, WHO-GBD comparison, IQVIA procurement analysis, HIV source join, and HIV negative-control models each have their own input requirements. A section may be skipped when its complete input set is absent. A partial but malformed input set can stop with an error. Review the console message and logs before interpreting an output as newly generated.

## Outputs

Scripts write generated figures to figures/, tables to tables/, intermediate model results to results/, and input or runtime audits to logs/. The CSV files already in tables/ and the PNG assets already in figures/ are supplied reference outputs. Running a script can overwrite files with the same names.

Do not interpret a copied figure asset or an existing result table as evidence that the corresponding analysis was refitted during a new run. See the main README and FIGURE_PROVENANCE.md.

