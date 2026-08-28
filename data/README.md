# Data contract

Patient-level controlled-access data, access tokens, and large raw matrices must not be committed.

## Required processed files

Each cohort requires a gene-expression CSV and a clinical CSV. The first column of each file must be `sample_id`. Expression files contain one row per patient and one column per gene. Clinical files must contain `sample_id`, `OS_time`, and `OS_status`, where `OS_time` is a positive follow-up duration in days and `OS_status` is coded 0 for censored/alive and 1 for death. Sample identifiers must be unique.

Validation files use the same schema. Gene-expression units and preprocessing must match the training cohort. The scripts stop if IDs are duplicated, survival coding is invalid, required genes are absent, or the matched sample count is too small.

## Public sources

The supplied code names TCGA and TARGET only generically. Before publication, the authors must add the exact portal, project/cohort identifiers, data type, workflow/version, download URL, release or download date, access date, and preprocessing steps for every cohort. Do not guess these fields from cancer abbreviations alone.

## Provenance fields to complete

| Cohort | Database/project | Data product | Version/date | URL | Access date | Preprocessing |
|---|---|---|---|---|---|---|
| LGG | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED |
| LUAD | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED |
| LIHC | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED |
| PAAD | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED |
| LAML | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED | AUTHOR_INPUT_NEEDED |

