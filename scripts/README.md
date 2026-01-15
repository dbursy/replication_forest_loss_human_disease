# Stata Replication Scripts

This repository contains Stata do-files that replicate and extend the results of:

Berazneva, Julia, and Tanya S. Byker. 2017. "Does Forest Loss Increase Human Disease? Evidence from Nigeria." American Economic Review 107 (5): 516–21. [Available here](https://doi.org/10.1257/aer.p20171132)

## `01_replication_original_with_annotations.do`

Replicates the original results using the authors’ AEA Replication Package (113539-V1).

- No substantive changes to the original analysis.
- Reproduces **Table 1**.
- Adds code to reproduce **Figure 2**.
- Code for **Figure 1** is provided in the Python notebook `5_Final_Dataframe.ipynb`.

Supporting materials, including a screenshot documenting female population of reproductive age used for updated weight de-normalization, are available in the `resources` directory.

## `02_replication_original_plus_extended_data.do`

Replicates and extends the main results using both the original replication dataset and a replicated, extended dataset constructed from raw data.

- Based on the authors’ original replication script.
- All extensions and deviations from the original analysis are documented in the code comments.