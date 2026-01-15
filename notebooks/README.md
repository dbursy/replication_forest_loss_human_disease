# Python Data Processing Notebooks

This directory contains Python notebooks used to process raw data from the sources cited in the original study and to construct the analysis datasets used in the replication and extension.

The *extended* datasets incorporate the additional **2018 DHS survey wave**, whereas the original study uses the **2008 and 2013 DHS waves only**.

## Overview of notebooks

- **`1_Preparation_Survey.ipynb`**  
  Imports and harmonizes DHS survey microdata used to measure child disease incidence and household and demographic controls.

- **`2_Preparation_Forest.ipynb`**  
  Processes forest cover and forest loss data and constructs measures consistent with the original study.

- **`3_Prepration_Lights.ipynb`**  
  Prepares night-time lights data and constructs covariates used as proxies for economic activity.

- **`4_Preparation_Soil.ipynb`**  
  Processes soil and environmental covariates used in the analysis.

- **`5_Final_Dataframe.ipynb`**  
  Merges the prepared survey and spatial covariates to create the final replication dataset based on the 2008 and 2013 DHS waves.

- **`6_Preparation_Forest_Extended.ipynb`**  
  Extends the forest data processing pipeline for use with the additional 2018 DHS survey wave.

- **`7_Final_Dataframe_Extended.ipynb`**  
  Constructs the extended analysis dataset by incorporating the 2018 DHS wave alongside the 2008 and 2013 waves.

- **`8_Descriptive_Statistics.ipynb`**  
  Produces descriptive statistics and summary checks for both the replication and extended datasets.

## Output

The outputs of this pipeline are:
- the replication analysis dataset (2008 and 2013 DHS waves), and  
- the extended analysis dataset (2008, 2013, and 2018 DHS waves),

which are subsequently used by the Stata do-files to reproduce tables and figures.
