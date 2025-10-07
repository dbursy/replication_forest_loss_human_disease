# Replication Exercise:  Forest Loss and Human Disease

The aim of this repository is to conduct a replication study of the following research paper:

Berazneva, Julia, and Tanya S. Byker. 2017. "Does Forest Loss Increase Human Disease? Evidence from Nigeria." American Economic Review 107 (5): 516–21. [Available here](https://doi.org/10.1257/aer.p20171132)

## 🚀 Prerequisites

### Virtual Environment

Create and activate the python virtual environment in the local project folder
 
```bash
python3 -m venv env
source env/bin/activate
```

If this step is successfull the comand line displays the prefix `(env)`. Continue by installing required Python packages
 
```bash
python -m pip install -r requirements.txt
```
 
Note: To create or update the requirements file, please run
 
```bash
pip freeze > requirements.txt
```

### Data Directory

The data directory is structured as follows

```bash
datassets
├── dhs_survey
│   ├── 2008
│   │   ├── NGGE52FL
│   │   └── NGKR53DT
│   ├── 2013
│   │   ├── NGGE6AFL
│   │   └── NGKR6ADT
│   └── 2018
│       ├── NGGE7BFL
│       └── NGKR7BDT
├── forest_change
│   ├── lossyear
│   │   ├── Hansen_GFC2015_lossyear_10N_000E.tif
│   │   ├── Hansen_GFC2015_lossyear_10N_010E.tif
│   │   ├── Hansen_GFC2015_lossyear_20N_000E.tif
│   │   └── Hansen_GFC2015_lossyear_20N_010E.tif
│   ├── lossyear_2019
│   │   ├── Hansen_GFC-2019-v1.7_lossyear_10N_000E.tif
│   │   ├── Hansen_GFC-2019-v1.7_lossyear_10N_010E.tif
│   │   ├── Hansen_GFC-2019-v1.7_lossyear_20N_000E.tif
│   │   └── Hansen_GFC-2019-v1.7_lossyear_20N_010E.tif
│   └── treecover
│       ├── Hansen_GFC-2019-v1.7_treecover2000_10N_000E.tif
│       ├── Hansen_GFC-2019-v1.7_treecover2000_10N_010E.tif
│       ├── Hansen_GFC-2019-v1.7_treecover2000_20N_000E.tif
│       └── Hansen_GFC-2019-v1.7_treecover2000_20N_010E.tif
── luminosity
│   ├── F162005.v4
│   ├── F162006.v4
│   ├── F162007.v4
│   ├── F162008.v4
│   ├── F182010.v4
│   ├── F182011.v4
│   ├── F182012.v4
│   └── F182013.v4
├── map_africa
│   ├── nga_admbnda_adm2_osgof_20190417.shp
│   └── ...
└── soil
    ├── af_CEC_T__M_sd1_250m.tif
    ├── af_CEC_T__M_sd2_250m.tif
    ├── af_ORCDRC_T__M_sd1_250m.tif
    ├── af_ORCDRC_T__M_sd2_250m.tif
    ├── af_PHIHOX_T__M_sd1_250m.tif
    └── af_PHIHOX_T__M_sd2_250m.tif
```
