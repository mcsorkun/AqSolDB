# Aqueous Solubility Data Curation

AqSolDB is created by the Autonomous Energy Materials Discovery [AMD] research group, consists of aqueous solubility values of 9,982 unique compounds curated from 9 different publicly available aqueous solubility datasets. This openly accessible dataset, which is the largest of its kind, and will not only serve as a useful reference source of measured solubility data, but also as a much improved and generalizable training data source for building data-driven models.


## Citation

If you use AqSolDB in your study, please cite the following paper.

Paper: [Nature Scientific Data](https://doi.org/10.1038/s41597-019-0151-1) - https://doi.org/10.1038/s41597-019-0151-1


## Other Platforms 

> Online reproducible code: [Code Ocean](https://doi.org/10.24433/CO.1992938.v1)

> Data science kernels: [Kaggle](https://www.kaggle.com/sorkun/aqsoldb-a-curated-aqueous-solubility-dataset)

> Easy search from AqSolDB: [AMD Website](https://www.amdlab.nl/database/AqSolDB/)

> Data Repository: [Harvard Dataverse](https://doi.org/10.7910/DVN/OVHAW8) 


## Overview

This repository has been developed in order to curate various aqueous solubility datasets into a broad and extensive dataset called AqSolDB.

The curation process in this work can be accomplished by executing two python scripts in the given sequence:
1. data-preprocess.py  - for pre-processing the raw data set to a standardized format
2. data-curation.py - for merging the standardized datasets, assigning reliability lables and adding 2D descriptors

These two python scripts call upon functions from other python modules that are defined in:
- preprocess.py
- merge.py
- descriptors.py

Further information about curation process can be found in the associated manuscript.

# Examples 

## data-preprocess.py 

This file converts 2 example sub-datasets (25 instances from raw forms of dataset-A[1] and dataset-H[6]) which are then converted into a standardized format. (This is an example how to preprocess datasets. The preporcessed data files already in the data folder.)

**inputs:**
1. raw-dataset-A.csv (various solubility metrics (g/L, mg/L..) with Name and CAS Number)
2. raw-dataset-H.csv (has solubility values(LogS) with SLN representations) 

**outputs:**
1. dataset-A.csv 
2. dataset-H.csv 

### Note

To apply this method to your own dataset, perform the following steps:
1. Check the available properties, representations, and solubility units of your dataset
2. Select the suitable preprocessing methods from the "preprocess.py" module.

## data-curation.py 

This file curates, i.e., merges datasets, selects most reliable values among multiple occurences, and adds 2D descriptors from 9 different standardized datasets that are obtained after the pre-processing step. 

**inputs:**
1. dataset-A.csv [1]
2. dataset-B.csv [2]
3. dataset-C.csv [3]
4. dataset-D.csv [4]
5. dataset-E.csv [5]
6. dataset-F.csv [6]
7. dataset-G.csv [7]
8. dataset-H.csv [6]
9. dataset-I.csv [8]

**outputs:**
1. dataset_curated.csv

### Note

To apply this method, your input dataset should be in the standardized format (output of preprocessing) having following columns:
- ID
- Name
- InChI
- InChIKey
- SMILES
- Solubility
- Prediction

# License and Copyright
MIT License

Copyright (c) 2019 Murat Cihan Sorkun

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



# References

[1] eChemPortal - The Global Portal to Information on Chemical Substances. https://www.echemportal.org/.

[2] Meylan, W. M.  Preliminary Report:  Water Solubility Estimation by Base Compound Modification.Environmental Science Center, Syracuse, NY (1995).

[3] Raevsky, O. A., Grigor’ev, V. Y., Polianczyk, D. E., Raevskaja, O. E.& Dearden, J. C. Calculation of aqueous solubility of crystalline un-ionized organic chemicals and drugs based on structural similarity and physicochemical descriptors.Journal of Chemical Information and Computer Sciences 54, 683–691 (2014).

[4] Meylan, W. M., Howard, P. H.  Upgrade of PCGEMS Water Solubility Estimation Method. Environmental Science Center, Syracuse, NY(1994)

[5] Huuskonen, J. Estimation of aqueous solubility for a diverse set of organic compounds based on molecular topology.Journal of Chemical Informationand Computer Sciences 40, 773–777 (2000).

[6] Wang, J., Hou, T. & Xu, X. Aqueous solubility prediction based on weighted atom type counts and solvent accessible surface areas. Journal of Chemical Information and Modeling 49, 571–581 (2009).

[7] Delaney, J. S. ESOL: estimating aqueous solubility directly from molecular structure. Journal of Chemical Information and Computer Sciences 44,1000–1005 (2004).

[8] Llinas, A., Glen, R. C. & Goodman, J. M. Solubility challenge: can you predict solubilities of 32 molecules using a database of 100 reliable measurements?.Journal of Chemical Information and Modeling 48, 1289–1303 (2008).
