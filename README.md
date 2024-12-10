# MCT_Fly_Val_public
 
Code and raw data to support Terry (2025) *Experimental validation of ecological coexistence theory to forecast extirpation under rising temperatures.*

# Contents

Top level analysis documents are held in `Analysis` folder:

1. Data preparation of raw data to more useful forms
2. Key summaries of data used in paper and analysis of extinction points
3. Fitting and comparison of suite of models
4. Inferences from best fit models
5. Inferences from alternative models fit using a Poisson Error
6. Simulation-based estimate of long-term impact of variability based on models and decomposition of impact of different processes
7. Tests using direct simulation
8. Predictions made not using certain experiments for both training and testing the models. 

Each are included in both .rmd and knited html form. 

`Scripts/` contains .R code specifying functions used to generate brms code to fit the suite of models. It also contains copies of the STAN code used for the best fit models. 

`AutoModels/` contains three folders:
-`Scripts/` contains the brms code generated specifying the models
-`Fits/` stores the large model fits generated by the scripts. For space reasons this is empty
-`KeyFits/` contains the best_fit models for each species as an accessible archive

`InitialDataStorage/` contains the raw data tables before any processing. 

`TubeDoc_CC3_ALL.csv` details each tube's metadata. The columns are:

- `CompTreat`  Competitor treatment 'COMP' or 'MONO'
- `TempTreat`  Temperature treatment 'STEADY' or VAR
- `VarBatch`   Which Trajectory, either 'S' for 'V01' to 'V12' for the twelve trajectories
- `Gen`        Generation the count came from. G01 to G10. 
- `ReplicateST`Replicate ID within steady temperatures (1:60)
- `ReplicateVA`Replicate ID within variable temperature trajectory (1:5)
- `Replicate`  Either of the preceding two columns as approproate. 
- `TempDiff`   Temperature Difference compared to mean
- `Incubator`  Which incubator ('Low', 'High' or 'Steady') the tube should be in 
- `MeanTemp`   The mean temperature of that generation
- `IncubTemp`  The Temperature of that incubator in that generation
- `LINE_ID`    Line ID, pasting together treatments and replicate values. 
- `TUBE_ID`    Unique tube ID, 1001:3400

`TubeDoc_CC3_ForDataEntry.csv` has columns detailing counts of each sex of each species `PAL_FEMALE`,`PAL_MALE`,`PAN_FEMALE`,`PAN_MALE`. Blank entries are zeros. 

 `TubeDoc_CC3_ALL.csv` relates to `TubeDoc_CC3_ForDataEntry.csv` via `Gen`,`IncubTemp`,`LINE_ID` & `TUBE_ID`. 

`Tables` and `Figures` both hold outputs from the analysis documents

`ProcessedData/` Includes dataframes with data wrangled for the model fitting. NB `_comp` columns are counts to be used in determining intraspecific competition, i.e. count-1, floored at zero, to avoid self-competition.  

# Package Versions

Code was originally built using R version 4.3.0, brms version 2.19 and tidyverse version 2.0.0. Revisions were made using R 4.4. All analysis documents include `sessionInfo()` in the htmls, which may help resolve any future package mystery.

# Licensing:

Code available for free reuse without restriction or warranty.
Data may be freely reused with attribution (ideally citing the paper: Terry (2025) "An experimental validation test of ecological coexistence theory to forecast extinction under rising temperatures" *Ecology Letters* . Preprint: https://doi.org/10.1101/2024.02.22.581553). 
