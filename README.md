# Boolean Network Inference using Minimum Description Length
Project ECEN 647

Usage:
1. `CreateRandNetwork.R`: Create a random boolean network using BoolNet package in R and simulate time series from that network. OR use provided .sbml file to generate time series from given Boolean Network using the same package.
2. `unix_processing.sh`: Process the generated csv files to matlab scripts using provided unix commands. Execute script from inside the folder containing csv files.
3. `MDL_inference.m`: Run the main matlab script after changing the directory loacation of data .m scripts generated in step 2. `Cml_calc.m` function file should be in the same folder as `MDL_inference.m`.

# Boolean Network used for testing:

 Involved genes:\
   Gene1 Gene2 Gene3 Gene4 Gene5 Gene6 Gene7 Gene8 Gene9 Gene10
 
 Transition functions:\
 Gene1 = ((!Gene2 & !Gene3) | Gene10)\
 Gene2 = (!Gene9 & Gene6 & !Gene8)\
 Gene3 = ((Gene9 & Gene10) | (Gene1 & Gene10))\
 Gene4 = (!Gene2 | (!Gene5 & !Gene6))\
 Gene5 = ((!Gene1 & Gene6) | (Gene1 & !Gene4 & !Gene6))\
 Gene6 = ((!Gene6 & !Gene10) | (Gene6 & !Gene1 & Gene10))\
 Gene7 = (!Gene8 | (Gene6 & !Gene7))\
 Gene8 = ((!Gene4 & !Gene7) | (Gene4 & Gene1 & Gene7))\
 Gene9 = ((!Gene9 & !Gene6) | (!Gene2 & !Gene6) | (Gene9 & Gene6))\
 Gene10 = ((Gene4 & Gene5) | (Gene7 & Gene5) | (Gene7 & Gene4))

# Results:

### Best predictor sets:
|Gene of Interest|Actual Predictors|Estimated Predictors Ts=10,Ns=10|Estimated Predictors Ts=10,Ns=100|Estimated Predictors Ts=100,Ns=100|Estimated Predictors Ts=100,Ns=1000|
| -------------- |:---------------:| :-----------------------------:| :------------------------------:|:------------------------------:|:------------------------------:|
|Gene1           |2,3,10           |3,10|2,3,10|2,3,10|2,3,10|
|Gene2           |6,8,9            |2|6,9|6,9|6,9|
|Gene3           |1,9,10           |10|9,10|9,10|9,10|
|Gene4           |2,5,6            |2|2,6|2,6|2,5|
|Gene5           |1,4,6            |1,4,6|1,4,6|1,4,6|1,4,6|
|Gene6           |1,6,10           |6,10|1,6,10|1,6,10|1,6,10|
|Gene7           |6,7,8            |8|2,8|6,8|6,8|
|Gene8           |1,4,7            |1,4,7|1,4,7|1,4,7|1,4,7|
|Gene9           |2,6,9            |2|2,6,9|2,6,9|2,6,9|
|Gene10          |4,5,7            |4,5,7|4,5,7|4,5,7|4,5,7|

### Rank of actual predictor set in results (from 175 total sets for each gene):
|Gene of Interest|Actual Predictors|Rank Ts=10,Ns=10|Rank Ts=10,Ns=100|Rank Ts=100,Ns=100|Rank Ts=100,Ns=1000|
| -------------- |:---------------:| :-----------------------------:| :------------------------------:|:------------------------------:|:------------------------------:|
|Gene1           |2,3,10           |2|1|1|1|
|Gene2           |6,8,9            |92|17|19|5|
|Gene3           |1,9,10           |9|10|10|10|
|Gene4           |2,5,6            |85|27|30|24|
|Gene5           |1,4,6            |1|1|1|1|
|Gene6           |1,6,10           |2|1|1|1|
|Gene7           |6,7,8            |43|34|18|22|
|Gene8           |1,4,7            |1|1|1|1|
|Gene9           |2,6,9            |25|1|1|1|
|Gene10          |4,5,7            |1|1|1|1|

# After accounting for dynamic errors:

# References:
1. [Inference of Gene Regulatory Networks Based on a Universal Minimum Description Length, John Dougherty, Ioan Tabus & Jaakko Astola](https://www.ncbi.nlm.nih.gov/pubmed/18437238)
2. [Normalized Maximum Likelihood Models for Boolean Regression with Application to Prediction and Classification in Genomics, Ioan Tabus, Jorma Rissanen, Jaakko Astola](https://link.springer.com/chapter/10.1007/0-306-47825-0_10)
