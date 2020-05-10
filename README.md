# MDL
Project ECEN 647

Usage:
1. Create a random boolean network using BoolNet package in R and simulate time series from that network. OR use provided .sbml file to generate time series from given Boolean Network using the same package.\
2. Process the generated csv files to matlab scripts using provided unix commands. Execute script from inside the folder containing csv files.\
3. Run the main matlab script after changing the directory loacation of data .m scripts generated in step 2.

# Boolean Network usind for testing:

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
 Gene10 = ((Gene4 & Gene5) | (Gene7 & Gene5) | (Gene7 & Gene4))\

# Results:

|Gene of Interest|Actual Predictors|Estimated Predictors Ts=10,Ns=10|Estimated Predictors Ts=10,Ns=100|
| -------------- |:---------------:| :-----------------------------:| :------------------------------:|
|Gene1           |2,3,10           |2,7                             |2,3,10                           |
Gene2
6,8,9
NULL
NULL
Gene3
1,9,10
10
9,10
Gene4
2,5,6
2
2
Gene5
1,4,6
1,2,4
1,4,6
Gene6
1,6,10
6,10
1,6,10
Gene7
6,7,8
8
2,8
Gene8
1,4,7
1,4,7
1,4,7
Gene9
2,6,9
2
2,6,9
Gene10
4,5,7
4,5,7
4,5,7
