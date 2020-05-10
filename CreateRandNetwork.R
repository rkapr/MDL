## R code for generating random boolean network with 10 nodes and 3 regulators per node:
library(BoolNet)
Rand_network1 = generateRandomNKNetwork(10,3)
# Boolean network with 10 genes
# 
# Involved genes:
#   Gene1 Gene2 Gene3 Gene4 Gene5 Gene6 Gene7 Gene8 Gene9 Gene10
# 
# Transition functions:
#   Gene1 = <f(Gene2,Gene3,Gene10){11010101}>
#   Gene2 = <f(Gene9,Gene6,Gene8){00100000}>
#   Gene3 = <f(Gene1,Gene9,Gene10){00010101}>
#   Gene4 = <f(Gene5,Gene2,Gene6){11101100}>
#   Gene5 = <f(Gene1,Gene4,Gene6){01011000}>
#   Gene6 = <f(Gene6,Gene1,Gene10){10100100}>
#   Gene7 = <f(Gene6,Gene7,Gene8){10101110}>
#   Gene8 = <f(Gene4,Gene1,Gene7){10100001}>
#   Gene9 = <f(Gene2,Gene9,Gene6){10111001}>
#   Gene10 = <f(Gene7,Gene4,Gene5){00010111}>

# Save file to load network in future

toSBML(Rand_network1,file="Rand_network1/Rand_network1.sbml", generateDNFs = FALSE,saveFixed = TRUE)
print(loadSBML("Rand_network1/Rand_network1.sbml"))
# Boolean network with 10 genes
# 
# Involved genes:
#   Gene1 Gene2 Gene3 Gene4 Gene5 Gene6 Gene7 Gene8 Gene9 Gene10
# 
# Transition functions:
#   Gene1 = ((!Gene2 & !Gene3) | Gene10)
# Gene2 = (!Gene9 & Gene6 & !Gene8)
# Gene3 = ((Gene9 & Gene10) | (Gene1 & Gene10))
# Gene4 = (!Gene2 | (!Gene5 & !Gene6))
# Gene5 = ((!Gene1 & Gene6) | (Gene1 & !Gene4 & !Gene6))
# Gene6 = ((!Gene6 & !Gene10) | (Gene6 & !Gene1 & Gene10))
# Gene7 = (!Gene8 | (Gene6 & !Gene7))
# Gene8 = ((!Gene4 & !Gene7) | (Gene4 & Gene1 & Gene7))
# Gene9 = ((!Gene9 & !Gene6) | (!Gene2 & !Gene6) | (Gene9 & Gene6))
# Gene10 = ((Gene4 & Gene5) | (Gene7 & Gene5) | (Gene7 & Gene4))

# Simulate 1000 datasets of 100 sample time series from the network
samples_RN1 = generateTimeSeries(Rand_network1,numSeries = 1000,numMeasurements=100);
for (i in 1:1000) 
{write.table(samples_RN1[[i]],paste0("Rand_network4/random_network4_",i,".csv"), quote = FALSE,col.names = FALSE,sep = ",")}
