# Combined Analysis of Individual SCN5A SP Experiments

# STEP 1: Combine each parameter from each experiment
# STEP 2: Analyze these combined datasets, perform outlier removal when appropriate
# STEP 3: Merge each parameter for final daataset used in manuscript and variant interpretation

library(stringr)
library(tidyverse)

########## STEP 1: Combine each parameter from each experiment ##########
# Note: Can edit working directory by adding workingDir='path/to/workingDir' when calling processMissing()

# STEP 1A: Combine SPs together (peak current and voltage of activation)
# STEP 1B: Combine SPs together (voltage of inactivation)
# STEP 1C: Combine SPs together (late current)
# STEP 1D: Combine SPs together (Recovery From Inactivation)
# STEP 1E: Combine SPs together (Inactivation time)
# STEP 1F: Combine SPs together (Peak current -90mV)

##### STEP 1A: Combine SPs together (peak current and voltage of activation) #####
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
# Excluding early SP runs from Calibration Manuscript 
#result=processMissing(workingDir,'SP1/Analysis/SP1_NaIV_processed3.csv','SP1/Analysis/SP1_missing.csv',result,'SP1',48)
#result=processMissing('SP2/Analysis/SP2_NaIV_processed3.csv','SP2/Analysis/SP2_missing.csv',result,'SP2',48)
#result=processMissing('SP3/Analysis/SP3_NaIV_processed3.csv','SP3/Analysis/SP3_missing.csv',result,'SP3',48)
#result=processMissing('SP4/Analysis/SP4_NaIV_processed3.csv','SP4/Analysis/SP4_missing.csv',result,'SP4',48)
#result=processMissing('SP5/Analysis/SP5_NaIV_processed3.csv','SP5/Analysis/SP5_missing.csv',result,'SP5',48)
result=processMissing('SP6+7/Analysis/SP6+7_NaIV_processed3.csv','SP6+7/Analysis/SP6+7_missing.csv',result,'SP6+7',48)
result=processMissing('SP8+9/Analysis/SP8+9_NaIV_processed3.csv','SP8+9/Analysis/SP8+9_missing.csv',result,'SP8+9',48)
result=processMissing('SP10/Analysis/SP10_NaIV_processed3.csv','SP10/Analysis/SP10_missing.csv',result,'SP10',48)
result=processMissing('SP11/Analysis/SP11_NaIV_processed3.csv','SP11/Analysis/SP11_missing.csv',result,'SP11',48)
result=processMissing('SP12/Analysis/SP12_NaIV_processed3.csv','SP12/Analysis/SP12_missing.csv',result,'SP12',48)
result=processMissing('SP13/Analysis/SP13_NaIV_processed3.csv','SP13/Analysis/SP13_missing.csv',result,'SP13',48)
result=processMissing('SP14/Analysis/SP14_NaIV_processed3.csv','SP14/Analysis/SP14_missing.csv',result,'SP14',48)
result=processMissing('SP15/Analysis/SP15_NaIV_processed3.csv','SP15/Analysis/SP15_missing.csv',result,'SP15',48)
result=processMissing('SP17/Analysis/SP17_NaIV_processed3.csv','SP17/Analysis/SP17_missing.csv',result,'SP17',48)
result=processMissing('SP18/Analysis/SP18_NaIV_processed3.csv','SP18/Analysis/SP18_missing.csv',result,'SP18',48)
result=processMissing('SP19/Analysis/SP19_NaIV_processed3.csv','SP19/Analysis/SP19_missing.csv',result,'SP19',48)
result=processMissing('SP20/Analysis/SP20_NaIV_processed3.csv','SP20/Analysis/SP20_missing.csv',result,'SP20',48)
result=processMissing('SP21/Analysis/SP21_NaIV_processed3.csv','SP21/Analysis/SP21_missing.csv',result,'SP21',48)
result=processMissing('SP22/Analysis/SP22_NaIV_processed3.csv','SP22/Analysis/SP22_missing.csv',result,'SP22',48)
result=processMissing('SP24/Analysis/SP24_NaIV_processed3.csv','SP24/Analysis/SP24_missing.csv',result,'SP24',48)
result=processMissing('SP26/Analysis/SP26_NaIV_processed3.csv','SP26/Analysis/SP26_missing.csv',result,'SP26',48)
result=processMissing('SP27/Analysis/SP27_NaIV_processed3.csv','SP27/Analysis/SP27_missing.csv',result,'SP27',48)
result=processMissing('SP28/Analysis/SP28_NaIV_processed3.csv','SP28/Analysis/SP28_missing.csv',result,'SP28',48)
result=processMissing('SP30/Analysis/SP30_NaIV_processed3.csv','SP30/Analysis/SP30_missing.csv',result,'SP30',48)
result=processMissing('SP31/Analysis/SP31_NaIV_processed3.csv','SP31/Analysis/SP31_missing.csv',result,'SP31',48)
result=processMissing('SP32/Analysis/SP32_NaIV_processed3.csv','SP32/Analysis/SP32_missing.csv',result,'SP32',48)
result=processMissing('SP35/Analysis/SP35_NaIV_processed3.csv','SP35/Analysis/SP35_missing.csv',result,'SP35',48)
result=processMissing('SP36/Analysis/SP36_NaIV_processed3.csv','SP36/Analysis/SP36_missing.csv',result,'SP36',48)
result=processMissing('SP37/Analysis/SP37_NaIV_processed3.csv','SP37/Analysis/SP37_missing.csv',result,'SP37',48)
result=processMissing('SP38/Analysis/SP38_NaIV_processed3.csv','SP38/Analysis/SP38_missing.csv',result,'SP38',48)
result=processMissing('SP39/Analysis/SP39_NaIV_processed3.csv','SP39/Analysis/SP39_missing.csv',result,'SP39',48)
result=processMissing('SP40/Analysis/SP40_NaIV_processed3.csv','SP40/Analysis/SP40_missing.csv',result,'SP40',48)
result=processMissing('SP41/Analysis/SP41_NaIV_processed3.csv','SP41/Analysis/SP41_missing.csv',result,'SP41',48)
result=processMissing('SP43/Analysis/SP43_NaIV_processed3.csv','SP43/Analysis/SP43_missing.csv',result,'SP43',48)
result=processMissing('SP44/Analysis/SP44_NaIV_processed3.csv','SP44/Analysis/SP44_missing.csv',result,'SP44',48)
result=processMissing('SP46/Analysis/SP46_NaIV_processed3.csv','SP46/Analysis/SP46_missing.csv',result,'SP46',48)
result=processMissing('SP47/Analysis/SP47_NaIV_processed3.csv','SP47/Analysis/SP47_missing.csv',result,'SP47',48)
result=processMissing('SP48/Analysis/SP48_NaIV_processed3.csv','SP48/Analysis/SP48_missing.csv',result,'SP48',48)
result=processMissing('SP49/Analysis/SP49_NaIV_processed3.csv','SP49/Analysis/SP49_missing.csv',result,'SP49',48)
result=processMissing('SP50/Analysis/SP50_NaIV_processed3.csv','SP50/Analysis/SP50_missing.csv',result,'SP50',48)
result=processMissing('SP51/Analysis/SP51_NaIV_processed3.csv','SP51/Analysis/SP51_missing.csv',result,'SP51',48)
result=processMissing('SP61/Analysis/SP61_NaIV_processed3.csv','SP61/Analysis/SP61_missing.csv',result,'SP61',48)
result=processMissing('SP63/Analysis/SP63_NaIV_processed3.csv','SP63/Analysis/SP63_missing.csv',result,'SP63',48)
result=processMissing('SP65/Analysis/SP65_NaIV_processed3.csv','SP65/Analysis/SP65_missing.csv',result,'SP65',48)
result=processMissing('SP67/Analysis/SP67_NaIV_processed3.csv','SP67/Analysis/SP67_missing.csv',result,'SP67',48)
result=processMissing('SP68/Analysis/SP68_NaIV_processed3.csv','SP68/Analysis/SP68_missing.csv',result,'SP68',48)
result=processMissing('SP69/Analysis/SP69_NaIV_processed3.csv','SP69/Analysis/SP69_missing.csv',result,'SP69',48)
result=processMissing('SP71/Analysis/SP71_NaIV_processed3.csv','SP71/Analysis/SP71_missing.csv',result,'SP71',48)
result=processMissing('SP72/Analysis/SP72_NaIV_processed3.csv','SP72/Analysis/SP72_missing.csv',result,'SP72',48)
result=processMissing('SP73/Analysis/SP73_NaIV_processed3.csv','SP73/Analysis/SP73_missing.csv',result,'SP73',48)
result=processMissing('SP74/Analysis/SP74_NaIV_processed3.csv','SP74/Analysis/SP74_missing.csv',result,'SP74',48)
result=processMissing('SP77/Analysis/SP77_NaIV_processed3.csv','SP77/Analysis/SP77_missing.csv',result,'SP77',48)
result=processMissing('SP78/Analysis/SP78_NaIV_processed3.csv','SP78/Analysis/SP78_missing.csv',result,'SP78',48)
result=processMissing('SP79/Analysis/SP79_NaIV_processed3.csv','SP79/Analysis/SP79_missing.csv',result,'SP79',48)
result=processMissing('SP80/Analysis/SP80_NaIV_processed3.csv','SP80/Analysis/SP80_missing.csv',result,'SP80',48)
result=processMissing('SP81/Analysis/SP81_NaIV_processed3.csv','SP81/Analysis/SP81_missing.csv',result,'SP81',48)

allAct=result[[1]]
empty=result[[2]]
allAct=allAct[!allAct$Type=='start',] #remove "start line"
allAct$PeakInclude=as.logical(allAct$PeakInclude)
allAct$VHalfActInclude=as.logical(allAct$VHalfActInclude)
allAct[is.na(allAct$Capacitance),'PeakInclude']=FALSE
allAct=fixSP6789names(allAct)
write.csv(allAct, 'CombinedAnalysis/SP5-81_allCells_act.csv',row.names=FALSE)
write.csv(empty,'CombinedAnalysis/SP5-81_missingCells.csv',row.names=FALSE)



##### STEP 1B: Combine SPs together (voltage of inactivation) #####
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
#result=processMissing('SP1/Analysis/SP1_Inac_processed2.csv','SP1/Analysis/SP1_missing.csv',result,'SP1',61)
#result=processMissing('SP2/Analysis/SP2_Inac_processed2.csv','SP2/Analysis/SP2_missing.csv',result,'SP2',61)
#result=processMissing('SP3/Analysis/SP3_Inac_processed2.csv','SP3/Analysis/SP3_missing.csv',result,'SP3',61)
#result=processMissing('SP4/Analysis/SP4_Inac_processed2.csv','SP4/Analysis/SP4_missing.csv',result,'SP4',61)
#result=processMissing('SP5/Analysis/SP5_Inac_processed2.csv','SP5/Analysis/SP5_missing.csv',result,'SP5',61)
result=processMissing('SP6+7/Analysis/SP6+7_Inac_processed2.csv','SP6+7/Analysis/SP6+7_missing.csv',result,'SP6+7',61)
result=processMissing('SP8+9/Analysis/SP8+9_Inac_processed2.csv','SP8+9/Analysis/SP8+9_missing.csv',result,'SP8+9',61)
result=processMissing('SP10/Analysis/SP10_Inac_processed2.csv','SP10/Analysis/SP10_missing.csv',result,'SP10',61)
result=processMissing('SP11/Analysis/SP11_Inac_processed2.csv','SP11/Analysis/SP11_missing.csv',result,'SP11',61)
result=processMissing('SP12/Analysis/SP12_Inac_processed2.csv','SP12/Analysis/SP12_missing.csv',result,'SP12',61)
result=processMissing('SP13/Analysis/SP13_Inac_processed2.csv','SP13/Analysis/SP13_missing.csv',result,'SP13',61)
result=processMissing('SP14/Analysis/SP14_Inac_processed2.csv','SP14/Analysis/SP14_missing.csv',result,'SP14',61)
result=processMissing('SP15/Analysis/SP15_Inac_processed2.csv','SP15/Analysis/SP15_missing.csv',result,'SP15',61)
result=processMissing('SP17/Analysis/SP17_Inac_processed2.csv','SP17/Analysis/SP17_missing.csv',result,'SP17',61)
result=processMissing('SP18/Analysis/SP18_Inac_processed2.csv','SP18/Analysis/SP18_missing.csv',result,'SP18',61)
result=processMissing('SP19/Analysis/SP19_Inac_processed2.csv','SP19/Analysis/SP19_missing.csv',result,'SP19',61)
result=processMissing('SP20/Analysis/SP20_Inac_processed2.csv','SP20/Analysis/SP20_missing.csv',result,'SP20',61)
result=processMissing('SP21/Analysis/SP21_Inac_processed2.csv','SP21/Analysis/SP21_missing.csv',result,'SP21',61)
result=processMissing('SP22/Analysis/SP22_Inac_processed2.csv','SP22/Analysis/SP22_missing.csv',result,'SP22',61)
result=processMissing('SP24/Analysis/SP24_Inac_processed2.csv','SP24/Analysis/SP24_missing.csv',result,'SP24',61)
result=processMissing('SP26/Analysis/SP26_Inac_processed2.csv','SP26/Analysis/SP26_missing.csv',result,'SP26',61)
result=processMissing('SP27/Analysis/SP27_Inac_processed2.csv','SP27/Analysis/SP27_missing.csv',result,'SP27',61)
result=processMissing('SP28/Analysis/SP28_Inac_processed2.csv','SP28/Analysis/SP28_missing.csv',result,'SP28',61)
result=processMissing('SP31/Analysis/SP31_Inac_processed2.csv','SP31/Analysis/SP31_missing.csv',result,'SP31',61)
result=processMissing('SP32/Analysis/SP32_Inac_processed2.csv','SP32/Analysis/SP32_missing.csv',result,'SP32',61)
result=processMissing('SP35/Analysis/SP35_Inac_processed2.csv','SP35/Analysis/SP35_missing.csv',result,'SP35',61)
result=processMissing('SP36/Analysis/SP36_Inac_processed2.csv','SP36/Analysis/SP36_missing.csv',result,'SP36',61)
result=processMissing('SP37/Analysis/SP37_Inac_processed2.csv','SP37/Analysis/SP37_missing.csv',result,'SP37',61)
result=processMissing('SP38/Analysis/SP38_Inac_processed2.csv','SP38/Analysis/SP38_missing.csv',result,'SP38',61)
result=processMissing('SP39/Analysis/SP39_Inac_processed2.csv','SP39/Analysis/SP39_missing.csv',result,'SP39',61)
result=processMissing('SP41/Analysis/SP41_Inac_processed2.csv','SP41/Analysis/SP41_missing.csv',result,'SP41',61)
result=processMissing('SP43/Analysis/SP43_Inac_processed2.csv','SP43/Analysis/SP43_missing.csv',result,'SP43',61)
result=processMissing('SP44/Analysis/SP44_Inac_processed2.csv','SP44/Analysis/SP44_missing.csv',result,'SP44',61)
result=processMissing('SP46/Analysis/SP46_Inac_processed2.csv','SP46/Analysis/SP46_missing.csv',result,'SP46',61)
result=processMissing('SP48/Analysis/SP48_Inac_processed2.csv','SP48/Analysis/SP48_missing.csv',result,'SP48',61)
result=processMissing('SP49/Analysis/SP49_Inac_processed2.csv','SP49/Analysis/SP49_missing.csv',result,'SP49',61)
result=processMissing('SP50/Analysis/SP50_Inac_processed2.csv','SP50/Analysis/SP50_missing.csv',result,'SP50',61)
result=processMissing('SP51/Analysis/SP51_Inac_processed2.csv','SP51/Analysis/SP51_missing.csv',result,'SP51',61)
result=processMissing('SP61/Analysis/SP61_Inac_processed2.csv','SP61/Analysis/SP61_missing.csv',result,'SP61',61)
result=processMissing('SP63/Analysis/SP63_Inac_processed2.csv','SP63/Analysis/SP63_missing.csv',result,'SP63',61)
result=processMissing('SP65/Analysis/SP65_Inac_processed2.csv','SP65/Analysis/SP65_missing.csv',result,'SP65',61)
result=processMissing('SP67/Analysis/SP67_Inac_processed2.csv','SP67/Analysis/SP67_missing.csv',result,'SP67',61)
result=processMissing('SP68/Analysis/SP68_Inac_processed2.csv','SP68/Analysis/SP68_missing.csv',result,'SP68',61)
result=processMissing('SP69/Analysis/SP69_Inac_processed2.csv','SP69/Analysis/SP69_missing.csv',result,'SP69',61)
result=processMissing('SP71/Analysis/SP71_Inac_processed2.csv','SP71/Analysis/SP71_missing.csv',result,'SP71',61)
result=processMissing('SP72/Analysis/SP72_Inac_processed2.csv','SP72/Analysis/SP72_missing.csv',result,'SP72',61)
result=processMissing('SP73/Analysis/SP73_Inac_processed2.csv','SP73/Analysis/SP73_missing.csv',result,'SP73',61)
result=processMissing('SP74/Analysis/SP74_Inac_processed2.csv','SP74/Analysis/SP74_missing.csv',result,'SP74',61)
result=processMissing('SP77/Analysis/SP77_Inac_processed2.csv','SP77/Analysis/SP77_missing.csv',result,'SP77',61)
result=processMissing('SP78/Analysis/SP78_Inac_processed2.csv','SP78/Analysis/SP78_missing.csv',result,'SP78',61)
result=processMissing('SP79/Analysis/SP79_Inac_processed2.csv','SP79/Analysis/SP79_missing.csv',result,'SP79',61)
result=processMissing('SP80/Analysis/SP80_Inac_processed2.csv','SP80/Analysis/SP80_missing.csv',result,'SP80',61)
result=processMissing('SP81/Analysis/SP81_Inac_processed2.csv','SP81/Analysis/SP81_missing.csv',result,'SP81',61)


allInact=result[[1]]
empty=result[[2]]
allInact=allInact[!allInact$Type=='start',] #remove "start line"
allInact$PeakInclude=as.logical(allInact$PeakInclude)
allInact[is.na(allInact$Capacitance),'PeakInclude']=FALSE
allInact$VHalfInactInclude=as.logical(allInact$VHalfInactInclude)
allInact=fixSP6789names(allInact)
write.csv(allInact, 'CombinedAnalysis/SP5-81_allCells_inact.csv',row.names=FALSE,quote=FALSE)


##### STEP 1C: Combine SPs together (late current) #####

empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
#result=processMissing('SP1/Analysis/SP1_Late_processed2.csv','SP1/Analysis/SP1_missing.csv',result,'SP1',45)
#result=processMissing('SP2/Analysis/SP2_Late_processed2.csv','SP2/Analysis/SP2_missing.csv',result,'SP2',45)
#result=processMissing('SP3/Analysis/SP3_Late_processed2.csv','SP3/Analysis/SP3_missing.csv',result,'SP3',45)
#result=processMissing('SP4/Analysis/SP4_Late_processed2.csv','SP4/Analysis/SP4_missing.csv',result,'SP4',45)
#result=processMissing('SP5/Analysis/SP5_Late_processed2.csv','SP5/Analysis/SP5_missing.csv',result,'SP5',45)
result=processMissing('SP6+7/Analysis/SP6+7_Late_processed2.csv','SP6+7/Analysis/SP6+7_missing.csv',result,'SP6+7',45)
result=processMissing('SP8+9/Analysis/SP8+9_Late_processed2.csv','SP8+9/Analysis/SP8+9_missing.csv',result,'SP8+9',45)
result=processMissing('SP10/Analysis/SP10_Late_processed2.csv','SP10/Analysis/SP10_missing.csv',result,'SP10',45)
result=processMissing('SP11/Analysis/SP11_Late_processed2.csv','SP11/Analysis/SP11_missing.csv',result,'SP11',45)
result=processMissing('SP12/Analysis/SP12_Late_processed2.csv','SP12/Analysis/SP12_missing.csv',result,'SP12',45)
result=processMissing('SP13/Analysis/SP13_Late_processed2.csv','SP13/Analysis/SP13_missing.csv',result,'SP13',45)
result=processMissing('SP14/Analysis/SP14_Late_processed2.csv','SP14/Analysis/SP14_missing.csv',result,'SP14',45)
result=processMissing('SP15/Analysis/SP15_Late_processed2.csv','SP15/Analysis/SP15_missing.csv',result,'SP15',45)
result=processMissing('SP17/Analysis/SP17_Late_processed2.csv','SP17/Analysis/SP17_missing.csv',result,'SP17',45)
result=processMissing('SP18/Analysis/SP18_Late_processed2.csv','SP18/Analysis/SP18_missing.csv',result,'SP18',45)
result=processMissing('SP19/Analysis/SP19_Late_processed2.csv','SP19/Analysis/SP19_missing.csv',result,'SP19',45)
result=processMissing('SP20/Analysis/SP20_Late_processed2.csv','SP20/Analysis/SP20_missing.csv',result,'SP20',45)
result=processMissing('SP21/Analysis/SP21_Late_processed2.csv','SP21/Analysis/SP21_missing.csv',result,'SP21',45)
result=processMissing('SP22/Analysis/SP22_Late_processed2.csv','SP22/Analysis/SP22_missing.csv',result,'SP22',45)
result=processMissing('SP24/Analysis/SP24_Late_processed2.csv','SP24/Analysis/SP24_missing.csv',result,'SP24',45)
result=processMissing('SP26/Analysis/SP26_Late_processed2.csv','SP26/Analysis/SP26_missing.csv',result,'SP26',45)
result=processMissing('SP27/Analysis/SP27_Late_processed2.csv','SP27/Analysis/SP27_missing.csv',result,'SP27',45)
result=processMissing('SP28/Analysis/SP28_Late_processed2.csv','SP28/Analysis/SP28_missing.csv',result,'SP28',45)
#result=processMissing('SP30/Analysis/SP30_Late_processed2.csv','SP30/Analysis/SP30_missing.csv',result,'SP30',45)
result=processMissing('SP31/Analysis/SP31_Late_processed2.csv','SP31/Analysis/SP31_missing.csv',result,'SP31',45)
result=processMissing('SP32/Analysis/SP32_Late_processed2.csv','SP32/Analysis/SP32_missing.csv',result,'SP32',45)
result=processMissing('SP35/Analysis/SP35_Late_processed2.csv','SP35/Analysis/SP35_missing.csv',result,'SP35',45)
result=processMissing('SP36/Analysis/SP36_Late_processed2.csv','SP36/Analysis/SP36_missing.csv',result,'SP36',45)
result=processMissing('SP37/Analysis/SP37_Late_processed2.csv','SP37/Analysis/SP37_missing.csv',result,'SP37',45)
result=processMissing('SP38/Analysis/SP38_Late_processed2.csv','SP38/Analysis/SP38_missing.csv',result,'SP38',45)
result=processMissing('SP39/Analysis/SP39_Late_processed2.csv','SP39/Analysis/SP39_missing.csv',result,'SP39',45)
result=processMissing('SP41/Analysis/SP41_Late_processed2.csv','SP41/Analysis/SP41_missing.csv',result,'SP41',45)
result=processMissing('SP43/Analysis/SP43_Late_processed2.csv','SP43/Analysis/SP43_missing.csv',result,'SP43',45)
result=processMissing('SP44/Analysis/SP44_Late_processed2.csv','SP44/Analysis/SP44_missing.csv',result,'SP44',45)
result=processMissing('SP46/Analysis/SP46_Late_processed2.csv','SP46/Analysis/SP46_missing.csv',result,'SP46',45)
result=processMissing('SP48/Analysis/SP48_Late_processed2.csv','SP48/Analysis/SP48_missing.csv',result,'SP48',45)
result=processMissing('SP49/Analysis/SP49_Late_processed2.csv','SP49/Analysis/SP49_missing.csv',result,'SP49',45)
result=processMissing('SP50/Analysis/SP50_Late_processed2.csv','SP50/Analysis/SP50_missing.csv',result,'SP50',45)
result=processMissing('SP51/Analysis/SP51_Late_processed2.csv','SP51/Analysis/SP51_missing.csv',result,'SP51',45)
result=processMissing('SP61/Analysis/SP61_Late_processed2.csv','SP61/Analysis/SP61_missing.csv',result,'SP61',45)
result=processMissing('SP63/Analysis/SP63_Late_processed2.csv','SP63/Analysis/SP63_missing.csv',result,'SP63',45)
result=processMissing('SP72/Analysis/SP72_Late_processed2.csv','SP72/Analysis/SP72_missing.csv',result,'SP72',45)
result=processMissing('SP73/Analysis/SP73_Late_processed2.csv','SP73/Analysis/SP73_missing.csv',result,'SP73',45)
result=processMissing('SP74/Analysis/SP74_Late_processed2.csv','SP74/Analysis/SP74_missing.csv',result,'SP74',45)
result=processMissing('SP77/Analysis/SP77_Late_processed2.csv','SP77/Analysis/SP77_missing.csv',result,'SP77',45)
result=processMissing('SP78/Analysis/SP78_Late_processed2.csv','SP78/Analysis/SP78_missing.csv',result,'SP78',45)
result=processMissing('SP79/Analysis/SP79_Late_processed2.csv','SP79/Analysis/SP79_missing.csv',result,'SP79',45)
result=processMissing('SP80/Analysis/SP80_Late_processed2.csv','SP80/Analysis/SP80_missing.csv',result,'SP80',45)
result=processMissing('SP81/Analysis/SP81_Late_processed2.csv','SP81/Analysis/SP81_missing.csv',result,'SP81',45)

# New late current protocol - first load old ones together, then insert proxy columns, then merge with V2 protocol runs and save
# Create a proxy data frame to align names correctly - will disregard incomplete results as NAs 

allLate=result[[1]]
empty=result[[2]]
allLate=allLate[!allLate$Type=='start',] #remove "start line"
proxy <- data.frame(matrix(ncol = 11, nrow = nrow(allLate)))
colnames(proxy) <- c("Peak_30Pre_Late", "Peak_30Post_Late", "Peak_30_Late", "Late200_30Pre_Late", "Late200_30Post_Late", "Late200_30_Late", "Late200_30Ratio_Late", "Late50_30Pre_Late", "Late50_30Post_Late", "Late50_30_Late", "Late50_30Ratio_Late")
allLate2 <- cbind(allLate[,1:45], proxy, Include = allLate[,46])
allLate2$PeakInclude=as.logical(allLate2$PeakInclude)


# Now read in SP late protocols from those using the -15 and the -30 pulses

empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
result=processMissing('SP65/Analysis/SP65_Late_processed2.csv','SP65/Analysis/SP65_missing.csv',result,'SP65',56)
result=processMissing('SP67/Analysis/SP67_Late_processed2.csv','SP67/Analysis/SP67_missing.csv',result,'SP67',56)
result=processMissing('SP68/Analysis/SP68_Late_processed2.csv','SP68/Analysis/SP68_missing.csv',result,'SP68',56)
result=processMissing('SP69/Analysis/SP69_Late_processed2.csv','SP69/Analysis/SP69_missing.csv',result,'SP69',56)
result=processMissing('SP71/Analysis/SP71_Late_processed2.csv','SP71/Analysis/SP71_missing.csv',result,'SP71',56)

# Now merge with the old protocol

allLate3=result[[1]]
allLate3=allLate3[!allLate3$Type=='start',] #remove "start line"
allLate3$PeakInclude=as.logical(allLate3$PeakInclude)
allLate4 <- rbind(allLate2, allLate3)
allLate4[is.na(allLate4$Capacitance),'PeakInclude']=FALSE
allLate4=fixSP6789names(allLate4)
write.csv(allLate4, 'CombinedAnalysis/SP5-81_allCells_late.csv',row.names=FALSE,quote=FALSE)


##### STEP 1D: Combine SPs together (Recovery From Inactivation) ##### 
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
#result=processMissing('SP1/Analysis/SP1_RecInact_processed2.csv','SP1/Analysis/SP1_missing.csv',result,'SP1',95)
#result=processMissing('SP2/Analysis/SP2_RecInact_processed2.csv','SP2/Analysis/SP2_missing.csv',result,'SP2',95)
#result=processMissing('SP3/Analysis/SP3_RecInact_processed2.csv','SP3/Analysis/SP3_missing.csv',result,'SP3',95)
#result=processMissing('SP4/Analysis/SP4_RecInact_processed2.csv','SP4/Analysis/SP4_missing.csv',result,'SP4',95)
#result=processMissing('SP5/Analysis/SP5_RecInact_processed2.csv','SP5/Analysis/SP5_missing.csv',result,'SP5',95)
result=processMissing('SP6+7/Analysis/SP6+7_RecInact_processed2.csv','SP6+7/Analysis/SP6+7_missing.csv',result,'SP6+7',95)
result=processMissing('SP8+9/Analysis/SP8+9_RecInact_processed2.csv','SP8+9/Analysis/SP8+9_missing.csv',result,'SP8+9',95)
result=processMissing('SP10/Analysis/SP10_RecInact_processed2.csv','SP10/Analysis/SP10_missing.csv',result,'SP10',95)
result=processMissing('SP11/Analysis/SP11_RecInact_processed2.csv','SP11/Analysis/SP11_missing.csv',result,'SP11',95)
result=processMissing('SP12/Analysis/SP12_RecInact_processed2.csv','SP12/Analysis/SP12_missing.csv',result,'SP12',95)
result=processMissing('SP13/Analysis/SP13_RecInact_processed2.csv','SP13/Analysis/SP13_missing.csv',result,'SP13',95)
result=processMissing('SP14/Analysis/SP14_RecInact_processed2.csv','SP14/Analysis/SP14_missing.csv',result,'SP14',95)
result=processMissing('SP15/Analysis/SP15_RecInact_processed2.csv','SP15/Analysis/SP15_missing.csv',result,'SP15',95)
result=processMissing('SP17/Analysis/SP17_RecInact_processed2.csv','SP17/Analysis/SP17_missing.csv',result,'SP17',95)
result=processMissing('SP18/Analysis/SP18_RecInact_processed2.csv','SP18/Analysis/SP18_missing.csv',result,'SP18',95)
result=processMissing('SP19/Analysis/SP19_RecInact_processed2.csv','SP19/Analysis/SP19_missing.csv',result,'SP19',95)
result=processMissing('SP20/Analysis/SP20_RecInact_processed2.csv','SP20/Analysis/SP20_missing.csv',result,'SP20',95)
result=processMissing('SP21/Analysis/SP21_RecInact_processed2.csv','SP21/Analysis/SP21_missing.csv',result,'SP21',95)
result=processMissing('SP22/Analysis/SP22_RecInact_processed2.csv','SP22/Analysis/SP22_missing.csv',result,'SP22',95)
result=processMissing('SP24/Analysis/SP24_RecInact_processed2.csv','SP24/Analysis/SP24_missing.csv',result,'SP24',95)
result=processMissing('SP26/Analysis/SP26_RecInact_processed2.csv','SP26/Analysis/SP26_missing.csv',result,'SP26',95)
result=processMissing('SP27/Analysis/SP27_RecInact_processed2.csv','SP27/Analysis/SP27_missing.csv',result,'SP27',95)
result=processMissing('SP28/Analysis/SP28_RecInact_processed2.csv','SP28/Analysis/SP28_missing.csv',result,'SP28',95)
result=processMissing('SP31/Analysis/SP31_RecInact_processed2.csv','SP31/Analysis/SP31_missing.csv',result,'SP31',95)
result=processMissing('SP32/Analysis/SP32_RecInact_processed2.csv','SP32/Analysis/SP32_missing.csv',result,'SP32',95)
result=processMissing('SP35/Analysis/SP35_RecInact_processed2.csv','SP35/Analysis/SP35_missing.csv',result,'SP35',95)
result=processMissing('SP36/Analysis/SP36_RecInact_processed2.csv','SP36/Analysis/SP36_missing.csv',result,'SP36',95)
result=processMissing('SP37/Analysis/SP37_RecInact_processed2.csv','SP37/Analysis/SP37_missing.csv',result,'SP37',95)
result=processMissing('SP38/Analysis/SP38_RecInact_processed2.csv','SP38/Analysis/SP38_missing.csv',result,'SP38',95)
result=processMissing('SP39/Analysis/SP39_RecInact_processed2.csv','SP39/Analysis/SP39_missing.csv',result,'SP39',95)
result=processMissing('SP41/Analysis/SP41_RecInact_processed2.csv','SP41/Analysis/SP41_missing.csv',result,'SP41',95)
result=processMissing('SP43/Analysis/SP43_RecInact_processed2.csv','SP43/Analysis/SP43_missing.csv',result,'SP43',95)
result=processMissing('SP44/Analysis/SP44_RecInact_processed2.csv','SP44/Analysis/SP44_missing.csv',result,'SP44',95)
result=processMissing('SP46/Analysis/SP46_RecInact_processed2.csv','SP46/Analysis/SP46_missing.csv',result,'SP46',95)
result=processMissing('SP48/Analysis/SP48_RecInact_processed2.csv','SP48/Analysis/SP48_missing.csv',result,'SP48',95)
result=processMissing('SP49/Analysis/SP49_RecInact_processed2.csv','SP49/Analysis/SP49_missing.csv',result,'SP49',95)
result=processMissing('SP50/Analysis/SP50_RecInact_processed2.csv','SP50/Analysis/SP50_missing.csv',result,'SP50',95)
result=processMissing('SP51/Analysis/SP51_RecInact_processed2.csv','SP51/Analysis/SP51_missing.csv',result,'SP51',95)
result=processMissing('SP61/Analysis/SP61_RecInact_processed2.csv','SP61/Analysis/SP61_missing.csv',result,'SP61',95)
result=processMissing('SP63/Analysis/SP63_RecInact_processed2.csv','SP63/Analysis/SP63_missing.csv',result,'SP63',95)
result=processMissing('SP65/Analysis/SP65_RecInact_processed2.csv','SP65/Analysis/SP65_missing.csv',result,'SP65',95)
result=processMissing('SP67/Analysis/SP67_RecInact_processed2.csv','SP67/Analysis/SP67_missing.csv',result,'SP67',95)
result=processMissing('SP68/Analysis/SP68_RecInact_processed2.csv','SP68/Analysis/SP68_missing.csv',result,'SP68',95)
result=processMissing('SP69/Analysis/SP69_RecInact_processed2.csv','SP69/Analysis/SP69_missing.csv',result,'SP69',95)
result=processMissing('SP71/Analysis/SP71_RecInact_processed2.csv','SP71/Analysis/SP71_missing.csv',result,'SP71',95)
result=processMissing('SP72/Analysis/SP72_RecInact_processed2.csv','SP72/Analysis/SP72_missing.csv',result,'SP72',95)
result=processMissing('SP73/Analysis/SP73_RecInact_processed2.csv','SP73/Analysis/SP73_missing.csv',result,'SP73',95)
result=processMissing('SP74/Analysis/SP74_RecInact_processed2.csv','SP74/Analysis/SP74_missing.csv',result,'SP74',95)
result=processMissing('SP77/Analysis/SP77_RecInact_processed2.csv','SP77/Analysis/SP77_missing.csv',result,'SP77',95)
result=processMissing('SP78/Analysis/SP78_RecInact_processed2.csv','SP78/Analysis/SP78_missing.csv',result,'SP78',95)
result=processMissing('SP79/Analysis/SP79_RecInact_processed2.csv','SP79/Analysis/SP79_missing.csv',result,'SP79',95)
result=processMissing('SP80/Analysis/SP80_RecInact_processed2.csv','SP80/Analysis/SP80_missing.csv',result,'SP80',95)
result=processMissing('SP81/Analysis/SP81_RecInact_processed2.csv','SP81/Analysis/SP81_missing.csv',result,'SP81',95)


allRFI=result[[1]]
empty=result[[2]]
allRFI=allRFI[!allRFI$Type=='start',] #remove "start line"
allRFI$PeakInclude=as.logical(allRFI$PeakInclude)
allRFI[is.na(allRFI$Capacitance),'PeakInclude']=FALSE
allRFI=fixSP6789names(allRFI)
write.csv(allRFI, 'CombinedAnalysis/SP5-81_allCells_rfi.csv',row.names=FALSE,quote=FALSE)

##### STEP 1E: Combine SPs together (Inactivation time) #####
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
#result=processMissing('SP1/Analysis/SP1_Inacttime_processed.csv','SP1/Analysis/SP1_missing.csv',result,'SP1',53)
#result=processMissing('SP2/Analysis/SP2_Inacttime_processed.csv','SP2/Analysis/SP2_missing.csv',result,'SP2',53)
#result=processMissing('SP3/Analysis/SP3_Inacttime_processed.csv','SP3/Analysis/SP3_missing.csv',result,'SP3',53)
#result=processMissing('SP4/Analysis/SP4_Inacttime_processed.csv','SP4/Analysis/SP4_missing.csv',result,'SP4',53)
#result=processMissing('SP5/Analysis/SP5_Inacttime_processed.csv','SP5/Analysis/SP5_missing.csv',result,'SP5',53)
result=processMissing('SP6+7/Analysis/SP6+7_Inacttime_processed.csv','SP6+7/Analysis/SP6+7_missing.csv',result,'SP6+7',53)
result=processMissing('SP8+9/Analysis/SP8+9_Inacttime_processed.csv','SP8+9/Analysis/SP8+9_missing.csv',result,'SP8+9',53)
result=processMissing('SP10/Analysis/SP10_Inacttime_processed.csv','SP10/Analysis/SP10_missing.csv',result,'SP10',53)
result=processMissing('SP11/Analysis/SP11_Inacttime_processed.csv','SP11/Analysis/SP11_missing.csv',result,'SP11',53)
result=processMissing('SP12/Analysis/SP12_Inacttime_processed.csv','SP12/Analysis/SP12_missing.csv',result,'SP12',53)
result=processMissing('SP13/Analysis/SP13_Inacttime_processed.csv','SP13/Analysis/SP13_missing.csv',result,'SP13',53)
result=processMissing('SP14/Analysis/SP14_Inacttime_processed.csv','SP14/Analysis/SP14_missing.csv',result,'SP14',53)
result=processMissing('SP15/Analysis/SP15_Inacttime_processed.csv','SP15/Analysis/SP15_missing.csv',result,'SP15',53)
result=processMissing('SP17/Analysis/SP17_Inacttime_processed.csv','SP17/Analysis/SP17_missing.csv',result,'SP17',53)
result=processMissing('SP18/Analysis/SP18_Inacttime_processed.csv','SP18/Analysis/SP18_missing.csv',result,'SP18',53)
result=processMissing('SP19/Analysis/SP19_Inacttime_processed.csv','SP19/Analysis/SP19_missing.csv',result,'SP19',53)
result=processMissing('SP20/Analysis/SP20_Inacttime_processed.csv','SP20/Analysis/SP20_missing.csv',result,'SP20',53)
result=processMissing('SP21/Analysis/SP21_Inacttime_processed.csv','SP21/Analysis/SP21_missing.csv',result,'SP21',53)
result=processMissing('SP22/Analysis/SP22_Inacttime_processed.csv','SP22/Analysis/SP22_missing.csv',result,'SP22',53)
result=processMissing('SP24/Analysis/SP24_Inacttime_processed.csv','SP24/Analysis/SP24_missing.csv',result,'SP24',53)
result=processMissing('SP27/Analysis/SP27_Inacttime_processed.csv','SP27/Analysis/SP27_missing.csv',result,'SP27',53)
result=processMissing('SP28/Analysis/SP28_Inacttime_processed.csv','SP28/Analysis/SP28_missing.csv',result,'SP28',53)
result=processMissing('SP31/Analysis/SP31_Inacttime_processed.csv','SP31/Analysis/SP31_missing.csv',result,'SP31',53)
result=processMissing('SP32/Analysis/SP32_Inacttime_processed.csv','SP32/Analysis/SP32_missing.csv',result,'SP32',53)
result=processMissing('SP35/Analysis/SP35_Inacttime_processed.csv','SP35/Analysis/SP35_missing.csv',result,'SP35',53)
result=processMissing('SP36/Analysis/SP36_Inacttime_processed.csv','SP36/Analysis/SP36_missing.csv',result,'SP36',53)
result=processMissing('SP37/Analysis/SP37_Inacttime_processed.csv','SP37/Analysis/SP37_missing.csv',result,'SP37',53)
result=processMissing('SP38/Analysis/SP38_Inacttime_processed.csv','SP38/Analysis/SP38_missing.csv',result,'SP38',53)
result=processMissing('SP39/Analysis/SP39_Inacttime_processed.csv','SP39/Analysis/SP39_missing.csv',result,'SP39',53)
result=processMissing('SP41/Analysis/SP41_Inacttime_processed.csv','SP41/Analysis/SP41_missing.csv',result,'SP41',53)
result=processMissing('SP43/Analysis/SP43_Inacttime_processed.csv','SP43/Analysis/SP43_missing.csv',result,'SP43',53)
result=processMissing('SP46/Analysis/SP46_Inacttime_processed.csv','SP46/Analysis/SP46_missing.csv',result,'SP46',53)
result=processMissing('SP48/Analysis/SP48_Inacttime_processed.csv','SP48/Analysis/SP48_missing.csv',result,'SP48',53)
result=processMissing('SP49/Analysis/SP49_Inacttime_processed.csv','SP49/Analysis/SP49_missing.csv',result,'SP49',53)
result=processMissing('SP50/Analysis/SP50_Inacttime_processed.csv','SP50/Analysis/SP50_missing.csv',result,'SP50',53)
result=processMissing('SP51/Analysis/SP51_Inacttime_processed.csv','SP51/Analysis/SP51_missing.csv',result,'SP51',53)
result=processMissing('SP61/Analysis/SP61_Inacttime_processed.csv','SP61/Analysis/SP61_missing.csv',result,'SP61',53)
result=processMissing('SP63/Analysis/SP63_Inacttime_processed.csv','SP63/Analysis/SP63_missing.csv',result,'SP63',53)
result=processMissing('SP65/Analysis/SP65_Inacttime_processed.csv','SP65/Analysis/SP65_missing.csv',result,'SP65',53)
result=processMissing('SP67/Analysis/SP67_Inacttime_processed.csv','SP67/Analysis/SP67_missing.csv',result,'SP67',53)
result=processMissing('SP68/Analysis/SP68_Inacttime_processed.csv','SP68/Analysis/SP68_missing.csv',result,'SP68',53)
result=processMissing('SP69/Analysis/SP69_Inacttime_processed.csv','SP69/Analysis/SP69_missing.csv',result,'SP69',53)
result=processMissing('SP71/Analysis/SP71_Inacttime_processed.csv','SP71/Analysis/SP71_missing.csv',result,'SP71',53)
result=processMissing('SP72/Analysis/SP72_Inacttime_processed.csv','SP72/Analysis/SP72_missing.csv',result,'SP72',53)
result=processMissing('SP73/Analysis/SP73_Inacttime_processed.csv','SP73/Analysis/SP73_missing.csv',result,'SP73',53)
result=processMissing('SP74/Analysis/SP74_Inacttime_processed.csv','SP74/Analysis/SP74_missing.csv',result,'SP74',53)
result=processMissing('SP77/Analysis/SP77_Inacttime_processed.csv','SP77/Analysis/SP77_missing.csv',result,'SP77',53)
result=processMissing('SP78/Analysis/SP78_Inacttime_processed.csv','SP78/Analysis/SP78_missing.csv',result,'SP78',53)
result=processMissing('SP79/Analysis/SP79_Inacttime_processed.csv','SP79/Analysis/SP79_missing.csv',result,'SP79',53)
result=processMissing('SP80/Analysis/SP80_Inacttime_processed.csv','SP80/Analysis/SP80_missing.csv',result,'SP80',53)
result=processMissing('SP81/Analysis/SP81_Inacttime_processed.csv','SP81/Analysis/SP81_missing.csv',result,'SP81',53)

allInacttime=result[[1]]
empty=result[[2]]
allInacttime=allInacttime[!allInacttime$Type=='start',] #remove "start line"
allInacttime$PeakInclude=as.logical(allInacttime$PeakInclude)
allInacttime[is.na(allInacttime$Capacitance),'PeakInclude']=FALSE
allInacttime=fixSP6789names(allInacttime)
write.csv(allInacttime, 'CombinedAnalysis/SP5-81_allCells_inacttime.csv',row.names=FALSE,quote=FALSE)


##### STEP 1F: Combine SPs together (Peak current -90mV) #####

empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)

#result=processMissing('SP3/Analysis/SP3_Peak90_processed2.csv','SP3/Analysis/SP3_missing2.csv',result,'S3',61)
result=processMissing('SP8+9/Analysis/SP8+9_Peak90_processed2.csv','SP8+9/Analysis/SP8+9_missing2.csv',result,'SP8+9',61)
result=processMissing('SP11/Analysis/SP11_Peak90_processed2.csv','SP11/Analysis/SP11_missing2.csv',result,'SP11',61)
result=processMissing('SP14/Analysis/SP14_Peak90_processed2.csv','SP14/Analysis/SP14_missing2.csv',result,'SP14',61)
result=processMissing('SP15/Analysis/SP15_Peak90_processed2.csv','SP15/Analysis/SP15_missing2.csv',result,'SP15',61)
result=processMissing('SP17/Analysis/SP17_Peak90_processed2.csv','SP17/Analysis/SP17_missing2.csv',result,'SP17',61)
result=processMissing('SP19/Analysis/SP19_Peak90_processed2.csv','SP19/Analysis/SP19_missing2.csv',result,'SP19',61)
result=processMissing('SP36/Analysis/SP36_Peak90_processed2.csv','SP36/Analysis/SP36_missing2.csv',result,'SP36',61)
result=processMissing('SP37/Analysis/SP37_Peak90_processed2.csv','SP37/Analysis/SP37_missing2.csv',result,'SP37',61)
result=processMissing('SP38/Analysis/SP38_Peak90_processed2.csv','SP38/Analysis/SP38_missing2.csv',result,'SP38',61)
result=processMissing('SP39/Analysis/SP39_Peak90_processed2.csv','SP39/Analysis/SP39_missing2.csv',result,'SP39',61)
result=processMissing('SP41/Analysis/SP41_Peak90_processed2.csv','SP41/Analysis/SP41_missing2.csv',result,'SP41',61)
result=processMissing('SP43/Analysis/SP43_Peak90_processed2.csv','SP43/Analysis/SP43_missing2.csv',result,'SP43',61)
result=processMissing('SP46/Analysis/SP46_Peak90_processed2.csv','SP46/Analysis/SP46_missing2.csv',result,'SP46',61)

result=processMissing('SP48/Analysis/SP48_Peak90_processed2.csv','SP48/Analysis/SP48_missing2.csv',result,'SP48',61)
result=processMissing('SP49/Analysis/SP49_Peak90_processed2.csv','SP49/Analysis/SP49_missing2.csv',result,'SP49',61)
result=processMissing('SP51/Analysis/SP51_Peak90_processed2.csv','SP51/Analysis/SP51_missing2.csv',result,'SP51',61)
result=processMissing('SP61/Analysis/SP61_Peak90_processed2.csv','SP61/Analysis/SP61_missing2.csv',result,'SP61',61)
result=processMissing('SP63/Analysis/SP63_Peak90_processed2.csv','SP63/Analysis/SP63_missing2.csv',result,'SP63',61)
result=processMissing('SP65/Analysis/SP65_Peak90_processed2.csv','SP65/Analysis/SP65_missing2.csv',result,'SP65',61)
result=processMissing('SP67/Analysis/SP67_Peak90_processed2.csv','SP67/Analysis/SP67_missing2.csv',result,'SP67',61)
result=processMissing('SP68/Analysis/SP68_Peak90_processed2.csv','SP68/Analysis/SP68_missing2.csv',result,'SP68',61)
result=processMissing('SP69/Analysis/SP69_Peak90_processed2.csv','SP69/Analysis/SP69_missing2.csv',result,'SP69',61)
result=processMissing('SP71/Analysis/SP71_Peak90_processed2.csv','SP71/Analysis/SP71_missing2.csv',result,'SP71',61)
result=processMissing('SP72/Analysis/SP72_Peak90_processed2.csv','SP72/Analysis/SP72_missing2.csv',result,'SP72',61)
result=processMissing('SP73/Analysis/SP73_Peak90_processed2.csv','SP73/Analysis/SP73_missing2.csv',result,'SP73',61)
result=processMissing('SP74/Analysis/SP74_Peak90_processed2.csv','SP74/Analysis/SP74_missing2.csv',result,'SP74',61)
result=processMissing('SP77/Analysis/SP77_Peak90_processed2.csv','SP77/Analysis/SP77_missing2.csv',result,'SP77',61)
result=processMissing('SP78/Analysis/SP78_Peak90_processed2.csv','SP78/Analysis/SP78_missing2.csv',result,'SP78',61)
result=processMissing('SP79/Analysis/SP79_Peak90_processed2.csv','SP79/Analysis/SP79_missing2.csv',result,'SP79',61)
result=processMissing('SP80/Analysis/SP80_Peak90_processed2.csv','SP80/Analysis/SP80_missing2.csv',result,'SP80',61)
result=processMissing('SP81/Analysis/SP81_Peak90_processed2.csv','SP81/Analysis/SP81_missing2.csv',result,'SP81',61)


allPeak90=result[[1]]
empty=result[[2]]
allPeak90=allPeak90[!allPeak90$Type=='start',] #remove "start line"
allPeak90$PeakInclude=as.logical(allPeak90$PeakInclude)
allPeak90[is.na(allPeak90$Capacitance),'PeakInclude']=FALSE
allPeak90=fixSP6789names(allPeak90)

write.csv(allPeak90, '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/E1225K_Peak90.csv',row.names=FALSE)
write.csv(empty,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/E1225K_missingCells2.csv',row.names=FALSE)

#write.csv(allPeak90, '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_Peak90.csv',row.names=FALSE)
#write.csv(empty,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_missingCells2.csv',row.names=FALSE)




########## STEP 2: Analyze these combined datasets, perform outlier removal when appropriate ##########


# STEP 2A) Calculate -120 mV peak current (no sqrt transformation)
# STEP 2B) Calculate -120 mV peak current (with sqrt transformation)
# STEP 2C) Calculate -90 mV peak current (no sqrt transformation)
# STEP 2D) Calculate -90 mV peak current (with sqrt transformation)
# STEP 2E) Calculate Recovery from Inactivation
# STEP 2F) Calculate Inactivation time
# STEP 2G) Calculate Voltage of Activation
# STEP 2H) Calculate Voltage of Inactivation
# STEP 2I) Calculate Late current

##### STEP 2A: Calculate -120 mV peak current (no sqrt transformation) #####
i=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_act.csv',stringsAsFactors=FALSE)
empty=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_missingCells.csv',stringsAsFactors = FALSE)
peakSummary=summarizeByMutationGeneral(i,'PeakDensityNorm','PeakInclude',missingFrame=empty)
peakSummary2=orderMutFrame(peakSummary,'PeakDensityNormMean','WT',decreasing=TRUE)
peakSummary3=getMutDescriptionAndFilter(peakSummary2)
write.csv(peakSummary3,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peakNormData.csv',row.names=FALSE)

##### STEP 2B: Calculate -120 mV peak current (with sqrt transformation) #####
i=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_act.csv',stringsAsFactors=FALSE)
empty=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_missingCells.csv',stringsAsFactors = FALSE)
j = sqrt_transform_120(i) #creates 2 variables, PeakDensitySQRTraw and PeakDensityNormSQRT
exportSingleCell(j, '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/control_ACMG_variants.csv', "PeakInclude", '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peakNormDataSingle120_v2.csv', empty)
peakSummary=summarizeByMutationGeneral(j,'PeakDensityNormSQRT','PeakInclude',missingFrame=empty)
peakSummary2=orderMutFrame(peakSummary,'PeakDensityNormSQRTMean','WT',decreasing=TRUE)
peakSummary3=getMutDescriptionAndFilter(peakSummary2)
write.csv(peakSummary3,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peakNormSQRTData.csv',row.names=FALSE)

##### STEP 2C: Calculate -90 mV peak current (no sqrt transformation) [aka "78"] #####
i=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_Peak90.csv',stringsAsFactors=FALSE)
empty=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_missingCells2.csv',stringsAsFactors = FALSE)
peakSummary=summarizeByMutationGeneral(i,'Peak78DensityNorm','Peak78Include',missingFrame=empty)
peakSummary2=orderMutFrame(peakSummary,'Peak78DensityNormMean','WT',decreasing=TRUE)
peakSummary3=getMutDescriptionAndFilter(peakSummary2)
write.csv(peakSummary3,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peak90NormData.csv',row.names=FALSE)

##### STEP 2D: Calculate -90 mV peak current (with sqrt transformation) #####
i=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_Peak90.csv',stringsAsFactors=FALSE)
empty=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_missingCells2.csv',stringsAsFactors = FALSE)
j = sqrt_transform_78(i,Positive2Zero = TRUE)
exportSingleCell78(j, '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/control_ACMG_variants.csv', "Peak78Include", '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peakNormDataSingle90_v5.csv', empty)
peakSummary=summarizeByMutationGeneral(j,'Peak78DensityNormSQRT','Peak78Include',missingFrame=empty)
peakSummary2=orderMutFrame(peakSummary,'Peak78DensityNormSQRTMean','WT',decreasing=TRUE)
peakSummary3=getMutDescriptionAndFilter(peakSummary2)
write.csv(peakSummary3,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peak78NormSQRTData.csv',row.names=FALSE)

##### STEP 2E: Calculate Recovery from Inactivation  #####
i=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_rfi.csv',stringsAsFactors=FALSE)
i[i$Mutation=='WT2','Mutation']='WT'
i[i$Mutation=='WT3','Mutation']='WT'
i_wt_act=i[i$Mutation=='WT' & i$VHalfActInclude,]
i_wt_act <- filter(i,!str_detect(Plate,"^SP26") & !str_detect(Plate,"^SP40") & !str_detect(Plate,"^SP44") & Mutation=="WT" & VHalfActInclude == TRUE)
# Outlier removal
i_good=i[i$rfiInclude,]
i_wtgood=i_good[i_good$Mutation=='WT',]
i2_good=filterOutliers(i_good,'rfi50',i_wtgood$rfi50,3)
counts2 = filter_stats(i_good,i2_good, "rfi50")
i3_good <- fancy_outlier_removal(i_good,counts2, 'rfi50',3, i_wtgood)
exportSingleCell(i3_good, '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/control_ACMG_variants.csv', "rfiInclude", '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_RFI_all_singlewell.csv')
RFISummary=summarizeByMutationGeneral(i3_good,'rfi50','rfiInclude')
RFISummary2=orderMutFrame(RFISummary,'rfi50Mean','WT',decreasing=TRUE)
RFISummary3=RFISummary2[RFISummary2$numCellsrfi50>=5,]
RFISummary4=getMutDescriptionAndFilter(RFISummary3,reorderCols=FALSE)
RFISummary4=RFISummary4[,c('Mutation','MutDescription','rfi50Mean','rfi50SE','numCellsrfi50')]
RFISummary5=RFISummary4[with(RFISummary4,order(MutDescription,rfi50Mean,decreasing=TRUE)),]
RFISummary5$rfi50NormMean=RFISummary5$rfi50Mean-RFISummary5[RFISummary5$Mutation=='WT','rfi50Mean']
write.csv(RFISummary5,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_RFI_all.csv',row.names=FALSE)

##### STEP 2F: Calculate Inactivation time #####
i=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_inacttime.csv',stringsAsFactors=FALSE)
i[i$Mutation=='WT2','Mutation']='WT'
i[i$Mutation=='WT3','Mutation']='WT'
# Outlier removal 
i_good=i[!is.na(i$inacttimeError) & i$inacttimeError<.1 & i$inacttimeTau>0,]
i_good=i_good[!is.na(i_good$Mutation),]
i_wtgood=i_good[i_good$Mutation=='WT',]
i2_good=filterOutliers(i_good,'inacttimeTau',i_wtgood$inacttimeTau,3)
counts2 = filter_stats(i_good,i2_good, "inacttimeTau")
i3_good <- fancy_outlier_removal(i_good,counts2, 'inacttimeTau', 3, i_wtgood)
exportSingleCell(i3_good, '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/control_ACMG_variants.csv', "inacttimeInclude", '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_inacttime_all_singlewell.csv')
InacttimeSummary=summarizeByMutationGeneral(i3_good,'inacttimeTau','inacttimeInclude',roundUnits=3)
InacttimeSummary2=orderMutFrame(InacttimeSummary,'inacttimeTauMean','WT',decreasing=TRUE)
InacttimeSummary3=InacttimeSummary2[InacttimeSummary2$numCellsinacttimeTau>=5,]
InacttimeSummary4=getMutDescriptionAndFilter(InacttimeSummary3,reorderCols=FALSE)
InacttimeSummary4=InacttimeSummary4[,c('Mutation','MutDescription','inacttimeTauMean','inacttimeTauSE','numCellsinacttimeTau')]
InacttimeSummary5=InacttimeSummary4[with(InacttimeSummary4,order(MutDescription,inacttimeTauMean,decreasing=TRUE)),]
row.names(InacttimeSummary5)=1:nrow(InacttimeSummary5)
InacttimeSummary5$inacttimeTauNormMean=InacttimeSummary5$inacttimeTauMean-InacttimeSummary5[InacttimeSummary5$Mutation=='WT','inacttimeTauMean']
write.csv(InacttimeSummary5,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_inacttime_all.csv',row.names=FALSE)

##### STEP 2G: Calculate Voltage of Activation #####
i=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_act.csv',stringsAsFactors=FALSE)
i[i$Mutation=='WT2','Mutation']='WT'
i[i$Mutation=='WT3','Mutation']='WT'
i_good= filter(i,!str_detect(Plate,"^SP26") & !str_detect(Plate,"^SP40") & !str_detect(Plate,"^SP44") & VHalfActInclude == TRUE
                  & Peak > (-3e-09))  # Only including cells with peak current <3pA
i_wtgood=i_good[i_good$Mutation=='WT',]
i2_good=filterOutliers(i_good,'VHalfAct',i_wtgood$VHalfAct,3)
counts2 = filter_stats(i_good,i2_good, "VHalfAct")
i3_good <- fancy_outlier_removal(i_good,counts2, 'VHalfAct',3, i_wtgood)
exportSingleCell(i3_good, '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/control_ACMG_variants.csv', "VHalfActInclude", '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_act_all_singlewell.csv')
actSummary=summarizeByMutationGeneral(i3_good,'VHalfAct','VHalfActInclude')
actSummary2=orderMutFrame(actSummary,'VHalfActMean','WT',decreasing=TRUE) #500
actSummary3=actSummary2[actSummary2$numCellsVHalfAct>=5,] #438
actSummary4=getMutDescriptionAndFilter(actSummary3,reorderCols=FALSE)
actSummary4=actSummary4[,c('Mutation','MutDescription','VHalfActMean','VHalfActSE','numCellsVHalfAct')]
actSummary5=actSummary4[with(actSummary4,order(MutDescription,VHalfActMean,decreasing=TRUE)),]
row.names(actSummary5)=1:nrow(actSummary5)
actSummary5$VHalfActNormMean=actSummary5$VHalfActMean-actSummary5[actSummary5$Mutation=='WT','VHalfActMean']
write.csv(actSummary5,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_act_all.csv',row.names=FALSE)

##### STEP 2H: Calculate Voltage of Inactivation #####
i=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_inact.csv',stringsAsFactors=FALSE)
i[i$Mutation=='WT2','Mutation']='WT'
i[i$Mutation=='WT3','Mutation']='WT'
i_good=i[i$VHalfInactInclude,]
i_good=i_good[!i_good$Plate == "SP26A" & !i_good$Plate == "SP26B" & !i_good$Plate == "SP26C" & !i_good$Plate == "SP26D" &
                !i_good$Plate == "SP26E" & !i_good$Plate == "SP44A" & !i_good$Plate == "SP44B" & !i_good$Plate == "SP44C" &
                !i_good$Plate == "SP44D", ]
i_wtgood=i_good[i_good$Mutation=='WT',]
i2_good=filterOutliers(i_good,'VHalfInact',i_wtgood$VHalfInact,3)
counts2 = filter_stats(i_good,i2_good, "VHalfInact")
i3_good <- fancy_outlier_removal(i_good,counts2, 'VHalfInact',3, i_wtgood)
exportSingleCell(i3_good, '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/control_ACMG_variants.csv', "VHalfInactInclude", '~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_inact_all_singlewell.csv')
inactSummary=summarizeByMutationGeneral(i3_good,'VHalfInact','VHalfInactInclude')
inactSummary2=orderMutFrame(inactSummary,'VHalfInactMean','WT',decreasing=TRUE)
inactSummary3=inactSummary2[inactSummary2$numCellsVHalfInact>=5,]
inactSummary4=getMutDescriptionAndFilter(inactSummary3,reorderCols=FALSE)
inactSummary4=inactSummary4[,c('Mutation','MutDescription','VHalfInactMean','VHalfInactSE','numCellsVHalfInact')]
inactSummary5=inactSummary4[with(inactSummary4,order(MutDescription,VHalfInactMean,decreasing=TRUE)),]
row.names(inactSummary5)=1:nrow(inactSummary5)
inactSummary5$VHalfInactNormMean=inactSummary5$VHalfInactMean-inactSummary5[inactSummary5$Mutation=='WT','VHalfInactMean']
write.csv(inactSummary5,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_inact_all.csv',row.names=FALSE)


##### STEP 2I: Calculate late current #####
i=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allCells_late.csv',stringsAsFactors=FALSE)
i[i$Mutation=='WT2','Mutation']='WT'
i[i$Mutation=='WT3','Mutation']='WT'
i=i[i$GoodSeal_Late,]
i=i[!i$Plate == "SP26A" & !i$Plate == "SP26B" & !i$Plate == "SP26C" & !i$Plate == "SP26D" &
                !i$Plate == "SP26E" & !i$Plate == "SP44A" & !i$Plate == "SP44B" & !i$Plate == "SP44C" &
                !i$Plate == "SP44D" & !i$Plate == "SP36A" & !i$Plate == "SP36B", ]
i=i[!is.na(i$GoodSeal_Late), ]
i_wt_late=i[i$Mutation=='WT' & i$GoodSeal_Late,]
i_wt_late$Peak_Late2=10^12*i_wt_late$Peak_Late

# Plot Late current at 50 ms and 200 ms
plot(i_wt_late$Peak_Late2,i_wt_late$Late200Ratio_Late,xlab='Peak current (pA)',ylab='late ratio',ylim=c(-.05,.05),pch='.')
abline(h=0)
abline(v=-500)
plot(i_wt_late$Peak_Late2,i_wt_late$Late50Ratio_Late,xlab='Peak current (pA)',ylab='late ratio',ylim=c(-.05,.05),pch='.')
abline(h=0)
abline(v=-500)

i_wt_late=i_wt_late[i_wt_late$GoodSeal_Late & i_wt_late$Peak_Late2<(-500),]
i$LateInclude=i$GoodSeal_Late & i$Peak_Late<=(-500)*10^-12 & i$PeakPost_Late/i$PeakPre_Late<0.1 # the pre and post indicate pre/post tetracaine treatment  
i_good=i[i$LateInclude,]
i_wtgood=i_good[i_good$Mutation=='WT',]
i2_good=filterOutliers(i_good,'Late200Ratio_Late',i_wtgood$Late200Ratio_Late,3)
counts2 = filter_stats(i_good,i2_good, "Late200Ratio_Late")
i3_good <- fancy_outlier_removal(i_good,counts2, 'Late200Ratio_Late',3, i_wtgood)
lateSummary=summarizeByMutationGeneral(i3_good,'Late50Ratio_Late','LateInclude',roundUnits=4)
lateSummaryB=summarizeByMutationGeneral(i3_good,'Late200Ratio_Late','LateInclude',roundUnits=4)
lateSummaryC=summarizeByMutationGeneral(i3_good,'RampMinRatio_Late','LateInclude',roundUnits=4)
lateSummary2 <- lateSummary %>% merge(lateSummaryB,all=TRUE) %>% merge(lateSummaryC,all=TRUE) %>%
  orderMutFrame(.,'Late50Ratio_LateMean','WT',decreasing=TRUE)
lateSummary3=lateSummary2[lateSummary2$numCellsLate50Ratio_Late>=5,]
lateSummary4=getMutDescriptionAndFilter(lateSummary3,reorderCols=FALSE)
lateSummary4=lateSummary4[,c('Mutation','MutDescription','Late50Ratio_LateMean','Late50Ratio_LateSE','numCellsLate50Ratio_Late','RampMinRatio_LateMean','RampMinRatio_LateSE','numCellsRampMinRatio_Late','Late200Ratio_LateMean','Late200Ratio_LateSE','numCellsLate200Ratio_Late')]
lateSummary5=lateSummary4[with(lateSummary4,order(MutDescription,Late50Ratio_LateMean,decreasing=TRUE)),]
row.names(lateSummary5)=1:nrow(lateSummary5)
lateSummary5$Late50Ratio_LateNormMean=lateSummary5$Late50Ratio_LateMean-lateSummary5[lateSummary5$Mutation=='WT','Late50Ratio_LateMean']
lateSummary5$Late200Ratio_LateNormMean=lateSummary5$Late200Ratio_LateMean-lateSummary5[lateSummary5$Mutation=='WT','Late200Ratio_LateMean']
write.csv(lateSummary5,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_late_all2.csv',row.names=FALSE)


########### STEP 3: Make merged file of summary data ##########
peakSummaryX=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peakNormData.csv',stringsAsFactors = FALSE)
peakSummarySQRTX=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peakNormSQRTData.csv',stringsAsFactors = FALSE)
peakSummarySQRTX=peakSummarySQRTX[,c(1,5,6)]
peakSummary90X=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peak90NormData.csv',stringsAsFactors = FALSE)
peakSummary90X=peakSummary90X[,-2]
peakSummary78SQRTX=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_NaIV_peak78NormSQRTData.csv',stringsAsFactors = FALSE)
peakSummary78SQRTX=peakSummary78SQRTX[,c(1,5,6)]
actSummaryX=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_act_all.csv',stringsAsFactors = FALSE)
actSummaryX=actSummaryX[,-2]
InacttimeSummaryX=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_inacttime_all.csv',stringsAsFactors=FALSE)
InacttimeSummaryX=InacttimeSummaryX[,-2]
inactSummaryX=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_inact_all.csv',stringsAsFactors = FALSE)
inactSummaryX=inactSummaryX[,-2]
RFISummaryX=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_RFI_all.csv',stringsAsFactors =FALSE)
RFISummaryX=RFISummaryX[,-2]
lateSummaryX=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_late_all2.csv',stringsAsFactors=FALSE)
lateSummaryX=lateSummaryX[,-2]

aa=merge(peakSummaryX, peakSummarySQRTX, all = TRUE)
aaa=merge(aa, peakSummary90X, all = TRUE)
aaaa = merge(aaa, peakSummary78SQRTX, all = TRUE)
ab=merge(aaaa,actSummaryX,all=TRUE)
abc=merge(ab,InacttimeSummaryX,all=TRUE)
abcd=merge(abc,inactSummaryX,all=TRUE)
abcde=merge(abcd,RFISummaryX,all=TRUE)
abcdef=merge(abcde,lateSummaryX,all=TRUE)
abcdef[,c('Mutation','MutDescription')]
allEP=abcdef

# Add in Z-scores based on benign distribution from VCCRI/VUMC collaboration
benign_mean <- 94.2
benign_sd <- 10.1 
allEP$z_score_120peak <- (allEP$PeakDensityNormSQRTMean-benign_mean) / benign_sd
allEP$z_se_120peak <- (allEP$PeakDensityNormSQRTSE-benign_sd) / benign_sd
allEP$ACMG_120peak <- convert_z_acmg(allEP$z_score_120peak)
summary(as.factor(allEP$MutDescription))
write.csv(allEP,'~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allEP_v3.csv',row.names=FALSE)

# Make summary of ACMG criteria for previously published variants
allEP <- read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-81_allEP_v3.csv')
variants <- read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SCN5A_mutant_ACMG_update.csv', header = TRUE)
variant_include <- as.vector(variants$Mutant)
allEP_subset <- allEP[allEP$Mutation %in% variant_include, ]
table <- allEP_subset[, c("Mutation", "PeakDensityNormSQRTMean", "PeakDensityNormSQRTSE",
                          "z_score_120peak", "z_se_120peak", "ACMG_120peak")]
write.csv(table, 'ACMG_scaled_criteria.csv')


