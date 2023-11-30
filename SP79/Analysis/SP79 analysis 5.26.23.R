# SP79 analysis on 5/26/2023 afternoon of experiment 
# LPneg line. Bunch of mutants, usual method (except lipofectamine).
# 4 plates A-D (SCN5A)


######### NaIV #########
# Read in and process NaIV (DONE)
a=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_raw.csv',skip=2,stringsAsFactors=FALSE)
a2=processInitialCsv(a,numSweeps=29,seqBy=5,firstVarStart=3,NA,capacitanceCol=7,capMin=5e-12,capMax=30e-12,sealCol=6,sealMin=0.5e9,sealMax=10e9,'SP79',averageCap=TRUE)
a2[1:384,1]='SP79A'
a2[385:768,1]='SP79B'
a2[769:1152,1]='SP79C'
a2[1153:1536,1]='SP79D'

# add in peak current and current density (DONE)
a2b=addInPeak(a2,6,6,5+29,'Peak')
a2c=addInPeakDensity(a2b,4,6,7)
b=a2c[a2c$GoodSeal,] #530/(384*4)=35%, pretty bad
plot(a2c$Seal/10^9,10^12*a2c$Capacitance,pch='.',xlab='Seal (GOhm)',ylab='Capacitance (pF)',xlim=c(0,10),ylim=c(0,30))
abline(v=0.5)
abline(v=10)
abline(h=5)
abline(h=30)
write.csv(a2c,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_processed.csv',row.names=FALSE)

# Make IV plot of all goodSeal traces (DONE)
sweepList=1:29
IV_voltages=seq(-80,60,5)
a2c=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_processed.csv',stringsAsFactors=FALSE)
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_traces.pdf'
makeIVplots(a2c,pdfname,IV_voltages,8,7+length(sweepList))

# Manually looked through .pdf and scored as none (<50pA), good, low (current but <50) or jump, 
# Also annotating as Type=badSeal sweeps with <500 MOhm seal or capacitance <5 pA
# saved as SP79_NaIV_processed2.csv
# First columns are: 
# Count, Plate, Well, Mutation, MutDescription, GoodSeal, Type, Capacitance, Seal, Peak, PeakDensity, PeakInclude

# Get initial stats and calculate peak density norm WTIHOUT missing adjustment (START)
a3=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_processed2.csv',stringsAsFactors=FALSE)
wtMean=mean(a3[a3$PeakInclude & a3$Mutation=='WT','PeakDensity'],na.rm=TRUE)

a3a=addInPeakDensityNorm(a3,wtMean,11,12,'PeakDensityNorm')
a3b=addInPeak(a3a,13,26,26,'Peak20')
a3c=addInPeakDensity(a3b,8,13,14,'Peak20Density')
wtMean=mean(a3c[a3c$PeakInclude & a3c$Mutation=='WT','Peak20Density'],na.rm=TRUE)
a3d=addInPeakDensityNorm(a3c,wtMean,14,15,'Peak20DensityNorm')
write.csv(a3d,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_processed2b.csv',row.names=FALSE)
mutFrameStats=getStats(a3d)
# numCellsPeak is all the cells that are in PeakInclude 
# observed is numCellsPeak / GoodSeal 

mutFramePeakDensityNorm=summarizeByMutationGeneral(a3d,'PeakDensityNorm','PeakInclude')
mutFrameStats2=merge(mutFrameStats,mutFramePeakDensityNorm[,c(1,3)])
mutFrameStats3=orderMutFrame(mutFrameStats2,'PeakDensityNormMean','WT',decreasing=TRUE)
mutFrameStats3


# Edit SP79_missing.csv repeat now with updated variants included 

# to calculate missingCells in the missing.csv file, multiply the columns
# Difference and GoodSeal together, and divide by 100 (since it is a percentage)

# Summarize peak by mutation WITH missing cells ###  (DONE)
missing=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_missing.csv',stringsAsFactors=FALSE,header=TRUE)
missing2=missing[,c('Mutation','MissingCells')]
mutFramePeak2=summarizeByMutationGeneral(a3d,'Peak','PeakInclude',scaleBy=10^12,missingFrame=missing2)
mutFramePeakDensity2=summarizeByMutationGeneral(a3d,'PeakDensity','PeakInclude',missingFrame=missing2)
mutFramePeakDensityNorm2=summarizeByMutationGeneral(a3d,'PeakDensityNorm','PeakInclude',missingFrame=missing2)
mutFramePeak20=summarizeByMutationGeneral(a3d,'Peak20','PeakInclude',scaleBy=10^12,missingFrame=missing2)
mutFramePeak20Density=summarizeByMutationGeneral(a3d,'Peak20Density','PeakInclude',missingFrame=missing2)
mutFramePeak20DensityNorm=summarizeByMutationGeneral(a3d,'Peak20DensityNorm','PeakInclude',missingFrame=missing2)
mutFramePeakAll=cbind(mutFramePeak2,mutFramePeakDensity2[,4:5],mutFramePeakDensityNorm2[,4:5],mutFramePeak20Density[,4:5],mutFramePeak20DensityNorm[,4:5])
mutFrame1=merge(mutFrameStats3[,1:3],mutFramePeakAll)
mutFrame1b=orderMutFrame(mutFrame1,'PeakDensityNormMean',"WT",decreasing=TRUE)
mutFrame1b[,c('Mutation','MutDescription','GoodSeal','numCellsPeak','numCellsPeakWithMissing','PeakMean','PeakDensityNormMean','PeakDensityNormSE')]
write.csv(mutFrame1b,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary.csv',row.names = FALSE)



######### ACTIVATION #########
# Look at IV, for each voltage calculate peak current (NEW)
# then divide by expected current with reversal potential of 45.3 mV, then fit Boltzmann
# Resave without A-E.
a3d=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_processed2b.csv',stringsAsFactors=FALSE)
wtMean=mean(a3d[a3d$PeakInclude & a3d$Mutation=='WT','PeakDensity'],na.rm=TRUE)
outfile='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_VHalfAct.pdf'
reversalPotential=45.26055 #60 outside, 10 inside, temp 20C physiologyweb.com/calculators/nernst_potential_calculator.html
a3e=calcActivation(a3d,minPeakForAct=-100e-12,errorCutoff=0.1,IV_voltages=seq(-80,15,5),reversalPotential=45.26055,firstCol=17,peakCol='Peak',outfile=outfile)
write.csv(a3e,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_processed3.csv',row.names=FALSE)

# Summarize act
mutFrame1b=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary.csv',stringsAsFactors=FALSE)
a3e=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_processed3.csv',stringsAsFactors=FALSE)
mutFrameAct=summarizeByMutationGeneral(a3e,'VHalfAct','VHalfActInclude')
mutFrame2=merge(mutFrame1b,mutFrameAct)
mutFrame2b=orderMutFrame(mutFrame2,'VHalfActMean','WT',decreasing=TRUE)
mutFrame2b
write.csv(mutFrame2b,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary2.csv',row.names = FALSE)

###### INACTIVATION  #########
# Read in and process inact file and merge with act data
act=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaIV_processed3.csv',stringsAsFactors=FALSE)
inact=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaInact_raw.csv',skip=2,stringsAsFactors=FALSE)
inact2=processInitialCsv(inact,numSweeps=35,seqBy=5,firstVarStart=3,NA,capacitanceCol=7,capMin=5e-12,capMax=30e-12,sealCol=6,sealMin=0.5e9,sealMax=10e9,'SP79',averageCap=TRUE)
inact2$GoodSeal=inact2$GoodSeal & act$GoodSeal
names(inact2)[3:5]=c('GoodSeal_Inact','Capacitance_Inact','Seal_Inact')
numSweeps=35
inact2b=addInPeak(inact2,6,6,5+numSweeps,'Peak_Inact')
inact3=cbind(act[,1:19],inact2b[,3:ncol(inact2b)])
write.csv(inact3,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Inac_processed.csv',row.names=FALSE)

# Make inact IV plot of all goodSeal traces
sweepList=1:35
Inact_voltages=seq(-160,10,5)
inact3=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Inac_processed.csv',stringsAsFactors=FALSE)
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Inac_traces.pdf'
makeIVplots(inact3,pdfname,Inact_voltages,24,23+length(sweepList))

# region of interest
inact4=calcInactivation(inact3,-100e-12,0.1,Inact_voltages,24,'Peak_Inact','~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_VHalfInact.pdf')
b=inact4[inact4$VHalfInactInclude,]
plot(b$VHalfInactError,b$VHalfInact)
write.csv(inact4,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Inac_processed2.csv',row.names=FALSE)

# Summarize inact by mutation
inact4=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Inac_processed2.csv',stringsAsFactors=FALSE)
mutFrame2b=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary2.csv',stringsAsFactors = FALSE)
mutFrameInact=summarizeByMutationGeneral(inact4,'VHalfInact','VHalfInactInclude')
mutFrame3=merge(mutFrame2b,mutFrameInact)[, union(names(mutFrame2b), names(mutFrameInact))]
mutFrame3=orderMutFrame(mutFrame3,'PeakDensityNormMean','WT',decreasing=TRUE)
mutFrame3[,c('Mutation','MutDescription','PeakDensityNormMean','VHalfActMean','VHalfInactMean')]
write.csv(mutFrame3,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary3.csv',row.names=FALSE)

#Mutation MutDescription PeakDensityNormMean VHalfActMean VHalfInactMean
#21                         WT             WT               100.0        -57.2         -112.5
#8                       M923T            VUS               103.0        -58.3         -102.3
#18                     T1131I          parse               100.2        -48.0         -100.3
#12                      Q692K        Control                93.3        -61.2         -115.2
#3                      G1262S            VUS                91.4        -61.7         -115.0
#17                      S524Y        Control                81.2        -62.4         -119.0
#11                      P656L        Control                69.1        -65.5         -115.0
#2                      F1486Q            VUS                62.4        -48.0          -78.8
#1                      E1225K        Control                57.8        -42.3          -90.2
#9                       N406S            VUS                55.0        -49.8         -108.1
#14                     R1195C            VUS                48.1        -47.6         -105.6
#19                      T220I        Control                36.4        -63.2         -121.2
#20                       toss           toss                27.7        -59.1         -115.5
#7                      L1222R        Control                17.5        -39.0          -90.5
#13                        QQQ        Control                11.3        -49.1          -69.1
#16                     S1904L        Control                 9.3        -52.4         -117.8
#10 p.Asn1379_Lys1380insValPhe            VUS                 0.6        -29.1             NA
#6                       H886Q            VUS                 0.5           NA             NA
#4                      G1408R        Control                 0.4        -38.0          -88.7
#5                      G1743R        Control                 0.4           NA             NA
#15                      R878H        Control                 0.3           NA             NA





###### LATE CURRENT  V2 #########


# Read in files
inact4=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Inac_processed2.csv',stringsAsFactors=FALSE)
late=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_noLeakCorr.csv',skip=2,stringsAsFactors=FALSE)
late2_peak200=processInitialCsv(late,numSweeps=51,seqBy=12,firstVarStart=6,secondVarStart=8,capacitanceCol=14,capMin=5e-12,capMax=30e-12,sealCol=13,sealMin=0.5e9,sealMax=10e9,'SP79',averageCap = FALSE)
late2_50ramp=processInitialCsv(late,numSweeps=51,seqBy=12,firstVarStart=7,secondVarStart=9,capacitanceCol=14,capMin=5e-12,capMax=30e-12,sealCol=13,sealMin=0.5e9,sealMax=10e9,'SP79',firstVarLetter='C',secondVarLetter='D',averageCap = FALSE)
late2_peak200_30=processInitialCsv(late,numSweeps=51,seqBy=12,firstVarStart=3,secondVarStart=5,capacitanceCol=14,capMin=5e-12,capMax=30e-12,sealCol=13,sealMin=0.5e9,sealMax=10e9,'SP79', firstVarLetter='E',secondVarLetter='F', averageCap = FALSE)
late2_50=processInitialCsv(late,numSweeps=51,seqBy=12,firstVarStart=4, secondVarStart=NA, capacitanceCol=14,capMin=5e-12,capMax=30e-12,sealCol=13,sealMin=0.5e9,sealMax=10e9,'SP79',firstVarLetter='G', secondVarLetter=NA, averageCap = FALSE)

late2=cbind(late2_peak200,late2_50ramp[,6:ncol(late2_50ramp)],late2_peak200_30[,6:ncol(late2_peak200_30)],late2_50[,6:ncol(late2_50)]) #A=Peak (-15), B=Late200 (-15), C=Late50 (-15), D=RampMin, E=Peak_30, F=Late200_30; G=Late50_30

# Edit GoodSeal based on pre and post-drug seals, capacitances
# 12-16=preRange, 22-26=postRange
sealPre=as.numeric(late[,145]) #sweep 12 seal
sealPost=as.numeric(late[,313]) #sweep 26 seal
plot(log10(sealPre),log10(sealPost))
capPre=as.numeric(late[,146]) #sweep 12 cap
capPost=as.numeric(late[,314]) #sweep 26 cap
plot(log10(capPre),log10(capPost))
late2$GoodSeal=late2$GoodSeal & !is.nan(sealPre) & !is.nan(sealPost) & sealPre>0.5e9 & sealPost>0.5e9 & sealPre<10e9 & sealPost<10e9 & sealPost/sealPre<1.5 & sealPost/sealPre>(1/1.5) & !is.nan(capPre) & !is.nan(capPost) & capPre>5e-12 & capPre<30e-12 & capPost/capPre<1.5 & capPost/capPre>(1/1.5)
names(late2)[3:5]=c('GoodSeal_Late','Capacitance_Late','Seal_Late')
late2$Capacitance_Late=capPre
late2$Seal_Late=sealPre
#d=late2[late2$GoodSeal_Late,]
#plot(sealPre[late2$GoodSeal_Late],sealPost[late2$GoodSeal_Late])
#plot(sealPre[late2$GoodSeal_Late],(sealPost/sealPre)[late2$GoodSeal_Late])
#plot(capPre[late2$GoodSeal_Late],capPost[late2$GoodSeal_Late])
#plot(capPre[late2$GoodSeal_Late],(capPost/capPre)[late2$GoodSeal_Late])
late3=cbind(inact4[,1:26],late2[,3:ncol(late2)])
write.csv(late3,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Late_processed.csv',row.names=FALSE)

# Set range variables (DONE)
late3[1,] # Drug added at sweep 17-21, pre is 12-16, post is 22-26
numSweeps=51
drugStart=16
drugStop=22

# Plot capacitances and seals over time (DONE)
late=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_noLeakCorr.csv',skip=2,stringsAsFactors=FALSE)
late2_capSeal=processInitialCsv(late,numSweeps=51,seqBy=12,firstVarStart=14,secondVarStart=13,capacitanceCol=14,capMin=5e-12,capMax=30e-12,sealCol=13,sealMin=0.5e9,sealMax=10e9,'SP79',averageCap = FALSE)
late2_capSeal=late2_capSeal[!is.na(late2_capSeal$GoodSeal),]
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_Capacitance_traces.pdf'
makeIVplots(late2_capSeal,pdfname,1:numSweeps,6,5+numSweeps,xlab='Sweep Number',keepCriteria = 'GoodSeal')
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_SealRes_traces.pdf'
makeIVplots(late2_capSeal,pdfname,1:numSweeps,57,107,scaleFactor=10^-6,xlab='Sweep Number',keepCriteria = 'GoodSeal')

#A=Peak (-15), B=Late200 (-15), C=Late50 (-15), D=RampMin, E=Peak_30, F=Late200_30; G=Late50_30

# Process Peak
varName='Peak'
names(late3)
varStart=30
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_Peak_traces.pdf'
late3_forplot=late3[!is.na(late3$GoodSeal) & !is.na(late3$Seal_Late),]
makeIVplots(late3_forplot,pdfname,1:numSweeps,varStart,varStart-1+numSweeps,keepCriteria = 'GoodSeal_Late',xlab='Sweep Number')
late3b=addInMean(late3,ncol(late3)+1,varStart+drugStart-6,varStart+drugStart-2,paste0(varName,'Pre_Late'))
late3b=addInMean(late3b,ncol(late3b)+1,varStart+drugStop,varStart+drugStop+4,paste0(varName,'Post_Late'))
diff=late3b[,ncol(late3b)-1]-late3b[,ncol(late3b)]
plot(diff,late3b$PeakPre_Late)
late3b=cbind(late3b,diff)
names(late3b)[ncol(late3b)]=paste0(varName,'_Late')
b=late3b[!is.na(late3b$GoodSeal_Late) & late3b$GoodSeal_Late,]
plot(10^12*b$PeakPre_Late,10^12*b$PeakPost_Late,xlab=paste(varName,'Pre (nA)'),ylab=paste(varName,'Post (nA)'),ylim=c(min(10^12*b$PeakPre_Late),500))
abline(0,1)
plot(10^12*b$PeakPre_Late,10^12*b$Peak)

# Process Late200 - now I get an error - means are NA 

varName='Late200'
names(late3b)
varStart=81
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_Late200_traces.pdf'
makeIVplots(late3b,pdfname,1:numSweeps,varStart,varStart-1+numSweeps,keepCriteria = 'GoodSeal_Late',xlab='Sweep Number')
late3c=addInMean(late3b,ncol(late3b)+1,varStart+drugStart-6,varStart+drugStart-2,paste0(varName,'Pre_Late'))
late3c=addInMean(late3c,ncol(late3c)+1,varStart+drugStop,varStart+drugStop+4,paste0(varName,'Post_Late'))
diff=late3c[,ncol(late3c)-1]-late3c[,ncol(late3c)]
late3c=cbind(late3c,diff)
names(late3c)[ncol(late3c)]=paste0(varName,'_Late')
late3c$Late200Ratio_Late=late3c$Late200_Late/late3c$Peak_Late
b=late3c[late3c$GoodSeal_Late,]
plot(10^12*b$Late200Pre_Late,10^12*b$Late200Post_Late,xlab=paste(varName,'Pre (nA)'),ylab=paste(varName,'Post (nA)'))
abline(0,1)
plot(10^12*b$Peak_Late,b$Late200Ratio_Late,xlab='Peak (nA)',ylab=paste(varName,'(nA)'),ylim=c(-.05,.05))
abline(h=0)
abline(v=-1000)
abline(v=-500)
d=b[10^12*b$Peak_Late<(-500) & 10^12*b$Peak_Late>(-1000),] #36
# need to look at QQQy and C341Y type variants that have teeny peak and large late
mean(d$Late200Ratio_Late) #.017
d=b[10^12*b$Peak_Late<(-1000),] 
mean(d$Late200Ratio_Late) 
sd(d$Late200Ratio_Late)

# Process Late50 
varName='Late50'
names(late3c)
varStart=132
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_Late50_traces.pdf'
makeIVplots(late3,pdfname,1:numSweeps,varStart,varStart-1+numSweeps,keepCriteria = 'GoodSeal_Late',xlab='Sweep Number')
late3d=addInMean(late3c,ncol(late3c)+1,varStart+drugStart-6,varStart+drugStart-2,paste0(varName,'Pre_Late'))
late3d=addInMean(late3d,ncol(late3d)+1,varStart+drugStop,varStart+drugStop+4,paste0(varName,'Post_Late'))
diff=late3d[,ncol(late3d)-1]-late3d[,ncol(late3d)]
late3d=cbind(late3d,diff)
names(late3d)[ncol(late3d)]=paste0(varName,'_Late')
late3d$Late50Ratio_Late=late3d$Late50_Late/late3d$Peak_Late
b=late3d[late3d$GoodSeal_Late,]
plot(10^12*b$Late50Pre_Late,10^12*b$Late50Post_Late,xlab=paste(varName,'Pre (nA)'),ylab=paste(varName,'Post (nA)'))
abline(0,1)
plot(10^12*b$Peak_Late,b$Late50Ratio_Late,xlab='Peak (nA)',ylab=paste(varName,'(nA)'),ylim=c(-.05,.05))
abline(h=0)
abline(v=-1000)
abline(v=-500)
d=b[10^12*b$Peak_Late<(-500) & 10^12*b$Peak_Late>(-1000),]
mean(d$Late50Ratio_Late,na.rm=TRUE)
d=b[10^12*b$Peak_Late<(-1000),] 
mean(d$Late50Ratio_Late,na.rm=TRUE) 
sd(d$Late50Ratio_Late) 

# Process RampMin
varName='RampMin'
names(late3d)
varStart=183
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_RampMin_traces.pdf'
makeIVplots(late3,pdfname,1:numSweeps,varStart,varStart-1+numSweeps,keepCriteria = 'GoodSeal_Late',xlab='Sweep Number')
late3e=addInMean(late3d,ncol(late3d)+1,varStart+drugStart-6,varStart+drugStart-2,paste0(varName,'Pre_Late'))
late3e=addInMean(late3e,ncol(late3e)+1,varStart+drugStop,varStart+drugStop+4,paste0(varName,'Post_Late'))
diff=late3e[,ncol(late3e)-1]-late3e[,ncol(late3e)]
late3e=cbind(late3e,diff)
names(late3e)[ncol(late3e)]=paste0(varName,'_Late')
late3e$RampMinRatio_Late=late3e$RampMin_Late/late3d$Peak_Late
b=late3e[late3e$GoodSeal_Late,]
plot(10^12*b$RampMinPre_Late,10^12*b$RampMinPost_Late,xlab=paste(varName,'Pre (nA)'),ylab=paste(varName,'Post (nA)'))
abline(0,1)
plot(10^12*b$Peak_Late,b$RampMinRatio_Late,xlab='Peak (nA)',ylab=paste(varName,'(nA)'),ylim=c(-.05,.05))
abline(h=0)
abline(v=-1000)
abline(v=-500)
d=b[10^12*b$Peak_Late<(-500) & 10^12*b$Peak_Late>(-1000),] 
mean(d$RampMinRatio_Late,na.rm=TRUE) 
d=b[10^12*b$Peak_Late<(-1000),] 
mean(d$RampMinRatio_Late,na.rm=TRUE) 
sd(d$RampMinRatio_Late) 


# Process Peak_30 - variable E, late3
varName='Peak_30'
names(late3)
varStart=234
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_Peak_30_traces.pdf'
late3_forplot=late3[!is.na(late3$GoodSeal) & !is.na(late3$Seal_Late),]
makeIVplots(late3_forplot,pdfname,1:numSweeps,varStart,varStart-1+numSweeps,keepCriteria = 'GoodSeal_Late',xlab='Sweep Number')
late3f=addInMean(late3e,ncol(late3e)+1,varStart+drugStart-6,varStart+drugStart-2,paste0(varName,'Pre_Late'))
late3f=addInMean(late3f,ncol(late3f)+1,varStart+drugStop,varStart+drugStop+4,paste0(varName,'Post_Late'))
diff=late3f[,ncol(late3f)-1]-late3f[,ncol(late3f)]
plot(diff,late3f$Peak_30Pre_Late)
late3f=cbind(late3f,diff)
names(late3f)[ncol(late3f)]=paste0(varName,'_Late')


# Process Late200_30 - variable F, late3g
varName='Late200_30'
names(late3b)
varStart=285
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_Late200_30_traces.pdf'
makeIVplots(late3,pdfname,1:numSweeps,varStart,varStart-1+numSweeps,keepCriteria = 'GoodSeal_Late',xlab='Sweep Number')
late3g=addInMean(late3f,ncol(late3f)+1,varStart+drugStart-6,varStart+drugStart-2,paste0(varName,'Pre_Late'))
late3g=addInMean(late3g,ncol(late3g)+1,varStart+drugStop,varStart+drugStop+4,paste0(varName,'Post_Late'))
diff=late3g[,ncol(late3g)-1]-late3g[,ncol(late3g)]
late3g=cbind(late3g,diff)
names(late3g)[ncol(late3g)]=paste0(varName,'_Late')
late3g$Late200_30Ratio_Late=late3g$Late200_30_Late/late3g$Peak_30_Late
b=late3g[late3g$GoodSeal_Late,]
plot(10^12*b$Late200_30Pre_Late,10^12*b$Late200_30Post_Late,xlab=paste(varName,'Pre (nA)'),ylab=paste(varName,'Post (nA)'))
abline(0,1)
plot(10^12*b$Peak_30_Late,b$Late200_30Ratio_Late,xlab='Peak (nA)',ylab=paste(varName,'(nA)'),ylim=c(-.05,.05))
abline(h=0)
abline(v=-1000)
abline(v=-500)


# Process Late50_30 - variable G, late3h
varName='Late50_30'
names(late3c)
varStart=336
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_CipaLate_Late50_traces.pdf'
makeIVplots(late3,pdfname,1:numSweeps,varStart,varStart-1+numSweeps,keepCriteria = 'GoodSeal_Late',xlab='Sweep Number')
late3h=addInMean(late3g,ncol(late3g)+1,varStart+drugStart-6,varStart+drugStart-2,paste0(varName,'Pre_Late'))
late3h=addInMean(late3h,ncol(late3h)+1,varStart+drugStop,varStart+drugStop+4,paste0(varName,'Post_Late'))
diff=late3h[,ncol(late3h)-1]-late3h[,ncol(late3h)]
late3h=cbind(late3h,diff)
names(late3h)[ncol(late3h)]=paste0(varName,'_Late')
late3h$Late50_30Ratio_Late=late3h$Late50_30_Late/late3h$Peak_30_Late
b=late3h[late3h$GoodSeal_Late,]
plot(10^12*b$Late50_30Pre_Late,10^12*b$Late50_30Post_Late,xlab=paste(varName,'Pre (nA)'),ylab=paste(varName,'Post (nA)'))
abline(0,1)
plot(10^12*b$Peak_30_Late,b$Late50_30Ratio_Late,xlab='Peak (nA)',ylab=paste(varName,'(nA)'),ylim=c(-.05,.05))
abline(h=0)
abline(v=-1000)
abline(v=-500)


# now make i

late3h$LateInclude=late3h$GoodSeal_Late & late3h$PeakPre_Late<(-1000e-12)
# recalculate the late include criteria in the combined analyses 
#.44, 18 cells
names(late3h)
late3i=cbind(late3h[,1:29],LateInclude=late3h[,413],late3h[,387:412],late3h[30:386])
write.csv(late3i,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Late_processed2.csv',row.names=FALSE)

b=late3i[late3i$LateInclude,]
plot(100*b$Late200Ratio_Late,100*b$RampMinRatio_Late,xlim=c(-3,3),ylim=c(-3,3))
model=lm(b$RampMinRatio_Late~b$Late200Ratio_Late)
summary(model) #0.03 (3%)
plot(100*b$Late200Ratio_Late,100*b$Late50Ratio_Late,xlim=c(-3,3),ylim=c(-3,3))
model=lm(b$Late50Ratio_Late~b$Late200Ratio_Late)
summary(model) #0.9
abline(0,1)
abline(0,2)

# Summarize late by mutation
late4=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Late_processed2.csv',stringsAsFactors=FALSE)
mutFrame4=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary3.csv',stringsAsFactors = FALSE)
mutFrameLate=summarizeByMutationGeneral(late4,'Late200Ratio_Late','LateInclude',roundUnits=5)
mutFrameLate2=summarizeByMutationGeneral(late4,'Late50Ratio_Late','LateInclude',roundUnits=5)
mutFrameLate3=summarizeByMutationGeneral(late4,'RampMinRatio_Late','LateInclude',roundUnits=5)
mutFrameLateAll=cbind(mutFrameLate,mutFrameLate2[,3:4],mutFrameLate3[,3:4])

mutFrame5=merge(mutFrame4,mutFrameLateAll)[, union(names(mutFrame4), names(mutFrameLateAll))]
mutFrame5=orderMutFrame(mutFrame5,'Late200Ratio_LateMean','WT',decreasing=TRUE)
mutFrame5[,c('Mutation','MutDescription','numCellsLate200Ratio_Late','Late200Ratio_LateMean','Late50Ratio_LateMean')]
mutFrame5[mutFrame5$numCellsLate200Ratio_Late>0,c('Mutation','MutDescription','numCellsLate200Ratio_Late','PeakDensityNormMean','Late200Ratio_LateMean','Late50Ratio_LateMean','RampMinRatio_LateMean')]
write.csv(mutFrame5,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary4.csv',row.names=FALSE)

# Complete on 5/26/2023 - update on 6/21/2023 w/ toss

###### RECOVERY FROM INACTIVATION ###### 
# Read in and process inact file and merge with act data (done for RecInac1)
late=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Late_processed2.csv',stringsAsFactors=FALSE)
rfi1=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaRecInac1.csv',skip=2,stringsAsFactors=FALSE)
rfi2=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaRecInac2.csv',skip=2,stringsAsFactors=FALSE)
rfi3=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_NaRecInac3.csv',skip=2,stringsAsFactors=FALSE)
RFI_times=c(1:10,seq(15,50,5),75,100,200)
#late=late[769:1920,]
#rfi1=rfi1[769:1920,]
#rfi2=rfi2[769:1920,]
#rfi3=rfi3[769:1920,]

rfi=processRFI(rfi1,rfi2,rfi3,late,-100e-12,RFI_times,3,4:10,4:10,4:10,1:45)
rfi[is.na(rfi$goodRFI),'goodRFI']=FALSE
write.csv(rfi,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_RecInac_processed.csv',row.names=FALSE)

# Make plots of raw and normalized RFI traces
sweepList=1:24
RFI_times=c(1:10,seq(15,50,5),75,100,200)
RFI_times_fake=c(0,0,0,RFI_times)
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_RFI_traces.pdf'
makeIVplots(rfi,pdfname,RFI_times_fake,47,46+length(sweepList),'goodRFI',xlab='Time(ms)')
pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_RFI_traces2.pdf'
makeIVplots(rfi,pdfname,RFI_times,71,70+length(sweepList)-3,'goodRFI',xlab='Time(ms)',scaleFactor=1)

# Calculate tau and error for all RFIs
rfi=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_RecInac_processed.csv',stringsAsFactors=FALSE)
rfi2=calcRFI(rfi,firstACol=47,Brange=71:91,pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_RFI_traces3.pdf',errorCutoff=.1,RFI_times=RFI_times)
write.csv(rfi2,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_RecInact_processed2.csv',row.names=FALSE)
plot(abs(10^12*rfi2$Peak),rfi2$rfiTau,ylim=c(0,50),xlim=c(0,5000))

plot(rfi2$rfi50,rfi2$rfiTau) # ln(2)*tau=halflife. tau=halflife/ln(2)
points(rfi2$rfi50,rfi2$rfi50/.693,col='red') # ok this is correct. I think half life is more intuitive #

# Summarize RFI by mutation
rfi2=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_RecInact_processed2.csv',stringsAsFactors=FALSE)
mutFrame5=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary4.csv',stringsAsFactors = FALSE)
mutFrameRFI=summarizeByMutationGeneral(rfi2,'rfi50','rfiInclude')
mutFrame6=merge(mutFrame5,mutFrameRFI,all.x=TRUE)
mutFrame6[mutFrame6$numCellsrfi50>1,c('Mutation','numCellsrfi50','rfi50Mean','rfi50SE')]
write.csv(mutFrame6,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary5.csv',row.names=FALSE)


#   Mutation numCellsrfi50 rfi50Mean rfi50SE
#3    G1262S             5       4.7     0.5
#9     N406S             3       5.0     0.5
#11    P656L             5      22.4     2.4
#12    Q692K             7      24.6     3.2
#13      QQQ             2       4.0     1.8
#16   S1787N             5      18.4     7.6
#17   S1904L            11      30.4     6.0
#18    S524Y             3      17.3     4.3
#21       WT             9      10.0     2.2



###### INACTIVATION TIME ######
rfi2=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_RecInact_processed2.csv',stringsAsFactors=FALSE)
inacttime=rfi2[,1:50]
summary(as.factor(inacttime$Plate))
inacttime_SP79A=inacttime[inacttime$Plate=='SP79A',]
inacttime_SP79B=inacttime[inacttime$Plate=='SP79B',]
inacttime_SP79C=inacttime[inacttime$Plate=='SP79C',]
inacttime_SP79D=inacttime[inacttime$Plate=='SP79D',]
inacttime2_SP79A=calcInactivationTime(inacttime_SP79A,folder='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/RawFiles/SP79A_rawTrace/',filePrefix='SP79A_NaIV_5mv_PP_09.30.05_',fileSuffix='_.csv',errorCutoff=0.1,pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79A_inactTime.pdf',makeFullPlot = TRUE,stimulusIncluded=TRUE,xlim=c(510,521), addTwo=FALSE)
inacttime2_SP79B=calcInactivationTime(inacttime_SP79B,folder='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/RawFiles/SP79B_rawTrace/',filePrefix='SP79B_NaIV_5mv_PP_10.08.05_',fileSuffix='_.csv',errorCutoff=0.1,pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79B_inactTime.pdf',makeFullPlot = TRUE,stimulusIncluded=TRUE,xlim=c(510,521), addTwo=FALSE)
inacttime2_SP79C=calcInactivationTime(inacttime_SP79C,folder='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/RawFiles/SP79C_rawTrace/',filePrefix='SP79C_NaIV_5mv_PP_10.42.55_',fileSuffix='_.csv',errorCutoff=0.1,pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79C_inactTime.pdf',makeFullPlot = TRUE,stimulusIncluded=TRUE,xlim=c(510,521), addTwo=FALSE)
inacttime2_SP79D=calcInactivationTime(inacttime_SP79D,folder='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/RawFiles/SP79D_rawTrace/',filePrefix='SP79D_NaIV_5mv_PP_11.20.20_',fileSuffix='_.csv',errorCutoff=0.1,pdfname='~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79D_inactTime.pdf',makeFullPlot = TRUE,stimulusIncluded=TRUE,xlim=c(510,521), addTwo=FALSE)

inacttime2=rbind(inacttime2_SP79A,inacttime2_SP79B,inacttime2_SP79C,inacttime2_SP79D)
sum(inacttime2$inacttimeInclude) #409 (27% of total data)
plot(inacttime2$inacttimeError,inacttime2$inacttimeTau,xlab='Inacttime error',ylab='Inacttime Tau')
points(inacttime2$inacttimeError[inacttime2$inacttimeInclude],inacttime2$inacttimeTau[inacttime2$inacttimeInclude],col='blue')
abline(v=0.1)
plot(inacttime2$inacttimeError[inacttime2$inacttimeInclude],inacttime2$inacttimeTau[inacttime2$inacttimeInclude],xlab='Inacttime error',ylab='Inacttime Tau')
# zoom in from 0-10% inacttime error
write.csv(inacttime2,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_InactTime_processed.csv',row.names=FALSE)

# Summarize Inacttime by mutation
inacttime2=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_InactTime_processed.csv',stringsAsFactors=FALSE)
mutFrame6=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary5.csv',stringsAsFactors = FALSE)
mutFrameInactTime=summarizeByMutationGeneral(inacttime2,'inacttimeTau','inacttimeInclude')
mutFrame7=merge(mutFrame6,mutFrameInactTime,all.x=TRUE)
mutFrame7[,c('Mutation','numCellsinacttimeTau','inacttimeTauMean','inacttimeTauSE')]
write.csv(mutFrame7,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary6.csv',row.names=FALSE)

pdf('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_inacttime_plots.pdf')
mutList=mutFrame7$Mutation
for(mutation in mutList){
  inacttime2_wt=inacttime2[inacttime2$Mutation=='WT' & inacttime2$inacttimeInclude,]
  inacttime2_mut=inacttime2[inacttime2$Mutation==mutation & inacttime2$inacttimeInclude,]
  if(nrow(inacttime2_mut)>0){
    plot(inacttime2_mut$inacttimeError,inacttime2_mut$inacttimeTau,col='red',xlim=c(0,0.1),ylim=c(min(inacttime2_wt$inacttimeTau),max(inacttime2_mut$inacttimeTau)),ylab='Inactivation time Tau',xlab='Inactivation time error',main=paste(mutation,mutFrame7[mutFrame7$Mutation==mutation,'inacttimeTauMean']))
    points(inacttime2_wt$inacttimeError,inacttime2_wt$inacttimeTau)
  }
}
dev.off()


###### Peak current from inact protocol #####
# Analyze from inactivation protocol the peak current from the -90 sweep (sweep 15)
# 3.2.23 example to use for calculating peak current from the -90 sweep. Missing cells2 now semi-automated
# Read in and process inact time file as a starting point
peak90=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_InactTime_processed.csv',stringsAsFactors =FALSE)
inact=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Inac_processed2.csv',stringsAsFactors =FALSE)
names(inact) # look for sweep 15
peak90$Peak90Include=peak90$PeakInclude & peak90$GoodSeal_Inact
# SKIP Error cutoff < 0.1
# SKIP Peak cutoff < -100e-12
# Allow good or jump (old peak include rules)
# Get initial stats and calculate peak90 density norm WTIHOUT missing adjustment (START)
peak90$Peak90=inact[,41] #sweep 15
peak90b=addInPeakDensity(peak90,8,55,56,'Peak90Density')
wtMean=mean(peak90b[peak90b$Peak90Include & peak90b$Mutation=='WT','Peak90Density'],na.rm=TRUE)
peak90c=addInPeakDensityNorm(peak90b,wtMean,56,57,'Peak90DensityNorm')

peak90c$Peak78Include=peak90$PeakInclude & peak90$GoodSeal_Inact
peak90c$Peak78=(inact[,43]*3+inact[,44]*2)/5 # get -78 mV data for Sydney 
peak90d=addInPeakDensity(peak90c,8,59,60,'Peak78Density')
wtMean=mean(peak90d[peak90d$Peak78Include & peak90d$Mutation=='WT','Peak78Density'],na.rm=TRUE)
peak90e=addInPeakDensityNorm(peak90d,wtMean,60,61,'Peak78DensityNorm')

write.csv(peak90e,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_Peak90_processed2.csv',row.names=FALSE)
mutFrameStats=getStats(peak90c,varName="Peak90",goodSealCol = "GoodSeal_Inact")
# numCellsPeak is all the cells that are in PeakInclude 
# observed is numCellsPeak / GoodSeal 

mutFramePeak90DensityNorm=summarizeByMutationGeneral(peak90c,'Peak90DensityNorm','Peak90Include')
mutFrameStats2=merge(mutFrameStats,mutFramePeak90DensityNorm[,c(1,3)])
mutFrameStats3=orderMutFrame(mutFrameStats2,'Peak90DensityNormMean','WT',decreasing=TRUE)
names(mutFrameStats3)=c('Mutation','MutDescription','GoodSeal','numCellsPeak','Observed','PeakDensityNormMean')
missingFrame=mutFrameStats3[,c(1,3,4,5,2,6)]
oldMissingFrame=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_missing.csv',stringsAsFactors=FALSE)
oldMissingFrame=oldMissingFrame[,c(1,5)]
missingFrame2=merge(missingFrame,oldMissingFrame)
missingFrame2=missingFrame2[order(missingFrame2$PeakDensityNormMean,decreasing=TRUE),]
missingFrame2$Difference=missingFrame2$Expected-missingFrame2$Observed
missingFrame2$MissingCells=(missingFrame2$Difference/100)*missingFrame2$GoodSeal
missingFrame2$MissingCells[missingFrame2$Difference < 10]=0
missingFrame2$PeakAdjustment=NA
missingFrame2$Include=TRUE
missingFrame3=missingFrame2[,c(1,2,3,4,7,8,9,10,11,5,6)]
write.csv(missingFrame3,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_missing2.csv',row.names=FALSE)
# Edit SP79_missing2.csv manually if necessary

# Summarize peak by mutation WITH missing cells ###  (DONE)
missing=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_missing2.csv',stringsAsFactors=FALSE,header=TRUE)
missing2=missing[,c('Mutation','MissingCells')]
mutFramePeak90=summarizeByMutationGeneral(peak90c,'Peak90','Peak90Include',scaleBy=10^12,missingFrame=missing2)
mutFramePeak90Density=summarizeByMutationGeneral(peak90c,'Peak90Density','Peak90Include',missingFrame=missing2)
mutFramePeak90DensityNorm=summarizeByMutationGeneral(peak90c,'Peak90DensityNorm','Peak90Include',missingFrame=missing2)
mutFramePeakAll=cbind(mutFramePeak90,mutFramePeak90Density[,4:5],mutFramePeak90DensityNorm[,4:5])
mutFrame1=merge(mutFrameStats3[,1:3],mutFramePeakAll)
mutFrame1b=orderMutFrame(mutFrame1,'Peak90DensityNormMean',"WT",decreasing=TRUE)
names(mutFrame1b)[3]='GoodSeal90'
mutFrame1b[,c('Mutation','MutDescription','GoodSeal90','numCellsPeak90','numCellsPeak90WithMissing','Peak90Mean','Peak90DensityNormMean','Peak90DensityNormSE')]
mutFrame6=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary6.csv',stringsAsFactors=FALSE)
mutFrame7=merge(mutFrame6,mutFrame1b,all.x=TRUE)
write.csv(mutFrame7,'~/Dropbox/Andrew-Matthew/Syncropatch/SP79/Analysis/SP79_mutationSummary7.csv',row.names = FALSE)

plot(mutFrame7$PeakDensityNormMean,mutFrame7$Peak90DensityNormMean,xlim=c(0,350),ylim=c(0,350))
abline(0,1)

b=peak90c[peak90c$Peak90Include,]
plot(b$PeakDensity,b$Peak90Density,xlim=c(-800,0),ylim=c(-800,0))
sum(peak90c$PeakInclude) #448
sum(peak90c$Peak90Include) #132
model=lm(b$Peak90Density~b$PeakDensity)
abline(model)
summary(model) # estimate = 0.4, so 40% currents at 
abline(0,1,lty=2)



