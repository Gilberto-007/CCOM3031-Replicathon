# Import files from the Replicaton GitHub

  #PR2017Replicaton:
! wget https://raw.githubusercontent.com/areyesq89/PR2017replicaton/master/README.md
  
  # Main Analysis Template R Markdown:

! wget https://raw.githubusercontent.com/areyesq89/PR2017replicaton/master/analysis_template.Rmd

  # Datasets:

! wget https://raw.githubusercontent.com/areyesq89/PR2017replicaton/master/rawPharmacoData.csv
! wget https://raw.githubusercontent.com/areyesq89/PR2017replicaton/master/summarizedPharmacoData.csv

  # Tutorials:

! wget https://raw.githubusercontent.com/areyesq89/PR2017replicaton/master/R_basics.Rmd
! wget https://raw.githubusercontent.com/areyesq89/PR2017replicaton/master/exploration.Rmd
! wget https://raw.githubusercontent.com/areyesq89/PR2017replicaton/master/pharmaco_correlation.Rmd
! wget https://raw.githubusercontent.com/areyesq89/PR2017replicaton/master/targeted_therapies.Rmd
! wget https://raw.githubusercontent.com/areyesq89/PR2017replicaton/master/regression.Rmd

#Setup R enviorment
%load_ext rpy2.ipython


#The packages we need for the project
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("plyr")


#Setup the file that contains the data we need
rawFile <- "rawPharmacoData.csv"
pharmacoData <- read.csv(rawFile)


NROW(pharmacoData[["cellLine"]])

length(unique(factor(pharmacoData$drug)))

library(ggplot2)
ggplot( pharmacoData, aes(viability, group=study, colour=study)) + geom_histogram(fill = "green", colour="black")



range( pharmacoData$viability )
sum( pharmacoData$viability > 0 )
sum( pharmacoData$viability < 100 )


sumFile <- "summarizedPharmacoData.csv"
sumpharmacoData <- read.csv(sumFile)
head( sumpharmacoData )
str( sumpharmacoData )


ggplot( pharmacoData, aes(viability, group=study)) + geom_histogram(fill = "green", colour="black") + facet_wrap(~doseID)
    
ggplot( aes(x=auc_GDSC, y=auc_CCLE), data=subset(sumpharmacoData) ) + geom_point() + facet_wrap(facets = ~drug)  


library(dplyr)

drugsCorrs <- sumpharmacoData %>%
  group_by(drug) %>%
  summarise( Pearson_auc = cor(auc_GDSC, auc_CCLE, method="pearson"),
             Spearman_auc = cor(auc_GDSC, auc_CCLE, method="spearman") )
  
drugsCorrs

cor.test(sumpharmacoData[["auc_CCLE"]], sumpharmacoData[["auc_GDSC"]]) 

cor.test(sumpharmacoData[["ic50_CCLE"]], sumpharmacoData[["ic50_GDSC"]])

cor.test(sumpharmacoData[["auc_CCLE"]], sumpharmacoData[["auc_GDSC"]], method = "pearson") 

cor.test(sumpharmacoData[["ic50_CCLE"]], sumpharmacoData[["ic50_GDSC"]], method = "pearson")

cor.test(sumpharmacoData[["auc_CCLE"]], sumpharmacoData[["auc_GDSC"]], method = "spearman")

cor.test(sumpharmacoData[["ic50_CCLE"]], sumpharmacoData[["ic50_GDSC"]], method = "spearman")


AUC_study1 <- rbeta(200, 1, 5)
AUC_study2 <- rbeta(200, 1, 5)
resistant <- data.frame(AUC_study1, AUC_study2)

resistant


AUC_study1 <- c(rbeta(100, 1, 5), rbeta(100, 4, 2))
AUC_study2 <- c(rbeta(100, 1, 5), rbeta(100, 4, 2))
resistant <- data.frame(AUC_study1, AUC_study2, 
                        CellLine=c(rep("Resistant", 100), rep("Sensitive", 100)))
ggplot(resistant, aes( y=AUC_study2, x=AUC_study1, colour=CellLine) ) +
    geom_point() + ggtitle("Simulated AUC with half sensitive and half resistant cell lines") +
    xlim(0,1) + ylim(0,1)


cellLinesSummary <- read.csv("summarizedPharmacoData.csv", header=TRUE)

drugAvg <- cellLinesSummary %>% 
              group_by(cellLine) %>%
              summarise(mean_ic50_CCLE = mean(-log10(ic50_CCLE/10^6)), 
                        mean_ic50_GDSC = mean(-log10(ic50_GDSC/10^6)),
                        mean_auc_CCLE = mean(auc_CCLE),
                        mean_auc_GDSC = mean(auc_GDSC)) 

ggplot(drugAvg, aes(x=mean_ic50_GDSC, y=mean_ic50_CCLE)) +
    geom_point(alpha=0.6) +
    ggtitle("Average IC50 value by cell line (averaged over drugs)")
    
    
cellLinesSummary <- cellLinesSummary %>% 
              mutate(cutoff = ifelse(drug=="paclitaxel", 0.4, 0.1)) %>%
              mutate(sensitivity_GDSC = factor(ifelse( auc_GDSC < cutoff, "Resistant", "Sensitive")), 
                     sensitivity_CCLE = factor(ifelse( auc_CCLE < cutoff, "Resistant", "Sensitive"))) 

table("GDSC"=cellLinesSummary$sensitivity_GDSC, "CCLE"=cellLinesSummary$sensitivity_CCLE)
    

mcc <- function (study1, study2)
{
  BS <- sum(study1 == "Sensitive" & study2 == "Sensitive") 
  BR <- sum(study1 == "Resistant" & study2 == "Resistant") 
  SR <- sum(study1 == "Sensitive" & study2 == "Resistant") 
  RS <- sum(study1 == "Resistant" & study2 == "Sensitive") 
  
  if (BS+SR == 0 | BS+RS == 0 | BR+SR == 0 |  BR+RS ==0){
    mcc <- ((BS*BR)-(SR*RS)) 
  }else{
    mcc <- ((BS*BR)-(SR*RS)) / sqrt(exp((log(BS+SR)+log(BS+RS)+log(BR+SR)+log(BR+RS))))
  }
  return(mcc)
}

No_Effect <- c("Sorafenib", "Erlotinib", "PHA-665752")

Narrow_Effect <- c("Nilotinib","Lapatinib", "Nutlin-3","PLX44720", "Crizotinib","PD-0332991", "AZD0530", "TAE684")

cellLinesSummary <- cellLinesSummary %>% 
  group_by(drug) %>% 
  mutate(drugClass=ifelse(drug==Narrow_Effect, "Narrow Effect", ifelse(drug==No_Effect, "No Effect","Broad Effect")))

cellLinesSummaryEffects <- cellLinesSummary %>% 
  group_by(drugClass) %>% 
  summarise(matthews_corr=mcc(sensitivity_GDSC, sensitivity_CCLE)) 
  cellLinesSummaryEffects
  
  
#Using the total raw data we eliminate the values below 50% viability and 
#substract them from the total raw data
rawData = read.csv("rawPharmacoData.csv")


Raw_amount_below50 <- length(unique(rawData$cellLine[rawData$viability < 0.5]))

Raw_amount <- length(unique(rawData$cellLine))

Potential_amount <- Raw_amount - Raw_amount_below50
Potential_amount


library(ggplot2)
library(dplyr)

plotResponse <- function(drugA, cellLineA, addPublishedIC50=TRUE ){
  pharSub <- filter( pharmacoData, drug==drugA, cellLine==cellLineA )
  sumSub <- filter( summarizedData, drug==drugA, cellLine==cellLineA )
  p <- ggplot( pharSub, aes( log10(concentration), viability, col=study)) +
      geom_point(size=2.1) + geom_line(lwd=1.1) + ylim(0, 150)
  if( addPublishedIC50 ){
      p <- p + geom_vline( sumSub, xintercept=log10( sumSub[,"ic50_CCLE"] ), col="#d95f02", linetype="longdash") +
          geom_vline( xintercept=log10( sumSub[,"ic50_GDSC"]), col="#1b9e77", linetype="longdash") +
          geom_hline( yintercept=50, col="#00000050", linetype="longdash")
  }
  p <- p + scale_colour_manual( values = c("CCLE" = "#d95f02", "GDSC" = "#1b9e77" ))
  xlims <- xlim( range(log10(c(pharSub$concentration, sumSub$ic50_CCLE, sumSub$ic50_GDSC ) ) ) )
  p + xlims
}

plotResponse(drugA = "17-AAG", cellLineA = "H4", TRUE)

library(ggplot2)
library(dplyr)

plotResponse <- function(drugA, cellLineA, addPublishedIC50=TRUE ){
  pharSub <- filter( pharmacoData, drug==drugA, cellLine==cellLineA )
  sumSub <- filter( summarizedData, drug==drugA, cellLine==cellLineA )
  p <- ggplot( pharSub, aes( log10(concentration), viability, col=study)) +
      geom_point(size=2.1) + geom_line(lwd=1.1) + ylim(0, 150)
  if( addPublishedIC50 ){
      p <- p + geom_vline( sumSub, xintercept=log10( sumSub[,"ic50_CCLE"] ), col="#d95f02", linetype="longdash") +
          geom_vline( xintercept=log10( sumSub[,"ic50_GDSC"]), col="#1b9e77", linetype="longdash") +
          geom_hline( yintercept=50, col="#00000050", linetype="longdash")
  }
  p <- p + scale_colour_manual( values = c("CCLE" = "#d95f02", "GDSC" = "#1b9e77" ))
  xlims <- xlim( range(log10(c(pharSub$concentration, sumSub$ic50_CCLE, sumSub$ic50_GDSC ) ) ) )
  p + xlims
}

plotResponse(drugA = "Nilotinib", cellLineA = "22RV1")


drugAvg = summarizedData %>% 
              group_by(cellLine) %>%
              summarise(mean_ic50_CCLE = mean(-log10(ic50_CCLE/10^6)), 
                        mean_ic50_GDSC = mean(-log10(ic50_GDSC/10^6)),
                        mean_auc_CCLE = mean(auc_CCLE),
                        mean_auc_GDSC = mean(auc_GDSC)) 

ggplot(drugAvg, aes(x = mean_ic50_GDSC, y = mean_ic50_CCLE)) +
    geom_point(alpha = 0.6) +
    ggtitle("Average IC50 value by cell line (averaged over drugs)")
    
    
linearModel = function(drugA, cellLineA, studyA){
    pharSub = filter(pharmacoData, drug == drugA, 
                      cellLine == cellLineA, study == studyA)
    pharSub$concentration = log10(pharSub$concentration)
    fit = lm(viability~ concentration, pharSub)
    coefficients(fit)["concentration"]
}

slopesCCLE = numeric(0)
slopesGDSC = numeric(0)
for (i in 1:length(summarizedData$drug)) {
  slopesCCLE = c(slopesCCLE, linearModel(summarizedData$drug[i],
                                         summarizedData$cellLine[i], "CCLE"))
  slopesGDSC = c(slopesGDSC, linearModel(summarizedData$drug[i], 
                                         summarizedData$cellLine[i], "GDSC"))
}

dataWithSlopes = cbind(summarizedData, slopesCCLE, slopesGDSC)


ggplot(aes(x = slopesCCLE, y = slopesGDSC), data = dataWithSlopes) +
    geom_point(cex=0.5) + 
    facet_wrap(facets=~drug) 


