library(readxl)
library(jtools)
library(ggplot2)

# Read excel with all the data of CogTE and ESB #
path_excelfile <- '.../Desktop/Complete Data CCA.xlsx'

sampleO <- read_xlsx(path_excelfile, 1)[,-1]
sampleR <- read_xlsx(path_excelfile, 2)[,-1]
sampleV_Nback <- read_xlsx(path_excelfile, 3)[,-1]
sampleV_Stroop <- read_xlsx(path_excelfile, 4)[,-1]

# Get the features to use: RMSSD, FC, PICs, Group
sampleO_RMSSD <- scale(sampleO[,1:15])
sampleO_betweenFC <- scale(sampleO[,16:36])
sampleO_withinFC <- scale(sampleO[,37:43])
sampleO_internalPIC <- sampleO[,44]
sampleO_externalPIC <- sampleO[,45]

sampleR_RMSSD <- scale(sampleR[,1:15])
sampleR_group <- sampleR[,22]
sampleR_internalPIC <- sampleR[,16]

sampleV_Nback_RMSSD <- scale(sampleV_Nback[,1:15])
sampleV_Stroop_RMSSD <- scale(sampleV_Stroop[,1:15])
sampleV_internalPIC<- sampleV_Nback[,16]
sampleV_MOCA <- sampleV_Nback[,21]

# Calculate CCA with SampleO
cca_Between <- cancor(sampleO_RMSSD, sampleO_betweenFC)
cca_Within <- cancor(sampleO_RMSSD, sampleO_withinFC)

# Obtain the canonical variances for RMSSD
varO_RMSSD_between  <- as.matrix(sampleO_RMSSD) %*% cca_Between$xcoef[, 1]
varO_RMSSD_within  <- as.matrix(sampleO_RMSSD) %*% cca_Within$xcoef[, 1]

varR_RMSSD_between <- as.matrix(sampleR_RMSSD) %*% cca_Between$xcoef[, 1]

varV_RMSSD_between_Stroop <- as.matrix(na.omit(sampleV_Stroop_RMSSD))%*% cca_Between$xcoef[, 1]
varV_RMSSD_between_Nback <- as.matrix(na.omit(sampleV_Nback_RMSSD))%*% cca_Between$xcoef[, 1]


# Since some subjects in ESB did only one of the tasks, we need to redo the matrices.
tempStroop <- c()
tempNback <- c()

i = 1
for(j in 1:nrow(sampleV_Stroop)){
  tempStroop[j] <- varV_RMSSD_between_Stroop[i]
  if (j == 13 || j == 23){ # Participants that didn't do Stroop
    tempStroop[j] <- NA
  }
  else{
    i = i + 1
  }
}

i = 1
for(j in 1:nrow(sampleV_Nback)){
  tempNback[j] <- varV_RMSSD_between_Nback[i]
  if (j == 8 || j == 25 || j == 29){# Participants that didn't do Nback
    tempNback[j] <- NA
  }
  else{
    i = i + 1
  }
}

varV_RMSSD_between_Stroop <- as.matrix(tempStroop)
varV_RMSSD_between_Nback <- as.matrix(tempNback)

# Dataframe of canonical variances and internal PIC
sampleO_Dataframe <- cbind(varO_RMSSD_between,varO_RMSSD_within,sampleO_internalPIC,sampleO_externalPIC)
sampleR_Dataframe <- cbind(varR_RMSSD_between,sampleR_internalPIC,sampleR_group)
sampleV_Dataframe <- na.omit(cbind(varV_RMSSD_between_Stroop,varV_RMSSD_between_Nback,sampleV_internalPIC,sampleV_MOCA))

# Linear regression
lm_SampleO_Between_InternalPIC <- lm(Internal_PIC~varO_RMSSD_between, data=as.data.frame(sampleO_Dataframe))
lm_SampleO_Within_InternalPIC <- lm(Internal_PIC~varO_RMSSD_within, data=as.data.frame(sampleO_Dataframe))
lm_SampleO_Between_ExternalPIC <- lm(External_PIC~varO_RMSSD_between, data=as.data.frame(sampleO_Dataframe))
lm_SampleO_Within_ExternalPIC <- lm(External_PIC~varO_RMSSD_within, data=as.data.frame(sampleO_Dataframe))

lm_SampleR_Between <- lm(Internal_PIC~varR_RMSSD_between+Group, data=as.data.frame(sampleR_Dataframe))

lm_SampleV_Between_Nback <- lm(Internal_PIC~varV_RMSSD_between_Nback+MOCA, data=as.data.frame(sampleV_Dataframe))
lm_SampleV_Between_Stroop <- lm(Internal_PIC~varV_RMSSD_between_Stroop+MOCA, data=as.data.frame(sampleV_Dataframe))

# Summary of Linear Regression
summary(lm_SampleO_Between_InternalPIC)
summary(lm_SampleO_Within_InternalPIC)
summary(lm_SampleO_Between_ExternalPIC)
summary(lm_SampleO_Within_ExternalPIC)

summary(lm_SampleR_Between)

summary(lm_SampleV_Between_Nback)
summary(lm_SampleV_Between_Stroop)

# Plots
sampleO_title=expression(Sample[O])
sampleR_title=expression(Sample[R])
sampleV_title=expression(Sample[V])

rmssdbet=expression(paste("RMSSD-betweenFC-", var[O]))
rmssdbetR=expression(paste("RMSSD-betweenFC-", var[R]))
rmssdbetV=expression(paste("RMSSD-betweenFC-", var[V]))

rmssdwit=expression(paste("RMSSD-withinFC-", var[O]))

loccog = expression(paste("LOC-", Cognition[O]))
loccogR = expression(paste("LOC-", Cognition[R]))
loccogV = expression(paste("LOC-", Cognition[V]))


lm_SampleO_Between_plot <- effect_plot(lm_SampleO_Between_InternalPIC, pred = varO_RMSSD_between, interval = TRUE, 
                          plot.points = TRUE)+ ggtitle(sampleO_title) +xlab(rmssdbet) +ylab(loccog)+xlim(-1,1) + 
                          ylim(1,6.5) + theme( panel.border = element_rect(colour = "black", fill=NA, size=1))

lm_SampleO_Within_plot <- effect_plot(lm_SampleO_Within_InternalPIC, pred = varO_RMSSD_within, interval = TRUE, 
                           plot.points = TRUE)+ ggtitle(sampleO_title) +xlab(rmssdwit) +ylab(loccog)+xlim(-1,1) + 
                            ylim(1,6.5) + theme( panel.border = element_rect(colour = "black", fill=NA, size=1))

lm_SampleR_Between_plot <- effect_plot(lm_SampleR_Between, pred = varR_RMSSD_between, interval = TRUE, 
                          plot.points = TRUE)+ ggtitle(sampleR_title) +xlab(rmssdbetR) +ylab(loccogR)+xlim(-1,1) + 
                          ylim(1,6.5) + theme( panel.border = element_rect(colour = "black", fill=NA, size=1))

lm_SampleV_Nback_plot <- effect_plot(lm_SampleV_Between_Nback, pred = varV_RMSSD_between_Nback, interval = TRUE, 
                          plot.points = TRUE)+ ggtitle(sampleV_title) +xlab(rmssdbetV) +ylab(loccogV)+xlim(-1,1) + 
                          ylim(1,6.5) + theme( panel.border = element_rect(colour = "black", fill=NA, size=1))

