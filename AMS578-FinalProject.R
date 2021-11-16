### Code for AMS578 Project ###
library(corrplot)
library(leaps)
library(knitr)
library(regclass)
library(Amelia)
library(MASS)

# Reading in the 3 data files.
print(getwd())
setwd("/Users/thomasconnolly/Downloads")
print(getwd())
IDE <- read.csv(("IDEgroup477375.csv"), header =  TRUE)[,-1]
IDG <- read.csv(("IDGgroup477375.csv"), header =  TRUE)[,-1]
IDY <- read.csv(("IDYgroup477375.csv"), header =  TRUE)[,-1]

# Merging the 3 data variables by ID number.
data1 <- merge(IDE, IDG, by="ID")
data2 <- merge(data1, IDY, by="ID")
data2

# Computing Summary Statistics for the merged data, before imputation.
summary_before <-summary(data2[,-1])
summary_before

apply(data2[,-1], MARGIN = 2,FUN =  sd, na.rm = TRUE)
apply(data2, MARGIN = 2, FUN = length)
cor <- cor(data2[,-1])
corr_NA <- cor(na.omit(data2[,-1]))
corr_NA2 <- cor(data2[,-1],use = "complete.obs")

# Fixing data by filling in missing data using the Amelia package.
set.seed(123)
a.out <- amelia(data2,noms=c('R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16','R17','R18','R19','R20','R21','R22','R23','R24','R25'),idvars = 'ID',m=1)
summary(a.out)
a.out$imputations$imp1
data <- a.out$imputations$imp1[,-1]

#Box-Cox transformation
(boxcox(lm(Y-min(Y)+1 ~ ., data=data))) #lambda=1 so most likely no transformation needed

# Correlations
abs(cor(data$E1, data$Y))  #.1098 ####
abs(cor(data$E2, data$Y))  #.2512 ####
abs(cor(data$E3, data$Y))  #.3085 ####
abs(cor(data$E4, data$Y))  #.4322 ####
abs(cor(data$E5, data$Y))  #.0120
abs(cor(data$E6, data$Y))  #.0258
abs(cor(data$R1, data$Y))  #.0287
abs(cor(data$R2, data$Y))  #.0049
abs(cor(data$R3, data$Y))  #.0285
abs(cor(data$R4, data$Y))  #.0013
abs(cor(data$R5, data$Y))  #.0074
abs(cor(data$R6, data$Y))  #.0014
abs(cor(data$R7, data$Y))  #.0090
abs(cor(data$R8, data$Y))  #.0307
abs(cor(data$R9, data$Y))  #.0318
abs(cor(data$R10, data$Y)) #.0021
abs(cor(data$R11, data$Y)) #.0057
abs(cor(data$R12, data$Y)) #.0397
abs(cor(data$R13, data$Y)) #.0080
abs(cor(data$R14, data$Y)) #.0296
abs(cor(data$R15, data$Y)) #.0351
abs(cor(data$R16, data$Y)) #.0028
abs(cor(data$R17, data$Y)) #.0142
abs(cor(data$R18, data$Y)) #.0057
abs(cor(data$R19, data$Y)) #.0171
abs(cor(data$R20, data$Y)) #.0498
abs(cor(data$R21, data$Y)) #.0111
abs(cor(data$R22, data$Y)) #.0008
abs(cor(data$R23, data$Y)) #.0230
abs(cor(data$R24, data$Y)) #.0182
abs(cor(data$R25, data$Y)) #.0390
abs(cor(data$Y, data$Y))   #1.0
corrplot(cor(data))

### Model Fitting
# Model with only intercept 
MI <- lm(Y~ 1, data=data)
summary(MI)

# Environment Only Model
M_E <- lm(Y~ E1+E2+E3+E4+E5+E6, data=data)
summary(M_E)
summary(M_E)$adj.r.squared #.383
BIC(M_E) #77516
AIC(M_E)#77473

# Full Linear Model
M1 <- lm(Y~(.), data = data)
summary(M1)
summary(M1)$coefficients[,4][summary(M1)$coefficients[,4]<0.01]
summary(M1)$adj.r.squared #.382
BIC(M1) #77678
AIC(M1) #77500
# The addition of R terms don't seem to help my model

# Full Model, Including All Possible Two-way Iteractions
M2 <- lm(Y~(.)^2, data = data)
summary(M2)
summary(M2)$coefficients[,4][summary(M2)$coefficients[,4]<0.001]
summary(M2)$adj.r.squared #.410
BIC(M2) #80476
AIC(M2) #77793

# Backwards Stepwise Regression On Full Linear Model
step(M1, direction='backward', scope=formula(MI), trace=0) #Outputs the following variables: E1,E2,E3,E4,R16,R19,R20

S1 <- lm(Y~(E1+E2+E3+E4+R16+R19+R20), data=data)
summary(S1)
summary(S1)$adj.r.squared #.386
summary(S1)$coefficients[,4][summary(S1)$coefficients[,4]<0.001]
BIC(S1) #77514
AIC(S1) #77466

# Only E1,E2,E3,E4 seem to be significant
# Model of all four way interactions between E1,E2,E3,E4
S2 <-lm(Y~(E1+E2+E3+E4)^4+poly(E1,4)+poly(E2,4)+poly(E3,4)+poly(E4,4), data=data)
summary(S2)
summary(S2)$adj.r.squared #.400
summary(S2)$coefficients[,4][summary(S2)$coefficients[,4]<0.001]

# Function below is used for creating best subset of variables using forward stepwise regression.
# I ran this on the model consisting of all interactions up to four-way between E1,E2,E3,E4 and
# choose the best model by BIC and adjusted R-squared to be Y = E1 + E2:E3:E4
var <- colnames(model.matrix(S2))
sets <- regsubsets(model.matrix(S2)[,-1],data$Y,nbest=1,nvmax=8,method='forward',intercept=TRUE)
best_subsets <- summary(sets)
subset_selection <-apply(best_subsets$which,1,function(x) paste0(var[x],collapse='+'))
subsets <-kable(data.frame(cbind(model=subset_selection ,adjR2=best_subsets$adjr2,BIC=best_subsets$bic)),caption='Model Summary')
subsets
summary(regsubsets(Y~(E1+E2+E3+E4)^4+poly(E1,4)+poly(E2,4)+poly(E3,4)+poly(E4,4), data=data, nvmax=5,method="forward"))

# Final Model and Results
final_model <- lm(Y~E1+E2:E3:E4, data = data)
summary(final_model)
BIC(final_model) #77461
AIC(final_model) #77440
summary(final_model)$adj.r.squared #.394
summary(final_model)$coefficients[,4][summary(final_model)$coefficients[,4]<0.001]

par(mfrow=c(2,2))
plot(final_model)
anova(final_model)
VIF(final_model)
