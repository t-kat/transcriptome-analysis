### Machine learning Project ###
#
# Analysis of variation of multi-tissue brain-gene regulation
# Discovery of underlying patterns
#
# 28.02.2020
###############################

# Read in data
data = read.table("expr4T.dat")
attach(data)
set.seed(1)

#---------------Project day 1-----------------------#
########################################## QUESTION 1

# finding the amount of NAs in our data to decide on the method to get rid of them
sum(is.na(data)) # 0 NAs in data

# filter genes based on minimum value and standard deviation/mean (coefficient of variation)
min.sm = 1
min.mean = 10

# mean and standard deviation calculations
nr.col = ncol(data) -1 # last column contains tissue label
mean = sapply(data[,1:nr.col], mean)
stand.dev = sapply(data[,1:nr.col], sd)
coef.var = stand.dev / mean

plot(coef.var, xlab = "Genes", ylab = "Coefficient of variation", main = "Coefficient of variation values for genes")
abline(1,0, col=2)
legend("topleft", legend = c("Minimum CV value"), col = c(2), lty = c(1))

# keep tissue label
mean = c(mean, min.mean+1)
coef.var = c(coef.var, min.sm+1)

dim(data)
length(which(coef.var > min.sm & mean > min.mean))
data2 = data[, which(coef.var > min.sm & mean > min.mean)]
dim(data2)

# check if gene "NSG00000271043.1_MTRNR2L2" is present
any(colnames(data) == "ENSG00000271043.1_MTRNR2L2")

########################################## QUESTION 2

# linear regression model using expression data from all tissues, ENSG00000271043.1_MTRNR2L2 as dependent variable

## MODEL 1
# split data into test and training set randomly
nr.rows <- nrow(data2)
trainindex <- sample(1:nr.rows, size=round(0.8*nr.rows), replace=FALSE)
train <- data2[trainindex, ]
test <- data2[-trainindex, ]

# fit model
lm_all <- lm(ENSG00000271043.1_MTRNR2L2 ~ . -tissue, data = train)  # MODEL 1 ALL TISSUES
summary(lm_all)

# how many variables are significant
sum(summary(lm_all)$coefficients[-1,4]<0.05)

par(mfrow = c(2, 2))  # change the number of plots displayed
par(mfrow = c(1, 1))
plot(lm_all)

# predictive value and true value comparison
lmall.pred <- predict(lm_all, test, type='response')

# plot Y values against predictions
plot(test$ENSG00000271043.1_MTRNR2L2, lmall.pred, xlab = "True values of gene of interest", ylab = "Predictive values",
     main = "MTRNR2L2 gene and all tissues")
abline(0,1, col=2) # 0 is intercept, 1 is slope
legend("topleft", legend = c("perfect linear line"), col = c(2), lty = c(1))

# root mean square error
value = (test$ENSG00000271043.1_MTRNR2L2 - lmall.pred)^2
average = mean(value)
sqrt(average)


## MODEL 2
# data from one specific tissue, ENSG00000271043.1_MTRNR2L2 as dependent variable

# choosing the tissue of interest
plot(data2$tissue, log(data2$ENSG00000271043.1_MTRNR2L2), 
     ylab = "log of gene expression levels", main = "MTRNR2L2 gene expression levels in all tissues", las = 2, cex.axis=0.70)

# selecting only samples from tissue of interest
caudate_tissue = "brain_spinalcord"
interest_data = data2[which(data2$tissue == caudate_tissue),]
interest_data = droplevels(interest_data)

# reducing nr of variables to get it to a level that's smaller than nr of observations
small = interest_data[, summary(lm_all)$coeff[-1,4] <0.05]

nr.rows2 <- nrow(small)
trainindex2 <- sample(1:nr.rows2, size=round(0.8*nr.rows2), replace=FALSE)
train2 <- small[trainindex2, ] # is small_set now train2?
test2 <- small[-trainindex2, ]

# new training set to reduce nr of variables
gene_interest <- interest_data["ENSG00000271043.1_MTRNR2L2"]
# which(names(train2) == "ENSG00000271043.1_MTRNR2L2")
small_set <- cbind(train2, ENSG00000271043.1_MTRNR2L2 = gene_interest[trainindex2,])#as.matrix(interest_data[,gene_interest])) # add missing gene row to the small

# new test set 
#gene_interest <- which(names(test2) == "ENSG00000271043.1_MTRNR2L2")
small_test <- cbind(test2, ENSG00000271043.1_MTRNR2L2 = gene_interest[-trainindex2,]) # add missing gene row to the small
#small_set$ENSG00000271043.1_MTRNR2L2 = interest_data$ENSG00000271043.1_MTRNR2L2 
#which(names(small_set) == "tissue") # to get index
#small = small[, -46] # to remove the tissue from data set

# fit model
lm_specific <- lm(ENSG00000271043.1_MTRNR2L2 ~ ., data = small_set) # MODEL 2 OF SELECTED TISSUE
summary(lm_specific)

# how many variables are significant
sum(summary(lm_specific)$coefficients[-1,4]<0.05)

par(mfrow = c(1, 1))
plot(lm_specific)

# predictive value and true value comparison
lmspecific.pred <- predict(lm_specific, small_test, type='response')

# plot Y values against predictions
plot(small_test$ENSG00000271043.1_MTRNR2L2, lmspecific.pred, xlab = "True values of gene of interest", ylab = "Predictive values",
     main = "MTRNR2L2 gene and spinal cord tissue")
abline(0,1, col=2) # 0 is intercept, 1 is slope
legend("topleft", legend = c("perfect linear line"), col = c(2), lty = c(1))

# root mean square error
value1 = (small_test$ENSG00000271043.1_MTRNR2L2 - lmspecific.pred)^2
average1 = mean(value1)
sqrt(average1)


## MODELS WITH OTHER GENE OF INTEREST
# picking other gene of interest
nr.col2 = ncol(data2) -1 # last column contains tissue label
mean2 = sapply(data2[,1:nr.col2], mean)
stand.dev2 = sapply(data2[,1:nr.col2], sd)
coef.var2 = stand.dev2 / mean2

# keep tissue label
mean2 = c(mean2, min.mean+1)
coef.var2 = c(coef.var2, min.sm+1)

five.max2 <-tail(sort(coef.var2),5) # ENSG00000161610.1_HCRT is selected


## MODEL 3
# all tissues, ENSG00000161610.1_HCRT as dependent variable

# fit model
lm_all1 <- lm(ENSG00000161610.1_HCRT ~ . -tissue, data = train)   # MODEL 3 ALL TISSUES
summary(lm_all1)

# how many variables are significant
sum(summary(lm_all1)$coefficients[-1,4]<0.05)

par(mfrow = c(1, 1))  # change the number of plots displayed
plot(lm_all1)

# predictive value and true value comparison
lmall.pred1 <- predict(lm_all1, test, type='response')

# plot Y values agains predictions
plot(test$ENSG00000161610.1_HCRT, lmall.pred1, xlab = "True values for gene of interest", ylab = "Predictive values",
     main = "HCRT gene and all tissues")
abline(0,1, col=2) # 0 is intercept, 1 is slope
legend("topright", legend = c("perfect linear line"), col = c(2), lty = c(1))

# root mean square error
value2 = (test$ENSG00000161610.1_HCRT - lmall.pred1)^2
average2 = mean(value2)
sqrt(average2)

## MODEL 4
# choosing the second tissue of interest
plot(data2$tissue, log(data2$ENSG00000161610.1_HCRT), ylab = "log of gene expression levels", las = 2,
     main = "HCRT gene expression levels in all tissues") # boxplot to pick tissue

# selecting only samples from tissue of interest
hypo_tissue = "brain_hypothalamus"
interest_data2 = data2[which(data2$tissue == hypo_tissue),]
interest_data2 = droplevels(interest_data2)

small1 = interest_data2[, summary(lm_all)$coeff[-1,4] <0.05]

nr.rows3 <- nrow(small1)
trainindex3 <- sample(1:nr.rows3, size=round(0.8*nr.rows3), replace=FALSE)
train3 <- small1[trainindex3, ]
test3 <- small1[-trainindex3, ]

# new training set
gene_interest1 <- interest_data2["ENSG00000161610.1_HCRT"]
small_set2 <- cbind(train3, ENSG00000161610.1_HCRT = gene_interest1[trainindex3,]) # add missing gene row to the small

# new test set 
small_test2 <- cbind(test3, ENSG00000161610.1_HCRT = gene_interest1[-trainindex3,])

lm_specific1 <- lm(ENSG00000161610.1_HCRT ~ ., data = small_set2) # MODEL 4 USING SELECTED TISSUE
summary(lm_specific1)

# how many variables are significant
sum(summary(lm_specific1)$coefficients[-1,4]<0.05)

par(mfrow = c(1, 1))  # change the number of plots displayed
plot(lm_specific1)

# predictive value and true value comparison
lspecific.pred1 <- predict(lm_specific1, small_test2, type='response')

# plot Y values agains predictions
plot(small_test2$ENSG00000161610.1_HCRT, lspecific.pred1, xlab = "True values for gene of interest", ylab = "Predictive values",
     main = "HCRT gene and hypothalamus tissue")
abline(0,1, col=2) # 0 is intercept, 1 is slope
legend("topright", legend = c("perfect linear line"), col = c(2), lty = c(1))

# root mean square error
value3 = (small_test2$ENSG00000161610.1_HCRT - lspecific.pred1)^2
average3 = mean(value3)
sqrt(average3)


# studentized residual plots without log-transforming the 1st model
student_res <- rstudent(lm_all)
lmall.pred.train <- predict(lm_all, train, type='response') # making predictions based on training data because of studentized residuals 

# plot of predictive values and studentized residuals
plot(lmall.pred.train, student_res, xlab = "Predictive values",
     ylab = "Studentized residuals", main = "Studentized residual plot without log transformation")
abline(3,0, col=2) # everything over 3 is an outlier
legend("bottomright", legend = c("Separation of outliers"), col = c(2), lty = c(1))

# log transform the data for 1st model
# split data into test and training set randomly
data4 <- data2[,-400] + 0.000001 # to get rid of tissue
log_data <- log(data4)

nr.rows4 <- nrow(log_data)
trainindex4 <- sample(1:nr.rows4, size=round(0.8*nr.rows4), replace=FALSE)
train4 <- log_data[trainindex4, ]
test4 <- log_data[-trainindex4, ]

# fit model
lm_all2 <- lm(ENSG00000271043.1_MTRNR2L2 ~ . , data = train4)  # MODEL 5 ALL TISSUES
summary(lm_all2)

# studentized residuals
student_res1 <- rstudent(lm_all2)
lmall.pred.train1 <- predict(lm_all2, train4, type='response') # making predictions based on training data because of studentized residuals 
# plot of predictive values and studentized residuals
plot(lmall.pred.train1, student_res1, xlab = "Predictive values",
     ylab = "Studentized residuals", main = "Studentized residual plot with log transformation")  # log transformation made the x-axis range smaller and separated expressed from not-expressed
abline(3,0, col=2) # everything over 3 is an outlier
legend("bottomright", legend = c("Separation of outliers"), col = c(2), lty = c(1))


########################################## QUESTION 3
## LOGISTIC REGRESSION

## brain_amygdala vs. brain_hippocampus
# get subset of 2 tissues
library(MASS)
tissue1 <- "brain_amygdala" #tissue 0
tissue2 <- "brain_hippocampus"  #tissue 1
mydata <- data2[which(data2$tissue==tissue1|data2$tissue==tissue2),]
mydata <- droplevels(mydata)

# new coef.var calculation for variable reduction purposes
nr.col5 = ncol(mydata) -1 # last column contains tissue label
mean5 = sapply(mydata[,1:nr.col5], mean)
stand.dev5 = sapply(mydata[,1:nr.col5], sd)
coef.var5 = stand.dev5 / mean5

small2 = mydata[, which(coef.var5 <0.6)]  # threshold selected based on number of variables
small2$tissue = mydata$tissue # to add tissue back to small2

# get test and training data
nr.rows5 <- nrow(small2)
trainindex5 <- sample(1:nr.rows5, size=round(0.8*nr.rows5), replace=FALSE)
train5 <- small2[trainindex5, ]
test5 <- small2[-trainindex5, ]

# fit the model
log_fit1 <- glm(tissue ~ . -tissue, data = train5, family = binomial)
summary(log_fit1)

# root mse
log_prob1.1 = predict(log_fit1, test5, type="response")
log_pred1[log_prob1 > .5]= 0
log_pred1[log_prob1 < .5]= 1


# calculate probabilities for specificity & sensitivity plot
log_prob1.1 = predict(log_fit1, test5, type="response")
log_pred1.1 = rep(tissue1, length(log_prob1.1))
log_pred1.1[log_prob1.1 > .25]= tissue2

matrix1.1 <- table(log_pred1.1, test5$tissue)

# calculate probabilities for specificity & sensitivity plot
log_prob1.2 = predict(log_fit1, test5, type="response")
log_pred1.2 = rep(tissue1, length(log_prob1.2))

log_pred1.2[log_prob1.2 > .75]= tissue2

matrix1.2 <- table(log_pred1.2, test5$tissue)

# calculate probabilities
log_prob1 = predict(log_fit1, test5, type="response")
log_pred1 = rep(tissue1, length(log_prob1))
log_pred1[log_prob1 > .5]= tissue2

# confusion matrix
matrix1 <-table(log_pred1, test5$tissue) 

# overall fraction of correct predictions
mean(log_pred1 == test5$tissue)

# test error rate
mean(log_pred1.1 != test5$tissue)
mean(log_pred1 != test5$tissue)

mean(log_pred1.2 != test5$tissue)


## brain_hippocampus vs. brain_nucleusaccumbens
# get subset of 2 tissues
tissue1 <- "brain_hippocampus"
tissue2 <- "brain_nucleusaccumbens"
mydata <- data2[which(data2$tissue==tissue1|data2$tissue==tissue2),]
mydata <- droplevels(mydata)

# new coef.var calculation for variable reduction purposes
nr.col3 = ncol(mydata) -1 # last column contains tissue label
mean3 = sapply(mydata[,1:nr.col3], mean)
stand.dev3 = sapply(mydata[,1:nr.col3], sd)
coef.var3 = stand.dev3 / mean3

small2 = mydata[, which(coef.var3 <0.6)]  # threshold selected based on number of variables
small2$tissue = mydata$tissue # to add tissue back to small2

# get test and training data
nr.rows5 <- nrow(small2)
trainindex5 <- sample(1:nr.rows5, size=round(0.8*nr.rows5), replace=FALSE)
train5 <- small2[trainindex5, ]
test5 <- small2[-trainindex5, ]
# fit the model
log_fit2 <- glm(tissue ~ . -tissue, data = train5, family = binomial)
summary(log_fit2)

# calculate probabilities for sensitivity and specificity plot 
log_prob2.1 = predict(log_fit2, test5, type="response")
log_pred2.1 = rep(tissue1, length(log_prob2.1))
log_pred2.1[log_prob2.1 > .25]= tissue2

matrix2.1 = table(log_pred2.1, test5$tissue) 

# calculate probabilities for sensitivity and specificity plot 
log_prob2 = predict(log_fit2, test5, type="response")
log_pred2 = rep(tissue1, length(log_prob2))
log_pred2[log_prob2 > .5]= tissue2

# calculate probabilities
log_prob2.2 = predict(log_fit2, test5, type="response")
log_pred2.2 = rep(tissue1, length(log_prob2.2))
log_pred2.2[log_prob2.2 > .25]= tissue2

matrix2.2 = table(log_pred2.2, test5$tissue)

# confusion matrix
matrix2 = table(log_pred2, test5$tissue) 

# overall fraction of correct predictions
mean(log_pred2 == test5$tissue)

# test error rate
mean(log_pred2.1 != test5$tissue)
mean(log_pred2 != test5$tissue)
mean(log_pred2.2 != test5$tissue)


## brain_spinalcord vs. brain_substantianigra
# get subset of 2 tissues
tissue1 <- "brain_spinalcord"
tissue2 <- "brain_substantianigra"
mydata <- data2[which(data2$tissue==tissue1|data2$tissue==tissue2),]
mydata <- droplevels(mydata)

# new coef.var calculation for variable reduction purposes
nr.col3 = ncol(mydata) -1 # last column contains tissue label
mean3 = sapply(mydata[,1:nr.col3], mean)
stand.dev3 = sapply(mydata[,1:nr.col3], sd)
coef.var3 = stand.dev3 / mean3

small2 = mydata[, which(coef.var3 <0.6)]  # threshold selected based on number of variables
small2$tissue = mydata$tissue # to add tissue back to small2

# get test and training data
nr.rows5 <- nrow(small2)
trainindex5 <- sample(1:nr.rows5, size=round(0.8*nr.rows5), replace=FALSE)
train5 <- small2[trainindex5, ]
test5 <- small2[-trainindex5, ]
# fit the model
log_fit3 <- glm(tissue ~ . -tissue, data = train5, family = binomial)
summary(log_fit3)

# calculate probabilities for sensitivity and specificity plot 
log_prob3.1 = predict(log_fit3, test5, type="response")
log_pred3.1 = rep(tissue1, length(log_prob3.1))
log_pred3.1[log_prob3.1 > .25]= tissue2

matrix3.1 = table(log_pred3.1, test5$tissue) 

# calculate probabilities for sensitivity and specificity plot 
log_prob3.2 = predict(log_fit3, test5, type="response")
log_pred3.2 = rep(tissue1, length(log_prob3.2))
log_pred3.2[log_prob3.2 > .75]= tissue2

matrix3.2 = table(log_pred3.2, test5$tissue) 

# calculate probabilities
log_prob3 = predict(log_fit3, test5, type="response")
log_pred3 = rep(tissue1, length(log_prob3))
log_pred3[log_prob3 > .5]= tissue2

# confusion matrix
matrix3 = table(log_pred3, test5$tissue) 

# overall fraction of correct predictions
mean(log_pred3 == test5$tissue)

# test error rate
mean(log_pred3.1 != test5$tissue)
mean(log_pred3 != test5$tissue) # cutoff value we used
mean(log_pred3.2 != test5$tissue)


## brain_hippocampus vs. brain_nucleusaccumbens
# get subset of 2 tissues
tissue1 <- "brain_hippocampus"
tissue2 <- "brain_nucleusaccumbens"
mydata <- data2[which(data2$tissue==tissue1|data2$tissue==tissue2),]
mydata <- droplevels(mydata)

# new coef.var calculation for variable reduction purposes
nr.col3 = ncol(mydata) -1 # last column contains tissue label
mean3 = sapply(mydata[,1:nr.col3], mean)
stand.dev3 = sapply(mydata[,1:nr.col3], sd)
coef.var3 = stand.dev3 / mean3

small2 = mydata[, which(coef.var3 <0.6)]  # threshold selected based on number of variables
small2$tissue = mydata$tissue # to add tissue back to small2

# get test and training data
nr.rows5 <- nrow(small2)
trainindex5 <- sample(1:nr.rows5, size=round(0.8*nr.rows5), replace=FALSE)
train5 <- small2[trainindex5, ]
test5 <- small2[-trainindex5, ]
# fit the model
log_fit4 <- glm(tissue ~ . -tissue, data = train5, family = binomial)
summary(log_fit4)

# calculate probabilities for sensitivity and specificity plot 
log_prob4.1 = predict(log_fit4, test5, type="response")
log_pred4.1 = rep(tissue1, length(log_prob4.1))
log_pred4.1[log_prob4.1 > .25]= tissue2

matrix4.1 = table(log_pred4.1, test5$tissue) 

# calculate probabilities for sensitivity and specificity plot 
log_prob4.2 = predict(log_fit4, test5, type="response")
log_pred4.2 = rep(tissue1, length(log_prob4.2))
log_pred4.2[log_prob4.2 > .75]= tissue2

matrix4.2 = table(log_pred4.1, test5$tissue) 

# calculate probabilities
log_prob4 = predict(log_fit4, test5, type="response")
log_pred4 = rep(tissue1, length(log_prob4))
log_pred4[log_prob4 > .5]= tissue2

# confusion matrix
matrix4 = table(log_pred4, test5$tissue) 

# overall fraction of correct predictions
mean(log_pred4 == test5$tissue)

# test error rate
mean(log_pred4.1 != test5$tissue)
mean(log_pred4 != test5$tissue) # cutoff value we used
mean(log_pred4.2 != test5$tissue)



## specificity & sensitivity plot
# create logistic regressions with different probabilities in predictions

# brain_amygdala vs. brain_hippocampus
install.packages("caret")
library(caret)

# cutoff value 0.25
spec1 = specificity(matrix1.1)
sens1 = sensitivity(matrix1.1)

# cutoff value 0.5
spec2 = specificity(matrix1)
sens2 = sensitivity(matrix1)

# cutoff value 0.75
spec3 = specificity(matrix1.2)
sens3 = sensitivity(matrix1.2)

plot(c(0.25, 0.5, 0.75), c(sens1, sens2, sens3), type="b", col='red', ylim = (0:1), 
     main = "Sensitivity and specificity for amygdala vs hippocampus", xlab = "cut-off values", ylab = "score")
lines(c(0.25, 0.5, 0.75), c(spec1, spec2, spec3), type="b",col = 'blue')
lines(c(0.25, 0.5, 0.75), c(0.3939394, 0.3939394, 0.3636364), type="b",col = 'black')
legend("topright", legend = c("sensitivity", "specificity", "error"), col = c("red", "blue", "black"), lty = c(1))

## brain_hippocampus vs. brain_nucleusaccumbens

# cutoff value 0.25
spec1.1 = specificity(matrix2.1)
sens1.1 = sensitivity(matrix2.1)

# cutoff value 0.5
spec2.1 = specificity(matrix2)
sens2.1 = sensitivity(matrix2)

# cutoff value 0.75
spec3.1 = specificity(matrix2.2)
sens3.1 = sensitivity(matrix2.2)

plot(c(0.25, 0.5, 0.75), c(sens1.1, sens2.1, sens3.1), type="b", col='red', ylim = (0:1), 
     main = "Sensitivity and specificity for hippocampus vs nucleusaccumbens", xlab = "cut-off values", ylab = "score")
lines(c(0.25, 0.5, 0.75), c(spec1.1, spec2.1, spec3.1), type="b",col = 'blue')
lines(c(0.25, 0.5, 0.75), c(0.1463415, 0.1707317, 0.1463415), type="b",col = 'black')
legend("bottomright", legend = c("sensitivity", "specificity", "error"), col = c("red", "blue", "black"), lty = c(1))

## brain_spinalcord vs. brain_substantianigra

# cutoff value 0.25
spec1.2 = specificity(matrix3.1)
sens1.2 = sensitivity(matrix3.1)

# cutoff value 0.5
spec2.2 = specificity(matrix3)
sens2.2 = sensitivity(matrix3)

# cutoff value 0.75
spec3.2 = specificity(matrix3.2)
sens3.2 = sensitivity(matrix3.2)

plot(c(0.25, 0.5, 0.75), c(sens1.2, sens2.2, sens3.2), type="b", col='red', ylim = (0:1), 
     main = "Sensitivity and specificity for spinalcord vs substantianigra", xlab = "cut-off values", ylab = "score")
lines(c(0.25, 0.5, 0.75), c(spec1.2, spec2.2, spec3.2), type="b",col = 'blue')
lines(c(0.25, 0.5, 0.75), c(0.3333333, 0.3333333, 0.3333333), type="b",col = 'black')
legend("bottomright", legend = c("sensitivity", "specificity", "error"), col = c("red", "blue", "black"), lty = c(1))

## brain_cerebellarhemisphere vs. brain_cerebellum

# cutoff value 0.25
spec1.3 = specificity(matrix4.1)
sens1.3 = sensitivity(matrix4.1)

# cutoff value 0.5
spec2.3 = specificity(matrix4)
sens2.3 = sensitivity(matrix4)

# cutoff value 0.75
spec3.3 = specificity(matrix4.2)
sens3.3 = sensitivity(matrix4.2)

plot(c(0.25, 0.5, 0.75), c(sens1.3, sens2.3, sens3.3), type="b", col='red', ylim = (0:1), 
     main = "Sensitivity and specificity for cerebellarhemisphere vs cerebellum", xlab = "cut-off values", ylab = "score")
lines(c(0.25, 0.5, 0.75), c(spec1.3, spec2.3, spec3.3), type="b",col = 'blue')
lines(c(0.25, 0.5, 0.75), c(0.1707317, 0.1707317, 0.1707317), type="b",col = 'black')
legend("bottomright", legend = c("sensitivity", "specificity", "error"), col = c("red", "blue", "black"), lty = c(1))



## KNN 

## brain_amygdala vs brain hippocampus
tissue1 <- "brain_amygdala"
tissue2 <- "brain_hippocampus"
mydata <- data2[which(data2$tissue==tissue1|data2$tissue==tissue2),]
mydata <- droplevels(mydata)

# new coef.var calculation for variable reduction purposes
nr.col3 = ncol(mydata) -1 # last column contains tissue label
mean3 = sapply(mydata[,1:nr.col3], mean)
stand.dev3 = sapply(mydata[,1:nr.col3], sd)
coef.var3 = stand.dev3 / mean3

small2 = mydata[, which(coef.var3 <0.6)]  # threshold selected based on number of variables
small2$tissue = mydata$tissue # to add tissue back to small2

# get test and training data
nr.rows5 <- nrow(small2)
trainindex5 <- sample(1:nr.rows5, size=round(0.8*nr.rows5), replace=FALSE)
train5 <- small2[trainindex5, ]
test5 <- small2[-trainindex5, ]

# get rid of tissue columns
which(names(train5) == "tissue") # to get index

# Make predictions
knn1 = knn(train = train5[,-118], test =test5[,-118], cl = train5$tissue, k=10)

# Print a confusion table
table(knn1, as.factor(test5$tissue))   

# Calculate correct predictions
mean(knn1 == test5$tissue)

# Calculate error
1 - mean(knn1 == test5$tissue)

## brain_hippocampus vs brain_nucleusaccumbens
tissue1 <- "brain_hippocampus"
tissue2 <- "brain_nucleusaccumbens"
mydata <- data2[which(data2$tissue==tissue1|data2$tissue==tissue2),]
mydata <- droplevels(mydata)

# new coef.var calculation for variable reduction purposes
nr.col3 = ncol(mydata) -1 # last column contains tissue label
mean3 = sapply(mydata[,1:nr.col3], mean)
stand.dev3 = sapply(mydata[,1:nr.col3], sd)
coef.var3 = stand.dev3 / mean3

small2 = mydata[, which(coef.var3 <0.6)]  # threshold selected based on number of variables
small2$tissue = mydata$tissue # to add tissue back to small2

# get test and training data
nr.rows5 <- nrow(small2)
trainindex5 <- sample(1:nr.rows5, size=round(0.8*nr.rows5), replace=FALSE)
train5 <- small2[trainindex5, ]
test5 <- small2[-trainindex5, ]

# get rid of tissue columns
which(names(train5) == "tissue") # to get index

# Make predictions
knn2 = knn(train = train5[,-97], test = test5[,-97], cl = train5$tissue, k=10)

# Print a confusion table
table(knn2, as.factor(test5$tissue)) 

# Calculate correct predictions
mean(knn2 == test5$tissue)

# Calculate error
1 - mean(knn2 == test5$tissue)

## brain_spinalcord vs brain_substantianigra
tissue1 <- "brain_spinalcord"
tissue2 <- "brain_substantianigra"
mydata <- data2[which(data2$tissue==tissue1|data2$tissue==tissue2),]
mydata <- droplevels(mydata)

# new coef.var calculation for variable reduction purposes
nr.col3 = ncol(mydata) -1 # last column contains tissue label
mean3 = sapply(mydata[,1:nr.col3], mean)
stand.dev3 = sapply(mydata[,1:nr.col3], sd)
coef.var3 = stand.dev3 / mean3

small2 = mydata[, which(coef.var3 <0.6)]  # threshold selected based on number of variables
small2$tissue = mydata$tissue # to add tissue back to small2

# get test and training data
nr.rows5 <- nrow(small2)
trainindex5 <- sample(1:nr.rows5, size=round(0.8*nr.rows5), replace=FALSE)
train5 <- small2[trainindex5, ]
test5 <- small2[-trainindex5, ]

# get rid of tissue columns
which(names(train5) == "tissue") # to get index

# Make predictions
knn3 = knn(train = train5[,-86], test =test5[,-86], cl = train5$tissue, k=10)

# Print a confusion table
table(knn3, as.factor(test5$tissue))   
# Calculate correct predictions
mean(knn3 == test5$tissue)

# Calculate error
1 - mean(knn3 == test5$tissue)


## brain_cerebellarhemisphere vs brain_cerebellum
tissue1 <- "brain_cerebellarhemisphere"
tissue2 <- "brain_cerebellum"
mydata <- data2[which(data2$tissue==tissue1|data2$tissue==tissue2),]
mydata <- droplevels(mydata)

# new coef.var calculation for variable reduction purposes
nr.col3 = ncol(mydata) -1 # last column contains tissue label
mean3 = sapply(mydata[,1:nr.col3], mean)
stand.dev3 = sapply(mydata[,1:nr.col3], sd)
coef.var3 = stand.dev3 / mean3

small2 = mydata[, which(coef.var3 <0.6)]  # threshold selected based on number of variables
small2$tissue = mydata$tissue # to add tissue back to small2

# get test and training data
nr.rows5 <- nrow(small2)
trainindex5 <- sample(1:nr.rows5, size=round(0.8*nr.rows5), replace=FALSE)
train5 <- small2[trainindex5, ]
test5 <- small2[-trainindex5, ]

# get rid of tissue columns
which(names(train5) == "tissue") # to get index

# Make predictions
knn4 = knn(train = train5[,-225], test =test5[,-225], cl = train5$tissue, k=10)

# Print a confusion table
table(knn4, as.factor(test5$tissue))   
# Calculate correct predictions
mean(knn4 == test5$tissue)

# Calculate error
1 - mean(knn4 == test5$tissue)


#---------------Project day 2-----------------------#
########################################## QUESTION 1
library(randomForest)
library(e1071)
library(gbm) 

#-----------Regression------------------------------------------------------------#
#predict expression levels for gene of interest: ENSG00000271043.1_MTRNR2L2
#three different methods will be used: boosting, randomforest and SVM
#get rid of tissue column
data_notissue = data2[,-400]

#BOOSTING---------------------------------------------
#Parameters chosen: nlambdas = 10, cv fold = 10, ntrees = 100
pows = seq(-10, -0.2, by = 1)
lambdas = 10^pows
length.lambdas = length(lambdas)
train.errors = rep(NA, length.lambdas)
test.errors = rep(NA, length.lambdas)
for (i in 1:length.lambdas) {gbm_model = gbm(data_notissue$ENSG00000271043.1_MTRNR2L2 ~ ., 
                                             data = data_notissue, distribution = "gaussian", n.trees = 100, cv.folds = 10, shrinkage = lambdas[i])
train.errors[i] = mean(gbm_model$train.error)
test.errors[i]  = mean(gbm_model$cv.error)}

#plot of test error
plot(lambdas, test.errors, type = "b", xlab = "Lambda", ylab = "Test error", main='Regression - Errors for Boosting classifier', col = "red", pch = 20)
#find minimum error
min(test.errors)
#corresponding lambda of minimum error
lambdas[which.min(test.errors)]

#Randomforest---------------------------------------------
# optimal mtry for regression: (p/3)-1 
p = dim(data_notissue)[2]
which(colnames(data_notissue) == "ENSG00000271043.1_MTRNR2L2") #check which column the MTRNR2L2 gene is in
data_nogene = data_notissue[,-82] #delete gene of interest from dataset

#create a function for input parameter
parameter = function(p){result = round((p/3)); return(max(1,result))} #created to obtain optimal mtry for regression
#create randomforest model, immediately apply CV
reg_randomforest = rfcv(trainy = data_notissue$ENSG00000271043.1_MTRNR2L2, trainx = data_nogene, mtry = parameter, cv.fold = 10, ntree = 100)
summary(reg_randomforest)
plot(reg_randomforest$n.var, reg_randomforest$error.cv, col = "purple", xlab = "Number of variables", ylab = "Test error", main = 'Regression - Errors for randomforest classifier')
#find minimum error
min(reg_randomforest$error.cv)

#SVM-------------------------------------------------
#costs 0.001,0.01,0.1,1
#linear SVM
tune.out1 = tune(svm, ENSG00000271043.1_MTRNR2L2 ~ ., data=data_notissue, kernel="linear", ranges=list(cost=c(0.001,0.01,0.1,1)))
svm.fit1 = tune.out1$best.model
summary(svm.fit1)
#polynomial SVM
tune.out2 = tune(svm, ENSG00000271043.1_MTRNR2L2 ~ ., data=data_notissue, kernel="polynomial", ranges=list(cost=c(0.001,0.01,0.1,1)))
svm.fit2 = tune.out2$best.model
summary(svm.fit2)
#radial SVM
tune.out3 = tune(svm, ENSG00000271043.1_MTRNR2L2 ~ ., data=data_notissue, kernel="radial", ranges=list(cost=c(0.001,0.01,0.1,1)))
svm.fit3 = tune.out3$best.model
summary(svm.fit3)
#obtain minimal erors of 3 models
min(tune.out1$performances$error)
min(tune.out2$performances$error)
min(tune.out3$performances$error)

tune.out1$performances$cost[which.min(tune.out1$performances$error)]
tune.out2$performances$cost[which.min(tune.out2$performances$error)]
tune.out3$performances$cost[which.min(tune.out3$performances$error)]

#create SVM plot
plot(log10(tune.out1$performances$cost), tune.out1$performances$error, type = 'b', col= 'red', xlab = 'Cost', ylab='Test error', main = 'Regression - Errors for SVM classifiers', ylim = c(10000,400000))
lines(log10(tune.out2$performances$cost), tune.out2$performances$error, type = 'b', col = 'blue')
lines(log10(tune.out3$performances$cost), tune.out3$performances$error, type = 'b', col= 'green')
legend( 'topright',  col = c('red', 'blue', 'green'), legend = c('Linear', 'Polynomial', 'Radial'), lty = c(1))


#-----------Classification------------------------------------------------------------------#
## get dataset containing only two tissuetypes: brain_hippocampus and brain_nucleusaccumbens
tissue1 <- "brain_hippocampus"
tissue2 <- "brain_nucleusaccumbens"
tissuedata <- data2[which(data2$tissue==tissue1|data2$tissue==tissue2),]
tissuedata <- droplevels(tissuedata)


#BOOSTING---------------------------------------------
#Parameters chosen: nlambdas = 10, cv fold = 10, ntrees = 100
pows = seq(-10, -0.2, by = 1)
lambdas = 10^pows
length.lambdas = length(lambdas)
train.errors_tissue = rep(NA, length.lambdas)
test.errors_tissue = rep(NA, length.lambdas)
for (i in 1:length.lambdas) {gbm_model_tissue = gbm(tissue ~ ., data = tissuedata, 
                                                    distribution = "gaussian",  n.trees = 100, cv.folds = 10, shrinkage = lambdas[i])
train.errors_tissue[i] = mean(gbm_model_tissue$train.error)
test.errors_tissue[i]  = mean(gbm_model_tissue$cv.error)}

#plot of test error
plot(lambdas, test.errors_tissue, type = "b", xlab = "Lambda", ylab = "Test error", col = "red", pch = 20, main = 'Classification - Errors for Boosting classifier')
#find minimum of test error
min(test.errors_tissue)
#most optimal lambda corresponding to minimal test error
lambdas[which.min(test.errors_tissue)]


#Randomforest---------------------------------------------
# optimal mtry for classification: sqrt(p)

#create a function for input parameter
parameter2 = function(p){result = round(sqrt(p)); return(max(1,result))}
#create randomforest model with CV
clas_randomforest = rfcv(trainy = tissuedata$tissue, trainx = tissuedata[,-400], mtry = parameter2, cv.fold = 10, ntree = 100, importanceSD = TRUE)
summary(clas_randomforest)
length(clas_randomforest$error.cv)
plot(clas_randomforest$n.var, clas_randomforest$error.cv, col = "purple", xlab = "Number of variables", ylab = "Test error", main = 'Classification - Errors for randomforest classifier')
min(clas_randomforest$error.cv)

#SVM-------------------------------------------------

#linear SVM
tune.out4 = tune(svm, tissue ~ ., data=tissuedata, kernel="linear", ranges=list(cost=c(0.001,0.01,0.1,1)))
svm.fit4 = tune.out4$best.model
summary(svm.fit4)
#polynomial SVM
tune.out5 = tune(svm, tissue ~ ., data=tissuedata, kernel="polynomial", ranges=list(cost=c(0.001,0.01,0.1,1)))
svm.fit5 = tune.out5$best.model
summary(svm.fit5)
#radial SVM
tune.out6 = tune(svm, tissue ~ ., data=tissuedata, kernel="radial", ranges=list(cost=c(0.001,0.01,0.1,1)))
svm.fit6 = tune.out6$best.model
summary(svm.fit6)
#obtain minimal erors of 3 models
min(tune.out4$performances$error)
min(tune.out5$performances$error)
min(tune.out6$performances$error)
tune.out4$performances$cost[which.min(tune.out4$performances$error)]
tune.out5$performances$cost[which.min(tune.out5$performances$error)]
tune.out6$performances$cost[which.min(tune.out6$performances$error)]

#create combined SVM plot
plot(log10(tune.out4$performances$cost), tune.out4$performances$error, type = 'b', col= 'red', xlab = 'Cost', ylab='Test error', main = 'Classification - Errors for SVM classifiers', ylim = c(0, 1))
lines(log10(tune.out5$performances$cost), tune.out5$performances$error, type = 'b', col = 'blue')
lines(log10(tune.out6$performances$cost), tune.out6$performances$error, type = 'b', col= 'green')
legend( 'topright',  col = c('red', 'blue', 'green'), legend = c('Linear', 'Polynomial', 'Radial'), lty = c(1))


########################################## QUESTION 2
#randomForest 
random_forest_clas = randomForest(y = tissuedata$tissue, x = tissuedata[,-400], nvar=6, importance= TRUE)
table_imp =importance(random_forest_clas)
table_imp = as.data.frame(table_imp)
table_imp$genenames = rownames(table_imp)
ordered_gini = table_imp$MeanDecreaseGini[order(table_imp$MeanDecreaseGini, decreasing = TRUE)]
ordered_genes = table_imp$genenames[order(table_imp$MeanDecreaseGini, decreasing = TRUE) ]
cbind(ordered_genes, ordered_gini)        

#Boosting
summary(gbm_model_tissue)


########################################## QUESTION 3 Boosting
#genes for new boosting model: #"ENSG00000171532.4_NEUROD2" and #"ENSG00000104888.5_SLC17A7"
twogenes_data = as.data.frame(cbind(tissuedata$ENSG00000171532.4_NEUROD2, tissuedata$ENSG00000104888.5_SLC17A7, tissuedata$tissue))

pows = seq(-10, -0.2, by = 1)
lambdas = 10^pows
length.lambdas = length(lambdas)
train.errors_twogenes = rep(NA, length.lambdas)
test.errors_twogenes = rep(NA, length.lambdas)
for (i in 1:length.lambdas) {
  gbm_model_twogenes = gbm(tissuedata$tissue ~ ., data = twogenes_data, distribution = "gaussian", 
                           n.trees = 100, cv.folds = 10, shrinkage = lambdas[i])
  train.errors_twogenes[i] = mean(gbm_model_twogenes$train.error)
  test.errors_twogenes[i]  = mean(gbm_model_twogenes$cv.error)}

#plot of test mse
plot(lambdas, test.errors_twogenes, type = "b", xlab = "Lambda", ylab = "Test error", col = "red", pch = 20, main = 'SLC17A7 and NEUROD2 classification - Error for boosting classifier')
#find minimum of MSE
min(test.errors_twogenes)
#most optimal lambda with minimal MSE
lambdas[which.min(test.errors_twogenes)]



########################################## QUESTION 4
#create dataset without two genes of interest
#first get the column numbers of tissues of interest
which(colnames(tissuedata) == "ENSG00000171532.4_NEUROD2")
data_no2genes = tissuedata[,-220]
which(colnames(data_no2genes) == "ENSG00000104888.5_SLC17A7")
data_no2genes = data_no2genes[,-259]
dim(data_no2genes)

pows = seq(-10, -0.2, by = 1)
lambdas = 10^pows
length.lambdas = length(lambdas)
train.errors_no2genes = rep(NA, length.lambdas)
test.errors_no2genes = rep(NA, length.lambdas)
for (i in 1:length.lambdas) { gbm_model_no2genes = gbm(tissuedata$tissue ~ ., 
                                                       data = data_no2genes, distribution = "gaussian", n.trees = 100, cv.folds = 10, shrinkage = lambdas[i])
train.errors_no2genes[i] = mean(gbm_model_no2genes$train.error)
test.errors_no2genes[i]  = mean(gbm_model_no2genes$cv.error)}

#plot of test error
plot(lambdas, test.errors_no2genes, type = "b", xlab = "Lambda", ylab = "Test error", col = "red", pch = 20, main = 'classification without SLC17A7 and NEUROD2 - Error for boosting classifier')
#find minimum test error
min(test.errors_no2genes)
#find corresponding lambda for minimal test error
lambdas[which.min(test.errors_no2genes)]



#---------------Project day 3-----------------------#
########################################## QUESTION 1
#create a new palette with nice colours
palette('default')
palette <- palette(c('pink', 'blue', 'yellow', 'black', 'orange', 'cyan', 'purple', 'brown', 'grey', 'coral', 'magenta', 'red', 'blueviolet' ))

###SCALING or NO SCALING?---------------------------
#PCA no scaling
PCA = prcomp(data_notissue)
PCA_var = PCA$sdev^2
PVE = PCA_var / sum(PCA_var) #calculating proportion of variance
#variance plot
plot(PVE[1:10], xlab = 'Principal Component', ylab = 'Proportion of variance explained', ylim=c(0,1), type = 'b', main = 'Scree plot nonscaled data ')
# create PCA plot ---- PC1 PC2
plot(PCA$x[,1:2], col= data2$tissue, pch = 20, main = 'PCA nonscaled data', xlab ='PC1 (85.37%)', ylab = 'PC2 (6.8%)')
legend('topleft',  col = 1:13, legend = c('Amygdala', 'Anteriorcortex', 'Caudate', 'Cerebellarhemisphere', 'Cerebellum', 'Cortex', 'Frontalcortex', 'Hippocampus', 'Hypothalamus', 'Nucleusaccumbens', 'Putamen', 'Spinalcord', 'Substantianigra' ), pch=16, cex=0.9)
# create Biplot
autoplot(PCA, data=data2, colour = 'tissue', loadings=TRUE)

#PCA scaling
PCA_scale = prcomp(data_notissue, scale = TRUE)
PCA_var_scale = PCA_scale$sdev^2
PVE_scale = PCA_var_scale / sum(PCA_var_scale) #calculating proportion of variance
#variance plot
plot(PVE_scale[1:10], xlab = 'Principal Component', ylab = 'Proportion of variance explained', ylim=c(0,1), type = 'b', main = 'Scree plot scaled data')
# create PCA plot ---- PC1 PC2
plot(PCA_scale$x[,1:2], col= data2$tissue, pch = 20, main = 'PCA scaled data', xlab ='PC1 (37.84%)', ylab = 'PC2 (14.08%)')
# create Biplot
autoplot(PCA_scale, data=data2, colour = 'tissue', loadings=TRUE)

#----------obtain genes that explain PCS----using scaled data --------
#what is the range of loadings
range(PCA_scale$rotation) #since the range is negative, we will use squared loadings from now on 

#Best loadings of scaled PCA
#PC1
loadings_sq_1 <- sort((PCA_scale$rotation[,1])^2, decreasing = TRUE)
loadings_sq_1[1:10]
#PC2
loadings_sq_2 <- sort((PCA_scale$rotation[,2])^2, decreasing = TRUE)
loadings_sq_2[1:10]


########################################## QUESTION 2
#Hierarcial clustering-----------------------------
# choose the most optimal type of linkage
#single  linkage
hc_complete1 = hclust(dist(data_notissue), method="single")
Single_linkage = table(cutree(hc_complete1, 13))
#average linkage
hc_complete2 <- hclust(dist(data_notissue), method="average")
Average_linkage = table(cutree(hc_complete2, 13))
#complete linkage
hc_complete3 = hclust(dist(data_notissue), method="complete")
Complete_linkage = table(cutree(hc_complete3, 13))

table_hclust = cbind(Single_linkage, Average_linkage, Complete_linkage)

#based on the table, we chose complete linkage
table(cutree(hc_complete3, 13), data2$tissue)
barplot(table(data2$tissue, cutree(hc_complete3, 13)), col = 1:13, ylim = c(0,1000), main = 'Tissue type division for H-clustering')
plot(PCA_scale$x[,1:2], col=cutree(hc_complete3, 13), main = 'PCA H-clustering based on complete linkage',xlab =  'PC1 (37.84%)', ylab = 'PC2 (14.08%)', pch=20)

#Kmeans clustering---------------------------------
k_means_result = kmeans(data_notissue, 13, nstart=20)
plot(PCA_scale$x[,1:2], col=k_means_result$cluster, main = 'PCA clustering based on 13-means',xlab =  'PC1 (37.84%)', ylab = 'PC2 (14.08%)', pch=20)
barplot(table(data2$tissue, k_means_result$cluster), col = 1:13, main = 'Tissue type division for 13-means clustering')

#create plot to visualize the amount of samples for each tissue type
plotje = data.frame(tissue = c('Amygdala', 'Anteriorcortex', 'Caudate', 'Cerebellarhemisphere', 'Cerebellum', 'Cortex', 'Frontalcortex', 'Hippocampus', 'Hypothalamus', 'Nucleusaccumbens', 'Putamen', 'Spinalcord', 'Substantianigra' ), count = c(72,84,117,105,125,114,108,94,96,113,97,71,63))    
barplot(plotje[,2], names.arg = plotje[,1], col = 1:13, main = 'Tissue type distribution in reduced dataset', las=2, cex.names=0.8)

########################################## QUESTION 3
#get data of 4PCs based on scree plot of the beginning
pca_data = as.data.frame(PCA$x[,1:4])
pca_data$tissue <- data2$tissue

tissue1 <- "brain_hippocampus"
tissue2 <- "brain_nucleusaccumbens"
tissuedata <- pca_data[which(pca_data$tissue==tissue1|pca_data$tissue==tissue2),]
tissuedata <- droplevels(tissuedata)

#apply boosting
#Parameters chosen: nlambdas = 10, cv fold = 10, ntrees = 100
pows = seq(-10, -0.2, by = 1)
lambdas = 10^pows
length.lambdas = length(lambdas)
train.errors_tissue = rep(NA, length.lambdas)
test.errors_tissue = rep(NA, length.lambdas)
for (i in 1:length.lambdas) {gbm_model_tissue = gbm(tissue ~ ., data = tissuedata, 
                                                    distribution = "gaussian",  n.trees = 100, cv.folds = 10, shrinkage = lambdas[i])
train.errors_tissue[i] = mean(gbm_model_tissue$train.error)
test.errors_tissue[i]  = mean(gbm_model_tissue$cv.error)}

#plot of test error
plot(lambdas, test.errors_tissue, type = "b", xlab = "Lambda", ylab = "Test error", col = "red", pch = 20, main = 'Classification of preprocessed data -  \n Errors for Boosting classifiernon two lines')
#find minimum of test error
min(test.errors_tissue)
#most optimal lambda corresponding to minimal test error
lambdas[which.min(test.errors_tissue)]


