library(Hmisc) # used to impute median
library(rattle) # used for fancy plotting 
library(dplyr)
library(rpart)
library(party)
library(forcats)
library(MLmetrics)
library(earth)
library(caret)
library(plyr)
library(xgboost)
library(MLeval)
install.packages("gbm")
library(ROCR)
library(pROC)
library(auprc)
library(InformationValue)
library(gbm)
# reading the train data
train_data <- read.csv(file = "D:/IDA/HW7/Train.csv")
str(train_data)

#reading test data
test_data <- read.csv(file = "D:/IDA/HW7/Test.csv")
str(test_data)

#pushing readmitted values to NA in test data
test_data$readmitted <- NA

# combining both train and test data
full_data <- train_data %>% rbind(test_data)

# If we have any blank spaces, replacing them with NA to the full data
full_data[full_data== ""] <- NA

# fetching the column names from 22 to 44
col<-full_data[,c(22:44)]

# removing those columns
full_data <- full_data[,-c(22:44)]

# removing medical_speciality and payer_code from full_data
full_data <- select(full_data, -c(medical_specialty, payer_code)) 

str(full_data)
sapply(full_data, function(x) sum(is.na(x)))

# Imputing the values for race,diagnosis,gender
full_data$race<- fct_explicit_na(full_data$race, na_level = mode(full_data$race))
full_data$diagnosis<- fct_explicit_na(full_data$diagnosis, na_level = mode(full_data$diagnosis))
full_data$gender<- fct_explicit_na(full_data$gender, na_level = mode(full_data$gender))

#Handling race
full_data$race <- ifelse (full_data$race == "AfricanAmerican", 1, full_data$race)
full_data$race <- ifelse (full_data$race == "Asian", 2, full_data$race)
full_data$race <- ifelse (full_data$race == "Caucasian", 3, full_data$race)
full_data$race <- ifelse (full_data$race == "Hispanic", 4, full_data$race)
full_data$race <- ifelse (full_data$race == "Other", 5, full_data$race)
full_data$race <- ifelse (full_data$race == "numeric", 6, full_data$race)

#performing Unclass to convert to numerical from factor
full_data$diagnosis <- unclass(full_data$diagnosis)
# changing the values of diagnosis to some levels like 1,2,3 till 17 based on the range of the value
full_data$diagnosis <- ifelse(full_data$diagnosis > 0 & full_data$diagnosis < 140 ,  1, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 139 & full_data$diagnosis < 240 ,  2, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 239 & full_data$diagnosis < 280 ,  3, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 279 & full_data$diagnosis < 290 ,  4, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 289 & full_data$diagnosis < 320 ,  5, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 319 & full_data$diagnosis < 390 ,  6, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 389 & full_data$diagnosis < 460 ,  7, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 459 & full_data$diagnosis < 520 ,  8, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 519 & full_data$diagnosis < 580 ,  9, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 579 & full_data$diagnosis < 630 ,  10, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 629 & full_data$diagnosis < 680 ,  11, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 679 & full_data$diagnosis < 710 ,  12, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 719 & full_data$diagnosis < 740 ,  13, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 739 & full_data$diagnosis < 760 ,  14, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 759 & full_data$diagnosis < 780 ,  15, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 779 & full_data$diagnosis < 800 ,  16, full_data$diagnosis)
full_data$diagnosis <- ifelse(full_data$diagnosis > 799 & full_data$diagnosis < 1000 ,  17, full_data$diagnosis)
#table(full_data$diagnosis)

# Handling the age 
# Replacing the value of age with a value which is maximum in the specified range
full_data$age <- ifelse(full_data$age == "[0-10)",  0, full_data$age)
full_data$age <- ifelse(full_data$age == "[10-20)", 10, full_data$age)
full_data$age <- ifelse(full_data$age == "[20-30)", 20, full_data$age)
full_data$age <- ifelse(full_data$age == "[30-40)", 30, full_data$age)
full_data$age <- ifelse(full_data$age == "[40-50)", 40, full_data$age)
full_data$age <- ifelse(full_data$age == "[50-60)", 50, full_data$age)
full_data$age <- ifelse(full_data$age == "[60-70)", 60, full_data$age)
full_data$age <- ifelse(full_data$age == "[70-80)", 70, full_data$age)
full_data$age <- ifelse(full_data$age == "[80-90)", 80, full_data$age)
full_data$age <- ifelse(full_data$age == "[90-100)", 90, full_data$age)

# Imputing the values for age
full_data$age<- fct_explicit_na(full_data$age, na_level = mode(full_data$age))

# Handling the max_glu_serum
full_data$max_glu_serum <- ifelse(full_data$max_glu_serum == "None",  0, full_data$max_glu_serum)
full_data$max_glu_serum <- ifelse(full_data$max_glu_serum == "Norm",  100, full_data$max_glu_serum)
full_data$max_glu_serum <- ifelse(full_data$max_glu_serum == ">200",  200, full_data$max_glu_serum)
full_data$max_glu_serum <- ifelse(full_data$max_glu_serum == ">300",  300, full_data$max_glu_serum)

#Handling the A1Cresult
full_data$A1Cresult <- ifelse(full_data$A1Cresult == "None",  0, full_data$A1Cresult)
full_data$A1Cresult <- ifelse(full_data$A1Cresult == "Norm",  5, full_data$A1Cresult)
full_data$A1Cresult <- ifelse(full_data$A1Cresult == ">7",    7, full_data$A1Cresult)
full_data$A1Cresult <- ifelse(full_data$A1Cresult == ">8",    8, full_data$A1Cresult)
#sapply(full_data, function(x) sum(is.na(x)))

# handling indicator_level
full_data$indicator_level[is.na(full_data$indicator_level)] <- 0
full_data$indicator_level <- ifelse(full_data$indicator_level < 0 ,  1, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 0 & full_data$indicator_level < 10 ,  1, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 10 & full_data$indicator_level < 20 ,  2, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 20 & full_data$indicator_level < 30 ,  3, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 30 & full_data$indicator_level < 40 ,  4, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 40 & full_data$indicator_level < 50 ,  5, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 50 & full_data$indicator_level < 60 ,  6, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 60 & full_data$indicator_level < 70 ,  7, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 70 & full_data$indicator_level < 80 ,  8, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 80 & full_data$indicator_level < 90 ,  9, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 90 & full_data$indicator_level < 100 ,  10, full_data$indicator_level)
full_data$indicator_level <- ifelse(full_data$indicator_level >= 100,  10, full_data$indicator_level)

# handling missing values in the column time_in_hospital and num_lab_procedures.
# Replacing those NA values with the mode of that column
full_data$time_in_hospital[is.na(full_data$time_in_hospital)] <- 1
full_data$num_lab_procedures[is.na(full_data$num_lab_procedures)] <- 1

#Converting numerical to factor
full_data$admission_type <- as.factor(full_data$admission_type)
full_data$discharge_disposition <- as.factor(full_data$discharge_disposition)
full_data$admission_source <- as.factor(full_data$admission_source)
full_data$diagnosis <-as.factor(full_data$diagnosis)
full_data$A1Cresult <- as.factor(full_data$A1Cresult)
full_data$diabetesMed <- as.factor(full_data$diabetesMed)
full_data$indicator_level <-as.factor(full_data$indicator_level)
full_data$race <-as.factor(full_data$race)
#glimpse(full_data)
#table(full_data$num_lab_procedures)

#splitting the full_data into Train and Test w.r.t readmitted column
Train <- full_data %>% filter(!is.na(readmitted))
Test <- full_data %>% filter(is.na(readmitted))

#--------------------------Modelling-----------------------------------------------------------

#--------------------------Decision Tree--------------------------------------------------------

fit<-rpart(Train$readmitted~.,data=Train,control=rpart.control(cp=0.001),xval=20)
cp.value <- fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]
pfit<-prune(fit,cp=cp.value)  #and we can prune to this level
summary(pfit)
fancyRpartPlot(pfit)

pred<-predict(pfit, newdata=Train)
LogLoss(y_pred=pred, y_true=Train$readmitted)
# logloss=0.6479

decision_Pred_Test<-predict(pfit, newdata=Test)     
TestResult<-tibble(patientID=Test$patientID, predReadmit=decision_Pred_Test)
write.csv(TestResult, 'submission1.csv', row.names = F)

# Confusion Matrix
threshold <- 0.5
confusionMatrix(factor(pred>threshold), factor(Train$readmitted==1), positive="TRUE")

varImp(fit)
#-----------------------------------------------------------------------------------------------------

#-----------------------------------------gbm---------------------------------------------------------

traincontrol <- trainControl(method = "repeatedcv", number = 5, repeats = 5, 
                             verboseIter = FALSE, allowParallel = TRUE)

gbmFit <- train(readmitted ~ ., method = "gbm", maximize = FALSE, 
                trControl = traincontrol, tuneGrid = expand.grid(n.trees = 250,
                                                                 interaction.depth = c(15), shrinkage = c(0.05), n.minobsinnode = c(10)), 
                data = Train, verbose = FALSE)

summary(gbmFit)
#model= gbmFit
#saveRDS(model, "GBM_model") 
#gbm1 <-readRDS("GBM_model") # load the model anytime/anywhere, no need to run it again
gbmpredict <- predict(gbm1, Train)

############ evaluation metrics #################################################

#Confusion matrix
threshold <- 0.5
cm<-confusionMatrix(factor(gbmpredict>threshold), factor(Train$readmitted==1), positive="TRUE")

#Precision
cm$byClass["Precision"]
#precision=0.6589

#Recall
cm$byClass["Recall"]
#recall=0.5599

#F1 score
cm$byClass["F1"]
#F1=0.6053

# precision recall curve
pred <- prediction(dummyTrain$gbmpredict,dummyTrain$readmitted)
perf <- performance(pred,"prec","rec")
plot(perf)

# AUC-ROC curve
test_roc = roc(Train$readmitted ~ gbmpredict, plot = TRUE, print.auc = TRUE)
as.numeric(test_roc$auc)
#0.717057

#ROC curve
roc_score = roc(Train$readmitted, gbmpredict)  # AUC score
plot(roc_score, main="ROC curve -- GBM ")

#Kolmogorov-Smirnov chart
ks_stat(Train$readmitted, gbmpredict)
ks_stat(Train$readmitted, gbmpredict, returnKSTable = T)
ks_plot(Train$readmitted, gbmpredict)
####################################################################################################

gbmpredict <- predict(gbmFit, Train)
gbmpredict <- ifelse(gbmpredict<0, 0, gbmpredict)
gbmpredict <- ifelse(gbmpredict>1, 1, gbmpredict)

LogLoss(y_pred=gbmpredict, y_true=Train$readmitted)
#logloss=0.6156328

gbm_Pred_Test <- predict(gbmFit, Test)
gbm_Pred_Test <- ifelse(gbm_Pred_Test<0, 0, gbm_Pred_Test)
gbm_Pred_Test <- ifelse(gbm_Pred_Test>1, 1, gbm_Pred_Test)

TestResult<-tibble(patientID=Test$patientID, predReadmit=gbm_Pred_Test)
write.csv(TestResult, 'submission3_gbm.csv', row.names = F)

################################################ PLOTS ###############################################
# rbinding the predicted readmitted values to the train dataset
dummyTrain <- Train
dummyTrain <- cbind(dummyTrain, gbmpredict)

#varimp.gbm(gbm1)

#1 Number of diagnosis vs readmitted values 
ggplot(dummyTrain, aes(y= gbmpredict, x = number_diagnoses), scale= TRUE) +
  geom_point(alpha = 0.7) 

#2 Age vs readmitted values 
ggplot(dummyTrain, aes(x = dummyTrain$num_medications , y = gbmpredict)) + 
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=FALSE)+
  theme_ipsum()

# 3 No of medications vs readmitted values 
dummyTrain %>%
  filter( Train$num_medications<35 ) %>%
  ggplot( aes(x=gbmpredict)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)


#---------------------------------------------Logistic------------------------------------------------

fitControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
lr_model <- train(readmitted ~ .,
                  data = Train,
                  method = 'glm', 
                  family = binomial(),
                  trControl = fitControl)

lrpredict_Train <- predict(lr_model, Train)
lrpredict_Train
LogLoss(y_pred=lrpredict_Train, y_true=Train$readmitted)
#logloss=0.6402

lrpredict_Test<-predict(lr_model, newdata=Test)     
lrpredict_Test

TestResult<-tibble(patientID=Test$patientID, predReadmit=lrpredict_Test)
write.csv(TestResult, 'submission2.csv', row.names = F)

#confusion matrix
threshold <- 0.5
confusionMatrix(factor(lrpredict_Train>threshold), factor(Train$readmitted==1), positive="TRUE")

#-------------------------------------------Mars------------------------------------------------------

marsModel <- earth(readmitted ~. ,data = Train,
                   degree=8,
                   nk=80,
                   pmethod="cv",
                   nfold=8,
                   ncross=8)

plotmo(marsModel)
summary(marsModel)

marsModel$fitted.values
marsModel$fitted.values<- ifelse(marsModel$fitted.values>1, 1, marsModel$fitted.values)
marsModel$fitted.values<- ifelse(marsModel$fitted.values<0, 0, marsModel$fitted.values)
mars_pred_Train <- predict(marsModel, newdata = Train)

LogLoss(y_pred = marsModel$fitted.values, y_true = Train$readmitted)
#logloss=0.6445


mars_Pred_Test<-predict(marsModel, newdata=Test)  

TestResult<-tibble(patientID=Test$patientID, predReadmit=mars_Pred_Test)
write.csv(TestResult, 'submission3.csv', row.names = F)

#Confusion matrix
threshold <- 0.5
confusionMatrix(factor(mars_pred_Train>threshold), factor(Train$readmitted==1), positive="TRUE")

#------------------------------------xgboost using caret----------------------------------------------

xgb_Train <- Train
t <- xgb_Train
xgb_Test <- Test[,-c(21)]
xgb_Train[] <- lapply(xgb_Train, as.numeric);
xgb_Test[] <- lapply(xgb_Test, as.numeric);
xbg <- t[,-20]
s <- t[,20]
xgb_D_Train <- xgb.DMatrix(as.matrix(xgb_Train[,-21]), label = xgb_Train[,21]);
xgb_D_Test <- xgb.DMatrix(as.matrix(xgb_Test))

watchlist <- list(xgb_Train = xgb_D_Train);

param <- list(
  objective           = "reg:logistic",
  booster             = "gbtree",
  eta                 = 0.03,
  max_depth           = 5,
  eval_metric         = "auc",
  min_child_weight    = 150,
  alpha               = 0.00,
  subsample           = 0.70,
  colsample_bytree    = 0.70
);

set.seed(1981);

xgb <- xgb.train(params                = param,
                 data                  = xgb_D_Train,
                 nrounds               = 20000,
                 verbose               = 1,
                 watchlist             = watchlist,
                 maximize              = TRUE,
                 nfold                 = 5,
                 nthread               = 4,
                 print_every_n         = 50,
                 stratified            = TRUE,
                 early_stopping_rounds = 10)
x<- xgb_Train
x <- x[-22]
x <-x[-21]
xgb_Train<- xgb_Train[-22]
xgb_Pred <- predict(xgb, data.matrix(xgb_Train[-21]))

LogLoss(y_pred=xgb_Pred, y_true=xgb_Train$readmitted)
#logloss=0.56674

xgb_pred_Test <- predict(xgb, data.matrix(xgb_Test[]))

xgb_Test$patientID
TestResult<-tibble(patientID=xgb_Test$patientID, predReadmit=xgb_pred_Test)
write.csv(TestResult, 'submission5.csv', row.names = F)

#confusion Matrix
threshold <- 0.5
confusionMatrix(factor(xgb_Pred>threshold), factor(Train$readmitted==1), positive="TRUE")

#--------------------------------------------End--------------------------------------------------------
