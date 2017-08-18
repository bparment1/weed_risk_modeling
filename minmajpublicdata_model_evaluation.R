setwd("/Users/chinchuharris/modeling_weed_risk/data")
data=read.csv(file="publicavailableaphisdatsetforbenoit.csv") #data for prediction
#MAJ = 1 Major invader 59 instances; MIN = 0 Minor invaders 35 instances
library(randomForest) #for random forests
library(lattice)
library(ggplot2)
library(caret) #for CV folds and data splitting
library(ROCR) #for diagnostics and ROC plots/stats
library(gplots)
library(pROC) #same as ROCR
###Using 43 variables specified in the firstmodel data Appendix###
data.full<-data[,c("invasion.status", 
                   "es1", "es2", "es3","es4", "es5","es6","es7",
                   "es8","es9","es10","es11","es12","es13","es14",
                   "es15","es16","es17","es18","es19","es20","es21",
                   "es22","es23", "impg1", "impg2", "impn1", "impn2", 
                   "impn3", "impn4", "impn5", "impn6", "impa1", "impa2",
                   "impa3", "impa4", "impp1", "impp2", "impp3", "impp4", 
                   "impp5", "impp6")]
###Converting DV into Factor with names for Caret Library###
data.full$invasion.status<-factor(data.full$invasion.status, 
                                  levels=c(0,1), 
                                  labels=c("MIN", "MAJ"))
set.seed(100) #the most metal seed for CV
#This method of data slicing-or CV will be used for logistic model
#in the trainControl function, the resampling method is "repeatedcv" (repeated cross-validation)
#number = 10 indicates that there are 10 folds in K-fold cross-validation
#repeats = 10 indicates that there are ten separate 10-fold cross-validations used as the resampling scheme
#verboseIter is a logical for printing a training log
#returnData is a logical for saving the data into a slot called trainingData
tc<-trainControl(method="repeatedcv", number=10, 
                 repeats=10, summaryFunction=twoClassSummary, 
                 verboseIter = T, returnData = T, 
                 classProb=T, savePredictions = T)

#creates 10 CV folds for this data, summaryFunction provides a ROC summary stat in call to model#
###Logistic Regression Model Specification###
minmaj.model.lr.1<-train(as.factor(invasion.status)~es1+es2+es3+
                          es4+es5+es6+es7+es8+es9+es10+es11+es12+
                          es13+es14+es15+es16+es17+es18+es19+es20+
                          es21+es22+es23+impg1+impg2+impn1+impn2+impn3+
                          impn4+impn5+impn6+impa1+impa2+impa3+impa4+impp1+
                          impp2+impp3+impp4+impp5+impp6, 
                        metric="ROC", method="glm", 
                        family="binomial", trControl=tc, 
                        data=data.full)

minmaj.model.lr.1$finalModel #gives model coefficients for each predictor variable
minmaj.model.lr.1$trainingData #gives the training log saved from the returnData argument = TRUE but can't access the entire training log
out_of_fold_predictions<-minmaj.model.lr.1$pred
out_of_fold_predictions #gives the prediction values for the glm model 
summary(minmaj.model.lr.1) #provides coefficients & traditional R model output
minmaj.model.lr.1 #provides CV summary stats #keep in mind caret takes minmaj class (here the minmaj class is O) ##Generalized Linear Model 
#as reference class so sensitivity and specificity are backwards
confusionMatrix(minmaj.model.lr.1, norm="average")
###Implementing RF (withCV) on entirety of data###
minmaj.model.rf<-train(as.factor(invasion.status)~., 
                      metric="ROC", method="rf", 
                      importance=T, proximity=F, 
                      ntree=1000, trControl=tc, 
                      data=data.full)
#Variable importance measures are retained when importance=T
#Number of trees grown is indicated on ntree=1000
minmaj.model.rf$finalModel
minmaj.model.rf$trainingData
minmaj.rf_out_of_fold_predictions<-minmaj.model.rf$pred
minmaj.rf_out_of_fold_predictions #gives the prediction values for the rf model
minmaj.model.rf
confusionMatrix(minmaj.model.rf, norm="average")
library(ROCR)
attach(data.full) #have to attach the data in order to get probabilities
###Gathering info for ROC Plots###
minmaj.LR.1.pred<-predict(minmaj.model.lr.1, data.full$invasion.status, type="prob")
minmaj.RF.1.pred<-predict(minmaj.model.rf, data.full$invasion.status, type="prob")
minmaj.pred.LR<-prediction(minmaj.LR.1.pred$MAJ, data$invasion.status)
minmaj.perf.LR<-performance(minmaj.pred.LR, "tpr", "fpr")
minmaj.pred.RF.1<-prediction(minmaj.RF.1.pred$MAJ, data$invasion.status)
minmaj.perf.RF.1<-performance(minmaj.pred.RF.1, "tpr", "fpr")
par(mar = c(7.5, 9.5, 1.5, 3.5), mgp = c(5, 1, 0))
plot(minmaj.perf.LR, 
     main="Logistic Regression and Random Forests", 
     lty=1, col=578)
plot(minmaj.perf.RF.1, xlab = "False Positive Rate", 
     ylab = "True Positive Rate", 
     add=T, lty=3, col=84)
legend(0.60, 0.20, c("Logistic Regression 0.841", 
                     "Random Forest  0.948"), 
       lty = c(1,3), col = c(578,84), 
       bty="o", cex=.60)
dev.print(tiff, "minmajpublicdata08-18-17.tiff", res=600, height=5, width=7, units="in")


###Separation Plots###
library(separationplot)
##Have to transform DV back to 0,1 values for separation plots
data.full$invasion.status<-factor(
  data.full$invasion.status, levels=c("MIN","MAJ"), 
  levels=c(0, 1))

Invasion.status<-as.vector(
  data.full$invasion.status) 
#transforming actual observations into vector
separationplot(LR.1.pred$MAJ, Invasion.status, type="line", line=T, lwd2=1, show.expected=T, heading="Logistic Regression", height=1.5, col0="white", col1="black")
separationplot(RF.1.pred$MAJ, Invasion.status, type="line", line=T, lwd2=1, show.expected=T, heading="Random Forests", height=2.5, col0="white", col1="black")
