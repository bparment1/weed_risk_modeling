############### SESYNC Research Support: weed risk ########## 
## Performing ROC on data for model assessment.
## 
## DATE CREATED: 06/15/2017
## DATE MODIFIED: 08/17/2017
## AUTHORS: Chinchu Harris 
## PROJECT: weed risk Chinchu Harris
## ISSUE: 
## TO DO:
##
## COMMIT: Benoit modifying and testing earlier code with caret 
##
## Links to investigate:
##https://stats.idre.ucla.edu/r/dae/logit-regression/
#
#
###################################################
#


###### Library used

library(randomForest) #for random forests
library(lattice)
library(ggplot2)
library(caret) #for CV folds and data splitting
library(ROCR) #for diagnostics and ROC plots/stats
library(gplots)
library(pROC) #same as ROCR


###### Functions used in this script and sourced from other files

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

### Other functions ####

function_sampling <- "sampling_function_06292017b.R" #PARAM 1
#function_modeling <- "CH07-26-2017roc_weed_risk_functions_06292017c.R" #PARAM 1 #changed this to another file
#function_modeling <- "CH07-19-2017roc_weed_risk_functions_06292017c.R" #PARAM 1 #changed this to another file
function_modeling <- "roc_weed_risk_functions_08112017d.R" #PARAM 1 #changed this to another file

#script_path <- "/Users/chinchuharris/modeling_weed_risk/scripts" #path to script #PARAM 
script_path <- "/nfs/bparmentier-data/Data/projects/modeling_weed_risk/scripts"

source(file.path(script_path,function_sampling)) #source all functions used in this script 1.
source(file.path(script_path,function_modeling)) #source all functions used in this script 1.

############################################################################
####################  Parameters and argument set up ###########

in_dir <- "/nfs/bparmentier-data/Data/projects/modeling_weed_risk/data"
out_dir <- "/nfs/bparmentier-data/Data/projects/modeling_weed_risk/outputs"

#in_dir <- "/Users/chinchuharris/modeling_weed_risk/data"
#setwd("/Users/chinchuharris/modeling_weed_risk/data")


#Chinchu data
#in_dir <- "/Users/chinchuharris/modeling_weed_risk/data" #local bpy50 , param 1
#out_dir <- "/Users/chinchuharris/modeling_weed_risk/outputs" #param 2

num_cores <- 2 #param 8 #normally I use only 1 core for this dataset but since I want to use the mclappy function the number of cores is changed to 2. If it was 1 then mclappy will be reverted back to the lapply function
create_out_dir_param=TRUE # param 9

out_suffix <-"roc_experiment_08172017" #output suffix for the files and ouptut folder #param 12

#MAJ=major invader=1 #MIN=minor invader=0
#68 instances of 1 #68 instances of 0

#data=read.csv(file="secondmodel.csv") #data for prediction
infile_data <- "publicavailableaphisdatsetforbenoit.csv"

model_names <- c("logistic","randomForest")

##############################  START SCRIPT  ############################

######### PART 0: Set up the output dir ################

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
  
}
#out_dir <- in_dir #output will be created in the input dir

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

options(scipen=999)  #remove scientific writing


### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

data <- read.csv(file.path(in_dir,infile_data))
dim(data)
View(data)
#> dim(data)
#[1] 94 42

###### PART 1: Use the whole dataset and set up equation and model predictions  #####

###Using 43 variables specified in the firstmodel data Appendix###
#data.full<-data[,c("invasion.status", "es1", "es2", "es3","es4", "es5","es6","es7",
#                   "es8","es9","es10","es11","es12","es13","es14","es15","es16","es17",
#                   "es18","es19","es20","es21","es22","es23", "impg1", 
#                   "impg2", "impn1", "impn2", "impn3", "impn4", "impn5", "impn6", 
#                   "impa1", "impa2", "impa3", "impa4", "impp1", "impp2", "impp3",
#                   "impp4", "impp5", "impp6")]


###Converting DV into Factor with names for Caret Library###

data.full <- data
data.full$invasion.status<-factor(data.full$invasion.status, levels=c(0,1), labels=c("ABS", "PRES"))

set.seed(100) #the most metal seed for CV

index_testing <- createDataPartition(data.full$invasion.status,
                    p=0.3,list=F)

data_testing <- data.full[index_testing,]
data_training <- data.full[-index_testing,]


#This method of data slicing-or CV will be used for logistic model
#https://www.analyticsvidhya.com/blog/2016/12/practical-guide-to-implement-machine-learning-with-caret-package-in-r-with-practice-problem/
  
tc <- trainControl(method="repeatedcv", 
                  number=10, 
                  repeats = 10, #note p is default value 0.75, i.e. 75% training
                  summaryFunction=twoClassSummary, 
                  classProb=T)

#https://stackoverflow.com/questions/31138751/roc-curve-from-training-data-in-caret


#creates 10 CV folds for this data, summaryFunction provides a ROC summary stat in call to model#
###Logistic Regression Model Specification###

model.lr.1<-train(as.factor(invasion.status)~es1+es2+es3+es4+es5+es6+es7+es8+es9+es10+es11+es12+es13+es14+es15+es16+es17+es18+es19+es20+es21+es22+es23+impg1+impg2+impn1+impn2+impn3+impn4+impn5+impn6+impa1+impa2+impa3+impa4+impp1+impp2+impp3+impp4+impp5+impp6, metric="ROC", method="glm", family="binomial", trControl=tc, data=data.full)
summary(model.lr.1) #provides coefficients & traditional R model output
model.lr.1 #provides CV summary stats #keep in mind caret takes first class (here the first class is O) ##Generalized Linear Model 

#as reference class so sensitivity and specificity are backwards
confusionMatrix(model.lr.1, norm="average")
###Implementing RF (withCV) on entirety of data###
model.rf<-train(as.factor(invasion.status)~., metric="ROC", method="rf", importance=T, proximity=F, ntree=1000, trControl=tc, data=data.full)
#Variable importance measures are retained when importance=T
#Number of trees grown is indicated on ntree=1000
model.rf
confusionMatrix(model.rf, norm="average")


library(ROCR)
attach(data.full) #have to attach the data in order to get probabilities
###Gathering info for ROC Plots###
LR.1.pred<-predict(model.lr.1, data.full$invasion.status, type="prob")
RF.1.pred<-predict(model.rf, data.full$invasion.status, type="prob")
pred.LR<-prediction(LR.1.pred$MAJ, data$invasion.status)
perf.LR<-performance(pred.LR, "tpr", "fpr")
pred.RF.1<-prediction(RF.1.pred$MAJ, data$invasion.status)
perf.RF.1<-performance(pred.RF.1, "tpr", "fpr")
#plot(perf.LR, main="Logistic Regression and Random Forests")
#plot(perf.RF.1, add=T, lty=2)
#legend(0.32, 0.25, c("Logistic Regression 0.723", "Random Forest  0.885"), lty=c(1,2), bty="n", cex=.75)
par(mar = c(6.5, 6.5, 1.5, 0.5), mgp = c(5, 1, 0))
plot(perf.LR, main="Logistic Regression and Random Forests", lty=1, col=578)
plot(perf.RF.1, xlab = "False Positive Rate", ylab = "True Positive Rate", add=T, lty=3, col=84)
legend(0.50, 0.50, c("Logistic Regression 0.729", "Random Forest  0.877"), lty = c(1,3), col = c(578,84), bty="o", cex=.3)
dev.print(tiff, "secondmodel08-15-2017.tiff", res=600, height=5, width=7, units="in")
###Separation Plots###


library(separationplot)
##Have to transform DV back to 0,1 values for separation plots
data.full$invasion.status<-factor(data.full$invasion.status, levels=c(0,1), labels=c("MIN","MAJ"))
invasion.status<-as.vector(data.full$invasion.status) #transforming actual observations into vector
separationplot(LR.1.pred$MAJ, invasion.status, type="line", line=T, lwd2=1, show.expected=T, heading="Logistic Regression", height=1.5, col0="white", col1="black")
separationplot(RF.1.pred$MAJ, invasion.status, type="line", line=T, lwd2=1, show.expected=T, heading="Random Forests", height=2.5, col0="white", col1="black")
