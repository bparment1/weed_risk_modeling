############### SESYNC Research Support: weed risk ########## 
## Checking ROC output and identification of CARET training and testing inputs/outputs
## 
## DATE CREATED: 08/16/2017
## DATE MODIFIED: 09/25/2017
## AUTHORS: Benoit Parmentier and Chinchu Harris
## PROJECT: weed risk Chinchu Harris
## ISSUE: 
## TO DO:
##
## COMMIT:  
##
## Links to investigate:

###################################################
#

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gdata)                               # various tools with xls reading, cbindX
library(rasterVis)                           # Raster plotting functions
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(plyr)                                # Various tools including rbind.fill
library(spgwr)                               # GWR method
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
library(ggplot2)
library(lubridate)
library(dplyr)
library(ROCR) #for diagnostics and ROC plots/stats
library(pROC)
library(TOC)
library(randomForest) #for random forests
library(lattice)
library(caret) #for CV folds and data splitting
library(gplots)

#
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

#Chinchu data
#in_dir <- "/Users/chinchuharris/modeling_weed_risk/data" #local bpy50 , param 1
#out_dir <- "/Users/chinchuharris/modeling_weed_risk/outputs" #param 2

num_cores <- 2 #param 8 #normally I use only 1 core for this dataset but since I want to use the mclappy function the number of cores is changed to 2. If it was 1 then mclappy will be reverted back to the lapply function
create_out_dir_param=TRUE # param 9

out_suffix <-"roc_experiment_09212017" #output suffix for the files and ouptut folder #param 12

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
#dim(data)
#View(data)
#> dim(data)
#[1] 94 42

###### PART 1: Use the whole dataset and set up equation and model predictions  #####

#head(data)

#data$invasion.status <- as.factor(data$invasion.status)
#y_var <- data$invasion.status
#explanatory_variables <- names(data)[-1]

#setwd("/Users/chinchuharris/modeling_weed_risk/data")

#in_dir <- "/nfs/bparmentier-data/Data/projects/modeling_weed_risk/data"
#setwd(in_dir)

#data=read.csv(file="publicavailableaphisdatsetforbenoit.csv") #data for prediction
#MAJ = 1 Major invader 59 instances; MIN = 0 Minor invaders 35 instances

###Using 43 variables specified in the firstmodel data Appendix###
#data.full<-data[,c("invasion.status", 
#                   "es1", "es2", "es3","es4", "es5","es6","es7",
#                   "es8","es9","es10","es11","es12","es13","es14",
#                   "es15","es16","es17","es18","es19","es20","es21",
#                   "es22","es23", "impg1", "impg2", "impn1", "impn2", 
#                   "impn3", "impn4", "impn5", "impn6", "impa1", "impa2",
#                   "impa3", "impa4", "impp1", "impp2", "impp3", "impp4", 
#                   "impp5", "impp6")]

data.full <- data
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

tc<-trainControl(method="repeatedcv", 
                 number=10, 
                 repeats=10, summaryFunction=twoClassSummary, 
                 verboseIter = T, returnData = T, 
                 classProb=T, savePredictions = T)

y_var <- "invasion.status"
#data[[y_var]] <- as.factor(data[[y_var]]) #this is needed for randomForest to get a classification

explanatory_variables <- names(data)[-1] #drop the first column

right_side_formula <- paste(explanatory_variables,collapse = " + ")
model_formula_str <- paste0(y_var," ~ ",right_side_formula)

#Check the option binomial effect on predictions
#mod_glm <- glm(model_formula_str,data = data,family=binomial())
#mod_glm2 <- glm(model_formula_str,data = data, family="binomial") #same as above
#predicted_val_glm <- predict(mod_glm,data=data)

#creates 10 CV folds for this data, summaryFunction provides a ROC summary stat in call to model#
###Logistic Regression Model Specification###
minmaj.model.lr.1 <-train(as.formula(model_formula_str), 
                         metric="ROC", 
                         method="glm", 
                         family="binomial", 
                         trControl=tc, 
                         data=data.full)

minmaj.model.lr.1$finalModel #gives model coefficients for each predictor variable
minmaj.model.lr.1$trainingData #gives the training log saved from the returnData argument = TRUE but can't access the entire training log
out_of_fold_predictions<-minmaj.model.lr.1$pred
dim(out_of_fold_predictions)

###
minmaj.model.lr.1$control$index #this is a list of training index row number

##### Let's selected each training dataset used in 10 folds, replication 10

#For fold 01, replication 01

index_selected <- minmaj.model.lr.1$control$index$Fold01.Rep01
index_selected <- minmaj.model.lr.1$control$index$Fold10.Rep10

data_training_01_01 <- data.full[index_selected,]
dim(data_training_01_01)

#Check the option binomial effect on predictions
mod_glm <- glm(model_formula_str,
               data = data_training_01_01,
               family="binomial")

index_selected_out <- minmaj.model.lr.1$control$indexOut$Resample001
data_testing_01_01 <- data.full[index_selected_out,]

#mod_glm2 <- glm(model_formula_str,data = data, family="binomial") #same as above
#predicted_val_glm <- predict(mod_glm,data=data)

attach(data.full) #have to attach the data in order to get probabilities

###Gathering info for ROC Plots###
#get porbability of prediction
minmaj.LR.1.pred_testing_01_01 <-predict(mod_glm, 
                                         data_testing_01_01, 
                                         type="prob")

minmaj.LR.1.pred_testing <-predict(mod_glm, 
                                         data_testing_01_01)#, 
                                         #type="prob")

minmaj.LR.1.pred_testing_01_01 <-predict(mod_glm, 
              data_testing_01_01, 
              type="response")

#minmaj.RF.1.pred<-predict(minmaj.model.rf, data.full$invasion.status, type="prob")
#minmaj.pred.LR<-prediction(minmaj.LR.1.pred$MAJ, data$invasion.status)

#minmaj.pred.LR.1_testing_01_01 <-prediction(minmaj.LR.1.pred_testing_01_01$PRES, 
#                             data_testing_01_01$invasion.status)

minmaj.pred.LR.1_testing_01_01 <-prediction(minmaj.LR.1.pred_testing_01_01, 
                                            data_testing_01_01$invasion.status)

minmaj.perf.LR<-performance(minmaj.pred.LR.1_testing_01_01,
                            measure="tpr",
                            x.measure = "fpr")
auc.perf = performance(minmaj.pred.LR.1_testing_01_01, measure = "auc")

par(mar = c(7.5, 9.5, 1.5, 3.5), mgp = c(5, 1, 0))
plot(minmaj.perf.LR, 
     main="Logistic Regression and Random Forests", 
     lty=1, col=578)
#plot(minmaj.perf.RF.1, xlab = "False Positive Rate", 
#     ylab = "True Positive Rate", 
#     add=T, lty=3, col=84)
legend(0.60, 0.20, c("Logistic Regression 0.841", 
                     "Random Forest  0.948"), 
       lty = c(1,3), col = c(578,84), 
       bty="o", cex=.60)
dev.print(tiff, "minmajpublicdata08-18-17.tiff", res=600, height=5, width=7, units="in")
#??performance

## USE TOC package ROC:

#breaks_val <- seq(0,1,0.1)

#hist(mod_glm$fitted.values,breaks=breaks_val)
#barplot(table(mod_glm$fitted.values))       

#mask_val <- 1:nrow(data)
#rocd2 <- ROC(index=modrf, boolean=data[[y_var]], mask=mask_val, nthres=100)
## Select column 2 of predicted matrix of probabilities

y_ref <- as.numeric(data_testing_01_01[[y_var]]) #boolean reference values
y_ref[y_ref==2]<- 1
y_ref[y_ref==1]<- 0

index_val <- minmaj.LR.1.pred_testing

mask_val <- rep(1,length=length(index_val))
rocd2_rf <- ROC(index=index_val, 
                boolean=y_ref, 
                mask=mask_val,
                nthres=100)

slot(rocd2_rf,"AUC") #this is your AUC from the logistic modeling
#Plot ROC curve:
plot(rocd2_rf,
     main=model_names[2])

names(rocd2_rf)
str(rocd2_rf)

#Access table: 
roc_table_rf <- slot(rocd2_rf,"table")

### Save ROC plot in a PNG file for glm model

mod_glm$fitted.values #these are the probability values from ROC

out_suffix_str <- paste0("full_data_",model_names[1],out_suffix)

res_pix<-480 #set as function argument...
col_mfrow<-1
#row_mfrow<-2
row_mfrow<-1

png_file_name<- paste("Figure_","ROC_plot_",out_suffix_str,".png", sep="")

png(filename=file.path(out_dir,png_file_name),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(row_mfrow,col_mfrow))

mask_val <- 1:nrow(data)
index_val <- mod_glm$fitted.values

rocd2_glm <- ROC(index= index_val, 
                 boolean=y_ref, 
                 mask=mask_val, 
                 nthres = 100)

slot(rocd2_glm,"AUC") #this is your AUC from the logistic modeling
#Plot ROC curve:
plot(rocd2_glm)

dev.off()

###
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



##################### End of script ##############################