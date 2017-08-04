############### SESYNC Research Support: weed risk ########## 
## Performing ROC on data for model assessment.
## 
## DATE CREATED: 06/15/2017
## DATE MODIFIED: 08/04/2017
## AUTHORS: Benoit Parmentier, Chinchu Harris 
## PROJECT: weed risk Chinchu Harris
## ISSUE: 
## TO DO:
##
## COMMIT:testing ROC and models of different samples for accuracy assessment
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
library(ROCR)
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
function_modeling <- "CH07-26-2017roc_weed_risk_functions_06292017c.R" #PARAM 1 #changed this to another file
#function_modeling <- "CH07-19-2017roc_weed_risk_functions_06292017c.R" #PARAM 1 #changed this to another file

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

out_suffix <-"roc_experiment_08032017" #output suffix for the files and ouptut folder #param 12

infile_data <- "publicavailableaphisdatsetforbenoit.csv"
#infile_genes_identity <- "genes_identity.csv"

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

head(data)

#data$invasion.status <- as.factor(data$invasion.status)
#y_var <- data$invasion.status
#explanatory_variables <- names(data)[-1]

y_var <- "invasion.status"
data[[y_var]] <- as.factor(data[[y_var]])

explanatory_variables <- names(data)[-1]

right_side_formula <- paste(explanatory_variables,collapse = " + ")
model_formula_str <- paste0(y_var," ~ ",right_side_formula)

#Check the option binomial effect on predictions
mod_glm <- glm(model_formula_str,data = data,family=binomial())
#mod_glm <- glm(model_formula_str,data = data)

breaks_val <- seq(0,1,0.1)

hist(mod_glm$fitted.values,breaks=breaks_val)
       
######Random Forest model##################

#This does not work:
#set.seed(100)

mod_rf <- randomForest(as.formula(model_formula_str),
                  type="classification",
                  data=data,
                  importance = TRUE, 
                  ntree = 10001, 
                  proximity=TRUE) 

#modrf <- predict(mod, data=data) #need to ask for probability otherwise the output is {0,1}

predicted_rf_mat <- predict(mod, data=data, type="prob")

#importance is set to True in order to assess the importance of the predictors (for downstream evaluations like variable importance plots)
#proximity is set to False because we are not determining the distance of the observations to one another
#confusion is set to False because we are using cross validation not OOB error rates to determine accuracy prediction
#err.rate is set to False because we are using cross validation not OOB measures

######## Random Forest model end#############


mask_val <- 1:nrow(data)
#rocd2 <- ROC(index=modrf, boolean=data[[y_var]], mask=mask_val, nthres=100)
## Select column 2 of predicted matrix of probabilities

y_ref <- as.numeric(as.character(data[[y_var]])) #boolean reference values
index_val <- predicted_rf_mat[,2] #probabilities
rocd2_rf <- ROC(index=index_val, 
                boolean=y_ref, 
                mask=mask_val, nthres=100)

#modrf[,1]
#rocd2 <- ROC(index=mod$fitted.values, boolean=data[[y_var]], mask=mask_val, nthres = 100)

slot(rocd2_rf,"AUC") #this is your AUC from the logistic modeling
#Plot ROC curve:
plot(rocd2_rf,main="Random Forest model")

names(rocd2_rf)
str(rocd2_rf)

#Access table: 
roc_table_rf <- slot(rocd2_rf,"table")

### Save ROC plot in a PNG file for glm model

mod_glm$fitted.values #these are the probability values from ROC


out_suffix_str <- paste0("full_data_",out_suffix)
res_pix<-480 #set as function argument...
col_mfrow<-1
#row_mfrow<-2
row_mfrow<-1

png_file_name<- paste("Figure_","ROC_plot_",out_suffix_str,".png", sep="")

png(filename=file.path(out_dir,png_file_name),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(row_mfrow,col_mfrow))

mask_val <- 1:nrow(data)
rocd2_glm <- ROC(index=mod_glm$fitted.values, 
             boolean=y_ref, mask=mask_val, nthres = 100)

slot(rocd2_glm,"AUC") #this is your AUC from the logistic modeling
#Plot ROC curve:
plot(rocd2_glm)

dev.off()


#### Generating plots with multiple models

#Access table: 
roc_table_rf <- slot(rocd2_rf,"table")
#Access table: 
roc_table_glm <- slot(rocd2_glm,"table")

plot(roc_table_rf$falseAlarms1,roc_table_rf$Model1,type="b",pch=7,col="red")
lines(roc_table_glm$falseAlarms1,roc_table_glm$Model1,type="b",pch=7,col="blue")
title("Model comparison Random Forest and GLM")

############# PART 2: Conduct modeling with training and testing data ##############

### First get the samples

seed_number <- 100
nb_sample <- 10
step <- 0.1
prop_minmax <- c(0.3,0.5)
out_suffix
out_dir

#sampling_training_testing(seed_number,nb_sample,step,prop_minmax,data_df,out_suffix,out_dir)
#debug(sampling_training_testing)
sampled_data_obj <- sampling_training_testing(data,nb_sample,step,prop_minmax,obs_id=NULL,seed_number=100,out_suffix="",out_dir=".")

names(sampled_data_obj)

sampled_data_obj$sampling_dat #sampling run summary data.frame with ID and settings
length(sampled_data_obj$data_training)
sampled_data_obj$data_training[[1]] # testing df for run sample ID 1
sampled_data_obj$data_testing[[1]]  # training df for run sampling 2

list_data_training <- sampled_data_obj$data_training
list_data_testing <- sampled_data_obj$data_testing

#paste(data_v)
#names()
###### Can repeat the logistic model for each sample of training!!!
#lf <- run_model_fun(data_df=data,model_formula_str = model_formula_str,model_opt="logistic",data_testing=NULL,save_fig=T,out_dir=".",out_suffix="")
#lf <- run_model_fun(data_df=data,model_formula_str = model_formula_str,model_opt="logistic",data_testing=NULL,save_fig=T,out_dir=".",out_suffix="")

test <- run_model_fun(data_df=list_data_training,
                      model_formula_str = model_formula_str,
                      model_opt="logistic",
                      data_testing=list_data_testing,
                      num_cores=1,
                      out_dir=".",
                      out_suffix="")

test <- run_model_fun(data_df=list_data_training,
                      model_formula_str = model_formula_str,
                      model_opt="randomForest",
                      data_testing=list_data_testing,
                      num_cores=1,
                      out_dir=".",
                      out_suffix="")
list_predicted_val <- test$predicted_val

#rep(data===)
names(list_data_testing) <- paste0("data_testing_",1:length(list_data_testing)) 
#debug(ROC_evaluation_fun)

roc_obj <- ROC_evaluation_fun(i=1,
                              list_data=list_data_testing,
                              y_var=y_var,
                              predicted_val=list_predicted_val,
                              save_fig=T,
                              out_suffix=out_suffix,
                              out_dir=out_dir)

list_roc_obj <- lapply(1:length(list_predicted_val),
                       FUN=ROC_evaluation_fun,
                       list_data=list_data_testing,
                       y_var=y_var,
                       predicted_val=list_predicted_val,
                       save_fig=T,
                       out_suffix=out_suffix,
                       out_dir=out_dir)

list_roc_obj[[2]] #gives ROC table values and overall ROC AUC value
list_roc_obj[[1]]
list_roc_obj[[3]]
roc_obj[[1]] #gives only ROC AUC value for specific object indicated
list_data_testing[[1]] #gives the data and obs id used to test roc_obj 1
list_data_training[[1]] #gives the data and obs id used to train roc_obj 1
################################ END OF SCRIPT ###################