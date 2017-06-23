############### SESYNC Research Support: weed risk ########## 
## Performing ROC on data for model assessment.
## 
## DATE CREATED: 06/15/2017
## DATE MODIFIED: 06/23/2017
## AUTHORS: Benoit Parmentier 
## PROJECT: weed risk Chinchu Harris
## ISSUE: 
## TO DO:
##
## COMMIT:ROC testing on publicly available data
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

#function_pca_eof <- "pca_eof_functions_06152017.R" #PARAM 1
#script_path <- "/nfs/bparmentier-data/Data/projects/neighborhood_regression_spatial_demography/scripts" #path to script #PARAM 
#source(file.path(script_path,function_pca_eof)) #source all functions used in this script 1.

############################################################################
####################  Parameters and argument set up ###########

in_dir <- "/nfs/bparmentier-data/Data/projects/modeling_weed_risk/data" #local bpy50 , param 1
out_dir <- "/nfs/bparmentier-data/Data/projects/modeling_weed_risk/outputs" #param 2

num_cores <- 2 #param 8
create_out_dir_param=TRUE # param 9

out_suffix <-"roc_experiment_06232017" #output suffix for the files and ouptut folder #param 12

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

options(scipen=999)


### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

data <- read.csv(file.path(in_dir,infile_data))
dim(data)
View(data)
#> dim(data)
#[1] 94 42

###### SET UP EQUATION AND MODEL PREDICIONS #####

#NON=non invader=0 #OTH=major invader + minor invader=1
#67 instances of 0 #136 instances of 1
#data=read.csv(file="firstmodelwithnewformatting.csv")
head(data)
y_var <- "invasion.status"

explanatory_variables <- names(data)[-1]

right_side_formula <- paste(explanatory_variables,collapse = " + ")
model_formula_str <- paste0(y_var," ~ ",right_side_formula,sep="")

mod <- glm(model_formula_str,data = data)

mod$fitted.values
data[,y_var]
table(data[,y_var])

mask_val <- 1:nrow(data)
rocd2 <- ROC(index=mod$fitted.values, boolean=data[[y_var]], mask=mask_val, nthres = 100)

slot(rocd2,"AUC") #this is your AUC
#Plot ROC curve:
plot(rocd2)

##### Add code for training and testing

#This function creates testing and training list for input sation data based on a list of dates. 
#This function works for montly time scale if dates are provided as mid-months or other forms of for monthly records.
#It requires 6 inputs:                                           
# 1) seed_number: allow comparison across runs, if seed zero then no seed number is used
# 2) nb_sample: number of time random sampling must be repeated for every hold out proportion 
# 3) step : step for proportion range 
# 4) constant: if value 1 then use the same samples as date one for the all set of dates
# 5) prop_minmax: if prop_min=prop_max and step=0 then predicitons are done for the number of dates...
# 6) dates: list of dates for prediction and subsetting                              
# 7) ghcn: station data as data.frame -daily input                                                                                
# 11) out_prefix: output suffix added to output names--it is the same in the interpolation script
#
#The output is a list of four shapefile names produced by the function:
# 1) sampling_dat: sampling information for every run by date and sampling combintation
# 2) sampling_index: list of indexes for training and testing for every dates
# 3) sampling_stat_id: list of station ID for training and testing for every dates
# 4) ghcn_data: ghcn subsets by date, can be monthly or daily with mulitple sampling

seed_number <-list_param_sampling$seed_number
nb_sample <- list_param_sampling$nb_sample
step<-list_param_sampling$step #if seed zero then no seed?     
constant <- list_param_sampling$constant
prop_minmax<-list_param_sampling$prop_minmax
dates<-list_param_sampling$dates
#ghcn_name<-list_param_sampling$ghcn_name
ghcn<-list_param_sampling$ghcn #can be daily or monthly!!
#ghcn<-get(ghcn_name) 


################################ END OF SCRIPT ###################