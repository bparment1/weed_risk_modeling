############### SESYNC Research Support: weed risk ########## 
## Performing ROC on data for model assessment.
## 
## DATE CREATED: 06/29/2017
## DATE MODIFIED: 06/29/2017
## AUTHORS: Benoit Parmentier 
## PROJECT: weed risk Chinchu Harris
## ISSUE: 
## TO DO:
##
## COMMIT:initial function script with model run function
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

run_model_fun <- function(data_df,model_formula_str,model_opt,data_testing=NULL,save_fig=F,out_dir=".",out_suffix=""){
  #data_df: input data.frame with data used in modeling
  #model_formula_str
  #model_opt: "logistic","randomForest"
  #out_dir
  #out_suffix
  
  if(model_opt=="logistic"){
    mod <- glm(model_formula_str,data = data_df) #this is the training data!!
  }
  
  ### add randomForest option
  
  mod$fitted.values #these are the probability values from ROC
  data_df[,y_var]
  table(data_df[,y_var])
  
  
  if(save_fig==TRUE){
    
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
    rocd2 <- ROC(index=mod$fitted.values, boolean=data[[y_var]], mask=mask_val, nthres = 100)
    
    slot(rocd2,"AUC") #this is your AUC from the logistic modeling
    #Plot ROC curve:
    plot(rocd2)
    
    dev.off()
  }
  
  return(png_file_name)
  
}


############################### END OF SCRIPT ################################