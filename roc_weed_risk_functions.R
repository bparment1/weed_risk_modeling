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

run_model_fun <- function(data_df,model_formula_str,model_opt,data_testing=NULL,num_cores=1,out_dir=".",out_suffix=""){
  #data_df: input data.frame with data used in modeling
  #model_formula_str
  #model_opt: "logistic","randomForest"
  #out_dir
  #out_suffix
  
  ### Make a list if data is not a list:
  if(class(data_df)!="list"){
    data_df <- list(data_df)
  }
  
  if(!is.null(data_testing)){
    if(class(data_testing)!="list"){
      data_testing <- list(data_testing)
    }
  }

  if(model_opt=="logistic"){
    
    #
    #if(class(data_df)!="list"){
    #  list_mod <- glm(model_formula_str,data = data_df,family=binomial()) #this is the training data!!
    #}else{
      
    list_mod <- mclapply(1:length(data_df),
                           FUN=function(i,model_formula_str,data_df){data_input<-data_df[[i]]; 
                           mod_glm <-  glm(model_formula_str, data = data_input);
                           return(mod_glm)},
                           model_formula_str=model_formula_str,
                           data_df=data_df,
                           mc.preschedule = FALSE,
                           mc.cores =num_cores)
      #mclapply used to return a list that is the same length as the x component (in this case x=1:length)
      #mc.preschedule=False because there aren't large numbers of x values 
    ### Apprend prediction to training data.frame!!!
    
    if(!is.null(data_testing)){
      #
      list_predicted_val <- mclapply(1:length(data_df),
                                     FUN=function(i,list_mod,data_testing){data_v<-data_testing[[i]]; 
                                     mod_glm <- list_mod[[i]];
                                     predicted_val <- predict(mod_glm,newdata=data_v,type='response');
                                     return(predicted_val)},
                                     list_mod=list_mod,
                                     data_testing=data_testing,
                                     mc.preschedule = FALSE,
                                     mc.cores =num_cores)
      
    }else{
      list_predicted_val <- NULL
    }
      
  }
  
  ##if(model=="RF)##
  if(model_opt=="randomForest"){
    list_mod <- mclapply(1:length(data_df),
                         FUN=function(i, model_formula_str, data_df){data_input<-data_df[[i]]; 
                         mod_rf <-randomForest(as.formula(model_formula_str), type = "classification",
                                            importance=TRUE,proximity=TRUE,ntree=1000,data= data_input);
                         return(mod_rf)},
                         model_formula_str=as.formula(model_formula_str),
                         data_df=data_df,
                         mc.preschedule = FALSE,
                         mc.cores =num_cores)
    
    if(!is.null(data_testing)){
      #
      list_predicted_val <- mclapply(1:length(data_df),
                                     FUN=function(i,list_mod,data_testing){data_v<-data_testing[[i]]; 
                                     mod_rf <- list_mod[[i]];
                                     predicted_val <- predict(mod_rf,newdata=data_v,type='response');
                                     return(predicted_val)},
                                     list_mod=list_mod,
                                     data_testing=data_testing,
                                     mc.preschedule = FALSE,
                                     mc.cores =num_cores)
    } else{
      list_predicted_val <- NULL
    }
    } 
    
  ################*bold bold *#####################
  model_obj <- list(list_mod,list_predicted_val,data_df,data_testing)
  names(model_obj) <- c("mod","predicted_val","data_training","data_testing")
  return(model_obj)
} 
  
  
ROC_evaluation_fun <- function(i,list_data,y_var,predicted_val,save_fig=T,out_suffix="",out_dir="."){

  ### Now do accuracy assessment:

  data_df <- list_data[[i]]
  dataset_name <- names(list_data)[i]
  #mod$fitted.values #these are the probability values from ROC
  #data_df[,y_var]
  #table(data_df[,y_var])

  mask_val <- 1:length(predicted_val[[i]])
  rocd2 <- ROC(index=predicted_val[[i]], 
               boolean=data_df[[y_var]], 
               mask=mask_val, 
               nthres = 100)
  
  AUC_val <- slot(rocd2,"AUC") #this is your AUC from the logistic modeling

  ROC_table <- slot(rocd2,"table")
  
  
  #Plot ROC curve:
  
  if(save_fig==TRUE){
    
    out_suffix_str <- paste0(dataset_name,out_suffix)
    
    res_pix<-480 #set as function argument...
    col_mfrow<-1
    #row_mfrow<-2
    row_mfrow<-1
    
    png_filename<- paste("Figure_","ROC_plot_",out_suffix_str,".png", sep="")
    
    png(filename=file.path(out_dir,png_filename),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    par(mfrow=c(row_mfrow,col_mfrow))
    
    plot(rocd2) #add option later...
    dev.off()
  }
  
  ###### prepare object ####
  roc_obj <- list(AUC_val,ROC_table,png_filename)
  names(roc_obj) <-c("AUC_val","ROC_table","png_filename")
  return(roc_obj)
  
}


############################### END OF SCRIPT ################################