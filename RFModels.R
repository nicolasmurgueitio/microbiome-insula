#RF Models

setwd()

library("randomForest")
library("plyr") 
library("rfUtilities") 
library("caret") 
library("mosaic")
library("sjPlot")
library(ggplot2)
library(jtools)
library(regclass)
library(jtools)

#Genus level bacteria

otu_table_brain <- read.table("", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")  
metadata_brain <- read.table("", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}
otu_table_brain <- remove_rare(table=otu_table_brain, cutoff_pro=0.50)


# Amygdala
otu_table_scaled_amygdala <- data.frame(t(otu_table_brain))  
otu_table_scaled_amygdala$amygdala <- metadata_brain[rownames(otu_table_scaled_amygdala), "amygdala"]  

set.seed(123)
RF_amygdala <- randomForest( x=otu_table_scaled_amygdala[,1:(ncol(otu_table_scaled_amygdala)-1)] , 
                             y=otu_table_scaled_amygdala$amygdala , ntree=501, importance=TRUE, proximities=TRUE )  
RF_amygdala


RF_amygdala_sig <- rf.significance( x=RF_amygdala ,  xdata=otu_table_scaled_amygdala[,1:(ncol(otu_table_scaled_amygdala)-1)] , nperm=1000 , ntree=501 )  
RF_amygdala_sig 

varImpPlot_thal<-varImpPlot(RF_amygdala)
importance(RF_thal)

RF_amygdala_imp <- as.data.frame(RF_amygdala$importance )
RF_amygdala_imp$features <- rownames( RF_amygdala_imp )
RF_amygdala_imp_sorted <- arrange( RF_amygdala_imp  , desc(`%IncMSE`)  )
RF_amygdala_imp_sorted 


# Insula
otu_table_scaled_insula <- data.frame(t(otu_table_brain))  
otu_table_scaled_insula$insula <- metadata_brain[rownames(otu_table_scaled_insula), "insula"]  

set.seed(123)
RF_insula <- randomForest( x=otu_table_scaled_insula[,1:(ncol(otu_table_scaled_insula)-1)] , 
                           y=otu_table_scaled_insula$insula , ntree=501, importance=TRUE, proximities=TRUE )  
RF_insula

RF_insula_sig <- rf.significance( x=RF_insula ,  xdata=otu_table_scaled_insula[,1:(ncol(otu_table_scaled_insula)-1)] , nperm=1000 , ntree=501 )  
RF_insula_sig 

varImpPlot_thal<-varImpPlot(RF_insula)
importance(RF_thal)

RF_insula_imp <- as.data.frame(RF_insula$importance )
RF_insula_imp$features <- rownames( RF_insula_imp )
RF_insula_imp_sorted <- arrange( RF_insula_imp  , desc(`%IncMSE`)  )
RF_insula_imp_sorted 

# Thalamus
otu_table_scaled_thalamus <- data.frame(t(otu_table_brain))  
otu_table_scaled_thalamus$thalamus <- metadata_brain[rownames(otu_table_scaled_thalamus), "thalamus"]  

set.seed(123)
RF_thalamus <- randomForest( x=otu_table_scaled_thalamus[,1:(ncol(otu_table_scaled_thalamus)-1)] , 
                             y=otu_table_scaled_thalamus$thalamus , ntree=501, importance=TRUE, proximities=TRUE )  
RF_thalamus


RF_thalamus_sig <- rf.significance( x=RF_thalamus ,  xdata=otu_table_scaled_thalamus[,1:(ncol(otu_table_scaled_thalamus)-1)] , nperm=1000 , ntree=501 )  
RF_thalamus_sig 

varImpPlot_thal<-varImpPlot(RF_thalamus)
importance(RF_thal)

RF_thalamus_imp <- as.data.frame(RF_thalamus$importance )
RF_thalamus_imp$features <- rownames( RF_thalamus_imp )
RF_thalamus_imp_sorted <- arrange( RF_thalamus_imp  , desc(`%IncMSE`)  )
RF_thalamus_imp_sorted 


# Anterior Cingulate
otu_table_scaled_acc <- data.frame(t(otu_table_brain))  
otu_table_scaled_acc$acc <- metadata_brain[rownames(otu_table_scaled_acc), "acc"]  

set.seed(123)
RF_acc <- randomForest( x=otu_table_scaled_acc[,1:(ncol(otu_table_scaled_acc)-1)] , 
                        y=otu_table_scaled_acc$acc , ntree=501, importance=TRUE, proximities=TRUE )  
RF_acc

RF_acc_sig <- rf.significance( x=RF_acc ,  xdata=otu_table_scaled_acc[,1:(ncol(otu_table_scaled_acc)-1)] , nperm=1000 , ntree=501 )  
RF_acc_sig 

#Gene pathways

otu_table_brain <- read.csv("",  header=T, row.names=1)  
metadata_brain <- read.table("", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

otu_table_brain <- data.frame(t(otu_table_brain))  


remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}
otu_table_brain <- remove_rare(table=otu_table_brain, cutoff_pro=0.50)
otu_table_brain <- data.frame(t(otu_table_brain))  

#explore data
dim(otu_table_brain) 

#explore metadata
dim(metadata_brain)

# Amygdala
otu_table_scaled_amygdala <- data.frame(t(otu_table_brain))  
otu_table_scaled_amygdala$amygdala <- metadata_brain[rownames(otu_table_scaled_amygdala), "amygdala"]  

set.seed(123)
RF_amygdala <- randomForest( x=otu_table_scaled_amygdala[,1:(ncol(otu_table_scaled_amygdala)-1)] , 
                             y=otu_table_scaled_amygdala$amygdala , ntree=501, importance=TRUE, proximities=TRUE )  
RF_amygdala


RF_amygdala_sig <- rf.significance( x=RF_amygdala ,  xdata=otu_table_scaled_amygdala[,1:(ncol(otu_table_scaled_amygdala)-1)] , nperm=1000 , ntree=501 )  
RF_amygdala_sig 


# Insula
otu_table_scaled_insula <- data.frame(t(otu_table_brain))  
otu_table_scaled_insula$insula <- metadata_brain[rownames(otu_table_scaled_insula), "insula"]  

set.seed(123)
RF_insula <- randomForest( x=otu_table_scaled_insula[,1:(ncol(otu_table_scaled_insula)-1)] , 
                           y=otu_table_scaled_insula$insula , ntree=501, importance=TRUE, proximities=TRUE )  
RF_insula

RF_insula_sig <- rf.significance( x=RF_insula ,  xdata=otu_table_scaled_insula[,1:(ncol(otu_table_scaled_insula)-1)] , nperm=1000 , ntree=501 )  
RF_insula_sig 

varImpPlot_thal<-varImpPlot(RF_insula)
importance(RF_thal)

RF_insula_imp <- as.data.frame(RF_insula$importance )
RF_insula_imp$features <- rownames( RF_insula_imp )
RF_insula_imp_sorted <- arrange( RF_insula_imp  , desc(`%IncMSE`)  )
RF_insula_imp_sorted 

# Thalamus
otu_table_scaled_thalamus <- data.frame(t(otu_table_brain))  
otu_table_scaled_thalamus$thalamus <- metadata_brain[rownames(otu_table_scaled_thalamus), "thalamus"]  

set.seed(123)
RF_thalamus <- randomForest( x=otu_table_scaled_thalamus[,1:(ncol(otu_table_scaled_thalamus)-1)] , 
                             y=otu_table_scaled_thalamus$thalamus , ntree=501, importance=TRUE, proximities=TRUE )  
RF_thalamus


RF_thalamus_sig <- rf.significance( x=RF_thalamus ,  xdata=otu_table_scaled_thalamus[,1:(ncol(otu_table_scaled_thalamus)-1)] , nperm=1000 , ntree=501 )  
RF_thalamus_sig 


# Anterior Cingulate
otu_table_scaled_acc <- data.frame(t(otu_table_brain))  
otu_table_scaled_acc$acc <- metadata_brain[rownames(otu_table_scaled_acc), "acc"]  

set.seed(123)
RF_acc <- randomForest( x=otu_table_scaled_acc[,1:(ncol(otu_table_scaled_acc)-1)] , 
                        y=otu_table_scaled_acc$acc , ntree=501, importance=TRUE, proximities=TRUE )  
RF_acc

RF_thalamus_sig <- rf.significance( x=RF_thalamus ,  xdata=otu_table_scaled_thalamus[,1:(ncol(otu_table_scaled_thalamus)-1)] , nperm=1000 , ntree=501 )  
RF_thalamus_sig 