setwd("")
data<-read.csv("")

library("lavaan")
library("interactions")
library("ggplot2")
library("stdmod")
library(ggiraphExtra)
library(dplyr)
library(regclass)
library(moments)

#Transform Data
data$Enterobacter<-(data$Enterobacter+1)
data$Enterobacter<-as.numeric(log10(data$Enterobacter))
data$Enterobacter<-as.numeric(scale(data$Enterobacter,center=TRUE, scale=TRUE))

data$Cutibacterium<-(data$Cutibacterium+1)
data$Cutibacterium<-as.numeric(log10(data$Cutibacterium))
data$Cutibacterium<-as.numeric(scale(data$Cutibacterium,center=TRUE, scale=TRUE))

data$Vibrio<-(data$Vibrio+1)
data$Vibrio<-as.numeric(log10(data$Vibrio))
data$Vibrio<-as.numeric(scale(data$Vibrio,center=TRUE, scale=TRUE))

data$Bifidobacterium<-(data$Bifidobacterium+1)
data$Bifidobacterium<-as.numeric(log10(data$Bifidobacterium))
data$Bifidobacterium<-as.numeric(scale(data$Bifidobacterium,center=TRUE, scale=TRUE))

data$Raoultella<-(data$Raoultella+1)
data$Raoultella<-as.numeric(log10(data$Raoultella))
data$Raoultella<-as.numeric(scale(data$Raoultella,center=TRUE, scale=TRUE))

data$Alistipes<-(data$Alistipes+1)
data$Alistipes<-as.numeric(log10(data$Alistipes))
data$Alistipes<-as.numeric(scale(data$Alistipes,center=TRUE, scale=TRUE))

data$Klebsiella<-(data$Klebsiella+1)
data$Klebsiella<-as.numeric(log10(data$Klebsiella))
data$Klebsiella<-as.numeric(scale(data$Klebsiella,center=TRUE, scale=TRUE))

data$Corynebacterium<-(data$Corynebacterium+1)
data$Corynebacterium<-as.numeric(log10(data$Corynebacterium))
data$Corynebacterium<-as.numeric(scale(data$Corynebacterium,center=TRUE, scale=TRUE))

data$Citrobacter<-(data$Citrobacter+1)
data$Citrobacter<-as.numeric(log10(data$Citrobacter))
data$Citrobacter<-as.numeric(scale(data$Citrobacter,center=TRUE, scale=TRUE))

data$Shigella<-(data$Shigella+1)
data$Shigella<-as.numeric(log10(data$Shigella))
data$Shigella<-as.numeric(scale(data$Shigella,center=TRUE, scale=TRUE))

data$Veillonella<-(data$Veillonella+1)
data$Veillonella<-as.numeric(log10(data$Veillonella))
data$Veillonella<-as.numeric(scale(data$Veillonella,center=TRUE, scale=TRUE))

data$Enterococcus<-(data$Enterococcus+1)
data$Enterococcus<-as.numeric(log10(data$Enterococcus))
data$Enterococcus<-as.numeric(scale(data$Enterococcus,center=TRUE, scale=TRUE))

#Amygdala and genus
fit<-(lm(amygdala~Alistipes+Klebsiella+Corynebacterium+Citrobacter+Shigella+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]

fit_amygdala<-(lm(amygdala~Alistipes+Klebsiella+Corynebacterium+Citrobacter+Shigella+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_amygdala)
VIF(fit_amygdala)

# Insula and Genus
fit<-(lm(insula~Veillonella+Enterococcus+Bifidobacterium+Vibrio+Corynebacterium+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]

fit_insula<-(lm(insula~Veillonella+Enterococcus+Bifidobacterium+Vibrio+Corynebacterium+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula)
VIF(fit_insula)

#Thalamus and genus
fit<-(lm(thalamus~Enterobacter+Cutibacterium+Vibrio+Bifidobacterium+Raoultella+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]
fit_thalamus<-(lm(thalamus~Enterobacter+Cutibacterium+Vibrio+Bifidobacterium+Raoultella+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_thalamus)
VIF(fit_thalamus)

#Corrected p-values and confidence intervals
models <- list(fit_amygdala, fit_insula,fit_thalamus)
result <- get_formatted_pvals_multiple_models(models)
print(result)

#Alpha diversity
fit<-(lm(thalamus~Shannon+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]
fit_thalamus<-(lm(thalamus~Shannon+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_thalamus)
VIF(fit_thalamus)

fit<-(lm(amygdala~Shannon+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]
fit_amygdala<-(lm(amygdala~Shannon+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_amygdala)
VIF(fit_amygdala)

fit<-(lm(insula~Shannon+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]
fit_insula<-(lm(insula~Shannon+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula)
VIF(fit_insula)

fit<-(lm(acc~Shannon+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]
fit_acc<-(lm(acc~Shannon+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_acc)
VIF(fit_acc)

#Corrected p-values and confidence intervals
models <- list(fit_amygdala, fit_insula,fit_thalamus,fit_acc)
result <- get_formatted_pvals_multiple_models(models)
print(result)


#Beta Diversity
fit<-(lm(thalamus~beta1+beta2+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]
fit_thalamus<-(lm(thalamus~beta1+beta2+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_thalamus)
VIF(fit_thalamus)

fit<-(lm(amygdala~beta1+beta2+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]
fit_amygdala<-(lm(amygdala~beta1+beta2+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_amygdala)
VIF(fit_amygdala)

fit<-(lm(insula~beta1+beta2+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]
fit_insula<-(lm(insula~beta1+beta2+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula)
VIF(fit_insula)

fit<-(lm(acc~beta1+beta2+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=data))
dof=df.residual(fit)
cooksThresh=4/dof
cdlist=cooks.distance(fit)[abs(cooks.distance(fit)) > cooksThresh]
ind.remove = which(row.names(data) %in% names(cdlist))
outdata=data[-ind.remove,]
fit_acc<-(lm(acc~beta1+beta2+GADays_T2+GADays_T2+tc_sex_T2+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_acc)
VIF(fit_acc)

#Corrected p-values and confidence intervals
models <- list(fit_amygdala, fit_insula,fit_thalamus,fit_acc)
result <- get_formatted_pvals_multiple_models(models)
print(result)
