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

#Functional pathways

path <- read.csv("")
path$PWY.5130<-(path$PWY.5130+1)
path$PWY.5130<-as.numeric(log10(path$PWY.5130))
path$PWY.5130<-as.numeric(scale(path$PWY.5130,center=TRUE, scale=TRUE))

path$PWY.6823<-(path$PWY.6823+1)
path$PWY.6823<-as.numeric(log10(path$PWY.6823))
path$PWY.6823<-as.numeric(scale(path$PWY.6823,center=TRUE, scale=TRUE))

path$PWY.5005<-(path$PWY.5005+1)
path$PWY.5005<-as.numeric(log10(path$PWY.5005))
path$PWY.5005<-as.numeric(scale(path$PWY.5005,center=TRUE, scale=TRUE))

path$PWY66.409<-(path$PWY66.409+1)
path$PWY66.409<-as.numeric(log10(path$PWY66.409))
path$PWY66.409<-as.numeric(scale(path$PWY66.409,center=TRUE, scale=TRUE))

path$ARO.PWY<-(path$ARO.PWY+1)
path$ARO.PWY<-as.numeric(log10(path$ARO.PWY))
path$ARO.PWY<-as.numeric(scale(path$ARO.PWY,center=TRUE, scale=TRUE))

path$PWY.6902<-(path$PWY.6902+1)
path$PWY.6902<-as.numeric(log10(path$PWY.6902))
path$PWY.6902<-as.numeric(scale(path$PWY.6902,center=TRUE, scale=TRUE))

path$P42.PWY<-(path$P42.PWY+1)
path$P42.PWY<-as.numeric(log10(path$P42.PWY))
path$P42.PWY<-as.numeric(scale(path$P42.PWY,center=TRUE, scale=TRUE))

path$P164.PWY<-(path$P164.PWY+1)
path$P164.PWY<-as.numeric(log10(path$P164.PWY))
path$P164.PWY<-as.numeric(scale(path$P164.PWY,center=TRUE, scale=TRUE))

path$PWY.5686<-(path$PWY.5686+1)
path$PWY.5686<-as.numeric(log10(path$PWY.5686))
path$PWY.5686<-as.numeric(scale(path$PWY.5686,center=TRUE, scale=TRUE))

path$WY.6606<-(path$WY.6606+1)
path$WY.6606<-as.numeric(log10(path$WY.6606))
path$WY.6606<-as.numeric(scale(path$WY.6606,center=TRUE, scale=TRUE))

fit_insula_1<-(lm(Insula~PWY.5130+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
cooksThresh=4/dof
cdlist=cooks.distance(fit_insula_1)[abs(cooks.distance(fit_insula_1)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_1<-(lm(Insula~PWY.5130+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_1)

fit_insula_2<-(lm(Insula~PWY.6823+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
cdlist=cooks.distance(fit_insula_2)[abs(cooks.distance(fit_insula_2)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_2<-(lm(Insula~PWY.6823+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_2)


fit_insula_3<-(lm(Insula~PWY66.409+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
summary(fit_insula_3)
cdlist=cooks.distance(fit_insula_3)[abs(cooks.distance(fit_insula_3)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_3<-(lm(Insula~PWY66.409+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_3)

fit_insula_4<-(lm(Insula~ARO.PWY+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
summary(fit_insula_4)
cdlist=cooks.distance(fit_insula_4)[abs(cooks.distance(fit_insula_4)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_4<-(lm(Insula~ARO.PWY+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_4)

fit_insula_5<-(lm(Insula~PWY.6902+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
summary(fit_insula_5)
cdlist=cooks.distance(fit_insula_5)[abs(cooks.distance(fit_insula_5)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_5<-(lm(Insula~PWY.6902+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_5)

fit_insula_6<-(lm(Insula~P42.PWY+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
summary(fit_insula_6)
cdlist=cooks.distance(fit_insula_6)[abs(cooks.distance(fit_insula_6)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_6<-(lm(Insula~P42.PWY+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_6)


fit_insula_7<-(lm(Insula~P164.PWY+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
summary(fit_insula_7)
cdlist=cooks.distance(fit_insula_7)[abs(cooks.distance(fit_insula_7)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_7<-(lm(Insula~P164.PWY+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_7)


fit_insula_8<-(lm(Insula~PWY.5686+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
summary(fit_insula_8)
cdlist=cooks.distance(fit_insula_8)[abs(cooks.distance(fit_insula_8)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_8<-(lm(Insula~PWY.5686+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_8)

fit_insula_9<-(lm(Insula~WY.6606+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
summary(fit_insula_9)
cdlist=cooks.distance(fit_insula_9)[abs(cooks.distance(fit_insula_9)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_9<-(lm(Insula~WY.6606+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_9)


fit_insula_10<-(lm(Insula~PWY.5005+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=path))
summary(fit_insula_10)
cdlist=cooks.distance(fit_insula_10)[abs(cooks.distance(fit_insula_10)) > cooksThresh]
ind.remove = which(row.names(path) %in% names(cdlist))
outdata=path[-ind.remove,]
fit_insula_10<-(lm(Insula~PWY.5005+GADays_T2+GADays_T2+tc_sex+ICV+fp_brestfed_T2+PH_cesarean_T2+ITN_adj_T0,data=outdata))
summary(fit_insula_10)

models <- list(fit_insula_1, fit_insula_2,fit_insula_3,fit_insula_4,fit_insula_5,fit_insula_6,fit_insula_7,fit_insula_8,fit_insula_9,fit_insula_10)
result <- get_formatted_pvals_multiple_models(models)
print(result)
