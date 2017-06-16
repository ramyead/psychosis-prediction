library(studerus)
library(mlr)
library(parallelMap)
library(beepr)
library(ggplot2)
library(plyr)
library(dplyr)

#####Loading the values creating by the script below

load('data/roc_values.Rda')

ggplot(roc, aes(x = FPR, y = TPR, group = set, color = set)) + geom_line(linetype= 1)
ggplot(roc[40:80,], aes(x = FPR, y = TPR, group = set, color = set)) + geom_line(linetype= 3)

################################creating the values######################

parallelStartSocket(3)

load('data/dat_the_one_and_only.Rda')
dat <- subset(dat, neurolep_now == 'no')
dat1 <- dat[dat$EEG_TDat < '2012-01-01', grep('group2|gamma|beta', names(dat))]


task1 <- makeClassifTask('ARMS.T_vs_ARMS.NT', droplevels(subset(dat1, group2 %in% c('ARMS-T', 'ARMS-NT'))),
                         target = 'group2', positive = 'ARMS-T')
task2 <- makeClassifTask('ARMS.T_vs_HC', droplevels(subset(dat1, group2 %in% c('ARMS-T', 'HC'))),
                         target = 'group2', positive = 'ARMS-T')
task3 <- makeClassifTask('ARMS.NT_vs_HC', droplevels(subset(dat1, group2 %in% c('ARMS-NT', 'HC'))),
                         target = 'group2', positive = 'HC')

tasks <- list(task1, task2, task3)
auc.train <- setAggregation(auc, aggr = train.mean)

meas <- list(bac, mlr::auc, auc.train, tpr, tnr, mmce, tp, fp, tn, fn, ppv, npv )

inner <- makeResampleDesc("RepCV", folds = 10, reps = 10, stratify = T)
outer <- makeResampleDesc("RepCV", folds = 10, reps = 10,  stratify = T, predict = 'both')
ctrl <- makeTuneControlGrid(resolution = 10)
ps <- makeParamSet(makeNumericParam("cost",lower =.1, upper = 15))

getWeight <- function(task){
  x <- table(getTaskData(task, target.extra = T)$target)
  max(x)/min(x)
}


# T VS NT
lrn2 <- makeLearner("classif.LiblineaRLogReg", id = 'Lasso', predict.type = "prob", type = 6)
weight.wrap <- makeWeightedClassesWrapper(lrn2, wcw.param = 'wi', wcw.weight =  getWeight(task1))
tune.wrap <- makeTuneWrapper(weight.wrap, inner, par.set = ps, control = ctrl, measures = ber)
set.seed(seed = 5)
out <- resample(tune.wrap, task1, outer, measures = meas, extract = getTuneResult)
beep()

##[Resample] Result: bac.test.mean=0.701,auc.test.mean=0.775,auc.train.mean=0.997,tpr.test.mean=0.63,tnr.test.mean=0.772
#mmce.test.mean=0.273,tp.test.mean=1.14,fp.test.mean=0.78,tn.test.mean=2.72,fn.test.mean=0.66

## PLOC THE ROC CURVE

AUC_coord <- function(results, thresholds = 20){
  
  dat = out$pred$data[,c(1,2,4,6,7)]
  names(dat)[3] = 'prob'
  pos = out$pred$task.desc$positive
  thres = seq(0, 1, len = thresholds)
  x <- ddply(dat, ~ set + iter, function(df){
    data.frame(thres,
               TPR = sapply(thres, function(i) with(df, sum(prob >= i & truth == pos)/sum(truth == pos))),
               FPR = sapply(thres, function(i) with(df, sum(prob >= i & truth != pos)/sum(truth != pos))))
  })
  ddply(x, ~ set + thres, summarise, TPR = mean(TPR), FPR = mean(FPR))
}

x <- AUC_coord(out)
x <- AUC_coord(results = out)
ggplot(x, aes(x = FPR, y = TPR, group = set, color = set)) + geom_line()

csd_plot=x
type= as.data.frame(rep("csd", 40))
names(type)= "type"
csd_plot= cbind(csd_plot, type)

load('data/coeff_lagged_Phase_Synchronicity.Rda')

dat <- subset(dat_c_phase, neurolep_now == 'no')
dat <- dat[dat$EEG_TDat < '2012-01-01', grep('group2|gamma|beta', names(dat))]

task1 <- makeClassifTask('ARMS.T_vs_ARMS.NT', droplevels(subset(dat, group2 %in% c('ARMS-T', 'ARMS-NT'))),
                         target = 'group2', positive = 'ARMS-T')
task2 <- makeClassifTask('ARMS.T_vs_HC', droplevels(subset(dat, group2 %in% c('ARMS-T', 'HC'))),
                         target = 'group2', positive = 'ARMS-T')
task3 <- makeClassifTask('ARMS.NT_vs_HC', droplevels(subset(dat, group2 %in% c('ARMS-NT', 'HC'))),
                         target = 'group2', positive = 'HC')

tasks <- list(task1, task2, task3)
auc.train <- setAggregation(auc, aggr = train.mean)

meas <- list(bac, mlr::auc, auc.train, tpr, tnr, mmce, tp, fp, tn, fn, ppv, npv )

inner <- makeResampleDesc("RepCV", folds = 10, reps = 10, stratify = T)
outer <- makeResampleDesc("RepCV", folds = 10, reps = 10,  stratify = T, predict = 'both')
ctrl <- makeTuneControlGrid(resolution = 10)
ps <- makeParamSet(makeNumericParam("cost",lower = .1, upper = 15))

getWeight <- function(task){
  x <- table(getTaskData(task, target.extra = T)$target)
  max(x)/min(x)
}

## T VS NT
lrn2 <- makeLearner("classif.LiblineaRLogReg", id = 'Lasso', predict.type = "prob", type = 6)
weight.wrap <- makeWeightedClassesWrapper(lrn2, wcw.param = 'wi', wcw.weight =  getWeight(task1))
tune.wrap <- makeTuneWrapper(weight.wrap, inner, par.set = ps, control = ctrl, measures = ber)
set.seed(seed = 5)
out <- resample(tune.wrap, task1, outer, measures = meas, extract = getTuneResult)
beep()

# [Resample] Result: bac.test.mean=0.52,auc.test.mean= 0.5,auc.train.mean=  NA,tpr.test.mean=0.445,tnr.test.mean=0.595,mmce.test.mean=0.484,tp.test.mean=0.73,
##fp.test.mean=1.47,tn.test.mean=2.03,fn.test.mean=1.07,ppv.test.mean= NaN,npv.test.mean= NaN


AUC_coord <- function(results, thresholds = 20){
  
  dat = out$pred$data[,c(1,2,4,6,7)]
  names(dat)[3] = 'prob'
  pos = out$pred$task.desc$positive
  thres = seq(0, 1, len = thresholds)
  x <- ddply(dat, ~ set + iter, function(df){
    data.frame(thres,
               TPR = sapply(thres, function(i) with(df, sum(prob >= i & truth == pos)/sum(truth == pos))),
               FPR = sapply(thres, function(i) with(df, sum(prob >= i & truth != pos)/sum(truth != pos))))
  })
  ddply(x, ~ set + thres, summarise, TPR = mean(TPR), FPR = mean(FPR))
}

x <- AUC_coord(out)
x <- AUC_coord(results = out)
ggplot(x, aes(x = FPR, y = TPR, group = set, color = set)) + geom_line()

lagged_plot=x
type= as.data.frame(rep("lagged", 40))
names(type)= "type"
lagged_plot= cbind(lagged_plot, type)

roc= rbind(csd_plot, lagged_plot)

save(roc, file = 'data/roc_values.Rda')
