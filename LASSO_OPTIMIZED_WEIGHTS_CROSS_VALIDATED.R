##LASSO TESTING WITH IDEAL WEIGHT
library(studerus)
library(mlr)
library(parallelMap)
library(beepr)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)


parallelStartSocket(3)

load('data/dat_the_one_and_only.Rda')
dat <- subset(dat, neurolep_now == 'no')
dat1 <- dat[dat$EEG_TDat < '2012-01-01', grep('group2|gamma|beta', names(dat))]


task1 <- makeClassifTask('ARMS.T_vs_ARMS.NT', droplevels(subset(dat1, group2 %in% c('ARMS-T', 'ARMS-NT'))),
                         target = 'group2', positive = 'ARMS-T')
task1 = normalizeFeatures(task1)
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
ggplot(x, aes(x = FPR, y = TPR, group = set, color = set)) +
  geom_path() + theme_bw() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = 'grey', size = 0.5) +
  coord_fixed(ratio = 1) +
  labs(y = 'Sensitivity', x = '1 - specificity',
       title = 'Receiver Operating Characteristic (ROC) curve\n',
       color = 'Data set') +
  theme(plot.margin = unit(c(0,0,0,0), "lines"))
ggsave('Figures/ROC_lasso_csd.png', width = 5.7, height = 5.4)


csd_plot=x
type= as.data.frame(rep("csd", 40))
names(type)= "type"
csd_plot= cbind(csd_plot, type)

#Result: bac.test.mean=0.698,auc.test.mean=0.788,auc.train.mean=0.997,tpr.test.mean=0.635,tnr.test.mean=0.762,mmce.test.mean=0.28,tp.test.mean=1.14,
#fp.test.mean=0.82,tn.test.mean=2.68,fn.test.mean=0.66,ppv.test.mean= NaN,npv.test.mean=0.828

















# #High freqs + BPRS
# dat1 <- dat[dat$EEG_TDat < '2012-01-01', grep('group2|gamma|beta|BPRS_Negative_symptoms|BPRS_Positive.Symptoms|timediff_EEG_TDat_bpTZPDat', names(dat))]
# # very High freqs only
# dat1 <- dat[dat$EEG_TDat < '2012-01-01', grep('group2|gamma|beta2', names(dat))]
# #All freqs + BPRS
# dat1 <- dat[dat$EEG_TDat < '2014-01-01', grep('group2|alpha|theta|delta|gamma|beta|BPRS_Negative_symptoms|BPRS_Positive.Symptoms|timediff_EEG_TDat_bpTZPDat', names(dat))]
# #BPRS only
# dat1 <- dat[dat$EEG_TDat < '2014-01-01', grep('group2|BPRS_Negative_symptoms|BPRS_Positive.Symptoms', names(dat))]
# #Low freqs only
# dat1 <- dat[dat$EEG_TDat < '2014-01-01', grep('group2|alpha|theta|delta', names(dat))]
# #low freqs + BPRS
# dat1 <- dat[dat$EEG_TDat < '2014-01-01', grep('group2|alpha|theta|delta|BPRS_Negative_symptoms|BPRS_Positive.Symptoms', names(dat))]
#dat1 = dat1[complete.cases(dat1$BPRS_Negative_symptoms),]

task1 <- makeClassifTask('ARMS.T_vs_ARMS.NT', droplevels(subset(dat1, group2 %in% c('ARMS-T', 'ARMS-NT'))),
                         target = 'group2', positive = 'ARMS-T')
task2 <- makeClassifTask('ARMS.T_vs_HC', droplevels(subset(dat1, group2 %in% c('ARMS-T', 'HC'))),
                         target = 'group2', positive = 'ARMS-T')
task3 <- makeClassifTask('ARMS.NT_vs_HC', droplevels(subset(dat1, group2 %in% c('ARMS-NT', 'HC'))),
                         target = 'group2', positive = 'HC')






##Filter Selection

lrn2 <- makeLearner("classif.LiblineaRLogReg", id = 'Lasso', predict.type = "prob", type = 6)
over.wrap <- makeWeightedClassesWrapper(lrn2, wcw.param = 'wi', wcw.weight =  getWeight(task1))
filt.wrap <- makeFilterWrapper(over.wrap, fw.method = 'information.gain', fw.val = 0.75)
inner <- makeResampleDesc("RepCV", folds = 8, reps = 10,  stratify = T)
outer <- makeResampleDesc("RepCV", folds = 8, reps = 10,  stratify = T)
ctrl <- makeTuneControlGrid(resolution = 10)
ps <- makeParamSet(makeNumericParam("cost",lower = .01, upper = 15))
tune.wrap <- makeTuneWrapper(filt.wrap, inner, par.set = ps, control = ctrl, measures = ber)
out <- resample(tune.wrap, task1, outer, measures = meas, extract = getTuneResult)
out$aggr
beep()


# bac.test.mean auc.test.mean tpr.test.mean tnr.test.mean 
# 0.4420573     0.4979167     0.3723958     0.5117188 


## T VS HC
lrn2 <- makeLearner("classif.LiblineaRLogReg", id = 'Lasso', predict.type = "prob", type = 6)
weight.wrap2 <- makeWeightedClassesWrapper(lrn2, wcw.param = 'wi', wcw.weight =  getWeight(task2))
tune.wrap <- makeTuneWrapper(weight.wrap, inner, par.set = ps, control = ctrl, measures = ber)
out <- resample(tune.wrap, task2, outer, measures = meas, extract = getTuneResult)
beep()

## NT VS HC
lrn2 <- makeLearner("classif.LiblineaRLogReg", id = 'Lasso', predict.type = "prob", type = 6)
weight.wrap2 <- makeWeightedClassesWrapper(lrn2, wcw.param = 'wi', wcw.weight =   getWeight(task3))
tune.wrap <- makeTuneWrapper(weight.wrap, inner, par.set = ps, control = ctrl, measures = ber)
system.time(out <- resample(tune.wrap, task3, outer, measures = meas, extract = getTuneResult))
beep()

parallelStop()


ggplot(data= out$measures.test, aes(x=iter, y= bac , group=1)) + geom_line() + geom_point()
