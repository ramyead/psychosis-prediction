library(studerus)
library(mlr)
library(parallelMap)
library(reshape)
library(ggplot2)
library(Cairo)
library(grid)
parallelStartSocket(4)
load('data/coeff_lagged_Phase_Synchronicity.Rda')
dat <- subset(dat_c_phase, neurolep_now == 'no')
dat1 <- dat[dat$EEG_TDat < '2012-01-01', grep('group2|gamma|beta', names(dat))]

task1 <- makeClassifTask('ARMS.T_vs_ARMS.NT', droplevels(subset(dat1, group2 %in% c('ARMS-T', 'ARMS-NT'))),
                         target = 'group2', positive = 'ARMS-T')
task1 = normalizeFeatures(task1)
task2 <- makeClassifTask('ARMS.T_vs_HC', droplevels(subset(dat1, group2 %in% c('ARMS-T', 'HC'))),
                         target = 'group2', positive = 'ARMS-T')
task3 <- makeClassifTask('ARMS.NT_vs_HC', droplevels(subset(dat1, group2 %in% c('ARMS-NT', 'HC'))),
                         target = 'group2', positive = 'HC')
tasks <- list(task1, task2, task3)
meas <- list(auc, bac, tpr, tnr)

lrn <- makeLearner("classif.LiblineaRLogReg", predict.type = "prob", type = 6)
outer <- makeResampleDesc("RepCV", folds = 10, reps = 20,  stratify = T)
ctrl <- makeTuneControlGrid(resolution = 10)
ps <- makeParamSet(makeNumericParam("cost",lower = 0.1, upper = 15))

coefs <- sapply(tasks, function(task){
  x <- as.vector(table(getTaskData(task, target.extra = T)$target))
  w <- max(x)/min(x)
  #print(w)
  weight.wrap <- makeWeightedClassesWrapper(lrn, wcw.param = 'wi', wcw.weight = w)
  out <- tuneParams(weight.wrap, task, outer, meas, ps, ctrl)
  opt.path <- data.frame(out$opt.path)
  opt.cost <- opt.path$cost[which.max(opt.path$auc.test.mean)]
  print(opt.cost)
  opt.lrn <- makeLearner("classif.LiblineaRLogReg", predict.type = "prob", type = 6, cost = opt.cost)
  opt.weight.wrap <- makeWeightedClassesWrapper(opt.lrn, wcw.param = 'wi', wcw.weight = w)
  train(opt.weight.wrap, task)$learner.model$next.model$learner.model$W
})

df <- data.frame(coefs[-nrow(coefs),])
names(df) <- gsub('\\.', '-', gsub('_', ' ', sapply(tasks, function(x) x$task.desc$id)))
df$freq = gsub('.*_', '', getTaskFeatureNames(task1))
df <- melt(df, measure.vars = 1:3)

ggplot(df, aes(fill = value, x = freq, y = freq)) + geom_tile(color = 'black')+
  scale_fill_gradient2(low = "blue", high = "red", mid = 'white',
                       guide = guide_colorbar(title = 'Lasso\nRegression\ncoefficient\n', 
                                              label.hjust = 1, barwidth = 1.2, barheight = 7)) +
  labs(x = 'Frequency', y = 'Region of Interest') +
  scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + 
  facet_grid(.~variable) + theme_bw()

coord <- read.table(textConnection("
-2.0;  6.0;  1
 2.0;  6.0;  2
-5.0;  3.5;  3
-2.5;  3.0;  4
 0.0;  3.0;  5
 2.5;  3.0;  6
 5.0;  3.5;  7
-6.0;  0.0;  8
-3.0;  0.0;  9
 0.0;  0.0; 10
 3.0;  0.0; 11
 6.0;  0.0; 12
-5.0; -3.5; 13
-2.5; -3.0; 14
 0.0; -3.0; 15
 2.5; -3.0; 16
 5.0; -3.5; 17
-2.0; -6.0; 18
 2.0; -6.0; 19
"), sep = ";", stringsAsFactors = F, strip.white = T, quote = "", col.names = c("x", "y", "ROI"))

df2 <- merge(df, coord)
df2$Sign <- ifelse(df2$value > 0, 'Positive', 'Negative')
df2$Sign[df2$value == 0] <- 'Neutral'


df2= subset(df2, subset =  variable == "ARMS-T vs ARMS-NT")
df2= df2[ order(-df2[,4], df2[,1]), ]
#df2$x= NULL
#df2$y = NULL



z =9
p <- ggplot(df2, aes(x = x, y = y)) +
  geom_point(aes(size = abs(value), color = Sign)) +
  scale_size(range = c(1, 10)) +
  facet_grid(freq~variable) +
  scale_color_manual(breaks = c('Positive', 'Neutral', 'Negative'),
                     values = c('blue', 'black', 'red')) + theme_bw()+
  labs(size = 'Size of\ncoefficient', x = NULL, y = NULL) +
  geom_point(x = 0, y = 0, size = 85, shape = 1, col = 'black')+
  geom_point(x = 0, y = 8.2, shape = 2, size = 6) +
  coord_fixed(xlim = c(-z, z), ylim = c(-z, z+1.5)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        plot.margin=unit(c(0,0,0,0),"lines"))
fn <- 'Figures/variable_importance.pdf'
ggsave(p, device = CairoPDF, filename = fn,
      width = 10, height = 10)
browseURL(fn)
  
  


