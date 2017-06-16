library(studerus)
library(compareGroups)

# data <- getBD(screening < '2012-01-01' & group %in% 'ARMS')
# load('L:/APA_fepsy/Datenbanken/EEG_daten/Misc/np.dat.Rda')
# med= dat[,c(1,156)]
# nomed <- getBD(screening < '2012-01-01' & group %in% 'ARMS')
# data= merge(med, nomed, by= "pcode")
# data = data[data$neurolep_ever == "no",]
# data$included <- ifelse(data$pcode %in% dat$pcode, 'included', 'not included')
# eeg.dat = getEEG(nphTZP == 0, c(pcode, EEG_TDat))
# eeg.dat = eeg.dat[!duplicated(eeg.dat$pcode),]
# data <- merge(data, eeg.dat, all.x = T)
# i = is.na(data$EEG_TDat)
# data[i,'EEG_TDat'] <- data[i, 'screening']
# data = mergeNearestDate(data, getBPRS(), date.x = 'EEG_TDat', date.y = 'bpTZPDat', maxdiff = 60)
# data = mergeNearestDate(data, getSANS(), date.x = 'EEG_TDat', date.y = 'SAN_TDat', maxdiff = 60)


load(file = "data/compare_included_notincluded.rda")
x <- compareGroups(included ~ sex + age + years.edu + BPRS_total +
                     BPRS_Positive_Symptoms_Yung + BPRS_Negative_symptoms, data)

createTable(x)

