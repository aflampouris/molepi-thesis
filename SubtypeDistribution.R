library(Hmisc)
library(rworldxtra)
library(rworldmap)
library(RColorBrewer)
library(plyr)
library(prevalence)
library(reporttools)

setwd("/home/andreas/MSc-Biostatistics/Thesis/Manuscript_new/figures/")

dfA<-read.csv2("~/MSc-Biostatistics/Thesis/Data/HIV-1_lanl.csv", sep = "\t", row.names=NULL)
dfA<-dfA[, c(6, 7, 8, 10)]
dfA<-dfA[(!is.na(dfA$Country) | !is.na(dfA$Georegion )) & !is.na(dfA$Subtype), ]
dfA<-dfA[(dfA$Country!='' | dfA$Georegion!='') & dfA$Subtype!='', ]

paletteCols<-c("#6CBB1D","#075E5E","#5E0707","#1DBBBB","#BB1D1D","#325E07","#6C1DBB","#666666")

dfA<-dfA[dfA[["Sampling.Year"]]  >1995,]
########################## WORLD #################################

mostPrev<-count(dfA, "Subtype")
mostPrev<-head(mostPrev[order(mostPrev$freq, decreasing=T), 1], 7)
levels(dfA[["Subtype"]])[!(levels(dfA[["Subtype"]]) %in% mostPrev)]<-"Other"

# tableNominal(vars = dfA[,1:2], cap ="Table of nominal variables.", lab = "tab: nominal",cumsum = FALSE,caption.placement = "top", longtable=F, print.pval = "fisher")   
# #     
# ci<-propCI(count(dfA$Subtype)[,2], dim(dfA)[1], method = 'cp', level = 0.95, sortby = "level")
# rownames(ci)<-count(dfA$Subtype)[,1]
# ci[,c(3,6:7)]<-sapply(ci[,c(3,6:7)],round,3)
# ci[,c(3,6:7)]*100


dfA<-count(dfA, c("Country", "Subtype", "Georegion"))

df3<-reshape(data = dfA[, c(1, 2, 4)], direction = "wide", timevar = "Subtype", idvar="Country")
df3[is.na(df3)]<-0
df3<-df3[c(1, 5, 3:4,7, 6,8,9, 2)]


dfMap<-joinCountryData2Map(dF = df3, joinCode = "NAME", "Country", nameCountryColumn = "Country", suggestForFailedCodes = T, verbose=T)
names(dfMap)<-gsub("freq.", "", names(dfMap))
dfMap$NAME<-gsub("N. Cyprus",NA,dfMap$NAME)
dfMap$NAME<-gsub("Macedonia","F.Y.R.O.M.",dfMap$NAME)
cairo_pdf("WorldSubB.pdf", height = 12, width = 26,family="Linux Biolinum" ,pointsize=20)
mapPies( dfMap, nameX="LON", nameY="LAT", maxZVal = 100, zColours = paletteCols, barOrient = "vert", landCol="gray90", 
 nameZs=names(dfMap)[c(51:58)], barWidth = 2, 
 symbolSize = 2, 
#  mapRegion='world', 
 ylim=c(-55, 75), xlim=c(-120,130),
 addSizeLegend = T, barRelative = F)
# labelCountries( col = "#2F2733", cex = 0.4)
dev.off()

########################## EUROPE #################################
dfA<-read.csv2("~/MSc-Biostatistics/Thesis/Data/HIV-1_lanl.csv", sep = "\t", row.names=NULL)
dfA<-dfA[, c(6, 7, 8, 10)]
dfA<-dfA[dfA[["Sampling.Year"]]>1995,]

dfA<-dfA[(!is.na(dfA$Country) | !is.na(dfA$Georegion )) & !is.na(dfA$Subtype), ]
dfA<-dfA[(dfA$Country!='' | dfA$Georegion!='') & dfA$Subtype!='', ]

dfA<-dfA[dfA[["Georegion"]]=='Europe', ]

mostPrev<-count(dfA, "Subtype")
mostPrev<-head(mostPrev[order(mostPrev$freq, decreasing=T), 1], 7)
levels(dfA[["Subtype"]])[!(levels(dfA[["Subtype"]]) %in% mostPrev)]<-"Other"

dfA<-count(dfA, c("Country", "Subtype", "Georegion"))

df3<-reshape(data = dfA[, c(1, 2, 4)], direction = "wide", timevar = "Subtype", idvar="Country")
df3[is.na(df3)]<-0
df3<-df3[c(1, 5, 3:4, 6:9, 2)]

dfMap<-joinCountryData2Map(dF = df3, joinCode = "NAME", "Country", nameCountryColumn = "Country", suggestForFailedCodes = T, verbose=T)
names(dfMap)<-gsub("freq.", "", names(dfMap))
dfMap$NAME<-gsub("N. Cyprus",NA,dfMap$NAME)
dfMap$NAME<-gsub("Macedonia","F.Y.R.O.M.",dfMap$NAME)
cairo_pdf("EuroSubB.pdf", height = 10, width = 16,family="Linux Biolinum" ,pointsize=20)
mapPies( dfMap, nameX="LON", nameY="LAT", zColours = paletteCols, 
 nameZs=names(dfMap)[c(51:58)], barWidth = 1, landCol="grey90",
 ylim=c(50,57), xlim=c(-49,35), maxZVal=1, symbolSize=2,
#  mapRegion='europe', 
 addSizeLegend = T, barRelative=F,lwdSymbols=0.1)
 labelCountries( col = "#2F2733", cex = 0.3)
dev.off()



########################## EUROPE HET#################################
dfHET<-read.csv2("~/MSc-Biostatistics/Thesis/Data/HET_EUR_subtypes.csv", sep = "\t", row.names=NULL)
dfHET<-dfHET[, c(6, 7, 8, 9, 10)]
dfHET<-dfHET[dfHET[["Sampling.Year"]]>1995,]
dfHET<-dfHET[(!is.na(dfHET$Country) | !is.na(dfHET$Georegion )) & !is.na(dfHET$Subtype), ]
dfHET<-dfHET[(dfHET$Country!='' | dfHET$Georegion!='') & dfHET$Subtype!='', ]

dfHET<-dfHET[dfHET[["Georegion"]]=='Europe', ]

mostPrev<-count(dfHET, "Subtype")
mostPrev<-head(mostPrev[order(mostPrev$freq, decreasing=T), 1], 7)
levels(dfHET[["Subtype"]])[!(levels(dfHET[["Subtype"]]) %in% mostPrev)]<-"Other"

# 
# tableNominal(vars = dfHET[,1:2], cap ="Table of nominal variables.", lab = "tab: nominal",cumsum = FALSE,caption.placement = "top", longtable=F, print.pval = "fisher")   
#     
    
dfHET<-count(dfHET, c("Country", "Subtype", "Georegion"))

df3<-reshape(data = dfHET[, c(1, 2, 4)], direction = "wide", timevar = "Subtype", idvar="Country")
df3[is.na(df3)]<-0
df3<-df3[c(1, 7, 2:5, 8:9, 6)]


dfMap<-joinCountryData2Map(dF = df3, joinCode = "NAME", "Country", nameCountryColumn = "Country", suggestForFailedCodes = T, verbose=T)

names(dfMap)<-gsub("freq.", "", names(dfMap))
dfMap$NAME<-gsub("N. Cyprus",NA,dfMap$NAME)
dfMap$NAME<-gsub("Macedonia","F.Y.R.O.M.",dfMap$NAME)

cairo_pdf("EuroHET.pdf", height = 10, width = 16,family="Linux Biolinum" ,pointsize=20)
mapPies( dfMap, nameX="LON", nameY="LAT", zColours = paletteCols, 
 nameZs=names(dfMap)[c(51:58)], barWidth = 1, landCol="grey90",
 ylim=c(50,57), xlim=c(-49,35), maxZVal=1, symbolSize=2,
#  mapRegion='europe', 
 addSizeLegend = T, barRelative=F,lwdSymbols=0.1)
 labelCountries( col = "#2F2733", cex = 0.3)
dev.off()

########################## EUROPE MSM#################################
dfMSM<-read.csv2("~/MSc-Biostatistics/Thesis/Data/MSM_EUR_subtypes.csv", sep = "\t", row.names=NULL)
dfMSM<-dfMSM[, c(6, 7, 8, 9, 10)]
dfMSM<-dfMSM[dfMSM[["Sampling.Year"]]>1995,]

dfMSM<-dfMSM[(!is.na(dfMSM$Country) | !is.na(dfMSM$Georegion )) & !is.na(dfMSM$Subtype), ]
dfMSM<-dfMSM[(dfMSM$Country!='' | dfMSM$Georegion!='') & dfMSM$Subtype!='', ]

dfMSM<-dfMSM[dfMSM[["Georegion"]]=='Europe', ]

mostPrev<-count(dfMSM, "Subtype")
mostPrev<-head(mostPrev[order(mostPrev$freq, decreasing=T), 1], 7)
levels(dfMSM[["Subtype"]])[!(levels(dfMSM[["Subtype"]]) %in% mostPrev)]<-"Other"

# 
# tableNominal(vars = dfMSM[,1:2], cap ="Table of nominal variables.", lab = "tab: nominal",cumsum = FALSE,caption.placement = "top", longtable=F, print.pval = "fisher")   
#     

dfMSM<-count(dfMSM, c("Country", "Subtype", "Georegion"))

df3<-reshape(data = dfMSM[, c(1, 2, 4)], direction = "wide", timevar = "Subtype", idvar="Country")
df3[is.na(df3)]<-0
df3<-df3[c(1:2, 4:9, 3)]


dfMap<-joinCountryData2Map(dF = df3, joinCode = "NAME", "Country", nameCountryColumn = "Country", suggestForFailedCodes = T, verbose=T)

names(dfMap)<-gsub("freq.", "", names(dfMap))
dfMap$NAME<-gsub("N. Cyprus",NA,dfMap$NAME)
dfMap$NAME<-gsub("Macedonia","F.Y.R.O.M.",dfMap$NAME)
cairo_pdf("EuroMSM.pdf", height = 10, width = 16,family="Linux Biolinum" ,pointsize=20)
mapPies( dfMap, nameX="LON", nameY="LAT", zColours = paletteCols, 
 nameZs=names(dfMap)[c(51:58)], barWidth = 1, landCol="grey90",
 ylim=c(50,57), xlim=c(-49,35), maxZVal=1, symbolSize=2,
#  mapRegion='europe', 
 addSizeLegend = T, barRelative=F,lwdSymbols=0.1)
 labelCountries( col = "#2F2733", cex = 0.3)
dev.off()


########################## EUROPE IDU#################################
dfIDU<-read.csv2("~/MSc-Biostatistics/Thesis/Data/IDU_EUR_subtypes.csv", sep = "\t", row.names=NULL)
dfIDU<-dfIDU[, c(6, 7, 8, 9, 10)]
dfIDU<-dfIDU[dfIDU[["Sampling.Year"]]>1995,]

dfIDU<-dfIDU[(!is.na(dfIDU$Country) | !is.na(dfIDU$Georegion )) & !is.na(dfIDU$Subtype), ]
dfIDU<-dfIDU[(dfIDU$Country!='' | dfIDU$Georegion!='') & dfIDU$Subtype!='', ]

dfIDU<-dfIDU[dfIDU[["Georegion"]]=='Europe', ]

mostPrev<-count(dfIDU, "Subtype")
mostPrev<-head(mostPrev[order(mostPrev$freq, decreasing=T), 1], 7)
levels(dfIDU[["Subtype"]])[!(levels(dfIDU[["Subtype"]]) %in% mostPrev)]<-"Other"

tableNominal(vars = dfIDU[,1:2], cap ="Table of nominal variables.", lab = "tab: nominal",cumsum = FALSE,caption.placement = "top", longtable=F, print.pval = "fisher")   


dfIDU<-count(dfIDU, c("Country", "Subtype", "Georegion"))

df3<-reshape(data = dfIDU[, c(1, 2, 4)], direction = "wide", timevar = "Subtype", idvar="Country")
df3[is.na(df3)]<-0
df3<-df3[c(1,2,5,7,4,6,8,9,3)]


dfMap<-joinCountryData2Map(dF = df3, joinCode = "NAME", "Country", nameCountryColumn = "Country", suggestForFailedCodes = T, verbose=T)
names(dfMap)<-gsub("freq.", "", names(dfMap))
dfMap$NAME<-gsub("N. Cyprus",NA,dfMap$NAME)
dfMap$NAME<-gsub("Macedonia","F.Y.R.O.M.",dfMap$NAME)
cairo_pdf("EuroIDU.pdf", height = 10, width = 16,family="Linux Biolinum" ,pointsize=20)
mapPies( dfMap, nameX="LON", nameY="LAT", zColours = paletteCols, 
 nameZs=names(dfMap)[c(51:58)], barWidth = 1, landCol="grey90",
 ylim=c(50,57), xlim=c(-49,35), maxZVal=1, symbolSize=2,
#  mapRegion='europe', 
 addSizeLegend = T, barRelative=F,lwdSymbols=0.1)
 labelCountries( col = "#2F2733", cex = 0.3)
dev.off()
