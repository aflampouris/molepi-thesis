
    list.of.packages <- c("RMySQL" ,"ggplot2" ,"Hmisc" ,"papeR" ,"scales" ,"reporttools","rworldmap","plyr")
         lapply(list.of.packages,require,character.only = T)
setwd("/home/andreas/MSc-Biostatistics/Thesis/Manuscript/figures/")




thesisdb.con <- dbConnect(MySQL(), user="R", password="1234", dbname="thesisdb", host="localhost")
dbData<- dbGetQuery(thesisdb.con, "SELECT * FROM DatasetFull WHERE region='Europe';")
dbDisconnect(thesisdb.con)

dbData <- dbData[,-c(1,2,10)]

################## percent from europe

thesisdb.con <- dbConnect(MySQL(), user="R", password="1234", dbname="thesisdb", host="localhost")
dbDataAll<- dbGetQuery(thesisdb.con, "SELECT * FROM Main ;")
dbDisconnect(thesisdb.con)

dbDataAll <- dbDataAll[,-c(1,2,10)]


x<-as.data.frame(table(dbData$submissionCountry))
p<-round(prop.table(x[,2])*100,2)
cbind(x,p)

library(countrycode)
euISO<-read.csv("~/whoEU.csv")
euISO<-countrycode(euISO[,1],"country.name" , "iso3c", warn = T)

from_eu<-x[["Var1"]] %in% euISO
from_eu_all<-dbDataAll[["submissionCountry"]] %in% euISO


table(dbData[dbDataAll[["submissionCountry"]] %in% euISO,7],useNA = "ifany")

n_eu_all<-table(from_eu_all)

eu_cl<-as.data.frame(table(dbData[(dbData[["submissionCountry"]] %in% euISO),1]))
all_cl<-as.data.frame(table(dbData[,1]))
monophyly<-round(eu_cl[,2]/all_cl[,2],4)*100
mono_table<-cbind(all_cl,eu_cl[,2],monophyly)

# latex(mono_table,rowname = NULL, label="EURclear",size = "normalsize",where="htpb", extracolsize='scriptsize'
# ,caption = paste("Statistical test results for monophyletic cluster No.")
# ,na.blank=F,col.just = c("c","c","c","c")
# ,file = "~/MSc-Biostatistics/Thesis/Manuscript_new/tables/monophyly.tex",ctable=F, booktabs=T
# ,dcolumn = T,colheads = c("Monophyletic","Total","European","monophyly"),extracolheads = c("cluster No.","Seqs","Seqs","\\%"))


# table(dbData[!(dbData[["submissionCountry"]] %in% eu[,2]),3])
# 
# summary(dbData[["testYear"]])


#################RISK GROUP DISTRIBUTION

x<-as.data.frame(table(dbData$submissionCountry))
p<-round(prop.table(x[,2])*100,2)
cbind(x,p)


##########Tabulation

library(Rz)
d<-crossTable(dbData[,1],dbData[,4],useNA = "ifany")
summary(d,latex=T)

# library(MASS)
# mytable <- xtabs(~cluster+iRsimple2, data=dbData, na.include)
# latex(mytable, booktabs=T,"bokktab.tex")
# tableNominal(vars = dbData, cap = "none", vertical = FALSE, lab = "tab: nominal1", longtable = FALSE)

# tableNominal(dbData[,3:4], weights = NA, subset = NA, 
#     group = NA, miss.cat = NA, print.pval = c("none", "fisher", 
#     "chi2"), pval.bound = 10^-4, fisher.B = 2000, vertical = F, 
#     cap = "", lab = "", col.tit.font = c("bf", "", "sf", "it", "rm"), 
#     font.size = "footnotesize", longtable = TRUE, nams = NA, 
#     cumsum = TRUE,file="~/MSc-Biostatistics/Thesis/Manuscript/tables/Descriptives.tex")

#######CODING IR

thesisdb.con <- dbConnect(MySQL(), user="R", password="1234", dbname="thesisdb", host="localhost")
dbData<- dbGetQuery(thesisdb.con, "SELECT * FROM DatasetFull WHERE region='Europe';")
irDesc<- dbGetQuery(thesisdb.con, "SELECT iRsimple2, iRdesc FROM thesisdb.iRmerge where iRdesc is not null;")
dbDisconnect(thesisdb.con)


dbData<-merge(irDesc,dbData[,c(4,6)],by.x="iRsimple2", by.y="iRsimple2", all.y=T)
mapping<-count(dbData)
mapping$percent <- round(mapping$freq/sum(mapping$freq)*100,2)
mapping[dim(mapping)[1]+1,4:5]<-round(sapply(mapping[,4:5],sum),0)
mapping<-mapping[,c(3,1,2,4,5)]
mapping[dim(mapping)[1],1]<-"Total"


latex(mapping,file="~/MSc-Biostatistics/Thesis/Manuscript_new/tables/iRmapping"
,rowname = NULL,labels = labels(results_tex),size = "scriptsize",where="htpb", extracolsize='tiny'
	,caption = "Initial risk group labelling, alongside final labelling for all sequences within European monophyletic clusters.Absolute and relative frequencies are displayed."
	      ,na.blank=F,insert.bottom=""
	      ,booktabs=T
	      ,dcolumn = F,extracolheads = c("Initial","Assigned","Description"),
	      colheads = c("Risk Group","Risk Group","Risk Group","N","\\%"))


t<-table(dbData[,4],useNA="ifany")
as.data.frame(prop.table(t))

t<-table(dbData[,2],useNA="ifany")
as.data.frame(prop.table(t))

tableNominal(vars = dbData[,c(4,6)], cap = "none", vertical = FALSE, lab = "tab: nominal1", longtable = T,font.size="tiny")

x<-
tableNominal(vars = dbData[,c(1,4)], group= dbData[,1],cap = "Table of nominal variables.", lab = "tab: nominal",vertical=F, cumsum=F,
longtable = F,font.size="tiny", print.pval="fisher")

tableNominal(vars = dbData[,c(2,4)],cap = "Table of nominal variables.", lab = "tab: nominal",vertical=F, cumsum=F,
longtable = T,font.size="tiny", print.pval="fisher")


ggplot(dbData, aes(x = iRsimple2, y = cluster)) + geom_point(size=2, aes(color = submissionCountry), position="dodge")

w<-xtabs(~iRsimple2 + cluster, data=dbData)
summary(w)

thesisdb.con <- dbConnect(MySQL(), user="R", password="1234", dbname="thesisdb", host="localhost")
iRmap<- dbGetQuery(thesisdb.con, "SELECT * FROM thesisdb.iRsummary where region = 'Europe';")
dbDisconnect(thesisdb.con)
iRmap<- iRmap[,-1]
iRmap<-count(iRmap)
iRmap<-iRmap[order(iRmap$iRsimple2),]


names(dbData)

#############MAP
paletteCols<-c("#BB1D1D","#5E0707","#6CBB1D","#6C1DBB","#325E07","#1DBBBB","#075E5E","#666666")


dfA<-count(dbData, c("submissionCountry", "iRsimple2"))
df3<-reshape(data = dfA, direction = "wide", timevar = "iRsimple2", idvar="submissionCountry",drop=NULL)
df3[is.na(df3)]<-0
df3<-df3[c(1, 5, 3:4, 6:9, 2)]


dfMap<-joinCountryData2Map(dF = df3, joinCode = "ISO3", "submissionCountry", nameCountryColumn = "submissionCountry", suggestForFailedCodes = T, verbose=T)
names(dfMap)<-gsub("freq.", "", names(dfMap))
dfMap$NAME<-gsub("N. Cyprus",NA,dfMap$NAME)
dfMap$NAME<-gsub("Macedonia","F.Y.R.O.M.",dfMap$NAME)


cairo_pdf("WorldSample.pdf", height = 10, width = 16,family="Linux Biolinum" ,pointsize=20)
mapPies( dfMap, nameX="LON", nameY="LAT", zColours = paletteCols, 
 nameZs=names(dfMap)[c(51:58)], barWidth = 1, landCol="grey90",
 ylim=c(50,57), xlim=c(-49,35), maxZVal=1, symbolSize=2,
#  mapRegion='europe', 
 addSizeLegend = T, barRelative=F,lwdSymbols=0.1)
 labelCountries( col = "#2F2733", cex = 0.3)
dev.off()