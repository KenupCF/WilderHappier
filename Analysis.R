## Loading packages
require(plyr)
require(dplyr)
require(magrittr)
require(reshape2)
require(devtools)
require(data.table)
require(AICcmodavg)
require(ggplot2)
require(combinat)

setwd() #Must be manually added!

# Importing helper functions
source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")

# Loading HDI+OLSI data
hdi.olsi<-read.csv('.\\HDI+OLSI.csv')

# Converting data to numeric
hdi.olsi$HDI<-force.numeric(gsub(x=hdi.olsi$HDI," ",""))
hdi.olsi$LEB<-force.numeric(gsub(x=hdi.olsi$LEB," ",""))
hdi.olsi$EYS<-force.numeric(gsub(x=hdi.olsi$EYS," ",""))
hdi.olsi$MYS<-force.numeric(gsub(x=hdi.olsi$MYS," ",""))
hdi.olsi$GNI<-force.numeric(gsub(x=hdi.olsi$GNI," ",""))
hdi.olsi$OLSI<-force.numeric(gsub(x=hdi.olsi$OLSI," ",""))

# Capping GNI e EYS
hdi.olsi$GNI.cap<-hdi.olsi$GNI
hdi.olsi$GNI.cap[hdi.olsi$GNI.cap>75000]<-75000
hdi.olsi$EYS.cap<-hdi.olsi$EYS
hdi.olsi$EYS.cap[hdi.olsi$EYS.cap>18]<-18

# Calculating socio-economic indexes
hdi.olsi%<>%mutate(Health.index=((LEB-20)/(85-20)))
hdi.olsi%<>%mutate(Income.index=((log(GNI.cap)-log(100))/(log(75000)-log(100))))
hdi.olsi%<>%mutate(
  MYS.index=(MYS-0)/(15-0),
  EYS.index=(EYS.cap-0)/(18-0))
hdi.olsi%<>%mutate(Education.index=(MYS.index+EYS.index)/2)

# Removing countries without OLSI data
hdi.olsi.bkp<-hdi.olsi
hdi.olsi%<>%filter(!is.na(OLSI))

# Loading PLC data
areas<-read.csv('.\\AREAS.csv')	##Total Land Area of Countries
doms<-read.csv('.\\DOMS.csv')	##Table containing Country Identifiers	

# Tidying PLC data frame
areas$PA_AREA<-force.numeric(gsub(x=areas$PA_AREA," ",""))
areas$TOTAL_AREA<-force.numeric(gsub(x=areas$TOTAL_AREA," ",""))
areas$COUNTRY_DOM%<>%tolower
PLC.table<-merge(areas,doms,by='COUNTRY_DOM')

# Generating PLC index
PLC.table$PLC.percent<-(PLC.table$PA_AREA*100)/(PLC.table$TOTAL_AREA)
PLC.table%<>%mutate(PLC.index=((PLC.percent-0)/(50-0)))

# Loading NLC data
lc<-read.csv('.\\LC.csv')	
pixels<-read.csv('.\\PIXELS.csv')
names<-read.csv('.\\NAMES.csv')

# Tidying data frames
names$COUNTRY_NAME%<>%trim
pixels.names<-merge(pixels,names,by='COUNTRY_ID')
lc.names<-merge(lc,names,by='COUNTRY_ID')

# Summarising NLC data
pixels_by_country<-lc.names%>%
	dplyr::group_by(COUNTRY_NAME)%>%
	dplyr::summarise(TOTAL.count=sum(PIXELS))
	
pixels_by_lc_by_country<-acast(
	data=lc.names,
	formula=COUNTRY_NAME~LC,
	value.var='PIXELS',
	fun.aggregate=sum)

# Simplifying Landcover data
lc.types<-colnames(pixels_by_lc_by_country)
ALC<-c('10','11','12','20','30','40','190') #Defining artificial landcovers
NLC<-lc.types[!lc.types%in%c(ALC)] 			#Defining natural landcovers

# Calculating total NLC
NLC.pixels<-apply(pixels_by_lc_by_country[,NLC],1,sum)
NLC.pixels<-data.frame(NLC.pixels=NLC.pixels,COUNTRY_NAME=names(NLC.pixels))
pixels_by_country%<>%merge(NLC.pixels,by='COUNTRY_NAME')
pixels_by_country$NLC.percent<-force.numeric(pixels_by_country$NLC.pixels)*100/force.numeric(pixels_by_country$TOTAL.count)
NLC.table<-merge(pixels_by_country,hdi.olsi,by='COUNTRY_NAME')

# Generating NLC index
NLC.table%<>%mutate(NLC.index=((NLC.percent-0)/(100-0)))

# Merging HDI.OLS + PLC.table + NLC.table
resu<-hdi.olsi.bkp%>%
	merge(PLC.table)%>%
	merge(NLC.table)%>%
	data.table(key=c('COUNTRY_NAME'))

# Generating different HDI indexes
#####Geometric HDI
resu$HDI.geo <- (resu$Education.index * resu$Income.index * resu$Health.index)^(1/3)
#####Arithmetic HDI
resu$HDI.ari <- (resu$Education.index + resu$Income.index + resu$Health.index)/3
#####Arithmetic HDI+PLC
resu$HDI.PLC <- (resu$Education.index + resu$Income.index + resu$Health.index + resu$PLC.index)/4
#####Arithmetic HDI+NLC
resu$HDI.NLC <- (resu$Education.index + resu$Income.index + resu$Health.index + resu$NLC.index)/4
#####Geometric HDI+NLC
resu$HDI.NLCgeo <- (resu$Education.index * resu$Income.index * resu$Health.index * resu$NLC.index)^(1/4)
#####Geometric HDI+PLC
resu$HDI.PLCgeo <- (resu$Education.index * resu$Income.index * resu$Health.index * resu$PLC.index)^(1/4)

# Calculating differences in HDI
resu$HDI.diff<-resu$HDI.NLC-resu$HDI.geo
resu$HDI.diffgeo<-resu$HDI.NLCgeo-resu$HDI.geo
resu%<>%
	dplyr::arrange(desc(HDI.geo))%>%
	dplyr::mutate(HDI.rank=1:n())%>%
	dplyr::arrange(desc(HDI.NLCgeo))%>%
	dplyr::mutate(
		HDI.NLC.rank=1:n(),
		change.in.rank=HDI.rank-HDI.NLC.rank)

###############
# Running GLMs#
###############

HDIgeo<-glm(OLSI~HDI.geo,Gamma(),data=resu)
HDIari<-glm(OLSI~HDI.ari,Gamma(),data=resu)
HEI<-glm(OLSI ~ Health.index + Education.index  + Income.index, Gamma(), data=resu)
HEIN<-glm(OLSI ~ Health.index + Education.index + Income.index  + NLC.index, Gamma(), data=resu)
HEIP<-glm(OLSI ~ Health.index + Education.index + Income.index + PLC.index, Gamma(), data=resu)

# Model selection
models<-list(HDIgeo,
	HDIari,HEI,HEIN,HEIP)
model.names<-c("HDIgeo",
	"HDIari","HEI","HEIN","HEIP")
res.table<-aictab(cand.set=models,modnames=model.names,second.ord=TRUE)
res.table

# Subsetting data by HDI level
VeryHigh<-resu$COUNTRY_NAME[which(resu$HDI>=0.8)]
High<-resu$COUNTRY_NAME[which (resu$HDI<0.8 & resu$HDI>=0.7)]
Middle<-resu$COUNTRY_NAME[which (resu$HDI<0.7 & resu$HDI>=0.55)]
Low<-resu$COUNTRY_NAME[which (resu$HDI<0.55)]

resu.VeryHigh<-resu%>%
  dplyr::filter(COUNTRY_NAME%in%VeryHigh)
resu.High<-resu%>%
  dplyr::filter(COUNTRY_NAME%in%High)
resu.Middle<-resu%>%
  dplyr::filter(COUNTRY_NAME%in%Middle)
resu.Low<-resu%>%
  dplyr::filter(COUNTRY_NAME%in%Low)

# Running GLMs with subsetted data
models.vhigh<-list(
  HEI=glm(OLSI ~ Health.index  +  Education.index  +  Income.index, Gamma(), data=resu.VeryHigh),
  HEIN=glm(OLSI ~ Health.index  +  Education.index  +  Income.index  +  NLC.index, Gamma(), data=resu.VeryHigh))
models.high<-list(
  HEI=glm(OLSI ~ Health.index  +  Education.index  +  Income.index, Gamma(), data=resu.High),
  HEIN=glm(OLSI ~ Health.index  +  Education.index  +  Income.index  +  NLC.index, Gamma(), data=resu.High))
models.mid<-list(
  HEI=glm(OLSI ~ Health.index  +  Education.index  +  Income.index, Gamma(), data=resu.Middle),
  HEIN=glm(OLSI ~ Health.index  +  Education.index  +  Income.index  +  NLC.index, Gamma(), data=resu.Middle))
models.low<-list(
  HEI=glm(OLSI ~ Health.index  +  Education.index  +  Income.index, Gamma(), data=resu.Low),
  HEIN=glm(OLSI ~ Health.index  +  Education.index  +  Income.index  +  NLC.index, Gamma(), data=resu.Low))

# Model selection
res.table.vhigh<-aictab(cand.set=models.vhigh,modnames=names(models.vhigh),second.ord=TRUE)
res.table.high<-aictab(cand.set=models.high,modnames=names(models.high),second.ord=TRUE)
res.table.mid<-aictab(cand.set=models.mid,modnames=names(models.mid),second.ord=TRUE)
res.table.low<-aictab(cand.set=models.low,modnames=names(models.low),second.ord=TRUE)

#######################
##Running null models##
#######################

temp.df<-resu%>%filter(!is.na(resu$OLSI)) #data.frame with only valid OLSI values

null.mod.list<-list() #list to hold null iterations
reposition<-TRUE #parameter defining if sampling is with replacement
n.iter<-1e4 #number of iterations run

#Running iterations to generate null data
for (i in 1:n.iter){
	null.mod.list[[i]]<-temp.df
	null.mod.list[[i]]$NLC.index<-NULL
	null.mod.list[[i]]$NLC.index<-sample(temp.df$NLC.index,replace=reposition)
}

#Running models from null ndata
null.mods<-lapply(null.mod.list,function(x){
	temp.resu<-glm(OLSI ~ Health.index  +  Education.index  +  Income.index + NLC.index, Gamma(), data=x)
	return(temp.resu)
	})

#Model table	
null.modsel<-aictab(cand.set=null.mods,modnames=paste0('null model ',1:n.iter),second.ord=TRUE)

#Proportion of null models with AICc lower than true AICc+2
null.p<-(sum(null.modsel$AICc<(res.table[1,'AICc'])+2))/n.iter

########################################
#Testing correlation between variables##
########################################
vars<-c('NLC.index','Health.index','Education.index','Income.index') #Variables
comb<-combn(vars,2) #combinations of variable pairs
corr.list<-lapply(1:ncol(comb),function(x){ #list of correlation tests
	cor.test(resu[,comb[1,x]],resu[,comb[2,x]])}) 
names(corr.list)<-apply(comb,2,paste0,collapse="-")
#creating dataframe of correlation tests
corr.p.df<-data.frame(pair=names(corr.list),
	p.value=sapply(corr.list,function(x){x$p.value}), 
	r=sapply(corr.list,function(x){x$estimate}))

###Analysing variance in HDI levels
{
resu.VeryHigh$HDI.lvl<-'Very High'
resu.High$HDI.lvl<-'High'
resu.Middle$HDI.lvl<-'Middle'
resu.Low$HDI.lvl<-'Low'

resu2<-rbind.fill(resu.VeryHigh,resu.High,resu.Middle,resu.Low)%>%

	
write.csv(resu2, file = "RESU2.csv", row.names = TRUE)
	
	
var.df<-rbind.fill(resu2%>%
	dplyr::group_by(HDI.lvl)%>%
	dplyr::summarise(
		Health.var=var(Health.index),
		Income.var=var(Income.index),
		NLC.var=var(NLC.index),
		Education.var=var(Education.index)),
	resu2%>%dplyr::summarise(
		Health.var=var(Health.index),
		Income.var=var(Income.index),
		NLC.var=var(NLC.index),
		Education.var=var(Education.index)))
var.df%<>%mutate(HDI.lvl= replace(HDI.lvl, is.na(HDI.lvl), 'All')) 

write.csv(var.df, "Variance Indexes.csv",row.names=TRUE)

}
