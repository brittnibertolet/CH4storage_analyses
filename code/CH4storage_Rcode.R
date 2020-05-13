rm(list=ls())
setwd("CH4storage_Analyses")

source("code/setPlotTheme.R")
library(reshape2)
library(RColorBrewer)
library(MASS)

####### Estimate CH4 storage for each lake-year #######

# Read in DIC data
dic=read.csv("data/hypoDIC_rawData.csv", stringsAsFactors = F)
# Average technical replicates of DIC concentrations 
dicM=aggregate(CH4_uM~dateSample+lakeID+year+doy, FUN=mean, data=dic)

# Read in lake area and volume data
lakeShape=read.csv("data/lakeVolumeArea.csv", stringsAsFactors = F)

# Get rates of CH4 storage from multiple regression 
fit=lm(CH4_uM~doy*lakeID*as.factor(year), data=dicM)
anova(fit)

coeff=as.data.frame(summary(fit)$coefficients)

B0=coeff$Estimate[1]
B1.doy=coeff$Estimate[2]
B2.lake=data.frame(lakeID=c("CR","EL", "HB", "MO", "WL"), estimate=c(0, coeff$Estimate[3:6]))
B3.year=data.frame(year=c("2014", "2015", "2016", "2018", "2019"), estimate=c(0, coeff$Estimate[7:10]))
B4.doy.lake=data.frame(lakeID=c("CR","EL", "HB", "MO", "WL"), estimate=c(0, coeff$Estimate[11:14]))
B5.doy.year=data.frame(year=c("2014", "2015", "2016", "2018", "2019"), estimate=c(0, coeff$Estimate[15:18]))
B6.lake.year=data.frame(lakeID=c(rep("CR", times=5), rep(c("EL", "HB", "MO", "WL"), times=5)), 
                        year=c(rep(c("2014", "2015", "2016", "2018", "2019"), each=1),rep(c("2014", "2015", "2016", "2018", "2019"), each=4)),   
                        estimate=c(rep(0, times=9), coeff$Estimate[19:30], 0, coeff$Estimate[31:33]))
B7.doy.lake.year=data.frame(lakeID=c(rep("CR", times=5), rep(c("EL", "HB", "MO", "WL"), times=5)), 
                            year=c(rep(c("2014", "2015", "2016", "2018", "2019"), each=1),rep(c("2014", "2015", "2016", "2018", "2019"), each=4)),
                            estimate=c(rep(0, times=9), coeff$Estimate[34:45], 0, coeff$Estimate[46:48]))

# Simulate line from regression coefficients 
lakes=as.character(unique(dicM$lakeID))
years=as.character(unique(dicM$year))
modelout=c()
modelcoeffs=c()
for(i in 1:length(lakes)){
  for(j in 1:length(years)){
    lakeID=as.character(lakes[i]); print(lakeID)
    year=as.character(years[j]); print(year)
    doy=125:260
    intercept= B0 + B2.lake$estimate[B2.lake$lakeID==lakes[i]] +
      B3.year$estimate[B3.year$year==years[j]] +
      B6.lake.year$estimate[B6.lake.year$lakeID==lakes[i] & B6.lake.year$year==years[j]]
    slope=     B1.doy+
      B4.doy.lake$estimate[B4.doy.lake$lakeID==lakes[i]] +
      B5.doy.year$estimate[B5.doy.year$year==years[j]] +
      B7.doy.lake.year$estimate[B7.doy.lake.year$lakeID==lakes[i] & B7.doy.lake.year$year==years[j]]
    
    #get lake volume and lake area from the lakeShape dataframe
    lakeV=lakeShape$hypoVol_m3[lakeShape$lakeID==lakes[i]]
    lakeA=lakeShape$sedArea_m2[lakeShape$lakeID==lakes[i]]
    cOut=data.frame(lakeID, year, intercept, slope=(slope*lakeV)/lakeA, stringsAsFactors = F)
    modelcoeffs=rbind(modelcoeffs, cOut)
    
    # get predicted line
    predCH4_uM = intercept + slope*doy
    predOut=data.frame(lakeID, year, doy, predCH4_uM, stringsAsFactors = F)
    modelout=rbind(modelout, predOut)
  }
}

# Get rid of model coefficients for EL 2019
modelcoeffs$intercept[modelcoeffs$year=="2019" & modelcoeffs$lakeID=="EL"]=NA
modelcoeffs$slope[modelcoeffs$year=="2019" & modelcoeffs$lakeID=="EL"]=NA
modelout$predCH4_uM[modelout$year=="2019" & modelout$lakeID=="EL"]=NA

# Plot Fig. 1 - CH4 time series
fig1=ggplot(dicM, aes(x=doy, y=CH4_uM))+
  facet_grid(lakeID~year, scales = "free")+
  stat_smooth(method="lm", fill="white", color="white")+
  xlab("Day of Year")+
  ylab(expression('Hypolimnion CH'['4']*' concentration (mmol L'^-1*')'))+
  geom_point()+
  geom_line(data=modelout, aes(x=doy, y=predCH4_uM))
#ggsave("Fig1.pdf", fig1, width=6.5, height=5)

####### Compare estimates of CH4 storage to other CH4 processes #######

# Read in West et al. 2016 data 
west=read.csv("data/CH4production_emissionRates_westetal2016.csv", stringsAsFactors = F, header=T)
west2=melt(west, id.vars = "lakeID")
CH4storage=modelcoeffs[,c("lakeID", "year", "slope")]

CH4budget=merge(west, CH4storage, by="lakeID", all.y = T)
CH4budget$CH4_storage=CH4budget$slope
CH4budget=CH4budget[,c("lakeID", "year", "CH4_prod", "CH4_diff", "CH4_tot_emiss",  "CH4_storage")]
CH4budgetagg=aggregate(.~lakeID, data=CH4budget[,c(1,3:6)], FUN = "mean", na.rm=T)
range(CH4budgetagg$CH4_storage/CH4budgetagg$CH4_prod, na.rm=T)

CH4budget=melt(CH4budget, id.vars = c("lakeID", "year"))

# Fix rate labels
CH4budget$variable=as.character(CH4budget$variable)
CH4budget$variable[CH4budget$variable=="CH4_prod"]="Production"
CH4budget$variable[CH4budget$variable=="CH4_diff"]="Diffusive emission"
CH4budget$variable[CH4budget$variable=="CH4_tot_emiss"]="Total emission"
CH4budget$variable[CH4budget$variable=="CH4_storage"]="Storage"

# Assign study
CH4budget$study=NA
CH4budget$study[CH4budget$variable%in%c("Production", "Diffusive emission", "Total emission")]="West et al. 2016"
CH4budget$study[CH4budget$variable%in%c("Storage")]="This study"

# Fix factors 
CH4budget$variable=factor(CH4budget$variable, levels=c("Production", "Diffusive emission", "Total emission", "Storage"))
CH4budget$study=factor(CH4budget$study, levels=c("West et al. 2016", "This study"))

# Plot Fig. 2 - CH4 budgets
fig2=ggplot(CH4budget[CH4budget$lakeID!="EL",], aes(x=variable, y=value))+
  stat_summary(geom="bar", fun.y="mean", aes(fill=study))+
  stat_summary(data=CH4budget[CH4budget$variable=="Storage" & CH4budget$lakeID!="EL",], geom="errorbar", fun.data="mean_se")+
  #geom_point(color="red")+
  facet_grid(.~lakeID, scales="free_x")+
  ylab(expression('Areal rate (mmol CH'['4']*' m'^-2*' d'^-1*')'))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  scale_fill_manual(values=c("grey70", "black"))+
  guides(fill=F)
#ggsave("Fig2.pdf", height=3.5, width=4)

####### Is CH4 storage related to GPP across all lakes and years? #######
# Read in GPP
gpp=read.csv("data/GPP_rawData.csv", stringsAsFactors = F)
# Convert to areal terms 
gpp$GPP_areal=gpp$GPP*gpp$zMix*1000/32 # mmol C m-2 d-1
# Calculate mean summer GPP
gpp.mean=aggregate(GPP_areal~lakeID+year, data=gpp, FUN=mean)

# Create data frame with summer CH4 storage and mean summer GPP
lakeVars=merge(modelcoeffs, gpp.mean, by=c("lakeID", "year"), all.x=T)

# Linear regression
fit1=lm(slope~GPP_areal, data=lakeVars) # CH4 storage ~ GPP
summary(fit1)

# Plot Fig. 3a - GPP-only model
fig3a=ggplot(lakeVars, aes(x=GPP_areal, y=slope))+
  stat_smooth(method="lm", color="black", se=F, size=0.5)+
  geom_point(aes(color=lakeID))+
  scale_y_continuous(limits=c(-1,9), breaks=c(0, 2, 4, 6, 8))+
  ylab(expression("CH"['4']* ' storage (mmol CH'['4']* ' m'^-2*' d'^-1*')'))+
  xlab(expression("mean GPP (mmol C m"^-2*' d'^-1*')'))+
  scale_color_brewer(name="Lake ID",type="qual", palette=6)+
  theme(legend.position = c(0.25, 0.8), 
        legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"))

####### Is CH4 storage related to GPP witin a lake? #######
# ANCOVA design regression
fit2=lm(slope~GPP_areal*lakeID, data=lakeVars) # CH4 storage ~ GPP * lakeID
anova(fit2)

# What is the best predictive model?
# Stepwise-regression:
step.model=stepAIC(fit2, direction="both", trace=F)
summary(step.model)
step.model$anova

# Get coefficients from lakeID-only model
fit3=lm(slope~lakeID, data=lakeVars)
coefficients(fit3)

# Plot Fig. 3b - intercept only model
pred=rbind(data.frame(lakeID=c("CR"), GPP_areal=14:32, intercept=0.035),
           data.frame(lakeID=c("EL"), GPP_areal=27:38, intercept=1.89),
           data.frame(lakeID=c("HB"), GPP_areal=24:34, intercept=1.92),
           data.frame(lakeID=c("MO"), GPP_areal=27:59, intercept=3.91),
           data.frame(lakeID=c("WL"), GPP_areal=20:29, intercept=1.84))
fig3b=ggplot(lakeVars, aes(x=GPP_areal, y=slope))+
  geom_line(data=pred, aes(x=GPP_areal, y=intercept, color=lakeID))+
  geom_point(aes(color=lakeID))+
  ylab(expression("CH"['4']* ' storage (mmol CH'['4']* ' m'^-2*' d'^-1*')'))+
  xlab(expression("mean GPP (mmol C m"^-2*' d'^-1*')'))+
  scale_color_brewer(name="Lake ID",type="qual", palette=6)+
  scale_y_continuous(limits=c(-1,9), breaks=c(0, 2, 4, 6, 8))+
  guides(color = F)

# Plot complete Fig. 3
fig3=plot_grid(fig3a, fig3b, labels = c("A", "B"))
#ggsave("Fig3.pdf", height=3.5, width=6)

# Stepwise regression
step.model=stepAIC(fit2, direction="both", trace=F)
summary(step.model)
step.model$anova

# Plot Fig 2 - GPP-only model and lakeID-only model

####### Simulation experiment #######


# linear 20 increase
GPP.pred=c(rep(20, times=49),seq(20,40, length.out = 51))
year=1:length(GPP.pred)
thought.pred=data.frame(year, GPP.pred, CH4stor1=NA, CH4stor2=NA, 
                        CH4stor3=NA, CH4stor4=NA, CH4stor5=NA, CH4stor6=NA)
thought.pred$CH4stor1=fit1$coefficients[1]+thought.pred$GPP.pred*fit1$coefficients[2]

delays=c(10, 20, 30, 40, 50)
for(j in 1:length(delays)){
  for(i in (delays[j]+1):length(year)){
    thought.pred[i,j+3]=fit1$coefficients[1]+mean(GPP.pred[(i-delays[j]):i])*fit1$coefficients[2]
  }
}
# change years to start at year 50
thought.pred$year=thought.pred$year-50

p1=ggplot(thought.pred[thought.pred$year>-10,], aes(year, GPP.pred))+geom_line()+
  xlab("Time")+ylab("Lake GPP")

thought.pred50=thought.pred[thought.pred$year>0,]
cumCH4=data.frame(lag=c(0, 10, 20, 30, 40, 50), CH4_sum=colSums(thought.pred50[,3:8], na.rm = T))
p3=ggplot(cumCH4, aes(x=lag, y=CH4_sum))+geom_point(aes(shape=as.factor(lag)))+
  xlab("Years represented in\nmethanogenic sediments (x)")+ylab(expression("Sum of summer CH"['4']*" storage"))+
  guides(shape=F)
range(cumCH4$CH4_sum)

thought.pred=melt(thought.pred, id.vars = "year")

p2=ggplot(thought.pred[thought.pred$variable!="GPP.pred" & thought.pred$year>-10,], aes(year, value))+
  geom_line(aes(linetype=variable))+
  xlab("Time")+ylab(expression("CH"['4']*" storage"))+
  scale_linetype_discrete(name="", labels=c(0, 10, 20, 30, 40, 50))+
  theme(legend.position = c(0, 1), 
        legend.justification = c(-0.1, 0.9),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.background = element_blank())
p2

plot_grid(p1, p2, nrow=1, labels=c("A", "B"))

linear20=plot_grid(plot_grid(p1, p2, nrow=2, labels=c("A", "B")),
                   p3, ncol=2, labels=c("", "C"))
linear20
#ggsave("Fig4.pdf", height=5, width=6.5)



