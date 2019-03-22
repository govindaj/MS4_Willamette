#August 13th, 2018

#TSS-THg analysis for MS4 Phase I stormwater data including Clean Water Services using Linear mixed-effects models

#Janani Govindan, Oregon Department of Environmental Quality


#Load libraries----


#To do ggplot

library(ggplot2)

#To do ggarrange

library(ggpubr)

#diagnostic plots and box-cox transformation
library(MASS)

#To organize data

library(plyr)

#To work with Linear Mixed Effects Models
library(lme4)

#-to get an R-Squared value for the mixed-effects linear model
library(MuMIn)

#To generate summaries of both models to compare-huxtable

library(huxtable)

#use summary table-export_summs

library(jtools)

#For Prediction intervals

library(merTools)


#Do diagnostic check for both LMEM's

library(nlme)
library(Matrix)
library(parallel)
library(gridExtra)
library(cowplot)

#broom
library(broom)


#Diagnostic Check
library(predictmeans)

# #To zoom into an out of a plot 
# library(zoom)

#lmerTest

library(lmerTest)

#to use the refernce grid function and emmeans

library(emmeans)

#to get standardized residuals

library(HLMdiag)

#To get the influential points

library(influence.ME)

#identify influential points

library(car)

#Zoom in function

library(zoom)

#To use methods to handle non-detects

library(NADA)
library(Rcmdr, quietly=TRUE)
library(EnvStats)


#attach the Clean Water Services data-----

#To read excel spreadsheet
X3Sites_CWS<-read.csv("3Sites_CWS.csv")

head(X3Sites_CWS)


#rename certain columns----

names(X3Sites_CWS)[13]<-"TSS_CONC_mg_L"
names(X3Sites_CWS)[10]<-"THg_CONC_mg_L"
names(X3Sites_CWS)[1]<-"ORGANIZATION_PROVIDING_DATA"
head(X3Sites_CWS)

#make TSS column from factor to numeric----
X3Sites_CWS$TSS_CONC_mg_L<-as.numeric(X3Sites_CWS$TSS_CONC_mg_L)

#Substitution of TSS non-detects (mg/L)----
#Already knew that there were 2 non-detect values for TSS and the location of the 2 non-detects within the dataset
#So just subsittiuted those values for TSS Concentrations (mg/L) at values half the detection limit


X3Sites_CWS$TSS_CONC_mg_L[50]<-1
X3Sites_CWS$TSS_CONC_mg_L[125]<-1 


#Subset both Gresham and CWs data to do Kaplan-Meier Method of handling non-detect values----
# DL_final_removed_3_HUC<-DL_final_removed[!(DL_final_removed$HUC_8=="17090002"|DL_final_removed$HUC_8=="17090004"|DL_final_removed$HUC_8=="17090012"|DL_final_removed$HUC_8=="17090005"|DL_final_removed$HUC_8=="17090010"|DL_final_removed$HUC_8=="17090001"),]

Gresham_THG<-X3Sites_CWS[(X3Sites_CWS$Place=="Gresham"),]
CWS_THG<-X3Sites_CWS[(X3Sites_CWS$Place=="Clean Water Services"),]

#Create a Gresham dataframe with just the THg concentrations, censored data, and detection limit data columns----


Gresham_THG$Numeric_THG<-Gresham_THG$THg_CONC_mg_L


CWS_THG$Numeric_THG<-CWS_THG$THg_CONC_mg_L




#Kaplan-Meier for Gresham data----
#with NADA package

nada_kmG_aveTHG<-cenfit(Gresham_THG$Numeric_THG, Gresham_THG$censor_data)
print(nada_kmG_aveTHG)


#NADA package gives a 5.93 x 10^-6 mg/L average for Gresham non-detect points----
#K_M with EnvStats package to compare (bootstrapping)

km_Gresham<-enparCensored(Gresham_THG$Numeric_THG, Gresham_THG$censor_data, ci=T, ci.method="bootstrap", n.bootstraps = 3000)

print(km_Gresham)


# MLE method for handling non-detects for Gresham

mle_summary<-enormCensored(Gresham_THG$Numeric_THG, 
                           Gresham_THG$censor_data)

print(mle_summary)

#5.88 x 10^-6 mg/L


#mle_Gresham<-enormCensored(Gresham_THG$Numeric_THG, Gresham_THG$censor_data, method="impute.w.mle", lb.impute=min(Gresham_THG$d))

#get the same mean (5.93 x 10 ^-6 mg/L) with bootstrapping
#Kaplan-Meier for CWS data--------

nada_kmCWS_aveTHG<-cenfit(CWS_THG$Numeric_THG, CWS_THG$censor_data)
print(nada_kmCWS_aveTHG)


#Nada package gave a 4.79 x 10^ -6 mg/ L mean for CWS non-detect points

km_CWS<-enparCensored(CWS_THG$Numeric_THG, CWS_THG$censor_data, ci=T, ci.method="bootstrap", n.bootstraps = 3000)


print(km_CWS)


#Bootstrapping gives 4.78* 10 ^-6 mg/L for CWS non-detect points


#MLE for CWS

mle_CWS<-enormCensored(CWS_THG$Numeric_THG, CWS_THG$censor_data)


#Do K-M for both Gresham and CWS


#combine the columns of both the Gresham and CWS dataframe into one dataframe----

Gresham_CWS_both_THG<-X3Sites_CWS[(X3Sites_CWS$Place=="Gresham"|X3Sites_CWS$Place=="Clean Water Services"),]
Gresham_CWS_both_THG$Numeric_THG<-Gresham_CWS_both_THG$THg_CONC_mg_L


Gresham_CWS_both_THG_nada<- cenfit(Gresham_CWS_both_THG$Numeric_THG, Gresham_CWS_both_THG$censor_data)

Gresham_CWS_both_THG$Numeric_THG
Gresham_CWS_both_THG$censor_data


print(Gresham_CWS_both_THG_nada)

#Average is 5.38 x 10^-6 mg/L


km_Gresh_CWS<-enparCensored(Gresham_CWS_both_THG$Numeric_THG, Gresham_CWS_both_THG$censor_data, ci=T, ci.method="bootstrap", n.bootstraps = 3000)


print(km_Gresh_CWS)


#Change the concentrations of THg to match the average conc. based on the Kaplan-meier method----

X3Sites_CWS$THg_CONC_mg_L[90]<-5.38*10^-6
X3Sites_CWS$THg_CONC_mg_L[91]<-5.38*10^-6
X3Sites_CWS$THg_CONC_mg_L[114]<-5.38*10^-6
X3Sites_CWS$THg_CONC_mg_L[125]<-5.38*10^-6
X3Sites_CWS$THg_CONC_mg_L[194]<-5.38*10^-6
X3Sites_CWS$THg_CONC_mg_L[195]<-5.38*10^-6
X3Sites_CWS$THg_CONC_mg_L[223]<-5.38*10^-6
X3Sites_CWS$THg_CONC_mg_L[225]<-5.38*10^-6


#Exploratory Analysis after accounting for non-detects----
#create variables for the untransformed TSS and THg concentrations and convert to numeric----

raw_THG_CWS<-as.numeric(X3Sites_CWS$THg_CONC_mg_L)
raw_TSS_CWS<-as.numeric(X3Sites_CWS$TSS_CONC_mg_L)


#Box-cox transformation of response variable

par(mfrow=c(1,2))

boxcox(raw_THG_CWS~raw_TSS_CWS, data = X3Sites_CWS)

boxcox(raw_THG_CWS~1, data =X3Sites_CWS)



#Box-cox transformation tells you that the response variable should be log-transformed so:----
#TSS and THg data were converted to numeric value, log transformed, and a column of log data was added to db

AllTSS_numeric<-as.numeric(X3Sites_CWS$TSS_CONC_mg_L)
AllTHG_numeric<-as.numeric(X3Sites_CWS$THg_CONC_mg_L)


#log transformation 

AllTSS_CWS<-log10(AllTSS_numeric)
AllTHG_CWS<-log10(AllTHG_numeric)


#Make into a column in data frame

X3Sites_CWS$All_TSS_permittees<-AllTSS_CWS
X3Sites_CWS$All_THG_permittees<-AllTHG_CWS


#Remove row 16 with the Undetected (TSS data)
#AllTSS_without_U<-AllTSS_CWS[-c(16)]
#AllTHG_without_U<-AllTHG_CWS[-c(16)]


#Summary Statistics for the raw (untransformed data)

#Summary of TSS concentrations
summary(raw_TSS_CWS)

#Summary of THg concentrations
summary (raw_THG_CWS)



#Make variables into a factor for visual analysis and for lmer models (i.e. ggplot-boxplot or scatterplot)----

#Organization as a factor
Organization_CWS<-as.factor(`ORGANIZATION PROVIDING DATA`)
Organization_CWS

#Nested Random effects of municipality and sites as a factor
#Add as columns to the data

X3Sites_CWS$CWS_city_sample<-as.factor(X3Sites_CWS$CITY)

X3Sites_CWS$CWS_site_sample<-as.factor(X3Sites_CWS$SITE)

#Make a explicitly nested variable as a factor

sample_CWS<-within(sample_CWS<-factor(CWS_site_sample:CWS_city_sample))



#Exploratory Analysis (Table 2)----

#Summary Statistics in mg/L

summary(raw_TSS_CWS)

summary(raw_THG_CWS)

#Visual exploratory analysis of each Phase I MS4 Permittee group (Figures 2 and 3)----
#Boxplot and violin plots (ggplot) for the Phase I MS4 Permittees

#For THG

theme_set((theme_gray(base_size = 20)))
box_violin1<-ggplot(X3Sites_CWS)+
  geom_violin(aes(x=Place, y=raw_THG_CWS, fill=City_name, color=City_name), scale = "width")+
  geom_boxplot(aes(x=Place, y=raw_THG_CWS), width=0.1)+
  xlab("Phase I MS4 Permittees")+ ylab("Untransformed THg Concentrations (mg/L)")+labs(fill="Phase I MS4 Permittees")+guides(color=FALSE)+
  theme(axis.title.x = element_blank())

#For TSS

box_violin2<-ggplot(X3Sites_CWS)+
  geom_violin(aes(x=Place, y=raw_TSS_CWS, fill=City_name, color=City_name), scale = "width")+
  geom_boxplot(aes(x=Place, y=raw_TSS_CWS), width=0.1)+
  xlab("Phase I MS4 Permittees")+ ylab("Untransformed TSS Concentrations (mg/L)")+labs(fill="Phase I MS4 Permittees")+guides(color=FALSE)


#combine both in grid arrange-can also do facet wrap

#grid.arrange(box_violin1, box_violin2)

ggarrange(box_violin1, box_violin2, ncol=1, nrow=2, common.legend=TRUE, legend = "bottom" )



#Scatterplot of TSS-THG concentrations in mg/L (Figure 2)----

#make the OLS model with no NA's

raw_TSS_CWS_no_NA<-as.numeric(X3Sites_CWS$TSS_CONC_mg_L[-c(16)])
raw_THG_CWS_no_NA<-as.numeric(X3Sites_CWS$THg_CONC_mg_L[-c(16)])

X3Sites_CWS_df<-X3Sites_CWS[-c(16),]
lm_CWS_full_no_NA<-lm(All_THG_no_NA_CWS~All_TSS_no_NA_CWS, data=X3Sites_CWS_df)
All_TSS_no_NA_CWS<-log10(raw_TSS_CWS_no_NA)
All_THG_no_NA_CWS<-log10(raw_THG_CWS_no_NA)
X3Sites_CWS_df$City_name_no_NA<-X3Sites_CWS_df$City_name

raw_TSS_CWS_no_NA

ggplot(lm_CWS_full_no_NA, aes(x=raw_TSS_CWS_no_NA, y=raw_THG_CWS_no_NA))+
  geom_point(aes(shape= X3Sites_CWS_df$City_name_no_NA, colour =  X3Sites_CWS_df$City_name_no_NA, size=3)) +
  guides(size=FALSE)+
  scale_fill_manual("Interval", values=c("grey33"))+
  scale_shape_manual(values=c(8,9,15,16,17,18,19,20))+xlab("TSS Concentrations (mg/L)") + 
  ylab("THg Concentrations (mg/L)")+
  scale_x_log10()+scale_y_log10()+
  labs(shape= "Phase I MS4 Permittees", color="Phase I MS4 Permittees")+
  guides(shape=guide_legend(override.aes = list(size=4)))



#Model Equations (Results Section and Tables 3,4,5, and 6)----

#OLS Model and summary with R-square values

lm_CWS_full<-lm(AllTHG_CWS~AllTSS_CWS, data=X3Sites_CWS)
summary(lm_CWS_full)

#LME Models

#Add columns to the table of both city and sites as a factor

X3Sites_CWS$CWS_city_sample<-as.factor(X3Sites_CWS$CITY)
X3Sites_CWS$CWS_site_sample<-as.factor(X3Sites_CWS$SITE)

#LME Model with Nested Sites as a Random Effect and it's summary results

mixed.lmer6_CWS<-lmer(AllTHG_CWS~AllTSS_CWS + (1|CWS_city_sample/CWS_site_sample), data=X3Sites_CWS)
mixed.lmer6_CWS
plot(mixed.lmer6_CWS)
summary(mixed.lmer6_CWS)

#R-Squared value for lmer 6

r.squaredGLMM(mixed.lmer6_CWS)

#Summary Results for lmer 6

summary (mixed.lmer6_CWS)


#LME Model with Permittees as a Random Effect

mixed.lmer7_CWS<-lmer(AllTHG_CWS~AllTSS_CWS + (1|CWS_city_sample), data=X3Sites_CWS)

#Summary results of lmer 7

summary(mixed.lmer7_CWS)


#R-squared value for lmer 7

r.squaredGLMM(mixed.lmer7_CWS)


#Diagnostic Plots for OLS and LME Models (Figures 7,8, and 9)----

#OLS Model Diagnostic Check (Figure 7)

par(mfrow=c(2,2))

plot(lm_CWS_full)
summary(lm_CWS_full)


#LME Model with nested sites as a random factor (Figure 9)-lmer 6

par(mfrow=c(2,2))

#Residuals vs fitted plot
plot(fitted(mixed.lmer6_CWS), residuals(mixed.lmer6_CWS), xlab="Fitted values", ylab="Residuals", main="Residuals vs Fitted")
abline(h=0)
lines(smooth.spline(fitted(mixed.lmer6_CWS), residuals(mixed.lmer6_CWS)), col="red")

#QQ-plot-do standardized residuals
qqnorm(AllTHG_CWS)
qqline(AllTHG_CWS)

#Scale-location Plot

MSE<-mean(resid(mixed.lmer6_CWS)^2)

Std_residuals<-as.numeric(resid(mixed.lmer6_CWS))/(sqrt(MSE))
Std_residuals
abs_Std_residuals<-abs(Std_residuals)
abs_Std_residuals
square_root_Std_residuals<-sqrt(abs_Std_residuals)

plot(fitted(mixed.lmer6_CWS), square_root_Std_residuals, xlab="Fitted values", ylab=expression(sqrt("Standardized residuals")), main="Scale-location")
lines(smooth.spline(fitted(mixed.lmer6_CWS), square_root_Std_residuals),col="red")


#Cook's Distance

cd<-cooks.distance(mixed.lmer6_CWS)
plot(cd,type="h", xlab="Num. of Observations", ylab="Cook's Distance", ylim=range(cd)*c(1,1.1), main="Cook's Distance Plot") #spread out the y-axis
threshold<-0.01680672 #threshold is 4/n(4/238)
lab<-cd>threshold
mylabel_sig_outliers<-1:length(cd) #label non-missing data (no NA data included)
text(which(lab), cd[lab], labels=mylabel_sig_outliers[lab], pos=3,cex=0.7,col="black")





#LME Model with Permittees as a Random Effect----

mixed.lmer7_CWS<-lmer(AllTHG_CWS~AllTSS_CWS + (1|CWS_city_sample), data=X3Sites_CWS)


#LME Model with Permittees as a random factor (Figure 8)-lmer 7----


par(mfrow=c(2,2))


#Residuals vs fitted plot
plot(fitted(mixed.lmer7_CWS), residuals(mixed.lmer7_CWS), xlab="Fitted values", ylab="Residuals", main="Residuals vs Fitted")
abline(h=0)
lines(smooth.spline(fitted(mixed.lmer7_CWS), residuals(mixed.lmer7_CWS)), col="red")

#QQ-plot-do standardized residuals
qqnorm(AllTHG_CWS)
qqline(AllTHG_CWS)

#Scale-location plot

MSE<-mean(resid(mixed.lmer7_CWS)^2)

Std_residuals<-as.numeric(resid(mixed.lmer7_CWS))/(sqrt(MSE))
Std_residuals
abs_Std_residuals<-abs(Std_residuals)
abs_Std_residuals
square_root_Std_residuals<-sqrt(abs_Std_residuals)

plot(fitted(mixed.lmer7_CWS), square_root_Std_residuals, xlab="Fitted values", ylab=expression(sqrt("Standardized residuals")), main="Scale-location")
lines(smooth.spline(fitted(mixed.lmer7_CWS), square_root_Std_residuals),col="red")


#Cook's Distance
cd<-cooks.distance(mixed.lmer7_CWS)
plot(cd,type="h", xlab="Num. of Observations", ylab="Cook's Distance", ylim=range(cd)*c(1,1.1), main="Cook's Distance Plot") #spread out the y-axis
threshold<-0.01680672 #threshold is 4/n(4/238)
lab<-cd>threshold
mylabel_sig_outliers<-1:length(cd) #label non-missing data (no NA data included)
text(which(lab), cd[lab], labels=mylabel_sig_outliers[lab], pos=3,cex=0.7,col="black")


#Prediction Intervals for OLS and LME Models (Figures 4,5, and 6)----

#Prediction interval plot for OLS Model (Figure 4)

plot(AllTSS_CWS, AllTHG_CWS, ylim=c(-6.5,-4.5), xlab ="Log10 Transformed TSS Concentrations", ylab="Log10 Transformed THg Concentrations",pch=21, bg="black")

x<-c(AllTSS_CWS, rev(AllTSS_CWS))
y<-c(pred_interval_lm_full[,2], rev(pred_interval_lm_full[,3]))

polygon(x,y, col="grey82", border=NA)
abline(lm_CWS_full, col="black", lwd=3)
points(AllTSS_CWS, AllTHG_CWS, pch=21,bg="black",cex=1.3)
zm()

#Prediction Interval Plot for LMEM 7 (Random factor-Permittes only) using jtools package (Figure 5)

pred_LMEM_CWS_7<-effect_plot(mixed.lmer7_CWS, pred=AllTSS_CWS, interval=T, plot.points = TRUE, level=0.95,x.label = "Log10 Transformed TSS Concentrations", y.label= "Log10 Transformed THg Concentrations")

pred_LMEM_CWS_7



#Prediction Interval Plot for LMEM 6 (Random factor-Nested Sites) using jtools package (Figure 6)

pred_LMEM_CWS_6<-effect_plot(mixed.lmer6_CWS, pred=AllTSS_CWS,interval=T, plot.points = TRUE, int.type = "prediction", x.label="Log10 Transformed TSS Concentrations", y.label="Log10 Transformed THg Concentrations")

pred_LMEM_CWS_6



#Alternative Application of LME Model-prediction intervals for each of the indivdual permittees (Table 8)----
#Prediction interval of permittee data


predictInterval(mixed.lmer7_CWS, X3Sites_CWS, include.resid.var=0, fix.intercept.variance = TRUE)



#Extract out the coefficients and the random effect intercepts for each random effect group to do alt.application----
#See DOcumentation for example calculation of each permittee's prediction interval

ranef(mixed.lmer7_CWS)
fixef(mixed.lmer7_CWS)
coef(mixed.lmer7_CWS)
summary(mixed.lmer7_CWS)
VarCorr(mixed.lmer7_CWS)
attributes(mixed.lmer7_CWS)
resid(mixed.lmer7_CWS)


#Calculate THg concentration values based on the lmer permittees-only equation----

R=100

Log10(Q) = 0.15970*log10(R)-5.30892



#Calculate % reduction of THg corresponding to % reduction in TSS (using permittee-only lmer equation)----
#From 100 mg/L

Z<-7.092184*10^(-6)

(1.024425*10^(-5)-Z)/(1.024425*10^(-5))*100

#From 60mg/L

H<-8.849709*10^(-6)

(9.441717*10^(-6)-H)/(1.024425*10^(-5))*100



#CDF-Cumulative Distribution Plots (Figure 10)---- 


#CDF Plot with lines (Figure 10)

plot(ecdf(raw_TSS_CWS), xlab="TSS Concentrations (mg/L)", main="", ylab="Cumulative Distribution Fraction" )



#Citations for R-for report----


citation()


sink("X3Sites_CWS_updated.txt")

print(citation())

#citation for packages

if(nchar(system.file(package="lattice"))) citation("ggplot2")
if(nchar(system.file(package="lattice"))) citation("lme4")
if(nchar(system.file(package="lattice"))) citation("MuMIn")
if(nchar(system.file(package="lattice"))) citation("gridExtra")
if(nchar(system.file(package="lattice"))) citation("nlme")
if(nchar(system.file(package="lattice"))) citation("jtools")
if(nchar(system.file(package="lattice"))) citation("predictmeans")
if(nchar(system.file(package="lattice"))) citation("emmeans")














