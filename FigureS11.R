# Figure S11 -------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(readxl)
library(tidyverse)
library(data.table)
library(growthrates)

samples = read.table("growth_curve_data.txt",header = T, check.names = F,sep = "\t")

# Figure 3A --------------------------------------------------------------------
# take the mean of the OD600 for technical replicates - both strains and blanks
# for each individual experiment, subtract the mean OD600 of the strain minus the mean OD600 of the blank, for each sweat/sebum concentration
# then, take the mean and standard deviation of all biological replicates (across multiple dates) for each sweat/sebum concentration and time point
samples_avg = samples %>% group_by(Strain, Concentration, Time) %>% dplyr::summarise(avg_OD600 = mean(OD600,na.rm = T),#
                                                                                      avg_SD=sd(OD600,na.rm=T)) #

samples_avg = samples_avg %>% filter(Strain != "NA")
samples_avg$Strain = factor(samples_avg$Strain, levels =c("HC1",	"HC4",	"HC5",	"HC6",	"HC10",	"HC12",	"HC14",	"HC16",	"HC18",	"HC20",
                                                           "HC21",	"HC23",	"HC25",	"HC27",	"HC30",	"HC101",	"AD31",	"AD33",	"AD34",	"AD35",	
                                                           "AD36",	"AD38",	"AD39",	"AD40",	"AD42",	"AD46",	"AD50",	"AD51",	"AD52",	"AD53",	"AD59",	"AD60",	
                                                           "Acne61",	"Acne63",	"Acne64",	"Acne69",	"Acne74",	"Acne79",	"Acne81",	"Acne82",	"Acne83",	"Acne85",	
                                                           "Acne86",	"Acne88",	"Acne91",	"Acne92",	"Acne93",	"Acne95"
))

#head(samples_avg)

# color palette
rainbow = c("#C92C31","#F16704","#FFB01F","#DECB36","#BCE64C","#63BD71","#0A9396","#005F73")

# Plot sebum growth curves
sebum_big = samples_avg %>%
  ggplot(aes(x=Time,y=avg_OD600,color=as.factor(Concentration))) +  
  theme_bw() +
  facet_wrap(~Strain) + 
  geom_ribbon(aes(x=Time,ymin=avg_OD600-avg_SD, ymax=avg_OD600+avg_SD, fill=as.factor(Concentration)),alpha=0.08, color=NA) + 
  geom_smooth(se = FALSE,aes(linetype=as.factor(Concentration)), size=0.5) +
  scale_color_manual(values=rainbow) + 
  labs(color="Sebum\nconcentration (%)",linetype="Sebum\nconcentration (%)") + 
  xlab("Time (hours)") + ylab("OD600") + scale_linetype_manual(values=c("0"=1,"0.0075"=2,"0.01"=1,"0.025"=2,"0.05"=1,"0.075"=2,"0.1"=1,"0.25"=2)) + 
  scale_fill_manual(values=rainbow, guide="none")+
  ylim(-0.3,1.3)

sebum_big

# Figure 3B --------------------------------------------------------------------
# Area Under the Curve
conc = c("0","0.0075","0.01","0.025","0.05","0.075","0.1","0.25")
strain = unique(DataAllStrains$Strain)
strain = strain[which(strain != "NA")]
strain = strain[which(strain != "Blank")]

head(samples_avg)
# Calculate area under the curve for each sebum concentration for each strain
AUC = lapply(strain, FUN=function(X){
  strain_filter = samples_avg %>% filter(Strain == X)
  areas = lapply(conc, FUN=function(Y){
    conc_filter = strain_filter %>% filter(Concentration == Y)
    fit = loess(conc_filter$avg_OD600~conc_filter$Time) # Fit LOESS model
    df = data.frame(x = conc_filter$Time)
    df = transform(df, y.pred = predict(fit, df))
    auc = pracma::trapz(df$x, df$y.pred)
  })
  names(areas) = conc
  return(areas)
})

head(AUC)
# combine all data for each strain and sebum concentration
aucauc = rbindlist(AUC)
aucauc$Strain = strain
head(aucauc)
auc_df = melt(aucauc,id.vars = "Strain")
colnames(auc_df) = c("Strain","Concentration","AUC")

#
auc_df$Strain = factor(auc_df$Strain, levels = c("HC1",	"HC4",	"HC5",	"HC6",	"HC10",	"HC12",	"HC14",	"HC16",	"HC18",	"HC20",
                                                  "HC21",	"HC23",	"HC25",	"HC27",	"HC30",	"HC101",	"AD31",	"AD33",	"AD34",	"AD35",	
                                                  "AD36",	"AD38",	"AD39",	"AD40",	"AD42",	"AD46",	"AD50",	"AD51",	"AD52",	"AD53",	"AD59",	"AD60",	
                                                  "Acne61",	"Acne63",	"Acne64",	"Acne69",	"Acne74",	"Acne79",	"Acne81",	"Acne82",	"Acne83",	"Acne85",	
                                                  "Acne86",	"Acne88",	"Acne91",	"Acne92",	"Acne93",	"Acne95"
))

head(auc_df)

# Plot sebum AUC
auc_plot = ggplot(auc_df, aes(x=Concentration, y=AUC, fill=Concentration)) + 
  theme_bw() + 
  geom_col(color="black") + 
  facet_wrap(~Strain) + 
  #facet_grid(~Strain) +
  scale_fill_manual(values = rainbow) + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none",
        text=element_text(size=8)) + xlab("Sebum concentration (%)") + 
  ylim(0,100)

auc_plot

# Linear Regression
# normalize AUC data out of 1
auc_df_min_max = auc_df %>% group_by(Strain) %>% summarise(minAUC = min(AUC), maxAUC = max(AUC))
auc_df = left_join(auc_df, auc_df_min_max)
auc_df$NormAUC = (auc_df$AUC - auc_df$minAUC) / (auc_df$maxAUC - auc_df$minAUC)

# linear regression with AUC data for each strain
fitted_models = auc_df %>% 
  group_by(Strain) %>% 
  do(model=lm(NormAUC~as.numeric(Concentration), data = .))

# extract p-values from fitted models
lmp = function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f = summary(modelobject)$fstatistic
  p = pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
lmp
# extract coefficients from fitted models
coef = function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f = summary(modelobject)$coefficients
  c = f[2]
  return(c)
}
coef
# create a data frame with extracted coefficients and p-values for each strain
auc_lm = apply(fitted_models, 1, FUN=function(X) {
  strain = X[[1]]
  coef = coef(X[[2]])
  pval = lmp(X[[2]])
  return(data.frame(Strain = strain, Coef = coef, Pval= pval))
})
auc_lm

auc_lm = rbindlist(auc_lm)
head(auc_lm)
write.csv(auc_lm,"auc_lm.csv");# add group information

# Figure 3C --------------------------------------------------------------------
auc_lm = read.csv("auc_lm.csv",row.names = 1,header = T) # with group information

#
auc_lm$Group = factor(auc_lm$Group, levels = c("HC","AD","Acne"))

# color palette
cols = c("#90A6BB","#FFD579","#9D4E3F")

# plotA
auc_point = ggplot(auc_lm, aes(x=Coef, y=-log10(Pval), color = Group)) + 
  theme_classic() + 
  geom_point(size=6,alpha=.8) +
  scale_color_manual(values=cols) +
  xlim(-0.3, 0.3) + 
  ylab("-log10(p-value)") +
  xlab("Slope coefficient") + 
  geom_hline(yintercept = -log10(0.05), color = "grey") + 
  geom_vline(xintercept = 0, color = "grey") +
  annotate("text", x=-0.16, y=4.5, label = "Low sebum") +
  annotate("text", x=0.16, y=4.5, label = "High sebum") + 
  annotate("text",x=-0.2, y=0.5, label = "No preference") + 
  theme(legend.position = "bottom")

auc_point

# classify sebum preferences
auc_lm$SebumPref = ifelse(auc_lm$Pval < 0.05 & auc_lm$Coef < 0, "Low sebum",
                           ifelse(auc_lm$Coef > 0 & auc_lm$Pval < 0.05, "High sebum", "No preference"))
auc_lm$SebumPref = factor(auc_lm$SebumPref, levels = c("No preference","High sebum","Low sebum"))

# plotB
pref_auc = ggplot(auc_lm, aes(x=SebumPref, fill=Group)) + geom_bar(color="black") + 
  theme_bw() +
  scale_fill_manual(values=cols) +
  scale_x_discrete(drop=FALSE) + 
  xlab("Sebum preference predicted by AUC") + 
  ylab("Number of strains") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.major.x = element_blank() , 
        panel.grid.major.y = element_line(size=0.1, color="black",) ,
        panel.grid.minor.y = element_blank()) + 
  scale_y_continuous(breaks = seq(0,50, by=5)) + 
  theme(legend.position = "none")
pref_auc





