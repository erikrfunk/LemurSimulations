# functions ----
harmonic_mean <- function(x) {
  1 / mean(1/x, na.rm = TRUE)
}
smooth_ribbon_mid95 = function(df,col) {
  dat = df[,c("Generation","Rep","Set",col)]
  res_names = c("Generation","Set",col)
  simMin = setNames(aggregate(get(col)~Generation+Set,data=dat,FUN = function(z) quantile(z,0.025)),res_names)
  simMinSmoothed = NULL
  for(i in 1:length(unique(dat$Set))){
    setMins = simMin[simMin$Set == unique(simMin$Set)[i],]
    loess_fit <- loess(setMins[,3] ~ setMins$Generation, span = 0.5)
    smoothed_min <- predict(loess_fit)
    smoothed_res = data.frame("Generation"=seq(1:max(dat$Generation)),
                              "Set" = unique(simMin$Set)[i],
                              "SmoothedMin" = smoothed_min)
    simMinSmoothed = rbind(simMinSmoothed,smoothed_res)
  }
  simMax = setNames(aggregate(get(col)~Generation+Set,data=dat,FUN = function(z) quantile(z,0.975)),res_names)
  simMaxSmoothed = NULL
  for(i in 1:length(unique(dat$Set))){
    setMaxes = simMax[simMax$Set == unique(simMax$Set)[i],]
    loess_fit <- loess(setMaxes[,3] ~ setMaxes$Generation,span = 0.5)
    smoothed_max <- predict(loess_fit)
    smoothed_res = data.frame("Generation"=seq(1:max(dat$Generation)),
                              "Set" = unique(simMax$Set)[i],
                              "SmoothedMax" = smoothed_max)
    simMaxSmoothed = rbind(simMaxSmoothed,smoothed_res)
  }
  ribbon = merge(simMinSmoothed,simMaxSmoothed,by=c("Generation","Set"))
  #names(ribbon) = c("Generation","Set","SmoothedMin","SmoothedMax")
  ribbon$Set = as.factor(ribbon$Set)
  return(ribbon)
}
smooth_ribbon_95CI = function(df,col) {
  # Calculate and smooth a 95% CI 
  dat = df[,c("Generation","Rep","Set",col)]
  n = max(dat$Rep)
  res_names = c("Generation","Set",col)
  simMin = setNames(aggregate(get(col)~Generation+Set,data=dat,
                              FUN = function(x) mean(x) - (1.96*(sd(x)/sqrt(n)))),res_names)
  simMinSmoothed = NULL
  for(i in 1:length(unique(dat$Set))){
    setMins = simMin[simMin$Set == unique(simMin$Set)[i],]
    loess_fit <- loess(setMins[,3] ~ setMins$Generation, span = 0.5)
    smoothed_min <- predict(loess_fit)
    smoothed_res = data.frame("Generation"=seq(1:max(dat$Generation)),
                              "Set" = unique(simMin$Set)[i],
                              "SmoothedMin" = smoothed_min)
    simMinSmoothed = rbind(simMinSmoothed,smoothed_res)
  }
  simMax = setNames(aggregate(get(col)~Generation+Set,data=dat,
                              FUN = function(x) mean(x) + (1.96*(sd(x)/sqrt(n)))),res_names)
  simMaxSmoothed = NULL
  for(i in 1:length(unique(dat$Set))){
    setMaxes = simMax[simMax$Set == unique(simMax$Set)[i],]
    loess_fit <- loess(setMaxes[,3] ~ setMaxes$Generation,span = 0.5)
    smoothed_max <- predict(loess_fit)
    smoothed_res = data.frame("Generation"=seq(1:max(dat$Generation)),
                              "Set" = unique(simMax$Set)[i],
                              "SmoothedMax" = smoothed_max)
    simMaxSmoothed = rbind(simMaxSmoothed,smoothed_res)
  }
  ribbon = merge(simMinSmoothed,simMaxSmoothed,by=c("Generation","Set"))
  #names(ribbon) = c("Generation","Set","SmoothedMin","SmoothedMax")
  ribbon$Set = as.factor(ribbon$Set)
  return(ribbon)
}


dir = "/home/centos/USS/erik/lemurs/Vrubra_Andranobe_simulations/super11_N100/"
dir = "/home/centos/USS/erik/lemurs/Vava_Ranomafana_simulations/super11_N150/"
fns = list.files(dir,pattern = "_stats.txt",full.names=T,recursive=T)
#fns = fns[grep("old",fns,invert=T)]

# Can run this for multiple dirs if you don't reset the res df in the next line
res = NULL
for(i in 1:length(fns)){
  rep = try(read.table(fns[i],sep=",",skip=22,header=T),silent = T)
  if(!(inherits(rep, "try-error"))){
    cat(paste0("Adding: ",fns[i],"\n"))
  rep$Rep = i
  rep$Set = basename(dir)
  res = rbind(res,rep)
  }
}
res$Set = as.factor(res$Set)
full_res = res # Can return here to retrieve an untrimmed version
res = full_res

# Plot Stats
# For Froh, check how much burn in is needed to equal the modern empirical estimate (0.02)
aggregate(Froh~Generation,data=res,FUN=mean)
res$Generation = res$Generation - 28
res = res[res$Generation > 0 & res$Generation <= 150,]
res_ribbon_mid95 = smooth_ribbon_mid95(res,"Froh")
res_ribbon_95CI = smooth_ribbon_95CI(res,"Froh")
ggplot(res,aes(x=Generation,y=Froh,group=Set,color=Set))+
  geom_ribbon(data = res_ribbon_95CI, inherit.aes = F, aes(x=Generation, ymin = SmoothedMin,ymax = SmoothedMax,fill=Set),alpha=0.3) +
  geom_ribbon(data = res_ribbon_mid95, inherit.aes = F, aes(x=Generation, ymin = SmoothedMin,ymax = SmoothedMax,fill=Set),alpha=0.15) +
  geom_smooth(span=1,alpha=0.2,se=F)+
  scale_color_manual(values=c("chocolate3","darkcyan"))+
  scale_fill_manual(values=c("chocolate3","darkcyan"))+
  labs(title = "Vrubra Froh",x="Years",y="Froh",color="Population Size")+
  theme_bw()

# Het reference lines
starting_het = 0.00156
starting_het = 0.000858
het_loss_gen = 50
p05 = (((1 - 0.05)*starting_het) - starting_het)/het_loss_gen
p10 = (((1 - 0.1)*starting_het) - starting_het)/het_loss_gen
p20 = (((1 - 0.2)*starting_het) - starting_het)/het_loss_gen

res_ribbon_95CI = smooth_ribbon_95CI(res,"ExpHet")
res_ribbon_mid95 = smooth_ribbon_mid95(res,"ExpHet")
ggplot(res,aes(x=Generation,y=ExpHet,group=Set,color=Set))+
  geom_ribbon(data = res_ribbon_95CI, inherit.aes = F, aes(x=Generation, ymin = SmoothedMin,ymax = SmoothedMax,fill=Set),alpha=0.3) +
  geom_ribbon(data = res_ribbon_mid95, inherit.aes = F, aes(x=Generation, ymin = SmoothedMin,ymax = SmoothedMax,fill=Set),alpha=0.15) +
  geom_smooth(span=1,alpha=0.2,se=F)+
  geom_abline(intercept=starting_het,slope=p05,alpha=0.5,linetype="dashed")+
  geom_abline(intercept=starting_het,slope=p10,alpha=0.5,linetype="dashed")+
  geom_abline(intercept=starting_het,slope=p20,alpha=0.5,linetype="dashed")+
  scale_color_manual(values=c("chocolate3","darkcyan"))+
  scale_fill_manual(values=c("chocolate3","darkcyan"))+
  labs(title = "Vrubra Het",x="Years",y="ExpHet",color="Population Size")+
  theme_bw()


#Statistics at a set year
year=50
#Froh - don't forget about burn in for this stat
avg = mean(res$Froh[res$Generation==year])
se = sd(res$Froh[res$Generation==year])/(10)
avg - 1.96*se
avg + 1.96*se

#Het as a percentage lost
starting_het = mean(res$ExpHet[res$Generation==1])
avg = mean(res$ExpHet[res$Generation==year])
se = sd(res$ExpHet[res$Generation==year])/(10) # divide by sqrt(n) = sqrt(100)
1-(avg/starting_het)
1-(avg + 1.96*se)/starting_het
1-(avg - 1.96*se)/starting_het

# Various Ne plots
# Read in the relevant gone results if plotting Ne
ne1 = read.table("/home/centos/USS/erik/lemurs/Vava_Ranomafana_simulations/super11_N150/AllGoneResults.txt", header=T)
ne1_harmean = aggregate(Ne_diploids~Iterations,data=ne1[ne1$Generation < 18 & ne1$Generation > 5,], FUN=harmonic_mean)
names(ne1_harmean)[1:2] = c("Rep","Size")
ne1_harmean$Set = "super11_N100"
ne1_harmean$SizeType = "NeGONE2"
ne2 = read.table("/home/centos/USS/erik/lemurs/Vava_Ranomafana_simulations/super11_N200/addedMortality/AllGoneResults.txt", header=T)
ne2_harmean = aggregate(Ne_diploids~Iterations,data=ne2[ne2$Generation < 18 & ne2$Generation > 5,], FUN=harmonic_mean)
names(ne2_harmean)[1:2] = c("Rep","Size")
ne2_harmean$Set = "super11_N150"
ne2_harmean$SizeType = "NeGONE2"
nes_gone2 = rbind(ne1_harmean,ne2_harmean)

ne_har_means = aggregate(NeTemporalAF~Rep,data=res[res$Generation>30,],FUN=harmonic_mean)
ggplot(res,aes(x=Generation,linetype=Set))+
  geom_smooth(span=1,alpha=0.25,aes(y=PopSize,color="census"))+
  geom_smooth(span=1,alpha=0.25,aes(y=NeTemporalAF,color="temporal"))+
  #geom_hline(yintercept=mean(ne_har_means$NeTemporalAF),linetype="dashed")+
  geom_smooth(span=1,alpha=0.25,aes(y=NeBreeders,color="breeders"))+
  scale_color_manual(labels = c(census = "Census", temporal="Ne: Temporal AF",breeders="Ne: Breeders"),
                     values = c(census = "#007193",temporal="#486C0D", breeders="#B59410"),
                     limits = c("census","temporal","breeders"))+
  labs(x="Years",y="Ne",color="Population Size")+
  theme_bw()

# A Boxplot visualization of Ne: First generate harmonic means for all
harmean_census = aggregate(PopSize~Rep+Set,data=res[res$Generation>30,],FUN=harmonic_mean)
names(harmean_census)[3] = "Size"
harmean_census$SizeType = "Census"
harmean_breeders = aggregate(NeBreeders~Rep+Set,data=res[res$Generation>30,],FUN=harmonic_mean)
names(harmean_breeders)[3] = "Size"
harmean_breeders$SizeType = "NeBreeders"
harmean_temporal = aggregate(NeTemporalAF~Rep+Set,data=res[res$Generation>30,],FUN=harmonic_mean)
names(harmean_temporal)[3] = "Size"
harmean_temporal$SizeType = "NeTemporal"

all_harmeans = rbind(harmean_census,harmean_breeders,harmean_temporal,nes_gone2[,c("Rep","Set","Size","SizeType")])
ggplot(all_harmeans,aes(x=SizeType,y=Size,fill=Set))+
  geom_boxplot()+
  scale_fill_manual(values = c("chocolate3","darkcyan"))+
  theme_bw()



