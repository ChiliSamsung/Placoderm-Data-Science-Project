

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Necessary files:
##A nexus file with all taxa
##A matrix with headers Taxon, Fauna, Lower, Upper, BA0, BA1, BA2, BA3, BA4, BA5, BA6. The latter are probabilities out of 1 for use with AncThresh. Max (Lower) and Min (Upper) ages cannot be the same, and cannot contain decimal points.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#install.packages("phytools")
#install.packages("ape")

library (phytools)
library(ape)
source("http://www.graemetlloyd.com/pubdata/functions_2.r")
#source("functions_2.R") 

##Load tree
tree <- read.nexus ("supplemented_supertree_new1.nex")
plot(tree)

##Load data
data <- read.table ("Placoderm_PCM_new.txt", header = TRUE)

##Get parameters
dimensions <- dim (data)
number.taxa <- dimensions [1]

##vectors used repeatedly in looped functions
fauna.names <- levels (data [, "Fauna"])

fauna.number <- length (fauna.names)

fauna.ages <- array (dim = number.taxa)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#STEP 1: Assign ages to taxa.  This is done by assuming ages of horizons 
#are uniformly distributed within the upper and lower age brackets. 
#All members of a given fauna are assigned the same age.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in 1:fauna.number) #loop through the faunas, and assign ages
  
{
  
  fauna.subset <- subset (data, Fauna == fauna.names[i])
  
  positions <- which (data [, "Fauna"] == fauna.names [i])
  
  new.age <- runif (1, min = fauna.subset [1, "Upper"], max = fauna.subset [1, "Lower"]) #Note that "Upper" = min and "Lower" = max due to syntax of function. Having decimal points produces NAs
  
  fauna.ages [positions] <- new.age 
  
  
}



tip.ages <- fauna.ages #assigns new ages to the a vector of tips

names (tip.ages) <- data [ , 1]

tip.ages

write.table(tip.ages, file="tipages.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#STEP 2: Timescale Tree
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ages<-read.table("tipages.txt", row.names=1)
ages

#source("functions_2.R") 
#options(error = recover)
scaled.tree <- date.phylo(tree, ages, rlen=3, method="equal")

write.nexus(scaled.tree, file="scaledtree.nex")

lengths<-scaled.tree$edge.lengths

write.table(lengths, file="branchlengths.txt")

plot(scaled.tree, edge.width=2, cex=0.25, label.offset=2); nodelabels(bg="white", cex=0.5)
plot(scaled.tree, edge.width=2, cex=0.25, label.offset=2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#STEP 2: Use AncThresh on scaled.tree
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##Subset Data for use with AncThresh
AT.data <- data [ , c("Taxon", "BA0", "BA1", "BA2", "BA3", "BA4", "BA5", "BA6")]
rownames(AT.data)<-AT.data[,1]
AT.data <-AT.data[,-1]
AT.data
scaled.tree <- read.nexus ("scaledtree.nex")

##Run ancThresh on BM
ancestralBM <- ancThresh(scaled.tree, AT.data, ngen=20000000, model="BM")

##Get results
write.csv(ancestralBM$par, file="ancThreshBMpar.csv")
write.csv(ancestralBM$liab, file="ancThreshBMliab.csv")
write.csv(ancestralBM$mcmc, file="ancThreshBMmcmc.csv")
write.csv(ancestralBM$ace, file="ancThreshBMace.csv")


##Run ancThresh on OU
ancestralOU <- ancThresh(scaled.tree, AT.data, ngen=20000000, model="OU")

##Get OU results
write.csv(ancestralOU$par, file="ancThreshOUpar.csv")
write.csv(ancestralOU$liab, file="ancThreshOUliab.csv")
write.csv(ancestralOU$mcmc, file="ancThreshOUmcmc.csv")
write.csv(ancestralOU$ace, file="ancThreshOUace.csv")

##Run ancThresh on Lambda
ancestralLambda <- ancThresh(scaled.tree, AT.data, ngen=20000000, model="lambda")

##Get Lambda results
write.csv(ancestralLambda$par, file="ancThreshLambdapar.csv")
write.csv(ancestralLambda$liab, file="ancThreshLambdaliab.csv")
write.csv(ancestralLambda$mcmc, file="ancThreshLambdamcmc.csv")
write.csv(ancestralLambda$ace, file="ancThreshLambdaace.csv")


##Get DIC scores
BMThreshDIC <-threshDIC(scaled.tree, AT.data, ancestralBM, burnin=200000)
LambdaThreshDIC <-threshDIC(scaled.tree, AT.data, ancestralLambda, burnin=200000)
OUThreshDIC <-threshDIC(scaled.tree, AT.data, ancestralOU, burnin=200000)

#Save results in Matrix
scores<-rbind(BMThreshDIC, OUThreshDIC, LambdaThreshDIC)
row.names(scores)<-c("BM", "OU", "Lambda")
write.csv(scores, file="ModelcomparisonDIC.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#STEP 4: Plot Thresholds and Parameters for Best Model ##Enter manually for now
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#for this run OU had the lowest DIC

#Load parameters

##debugging code
par <- read.csv(file="ancThreshOUpar.csv")
#NEED to chop off the 20% first
BA0 <- par[202:1001, "BA0"]
BA1 <- par[202:1001, "BA1"]
BA2 <- par[202:1001, "BA2"]
BA3 <- par[202:1001, "BA3"]
BA4 <- par[202:1001, "BA4"]
BA5 <- par[202:1001, "BA5"]
BA6 <- par[202:1001, "BA6"]

mean(BA0)
mean(BA1)
mean(BA2)
mean(BA3)
mean(BA4)
mean(BA5)
mean(BA6)
###

par<-ancestralOU$par
par = read.csv("ancThreshOUpar.csv")

#breaks
breaks<-c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)

#Get histograms of values for OU
breaks <- c(0, 1, 2, 3, 4, 5)
BA1<-hist(par[202:1001,"BA1"], breaks=breaks, ylim = c(0, 1), freq = FALSE) 

breaks <- c(0, 1, 2, 3, 4, 5)
BA2<-hist(par[202:1001,"BA2"], breaks=breaks, ylim = c(0, 1), freq = FALSE) 

breaks<-c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100)
BA3<-hist(par[202:1001,"BA3"], breaks=breaks, freq = FALSE, ylim = (c(0, 0.05)))

breaks<-c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160)
BA4<-hist(par[202:1001,"BA4"], breaks=breaks)

breaks<-c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)
BA5<-hist(par[202:1001,"BA5"], breaks=breaks, freq = FALSE)

#Plot
plot(BA1$mids,BA1$density,type="s",ylim=c(0, max(c(BA1$density,BA2$density))),xlab="liability", ylab="density",col="red")
lines(BA1$mids,BA2$density,type="s",col="green")
lines(BA1$mids,BA3$density,type="s",col="blue")
lines(BA1$mids,BA4$density,type="s",col="orange")
lines(BA1$mids,BA5$density,type="s",col="purple")

#Plot Alpha Distribution
burnin<-200000
ii<-which(par[,"gen"]==burnin)
ps.alpha<-par[ii:nrow(par),"alpha"]
mean(ps.alpha)
pd<-density(ps.alpha,bw=0.4)
plot(pd,xlab="alpha",main="posterior density of alpha", xlim = c(-1,1))
