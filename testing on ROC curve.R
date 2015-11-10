#test for cumulative distribution and normal calculation

a <- rnorm(120, 3, 1.45)
sum(a>5)/length(a)
pnorm(5,3,1.45, lower.tail = FALSE)


b <- cbind(total_pop,total_pop>5)
roc(b[,2],b[,1])


#roc example
data(aSAH)
# Build a ROC object and compute the AUC
roc(aSAH$outcome, aSAH$s100b)
roc(outcome ~ s100b, aSAH)
# Smooth ROC curve
roc(outcome ~ s100b, aSAH, smooth=TRUE)
# more options, CI and plotting
roc1 <- roc(aSAH$outcome,
            aSAH$s100b, percent=TRUE,
            # arguments for auc
            partial.auc=c(100, 90), partial.auc.correct=TRUE,
            partial.auc.focus="sens",
            # arguments for ci
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            # arguments for plot
            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
            print.auc=TRUE, show.thres=TRUE)
# Add to an existing plot. Beware of 'percent' specification!
roc2 <- roc(aSAH$outcome, aSAH$wfns,
            plot=TRUE, add=TRUE, percent=roc1$percent)


#testing ROC
nn <- 200
pr <- .1
cutoff <- 3.5

sen_pop <- rlnorm(nn*(1-pr),log(3),log(1.45))
res_pop <- rlnorm(nn*pr, log(4),log(1.22))
total_pop <- c(sen_pop,res_pop)
#b <- cbind(total_pop,c(rep(0,length(sen_pop)),rep(1,length(res_pop))))
b <- total_pop
#going back to base plotting system
hist(b, probability = TRUE, col="grey",lwd=2,ps=20,breaks=as.numeric(floor(min(b)):ceiling(max(b))), main="Histogram of Simulated Half-Lives", xlab="Half-life (hours)")
hist(res_pop, probability = TRUE, col="red", add=T, breaks=as.numeric(floor(min(b)):ceiling(max(b))))
lines(density(res_pop),lwd=5, col="red")
abline(v=cutoff, lwd=3, col="blue")


#roc(b[,2],b[,1])
#roc(b[,2],b[,1], smooth=TRUE)
popDF2 <- b

popDF2[popDF2[,2]==0,2] <- "Sensitive"
popDF2[popDF2[,2]==1,2] <- "Resistant"

popDF2 <- as.data.frame(popDF2)
popDF2[,1] <- as.numeric(as.character(popDF2[,1]))
names(popDF2) <- c("Half-life (hours)","Sensitivity")

ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour=Sensitivity)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks=as.numeric(floor(min(popDF2[,1])):ceiling(max(popDF2[,1]))))

ggplot(popDF2, aes(x=`Half-life (hours)`, colour=Sensitivity)) +
  geom_density() +
  scale_x_continuous(breaks=as.numeric(floor(min(popDF2[,1])):ceiling(max(popDF2[,1]))))



#NOT smoothed density
hist(b[,1], freq=FALSE)
lines(density(b[,1]))
lines(density(sen_pop))
lines(density(res_pop),col="red")



roc3 <- roc(b[,2], b[,1], percent=TRUE, partial.auc=c(100, 90), partial.auc.correct=TRUE, partial.auc.focus="sens",ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE, show.thres=TRUE)

roc3 <- roc(b[,2], b[,1],  partial.auc.correct=TRUE, partial.auc.focus="sens",ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE, show.thres=TRUE)


#pointing the TPR and FPR
FPR <- plnorm(cutoff, log(3), log(1.45), lower.tail=FALSE)
TPR <- plnorm(cutoff, log(4), log(1.22), lower.tail=FALSE)

TPR2 <- sum(res_pop>cutoff)/length(res_pop)
FPR2 <- sum(sen_pop>cutoff)/length(sen_pop)

points((1-FPR),TPR, col="red")


points((1-FPR2),TPR2)


#stacked histogram in R
library(ggplot2)
#ggplot doesn't deal with matrix so,
c <- as.data.frame(b)
c <- cbind(c, c(rep("sen",length(sen_pop)),rep("res",length(res_pop))))
names(c)[3] <- "lab"
ggplot(c, aes(x=total_pop, fill=lab),binwidth = 1) +
  geom_histogram()

library(mixtools)
mxmdl = normalmixEM(b)
plot(mxmdl, which=2)
lines(density(b, lty=2, lwd=2))

#example of normalized histogram from web
# Fake data (two normal distributions)
set.seed(20)
dat1 = data.frame(x=rnorm(1000, 100, 10), group="A")
dat2 = data.frame(x=rnorm(200, 120, 20), group="B")
dat = rbind(dat1, dat2)

ggplot(dat, aes(x, fill=group, colour=group)) +
  geom_histogram(breaks=seq(0,200,5), alpha=0.6, 
                 position="identity", lwd=0.2) +
  ggtitle("Unormalized")

ggplot(dat, aes(x, fill=group, colour=group)) +
  geom_histogram(aes(y=2*(..density..)/sum(..density..)), breaks=seq(0,200,5), alpha=0.6, 
                 position="identity", lwd=0.2) +
  ggtitle("Normalized")

ggplot(dat, aes(x, fill=group, colour=group)) +
  geom_histogram(aes(y=(..count..)/sum(..count..)), breaks=seq(0,200,5), alpha=0.6, 
                 position="identity", lwd=0.2) +
  ggtitle("Normalized")
