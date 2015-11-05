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
b <- cbind(total_pop,c(rep(0,length(sen_pop)),rep(1,length(res_pop))))
#roc(b[,2],b[,1])
#roc(b[,2],b[,1], smooth=TRUE)


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