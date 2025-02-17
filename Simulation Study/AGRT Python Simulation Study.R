# To install now archived "mdsdt" package:
install.packages(c("ellipse","polycor"))
install.packages("https://cran.r-project.org/src/contrib/Archive/mdsdt/mdsdt_1.2.tar.gz", type="source", repos=NULL)

# Load required packages
require(mnormt)
require(stringr)
require(mdsdt)
require(scales)

# Generalized subject model
sm.full <- function (x,y,delta,epsilon,sig,m,n,rho, R=NULL) {
  # m := 2 element vector specifying the mean spreading for each dimension
  # n := 2 element vector specifying the var spreading for each dimension
  # rho := 4 element vector specifying the correlation for each quadrant
  # sig := 2 element vector specifying the standard deviation for each dimension
  # delta := 2 element vector specifying the decision bound intercept for each dimension
  # epsilon := 2 element vector specifying the decision bound slopes for each dimension
  
  # calculate decision bounds based on where we are in stimulus space
  tau <- delta + epsilon * c(y,x)
  # create indicator for which quadrant we're in
  zeta <- c(if (x > tau[1]) 1 else -1, if (y > tau[2]) 1 else -1)
  # calculate mean vector
  mymean <- c(x,y) + rev(zeta) * m * c(x,y)
  # calculate sd vector
  mysd <- sig + rev(zeta) * n * sig
  # pick the rho that corresponds to our quadrant
  if (sum(zeta)==2) {
    myrho <- rho[4]
  } else if (sum(zeta)==-2) {
    myrho <- rho[1]
  } else if (zeta[1]>0) {
    myrho <- rho[3]
  } else {
    myrho <- rho[2]
  }
  # create the cov matrix
  mycov <- matrix(c(mysd[1]^2, myrho*mysd[1]*mysd[2], myrho*mysd[1]*mysd[2], mysd[2]^2), 2,2,T)
  c(
    sadmvn(lower=c(-Inf,-Inf), upper=tau, mean=mymean, varcov = mycov)[1],
    sadmvn(lower=c(-Inf,tau[2]), upper=c(tau[1],Inf), mean=mymean, varcov = mycov)[1],
    sadmvn(lower=c(tau[1],-Inf), upper=c(Inf,tau[2]), mean=mymean, varcov = mycov)[1],
    sadmvn(lower=tau, upper=c(Inf,Inf), mean=mymean, varcov = mycov)[1]
  )[if (is.null(R)) 1:4 else R]
}

# Subject model used to generate data
test.sm <- function (x,y,correct.resp.ind,sm) {
  switch(sm, 
         "Recovery" = .04 * .25 + (1-.04) * sm.full(x, y, c(600, 4), c(0, 0), c(100, .8), c(0, 0), c(0, 0), c(0, 0, 0, 0), R=correct.resp.ind+1),
         "PSM" = .04 * .25 + (1-.04) * sm.full(x, y, c(600, 4), c(0, 0), c(100, .8), c(.1, 0), c(0, 0), c(0, 0, 0, 0), R=correct.resp.ind+1),
         "PSV" = .04 * .25 + (1-.04) * sm.full(x, y, c(600, 4), c(0, 0), c(100, .8), c(0, 0), c(0, .4), c(0, 0, 0, 0), R=correct.resp.ind+1),
         "PSMV" = .04 * .25 + (1-.04) * sm.full(x, y, c(600, 4), c(0, 0), c(100, .8), c(.1, 0), c(0, .4), c(0, 0, 0, 0), R=correct.resp.ind+1),
         "PI" = .04 * .25 + (1-.04) * sm.full(x, y, c(600, 4), c(0, 0), c(100, .8), c(0, 0), c(0, 0), c(.9, 0, -.9, 0), R=correct.resp.ind+1),
         "PI.PSM" = .04 * .25 + (1-.04) * sm.full(x, y, c(600, 4), c(0, 0), c(100, .8), c(.1, 0), c(0, 0), c(.9, 0, -.9, 0), R=correct.resp.ind+1),
         "PI.PSV" = .04 * .25 + (1-.04) * sm.full(x, y, c(600, 4), c(0, 0), c(100, .8), c(0, 0), c(0, .4), c(.9, 0, -.9, 0), R=correct.resp.ind+1),
         "PI.PSMV" = .04 * .25 + (1-.04) * sm.full(x, y, c(600, 4), c(0, 0), c(100, .8), c(.1, 0), c(0, .4), c(.9, 0, -.9, 0), R=correct.resp.ind+1))
} 

# Import Data

# Without lapse
#mydata <- read.table("/Users/jjg367/Documents/Outside Research/Adaptive GRT/Psi Source/Test_Simulation/data/2023-09-15_21h12.50.080.txt", sep="\t", header=T) # 1 subject, 100 adaptive trials, 1000 GRT trials
#mydata <- read.table("/Users/jjg367/Documents/Outside Research/Adaptive GRT/Psi Source/Test_Simulation/data/2023-09-15_21h13.55.920.txt", sep="\t", header=T) # 100 subjects, 100 adaptive trials, 1000 GRT trials
#mydata <- read.table("/Users/jjg367/Documents/Outside Research/Adaptive GRT/Psi Source/Test_Simulation/data/2023-09-15_23h18.41.285.txt", sep="\t", header=T) # 100 subjects, 300 adaptive trials, 1000 GRT trials

# With lapse
#mydata <- read.table("/Users/jjg367/Documents/Outside Research/Adaptive GRT/Psi Source/Test_Simulation/data/2023-09-17_18h43.35.612.txt", sep="\t", header=T) # 1 subject, 300 adaptive trials, 1000 GRT trials
#mydata <- read.table("/Users/jjg367/Documents/Outside Research/Adaptive GRT/Psi Source/Test_Simulation/data/2023-09-17_18h52.47.370.txt", sep="\t", header=T) # 100 subjects, 300 adaptive trials, 1000 GRT trials
#mydata <- read.table("/Users/jjg367/Documents/Outside Research/Adaptive GRT/Psi Source/Test_Simulation/data/2023-10-02_17h52.55.027.txt", sep="\t", header=T) # 100 subjects, 300 adaptive trials, 1000 GRT trials

# All Subject Models
mydata <- read.table("/Users/jjg367/Documents/Outside Research/Adaptive GRT/Psi Source/Test_Simulation/data/All Models 300 Adaptive Trials 2023-10-02_18h31.13.453.txt", sep="\t", header=T) # 100 subjects, 300 adaptive trials, 1000 GRT trials


# Add Missing Columns
mydata$Dim1 <- sapply(mydata$Stimulus, function (x) as.numeric(unlist(strsplit(str_sub(x,2,-2), ","))[1]))
mydata$Dim2 <- sapply(mydata$Stimulus, function (x) as.numeric(unlist(strsplit(str_sub(x,2,-2), ","))[2]))
mydata$StimInd <- NA
for (sm in unique(mydata$Version)) {
  for (s in unique(mydata$Subject[mydata$Version==sm])) {
    # 00 01 10 11
    d1 <- range(mydata$Dim1[mydata$Block=="main" & mydata$Subject==s & mydata$Version==sm])
    d2 <- range(mydata$Dim2[mydata$Block=="main" & mydata$Subject==s & mydata$Version==sm])
    mydata$StimInd[mydata$Block=="main" & mydata$Subject==s & mydata$Version==sm] <- sapply(mydata$Stimulus[mydata$Block=="main" & mydata$Subject==s & mydata$Version==sm], function (x) {
      mystim <- eval(parse(text=paste0("c",x)))
      if (mystim[1] == d1[1] & mystim[2] == d2[1]) 0 else
        if (mystim[1] == d1[1] & mystim[2] == d2[2]) 1 else
          if (mystim[1] == d1[2] & mystim[2] == d2[1]) 2 else
            if (mystim[1] == d1[2] & mystim[2] == d2[2]) 3
    })
  } 
}

mydata$RespInd <- NA; mydata$RespInd[mydata$Response=="(0, 0)"] <- 0; mydata$RespInd[mydata$Response=="(0, 1)"] <- 1; mydata$RespInd[mydata$Response=="(1, 0)"] <- 2; mydata$RespInd[mydata$Response=="(1, 1)"] <- 3
mydata$Correct <- mydata$StimInd == mydata$RespInd
mydata$Alpha1 <- sapply(mydata$Lambda, function (x) as.numeric(unlist(strsplit(gsub("[()]","",x), ",")))[1])
mydata$Beta1 <- sapply(mydata$Lambda, function (x) as.numeric(unlist(strsplit(gsub("[()]","",x), ",")))[2])
mydata$Alpha2 <- sapply(mydata$Lambda, function (x) as.numeric(unlist(strsplit(gsub("[()]","",x), ",")))[3])
mydata$Beta2 <- sapply(mydata$Lambda, function (x) as.numeric(unlist(strsplit(gsub("[()]","",x), ",")))[4])
mydata$Pred.Dim1min <- mydata$Alpha1 - mydata$Beta1 * sqrt(2) * erfinv((2*sqrt(.75)-.04)/(1-.04)-1)
mydata$Pred.Dim1max <- mydata$Alpha1 - mydata$Beta1 * sqrt(2) * erfinv((2*(1-sqrt(.75))-.04)/(1-.04)-1)
mydata$Pred.Dim2min <- mydata$Alpha2 - mydata$Beta2 * sqrt(2) * erfinv((2*sqrt(.75)-.04)/(1-.04)-1)
mydata$Pred.Dim2max <- mydata$Alpha2 - mydata$Beta2 * sqrt(2) * erfinv((2*(1-sqrt(.75))-.04)/(1-.04)-1)
mydata$Pred.Acc <- NA
mydata$Pred.Acc[mydata$Block=="adapt"] <- mapply(function (x1,x2,y1,y2,z) {
  mean(c(test.sm(x1,y1,0,z), test.sm(x1,y2,1,z), test.sm(x2,y1,2,z), test.sm(x2,y2,3,z)))
}, x1=mydata$Pred.Dim1min[mydata$Block=="adapt"], x2=mydata$Pred.Dim1max[mydata$Block=="adapt"], 
y1=mydata$Pred.Dim2min[mydata$Block=="adapt"], y2=mydata$Pred.Dim2max[mydata$Block=="adapt"], z = mydata$Version[mydata$Block=="adapt"])

# For only one subject:
# par(mfrow=c(2,2))
# plot(mydata$Trial[mydata$Block=="adapt"], mydata$Alpha1[mydata$Block=="adapt"], type='l', main="Dimension 1", xlab="Trial", ylab="Alpha (Decision Threshold)"); abline(h=600,lty=2)
# plot(mydata$Trial[mydata$Block=="adapt"], mydata$Beta1[mydata$Block=="adapt"], type='l', main="Dimension 1", xlab="Trial", ylab="Beta (Stimulus Separation)"); abline(h=100,lty=2)
# plot(mydata$Trial[mydata$Block=="adapt"], mydata$Alpha2[mydata$Block=="adapt"], type='l', main="Dimension 2", xlab="Trial", ylab="Alpha (Decision Threshold)"); abline(h=4,lty=2)
# plot(mydata$Trial[mydata$Block=="adapt"], mydata$Beta2[mydata$Block=="adapt"], type='l', main="Dimension 2", xlab="Trial", ylab="Beta (Stimulus Separation)"); abline(h=.8,lty=2)
# plot(mydata$Trial[mydata$Block=="adapt"], mydata$Pred.Acc[mydata$Block=="adapt"], type='l', main="Predicted Accuracy", xlab="Trial", ylab="Overall Accuracy"); abline(h=.75, lty=2)
# par(mfrow=c(1,1))

# SINGLE SUBJECT MODEL, MULTIPLE SUBJECTS
# Plot Psi Parameters
# par(mfrow=c(2,2))
# plot(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Alpha1[mydata$Block=="adapt" & mydata$Trial==x])
#   }), type='l', main="Dimension 1", xlab="Trial", ylab="Alpha (Decision Threshold)", ylim=c(100,800)); abline(h=600,lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Alpha1[mydata$Block=="adapt" & mydata$Trial==x]) + sd(mydata$Alpha1[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Alpha1[mydata$Block=="adapt" & mydata$Trial==x]) - sd(mydata$Alpha1[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# plot(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Beta1[mydata$Block=="adapt" & mydata$Trial==x])
# }), type='l', main="Dimension 1", xlab="Trial", ylab="Beta (Stimulus Separation)", ylim=c(0,200)); abline(h=100,lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Beta1[mydata$Block=="adapt" & mydata$Trial==x]) + sd(mydata$Beta1[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Beta1[mydata$Block=="adapt" & mydata$Trial==x]) - sd(mydata$Beta1[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# plot(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Alpha2[mydata$Block=="adapt" & mydata$Trial==x])
# }), type='l', main="Dimension 2", xlab="Trial", ylab="Alpha (Decision Threshold)", ylim=c(1,12)); abline(h=4,lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Alpha2[mydata$Block=="adapt" & mydata$Trial==x]) + sd(mydata$Alpha2[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Alpha2[mydata$Block=="adapt" & mydata$Trial==x]) - sd(mydata$Alpha2[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# plot(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Beta2[mydata$Block=="adapt" & mydata$Trial==x])
# }), type='l', main="Dimension 2", xlab="Trial", ylab="Beta (Stimulus Separation)", ylim=c(0,2)); abline(h=.8,lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Beta2[mydata$Block=="adapt" & mydata$Trial==x]) + sd(mydata$Beta2[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Beta2[mydata$Block=="adapt" & mydata$Trial==x]) - sd(mydata$Beta2[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# par(mfrow=c(1,1))
# # Plot Predicted Accuracy
# plot(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#      mean(mydata$Pred.Acc[mydata$Block=="adapt" & mydata$Trial==x])
#   }), type='l', main="Predicted Accuracy", xlab="Trial", ylab="Overall Accuracy", ylim=c(0,1)); abline(h=.75, lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Pred.Acc[mydata$Block=="adapt" & mydata$Trial==x]) + sd(mydata$Pred.Acc[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# lines(unique(mydata$Trial[mydata$Block=="adapt"]), sapply(unique(mydata$Trial[mydata$Block=="adapt"]), function (x) {
#   mean(mydata$Pred.Acc[mydata$Block=="adapt" & mydata$Trial==x]) - sd(mydata$Pred.Acc[mydata$Block=="adapt" & mydata$Trial==x])
# }), lty=2)
# # Plot Actual Accuracy
# barplot(sapply(sort(unique(mydata$Subject)), function (x) mean(mydata$Correct[mydata$Block=="main" & mydata$Subject==x])),
#         main="Overall Accuracy", xlab="Simulated Subject", names.arg=sort(unique(mydata$Subject)), ylab="Accuracy", ylim=c(0,1)); abline(h=.75, lty=2)
# # Plot fitted GRT models for first 10 subjects
# for (s in 0:10) plot(fit.grt(table(mydata$StimInd[mydata$Block=="main" & mydata$Subject==s],mydata$RespInd[mydata$Block=="main" & mydata$Subject==s], exclude=NA)), main=s, marginals=T)


# ALL SUBJECT MODELS ALL SUBJECTS
# Use Shading
# x loc
# PSM MODELS IN RED
plot(0,0,type='n',xlim=c(0,300),ylim=c(400,800),xlab="Trial",ylab="Parameter Estimate",main="X Threshold")
abline(h=600,lty=1, lwd=2, xpd=F)
for (sm in unique(mydata$Version)) {
  conv.data <- subset(mydata, Block=="adapt" & Version==sm)
  lines(sort(unique(conv.data$Trial)), sapply(sort(unique(conv.data$Trial)), function (x) mean(conv.data$Alpha1[conv.data$Trial==x])), 
        col=if (sm %in% c("PSM","PSMV","PI.PSM","PI.PSMV")) 'red' else 'black')
  polygon(c(sort(unique(conv.data$Trial)),rev(sort(unique(conv.data$Trial)))), 
          c(sapply(sort(unique(conv.data$Trial)), 
                   function (x) mean(conv.data$Alpha1[conv.data$Trial==x]) + sd(conv.data$Alpha1[conv.data$Trial==x])),
            sapply(rev(sort(unique(conv.data$Trial))),
                   function (x) mean(conv.data$Alpha1[conv.data$Trial==x]) - sd(conv.data$Alpha1[conv.data$Trial==x]))),
          col=alpha(if (sm %in% c("PSM","PSMV","PI.PSM","PI.PSMV")) 'red' else 'black', alpha=.05), border=NA)
}
legend('top', legend=c("Violate Perceptual Separability", "Other"), fill=c('red','black'), horiz=T, bty='n', cex=.6)
# x slope
plot(0,0,type='n',xlim=c(0,300),ylim=c(0,200),xlab="Trial",ylab="Parameter Estimate",main="X SD")
abline(h=100,lty=1,lwd=2, xpd=F)
for (sm in unique(mydata$Version)) {
  conv.data <- subset(mydata, Block=="adapt" & Version==sm)
  lines(sort(unique(conv.data$Trial)), sapply(sort(unique(conv.data$Trial)), function (x) mean(conv.data$Beta1[conv.data$Trial==x])), 
        col=if (sm %in% c("PSM","PSMV","PI.PSM","PI.PSMV")) 'red' else 'black')
  polygon(c(sort(unique(conv.data$Trial)),rev(sort(unique(conv.data$Trial)))), 
          c(sapply(sort(unique(conv.data$Trial)), 
                   function (x) mean(conv.data$Beta1[conv.data$Trial==x]) + sd(conv.data$Beta1[conv.data$Trial==x])),
            sapply(rev(sort(unique(conv.data$Trial))),
                   function (x) mean(conv.data$Beta1[conv.data$Trial==x]) - sd(conv.data$Beta1[conv.data$Trial==x]))),
          col=alpha(if (sm %in% c("PSM","PSMV","PI.PSM","PI.PSMV")) 'red' else 'black', alpha=.05), border=NA)
}
legend('top', legend=c("Violate Perceptual Separability", "Other"), fill=c('red','black'), horiz=T, bty='n', cex=.6)
# y loc
plot(0,0,type='n',xlim=c(0,300),ylim=c(0,8),xlab="Trial",ylab="Parameter Estimate",main="Y Threshold")
abline(h=4,lty=1,lwd=2,xpd=F)
for (sm in unique(mydata$Version)) {
  conv.data <- subset(mydata, Block=="adapt" & Version==sm)
  lines(sort(unique(conv.data$Trial)), sapply(sort(unique(conv.data$Trial)), function (x) mean(conv.data$Alpha2[conv.data$Trial==x])), 
        col=if (sm %in% c()) 'red' else 'black')
  polygon(c(sort(unique(conv.data$Trial)),rev(sort(unique(conv.data$Trial)))), 
          c(sapply(sort(unique(conv.data$Trial)), 
                   function (x) mean(conv.data$Alpha2[conv.data$Trial==x]) + sd(conv.data$Alpha2[conv.data$Trial==x])),
            sapply(rev(sort(unique(conv.data$Trial))),
                   function (x) mean(conv.data$Alpha2[conv.data$Trial==x]) - sd(conv.data$Alpha2[conv.data$Trial==x]))),
          col=alpha(if (sm %in% c()) 'red' else 'black', alpha=.05), border=NA)
}
# y slope
plot(0,0,type='n',xlim=c(0,300),ylim=c(0,2),xlab="Trial",ylab="Parameter Estimate",main="Y SD")
abline(h=.8,lty=1,lwd=2, xpd=F)
for (sm in unique(mydata$Version)) {
  conv.data <- subset(mydata, Block=="adapt" & Version==sm)
  lines(sort(unique(conv.data$Trial)), sapply(sort(unique(conv.data$Trial)), function (x) mean(conv.data$Beta2[conv.data$Trial==x])), 
        col=if (sm %in% c()) 'red' else 'black')
  polygon(c(sort(unique(conv.data$Trial)),rev(sort(unique(conv.data$Trial)))), 
          c(sapply(sort(unique(conv.data$Trial)), 
                   function (x) mean(conv.data$Beta2[conv.data$Trial==x]) + sd(conv.data$Beta2[conv.data$Trial==x])),
            sapply(rev(sort(unique(conv.data$Trial))),
                   function (x) mean(conv.data$Beta2[conv.data$Trial==x]) - sd(conv.data$Beta2[conv.data$Trial==x]))),
          col=alpha(if (sm %in% c()) 'red' else 'black', alpha=.05), border=NA)
}
# acc
# PI MODELS IN RED
plot(0,0,type='n',xlim=c(0,300),ylim=c(.5,1),xlab="Trial",ylab="Accuracy Estimate")#,main="Predicted Accuracy")
abline(h=.75,lty=1,lwd=2, xpd=F)
for (sm in unique(mydata$Version)) {
  conv.data <- subset(mydata, Block=="adapt" & Version==sm)
  lines(sort(unique(conv.data$Trial)), sapply(sort(unique(conv.data$Trial)), function (x) mean(conv.data$Pred.Acc[conv.data$Trial==x])), 
        col=if (sm %in% c("PI","PI.PSM","PI.PSV","PI.PSMV")) 'red' else 'black')
  polygon(c(sort(unique(conv.data$Trial)),rev(sort(unique(conv.data$Trial)))), 
          c(sapply(sort(unique(conv.data$Trial)), 
                   function (x) mean(conv.data$Pred.Acc[conv.data$Trial==x]) + sd(conv.data$Pred.Acc[conv.data$Trial==x])),
            sapply(rev(sort(unique(conv.data$Trial))),
                   function (x) mean(conv.data$Pred.Acc[conv.data$Trial==x]) - sd(conv.data$Pred.Acc[conv.data$Trial==x]))),
          col=alpha(if (sm %in% c("PI","PI.PSM","PI.PSV","PI.PSMV")) 'red' else 'black', alpha=.05), border=NA)
}
legend('top', legend=c("Violate Perceptual Independence", "Other"), fill=c('red','black'), horiz=T, bty='n', cex=.6)

# Calculate assumption violation rates
grtResults <- matrix(nrow=2, ncol=32); m <- 0
for (sm in c("Recovery","PSM","PSV","PSMV","PI","PI.PSM","PI.PSV","PI.PSMV")) {
  m <- m + 1
  this.data <- subset(mydata, Version == sm)
  violations <- lapply(sort(unique(this.data$Subject)), function (s) c(1 * (riTest(table(this.data$StimInd[this.data$Block=="main" & this.data$Subject==s],this.data$RespInd[this.data$Block=="main" & this.data$Subject==s], exclude=NA))$p.value < .05),
                                                                    1 * (mriTest(table(this.data$StimInd[this.data$Block=="main" & this.data$Subject==s],this.data$RespInd[this.data$Block=="main" & this.data$Subject==s], exclude=NA))$p.value < .05)))
  violations <- sapply(1:8, function (x) sum(sapply(1:length(unique(this.data$Subject)), function (y) violations[[y]][x]), na.rm=T)/sum(sapply(1:length(unique(this.data$Subject)), function (y) !is.na(violations[[y]][x]))))
  grtResults[1, 1:4 + 4*(m-1)] <- violations[1:4]
  grtResults[2, 1:4 + 4*(m-1)] <- violations[5:8]
}
barplot(grtResults, beside=T, main="", xlab="Stimulus/Response", ylab=paste("Null Hypothesis Rejection Rate (",expression(alpha)," = .05)",sep=''),
        legend.text = c("Report Independence", "Marginal Response Invariance"), 
        args.legend = list(title="Test", bty='n', horiz=T, xpd=NA, x=48,y=-.2, xjust=.5,cex=.8, col='gray'))
abline(v=seq(12.5,by=12,length.out = length(unique(mydata$Version))-1), lty=2)





# Plot stimulus intensity estimates

       
       
       
