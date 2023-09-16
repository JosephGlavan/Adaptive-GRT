library(stringr)

# Import Data
mydata <- read.table("/Users/jjg367/Documents/Outside Research/Adaptive GRT/Psi Source/Test_Simulation/data/2023-09-15_21h12.50.080.txt", sep="\t", header=T)

# Add Missing Columns
mydata$Dim1 <- sapply(mydata$Stimulus, function (x) as.numeric(unlist(strsplit(str_sub(x,2,-2), ","))[1]))
mydata$Dim2 <- sapply(mydata$Stimulus, function (x) as.numeric(unlist(strsplit(str_sub(x,2,-2), ","))[2]))
mydata$StimInd <- NA
for (s in unique(mydata$Subject)) {
  # 00 01 10 11
  d1 <- range(mydata$Dim1[mydata$Block=="main" & mydata$Subject==s])
  d2 <- range(mydata$Dim2[mydata$Block=="main" & mydata$Subject==s])
  mydata$StimInd[mydata$Block=="main" & mydata$Subject==s] <- sapply(mydata$Stimulus[mydata$Block=="main" & mydata$Subject==s], function (x) {
    mystim <- eval(parse(text=paste0("c",x)))
    if (mystim[1] == d1[1] & mystim[2] == d2[1]) 0 else
      if (mystim[1] == d1[1] & mystim[2] == d2[2]) 1 else
        if (mystim[1] == d1[2] & mystim[2] == d2[1]) 2 else
          if (mystim[1] == d1[2] & mystim[2] == d2[2]) 3
  })
} 
mydata$RespInd <- NA; mydata$RespInd[mydata$Response=="(0, 0)"] <- 0; mydata$RespInd[mydata$Response=="(0, 1)"] <- 1; mydata$RespInd[mydata$Response=="(1, 0)"] <- 2; mydata$RespInd[mydata$Response=="(1, 1)"] <- 3
mydata$Correct <- mydata$StimInd == mydata$RespInd
mydata$Alpha1 <- sapply(mydata$Lambda, function (x) as.numeric(unlist(strsplit(gsub("[()]","",x), ",")))[1])
mydata$Beta1 <- sapply(mydata$Lambda, function (x) as.numeric(unlist(strsplit(gsub("[()]","",x), ",")))[2])
mydata$Alpha2 <- sapply(mydata$Lambda, function (x) as.numeric(unlist(strsplit(gsub("[()]","",x), ",")))[3])
mydata$Beta2 <- sapply(mydata$Lambda, function (x) as.numeric(unlist(strsplit(gsub("[()]","",x), ",")))[4])

# For only one subject:
# par(mfrow=c(2,2))
# plot(mydata$Trial[mydata$Block=="adapt"], mydata$Alpha1[mydata$Block=="adapt"], type='l', main="Dimension 1", xlab="Trial", ylab="Alpha (Decision Threshold)"); abline(h=600,lty=2)
# plot(mydata$Trial[mydata$Block=="adapt"], mydata$Beta1[mydata$Block=="adapt"], type='l', main="Dimension 1", xlab="Trial", ylab="Beta (Stimulus Separation)"); abline(h=100,lty=2)
# plot(mydata$Trial[mydata$Block=="adapt"], mydata$Alpha2[mydata$Block=="adapt"], type='l', main="Dimension 2", xlab="Trial", ylab="Alpha (Decision Threshold)"); abline(h=4,lty=2)
# plot(mydata$Trial[mydata$Block=="adapt"], mydata$Beta2[mydata$Block=="adapt"], type='l', main="Dimension 2", xlab="Trial", ylab="Beta (Stimulus Separation)"); abline(h=.8,lty=2)
# par(mfrow=c(1,1))



