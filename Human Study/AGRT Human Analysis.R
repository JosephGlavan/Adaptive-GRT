require(mdsdt)

ma <- function(x, n = 5) as.vector(filter(x, rep(1 / n, n), sides = 2))

##### PILOT-CONDITION #####

# "good" Subject IDs
good.subjects <- c(1, 5, 28, 29, 31, 32, 33, 34, 35, 37, # separable condition
                   19, 20, 21, 22, 23, 24, 25, 26, 30, 36) # integral condition

## Separable Subject 2 needs to reverse brightness and redo results


mydata <- do.call(rbind.data.frame, lapply(list.files(path="/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/agrt-2023-04-03", full.names=T, recursive=T), function (x) {
  if (as.integer(unlist(strsplit(basename(x), "_"))[1]) %in% good.subjects) read.table(x, header=T, stringsAsFactors = F)
}))

pract.acc.p <- list()
main.acc.p <- list()
mri.p <- list()
ri.p <- list()
for (exp in 1:length(unique(mydata$Experiment))) {
  for (subj in 1:length(unique(mydata$Subject[mydata$Experiment==unique(mydata$Experiment)[exp]]))) {
    this.data <- subset(mydata, Experiment==sort(unique(mydata$Experiment))[exp] & Subject==sort(unique(mydata$Subject[mydata$Experiment==sort(unique(mydata$Experiment))[exp]]))[subj])
    
    # Fix Separable Subject 2 who seems to have reversed the orientation instructions
    if (exp==2 & subj==2) {
      # Swap 0 and 1, 2 and 3
      this.data$Response[this.data$Response==0] <- 4
      this.data$Response[this.data$Response==1] <- 0
      this.data$Response[this.data$Response==4] <- 1
      this.data$Response[this.data$Response==2] <- 4
      this.data$Response[this.data$Response==3] <- 2
      this.data$Response[this.data$Response==4] <- 3
      this.data$Correct <- 1*(this.data$Response==this.data$Stimulus)
    }
    
    ### Accuracy by Experiment/Subject/Stimulus/Overall
    p.means <- sapply(0:3, function (x) mean(this.data$Correct[this.data$Block=="PRACT" & this.data$Stimulus==x]))
    p.sds <- sapply(0:3, function (x) sd(this.data$Correct[this.data$Block=="PRACT" & this.data$Stimulus==x]))
    pract.acc.p[[length(pract.acc.p)+1]] <- c(c("I","S")[exp], subj, rbind(p.means,p.sds))
    m.means <- sapply(0:3, function (x) mean(this.data$Correct[this.data$Block=="MAIN" & this.data$Stimulus==x]))
    m.sds <- sapply(0:3, function (x) sd(this.data$Correct[this.data$Block=="MAIN" & this.data$Stimulus==x]))
    main.acc.p[[length(main.acc.p)+1]] <- c(c("I","S")[exp], subj, rbind(m.means,m.sds))
    png(paste0("/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/Figures/PC Practice Accuracy MA ", c("Integral","Separable")[exp]," ", subj, ".png"), width=600, height=600)
    plot(1:160,ma(this.data$Correct[this.data$Block=="PRACT"], 16), type='l', ylim=c(0,1), main=paste0("Practice Block\n", c("Integral Subject ", "Seperable Subject ")[exp], subj), xlab="Trial", ylab="P(correct)")
    dev.off()
    png(paste0("/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/Figures/PC Main Accuracy MA ", c("Integral","Separable")[exp]," ", subj, ".png"), width=600, height=600)
    plot(1:1000,ma(this.data$Correct[this.data$Block=="MAIN"], 16), type='l', ylim=c(0,1), main=paste0("Main Block\n", c("Integral Subject ", "Seperable Subject ")[exp], subj), xlab="Trial", ylab="P(correct)")
    dev.off()
    
    ### GRT Results
    # Build mri / ri test results
    # Plots (equal-likelihood countours / plot.fit.grt)
    srm <- table(this.data$Stimulus[this.data$Block=="MAIN"],this.data$Response[this.data$Block=="MAIN"], exclude=-1)
    png(paste0("/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/Figures/PC GRT Model ", c("Integral","Separable")[exp]," ", subj, ".png"), width=600, height=600)
    plot(fit.grt(srm), marginals=T, xlab=c("Size", "Saturation")[exp], ylab=c("Orientation", "Brightness")[exp], main=paste0(c("Integral Subject ", "Seperable Subject ")[exp], subj))
    dev.off()
    mri <- mriTest(srm); ri <- riTest(srm)
    mri.p[[length(mri.p)+1]] <- c(c("I","S")[exp], subj, rbind(mri$z,mri$p.value))
    ri.p[[length(ri.p)+1]] <- c(c("I","S")[exp], subj, rbind(ri$chi.2,ri$p.value))
  }
}
pract.acc.p <- data.frame(do.call(rbind, pract.acc.p)); names(pract.acc.p) <- c("Experiment", "Subject", "A1B1.M", "A1B1.SD", "A1B2.M", "A1B2.SD", "A2B1.M", "A2B1.SD", "A2B2.M", "A2B2.SD"); for (i in 2:ncol(pract.acc.p)) pract.acc.p[,i] <- as.numeric(pract.acc.p[,i])
main.acc.p <- data.frame(do.call(rbind, main.acc.p)); names(main.acc.p) <- c("Experiment", "Subject", "A1B1.M", "A1B1.SD", "A1B2.M", "A1B2.SD", "A2B1.M", "A2B1.SD", "A2B2.M", "A2B2.SD"); for (i in 2:ncol(main.acc.p)) main.acc.p[,i] <- as.numeric(main.acc.p[,i])
mri.p <- data.frame(do.call(rbind, mri.p)); names(mri.p) <- c("Experiment", "Subject", "(A1,-).Z", "(A1,-).P", "(A2,-).Z", "(A2,-).P", "(-,B1).Z", "(-,B1).P", "(-,B2).Z", "(-,B2).P"); for (i in 2:ncol(mri.p)) mri.p[,i] <- as.numeric(mri.p[,i])
ri.p <- data.frame(do.call(rbind, ri.p)); names(ri.p) <- c("Experiment", "Subject", "A1B1.X2", "A1B1.P", "A1B2.X2", "A1B2.P", "A2B1.X2", "A2B1.P", "A2B2.X2", "A2B2.P"); for (i in 2:ncol(ri.p)) ri.p[,i] <- as.numeric(ri.p[,i])
# Print Tables and Group Figures
write.csv(pract.acc.p, file="/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/Practice Accuracy Pilot Condition.csv", row.names = F)
write.csv(main.acc.p, file="/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/Main Accuracy Pilot Condition.csv", row.names = F)
write.csv(mri.p, file="/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/MRI Test Results Pilot Condition.csv", row.names = F)
write.csv(ri.p, file="/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/RI Test Results Pilot Condition.csv", row.names = F)


pi.acc.m <- t(as.matrix(main.acc.p[main.acc.p$Experiment=="I", seq(3,by=2,length.out=4)]))
pi.acc.group.m <- apply(pi.acc.m, 2, mean); pi.acc.group.sd <- apply(pi.acc.m, 2, sd)
png(paste0("/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/Figures/PC Main Accuracy Stimulus Integral.png"), width=600, height=600)
barplot(pi.acc.m, beside=T, ylim=c(0,1), main="Main Block Accuracy", ylab="P(correct)", xlab="Subject", names.arg=1:10)#, legend.text=c("A1B1","A1B2","A2B1","A2B2"), args.legend=list(title="Stimulus") )
dev.off()
png(paste0("/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/Figures/PC Main Accuracy Total Integral.png"), width=600, height=600)
barplot(pi.acc.group.m, beside=F, ylim=c(0,1), main="Main Block Accuracy", ylab="P(correct)", xlab="Subject", names.arg=1:10)#, legend.text=c("A1B1","A1B2","A2B1","A2B2"), args.legend=list(title="Stimulus") )
abline(h=.75, lty=2)
dev.off()
ps.acc.m <- t(as.matrix(main.acc.p[main.acc.p$Experiment=="S", seq(3,by=2,length.out=4)]))
ps.acc.group.m <- apply(ps.acc.m, 2, mean); ps.acc.group.sd <- apply(ps.acc.m, 2, sd)
png(paste0("/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/Figures/PC Main Accuracy Stimulus Separable.png"), width=600, height=600)
barplot(ps.acc.m, beside=T, ylim=c(0,1), main="Main Block Accuracy", ylab="P(correct)", xlab="Subject", names.arg=1:10)#, legend.text=c("A1B1","A1B2","A2B1","A2B2"), args.legend=list(title="Stimulus") )
dev.off()
png(paste0("/Users/jjg367/Documents/Adaptive GRT/Human Study/Analysis/Figures/PC Main Accuracy Total Separable.png"), width=600, height=600)
barplot(ps.acc.group.m, beside=F, ylim=c(0,1), main="Main Block Accuracy", ylab="P(correct)", xlab="Subject", names.arg=1:10)#, legend.text=c("A1B1","A1B2","A2B1","A2B2"), args.legend=list(title="Stimulus") )
abline(h=.75, lty=2)
dev.off()



##### ADAPTIVE-CONDITION #####
