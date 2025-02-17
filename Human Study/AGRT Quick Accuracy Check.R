#### FILE PATH ####
myfile <- "/Users/joeglavan/Documents/Research Projects/Adaptive GRT/Human Study/Experiment Code/data/TEST_2021_Nov_12_1623.txt"
###################

mydata <- read.table(myfile, header=T, sep="\t", stringsAsFactors = F)
mydata <- subset(mydata, Block=="MAIN")
print("Accuracy by Stimulus",quote=F)
for (s in 0:3) {
  print(paste("Stimulus ",s," (",switch(s+1, "LL","LH","HL","HH"),"): ", round(mean(mydata$Correct[mydata$Stimulus==s]), 4), sep=''), quote=F)
}
print("Overall Accuracy (targeting 75%)", quote=F)
print(round(mean(sapply(0:3, function (x) mean(mydata$Correct[mydata$Stimulus==x]))), 4), quote=F)

if (chisq.test(table(mydata$Stimulus, mydata$Response))$p.value < .05) print("Responses are statistically non-random. ACCEPT SUBJECT.", quote=F) else print("Responses are statistically random. REJECT SUBJECT.", quote=F)
