pdf("CoexistEpiBurnSurvRatesPlotB.pdf", height = 4, width = 4)
par(mai=c(0.65,0.65,0.25,1.25))


if(1){
par(mfrow=c(2,1))
Data = read.csv("CoexistEpiBurnoutSurvivalRatesA.csv")
plot(log10(Data$Wr),log10(Data$Wi.Hi.25),type="l",col="RED",xlab="",ylab="");
lines(log10(Data$Wr),log10(Data$Wi.Hi.25),col="BLACK");
mtext(side=1,text="log10 W_r",padj=3);
mtext(side=2,text="log10 W_i",padj=-3.5);

Data = read.csv("CoexistEpiBurnoutSurvivalRatesC.csv")
plot(log10(Data$Wr),log10(Data$Wi.Hi.25),type="l",col="BLACK",xlab="",ylab="");
#lines(log10(Data$Wr),log10(Data$Wi.Hi.250),col="BLACK");
mtext(side=1,text="log10 W_r",padj=3);
mtext(side=2,text="log10 W_i",padj=-3.5);
}


if(0){
  
  par(mfrow=c(2,1))
  Data = read.csv("CoexistEpiBurnoutSurvivalRatesA.csv")

  plot((Data$Wr),(Data$Wi.Hi.25),type="l",axes="FALSE",xlab="",ylab="");

  axis(1,pretty(range((Data$Wr)),5),pos=0.0);
  mtext(side=1,text="log10 Resident Survival W_r",padj=2.5);

  axis(2,pretty(range((Data$Wi.Hi.25)),5),pos=0)
  mtext(side=2,text="log10 Invader Survival W_i",padj=-3);

Data = read.csv("CoexistEpiBurnoutSurvivalRatesC.csv")

plot(Data$Wr,Data$Wi.Hi.25,type="l",axes="FALSE",xlim=c(0,0.01),xlab="",ylab="");
axis(1,pretty(range(Data$Wr),5),pos=0.0);
mtext(side=1,text="log10 Resident Survival W_r",padj=2.5);

#mtext(side=1,text="Time",padj=1.85);
axis(2,pretty(range(Data$Wi.Hi.25),10),pos=0)
mtext(side=2,text="log10 Invader Survival W_i",padj=-2.5);
#abline(v=0.007)
}

dev.off()

