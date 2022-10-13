
pdf("BurnoutPlotInclR0Ver4.pdf", height = 4, width = 6)

Burnout<-function(R0,S0,I0,i){ #This function is a straight line

	return (1 - i - exp(-(R0)*(S0*i + I0)) )
  
}



R0 = 4.0;

par(mfrow=c(1,1));
par(mai=c(1.2,0.55,0.5,0.5));
LineWidths = c(3,2,1)

#Using the burnout approximation to calculate the fraction infected at the end of the epidemic

rootStor = numeric(); #Declare a vector to store the roots in
Low = 1e-10;
S0Range = seq(from=Low,to=1.0,by=0.01); #the range of values over which we will calculate the root

I0Range = c(Low,1e-3,1e-2);
#I0Range = c(1e-10,1e-2,1e-1);


k = 1;
for(I0 in I0Range){ 
	j = 1; #j keeps track of how many roots we have calculated
	for(S0 in S0Range){ #loop over the values of S0
		#z = uniroot(f=Burnout,interval=c(0.0,1.0),S0=S0,I0=I0);
		z = uniroot(f=Burnout,interval=c(Low,(1.0)),R0=R0,S0=S0,I0=I0);
		
		rootStor[j] = z$root; #store the root
		j = j + 1; #update j
	}

	if(k==1){
		plot(S0Range,rootStor,type="l",xlab="",ylab="",axes=FALSE,xlim=c(0,max(S0Range)),ylim=c(0,1),lwd=LineWidths[k]); #plot the roots versus the values of m
	}else{
		lines(S0Range,rootStor,lwd=LineWidths[k]);
	}
	k = k + 1;
}


#par(mai=c(0.45,0.45,0.0,0.0));

lgnd = c();
lgnd[3] = expression(I[r](0)~paste("=")~10^-10);
lgnd[2] = expression(I[r](0)~paste("=")~10^-3);
lgnd[1] = expression(I[r](0)~paste("=")~10^-2);
LineWidths = c(1,2,3)

#text(x=0.35,y=0.05,labels=lgnd[1],bty="n",cex=0.7)
#text(x=0.1,y=0.17,labels=lgnd[2],bty="n",cex=0.7)
#text(x=0.1,y=0.3,labels=lgnd[3],bty="n",cex=0.7)

text(x=0.725,y=0.6,labels="Initial Fraction Infected",bty="n",cex=0.8)
lbl = expression(Initial~Infected~I[r](0))
legend(x=0.5,y=0.6,legend=lgnd,bty="n",lwd = LineWidths,cex=0.8);
axis(1,pretty(range(S0Range),5),pos=0.0,padj=-0.5);

axis(2,pretty(range(rootStor),5),pos=0.0);
yaxislbl = expression(Cumulative~Fraction~Infected~Z[r](infinity))
mtext(side=2,text=yaxislbl,padj=-1.65);

xaxislbl = expression(Initial~Fraction~Susceptible~S[0]);
mtext(side=1,text=xaxislbl,padj=1.75);
xaxislbl = expression(Reproductive~Number~R[0][r])
mtext(side=1,text=xaxislbl,padj=5.75);

axis(side=1,at=c(0,0.25,0.5,0.75,1.0),labels=c(0,1,2,3,4),line=3,pos=-0.285,padj=-0.5);

dev.off()

