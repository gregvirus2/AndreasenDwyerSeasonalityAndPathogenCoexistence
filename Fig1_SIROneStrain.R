require(deSolve);

pdf("SIROneStrain.pdf", height = 4, width = 4)


SIR<-function(t,y,p){
	S = y[1];
	I = y[2];
	with(as.list(p),{
		dS.dt = -R0*S*I;
		dI.dt = R0*S*I - I;
		return(list(c(dS.dt,dI.dt)));
	})
}

par(mfrow=c(1,1));
#par=(mai=c(0.1,0.1,0.1,0.1))

par(mai=c(0.6,0.55,0.05,0.05))

R0 = 2.3;

p2 = list(R0=R0);

#t2 = c(1,5,10,20); #stupid
timeHi = 20;
t2 = seq(from=0,to=timeHi,by=0.1);

I0 = 0.001;
S0 = 1.0 - I0; 
N0 = c(S0,I0);
out = ode(y=N0,times=t2,func=SIR,parms=p2);

plot(out[,1],out[,2],type="l",xlab="Time",ylab="S and I",ylim=c(0,1),lwd=2,col=1,xlim=c(0,timeHi),axes=FALSE);
lines(out[,1],out[,3],col=2);
lgndLbl = c();
lgndLbl[1] = expression(Susceptible~S[r]);
lgndLbl[2] = expression(Infected~I[r]);
legend(x=8,y=0.6,legend=c(lgndLbl[1],lgndLbl[2]),bty="n",lwd = c(2,1),col=c("BLACK","RED"));

axis(1,pretty(range(out[,1]),5),pos=0.0);
mtext(side=1,text="Time",padj=1.75);
axis(2,pretty(range(out[,2]),5),pos=0);
y <- expression(Susceptible~Hosts~S[r]~and~Infected~Hosts~I[r])
mtext(side=2,text=y,padj=-2.0);

dev.off()






		