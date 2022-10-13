#pdf("TwoPathEpiVsXiVer3B.pdf", height = 4, width = 6)

require(deSolve);

########### Function definitions ###############

####### Burnout approximation
Burnout<-function(I0,R0,z){
  
  output = 1 - z - exp(-R0 * ((1-I0)*z + I0));
  #cat("i:",i,"Output:",output,"\n")
  return(output);
  #return(0);
  
}

######## ODE function for Xi, the approximation

XiODE<-function(t,y,p){
  Xi = y[1];
  
  with(as.list(p), {
    
    phiOft = (t + log(1-t))/(t^2);
    gamma = u*(R0i - 1)/(R0r - 1)
    numer = (u/(t-1)) - gamma*phiOft
    denom = R0r - 1 + t*phiOft;
    
    #dXh.dt = nuH - r*Thm*Ym*Xh - muH*Xh
    dXi.dt = (numer/denom)*Xi
    return(list(c(dXi.dt)));
  })
}

############## ODE function for the SIR model with 2 strains
TwoPathEpidVer2<-function(t,y,p){
  S = y[1];
  Ir = y[2];
  Ii = y[3];
  Zr = y[4];
  Zi = y[5];
  
  with(as.list(p), {
    dS.dt = - R0r*S*Ir - u*R0i*S*Ii;
    dIr.dt = R0r*S*Ir - Ir;
    dIi.dt = u*R0i*S*Ii - u*Ii;
    dZr.dt = Ir;
    dZi.dt = u*Ii;
    return(list(c(dS.dt,dIr.dt,dIi.dt,dZr.dt,dZi.dt)));
  })
  
}


R0r = 3; R0i = 1.05; u = 0.5;
Steps  = 1e-3; #It's convenient to define the number of steps, for use in spectrum
t = seq(from=0,to=50,by=Steps);
IrStart = 1e-3; IiStart = 1e-10;
SStart = 1.0-IrStart-IiStart;
eps = IrStart;

y0 = c(SStart,IrStart,IiStart,0.0,0.0);

p = list(R0r=R0r,R0i=R0i,u=u);
out = ode(y=y0,times=t,func=TwoPathEpidVer2,parms=p);


if(0){
  plot(out[,1],out[,2],type="l",ylim=c(0,1.0))
  lines(out[,1],out[,3],col="RED")
  lines(out[,1],out[,4],col="BLUE")
}
if(0){
par(mfrow=c(1,2))
plot(out[,1],out[,3],col="RED",type="l")
lines(out[,1],out[,4],col="BLUE")

plot(out[,1],out[,5],col="RED",type="l")
lines(out[,1],out[,6],col="BLUE")
}

end = length(out[,1])
FI = 1 - out[end,2]/out[1,2]


require(deSolve);


XiZero = 1;
y0 = c(XiZero);
StepSize = 1e-4;
p = list(R0r=R0r,R0i = R0i,u=u)

Low = 1e-10;
FractInf = uniroot(f=Burnout,interval=c(0.0,1.0),I0=Low,R0=R0r);

par(mfrow=c(1,1));
par(mai=c(0.75,0.75,0.25,0.25));
t = seq(from=Low,to=FractInf$root,by=StepSize);
out2 = ode(y=y0,times=t,func=XiODE,parms=p); 
plot(out2[,1],out2[,2],type="l",col="RED",lwd=3,axes="FALSE",xlab="",ylab="",xlim=c(0,1.0))
axis(1,pretty(range(out2[,1]),5),pos=0.0);
#axis(1,c(0,0.5,1.0),pos=0.0);

mtext(side=1,text="Cumulative Fraction Infected By Either Strain",padj=2.75);
axis(2,pretty(range(out2[,2]),5),pos=0);
yAxisLbl = expression(xi~or~Scaled~I[i]);

mtext(side=2,text=yAxisLbl,padj=-1.75);

lgnd0 = expression(I[r](0)~paste("=")~0.001);
lgnd1 = expression(xi);
lgnd2 = expression(Scaled~I[i])
lgnd=c(lgnd1,lgnd2);
text(x=0.35,y=0.55,labels=lgnd0,bty="n")
legend(x=0.25,y=0.5,legend=lgnd,bty="n",lwd = c(1,2),cex=1.0,col=c("BLACK","RED"));


gam = u*(R0i-1)/(R0r-1);
arg1 = R0r^gam;
omega = u*(R0i-1);
eps = IrStart/omega;
#arg2 = (1 +  (1-out[,2])/(IiStart*gam))^(-gam)
arg2 = (1 +  (1-out[,2])/(eps*gam))^(-gam)

#plot(1-out[,2],(1/IiStart)*arg1*arg2*out[,4],type="l")
lines(1-out[,2],(1/IiStart)*arg1*arg2*out[,4])

#dev.off();


