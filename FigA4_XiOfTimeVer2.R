pdf("XiOverTimeVer2.pdf", height = 4, width = 4)

require(deSolve);


Burnout<-function(I0,R0,z){
  
  output = 1 - z - exp(-R0 * ((1-I0)*z + I0));
  #cat("i:",i,"Output:",output,"\n")
  return(output);
  #return(0);
  
}


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
    dS.dt = - R0r*S*Ir -u*R0i*S*Ii - u*R0i*S*Ii;
    dIr.dt = R0r*S*Ir - Ir;
    dIi.dt = u*R0i*S*Ii - u*Ii;
    dZr.dt = Ir;
    dZi.dt = u*Ii;
    return(list(c(dS.dt,dIr.dt,dIi.dt,dZr.dt,dZi.dt)));
  })
  
}

R0r = 2; u = 0.7;
mui = 0.1; mur = 0.2; betar = 2.4;
Seqbeta = seq(from=0.101,to=0.15,by=0.001);
betai = 0.15; R0i = betai/mui;
betai2 = 0.105; R0i2 = betai2/mui
u = mui/mur;
gamma1 = u*(R0i-1)/(R0r-1);
gamma2 = u*(R0i2-1)/(R0r-1);

IrStart = 1e-5; IiStart = 1e-5;
SStart = 1.0-IrStart-IiStart;

y0 = c(SStart,IrStart,IiStart,0.0,0.0);
Steps  = 1e-3; 
EpidTime = 25;
t = seq(from=0,to=EpidTime,by=Steps);
p = list(R0r=R0r,R0i=R0i,u=u);
out = ode(y=y0,times=t,func=TwoPathEpidVer2,parms=p);

XiZero = 1;
y0 = c(XiZero);
StepSize = 1e-5; #5e-6;
Low = 1e-5;
p = list(R0r=R0r,R0i = R0i,u=u);
FractInf = uniroot(f=Burnout,interval=c(0.0,1.0),I0=Low,R0=R0r);
t = seq(from=Low,to=FractInf$root,by=StepSize);
out2 = ode(y=y0,times=t,func=XiODE,parms=p); 

#plot(out[,1],out[,2],type="l",ylim=c(0,1))
#lines(out[,1],out[,3])


#plot(1-out[,2],type="l")
#lines(out2[,1])

#plot(out2[,1],type="l")
#lines(1-out[,2])

TimeStor = c(); XiStor = c(); xStor = c();
oldTimeIndex = 0;
j = 1;
EndSimuln = nrow(out);
End = nrow(out2); 
i = 1; Stop = 0;
j = 1;
while((i<=End)&&(Stop!=1)){
#while(j<8316){

  
  #Value = 1 - out[i,2]; 
  #Index = which(out2[,1]>Value)
  Value = out2[i,1]; 
  Index = which((1-out[,2])>Value)
  IndexLength = length(Index)
  if(IndexLength<1) Stop = 1;
  if(Stop<1){
    TimeIndex = min(Index);
    if((TimeIndex>oldTimeIndex)&&(TimeIndex<EndSimuln)){
      TimeStor[j] = out[TimeIndex,1]; XiStor[j] = out2[i,2]; 
      xStor[j] = out2[i,1];
      j = j + 1;
    } else{
      if(TimeIndex>EndSimuln) Stop = 1;
    }
    oldTimeIndex = TimeIndex;
  } 
  i = i + 1;
}


par(mai=c(0.45,0.5,0.25,0.1))
arg = IrStart*R0r/(R0r-1);
yArg = (IiStart*((xStor/arg)^gamma1)*XiStor);
End = length(yArg);
TimeStorB = TimeStor[2:End]
yArgB = yArg[2:End]
plot(TimeStorB,yArgB,type="l",axes=FALSE,xlab="",ylab="",xlim=c(0,EpidTime),ylim=c(0,max(yArg)),col="RED",lwd=1);
#plot(TimeStor,yArg,type="l",axes=FALSE,xlab="",ylab="",xlim=c(0,25));

axis(1,pretty(range(TimeStor),5),pos=0.0);
mtext(side=1,text="Time",padj=1.75);
axis(2,pretty(range(yArg),5),pos=0);
yLabel = expression(Invader~Fraction~Infected~I[i](t));

mtext(side=2,text=yLabel,padj=-1.25);
legend(x=11,y=1e-4,legend=c("Full model","Approximate model"),bty="n",cex=0.8,col=c("BLACK","RED"),lwd=c(1,1));

par(new="TRUE")
#yaxisArg = (out[,3]^gamma)*out[,4];
yaxisArg = (out[,4]);
#plot(out[,1],out[,4],type="l",col="RED",xlim=c(0,25))
plot(out[,1],yaxisArg,col="BLACK",type="l",axes=FALSE,xlab="",ylab="",xlim=c(0,EpidTime),ylim=c(0,max(yArg)),lwd=1)

#plot(out[,1],yaxisArg,col="RED",type="l",axes=FALSE,xlab="",ylab="",xlim=c(0,25),ylim=c(-6,max(yArg)))

#axis(4,pretty(range(out[,4]),5),pos=EpidTime);
#mtext(side=4,text="Full model",padj=2.2,col="RED");

dev.off()


              

              