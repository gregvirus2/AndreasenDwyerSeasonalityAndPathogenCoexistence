pdf("DiVersusBetaOrMuVer8.pdf", height = 4, width = 6)

require(deSolve);

Burnout<-function(I0,R0,z){
  
  output = 1 - z - exp(-R0 * ((1-I0)*z + I0));
  #cat("i:",i,"Output:",output,"\n")
  return(output);
  #return(0);
  
}

IeqmCalc<-function(R0,W,I0){
  
  #arg = ((w+1)/w)*I0 - (I0^2)/w;
  arg = (1/W)*I0*(1-I0) + I0;
  output = 1 - (I0/W) - exp(-R0*arg);
  #cat("I:",I0,"Output:",output,"\n")
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


R0iStor = c();
DiStor = c();
gamStor = c();
betaStor = c();
muStor = c();
XiStor = c();

#par(mfrow=c(3,4))
par(mfrow=c(1,1))
#par(ask="TRUE")

par(mai=c(0.45,0.45,0.05,0.1))

R0r = 2.5;
u = 1;

LineWidths = c(3,2,1)


plotNum = 1;


XiZero = 1;
y0 = c(XiZero);
StepSize = 1e-4;
#R0r = 1.25; R0i = 2.4; u = 0.7;
#R0r = 4.0; R0i = 3; u = 0.5;
#R0r = 1.5;
#R0i = 3; 
#u = 0.5;

mui = 0.5; 

IrStart = INought = 1e-3;
Low = 1e-2;
FractInf = uniroot(f=Burnout,interval=c(0.0,1.0),I0=Low,R0=R0r);
t = seq(from=Low,to=FractInf$root,by=StepSize);

maxR0i = 1/(1-FractInf$root);
UpperThold = u*(maxR0i-1)/(R0r-1);
if(UpperThold>3) UpperThold = 3;
betaCount = 25;
Step = (UpperThold)/betaCount;
Start = Step; 
#Last = min((End-(5*Step)),3);
cat("u",u,"R0r:",R0r,"UpperThold:",UpperThold,"\n")
Seqgam = seq(from=Start,to=UpperThold-5*Step,by=Step)
cat("Seqgam:",Seqgam,"\n");


Seqmu = seq(from=0.01,to=0.09,by=0.001)

i = 1;

if(plotNum>1){
  
  rm(R0iStor);
  rm(DiStor);
  rm(gamStor);
  rm(betaStor);
  rm(muStor);
  rm(XiStor);
  
  R0iStor = c();
  DiStor = c();
  gamStor = c();
  betaStor = c();
  muStor = c();
  XiStor = c();
}
  


 
#for(betai in Seqbeta){

for(gamma in Seqgam){
  
  R0i = (gamma/u)*(R0r-1) + 1;
  
  Low = 1e-5;
  FractInf = uniroot(f=Burnout,interval=c(0.0,1.0),I0=Low,R0=R0r);
  t = seq(from=Low,to=(FractInf$root-Low),by=StepSize);
  p = list(R0r=R0r,R0i = R0i,u=u);
  out = ode(y=y0,times=t,func=XiODE,parms=p); 
  End = nrow(out);
  XiStor[i] = out[End,2];
  
  if(R0i*(1-FractInf$root)>1){
    cat("problem - R0i:",R0i,"1-FI:",1-FractInf$root,"Product:",R0i*(1-FractInf$root),"\n")
  }
  
  #gamma = u*(R0i-1)/(R0r-1);
  
  gamStor[i] = gamma;
  #betaStor[i] = betai;
  muStor[i] = mui;
  
  numer = (out[,1]^gamma)*out[,2];
  denom = (1-out[,1])*(R0r*out[,1] + log(1-out[,1]))
  arg = u*(((R0r-1)/R0r)^gamma);
  Di = StepSize*arg*sum(numer/denom)
  last = length(Di)
  
  numer1 = (out[1,1]^gamma)*out[1,2];
  denom1 = (1-out[1,1])*(R0r*out[1,1] + log(1-out[1,1]))
  FirstTerm = StepSize*arg*(numer1/denom1);
  numer1 = (out[last,1]^gamma)*out[last,2];
  denom1 = (1-out[last,1])*(R0r*out[last,1] + log(1-out[last,1]))
  LastTerm = StepSize*arg*(numer1/denom1);
  
  Di = Di - 0.5*(FirstTerm + LastTerm);
  
  R0iStor[i] = R0i;
  DiStor[i] = Di;
  #cat("gamma:",gamma,"R0i:",R0i,"Di:",Di,"\n");
  
  i = i + 1;
}


plot(gamStor,log10((INought^(-gamStor))*DiStor),type="l",xlab="",ylab="",axes="FALSE",lwd=LineWidths[1],ylim=c(-1,max(log10((INought^(-gamStor))*DiStor))));

axis(1,pretty(range(gamStor),5),pos=0);

#axis(2,pretty(range(log10((INought^(-gamStor))*DiStor)),5),pos=0.12);
#axis(2,pretty(range(log10(DiStor)),5),pos=0.12);
x = c(-0.5,0,1,2,3,4,5,6,7)
axis(2,pretty(range(x),10),pos=0.12);


#axis(2,pretty(-5:7,5),pos=0.12);

xaxisLbl = expression(Initial~relative~invader~fitness~gamma[0])

#yaxisLbl = expression(log[10]~Fitness~rise~Z[i](infinity)/I[i](0)~or~component~I[r](0)^-gamma[0]~or~D[i])
yaxisLbl = expression(log[10]~Fitness~and~its~components)

mtext(side=1,text=xaxisLbl,padj=-0.25,cex=1.25);
mtext(side=2,text=yaxisLbl,padj=-0.75,cex=1.25);


#par(new="TRUE")
#plot(gamStor,log10((INought^(-gamStor))),type="l",xlab="",lwd=LineWidths[2],col="BLUE",ylab="",axes="FALSE",ylim=c(0,max(log10((INought^(-gamStor))))));
lines(gamStor,log10((INought^(-gamStor))),type="l",lwd=LineWidths[2],col="BLUE");
lines(gamStor,log10(DiStor),type="l",col="RED",lwd=LineWidths[3])
#axis(4,pretty(range(log10((INought^(-gamStor)))),5),pos=max(gamStor));

lgnd0 = expression(log[10]~Z[i](infinity)/I[i](0));
lgnd1 = expression(log[10]~I[r](0)^-gamma[0]);
lgnd2 = expression(log[10]~D[i]);
lgnd=c(lgnd0,lgnd1,lgnd2);
legend(x=0.2,y=5.5,legend=lgnd,bty="n",lwd = c(3,2,1),cex=1.0,col=c("BLACK","BLUE","RED"));


dev.off();



