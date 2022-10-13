par(mfrow=c(1,5));


require(deSolve);

OneStrainEqm<-function(R0,W,i){
  
  output = 1 - i - exp(-R0 * ((1-W*i)*i + W*i));
  #cat("i:",i,"Output:",output,"\n")
  return(output);
  #return(0);
  
}

Burnout<-function(I0,R0,z){
  
  output = 1 - z - exp(-R0 * ((1-I0)*z + I0));
  #cat("i:",i,"Output:",output,"\n")
  return(output);
  #return(0);
  
}


############## ODE function for the SIR model with 2 strains
TwoPathEpidVer2<-function(t,y,p){
  S = y[1];
  Ir = y[2];
  Ii = y[3];
  Zr = y[4];
  Zi = y[5];
  
  with(as.list(p), {
    dS.dt = - betar*S*Ir -betai*S*Ii;
    dIr.dt = betar*S*Ir - mur*Ir;
    dIi.dt = betai*S*Ii - mui*Ii;
    dZr.dt = mur*Ir;
    dZi.dt = mui*Ii;
    return(list(c(dS.dt,dIr.dt,dIi.dt,dZr.dt,dZi.dt)));
  })
  
}


#par(mai=c(0.45,0.5,0.25,1.25));
#par(mai=c(0.45,0.5,0.25,0.25));
par(mai=c(0.3,0.3,0.0,0.05));



#Seqbeta = seq(from=0.101,to=0.15,by=0.001);

#Lower right

if(0){ # Number 1
  betar = 3.25; mur = 1.0; 
  betai1 = 4.0; mui1 = 1.5; 
  #Wr = 1e-3;
  #Wi = 2.4e-4; 
  #Wr = 1e-5;
  #Wi = 2.5e-6;
  Wr = 1e-7;
  Wi = 2e-8;
  LastGen = 1000;
  IiStart = 1e-10;
}


if(0){ # Number 2
  betar = 4.0; mur = 1.0; 
  betai1 = 3.0; mui1 = 1.0; 
  Wr = 1e-6;
  Wi = 8.5e-5; 
  LastGen = 1000;
  IiStart = 1e-8;
  EndTime = 15;

}

if(1){ # Number 2
  betar = 3.5; mur = 1.0; 
  betai1 = 4; mui1 = 1.3; 
  Wr = 1.1e-6;
  Wi = 4.05e-7; 
  LastGen = 5000;
  IiStart = 1e-8;
  EndTime = 20;

}


#Upper left



u1 = mui1/mur;  


InvadeTime = 1;


Steps  = 1e-2; 

t = seq(from=0,to=EndTime,by=Steps);
R0r = betar/mur; R0i1 = betai1/mui1; 
lambdai1 = betai1 - mui1; u1 = mui1/mur;
lambdar = betar - mur;

Low = 1e-10;
OneStrain = uniroot(f=OneStrainEqm,interval=c(Low,1.0),W=Wr,R0=R0r);
FractInf = uniroot(f=Burnout,interval=c(0.0,1.0),I0=Low,R0=R0r);
Ir = Wr*OneStrain$root;
Ii = 0;

SStart = 1.0-Ir-Ii;

Rptd = R0i1*(1-FractInf$root)*SStart;
cat("R0r:",R0r,"R0i1:",R0i1,"Rptd:",Rptd,"lambdar:",lambdar,"lambdai1:",lambdai1,"mur:",mur,"mui1:",mui1,"\n")
cat("betar:",betar,"mur:",mur,"betai1:",betai1,"mui1:",mui1,"\n")



IrStor = c();
IiStor = c();
IrStorPoints = c();
IiStorPoints = c();
genStorPoints = c();

IrStor[1] = OneStrain$root;
IiStor[1] = 0;


genPoints = 1;
for(gen in 1:LastGen){
  
  if(gen==InvadeTime){
    Ii = IiStart;
  }else{
    if(gen<InvadeTime) Ii = 0;
  }
  

 

  
y0 = c(SStart,Ir,Ii,0.0,0.0);

p = list(betar=betar,betai=betai1,mui=mui1,mur=mur);
out = ode(y=y0,times=t,func=TwoPathEpidVer2,parms=p);

Index = nrow(out);
if(gen>1){
  IrStor[gen] = out[Index,5];
  IiStor[gen] = out[Index,6];
}

#cat("gen:",gen,"\n")
#if(((gen==60)||(gen==80))||(gen==120)){
#if(((gen==100)||(gen==250))||(gen==500)){
if(((gen==1000)||(gen==2000))||(gen==3000)){

  if(genPoints==1){ pdf("EpiDynTransient1.pdf", height = 6, width = 6); 
    cat("genPoints:",genPoints,"\n")
  }
  if(genPoints==2){ 
    pdf("EpiDynTransient2.pdf", height = 6, width = 6);
    cat("genPoints:",genPoints,"\n")
  }
  if(genPoints==3){
    cat("genPoints:",genPoints,"\n")
    pdf("EpiDynTransient3.pdf", height = 6, width = 6);
  }
  
  plot(out[,1],out[,3],ylim=c(0,0.3),type="l",axes=FALSE,xlab="",ylab="",lwd=3.5);
  lines(out[,1],out[,4],col="RED",lwd=3.5) 
  axis(1,pretty(range(out[,1]),5),pos=0.0,cex=1.5);
  
  Ylbls = c(0,0.1,0.2,0.3)
  axis(2,at=Ylbls,pos=0.0,cex=1.5);
  #axis(2,pretty(range(out[,3]),5),pos=0.0);
  dev.off();
  par(ask="FALSE");
  IrStorPoints[genPoints] = IrStor[gen];
  IiStorPoints[genPoints] = IiStor[gen];
  genStorPoints[genPoints] = gen;
  genPoints = genPoints + 1;

}

Ir = Wr*out[Index,5];
Ii = Wi*out[Index,6];
if((gen%%50)==0) cat("gen:",gen,"Ir:",IrStor[gen],"Ii:",IiStor[gen],"\n")
SStart = 1 - Ir - Ii;

}

lgnd = c();
lgnd[1] = expression(beta[r]~paste("=")~3);
lgnd[2] = expression(mu[r]~paste("=")~1);
lgnd[3] = expression(R[0][i]~paste("=")~6);
lgnd[4] = expression(lambda[i]~paste("=")~0.05);
lgnd[5] = expression(R[0][i]~paste("=")~5);
lgnd[6] = expression(lambda[i]~paste("=")~.033);
lgnd[7] = expression(I[r](0)~paste("=")~0.025);
lgnd[8] = expression(I[i](0)~paste("=")~10^-5);

pdf("LongTermTransientb.pdf", height = 4, width = 6)
plot(1:LastGen,IrStor,ylim=c(0,1),type="l",axes=FALSE,xlab="",ylab="",lwd=2.25)
points(genStorPoints,IiStorPoints,col="RED",pch=19,cex=1.25)
points(genStorPoints,IrStorPoints,col="BLACK",pch=19,cex=1.25)
lines(1:LastGen,IiStor,type="l",col="RED",lwd=2.25);
axis(1,pretty(range(1:LastGen),5),pos=0.0);
Ylbls = c(0,0.5,1.0)
axis(2,at=Ylbls,pos=0.0,);
#axis(2,pretty(range(0:1),5),pos=0.0);
dev.off()

#text(x=100,y=0.7,labels="Resident",bty="n",cex=0.85,col="BLACK")
#text(x=100,y=0.6,labels=lgnd[1],bty="n",cex=0.85,col="BLACK");

#dev.off();
