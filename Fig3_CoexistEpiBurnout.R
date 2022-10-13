
par(mfrow=c(1,1));


require(deSolve);

OneStrainEqm<-function(R0,W,i){
  
  arg = (1 - (W*i))*i + W*i;
  output = 1 - i - exp(-R0 * arg);
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
par(mai=c(0.5,0.5,0.25,0.25));



#Seqbeta = seq(from=0.101,to=0.15,by=0.001);

#Lower right

if(1){ # Number 1
  betar = 1.25; mur = 1.0; 
  betai1 = 1.8; mui1 = 1.0; 
  Wr = 1e-3;
  Wi = 1e-9; #1.5e-10; 
  LastGen = 10;
  IiStart = 1e-15;
}


if(0){ # Number 2
  betar = 3.0; mur = 1.5; 
  betai1 = 4.0; mui1 = 1.0; 
  Wr = 1e-3; 
  Wi = 1e-6; 
  LastGen = 100;
  IiStart = 1e-10;
}

u1 = mui1/mur;  
EndTime = 100;

InvadeTime = 1;


Steps  = 1e-2; 

t = seq(from=0,to=EndTime,by=Steps);
R0r = betar/mur; R0i1 = betai1/mui1; 
lambdai1 = betai1 - mui1; u1 = mui1/mur;
lambdar = betar - mur;

Low = 1e-10;
OneStrain = uniroot(f=OneStrainEqm,interval=c(Low,1.0),W=Wr,R0=R0r);

#FractInf = uniroot(f=Burnout,interval=c(0.0,1.0),I0=Low,R0=R0r);
FractInf = uniroot(f=Burnout,interval=c(0.0,1.0),I0=Wr*OneStrain$root,R0=R0r);

IrEqm = Ir = Wr*OneStrain$root;
cat("IrEqm:",Ir,"\n")
Ii = 0;

SStart = 1.0-Ir-Ii;

#Rptd = R0i1*(1-OneStrain$root);
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
#if(gen>1){
  IrStor[gen] = out[Index,5];
  IiStor[gen] = out[Index,6];
#}

#cat("gen:",gen,"\n")

#if(((gen==10)||(gen==13))||(gen==17)){ #2
if(((gen==2)||(gen==4))||(gen==6)){ #1

  if(genPoints==1){ pdf("EpiDynBurnout1b.pdf", height = 6, width = 6); 
    cat("genPoints:",genPoints,"\n")
  }
  if(genPoints==2){ 
    pdf("EpiDynBurnout2b.pdf", height = 6, width = 6);
    cat("genPoints:",genPoints,"\n")
  }
  if(genPoints==3){
    cat("genPoints:",genPoints,"\n")
    pdf("EpiDynBurnout3b.pdf", height = 6, width = 6);
  }
  
  plot(out[,1],out[,3],ylim=c(0,0.025),type="l",axes=FALSE,xlab="",ylab="",lwd=3.5);
  lines(out[,1],out[,4],col="RED",lwd=2.25) 
  axis(1,pretty(range(out[,1]),5),pos=0.0,cex=3.5);
  
  Ylbls = c(0,0.01,0.02)
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
#cat("gen:",gen,"Ir:",IrStor[gen],"Ii:",IiStor[gen],"\n")
#if((gen%%50)==0) 
  cat("gen:",gen,"Ir:",IrStor[gen],"Ii:",IiStor[gen],"\n")

SStart = 1 - Ir - Ii;

}

pdf("LongTermBurnoutb.pdf", height = 4, width = 6)
plot(1:LastGen,IrStor,ylim=c(0,0.5),type="l",axes=FALSE,xlab="",ylab="",lwd=2.25)
#plot(1:LastGen,IrStor,type="l",ylim=c(1,0),axes=TRUE,xlab="",ylab="")

points(genStorPoints,IiStorPoints,col="RED",pch=19,cex=1.25)
points(genStorPoints,IrStorPoints,col="BLACK",pch=19,cex=1.25)
lines(1:LastGen,IiStor,type="l",col="RED",lwd=2.25);
axis(1,pretty(range(1:LastGen),5),pos=0.0);
Ylbls = c(0,0.25,0.5)
axis(2,at=Ylbls,pos=1.0);
dev.off()

#text(x=100,y=0.7,labels="Resident",bty="n",cex=0.85,col="BLACK")
#text(x=100,y=0.6,labels=lgnd[1],bty="n",cex=0.85,col="BLACK");

#dev.off();
Index = length(IrStor)
cat("IrEqm:",IrEqm,"Ir obsvd:",IrStor[Index],"\n")
