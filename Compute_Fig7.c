#include <stdio.h>
#include <math.h>

FILE *fopen(), *fp;
void main()
{ double rho1a,rho1b,rho2a,rho2b,dum1,dum2,res1a,res1b,dum3,dum4,dum2b;
  int flag, curve_type;
  /* curve_type = 1 means left curve (A)
     curve_type = 2 means right hand curve (B)*/
  flag =55;
  curve_type = 2;
  fp= fopen("MyLog11FIGextra2.dat","r");
  rho1b=20.0; rho2b=0.0;
  if(curve_type ==2) {res1b=1.0; dum2b=-1.0;}
  if(curve_type ==1) {res1b=-1.0; dum2b=-1.0;}
  while( (flag=fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",
		     &rho1a,&rho2a,&dum1,&dum2,&res1a,&dum3))!=EOF)
    { 
      if(flag==6 && curve_type ==2)
      {  
	if(res1b<0.0&&res1a>0.0&&rho1a==rho1b && fabs(dum2b+1.0)>0.00001)
	  printf("%f  %f %f\n",rho1a,(rho2a*res1b-rho2b*res1a)/(res1b-res1a),dum2b);
       res1b=res1a; rho1b=rho1a,rho2b=rho2a;
      }
      if(flag==6 && curve_type ==1)
      {  
	if(res1a<0.0&&res1b>0.0&&rho1a==rho1b && fabs(dum2b+1.0)>0.00001)
	  printf("%f  %f %f\n",rho1a,(rho2a*res1b-rho2b*res1a)/(res1b-res1a),dum2b);
       res1b=res1a; rho1b=rho1a,rho2b=rho2a;
      }



      dum2b=dum2;
    }
  fclose(fp);
}
