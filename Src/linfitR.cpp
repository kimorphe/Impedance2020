#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spline.h"
#include <complex.h>

using namespace std;

int argmax(double *y, int n){
	int i,imax=0;
	double ymax=y[0];

	for(i=1;i<n;i++){
		if(y[i]> ymax){
			imax=i;
			ymax=y[i];
		}	
	};
	return(imax);
};
int argmin(double *y, int n){
	int i,imin=0;
	double ymin=y[0];

	for(i=1;i<n;i++){
		if(y[i]<ymin){
			imin=i;
			ymin=y[i];
		}	
	};
	return(imin);
};

class IMP{
	public:
		int ndat;
		double *X;
		double *Y;
		double *freq;
		int imax,imin;
		void set_ylim();
		void load(char *fname);
		void Write_Bode(char *fname, bool flip);
		void  RC(double R0, double Rc, double C, double p);
		complex<double> get_Zw(int i);
		IMP();
		IMP(int n);
	private:
};
IMP::IMP(){
};
IMP::IMP(int n){
	ndat=n;
	X=(double *)malloc(sizeof(double)*ndat);
	Y=(double *)malloc(sizeof(double)*ndat);
	freq=(double *)malloc(sizeof(double)*ndat);
};
complex<double> IMP::get_Zw(int i){
	return(complex<double>(X[i],Y[i]));
};
void IMP::RC(double R0, double Rc, double C,double p){
	double RC=Rc*C;
	complex<double> zi,Z;

	zi=complex<double>(0.0,1.0);
	int i;
	double w;
	double PI=4.0*atan(1.0);
	//double p=0.825;
	for(i=0;i<ndat;i++){
		w=2.*PI*freq[i];
		//Z=pow(1.+zi*w*RC,p);
		Z=1.+pow(zi*w*RC,p);
		Z=Rc/Z;
		X[i]=Z.real()+R0;
		Y[i]=Z.imag();
	}
};

void IMP::load(char *fname){
	FILE *fp=fopen(fname,"r");
	
	if(fp==NULL){
		printf("File not found !!\n");
		printf("  --> abort pocess\n");
		printf("%s",fname);
		exit(-1);
	};
	char cbff[256];

	ndat=0;
	while(fgets(cbff,128,fp)!=NULL) ndat++; 
	ndat--;
	fclose(fp);

	X=(double *)malloc(sizeof(double)*ndat);
	Y=(double *)malloc(sizeof(double)*ndat);
	freq=(double *)malloc(sizeof(double)*ndat);

	fp=fopen(fname,"r");
	int i,num;
	fgets(cbff,128,fp);
	for(i=0;i<ndat;i++){
		fscanf(fp,"%d, %lf, %lf, %lf\n",&num,freq+i,X+i,Y+i);
		//Y[i]*=-1;
	};
	fclose(fp);
};
void IMP::Write_Bode(char *fname,bool flip){
	FILE *fp=fopen(fname,"w");
	int i;
	double A,phi;
	static double PI=4.0*atan(1.0);
	int isgn;
	for(i=0;i<ndat;i++){
		A=sqrt(X[i]*X[i]+Y[i]*Y[i]);
		phi=acos(X[i]/A);
		if(Y[i]<0.0) phi*=-1;
		phi/=PI;
		phi+=180.0;
		isgn=1;
		if(flip) isgn=-1;
		fprintf(fp,"%lf %lf %lf %lf %lf\n",log10(freq[i]),log10(A),phi,X[i],isgn*Y[i]);
	}
	fclose(fp);
};
void IMP::set_ylim(){
	imax=argmax(Y,ndat);
	imin=argmin(Y,ndat);
	//printf("imax=%d, imin=%d\n",imax,imin);
	//printf("fmax=%lf, fmin=%lf\n",freq[imax],freq[imin]);
	//imax*=0.5;
};

double L2( double *y1, double *y2, int i1, int i2){
	int i;
	double dy,res=0.0;
	for(i=i1; i<i2; i++){
		dy=y2[i]-y1[i];
		res+=(dy*dy);
	}
	return(res);
};

complex<double> gw(double freq,double p, double RC){
	complex<double> zi,zw;
	zi=complex<double>(0.0,1.0);
	static double PI=4.0*atan(1.0);
	double omg=2.0*PI*freq;
	zw=1.+pow(zi*omg*RC,p);
	return(zw);	
};

double p_search(IMP Z, double p1, double p2, int npnt, double RC, double prms[4],bool verb){
	double dp=(p2-p1)/(npnt-1);
	double a12,a22,b1,b2;
	double A[2][2],b[2];
	complex<double> zg,zw;
	double g2,det;
	double rmin,res;
	int nn,imin;

	double R0,Rc,C,p;
	double Rc_min,R0_min;
	static IMP Z_RC(Z.ndat);
	Z_RC.freq=Z.freq;

	int itr,i;
	Z.set_ylim();
	for(itr=0;itr<npnt;itr++){
		p=p1+dp*itr;
		//RC=RC1+dRC*itr;
		a12=0.0;
		a22=0.0;
		b1=0.0;
		b2=0.0;
		for(i=0;i<Z.imax+1;i++){
			zg=gw(Z.freq[i],p,RC);
			g2=abs(zg*zg);
			a22+=(1./g2);
			a12+=(zg.real()/g2);
			b1+=Z.X[i];
			zw=Z.get_Zw(i);
			b2+=(real((zw)*zg)/g2);
		}

		nn=Z.imax+1;
		A[0][0]=1.0;
		A[0][1]=a12/nn;
		A[1][0]=A[0][1];
		A[1][1]=a22/nn;
		b[0]=b1/nn;
		b[1]=b2/nn;

		det=A[0][0]*A[1][1]-A[0][1]*A[1][0];
		R0=( A[1][1]*b[0]-A[0][1]*b[1])/det;
		Rc=(-A[1][0]*b[0]+A[0][0]*b[1])/det;

		C=RC/Rc;
		Z_RC.RC(R0,Rc,C,p);
		res=sqrt(L2(Z.X, Z_RC.X,0,Z.imax)+L2(Z.Y,Z_RC.Y,0,Z.imax))/nn;
		if(verb) printf("R0,Rc,C,p,res=%lf %lf %lf %lf %lf\n",R0,Rc,C,p,res);
		if(itr==0){
			rmin=res;
			Rc_min=Rc;
			R0_min=R0;
			imin=0;
		}
		if(res <rmin){
			rmin=res; 
			Rc_min=Rc;
			R0_min=R0;
			imin=itr;
		}
	}
	p=p1+dp*imin;
	C=RC/Rc_min;
	Z_RC.RC(R0_min,Rc_min,C,p);
	prms[0]=R0_min;
	prms[1]=Rc_min;
	prms[2]=C;
	prms[3]=p;
	return(rmin);
};
double C_search(IMP Z, double RC1, double RC2, int npnt, double p, double prms[4], bool verb){
	double C;
	int itr,i;
	double RC;
	double dRC=(RC2-RC1)/(npnt-1);
	double a12,a22,b1,b2;
	double A[2][2],b[2];
	complex<double> zg,zw;
	double g2,det;
	double rmin,res;
	int nn,imin;

	double R0,Rc;
	double Rc_min,R0_min;
	static IMP Z_RC(Z.ndat);
	Z_RC.freq=Z.freq;

	double m1=log10(RC1);
	double m2=log10(RC2);
	double dm=(m2-m1)/(npnt-1);
	double mi,RCm;

	for(itr=0;itr<npnt;itr++){	
		mi=m1+dm*itr;
		RC=pow(10,mi);
		//RC=RC1+dRC*itr;
		a12=0.0;
		a22=0.0;
		b1=0.0;
		b2=0.0;
		for(i=0;i<Z.imax+1;i++){
			zg=gw(Z.freq[i],p,RC);
			g2=abs(zg*zg);
			a22+=(1./g2);
			a12+=(zg.real()/g2);
			b1+=Z.X[i];
			zw=Z.get_Zw(i);
			b2+=(real((zw)*zg)/g2);
		}
		
		nn=Z.imax+1;
		A[0][0]=1.0;
		A[0][1]=a12/nn;
		A[1][0]=A[0][1];
		A[1][1]=a22/nn;
		b[0]=b1/nn;
		b[1]=b2/nn;

		det=A[0][0]*A[1][1]-A[0][1]*A[1][0];
		R0=( A[1][1]*b[0]-A[0][1]*b[1])/det;
		Rc=(-A[1][0]*b[0]+A[0][0]*b[1])/det;

		C=RC/Rc;
		Z_RC.RC(R0,Rc,C,p);

		res=sqrt(L2(Z.X, Z_RC.X,0,Z.imax)+L2(Z.Y,Z_RC.Y,0,Z.imax))/nn;
		if(verb) printf("R0,Rc,C,p,res=%lf %lf %lf %lf %lf\n",R0,Rc,C,p,res);
		if(itr==0){
			imin=0;
			rmin=res;
			Rc_min=Rc;
			R0_min=R0;
		}
		if(res <rmin){
			rmin=res; 
			Rc_min=Rc;
			R0_min=R0;
			imin=itr;
		}
	}

	mi=m1+dm*imin;
	RC=pow(10,mi);
	//RC=RC1+dRC*imin;
	C=RC/Rc_min;
	//Z_RC.RC(R0_min,Rc_min,C,p);
	prms[0]=R0_min;
	prms[1]=Rc_min;
	prms[2]=C;
	prms[3]=p;
	return(rmin);
};

int main(){
	int itr,n_itr,j;
	double R0,Rc,C,p,RC;
	double Rc_init,C_init;
	double prms[4],prms_min[4];

	bool verb=false,flip=true;
	double p1,p2;
	double RC0,RC1,RC2; 
	double res,rmin;
	int n_try;
	double alph,beta=0.8,rt=0.95;
	double p0=0.75,pw=0.5,prt=0.9;
	//Rc=300; C=3.0e-02;
	//p1=0.5; p2=1.0; n_try=21;

	char fname[128];
	int ifile,nfile,nchar;
	char cbff[128];
	FILE *fin=fopen("linfitR.inp","r");
	FILE *fpr=fopen("params.dat","w");
	char fnm[128],fnbode[128],fnfit[128];

	fgets(cbff,128,fin);
		fscanf(fin,"%lf %lf\n",&Rc_init,&C_init);
		printf("Rc,C=%lf %lf <-- initial guess)\n",Rc_init,C_init);
	fgets(cbff,128,fin);
		fscanf(fin,"%d\n",&nfile);
		printf("nfile=%d\n",nfile);
	fgets(cbff,128,fin);

	for(ifile=0;ifile<nfile;ifile++){// START_FILE
		fscanf(fin,"%s\n",fname);

	puts(fname);
	nchar=int(strlen(fname));
	strncpy(fnm,fname,nchar-4);
	sprintf(fnbode,"%s%s",fnm,".bode");// output file name (measured)
	sprintf(fnfit,"%s%s",fnm,".fit");  // output file name (fitted) 
	puts(fnbode);
	puts(fnfit);


	IMP Z;
	Z.load(fname);		// load measured data
	Z.Write_Bode(fnbode,flip);	// write measured data
	Z.set_ylim();	// set cutoff frequency

	IMP Z_RC(Z.ndat);	// IMP class for curve fitting
	Z_RC.freq=Z.freq;	// share frequency axis

	//alph,beta=0.8,rt=0.95;
	p0=0.7,pw=0.6,prt=0.95;

	RC0=Rc_init*C_init;
	n_itr=101;
	n_try=51;
	RC0=1.0;
	//RC1=RC0*1.e-06;
	//RC2=RC0*1.e+06;
	double  mexp=4.0;
	for(itr=0;itr<n_itr;itr++){
		p1=p0-pw*0.5;
		p2=p0+pw*0.5;;
		res=p_search(Z,p1,p2,n_try,RC0,prms,verb);
		printf("  p-optimized: %lf %lf %lf %lf (r=%lf)\n",prms[0],prms[1],prms[2],prms[3],res);
		p0=prms[3];
		pw*=prt;
	
		//alph=1.0-beta; //printf("alph=%lf\n",alph);
		//RC1=RC0*alph;
		//RC2=RC0/alph;
		RC1=RC0*pow(10.,-mexp);
		RC2=RC0*pow(10., mexp);
		res=C_search(Z,RC1,RC2,n_try,p0,prms,verb);
		C=prms[3];
		printf("  C-optimized: %lf %lf %lf %lf (r=%lf)\n",prms[0],prms[1],prms[2],prms[3],res);
		RC0=prms[1]*prms[2]; // prms={R0,Rc,C,p}
		//beta*=rt;	
		mexp*=0.95;

		if(itr==0||res<rmin){
			rmin=res;
			for(j=0;j<4;j++) prms_min[j]=prms[j];
		};
	}
	puts("------OPTIMUM-----------");
	printf("%lf %lf %lf %lf (r=%lf)\n",prms_min[0],prms_min[1],prms_min[2],prms_min[3],rmin);
	fprintf(fpr,"%lf %lf %lf %lf %lf\n",prms_min[0],prms_min[1],prms_min[2],prms_min[3],rmin);
	puts("------------------------");
	puts("------------------------");

	Z_RC.RC(prms_min[0],prms_min[1],prms_min[2],prms_min[3]);
	Z_RC.Write_Bode(fnfit,flip); // write fitted curve

	} //END_FILE


	return(0);
};
