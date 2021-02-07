#include <stdio.h>
#include <stdlib.h>
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
		void Write_Bode(char *fname);
		void  RC(double R0, double Rc, double C);
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
void IMP::RC(double R0, double Rc, double C){
	double RC=Rc*C;
	complex<double> zi,Z;

	zi=complex<double>(0.0,1.0);
	int i;
	double w;
	double PI=4.0*atan(1.0);
	double p=0.825;
	for(i=0;i<ndat;i++){
		w=2.*PI*freq[i];
		Z=pow(1.+zi*w*RC,p);
		Z=Rc/Z;
		X[i]=Z.real()+R0;
		Y[i]=-Z.imag();
	}
};

void IMP::load(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[256];

	ndat=0;
	while(fgets(cbff,128,fp)!=NULL) ndat++; 
	ndat--;
	printf("# ndat=%d\n",ndat);
	fclose(fp);

	X=(double *)malloc(sizeof(double)*ndat);
	Y=(double *)malloc(sizeof(double)*ndat);
	freq=(double *)malloc(sizeof(double)*ndat);

	fp=fopen(fname,"r");
	int i,num;
	fgets(cbff,128,fp);
	for(i=0;i<ndat;i++){
		fscanf(fp,"%d, %lf, %lf, %lf\n",&num,freq+i,X+i,Y+i);
		Y[i]*=-1;
	};
	fclose(fp);
};
void IMP::Write_Bode(char *fname){
	FILE *fp=fopen(fname,"w");
	int i;
	double A,phi;
	static double PI=4.0*atan(1.0);
	for(i=0;i<ndat;i++){
		A=sqrt(X[i]*X[i]+Y[i]*Y[i]);
		phi=acos(X[i]/A);
		if(Y[i]<0.0) phi*=-1;
		phi/=PI;
		phi+=180.0;
		fprintf(fp,"%lf %lf %lf %lf %lf\n",log10(freq[i]),log10(A),phi,X[i],Y[i]);
	}
	fclose(fp);
};
void IMP::set_ylim(){
	imax=argmax(Y,ndat);
	imin=argmin(Y,ndat);
	printf("imax=%d, imin=%d\n",imax,imin);
	printf("fmax=%lf, fmin=%lf\n",freq[imax],freq[imin]);
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
int main(){
	char fname[128]="../0129/data16.csv";
	IMP Z;
	Z.load(fname);	
	sprintf(fname,"bode.dat");
	Z.Write_Bode(fname);

	IMP Z_RC(Z.ndat);
	Z_RC.freq=Z.freq;
	double R0=175;
	double Rc=330;
	double C=1.8e-02;


	Z.set_ylim();
	int i,imin=0;
	double rmin,res;
	for(i=0;i<20;i++){
		R0=175-20+5*i;
		Z_RC.RC(R0,Rc,C);
		res=L2(Z.X, Z_RC.X,0,Z.imin)+L2(Z.Y,Z_RC.Y,0,Z.imin);
		if(i==0) rmin=res;
		if(res <rmin){
			rmin=res; 
			imin=i;
		}
		printf("%lf\n", res);
	}
	R0=175-20+5*imin;
	Z_RC.RC(R0,Rc,C);
	sprintf(fname,"rc.dat");
	Z_RC.Write_Bode(fname);

	Curve2D cv0,cv1;
	cv0.init(Z.ndat);
	cv1.init(Z.ndat);

	cv0.set_xy(Z.X, Z.Y, Z.ndat);
	cv0.smooth(5);
	cv0.spline();

	//int imax=argmax(Z.Y,Z.ndat);
	//int imin=argmin(Z.Y,Z.ndat);

	double dx,dy,dz,th;
	double PI=atan(1.0)*4.0;
	FILE *fout=fopen("grad.out","w");
	for(i=0;i<Z.ndat;i++){
		dx=cv0.dxds(i);
		dy=cv0.dyds(i);
		dz=sqrt(dx*dx+dy*dy);
		th=acos(dx/dz);
		if(dy<0.0) th=2.*PI-th;
		//if(dy<0.0) th=-th;
		th/=PI;
		th*=180.0;
		fprintf(fout,"%lf %lf %lf %lf\n",log10(Z.freq[i]), dx,dy/dx,th);
		cv1.x[i]=dx;
		cv1.y[i]=dy;
	}

	return(0);
};
