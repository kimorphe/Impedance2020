import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import datafiles as dtf

class IMP:
    def load(self,fname):
        #dat=pd.read_csv(fname,skiprows=0)
        dat=pd.read_csv(fname,skiprows=0)
        self.R=np.array(dat['Rs'])
        self.X=np.array(dat['X'])
        self.f=np.array(dat['FREQUENCY(Hz)'])
        self.fname=fname
        self.icut=-1

    def plot(self,ax,name=""):
    	if name=="":
    		name=self.fname
    	ax.plot(self.R,-self.X,"o",linewidth=1,label=name)
    def Nyquist(self,ax,cut=False,name=""):
        imax=-1
        if cut:
            imax=self.icut
        ax.plot(self.R[0:imax],-self.X[0:imax],"-",linewidth=1,label=name)

    def fplot(self,ax1,ax2,name=""):
    	if name=="":
    		name=self.fname
    	ax1.semilogx(self.f,self.R,label=name)
    	ax2.semilogx(self.f,self.X,label=name)

    def diff(self):
        self.dX=np.diff(self.R)
        self.dY=np.diff(self.X)
    def Bode(self,ax1,ax2,name="",cut=False):
        if name=="":
            name=self.fname
        Z=self.R+1j*self.X
        phi=np.angle(Z)
        phi=np.degrees(phi)

        imax=len(self.f)
        #print(imax)
        if cut==True:
            imax=self.icut
        ax1.loglog(self.f[0:imax],np.abs(Z[0:imax]),"-",label=name)
        ax2.semilogx(self.f[0:imax],phi[0:imax],"-",label=name)
    def ftrim(self):
        self.diff()
        dYdX=self.dY/self.dX
        indx=dYdX < 0 
        INDX=np.cumsum(indx.astype(int))
        imax=np.argmax(INDX)
        #print(imax)
        #print(self.f[imax])
        self.fcut=self.f[imax]
        self.icut=imax


if __name__=="__main__":
    dir_name="./"
    #date="0128"; name="w15"  
    #date="0201"; name="w22" 
    #date="0126"; name="w20"
    #date="0202"; name="w17"
    date="0129"; name="w25"  

    # dry density-water content plot
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.grid(True)

    fsz=14
    # Nyquis Plot
    fig2=plt.figure();
    bx=fig2.add_subplot(111)
    bx.grid(True)
    bx.set_xlabel(r"Re{Z} [$\Omega$]",fontsize=fsz)
    bx.set_ylabel(r"-Im{Z} [$\Omega$]",fontsize=fsz)
    bx.set_aspect(1.0)

    # Bode  Plot
    fig3=plt.figure()
    bx1=fig3.add_subplot(211)
    bx2=fig3.add_subplot(212)
    bx1.grid(True)
    bx2.grid(True)

    
    Z=IMP() # Measured Impedance 

    fname="datafiles"+date+".csv"
    print(fname)
    FL=dtf.FILE_LIST(dir_name,fname)
    A=FL.shape()

    rho_dry=np.array(FL.get_col("rho_dry"))
    rho_dry=rho_dry.astype(float)
    wt=np.array(FL.get_col("water_content"))
    wt=wt.astype(float)
    ax.plot(wt,rho_dry,"o",markersize=8)
    data_dir=date
    USE=FL.get_col("use")
    FNAME=FL.get_col("name")
    k=0
    for use in USE: # for each data file
        if int(use)==1:
            fname=data_dir+"/"+FNAME[k]+".csv"
            print(fname)
            Z.load(fname)
            #Z.diff(); Z.ftrim(); imax=Z.icut
            Z.Bode(bx1,bx2,cut=True)
            Z.Nyquist(bx,cut=True,name=fname)
        k+=1

    #bx.legend()
    w=np.linspace(0.12,0.25)
    Sr=np.linspace(1.0,0.3,8)
    for sr in Sr:
        rho_d=dtf.Rho_Dry(sr,w)
        ax.plot(w,rho_d,"--k",label="Sr="+str(sr),linewidth=1)

    ax.set_xlabel("water content",fontsize=fsz)
    ax.set_ylabel("dry density [g/cm$^3$]",fontsize=fsz)
    ax.tick_params(labelsize=fsz)
    ax.set_xlim([0.12,0.25])
    ax.set_ylim([1.2,2.0])

    label=name+" series ($w$="+str(wt[0])+"%)"
    bx.set_title(label,fontsize=fsz)
    bx1.set_title(label,fontsize=fsz)
    #bx2.set_title(label,fontsize=fsz)
    bx2.set_xlabel("frequency [Hz]",fontsize=fsz)
    bx1.set_ylabel("|Z|",fontsize=fsz)
    bx2.set_ylabel(r"$\theta$[deg]",fontsize=fsz)
    bx.tick_params(labelsize=fsz)
    bx1.tick_params(labelsize=fsz)
    bx2.tick_params(labelsize=fsz)
    x1,x2=bx.get_xlim()
    bx.set_xlim([x1,x1+500])
    bx.set_ylim([0,300])

    
    fig2.savefig("Nyq"+date+".png",bbox_inches="tight")
    fig3.savefig("Bode"+date+".png",bbox_inches="tight")
    plt.show()

