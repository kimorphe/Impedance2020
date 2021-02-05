import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

class IMP:
    def load(self,fname):
        #dat=pd.read_csv(fname,skiprows=0)
        dat=pd.read_csv(fname,skiprows=0)
        self.R=np.array(dat['Rs'])
        self.X=np.array(dat['X'])
        self.f=np.array(dat['FREQUENCY(Hz)'])
        self.fname=fname

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
        print(imax)
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
        print(imax)
        print(self.f[imax])
        self.fcut=self.f[imax]
        self.icut=imax

def Rho_Dry(Sr,w,rho_s=2.7):
    rho_w=1.0;
    rho_d=Sr/rho_s+w/rho_w
    rho_d=Sr/rho_d
    return(rho_d)


class FILE_LIST:
    def __init__(self,dir_name,fname):
        self.dir_name=dir_name
        self.fname=fname

        with open(dir_name+"/"+fname) as f:
            reader=csv.reader(f)
            # Use comprehension notation
            D=[row for row in reader]   # reads "stack row for each element(row) in reader"
            keys=D[0]
            N=np.shape(D)
            Key_Num=dict(zip(keys,range(N[1])))
            del(D[0])   # remove header line
            del(D[0])   # remove header line
            del(D[-1])   # remove header line
            N=np.shape(D)
            nkeys=N[1]
            nfile=N[0]

        self.D=D
        self.nkeys=nkeys
        self.nfile=nfile
        self.Key_Num=Key_Num
    def shape(self):
        n1=self.nkeys
        n2=self.nfile
        print("Number of keys=",n1)
        print("Number of lines=",n2)
        return([n1,n2])
    def get_elem(self,key,line):
        icol=self.Key_Num[key]
        L=self.D[line]
        return(L[icol])
    def get_line(self,line):
        return(self.D[line])
    def get_col(self,key):
        icol=self.Key_Num[key]
        data=[]
        for k in range(self.nfile):
            L=self.D[k]
            data.append(L[icol])
        return(data)

if __name__=="__main__":
    dir_name="./"
    dates=["0126","0128","0129","0201","0202"]

    # dry density-water content plot
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.grid(True)

    # Nyquis Plot
    fig2=plt.figure();
    bx=fig2.add_subplot(111)
    bx.grid(True)
    bx.set_xlabel(r"Re{Z} [$\Omega$]",fontsize=12)
    bx.set_ylabel(r"Im{Z} [$\Omega$]",fontsize=12)
    bx.set_aspect(1.0)

    # Bode  Plot
    fig3=plt.figure()
    bx1=fig3.add_subplot(211)
    bx2=fig3.add_subplot(212)
    bx1.grid(True)
    bx2.grid(True)

    
    Z=IMP() # Measured Impedance 
    for date in dates:  # for each data folder
        fname="datafiles"+date+".csv"
        print(fname)
        FL=FILE_LIST(dir_name,fname)
        A=FL.shape()
        print(FL.get_col("ID"))

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
            if use:
                fname=data_dir+"/"+FNAME[k]+".csv"
                print(fname)
                Z.load(fname)
                Z.diff()
                Z.ftrim()
                Z.Bode(bx1,bx2,cut=True)
                Z.Nyquist(bx,cut=True,name=fname)
                imax=Z.icut
            k+=1

    #bx.legend()
    w=np.linspace(0.12,0.25)
    Sr=np.linspace(1.0,0.3,8)
    for sr in Sr:
        rho_d=Rho_Dry(sr,w)
        ax.plot(w,rho_d,"--k",label="Sr="+str(sr),linewidth=1)

    fsz=14
    ax.set_xlabel("water content",fontsize=fsz)
    ax.set_ylabel("dry density [g/cm$^3$]",fontsize=fsz)
    ax.tick_params(labelsize=fsz)
    ax.set_xlim([0.12,0.25])
    ax.set_ylim([1.2,2.0])

    fig2.savefig("Nyquist.png",bbox_inches="tight")
    fig3.savefig("Bode.png",bbox_inches="tight")
    plt.show()

