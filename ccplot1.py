import numpy as np
import matplotlib.pyplot as plt
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
    	ax.plot(self.R,-self.X,"-",linewidth=1,label=name)

    def fplot(self,ax1,ax2,name=""):
    	if name=="":
    		name=self.fname
    	ax1.semilogx(self.f,self.R,label=name)
    	ax2.semilogx(self.f,self.X,label=name)

    def diff(self):
        self.dX=np.diff(self.R)
        self.dY=np.diff(self.X)

if __name__=="__main__":


    fig=plt.figure();
    ax=fig.add_subplot(111)
    ax.grid(True)
    ax.set_xlabel(r"Re{Z} [$\Omega$]",fontsize=12)
    ax.set_ylabel(r"Im{Z} [$\Omega$]",fontsize=12)
    ax.set_aspect(1.0)

    fig2=plt.figure()
    bx1=fig2.add_subplot(211)
    bx2=fig2.add_subplot(212)
    bx1.grid(True)
    bx2.grid(True)

    Z=IMP()
    #dir="2kg1V";
    dir="0112/h615";
    n1=5;
    n2=10;
    nums=np.arange(n1,n2+1,1)
    nums=nums.astype(int)
    
    fig3=plt.figure()
    ax2=fig3.add_subplot(111)
    ax2.grid(True)
    for k in nums:
    	fname=dir+"/"+"data"+str(k)+".csv"
    	Z.load(fname)
    	Z.plot(ax,name=str(k))
    	Z.fplot(bx1,bx2);Z.diff();ax2.semilogx(Z.f[1:],Z.dY/Z.dX)

    ax.set_title(dir)
    ax.legend()
    bx1.legend()
    fig.savefig(dir+".png",bbox_inches="tight");
    fig2.savefig(dir+"f.png",bbox_inches="tight");

    plt.show()



