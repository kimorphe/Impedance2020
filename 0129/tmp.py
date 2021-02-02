import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#   CIRCUIT ELEMENTS
class Cp:   # Capacitor 
    def setup(self,C,p=1): #(capacitance, exponent of CPE:constant phase element) 
        self.C=C
        self.p=p
    def set_omg(self,f1,f2,nf=100):
        m1=np.log10(f1)
        m2=np.log10(f2)
        m=np.linspace(m1,m2,nf)
        self.freq=10**m;
        self.w=2.*np.pi*self.freq
        return(self.w)
    def eval(self):
        self.Z=self.C*(1j*self.w)**(self.p)
        self.Z=1./self.Z
    def eval_w(self,w):
        self.w=w
        self.freq=w/np.pi*0.5
        return(self.Z)
    def show(self,ax):
        X=np.real(self.Z)
        Y=np.imag(self.Z)
        ax.plot(X,-Y,"-s")

class Wb:   # Warburg Impedance of INFINITE diffusion layer
    def setup(self,sig):
        self.sig=sig
    def set_omg(self,f1,f2,nf=100):
        m1=np.log10(f1)
        m2=np.log10(f2)
        m=np.linspace(m1,m2,nf)
        self.freq=10**m;
        self.w=2.*np.pi*self.freq
        return(self.w)
    def eval(self):
        self.Z=(1-1j)/np.sqrt(self.w)*self.sig;
    def eval_w(self,w):
        self.w=w
        self.freq=w/np.pi*0.5
        self.eval()
    def show(self,ax):
        X=np.real(self.Z)
        Y=np.imag(self.Z)
        ax.plot(X,-Y,"-o")

class WbF: # Warburg impedance of FINITE diffusion layer
    def setup(self,A0,W):
        self.A0=A0  # Amplitude scaling factor
        self.W=W;   # Frequency scaling factor
    def set_omg(self,f1,f2,nf=100): # set frequencies
        m1=np.log10(f1)
        m2=np.log10(f2)
        m=np.linspace(m1,m2,nf)
        self.freq=10**m;
        self.w=2.*np.pi*self.freq
        return(self.w)
    def eval(self): # # Evaluate Impledance
        wb=np.sqrt(self.w*1j*self.W)
        self.Z=self.A0*np.tanh(wb)/wb;
    def eval_w(self,w): 
        self.w=w
        self.freq=w/np.pi*0.5
        self.eval()
    def show(self,ax):
        X=np.real(self.Z)
        Y=np.imag(self.Z)
        ax.plot(X,-Y,"-o")

class ZCD: # Cole-Davidson Impedance
    def setup(self,Rsol,R,C,p):
        self.Rsol=Rsol
        self.R=R
        self.C=C
        self.p=p
    def set_omg(self,f1,f2,nf=100): # set frequencies
        m1=np.log10(f1)
        m2=np.log10(f2)
        m=np.linspace(m1,m2,nf)
        self.freq=10**m;
        self.w=2.*np.pi*self.freq
        return(self.w)
    def eval(self): # # Evaluate Impledance
        #Z=(1+1j*self.R*self.C*self.w)**self.p
        Z=(1+1j*self.C*self.w)**self.p
        Z=self.R/Z+self.Rsol
        #Z=(1+(1j*self.R*self.C*self.w)**self.p)
        #Z=self.R/Z+self.Rsol
        self.Z=Z
    def eval_w(self,w):
        self.w=w
        self.freq=w/np.pi*0.5
        self.eval()
    def Nyquist(self,ax):
        X=np.real(self.Z)
        Y=np.imag(self.Z)
        ax.plot(X,-Y,"-")
    def Bode(self,ax1,ax2,name=""):
        phi=np.degrees(np.angle(self.Z))
        ax1.loglog(self.freq,np.abs(self.Z),label=name)
        ax2.semilogx(self.freq,phi,label=name)

def para(Z1,Z2):
    Z=1./Z1+1./Z2
    return(1/Z)

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
    dir="./";
    n1=18;
    n2=18;
    nums=np.arange(n1,n2+1,5)
    nums=nums.astype(int)

    nums=[6,11,14,20]
    

    for k in nums:
        fname=dir+"/"+"data"+str(k)+".csv"
        Z.load(fname)
        #Z.R-=175
        #Z.plot(ax,name=str(k)); 
        #Z.fplot(bx1,bx2);Z.diff();
        #ax2.semilogx(Z.f[1:],Z.dY/Z.dX)
        Z.diff()
        Z.ftrim()
        Z.Bode(bx1,bx2,cut=True);
        Z.Nyquist(ax,cut=True)
        imax=Z.icut

    plt.show()

