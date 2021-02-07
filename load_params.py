import numpy as np
import matplotlib.pyplot as plt

class PRMS:
    def __init__(self,fname):
        self.fname=fname
        fp=open(fname,"r")
        R0=[]
        Rc=[]
        C=[]
        p=[]
        res=[]
        ndat=0;
        for row in fp:
            data=row.strip()
            data=data.split(" ")
            ndat+=1;
            R0.append(float(data[0]))
            Rc.append(float(data[1]))
            C.append(float(data[2]))
            p.append(float(data[3]))
            res.append(float(data[4]))
        self.R0=np.array(R0)
        self.Rc=np.array(Rc)
        self.C=np.array(C)
        self.p=np.array(p)
        self.res=np.array(res)

    def plot(self,ax,idat):
        fsz=14
        if idat==0:
            ax.plot(self.R0)
            ax.set_ylabel("R0 [$\Omega$]")
        if idat==1:
            ax.plot(self.Rc)
            ax.set_ylabel("Rc [$\Omega$]")
        if idat==2:
            ax.plot(self.C)
            ax.set_ylabel("C [F]")
        if idat==3:
            ax.plot(self.p)
            ax.set_ylabel("p")
        if idat==4:
            ax.plot(self.res)
            ax.set_ylabel("residual")

if __name__=="__main__":
    prm=PRMS("0128/params.dat")
    fig=plt.figure()
    ax=fig.add_subplot(111)
    prm.plot(ax,2)
    ax.grid(True)

    plt.show()


                




