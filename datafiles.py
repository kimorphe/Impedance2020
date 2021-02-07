import numpy as np
import matplotlib.pyplot as plt
import csv


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
            k=0
            k2del=[]
            for line in D:
                if line[0]=="": 
                    k2del.append(k)
                k+=1
            #print(k2del)
            del(D[k2del[0]:k2del[-1]+1])
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
        #print("Number of keys=",n1)
        #print("Number of lines=",n2)
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


    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.grid(True)

    for date in dates:
        fname="datafiles"+date+".csv"
        print(fname)
        FL=FILE_LIST(dir_name,fname)
        A=FL.shape()
        #print(FL.get_elem("m_dry",5))
        #print(FL.get_elem("name",5))
        #print(FL.get_line(0))
        print(FL.get_col("ID"))

        rho_dry=np.array(FL.get_col("rho_dry"))
        rho_dry=rho_dry.astype(float)
        wt=np.array(FL.get_col("water_content"))
        wt=wt.astype(float)

        ax.plot(wt,rho_dry,"o",markersize=8)

    w=np.linspace(0.12,0.25)
    Sr=np.linspace(1.0,0.3,8)
    print(Sr)
    for sr in Sr:
        rho_d=Rho_Dry(sr,w)
        ax.plot(w,rho_d,"--k",label="Sr="+str(sr),linewidth=1)

    fsz=14
    ax.set_xlabel("water content",fontsize=fsz)
    ax.set_ylabel("dry density [g/cm$^3$]",fontsize=fsz)
    ax.tick_params(labelsize=fsz)
    ax.set_xlim([0.12,0.25])
    ax.set_ylim([1.2,2.0])

    fig.savefig("rhod.png",bbox_inches="tight")
    plt.show()

