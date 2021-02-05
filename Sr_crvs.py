import numpy as np
import matplotlib.pyplot as plt


def Rho_Dry(Sr,w,rho_s=2.7):
    rho_w=1.0;
    rho_d=Sr/rho_s+w/rho_w
    rho_d=Sr/rho_d
    return(rho_d)

if __name__=="__main__":

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.grid(True)

    w=np.linspace(0.15,0.25)
    Sr=np.linspace(1.0,0.5,6)
    for sr in Sr:
        rho_d=Rho_Dry(sr,w)
        ax.plot(w,rho_d,label="Sr="+str(sr))

    ax.set_xlim([w[0],w[-1]])
    ax.set_ylim([1.0,2.0])

    fsz=14
    ax.tick_params(labelsize=fsz)
    ax.set_xlabel("water content",fontsize=fsz)
    ax.set_ylabel("dry density [g/cm$^3$]",fontsize=fsz)
    ax.legend()
    plt.show()


