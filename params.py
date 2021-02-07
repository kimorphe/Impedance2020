import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import load_params as PRM
import datafiles as dtf

if __name__=="__main__":
    dir_name="./"
    dates=["0126","0128","0129","0201","0202"]

    # dry density-water content plot
    fig1=plt.figure()
    ax0=fig1.add_subplot(111)
    ax0.grid(True)

    fgsz=[9,5]
    # C-plot
    fig2=plt.figure(figsize=fgsz);
    ax1=fig2.add_subplot(121)
    ax2=fig2.add_subplot(122)
    ax1.grid(True)
    ax2.grid(True)

    # p-plot
    fig3=plt.figure(figsize=fgsz);
    bx1=fig3.add_subplot(121)
    bx2=fig3.add_subplot(122)
    bx1.grid(True)
    bx2.grid(True)

    # Rc-plot
    fig4=plt.figure(figsize=fgsz);
    cx1=fig4.add_subplot(121)
    cx2=fig4.add_subplot(122)
    cx1.grid(True)
    cx2.grid(True)

    # R0-plot
    fig5=plt.figure(figsize=fgsz);
    dx1=fig5.add_subplot(121)
    dx2=fig5.add_subplot(122)
    dx1.grid(True)
    dx2.grid(True)

    # residual plot
    fig6=plt.figure(figsize=fgsz);
    ex1=fig6.add_subplot(121)
    ex2=fig6.add_subplot(122)
    ex1.grid(True)
    ex2.grid(True)

    # T_CPE plot
    fig7=plt.figure(figsize=fgsz);
    fx1=fig7.add_subplot(121)
    fx2=fig7.add_subplot(122)
    fx1.grid(True)
    fx2.grid(True)


    for date in dates:  # for each data folder
        fname="datafiles"+date+".csv"
        print(fname)
        FL=dtf.FILE_LIST(dir_name,fname)
        A=FL.shape()
        #print(FL.get_col("ID"))

        rho_dry=np.array(FL.get_col("rho_dry"))
        rho_dry=rho_dry.astype(float)
        wt=np.array(FL.get_col("water_content"))
        wt=wt.astype(float)
        ax0.plot(wt*100,rho_dry,"o",markersize=10)
        sat=np.array(FL.get_col("saturation"))
        sat=sat.astype(float)
        diam=np.array(FL.get_col("diameter"))
        diam=diam.astype(float)

        data_dir=date
        USE=FL.get_col("use")
        FNAME=FL.get_col("name")
        k=0
        rhod=[];
        Sr=[];
        pexp=[]
        Cp=[]
        Rc=[]; R0=[]; res=[];
        wct=[]
        area=[]
        prm=PRM.PRMS(date+"/params.dat")
        for use in USE: # for each data file
            if int(use) ==1:
                rhod.append(rho_dry[k])
                Sr.append(sat[k]);
                wct.append(wt[k]);
                pexp.append(prm.p[k])
                Cp.append(prm.C[k])
                Rc.append(prm.Rc[k])
                R0.append(prm.R0[k])
                res.append(prm.res[k])
                D=diam[k]
                area.append(D*D*0.25*np.pi)
            k+=1

        M=1.e06; K=1.e03;
        area=np.array(area)
        wct=np.array(wct)*100
        ax1.semilogy(rhod,Cp/area*M,"o",markersize=12)
        ax2.semilogy(wct,Cp/area*M,"o",markersize=12)
        bx1.plot(rhod,pexp,"o",markersize=12)
        bx2.plot(wct,pexp,"o",markersize=12)
        cx1.plot(rhod,Rc/area,"o",markersize=12)
        cx2.plot(wct,Rc/area,"o",markersize=12)
        dx1.plot(rhod,R0/area,"o",markersize=12)
        dx2.plot(wct,R0/area,"o",markersize=12)
        ex1.plot(rhod,res,"o",markersize=12)
        ex2.plot(wct,res,"o",markersize=12)

        Rc=np.array(Rc)
        Cp=np.array(Cp)
        pexp=np.array(pexp)
        T_CPE=((Rc*Cp)**pexp)/Rc
        fx1.plot(rhod,T_CPE*K,"o",markersize=12)
        fx2.plot(wct,T_CPE*K,"o",markersize=12)


    #bx.legend()
    w=np.linspace(0.12,0.25)
    Sr=np.linspace(1.0,0.3,8)
    for sr in Sr:
        rho_d=dtf.Rho_Dry(sr,w)
        ax0.plot(w*100,rho_d,"--k",label="Sr="+str(sr),linewidth=1)

    fsz=14
    ax0.set_xlabel("water content [$\%$]",fontsize=fsz)
    ax0.set_ylabel("dry density [g/cm$^3$]",fontsize=fsz)
    ax0.set_xlim([12,25])
    ax0.set_ylim([1.2,2.0])

    ax1.set_xlabel("dry density [g/cm$^3$]",fontsize=fsz);
    ax2.set_xlabel("water content [$\%$]",fontsize=fsz);
    ax1.set_ylabel("Capacitance [$\mu$F/mm$^2$]",fontsize=fsz);
    ax1.set_ylim([2.,4.e02])
    ax2.set_ylim([2.,4.e02])

    bx1.set_xlabel("dry density [g/cm$^3$]",fontsize=fsz);
    bx2.set_xlabel("water content [$\%$]",fontsize=fsz);
    bx1.set_ylabel("exponent $p$",fontsize=fsz);

    cx1.set_xlabel("dry density [g/cm$^3$]",fontsize=fsz);
    cx2.set_xlabel("water content [$\%$]",fontsize=fsz);
    cx1.set_ylabel("Resistance $R_{ct}$[$\Omega$/mm$^2$]",fontsize=fsz);

    dx1.set_xlabel("dry density [g/cm$^3$]",fontsize=fsz);
    dx2.set_xlabel("water content [$\%$]",fontsize=fsz);
    dx1.set_ylabel("Resistance $R_0$[$\Omega$/mm$^2$]",fontsize=fsz);

    ex1.set_xlabel("dry density [g/cm$^3$]",fontsize=fsz);
    ex2.set_xlabel("water content [$\%$]",fontsize=fsz);
    ex1.set_ylabel("residual [$\Omega$/data point]",fontsize=fsz);

    fx1.set_xlabel("dry density [g/cm$^3$]",fontsize=fsz);
    fx2.set_xlabel("water content [$\%$]",fontsize=fsz);
    fx1.set_ylabel("CPE Constant [mF/s$^{1-p}$]",fontsize=fsz);

    ax0.tick_params(labelsize=fsz)
    ax1.tick_params(labelsize=fsz)
    ax2.tick_params(labelsize=fsz)
    bx1.tick_params(labelsize=fsz)
    bx2.tick_params(labelsize=fsz)
    cx1.tick_params(labelsize=fsz)
    cx2.tick_params(labelsize=fsz)
    dx1.tick_params(labelsize=fsz)
    dx2.tick_params(labelsize=fsz)
    ex1.tick_params(labelsize=fsz)
    ex2.tick_params(labelsize=fsz)
    fx1.tick_params(labelsize=fsz)
    fx2.tick_params(labelsize=fsz)

    fig1.savefig("rhod.png",bbox_inches="tight")
    fig2.savefig("Cp.png",bbox_inches="tight")
    fig3.savefig("pexp.png",bbox_inches="tight")
    fig4.savefig("Rc.png",bbox_inches="tight")
    fig5.savefig("R0.png",bbox_inches="tight")
    fig6.savefig("res.png",bbox_inches="tight")
    fig7.savefig("cpe.png",bbox_inches="tight")
    plt.show()

