#################################################################################################################################################################
##
######## Thiziri Amezza 23-02-2023 #############################################################################################################################
##
#################################################################################################################################################################
import numpy as np
import LoadData as ld
import pandas as pd
import LuminosityOptimization as lo
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import matplotlib.ticker as ticker
import time as t
#from lmfit import Model
from scipy import integrate
import scipy.optimize
from scipy.optimize import curve_fit
#_____________________________________________
plt.close("all")

#Font setting
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica",
  "font.size": 12
})

#defining the start time of the program
start=t.time()
# create an empty list to store the L_fit data for each fill number

year = 16

#Plotting Luminosity Evolution Graphs
plot=True

#model parameters 

#loading fill number
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()

if year==16:
        FillNumber=FillNumber16

elif year==17:
    FillNumber=FillNumber17
    FillNumber_Prev=FillNumber16
    previous_year=16
elif year==18:
    FillNumber=FillNumber18
    FillNumber_Prev=FillNumber17
    previous_year=17

#load turnaround times and fill times 
data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 = ld.loadFill()
data_ta16_sec = data_ta16*3600 
data_tf16_sec = data_tf16*3600  
data_ta17_sec = data_ta17*3600 
data_tf17_sec = data_tf17*3600
data_ta18_sec = data_ta18*3600 
data_tf18_sec = data_tf18*3600

# Fill to delete:
skip16=[5111,5112,5116,5117,5257,5258,5264,5406,5427,5439,5451]
skip17=[5837,5840,5856,6019,6055,6060,6082,6084,6089,6116,6152,6156,6158,6160,6167,6168,6169,6170,6193,6258,6268,6271,6272,6275,6283,6285,6287,6288,6291,6298,6303,6304,6305,6308,6311,6312,6314,6317,6324,6325,6337,6346,6356,6360,6362,6364,6371]
skip18=[6640,6654,6659,6688,6690,6693,6694,6696,6700,6706,6729,6741,6747,6752,6782,6778,6882,6890,6892,6901,6929,6939,6940,6961,7026,7033,7036,7039,7056,7069,7105,7118,7120,7135,7145,7217,7245,7259,7271,7324]
skip = skip16 + skip17 + skip18
for fill in skip:
    FillNumber=np.delete(FillNumber, (np.where(FillNumber==(fill))[0]))
for i in range(len(FillNumber)):
    if FillNumber[i] in skip:
        continue


for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
   
    with open('20{}/Coeff/{}_exp_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        A=[]
        B=[]
        C=[]
        D=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            A.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            C.append(float(x.split(' ')[4])) 
            D.append(float(x.split(' ')[5]))  
            R_squared.append(float(x.split(' ')[6]))

          
    Time_list = np.array(time)
    a_list = np.array(A)
    b_list = np.array(B)
    c_list = np.array(C)
    d_list = np.array(D)
    R_squared_list = np.array(R_squared)
    
    with open('20{}/data/{}_exp_Fitdata.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        Lfit=[]
        LL=[]

        

        for x in lines:
            time.append(float(x.split(' ')[0]))
            Lfit.append(float(x.split(' ')[1]))
            LL.append(float(x.split(' ')[2]))

            
          
    Timefit_list = np.array(time)
    L_list = np.array(Lfit)
    LL_list = np.array(LL)



    
    fig, ax= plt.subplots(num=1, clear=True)
    ax.plot(Timefit_list/3600, L_list/np.max(L_list), "b.", label='${L}/{L_i}$')
    #ax.plot(Timefit_list/3600, LL_list/np.max(LL_list), color = "coral", label='$Exp-Fit$')
    ax.plot(Time_list/3600, a_list/np.max(a_list), "y.",label='$a$')
    ax.plot(Time_list/3600, b_list/np.max(b_list), "r.",label='$b$')
    ax.plot(Time_list/3600, c_list/np.max(c_list), "c.",label='$c$')
    ax.plot(Time_list/3600, d_list/np.max(d_list), "m.",label='$d$')
    ax.set_xlim(-1, (np.max(Time_list/3600)+1))
    #ax.set_title('FillNumber: {}'.format(text))

    ax2 = ax.twinx()
    ax2.plot(Time_list/3600, R_squared_list/np.max(R_squared_list) ,"g.",label='${R_{adj}}$$^2$')
    #ax2.set_ylim(0.8, 1)
    ax.set_xlabel('Time [h]')
    ax.set_ylabel('$Normlised Values$')
    ax2.set_ylabel('$Normlised$-${R_{adj}}$$^2$')
    ax.set_title('FillNumber: {}'.format(text))
    #ax.xaxis.grid()
    for t in (Time_list/3600):
       plt.axvline((t), color='grey', linestyle='--',linewidth=0.5)

    #ax.grid(which='both', axis='x', linestyle=':', linewidth=0.5)  # Show grid lines on both major and minor ticks of x-axis
    #ax2.grid(False)    
    #ax.legend(loc='best')
    #ax2.legend(loc='best')
    # Create a combined legend for ax and ax2
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc='best')
    # Remove the individual legends from ax and ax2
    #ax.get_legend().remove()
    #ax2.get_legend().remove()
    plt.savefig('20{}/param/{}_DotParam_double.pdf'.format(str(year),text))
    #plt.show()
    plt.close("all")
    
    fig, ax= plt.subplots(num=1, clear=True)
    ax.plot(Timefit_list/3600, L_list/np.max(L_list), "b.", label='${L}/{L_i}$')
    #ax.plot(Timefit_list/3600, LL_list/np.max(LL_list), color = "coral", label='$Exp-Fit$')
    ax.plot(Time_list/3600, a_list/np.max(a_list), "y-",label='$a$')
    ax.plot(Time_list/3600, b_list/np.max(b_list), "r-",label='$b$')
    ax.plot(Time_list/3600, c_list/np.max(c_list), "c-",label='$c$')
    ax.plot(Time_list/3600, d_list/np.max(d_list), "m-",label='$d$')
    ax.set_xlim(-1, (np.max(Time_list/3600)+1))
    ax.set_xlabel('Time [h]')
    ax.set_ylabel('$Normlised Values$')
    ax.set_title('FillNumber: {}'.format(text))
    ax2 = ax.twinx()
    ax2.plot(Time_list/3600, R_squared_list/np.max(R_squared_list) ,"g-",label='${R_{adj}}$$^2$')
    #ax2.set_ylim(0.8, 1)
    ax2.set_ylabel('$Normlised$-${R_{adj}}$$^2$')

    #ax.xaxis.grid()
    for t in (Time_list/3600):
       plt.axvline((t), color='grey', linestyle='--',linewidth=0.5)
    # Create a combined legend for ax and ax2
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc='best')
    plt.savefig('20{}/param/{}_LineParam_double.pdf'.format(str(year),text))
    #plt.show()
    plt.close("all")

    
    fig2, ax2= plt.subplots(3,2,figsize=(10,8))
    #plt.subplots(5)

    ax2[0,0].plot(Timefit_list/3600, L_list/np.max(L_list), "b.")
    #ax2[0,0].plot(Timefit_list/3600, LL_list/np.max(LL_list), color = "coral", label='$Exp-Fit$')
    ax2[0,0].set_xlim(-1, (np.max(Timefit_list/3600)+1))
    ax2[0,0].set_ylabel('${L}/{L_i}$')
    ax2[0,0].set_title('FillNumber: {}'.format(text))
   
    ax2[1,0].plot(Time_list/3600, b_list/np.max(b_list), "r.")
    ax2[1,0].set_xlim(-1, (np.max(Time_list/3600)+1))
    ax2[1,0].set_ylabel('$b$')


    ax2[2,0].plot(Time_list/3600, d_list/np.max(d_list), "m.") 
    ax2[2,0].set_xlim(-1, (np.max(Time_list/3600)+1))
    ax2[2,0].set_ylabel('$d$')
    
    
    ax2[0,1].plot(Time_list/3600, R_squared_list/np.max(R_squared_list) ,"g.")
    #ax2[0,1].plot(Timefit_list/3600, LL_list/np.max(LL_list), color = "coral", label='$Exp-Fit$')
    ax2[0,1].set_xlim(-1, (np.max(Timefit_list/3600)+1))
    ax2[0,1].set_ylabel('$Normlised$-${R_{adj}}$$^2$')
    ax2[0,1].set_title('FillNumber: {}'.format(text))

    ax2[1,1].plot(Time_list/3600, a_list/np.max(a_list) ,"y.")
    ax2[1,1].set_xlim(-1, (np.max(Time_list/3600)+1))
    #ax2[1,1].set_xlabel('Time [h]')
    ax2[1,1].set_ylabel('$a$')
    
    ax2[2,1].plot(Time_list/3600, c_list/np.max(c_list), "c.")
    ax2[2,1].set_xlim(-1, (np.max(Time_list/3600)+1))
    ax2[2,1].set_xlabel('Time [h]')
    ax2[2,1].set_ylabel('$c$')

    
    plt.savefig('20{}/param/{}_Param_double.pdf'.format(str(year),text))
    #plt.show()
    plt.close("all")



