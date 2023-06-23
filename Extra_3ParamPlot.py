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

#Selecting Current Year
year=16


#plotting
plot=True

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
#print(FillNumber)

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

fill_nums = []
fill_nums1 = []
fill_nums2 = []
fill_nums3 = []

R = []
R1 = []
R2 = []
R3 = []

Ra = []
Ra1 = []
Ra2 = []
Ra3 = []

L = []
L1 = []
L2 = []
L3 = []

La = []
La1 = []
La2 = []
La3 = []



cmap = plt.get_cmap('jet', len(FillNumber))

fig, ax = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill

    delta_R = []
    delta_Lum = []
    
    with open('20{}/Extrapolation/EOF/{}_3param_Eof_R_extra.txt'.format(str(year), text), 'r') as f:


        lines=f.readlines()
        tm = []
        fn=[] #fill number
        L_Tau=[]
        R_adj=[]
        Lumi_evol=[]
        extra_R_adj=[]


        for x in lines:
            tm.append(float(x.split(' ')[0]))
            fn.append(float(x.split(' ')[1]))
            L_Tau.append(float(x.split(' ')[2]))
            R_adj.append(float(x.split(' ')[3]))
            Lumi_evol.append(float(x.split(' ')[4]))  
            extra_R_adj.append(float(x.split(' ')[5]))


    time_list = np.array(tm)  
    HR = time_list / 3600    
    #fn_list = np.array([text])
    fn_list = np.unique(fn)
    #fn_list = np.array(fn)
    L_Tau_list = np.array(L_Tau)
    R_adj_list = np.array(R_adj)
    Lumi_evol_list = np.array(Lumi_evol)
    extra_R_adj_list = np.array(extra_R_adj)

    

    delta_R_list = (extra_R_adj_list - R_adj_list).tolist()
    delta_Lumi_evol_list = (Lumi_evol_list - L_Tau_list).tolist()

    RelativeDiff_Lum = delta_Lumi_evol_list/L_Tau_list
    RelativeDiff_R = delta_Lumi_evol_list/R_adj_list


    

    #cax = ax.scatter(time_list, L_Tau_list, c=i, cmap='viridis', marker='.')
    ax.scatter(HR,RelativeDiff_Lum, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')


ax.set_xlabel('Time[h]')
ax.set_ylabel('${L/{L_i}}$')
ax.set_title('Relative-Difference on last-point for 20{}'.format(str(year)))
ax.grid(True)

#cbar.set_label('FillNumber')
norm = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm, ticks=np.linspace(5000, 5500, 6),label='FillNumber')


plt.savefig('20{}/Extrapolation/Relative_Diff_{}_3param_mod2_extrapolate_Lum.pdf'.format(str(year),str(year)))
plt.savefig('20{}/Extrapolation/Relative_Diff_{}_3param_mod2_extrapolate_Lum.png'.format(str(year),str(year)))

#plt.show()
plt.close("all")


cmap = plt.get_cmap('jet', len(FillNumber))

fig, ax = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill

    delta_R = []
    delta_Lum = []
    
    with open('20{}/Extrapolation/EOF/{}_3param_Eof_R_extra.txt'.format(str(year), text), 'r') as f:


        lines=f.readlines()
        tm = []
        fn=[] #fill number
        L_Tau=[]
        R_adj=[]
        Lumi_evol=[]
        extra_R_adj=[]


        for x in lines:
            tm.append(float(x.split(' ')[0]))
            fn.append(float(x.split(' ')[1]))
            L_Tau.append(float(x.split(' ')[2]))
            R_adj.append(float(x.split(' ')[3]))
            Lumi_evol.append(float(x.split(' ')[4]))  
            extra_R_adj.append(float(x.split(' ')[5]))


    time_list = np.array(tm)  
    HR = time_list / 3600    
    #fn_list = np.array([text])
    fn_list = np.unique(fn)
    #fn_list = np.array(fn)
    L_Tau_list = np.array(L_Tau)
    R_adj_list = np.array(R_adj)
    Lumi_evol_list = np.array(Lumi_evol)
    extra_R_adj_list = np.array(extra_R_adj)

    

    delta_R_list = (extra_R_adj_list - R_adj_list).tolist()
    delta_Lumi_evol_list = (Lumi_evol_list - L_Tau_list).tolist()

    RelativeDiff_Lum = delta_Lumi_evol_list/L_Tau_list
    RelativeDiff_R = delta_Lumi_evol_list/R_adj_list


    

    #cax = ax.scatter(time_list, L_Tau_list, c=i, cmap='viridis', marker='.')
    ax.scatter(HR,RelativeDiff_R, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')


ax.set_xlabel('Time[h]')
ax.set_ylabel(r'$\Delta {R_{adj}}^2 / {R_{adj}}^2$ $({\%})$')
ax.set_title('Relative-Difference on fit quality {R_{adj}}^2 for 20{}'.format(str(year)))
ax.grid(True)

#cbar.set_label('FillNumber')
norm = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm, ticks=np.linspace(5000, 5500, 6),label='FillNumber')


plt.savefig('20{}/Extrapolation/Relative_Diff_{}_3param_mod2_extrapolate_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/Extrapolation/Relative_Diff_{}_3param_mod2_extrapolate_R.png'.format(str(year),str(year)))

#plt.show()
plt.close("all")




fig, ax = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill

    delta_R = []
    delta_Lum = []
    
    with open('20{}/Extrapolation/EOF/{}_3param_Eof_R_extra.txt'.format(str(year), text), 'r') as f:


        lines=f.readlines()
        tm = []
        fn=[] 
        L_Tau=[]
        R_adj=[]
        Lumi_evol=[]
        extra_R_adj=[]


        for x in lines:
            tm.append(float(x.split(' ')[0]))
            fn.append(float(x.split(' ')[1]))
            L_Tau.append(float(x.split(' ')[2]))
            R_adj.append(float(x.split(' ')[3]))
            Lumi_evol.append(float(x.split(' ')[4]))  
            extra_R_adj.append(float(x.split(' ')[5]))


    time_list = np.array(tm)  
    HR = time_list / 3600    
    #fn_list = np.array([text])
    fn_list = np.unique(fn)
    #fn_list = np.array(fn)
    L_Tau_list = np.array(L_Tau)
    R_adj_list = np.array(R_adj)
    Lumi_evol_list = np.array(Lumi_evol)
    extra_R_adj_list = np.array(extra_R_adj)

    

    delta_R_list = (extra_R_adj_list - R_adj_list).tolist()
    delta_Lumi_evol_list = (Lumi_evol_list - L_Tau_list).tolist()

    RelativeDiff_Lum = delta_Lumi_evol_list/L_Tau_list
    RelativeDiff_R = delta_Lumi_evol_list/R_adj_list

    #cmap = plt.get_cmap('jet', len(HR))


    

    #cax = ax.scatter(time_list, L_Tau_list, c=i, cmap='viridis', marker='.')
    #ax.scatter(i,delta_Lumi_evol_list/L_Tau_list, color=cmap(HR), label='Hour {}'.format(HR),marker='.')
    #for j in range(len(HR)):

    for j in range(len(fn_list)):

        #ax.scatter(fn_list[j], delta_Lumi_evol_list[j] / L_Tau_list[j], c=[HR[j]], cmap='jet', marker='.')
        ax.scatter(fn_list[j], RelativeDiff_Lum[j], c=HR[j], cmap='jet', marker='.')


ax.set_xlabel('Time[h]')
ax.set_ylabel('${L/{L_i}}$')
ax.set_title('20{}'.format(str(year)))
ax.grid(True)

#cbar = fig.colorbar(cax)
#cbar.set_label('Time [h]')
cbar = fig.colorbar(ax.collections[0], ticks=np.arange(min(HR), max(HR) + 1), label='Time [h]')


plt.savefig('20{}/Extrapolation/R_{}_3param_mod2_extrapolate778.pdf'.format(str(year),str(year)))
plt.savefig('20{}/Extrapolation/R_{}_3param_mod2_extrapolate778.png'.format(str(year),str(year)))

#plt.show()
plt.close("all")









fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill

    delta_R = []
    delta_Lum = []
    
    with open('20{}/Extrapolation/EOF/{}_3param_Eof_R_extra.txt'.format(str(year), text), 'r') as f:


        lines=f.readlines()
        tm = []
        fn=[] 
        L_Tau=[]
        R_adj=[]
        Lumi_evol=[]
        extra_R_adj=[]


        for x in lines:
            tm.append(float(x.split(' ')[0]))
            fn.append(float(x.split(' ')[1]))
            L_Tau.append(float(x.split(' ')[2]))
            R_adj.append(float(x.split(' ')[3]))
            Lumi_evol.append(float(x.split(' ')[4]))  
            extra_R_adj.append(float(x.split(' ')[5]))


    time_list = np.array(tm)  
    HR = time_list / 3600    
    #fn_list = np.array([text])
    fn_list = np.unique(fn)
    #fn_list = np.array(fn)
    L_Tau_list = np.array(L_Tau)
    R_adj_list = np.array(R_adj)
    Lumi_evol_list = np.array(Lumi_evol)
    extra_R_adj_list = np.array(extra_R_adj)

       

    delta_R_list = (extra_R_adj_list - R_adj_list).tolist()
    delta_Lumi_evol_list = (Lumi_evol_list - L_Tau_list).tolist()

    RelativeDiff_Lum = delta_Lumi_evol_list/L_Tau_list
    RelativeDiff_R = delta_Lumi_evol_list/R_adj_list


    #cmap = plt.get_cmap('jet', len(HR))

    #for j in range(len(fn_list)):

        #ax.scatter(fn_list[j], delta_Lumi_evol_list[j] / L_Tau_list[j], c=HR[j], cmap='jet', marker='.')
    #ax.scatter(fn_list, delta_Lumi_evol_list/L_Tau_list,HR, c=HR, cmap='jet', marker='.')
    ax.scatter(fn_list,HR,RelativeDiff_Lum, marker='.')


ax.set_xlabel('FillNumber')
ax.set_zlabel('Relative Difference on {L/{L_i}}')
ax.set_ylabel('Time[h]')
ax.set_title('20{}'.format(str(year)))
ax.grid(True)

#cbar = fig.colorbar(cax)
#cbar.set_label('Time [h]')
#cbar = fig.colorbar(ax.collections[0], ticks=np.arange(min(HR), max(HR) + 1), label='Time [h]')


plt.savefig('20{}/Extrapolation/R_{}_3param_mod2_extrapolate779.pdf'.format(str(year),str(year)))
plt.savefig('20{}/Extrapolation/R_{}_3param_mod2_extrapolate779.png'.format(str(year),str(year)))

#plt.show()
plt.close("all")









