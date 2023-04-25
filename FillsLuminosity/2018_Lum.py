####################################################################
# @Giulia Faletti                                                  #
# Cutting and Fitting algorithm for the luminosity evolution model #
####################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import LoadData as ld
import LuminosityOptimization as lo
from lmfit import Model
import scipy.integrate as integrate
from matplotlib import cm

#Font setting
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica",
  "font.size": 12
})

# create an empty list to store the L_fit data for each fill number
L_fit_list = []
L_norm_list = []
Time_fit_list = []
Turn_fit_list = []
Fill_fit_list = []


#selecting current year
#year=16
# Selecting Current Year
#years = [16, 17, 18]
year = 18

#Plotting Luminosity Evolution Graphs
plot=True

#model parameters 
if year==16:
   n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2016()

elif year==17:
    n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2017()

elif year==18:
    n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2018()

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
def Cut_Fit(year, text):
    """Function that performs the necessary cut on the current fill

    Args:
        year (int): current year
        text (str): current fill

    Returns:
        L_fit: cutted data
        T_fit_real: times in second for the cutted data fit
        L_evol: raw luminosity data
        Times: raw Unix time
    """
    year=str(year)
    f=open('ATLAS/ATLAS_fill_20{}/{}_lumi_ATLAS.txt'.format(year, text),"r")
    lines=f.readlines()
    L_evolx=[]
    times=[]
    for x in lines:
        times.append(int(x.split(' ')[0]))  
        L_evolx.append(float(x.split(' ')[2]))
        
    f.close()
    Times = np.array(times)
    L_evol = np.array(L_evolx)

    #deleting the null values of the luminosity
    zero=np.where(L_evol<100)
    L_zero=np.delete(L_evol, zero)
    T_zero=np.delete(Times, zero)
        
    #check for enough points
    if len(L_zero)<10:
        zero=np.where(L_evol<5)
        L_zero=np.delete(L_evol, zero)
        T_zero=np.delete(Times, zero)

    #defining the derivative 
    dy = np.zeros(L_zero.shape)
    dy[0:-1] = np.diff(L_zero)/np.diff(T_zero)


    #start to slim down the fit interval       
    L_tofit=[]
    T_tofit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1.5:
            L_tofit.append(L_zero[idx])
            T_tofit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1.5:
            continue 

   
   #evaluating the differences between two subsequent points
    diff=np.diff(L_tofit)
        
    #deleting the discrepancies
    thr=np.max(abs(diff))*0.05
    idx_diff= np.where(abs(diff)>thr)[0]+1
        
    #new slim down of data
    L_tofit2=np.delete(L_tofit, idx_diff)
    T_tofit2=np.delete(T_tofit, idx_diff)
        
    #check for enough points
    if len(L_tofit2) < 30:
        L_tofit2=L_tofit
        T_tofit2=T_tofit
        
    L_fit=L_tofit2
    T_fit=T_tofit2     

    L_fit=np.array(L_fit)
    T_fit=np.array(T_fit) 

    #transforming the times from unix in seconds
    T_fit_real=T_fit-np.amin(T_fit)    
   
    return  L_fit, T_fit_real, L_evol, Times

def PlotLumiEvol(L_fit, T, text,year, L_ev, Tim):
    """Function that plots all the needed graphs.

    Args:
        L_fit: cutted data
        T_fit_real: times in second for the cutted data fit
        L_evol: raw luminosity data
        Times: raw Unix time
    """
    fig, ax= plt.subplots(num=1, clear=True)
    ax.plot(T/3600, L_fit, "c-", markersize=4)
    ax.set_title('{}'.format(text))
    ax.set_xlabel('Times [h]')
    ax.set_ylabel('Luminosity evolution [$\mathrm{Hz}/\mathrm{\mu b}$]')
    plt.legend(loc='best')
    plt.savefig('20{}_Graphs/{}_LuminosityEvolutionCutFit.pdf'.format(str(year),text))
    
    fig1, ax1= plt.subplots(num=1, clear=True)
    ax1.plot(Tim, L_ev, "k-")
    ax1.set_title('{}'.format(text))
    ax1.set_xlabel('Unix Times [s]')
    ax1.set_ylabel('Luminosity evolution [$\mathrm{Hz}/\mathrm{\mu b}$]')
    plt.savefig('20{}_Graphs/init/{}_LuminosityEvolution.pdf'.format(str(year),text))
    
    
    fig2, ax2= plt.subplots(num=1, clear=True)
    ax2.plot(T/3600, L_fit, "b.", markersize=4)
    #ax2.plot(T/3600, L, 'r-', label='Fit')
    #ax2.plot([], [], 'kx ', label=r'$\tilde{\chi}^2$='+'{:.2f}'.format(chi))
    ax2.set_yscale('log')
    ax2.set_title('{}'.format(text))
    ax2.set_xlabel('Times [h]')
    ax2.set_ylabel('Luminosity evolution [$\mathrm{Hz}/\mathrm{\mu b}$]')
    plt.legend(loc='best')
    plt.savefig('20{}_Graphs/log/{}_LuminosityEvolutionCutFit_log.pdf'.format(str(year),text))


    
# Define colormap
#cmap = plt.get_cmap('jet', len(FillNumber))
cmap = plt.get_cmap('jet', len(FillNumber))

fig, ax = plt.subplots()
for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
   
    #performing the cutting and fitting algorithm
    L_fit, T, L_ev, Tim= Cut_Fit(year, text)
    
    #defining the new time variable representing the number of turns
    Turn=[] 
    Turn=np.array(Turn)
    for el in T:
        tau=(f_rev*el+1)
        Turn=np.append(Turn, tau)

    L_min = min(L_fit)
    L_max = max(L_fit)
    #
    # Normlise L_fit
    L_Tau = L_fit/L_max 

    L_ma = np.array(L_max)
    L_norm = np.array(L_Tau)
        
    if plot==True:
       PlotLumiEvol(L_fit, T, text, year, L_ev, Tim)


    L_fit_list.append(L_fit) # append the L_fit data to the list
    L_norm_list.append(L_norm) # append the L_norm data to the list
    Time_fit_list.append(T)
    Turn_fit_list.append(Turn)
    Fill_fit_list.append(str(int(text)))

    fig, ax = plt.subplots()
        

    # loop over each L_fit data and plot it] ${L}/{L_i}$
    for i, L in enumerate(L_fit_list):

        print('i',i,'Year',year)
        #ax.plot(Time_fit_list[i]/3600, L*1e30/1e34, c=cmap(i))
        ax.plot(Time_fit_list[i]/3600, L*1e30/1e34, c=cmap(i))
         
    ax.set_xlabel(r'Times [$\mathrm{h}$]')
    #ax.set_ylabel(" Luminosity [$\mathrm{10}^{34}{cm}^{2}{s}^{-1}$]")
    ax.set_ylabel('Luminosity [$\mathrm{10}^{34}{cm}^{2}{s}^{-1}$]')
    # Normalizer
    norm = mpl.colors.Normalize(vmin=6638, vmax=7334) 
    # creating ScalarMappable
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
  
    plt.colorbar(sm, ticks=np.linspace(6600, 7400, 9))
    #ax.set_title('Luminosity evolution of fill {}'.format(text))
    fig.tight_layout()
    plt.savefig('{}_Time_Lumi.pdf'.format(str(year)))
                            

    plt.close("all")
    
    """fig1, ax1 = plt.subplots()
        

    # loop over each L_fit data and plot it]
    for i, L in enumerate(L_norm_list):

        ax1.plot(Turn_fit_list[i]/1e9, L_norm, c=cmap(i))
        # add a legend to the plot
    #ax1.legend()
    # set the axis labels and title
    ax1.set_xlabel(r'Time in Number of Turns 'r'$\tau$(10$^9$)')
    ax1.set_ylabel("Relative Luminosity [$\mathrm{cm}^{2}{s}^{-1}$]")
    # Normalizer
    norm = mpl.colors.Normalize(vmin=5000, vmax=5500) 
    # creating ScalarMappable
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
  
    plt.colorbar(sm, ticks=np.linspace(5000, 5500, 6))
    fig1.tight_layout()
    #ax.set_title('Luminosity evolution of fill {}'.format(text))
    plt.savefig('Turn_Norm_Lumi_{}.pdf'.format(str(year)))
    """

          