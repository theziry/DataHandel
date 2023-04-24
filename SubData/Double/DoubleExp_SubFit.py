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
from lmfit import Model
from scipy import integrate
import scipy.optimize
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
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

#defining the double-exponential decay fit function
def fit(x, a, b, c, d):
    return (a*np.exp((-b)*x))+(c*np.exp((-d)*x))
    #return (a*np.exp(-x/(b)))+(c*np.exp(-x/(d)))

model=Model(fit)

#defining the DA MODEL-4 fit function
def L_model(x, eps_ni, B, k):
        x = np.where(x>1, x, 2 )
        #x = np.where(x>1, x, 1.1) 
        np.nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None) 
        D = B / (np.log(x) ** k)
        #D = B *np.exp(-k*np.log((np.log(x))))
        L = (1/(1 + eps_ni * (x-1))**2) - ((1 + D**2)*np.exp(-(D**2))) * ((2 - (1 + D**2)*np.exp(-(D**2)/2)))
        return L

#####################################################  NOW WE DEFINE THE FUNCTION FOR DATA SELECTION  ###################################################################################
def Cut_Fit(year, text):
    """Function that performs the necessary cut on the current fill

    Args:
        year (int): current year
        text (str): current fill

    Returns:
        L_fit: cutted data
        T_fit_real: times in second for the cutted data fit
        Y: Luminosity evolution form the fit
        a: fitting parameter
        b: fitting parameter
        c: fitting parameter
        d: fitting parameter
        chi: reduced chi square of the fit
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

    ################## Normlise L_fit ############ L_norm = L_fit/L_fit[0] ##########################
    L_min = min(L_fit)
    L_max = max(L_fit)

    L_norm = L_fit/L_max 

    L_ma = np.array(L_max)
    L_Tau = np.array(L_norm)

    return T_fit_real,L_fit, L_ma, L_Tau, L_evol, Times

# create an empty list to store the L_fit data for each fill number
a_list = []
b_list = []
c_list = []
d_list = []
R_squared_list = []
Time_list = []
FillNbr_list = []

for i in range(len(FillNumber)):

    text = str(int(FillNumber[i])) #number of current fill
    #Times,L_fit,a,b,c,d=ld.CuttedData(year, text)
    T_fit_real,L_fit, L_ma, L_Tau, L_evol, Times = Cut_Fit(year, text)

    #### defining the Tau different time variable representing the number of turns #### 13-12-2022 ######
    Turn=[] 
    Turn=np.array(Turn)
    for el in T_fit_real:
       tau=(f_rev*el+1)
       Turn=np.append(Turn, tau) 

    Tau_norm_0= (Turn)/(np.max(Turn))
    Tau_norm_0= np.array(Tau_norm_0)


    ############################################ STORE COEFF ########################################################################################################## 
    with open('20{}/Coeff/{}_exp_FitCoeff.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 

    with open('20{}/data/{}_exp_Fitdata.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 

    with open('20{}/data/{}_Fitdata.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 
        

     
    p0 = [0.1,0.2e-9,0.1,0.1e-9]
    params, pcov = curve_fit(fit, Tau_norm_0,L_Tau,p0,bounds=(0, [1,500,1,500]),maxfev=5000)
    p0 = params 

    a,b,c,d = params

    print('p0', p0)

    LL = fit(Turn, a, b/(np.max(Turn)),c, d/(np.max(Turn)))

    residuals = L_Tau - LL
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((L_Tau - np.mean(L_Tau))**2)
    if ss_tot == 0:
        r_2 = 1
    else:
        r_2 = 1 - (ss_res / ss_tot) 
    # number of observations
    n = len(L_Tau)
    # number of predictors (not including the constant term)0.008*
    p = len(params) - 1
    # calculate the adjusted R-squared value
    Radj = 1 - (1 - r_2) * (n - 1) / (n - p - 1)

    for x in range(len(T_fit_real)):
        with open('20{}/data/{}_exp_Fitdata.txt'.format((str(year)),text), 'a') as f:  
            f.write(str(T_fit_real[x]))
            f.write(' ')
            f.write(str(L_fit[x]))
            f.write(' ')
            f.write(str(LL[x]))
            f.write('\n') 




    #################################################################################################################################


    start_time = T_fit_real[0]
    #start_time = 3600
    print('T_fit_real[0] start_time', T_fit_real[0])
    end_time = T_fit_real[-1]
    print('T_fit_real[-1] end_time', T_fit_real[-1])

    step_size = 3600  # Fixed step size
    
    subsets = []
    subsets2 = []

    #fig0, ax0= plt.subplots(num=1, clear=True)
        
    while start_time <= (end_time):
        # Select a subset of the data
        start_idx = int(start_time)
        end_idx = int(start_time + step_size)


        Time_subset = T_fit_real[T_fit_real <=end_idx]
        Time_subset[-1] = end_idx

        if end_idx >= T_fit_real[-1]:
            Time_subset[-1] = T_fit_real[-1]
        #L_Tau_subset = L_Tau[:len(Time_subset)]
        L_Tau_subset = L_Tau[T_fit_real <=end_idx]        

        Turn_subset=[] 
        Turn_subset=np.array(Turn_subset)
        for el in Time_subset:
           tau=(f_rev*el+1)
           Turn_subset=np.append(Turn_subset, tau)

        #print('length Turn', len(Turn))
        print('length Time', len(T_fit_real))
        print('length Time_subset', len(Time_subset))
        print('length Turn_subset', len(Turn_subset))
        
        tim = int(end_idx)
        tmp = int(Turn_subset[-1])

        # FIRST NORMLIZE THE NUMBER OF TURNS
        Tau_norm = (Turn_subset)/(np.max(Turn_subset))
        Tau_norm = np.array(Tau_norm)


        print('time time time ',(Time_subset[-1]))
        
        # define bounds as a percentage of p0
        lower_bound = 0.95 *np.array(p0)
        upper_bound = 1.05 *np.array(p0)
        bounds = (lower_bound, upper_bound)

        popt, pcov = curve_fit(fit, Tau_norm, L_Tau_subset,p0,bounds=([0.,0.,0.,0.], [1,500,1,500]),maxfev=5000) 
        #popt, pcov = curve_fit(fit, Tau_norm, L_Tau_subset,p0,bounds=bounds,maxfev=5000) 
             
        a,b,c,d = popt
        
        Y = fit(Turn_subset, a, b/(np.max(Turn_subset)),c, d/(np.max(Turn_subset)))

        residuals = L_Tau_subset - Y
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((L_Tau_subset - np.mean(L_Tau_subset))**2)

        if ss_tot == 0:
           r_2 = 1
        else:
           r_2 = 1 - (ss_res / ss_tot) 
        # number of observations
        n = len(L_Tau_subset)
        # number of predictors (not including the constant term)0.008*
        p = len(popt) - 1
        # calculate the adjusted R-squared value
        R_squared = 1 - (1 - r_2) * (n - 1) / (n - p - 1)

        print('FillNumber', text,'R_squared', R_squared,'ss_tot', ss_tot)
        print('popt p0 is',popt) 

        p0 = popt
        print('p0 is',p0)
        



        with open('20{}/Coeff/{}_exp_FitCoeff.txt'.format((str(year)),text), 'a') as f:
            f.write(str(tim))
            f.write(' ')
            f.write(str(tmp))
            f.write(' ')
            f.write(str(a))
            f.write(' ')
            f.write(str(b/(np.max(Turn_subset))))
            f.write(' ')
            f.write(str(c))
            f.write(' ')
            f.write(str(d/(np.max(Turn_subset))))
            f.write(' ')
            f.write(str(R_squared))
            f.write('\n')
        #start_idx += step
        start_time += step_size

        ################################################### plot ##################################################
        
        
        subsets.append((Time_subset, Y))

        fig0, ax0= plt.subplots(num=1, clear=True)
        ax0.plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')  
                  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, (Time_subset, Y) in enumerate(subsets):
            color = colors[i % len(colors)]
            label = f'Subset {i+1}'
            
            ax0.plot(Time_subset/3600, Y, color=color, label=label)
            
        ax0.set_title('FillNumber: {}'.format(text))
        ax0.set_xlabel(' Time in Number of Turns 'r'$\tau$(10$^9$)')
        ax0.set_ylabel('${L}/{L_i}$' ) 
        #plt.legend(loc='best')
        plt.savefig('20{}/exp/{}_Fit_double.pdf'.format(str(year),text))
        #plt.show()

        
        ################################################### 2 plot ##################################################

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

        fig, ax= plt.subplots(1,2, figsize=(10,6))        
        #fig, ax= plt.subplots(num=1, clear=True)
        #ax.plot(Timefit_list/3600, L_list/np.max(L_list), "b.", label='${L}/{L_i}$')
        ax[1].plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')
        #ax.plot(Timefit_list/3600, LL_list/np.max(LL_list), color = "coral", label='$Exp-Fit$')
        ax[1].plot(Time_list/3600, a_list/np.max(a_list), "y.",label='$a$')
        ax[1].plot(Time_list/3600, b_list/np.max(b_list), "r.",label='$b$')
        ax[1].plot(Time_list/3600, c_list/np.max(c_list), "c.",label='$c$')
        ax[1].plot(Time_list/3600, d_list/np.max(d_list), "m.",label='$d$')
        ax[1].set_xlim(-0.8, (np.max(Time_list/3600)+1))
        #ax.set_title('FillNumber: {}'.format(text))

        ax2 = ax[1].twinx()
        ax2.plot(Time_list/3600, R_squared_list/np.max(R_squared_list) ,"g.",label='${R_{adj}}$$^2$')
        #ax2.set_ylim(0.8, 1)
        ax[1].set_xlabel('Time [h]')
        ax[1].set_ylabel('$Normlised Values$')
        ax2.set_ylabel('$Normlised$ - ${R_{adj}}$$^2$')
        ax[1].set_title('FillNumber: {}'.format(text))
        #ax.xaxis.grid()
        for t in (Time_list/3600):
            plt.axvline((t), color='grey', linestyle='--',linewidth=0.5)
        lines, labels = ax[1].get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax[1].legend(lines + lines2, labels + labels2, loc='best')

        ax[0].plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')
        ax[0].set_xlim(-0.8, (np.max(Time_subset/3600)+1))
                  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, (Time_subset, Y) in enumerate(subsets):
            color = colors[i % len(colors)]
            label = f'Subset {i+1}'
            
            ax[0].plot(Time_subset/3600, Y, color=color, label=label)

            #for ti in (int(tim)/3600):
                #ax[0].axvline((ti), color='grey', linestyle='--',linewidth=0.5)
            
        ax[0].set_title('FillNumber: {}'.format(text))
        ax[0].set_xlabel(' Time in Number of Turns 'r'$\tau$(10$^9$)')
        ax[0].set_ylabel('${L}/{L_i}$' )
        plt.savefig('20{}/param/fitpar/{}_DotParam_double.pdf'.format(str(year),text))
        #plt.show()
        plt.close("all")

        #### A C ###########################

        fig, ax= plt.subplots(1,2, figsize=(10,6))        
        ax[1].plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')
        #ax.plot(Timefit_list/3600, LL_list/np.max(LL_list), color = "coral", label='$Exp-Fit$')
        ax[1].plot(Time_list/3600, a_list/np.max(a_list), "y.",label='$a$')
        #ax[1].plot(Time_list/3600, b_list/np.max(b_list), "r.",label='$b$')
        ax[1].plot(Time_list/3600, c_list/np.max(c_list), "c.",label='$c$')
        #ax[1].plot(Time_list/3600, d_list/np.max(d_list), "m.",label='$d$')
        ax[1].set_xlim(-0.8, (np.max(Time_list/3600)+1))
        #ax.set_title('FillNumber: {}'.format(text))

        ax2 = ax[1].twinx()
        ax2.plot(Time_list/3600, R_squared_list/np.max(R_squared_list) ,"g.",label='${R_{adj}}$$^2$')
        #ax2.set_ylim(0.8, 1)
        ax[1].set_xlabel('Time [h]')
        ax[1].set_ylabel('$Normlised Values$')
        ax2.set_ylabel('$Normlised$ - ${R_{adj}}$$^2$')
        ax[1].set_title('FillNumber: {}'.format(text))
        #ax.xaxis.grid()
        for t in (Time_list/3600):
            plt.axvline((t), color='grey', linestyle='--',linewidth=0.5)
        lines, labels = ax[1].get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax[1].legend(lines + lines2, labels + labels2, loc='best')

        ax[0].plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data') 

        ax[0].plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')
        ax[0].set_xlim(-0.8, (np.max(Time_subset/3600)+1))
                  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, (Time_subset, Y) in enumerate(subsets):
            color = colors[i % len(colors)]
            label = f'Subset {i+1}'
            
            ax[0].plot(Time_subset/3600, Y, color=color, label=label)

            #for ti in (int(tim)/3600):
                #ax[0].axvline((ti), color='grey', linestyle='--',linewidth=0.5)
            
        ax[0].set_title('FillNumber: {}'.format(text))
        ax[0].set_xlabel(' Time in Number of Turns 'r'$\tau$(10$^9$)')
        ax[0].set_ylabel('${L}/{L_i}$' )
        plt.savefig('20{}/param/ABCD/{}_A_C_DotParam_double.pdf'.format(str(year),text))
        #plt.show()
        plt.close("all")

        #### B and D #################

        fig, ax= plt.subplots(1,2, figsize=(10,6))        
        ax[1].plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')
        #ax.plot(Timefit_list/3600, LL_list/np.max(LL_list), color = "coral", label='$Exp-Fit$')
        #ax[1].plot(Time_list/3600, a_list/np.max(a_list), "y.",label='$a$')
        ax[1].plot(Time_list/3600, b_list/np.max(b_list), "r.",label='$b$')
        #ax[1].plot(Time_list/3600, c_list/np.max(c_list), "c.",label='$c$')
        ax[1].plot(Time_list/3600, d_list/np.max(d_list), "m.",label='$d$')
        ax[1].set_xlim(-0.8, (np.max(Time_list/3600)+1))
        #ax.set_title('FillNumber: {}'.format(text))

        ax2 = ax[1].twinx()
        ax2.plot(Time_list/3600, R_squared_list/np.max(R_squared_list) ,"g.",label='${R_{adj}}$$^2$')
        #ax2.set_ylim(0.8, 1)
        ax[1].set_xlabel('Time [h]')
        ax[1].set_ylabel('$Normlised Values$')
        ax2.set_ylabel('$Normlised$ - ${R_{adj}}$$^2$')
        ax[1].set_title('FillNumber: {}'.format(text))
        #ax.xaxis.grid()
        for t in (Time_list/3600):
            plt.axvline((t), color='grey', linestyle='--',linewidth=0.5)
        lines, labels = ax[1].get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax[1].legend(lines + lines2, labels + labels2, loc='best')

        ax[0].plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data') 

        ax[0].plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')
        ax[0].set_xlim(-0.8, (np.max(Time_subset/3600)+1))
                          
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, (Time_subset, Y) in enumerate(subsets):
            color = colors[i % len(colors)]
            label = f'Subset {i+1}'
                        
            ax[0].plot(Time_subset/3600, Y, color=color, label=label)

            #for ti in (Time_subset/3600):
                #ax[0].axvline((ti), color='grey', linestyle='--',linewidth=0.5)
            
        ax[0].set_title('FillNumber: {}'.format(text))
        ax[0].set_xlabel(' Time in Number of Turns 'r'$\tau$(10$^9$)')
        ax[0].set_ylabel('${L}/{L_i}$' )
        plt.savefig('20{}/param/ABCD/{}_B_D_DotParam_double.pdf'.format(str(year),text))
        #plt.show()
        plt.close("all")


        
       

   
#defining the stop time of the program      
stop=t.time()
#evaluating the run time 
runtime_seq=stop-start
print('The runtime is:', runtime_seq, '[s]')      
