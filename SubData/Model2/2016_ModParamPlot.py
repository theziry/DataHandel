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
        
FillNbr_list = []
  
# Define colormap
#cmap = plt.get_cmap('jet', len(FillNumber))
cmap = plt.get_cmap('jet', len(FillNumber))

fig, ax = plt.subplots()
for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    # Plot B_list vs Time_list

    ax.plot(Time_list/3600, R_squared_list*100, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
ax.set_ylim([-0.5, 101])
ax.set_yticks(np.linspace(0, 100, 11))
ax.set_xlabel('Time [h]')
ax.set_ylabel('${R_{adj}}$$^2${({\%})}')
ax.set_title('$Model$-$2$ ${R_{adj}}$$^2$ 20%s' % str(year))
# Normalizer
norm = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig.tight_layout()
plt.savefig('20{}/plot/R_{}_mod2.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")


fig1, ax1 = plt.subplots()
for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax1.plot(Time_list/3600, Eps_list*1e10, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

ax1.set_ylim([-0.5, 14])
ax1.set_yticks(np.linspace(0, 13, 14))
ax1.set_xlabel('Time [h]')
ax1.set_ylabel('$\epsilon$$*$$N_i$ {(10$^{-10}$)}')
ax1.set_title('$Model$-$2$ $\epsilon$$*$$N_i$ Parameter - 20%s' % str(year))
# Normalizer
norm1 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm1 = plt.cm.ScalarMappable(cmap=cmap, norm=norm1)
sm1.set_array([])
plt.colorbar(sm1, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig1.tight_layout()
plt.savefig('20{}/plot/EPS_{}_mod2.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")



fig2, ax2 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax2.plot(Time_list/3600, B_list, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')


ax2.set_ylim([-0.5, 301])
ax2.set_yticks(np.linspace(0, 300, 9))
#ax2.set_yticks(np.linspace(0, 140, 8))
ax2.set_xlabel('Time [h]')
ax2.set_ylabel(r'$\rho_\ast$')
ax2.set_title(r'$Model$-$2$ $\rho_\ast$ - 20{}'.format(str(year)))
# Normalizer
norm2 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm2 = plt.cm.ScalarMappable(cmap=cmap, norm=norm2)
sm2.set_array([])
plt.colorbar(sm2, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig2.tight_layout()
plt.savefig('20{}/plot/B_{}_mod2.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")


fig3, ax3 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax3.plot(Time_list/3600, k_list, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
ax3.set_ylim([0.45, 1.45])
ax3.set_yticks(np.linspace(0.5, 1.4, 10))
ax3.set_xlabel('Time [h]')
ax3.set_ylabel('$k$')
ax3.set_title('$Model$-$2$ $k$ Parameter - 20{}'.format(str(year)))
# Normalizer
norm3 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm3 = plt.cm.ScalarMappable(cmap=cmap, norm=norm3)
sm3.set_array([])
plt.colorbar(sm3, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig3.tight_layout()
plt.savefig('20{}/plot/K_{}_mod2.pdf'.format(str(year),str(year)))

plt.show()
plt.close("all")

#######################################################################################################
fig4, ax4 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax4.plot(Time_list/3600, B_list, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')


#ax4.set_ylim([-0.5, 301])
# ax4.set_yticks(np.linspace(0, 300, 9))
#ax2.set_yticks(np.linspace(0, 140, 8))
ax4.set_xlabel('Time [h]')
ax4.set_ylabel(r'$\rho_\ast$')
ax4.set_title(r'$Model$-$2$ $\rho_\ast$ - 20{}'.format(str(year)))
# Normalizer
norm4 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm4 = plt.cm.ScalarMappable(cmap=cmap, norm=norm4)
sm4.set_array([])
plt.colorbar(sm4, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig4.tight_layout()
plt.savefig('20{}/plot/B_{}.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")

##############################################################

fig5, ax5 = plt.subplots()
for i in range(len(FillNumber)):
    text = str(int(FillNumber[i]))

    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    # Plot B_list vs Time_list

    ax5.plot(Time_list/3600, R_squared_list*100, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
#ax5.set_ylim([-0.5, 101])
#ax5.set_yticks(np.linspace(0, 100, 11))
ax5.set_xlabel('Time [h]')
ax5.set_ylabel('${R_{adj}}$$^2${({\%})}')
ax5.set_title('$Model$-$2$ ${R_{adj}}$$^2$ 20%s' % str(year))
# Normalizer
norm5 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm5 = plt.cm.ScalarMappable(cmap=cmap, norm=norm5)
sm5.set_array([])
plt.colorbar(sm5, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig5.tight_layout()
plt.savefig('20{}/plot/R_{}.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")


############################################################################################################################
################################### Radj corr with parameters check #########################################################

fig7, ax7 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax7.plot(k_list, R_squared_list*100, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
#ax7.set_ylim([0, 1.5])
#ax7.set_yticks(np.linspace(0.1, 1.4, 14))  ('${R_{adj}}$$^2${({\%})}')
#ax.set_title('$Model$-$2$ ${R_{adj}}$$^2$ 20%s' % str(year))

ax7.set_xlabel('$k$')
ax7.set_ylabel('${R_{adj}}$$^2${({\%})}')
ax7.set_title('$Model$-$2$ $k$ and ${R_{adj}}$$^2$ 20%s' % str(year).format(str(year)))
# Normalizer
norm7 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm7 = plt.cm.ScalarMappable(cmap=cmap, norm=norm7)
sm7.set_array([])
plt.colorbar(sm7, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig7.tight_layout()
plt.savefig('20{}/plot/K_Radj_{}.pdf'.format(str(year),str(year)))

plt.show()
plt.close("all")


##############################################################


fig8, ax8 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax8.plot(B_list, R_squared_list*100, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
#ax3.set_ylim([0, 1.5])
#ax3.set_yticks(np.linspace(0.1, 1.4, 14)) ax7.set_title('$Model$-$2$ $k$ and ${R_{adj}}$$^2$ 20%s'

ax8.set_xlabel(r'$\rho_\ast$')
ax8.set_ylabel('${R_{adj}}$$^2${({\%})}')
ax8.set_title(r'$Model$-$2$ ${R_{adj}}$$^2$ and $\rho_\ast$ -  20%s' % str(year).format(str(year)))
# Normalizer
norm8 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm8 = plt.cm.ScalarMappable(cmap=cmap, norm=norm8)
sm8.set_array([])
plt.colorbar(sm8, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig8.tight_layout()
plt.savefig('20{}/plot/Rho_Radj_{}.pdf'.format(str(year),str(year)))

plt.show()
plt.close("all")


##############################################################

fig9, ax9 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax9.plot(Eps_list*1e10, R_squared_list*100, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
#ax3.set_ylim([0, 1.5])
#ax3.set_yticks(np.linspace(0.1, 1.4, 14)) ax7.set_title('$Model$-$2$ $k$ and ${R_{adj}}$$^2$ 20%s'

ax9.set_xlabel('$\epsilon$$*$$N_i$ {(10$^{-10}$)}')
ax9.set_ylabel('${R_{adj}}$$^2${({\%})}')
ax9.set_title(r'$Model$-$2$ ${R_{adj}}$$^2$ and $\epsilon$$*$$N_i$ -  20%s' % str(year).format(str(year)))
# Normalizer
norm9 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm9 = plt.cm.ScalarMappable(cmap=cmap, norm=norm9)
sm9.set_array([])
plt.colorbar(sm9, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig9.tight_layout()
plt.savefig('20{}/plot/Eps_Radj_{}.pdf'.format(str(year),str(year)))

plt.show()
plt.close("all")


##############################################################


fig10, ax10 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax10.plot(B_list, k_list, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
#ax3.set_ylim([0, 1.5])
#ax3.set_yticks(np.linspace(0.1, 1.4, 14))

ax10.set_xlabel(r'$\rho_\ast$')
ax10.set_ylabel('$k$')
ax10.set_title(r'$Model$-$2$ $k$ and $\rho_\ast$ - 20{}'.format(str(year)))
# Normalizer
norm10 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm10 = plt.cm.ScalarMappable(cmap=cmap, norm=norm10)
sm10.set_array([])
plt.colorbar(sm10, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig10.tight_layout()
plt.savefig('20{}/plot/K_Rho_{}.pdf'.format(str(year),str(year)))

plt.show()
plt.close("all")


##############################################################


fig11, ax11 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax11.plot(B_list, Eps_list*1e10, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
#ax3.set_ylim([0, 1.5])
#ax3.set_yticks(np.linspace(0.1, 1.4, 14))

ax11.set_xlabel(r'$\rho_\ast$')
ax11.set_ylabel('$\epsilon$$*$$N_i$ {(10$^{-10}$)}')
ax11.set_title(r'$Model$-$2$ $\epsilon$$*$$N_i$ and $\rho_\ast$ - 20{}'.format(str(year)))
# Normalizer
norm11 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm11 = plt.cm.ScalarMappable(cmap=cmap, norm=norm11)
sm11.set_array([])
plt.colorbar(sm11, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig11.tight_layout()
plt.savefig('20{}/plot/Rho_Eps_{}.pdf'.format(str(year),str(year)))

plt.show()
plt.close("all")


##############################################################


fig12, ax12 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    ############################################ READ COEFF ##########################################################################################################
    with open('20{}/Coeff/{}_mod2_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[2]))
            B.append(float(x.split(' ')[3]))
            k.append(float(x.split(' ')[4]))  
            R_squared.append(float(x.split(' ')[5]))

          
    Time_list = np.array(time)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_squared_list = np.array(R_squared)

    ax12.plot(Eps_list*1e10, k_list, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
#ax3.set_ylim([0, 1.5])
#ax3.set_yticks(np.linspace(0.1, 1.4, 14))

ax12.set_xlabel('$\epsilon$$*$$N_i$ {(10$^{-10}$)}')
ax12.set_ylabel('$k$')
ax12.set_title('$Model$-$2$ $k$ and $\epsilon$$*$$N_i$ - 20{}'.format(str(year)))
# Normalizer
norm12 = mpl.colors.Normalize(vmin=5017, vmax=5451) 
# creating ScalarMappable
sm12 = plt.cm.ScalarMappable(cmap=cmap, norm=norm12)
sm12.set_array([])
plt.colorbar(sm12, ticks=np.linspace(5000, 5500, 6),label='FillNumber')
fig12.tight_layout()
plt.savefig('20{}/plot/K_EPS_{}.pdf'.format(str(year),str(year)))

plt.show()
plt.close("all")


##############################################################