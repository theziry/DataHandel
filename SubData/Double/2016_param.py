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
cmap = plt.get_cmap('jet', len(FillNumber))

fig, ax = plt.subplots()


for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill

    with open('20{}/Coeff/{}_exp_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        a=[]
        b=[]
        c=[]
        d=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            a.append(float(x.split(' ')[2]))
            b.append(float(x.split(' ')[3]))
            c.append(float(x.split(' ')[4]))  
            d.append(float(x.split(' ')[5]))  
            R_squared.append(float(x.split(' ')[6]))

          
    Time_list = np.array(time)
    a_list = np.array(a)
    b_list = np.array(b)
    c_list = np.array(c)
    d_list = np.array(d)
    R_squared_list = np.array(R_squared)

    ax.plot(Time_list/3600, R_squared_list*100, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')
    
# Set plot parameters
ax.set_ylim([-0.5, 101.5])
ax.set_yticks(np.linspace(0, 100, 11))
ax.set_xlabel('Time [h]')
ax.set_ylabel('${R_{adj}}$$^2${({\%})}')
ax.set_title('$Double$ - $Exponential$ ${R_{adj}}$$^2$ 20%s' % str(year))
# Normalizer
norm = mpl.colors.Normalize(vmin=5017, vmax=5451) #2016
# creating ScalarMappable
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm, ticks=np.linspace(5000, 5500, 6),label='FillNumber') #2016
fig.tight_layout()
plt.savefig('20{}/plot/R_{}_Exp.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")

fig1, ax1 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    with open('20{}/Coeff/{}_exp_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        a=[]
        b=[]
        c=[]
        d=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            a.append(float(x.split(' ')[2]))
            b.append(float(x.split(' ')[3]))
            c.append(float(x.split(' ')[4]))  
            d.append(float(x.split(' ')[5]))  
            R_squared.append(float(x.split(' ')[6]))

          
    Time_list = np.array(time)
    a_list = np.array(a)
    b_list = np.array(b)
    c_list = np.array(c)
    d_list = np.array(d)
    R_squared_list = np.array(R_squared)

    ax1.plot(Time_list/3600, a_list, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
ax1.set_xlabel('Time [h]')
ax1.set_ylabel('$a$')
ax1.set_title('$a$- parameter 20{}'.format(str(year)))
# Normalizer
norm1 = mpl.colors.Normalize(vmin=5017, vmax=5451) #2016
# creating ScalarMappable
sm1 = plt.cm.ScalarMappable(cmap=cmap, norm=norm1)
sm1.set_array([])
plt.colorbar(sm1, ticks=np.linspace(5000, 5500, 6),label='FillNumber') #2016
fig1.tight_layout()
plt.savefig('20{}/plot/a_{}_Exp.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")


fig2, ax2 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    with open('20{}/Coeff/{}_exp_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        a=[]
        b=[]
        c=[]
        d=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            a.append(float(x.split(' ')[2]))
            b.append(float(x.split(' ')[3]))
            c.append(float(x.split(' ')[4]))  
            d.append(float(x.split(' ')[5]))  
            R_squared.append(float(x.split(' ')[6]))

          
    Time_list = np.array(time)
    a_list = np.array(a)
    b_list = np.array(b)
    c_list = np.array(c)
    d_list = np.array(d)
    R_squared_list = np.array(R_squared)

    ax2.plot(Time_list/3600, b_list, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
#ax2.set_ylim([-0.5, 30.5])
#ax2.set_yticks(np.linspace(0, 30, 10))
ax2.set_xlabel('Time [h]')
ax2.set_ylabel('$b$')
ax2.set_title('$b$- parameter 20{}'.format(str(year)))
# Normalizer
norm2 = mpl.colors.Normalize(vmin=5017, vmax=5451) #2016
# creating ScalarMappable
sm2 = plt.cm.ScalarMappable(cmap=cmap, norm=norm2)
sm2.set_array([])
plt.colorbar(sm2, ticks=np.linspace(5000, 5500, 6),label='FillNumber') #2016
fig2.tight_layout()
plt.savefig('20{}/plot/b_{}_Exp.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")

####################################################################################################################################

fig3, ax3 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    with open('20{}/Coeff/{}_exp_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        a=[]
        b=[]
        c=[]
        d=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            a.append(float(x.split(' ')[2]))
            b.append(float(x.split(' ')[3]))
            c.append(float(x.split(' ')[4]))  
            d.append(float(x.split(' ')[5]))  
            R_squared.append(float(x.split(' ')[6]))

          
    Time_list = np.array(time)
    a_list = np.array(a)
    b_list = np.array(b)
    c_list = np.array(c)
    d_list = np.array(d)
    R_squared_list = np.array(R_squared)

    ax3.plot(Time_list/3600, c_list, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')
    
# Set plot parameters
ax3.set_xlabel('Time [h]')
ax3.set_ylabel('$c$')
ax3.set_title('$c$- parameter 20{}'.format(str(year)))
# Normalizer
norm3 = mpl.colors.Normalize(vmin=5017, vmax=5451) #2016
# creating ScalarMappable
sm3 = plt.cm.ScalarMappable(cmap=cmap, norm=norm3)
sm3.set_array([])
plt.colorbar(sm3, ticks=np.linspace(5000, 5500, 6),label='FillNumber') #2016
fig3.tight_layout()
plt.savefig('20{}/plot/c_{}_Exp.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")


#############################################################################################################################################


fig4, ax4 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    with open('20{}/Coeff/{}_exp_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        a=[]
        b=[]
        c=[]
        d=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            a.append(float(x.split(' ')[2]))
            b.append(float(x.split(' ')[3]))
            c.append(float(x.split(' ')[4]))  
            d.append(float(x.split(' ')[5]))  
            R_squared.append(float(x.split(' ')[6]))

          
    Time_list = np.array(time)
    a_list = np.array(a)
    b_list = np.array(b)
    c_list = np.array(c)
    d_list = np.array(d)
    R_squared_list = np.array(R_squared)

    ax4.plot(Time_list/3600, d_list, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
ax4.set_xlabel('Time [h]')
ax4.set_ylabel('$d$')
ax4.set_title('$d$- parameter 20{}'.format(str(year)))
# Normalizer
norm4 = mpl.colors.Normalize(vmin=5017, vmax=5451) #2016
# creating ScalarMappable
sm4 = plt.cm.ScalarMappable(cmap=cmap, norm=norm4)
sm4.set_array([])
plt.colorbar(sm4, ticks=np.linspace(5000, 5500, 6),label='FillNumber') #2016
fig4.tight_layout()
plt.savefig('20{}/plot/d_{}_Exp.pdf'.format(str(year),str(year)))

plt.show()
plt.close("all")


fig5, ax5 = plt.subplots()

for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    with open('20{}/Coeff/{}_exp_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        time=[] #fill number
        a=[]
        b=[]
        c=[]
        d=[]
        R_squared=[]


        for x in lines:
            time.append(float(x.split(' ')[0]))
            a.append(float(x.split(' ')[2]))
            b.append(float(x.split(' ')[3]))
            c.append(float(x.split(' ')[4]))  
            d.append(float(x.split(' ')[5]))  
            R_squared.append(float(x.split(' ')[6]))

          
    Time_list = np.array(time)
    a_list = np.array(a)
    b_list = np.array(b)
    c_list = np.array(c)
    d_list = np.array(d)
    R_squared_list = np.array(R_squared)

    ax5.plot(Time_list/3600, R_squared_list*100, color=cmap(i), label='FillNumber {}'.format(FillNumber[i]),marker='.')

# Set plot parameters
ax5.set_xlabel('Time [h]')
ax5.set_ylabel('${R_{adj}}$$^2${({\%})}')
ax5.set_title('$Double$ - $Exponential$ ${R_{adj}}$$^2$ 20%s' % str(year))
# Normalizer
norm5 = mpl.colors.Normalize(vmin=5017, vmax=5451) #2016
# creating ScalarMappable
sm5 = plt.cm.ScalarMappable(cmap=cmap, norm=norm5)
sm5.set_array([])
plt.colorbar(sm5, ticks=np.linspace(5000, 5500, 6),label='FillNumber') #2016
fig5.tight_layout()
plt.savefig('20{}/plot/R_{}_Exp_all.pdf'.format(str(year),str(year)))
plt.show()
plt.close("all")









