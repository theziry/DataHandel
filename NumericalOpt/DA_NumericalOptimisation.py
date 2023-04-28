import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import quad
import LoadData as ld
import LuminosityOptimization as lo
#from autograd import jacobian

plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica",
  "font.size": 12
})

#selecting current year
year=17

#model parameters 
if year==16:
   n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2016()

elif year==17:
    n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2017()

elif year==18:
    n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2018()



#loading fill number
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()

#load turnaround times and fill times - cheking
data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 = ld.loadFill()
data_ta16_sec = data_ta16*3600 
data_tf16_sec = data_tf16*3600  
data_ta17_sec = data_ta17*3600 
data_tf17_sec = data_tf17*3600
data_ta18_sec = data_ta18*3600 
data_tf18_sec = data_tf18*3600


if year==16:
    FillNumber=FillNumber16
    ta=data_ta16_sec
    tf=data_tf16_sec
elif year==17:
    FillNumber=FillNumber17
    FillNumber_Prev=FillNumber16
    previous_year=16
    ta=data_ta17_sec
    tf=data_tf17_sec
elif year==18:
    FillNumber=FillNumber18
    FillNumber_Prev=FillNumber17
    previous_year=17
    ta=data_ta18_sec
    tf=data_tf18_sec
    
Lmes16,Lmes17,Lmes18=ld.MeasuredLuminosity()
    
if year==16:
   Lmes=Lmes16
elif year==17:
   Lmes=Lmes17
elif year==18:
   Lmes=Lmes18

#Fill to delete [6638,6666,6174,7061,7065,7087,7124,7127]
skip16=[5111,5112,5116,5117,5257,5258,5264,5406,5427,5439,5451]
skip17=[5837,5840,5856,6019,6055,6060,6082,6084,6089,6116,6152,6156,6158,6160,6167,6168,6169,6170,6193,6258,6268,6271,6272,6275,6283,6285,6287,6288,6291,6298,6303,6304,6305,6308,6311,6312,6314,6317,6324,6325,6337,6346,6356,6360,6362,6364,6371]
skip18=[6640,6654,6659,6688,6690,6693,6694,6696,6700,6706,6729,6741,6747,6752,6782,6778,6882,6890,6892,6901,6929,6939,6940,6961,7026,7033,7036,7039,7056,7069,7105,7118,7120,7135,7145,7217,7245,7259,7271,7324]
skip = skip16 + skip17 + skip18

for fill in skip:
    FillNumber=np.delete(FillNumber, (np.where(FillNumber==(fill))[0]))
for i in range(len(FillNumber)):
    if FillNumber[i] in skip:
        continue

f=open('MODEL2/20{}/Coeff/mod2_Coeff_20{}.txt'.format(str(year),str(year)),"r")
lines=f.readlines()
a=[] #eps
b=[] #B
c=[] #k
for x in lines:
    a.append(float(x.split(' ')[1]))
    b.append(float(x.split(' ')[2]))
    c.append(float(x.split(' ')[3])) 

f.close()


#model parameters and initial guesses
a=np.array(a)
b=np.array(b)
c=np.array(c)

# define the new objective function
def fun(t1): 
    result = np.empty(len(x0))  
    for i in range(len(x0)): 
        #D = b[i] / (np.log(t1[i]) ** c[i])
        D = b[i] * np.power(c[i] / (2 * np.exp(1)), c[i]) / (np.power(np.log(t1[i]), c[i]))
        lam = lambda x1: (1 / (1 + a[i] * (x1 - 1)) ** 2) - ((1 + D ** 2) * np.exp(-(D ** 2))) * ((2 - (1 + D ** 2) * np.exp(-(D**2) / 2)))
        result[i] = -quad(lam, 0, t1[i])[0] 
        
    result = np.sum(result) 
    return result

#constraint
def cons(t1): 
    res = np.sum(t1) - ((tot)) 
    
#jacobian of the objective function: modify the gradient function to use t

def jacb(t1):
    der = np.empty(len(x0))
    for i in range(len(x0)):
        D = b[i] * np.power(c[i] / (2 * np.exp(1)), c[i]) / (np.power(np.log(t1[i]), c[i]))
        der[i] = -2 * a[i] * (1 - 1 / (1 + a[i] * (t1[i] - 1)) ** 3) * np.exp(-(D ** 2)) * ((1 + D ** 2) * np.exp(-(D ** 2) / 2) - 2)
    return der 

for fill in skip:
    indices = np.where(FillNumber == fill)[0]
    FillNumber = np.delete(FillNumber, indices)
    tf = np.delete(tf, indices)
    Lmes = np.delete(Lmes, indices)
    print('length tf', len(tf))
for i in range(len(FillNumber)):
    if FillNumber[i] in skip:
        continue

#Initial guesses    
x0=(f_rev*tf+1)

#constraint determination
tot=sum(x0)

list=[[f_rev*1800+1, f_rev*28*3600+1]]
#for i in range(len(a) - 1):
for li in range(1, len(a)):
#for li in range(1, len(x0)):
    list=list+[[f_rev*1800+1, f_rev*28*3600+1]]
    #list=list+[[(1800 / f_rev) + 1, (28 * 3600 / f_rev) + 1]]
        
bnd=list

#optimization
res = minimize(fun, x0, options={'disp': True, 'maxiter':10000}, constraints={'type':'eq', 'fun': cons, 'jac': lambda x: np.ones(len(x0))}, jac=jacb, method='SLSQP', bounds=bnd) #

#saving optimized times
with open('NumericalOptimization/res_opt_20{}.txt'.format(str(year)), 'w') as f:
        f.write('')
        f.close()
for el in res.x:
    with open('NumericalOptimization/res_opt_20{}.txt'.format(str(year)), 'a') as f:
        f.write(str(el))
        f.write('\n')

plot=True

if plot==True:
    #comparison between real and Optimised times
    fig, ax1= plt.subplots()
    ax1.hist(x0/1e8, facecolor='steelblue', density=True, alpha=0.4, label="Real Fill Times" )
    ax1.hist(res.x/1e8, color='red', histtype='step', density=True, label="Optimised Fill Times")
    ax1.set_xlabel(r'Time in Turns 'r'$\tau$(10$^8$)')
    ax1.set_ylabel(r'Normalised Frequencies')
    ax1.set_title('20{}'.format(year))
    plt.legend(loc='best')
    plt.savefig('NumericalOptimization/20{}_times.pdf'.format(year))
    
    Lmes16,Lmes17,Lmes18=ld.MeasuredLuminosity()
    
    if year==16:
        Lmes=Lmes16
    elif year==17:
        Lmes=Lmes17
    elif year==18:
        Lmes=Lmes18
    
    r=res.x
    
    Lopt=[]
    for i in range(len(res.x)):
        D = b[i] / (np.log(r[i]) ** c[i])
        fit= lambda x: (1 / (1 + a[i] * (x - 1)) ** 2) - ((1 + D ** 2) * np.exp(-(D ** 2))) * ((2 - (1 + D ** 2) * np.exp(-(D**2) / 2)))
        Li=quad(fit, 0, r[i])[0]
        Lopt.append(Li)
    
    Lopt=np.array(Lopt)
    Lmes_tot = np.sum(Lmes)
    Lopt_tot = np.sum(Lopt)


     #comparison between Actual and Optimised integrated luminosity
    fig, ax1= plt.subplots()
    ax1.hist(Lmes,  facecolor="steelblue", alpha=0.4, density=True, label="Measured Integrated Luminosities" )
    ax1.hist( Lopt/1e9, histtype='step', density=True, color='red', label="Optimised Integrated Luminosity")
    plt.plot([],[], "k.", label=r'Initial $L_{\mathrm{tot}}$='+'{:.2f}'.format(Lmes_tot)+r' [$\mathrm{fb}^{-1}$]')
    plt.plot([],[], "k.", label=r'Optimized $L_{\mathrm{tot}}$='+'{:.2f}'.format(Lopt_tot/1e9)+r' [$\mathrm{fb}^{-1}$]')
    plt.legend(loc='upper left')
    ax1.set_xlabel(r'Integrated Luminosity [$\mathrm{fb}^{-1}$]')
    ax1.set_ylabel('Normalised Frequencies')
    ax1.set_title('20{}'.format(year))
    plt.legend(loc='upper left')
    plt.savefig('NumericalOptimization/20{}_lumi.pdf'.format(year))
