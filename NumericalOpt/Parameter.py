import numpy as np
import LoadData as ld
import matplotlib.pyplot as plt

# Selecting Current Year
years = [16, 17, 18]

# create empty arrays to store fill numbers and parameters
fill_nums = []
EPS = []
B = []
K = []
R = []

for year in years:

    
    with open('MODEL2/20{}/Coeff/mod2_Coeff_20{}.txt'.format(str(year),str(year)), 'r') as f:  

        lines=f.readlines()
        fn=[] #fill number
        eps_ni=[]
        B=[]
        k=[]
        Eps_err=[]
        B_err=[]
        K_err=[]
      
        for x in lines:
            fn.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[1]))
            B.append(float(x.split(' ')[2]))
            k.append(float(x.split(' ')[3]))  
            Eps_err.append(float(x.split(' ')[4]))
            B_err.append(float(x.split(' ')[5]))
            K_err.append(float(x.split(' ')[6]))
       
        fn = np.array(fn)
        eps = np.array(eps_ni)
        b = np.array(B)
        k = np.array(k)
        eps_err = np.array(Eps_err)
        b_err = np.array(B_err)
        k_err = np.array(K_err)
    
    
        with open('MODEL2/20{}/Coeff/mod2_Rajd_20{}.txt'.format(str(year),str(year)), 'r') as f:
            
            lines=f.readlines()
            fn1=[] #fill number
            Radj=[]
            RChi=[]

            for x in lines:
                fn1.append(float(x.split(' ')[0]))
                Radj.append(float(x.split(' ')[1]))
                RChi.append(float(x.split(' ')[2]))
          
        fn1 = np.array(fn1)
        Rs = np.array(Radj)
        chi = np.array(RChi)  
        
        # plot B vs. fill number for different years
    plt.figure(figsize=(8, 6))

    colors = ['r', 'g', 'b']  # define different colors for different years

    for i, year in enumerate(years):
        if year==16:
           fn16=fn[i]
        elif year==17:
           fn17=fn[i]
           previous_year=16
        elif year==18:
           fn18=fn[i]
           previous_year=17 
        plt.plot(fn[i], B[i], color=colors[i], label='20{}'.format(year))
        print('length 16 ', len(fn16))
        print('length 17 ', len(fn17))
        print('length 18 ', len(fn18))

    plt.xlabel('Fill Number')
    plt.ylabel('B')
    plt.title('B vs. Fill Number for Different Years')
    plt.legend()
    plt.show()
        
    
        
    fill_nums.append(fn)
    EPS.append(eps)
    B.append(b)
    K.append(k)
    #R.append(Rs)
    
    











