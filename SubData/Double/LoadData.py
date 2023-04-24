#########################################################
#      @Giulia Faletti                                  #
#   Loading run 2 data for the analysis                 #
#########################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#old functions
def Create_DataSet():
     """Reads the Excel file and creates the corret datasets.    
    Returns
    -------
    data1, data2, data3: DataFrame"""

     #Importing datas from the excel file
     data = pd.read_excel(r'TurnAroundData.xlsx')

     df1 = pd.DataFrame(data, columns=['sample16'])
     df2 = pd.DataFrame(data, columns=['sample17'])
     df3 = pd.DataFrame(data, columns=['sample18'])

     #Setting the correct different data sets
     data1 = df1.dropna()#loc[:156]
     data2 = df2.dropna()#loc[:194]
     data3 = df3.dropna()#loc[:217]

     return data1, data2, data3

def DataToLists(data1, data2, data3):
     """Transfroms Pandas Dataframes into Lists.
    
    Parameters
    ----------
    data1: dataframe
    data2: dataframe
    data3: dataframe
    
    
    Returns
    -------
    data_16, data_17, data_18: lists"""
     
     #Converting data into python lists
     data_16 = data1['sample16'].values.tolist()
     data_17 = data2['sample17'].values.tolist()
     data_18 = data3['sample18'].values.tolist()
     
     return data_16, data_17, data_18
 
def FromListsToArrays(data_16, data_17, data_18):
     """Transfroms Lists into arrays.
    
    Parameters
    ----------
    data_16: list
    data_17: list
    data_18: list
    
    
    Returns
    -------
    array16, array17, array18: arrays"""
     
     #to array
     array16 = np.array(data_16)
     array17 = np.array(data_17)
     array18 = np.array(data_18)
     
     return array16, array17, array18
 
def TotalDataSet(data_16, data_17, data_18):
     """Creates the total dataset and transforms the inital list into a dataframe and an array.
    
    Parameters
    ----------
    data_16: list
    data_17: list
    data_18: list
    
    
    Returns
    -------
    dataTot: dataframe
    array_tot: array"""
     
     #Summing lists for the total dataset 
     data_tot = data_16 + data_17 + data_18
     dataTot = pd.DataFrame(data_tot)
     array_tot = dataTot.to_numpy()
     
     return data_tot, dataTot, array_tot
 
def PartialDataSets(data_16, data_17, data_18):
      """Creates the partial dataseta and transforms the inital lists into dataframea and arrays.
    
    Parameters
    ----------
    data_16: list
    data_17: list
    data_18: list
    
    
    Returns
    -------
    data_tot_A, data_tot_B, data_tot_C: dataframes
    array_totA, array_totB, array_totC: arrays"""
      
      #Summing lists for the partial datasets 
      data_tot_A = data_16 + data_17
      data_tot_B = data_17 + data_18
      data_tot_C = data_16 + data_18
      array_totA = np.array(data_tot_A)
      array_totB = np.array(data_tot_B)
      array_totC = np.array(data_tot_C)
      
      return data_tot_A, data_tot_B, data_tot_C, array_totA, array_totB, array_totC
 
def Data():
     """Generate the whole set of data sample needed.

     Returns:
         data_16, data_17, data_18: list
         array16, array17, array18: array
     """
     #Importing datas from the excel file
     data = pd.read_excel(r'TurnAroundData.xlsx')

     df1 = pd.DataFrame(data, columns=['sample16'])
     df2 = pd.DataFrame(data, columns=['sample17'])
     df3 = pd.DataFrame(data, columns=['sample18'])

     #Setting the correct different data sets
     data1 = df1.dropna()#loc[:156]
     data2 = df2.dropna()#loc[:194]
     data3 = df3.dropna()#loc[:217]
     
     #Converting data into python lists
     data_16 = data1['sample16'].values.tolist()
     data_17 = data2['sample17'].values.tolist()
     data_18 = data3['sample18'].values.tolist()
     
     #to array
     array16 = np.array(data_16)
     array17 = np.array(data_17)
     array18 = np.array(data_18)
     
     return data_16, data_17, data_18, array16, array17, array18
 
#data needed for online strategy    
def loadFill():
     """Reads the Excel file and creates the corret dataset lists.    
     Returns
     -------
     data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 : list"""

     #Importing datas from the excel file
     data = pd.read_excel(r'FillData.xlsx')

     df1 = pd.DataFrame(data, columns=['ta16'])
     df2 = pd.DataFrame(data, columns=['tf16'])
     df3 = pd.DataFrame(data, columns=['ta17'])
     df4 = pd.DataFrame(data, columns=['tf17'])
     df5 = pd.DataFrame(data, columns=['ta18'])
     df6 = pd.DataFrame(data, columns=['tf18'])
     

     #Setting the correct different data sets
     data1 = df1.dropna() #df1.loc[:44]
     data2 = df2.dropna()#df2.loc[:44]
     data3 = df3.dropna()#loc[:83]
     data4 = df4.dropna()#loc[:83]
     data5 = df5.dropna()#loc[:94]
     data6 = df6.dropna()#loc[:94]
     
     #Converting data into python lists
     data_ta16 = np.array(data1['ta16'].values.tolist())
     data_tf16 = np.array(data2['tf16'].values.tolist())
     data_ta17 = np.array(data3['ta17'].values.tolist())
     data_tf17 = np.array(data4['tf17'].values.tolist())
     data_ta18 = np.array(data5['ta18'].values.tolist())
     data_tf18 = np.array(data6['tf18'].values.tolist())
     
     return data_ta16, data_tf16, data_ta17,\
          data_tf17, data_ta18, data_tf18 

def Data_sec(array16, data_ta16, data_tf16, array17, data_ta17, data_tf17, array18, data_ta18, data_tf18):
     """Transform the data from hours to seconds.

     Args:
         array16, array17, array18, data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18: array

     Returns:
         data_16_sec, data_ta16_sec, data_tf16_sec, data_17_sec, data_ta17_sec, data_tf17_sec, data_18_sec, data_ta18_sec, data_tf18_sec: array
     """
     data_16_sec = array16*3600
     data_ta16_sec = data_ta16*3600 
     data_tf16_sec = data_tf16*3600  
     data_17_sec = array17*3600
     data_ta17_sec = data_ta17*3600 
     data_tf17_sec = data_tf17*3600
     data_18_sec = array18*3600
     data_ta18_sec = data_ta18*3600 
     data_tf18_sec = data_tf18*3600
     return data_16_sec, data_ta16_sec, data_tf16_sec, data_17_sec, data_ta17_sec, data_tf17_sec, data_18_sec, data_ta18_sec, data_tf18_sec
     
def FillNumber():
     """Creates the arrays that contain the number of fills' list.

     Returns:
         FillNumber16, FillNumber17, FillNumber18: array
     """
     data = pd.read_excel(r'FillData.xlsx')

     df1 = pd.DataFrame(data, columns=['NrFill_2016'])
     df2 = pd.DataFrame(data, columns=['NrFill_2017'])
     df3 = pd.DataFrame(data, columns=['NrFill_2018'])
     
     
     data1 = df1.dropna()#df1.loc[:44]
     data2 = df2.dropna()#loc[:83]
     data3 = df3.dropna()#loc[:94]
     
     NrF_16 = data1['NrFill_2016'].values.tolist()
     NrF_17 = data2['NrFill_2017'].values.tolist()
     NrF_18 = data3['NrFill_2018'].values.tolist()
     
     
     FillNumber16 = np.array(NrF_16)
     FillNumber17 = np.array(NrF_17)
     FillNumber18 = np.array(NrF_18)
     return FillNumber16, FillNumber17, FillNumber18

#evaluating the measured total integrated luminosities
def MeasuredLuminosity():
     """Evaluates the measured luminosity from ATLAS data.

     Returns:
         L_mes16, L_mes17, L_mes18: Measured Luminosity in fb^-1
     """
     
     FillNumber16, FillNumber17, FillNumber18 = FillNumber()
     L_mes16=[]
     for i in FillNumber16:
          text = str(int(i))
          f=open('ATLAS/ATLAS_summary_2016/{}_summary_ATLAS.txt'.format(text),"r")
          lines=f.readlines()
          for x in lines:
             result= float(x.split(' ')[3])
    
          L_mes16.append(result)   
     f.close()
     L_mes16=np.array(L_mes16)/1e9 
     
     L_mes17=[]
     for i in FillNumber17:
          text = str(int(i))
          f=open('ATLAS/ATLAS_summary_2017/{}_summary_ATLAS.txt'.format(text),"r")
          lines=f.readlines()
          for x in lines:
             result= float(x.split(' ')[3])
    
          L_mes17.append(result)   
     f.close()
     L_mes17=np.array(L_mes17)/1e9 
     
     L_mes18=[]
     for i in FillNumber18:
          text = str(int(i))
          f=open('ATLAS/ATLAS_summary_2018/{}_summary_ATLAS.txt'.format(text),"r")
          lines=f.readlines()
          for x in lines:
             result= float(x.split(' ')[3])
    
          L_mes18.append(result)   
     f.close()
     L_mes18=np.array(L_mes18)/1e9  #from microbarn^-1 to fb^-1
     return L_mes16, L_mes17, L_mes18
         
#saving fill numbers in txt files --> needed for ATLAS Root files data extraction                                    
def savingFillNumber_txt(year, FillNumber):
     """saving fill number arrays into txt files.

     Args:
         year (int): current year
         FillNumber (array): fill numbers of the current year
     """
     #saving data in txt files
     with open('FillData{}.txt'.format(year), 'w') as f:
          f.write('')
          f.close()
     for el in FillNumber:
          with open('FillData{}.txt'.format(year), 'a') as f:
               f.write(str(el))
               f.write(' ')
               f.write('\n')
               f.close()
               
#Script to save the FillNumbers in txt files
#FillNumber16, FillNumber17, FillNumber18=FillNumber()
#savingFillNumber_txt(16, FillNumber16)
#savingFillNumber_txt(17, FillNumber17)
#savingFillNumber_txt(18, FillNumber18)
       
#Loading New Atlas Data from Root Files
def AtlasData(text, year, grl=True):
     """
     Function that sets the data extracted from Atlas ROOT files for them to be usable 
     in optimization algorithms, allowing all data or only those with positive flags to
     be selected in the good atlas run list.
     Args:
         text (str): current fill
         year (int): current year
         grl (bool, optional): flag for to reject bad atlas lumi block. Defaults to True.

     Returns:
         T (array): times in seconds
         Li (array): istantaneous luminosity
     """

     #exctraction of data from txt files
     f=open('LumiData/20{}Data/{}_Data.txt'.format(year, text),"r")
     lines=f.readlines()
     dt=[]
     L=[]
     G=[]
     for x in lines:
          g=float(x.split('\t')[3])
          t=float(x.split('\t')[1])
          l=float(x.split('\t')[4])
          dt.append(t)
          L.append(l)
          G.append(g)
          
     #Evaluating the time in seconds
     T=[]
     t=dt[0]
     T.append(t)
     for i in range(1, len(dt)):
          t=t+dt[i]
          T.append(t) 
     
     #Evaluating the istantaneous luminosity
     Li=[]
     for i in range(len(dt)):
          li=L[i]/dt[i]
          Li.append(li)
     #print(L, Li)
     
     T=np.array(T)
     Li=np.array(Li)
     #cutting bad atlas values
     G=np.array(G)
     if grl==True:
          ind=np.where(G==0)[0]
          T=np.delete(T, ind)
          Li=np.delete(Li, ind)
           
     f.close()
     
     return  T, Li
   
#Script to plot new root data  
#T, Li=AtlasData(5257, 17, grl=True)
#T1, Li1=AtlasData(5017, 17, grl=False)
#fig, ax=plt.subplots()
#ax.plot(T, Li, "r.", markersize=4)
#ax.plot(T1, Li1, "b.", markersize=3)
#plt.show()

def CuttedData(year, text):
     """Function that read the saved cutted and fitted data from txt files.

     Args:
         year (int): current year
         text (str): current fill

     Returns:
         Times (ndarray): times in second of the current fill
         L_fit (ndarray): Luminosity in second of the current fill
         a[0],b[0],c[0],d[0]: fit coefficients
     """
     year=str(year)
     f=open('Cutting_Fitting/20{}/{}.txt'.format(year, text),"r")
     lines=f.readlines()
     Times=[]
     L_fit=[] #I added this 
     for x in lines:
          Times.append(int(x.split(' ')[0]))
          L_fit.append(float(x.split(' ')[1]))  #I added this
          
     f.close()
     Times = np.array(Times)
     L_fit = np.array(L_fit) #I added this
     f1=open('Cutting_Fitting/Fit_Coeff/FitCoefficients{}.txt'.format(year),"r")
     lines1=f1.readlines()
     a=0
     b=0
     c=0
     d=0
     for x1 in lines1:
          if str(int(x1.split(' ')[0]))==text:
               a=float(x1.split(' ')[1])
               b=float(x1.split(' ')[2])
               c=float(x1.split(' ')[3])
               d=float(x1.split(' ')[4])
          
     f1.close()

     #return Times,a,b,c,d
     return Times,L_fit,a,b,c,d

#text=str(5017)
#Times, a,b,c,d=CuttedData(16,text)
#print(a,b,c,d)

def CuttedNewData(year, text):

     """Function that read the saved cutted and fitted new data from txt files.

     Args:
         year (int): current year
         text (str): current fill

     Returns:
         Times (ndarray): times in second of the current fill
         a[0],b[0],c[0],d[0]: fit coefficients
     """
     year=str(year)
     f=open('Cutting_FittingNew/20{}/{}.txt'.format(year, text),"r")
     lines=f.readlines()
     Times=[]
     for x in lines:
          Times.append(float(x.split(' ')[0]))  
          
     f.close()
     Times = np.array(Times)
     f1=open('Cutting_FittingNew/FitCoefficients{}.txt'.format(year),"r")
     lines1=f1.readlines()
     a=0
     b=0
     c=0
     d=0
     for x1 in lines1:
          if str(int(x1.split(' ')[0]))==text:
               a=float(x1.split(' ')[1])
               b=float(x1.split(' ')[2])
               c=float(x1.split(' ')[3])
               d=float(x1.split(' ')[4])
          
     f1.close()

     
     return Times, a,b,c,d


def OldCuttedData(year, text):
     """Function that read the saved cutted and fitted data from txt files.

     Args:
         year (int): current year
         text (str): current fill

     Returns:
         Times (ndarray): times in second of the current fill
         a[0],b[0],c[0],d[0]: fit coefficients
     """
     year=str(year)
     f=open('Cutting_Fitting/20{}/{}.txt'.format(year, text),"r")
     lines=f.readlines()
     Times=[]
     for x in lines:
          Times.append(int(x.split(' ')[0]))  
          
     f.close()
     Times = np.array(Times)
     f1=open('Cutting_Fitting/FitCoefficients{}.txt'.format(year),"r")
     lines1=f1.readlines()
     a=0
     b=0
     c=0
     d=0
     for x1 in lines1:
          if str(int(x1.split(' ')[0]))==text:
               a=float(x1.split(' ')[1])
               b=float(x1.split(' ')[2])
               c=float(x1.split(' ')[3])
               d=float(x1.split(' ')[4])
          
     f1.close()
     
     a=np.array(a)
     b=np.array(b)
     c=np.array(c)
     d=np.array(d)
     L=a*np.exp(-b*Times)+c*np.exp(-d*Times)
     
     return Times, L


def AtlasCuttedData(year, text):

     """Function that read the saved cutted and fitted new data from txt files.

     Args:
         year (int): current year
         text (str): current fill

     Returns:
         Times (ndarray): times in second of the current fill
         a[0],b[0],c[0],d[0]: fit coefficients
     """
     year=str(year)
     f=open('Cutting_FittingNew/20{}/{}.txt'.format(year, text),"r")
     lines=f.readlines()
     Times=[]
     for x in lines:
          Times.append(float(x.split(' ')[0]))  
          
     f.close()
     Times = np.array(Times)
     f1=open('Cutting_FittingNew/FitCoefficients{}.txt'.format(year),"r")
     lines1=f1.readlines()
     a=0
     b=0
     c=0
     d=0
     for x1 in lines1:
          if str(int(x1.split(' ')[0]))==text:
               a=float(x1.split(' ')[1])
               b=float(x1.split(' ')[2])
               c=float(x1.split(' ')[3])
               d=float(x1.split(' ')[4])
          
     f1.close()

     a=np.array(a)
     b=np.array(b)
     c=np.array(c)
     d=np.array(d)
     L=a*np.exp(-b*Times)+c*np.exp(-d*Times)
     
     return Times, L
