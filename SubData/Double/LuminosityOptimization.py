import numpy as np
from scipy.integrate import quad
c = 299792458 #[m/s]
C_LHC = 26658.883 #[m]
B_r = 0.9999999895876449 #beta_r
G_r = 6929.637526652453 #gamma_r
S_z = 1.35e-9*c/4#[m] sigma_z longitudinal RMS dimension 
S_int = 7.95e-30 #[m^2]interaction cross section
n_c = 2 #Number of collision points

def selection_sort(x):
    """Sort an array.

    Args:
       x : array

    Returns:
       x: array
    """
    for i in range(len(x)):
        swap = i + np.argmin(x[i:])
        (x[i], x[swap]) = (x[swap], x[i])
    return x

def Parameters2016():
  """Evaluates some parameters of the Luminosity model.
   
  Returns
  -------
  Xi, Eps, S_s, Fe: any
  n_c, k_b: int
  n_i, f_rev, B_s, E_s, B_r, G_r, S_int, N_i, T_hc, T_ph, S_z : float"""
  #Defining Constants and Variables
  n_i = 1.25e11 #Intensity of the beam (only time dependent beam parameter)
  B_s = 0.4 #[m] beta* -  value of the beta-function at the collision point
  E_s = 2.2e-6 #[m] epsilon* - RMS normalized transverse emittance
  T_hc = 1.85e-4 #[rad] theta_hc - half crossing angle 
  k_b = 2220  #number of colliding bunches
  N_i = k_b*n_i #Intensity of the beam (only time dependent beam parameter)
  T_ph = 12614400 #[s] total physics time (146 days)

   #Definition of sigma* - Tranverse RMS dimension
  S_s = np.sqrt((B_s*E_s)/(B_r*G_r))

  #Definition of F(theta_c, sigma_z, sigma*) -  that accounts for the reduction in volume overlap
  # between the colliding bunches due to the presence of a crossing angle
  Fe = 1/(np.sqrt(1+(((T_hc*S_z)/(S_s))**2)))


  #Definition of the revolution frequency
  f_rev = c/(C_LHC) #Hz

  #Definition of Xi
  Xi = ((G_r*f_rev)/(4*np.pi*E_s*B_s*k_b))*Fe
  
  #Epsilon Definition
  Eps = (S_int*n_c*Xi)/f_rev
  
  return n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps

def Model_L16(t_fill, t_a):
  """Evaluates the Instantaneous luminosity, the integrated luminosity and the Total Luminosity.
    
  Parameters
  ----------
  t_fill, t_a: any 
    
  Returns
  -------
  L_inst, L_int, L_tot: any"""
  #Parameters
  n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = Parameters2016()
  
  #Instantaneous Luminosity
  L_inst=(Xi*(N_i**2))
  
  #Integrated Luminosity
  L_int = ((N_i*Xi)/(f_rev*Eps))*((N_i*Eps*f_rev*t_fill)/(1+(Eps*N_i*f_rev*t_fill)))
  
  
  #Total Luminosity  
  L_tot = (T_ph/(t_a+t_fill))*L_int
  
  return   L_inst, L_int, L_tot

def L_optimal_16(t_a):
  """Evaluates the Optimized Total Luminosity in the simplest case: the t_fill is the one that maximize the total luminosity, 
  and the t_a is given.
    
  Parameters
  ----------
  t_a: any 
    
  Returns
  -------
  t_opt, L_tot_opt, L_int_opt, t_a16, t_opt16, L_int16, L_tot16: any"""
  n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = Parameters2016()
  
  #Definition of t_opt
  t_opt = np.sqrt(t_a)/(np.sqrt(N_i*n_c*Xi*S_int)) 
  
  L_inst_opt, L_int_opt, L_tot_opt = Model_L16(t_opt, t_a)
  
  t_a16 = 8.861574692646066*3600 #expectation value 2016
  t_opt16 = np.sqrt(t_a16)/(np.sqrt(N_i*n_c*Xi*S_int))
  L_inst16, L_int16, L_tot16 = Model_L16(t_opt16, t_a16)

  return t_opt, L_tot_opt, L_int_opt, t_a16, t_opt16, L_int16, L_tot16

def Parameters2017():
  """Evaluates some parameters of the Luminosity model.
   
  Returns
  -------
  Xi, Eps, S_s, Fe: any
  n_c, k_b: int
  n_i, f_rev, B_s, E_s, B_r, G_r, S_int, N_i, T_hc, T_ph, S_z : float"""
  #Defining Constants and Variables
  n_i = 1.25e11 #Intensity of the beam (only time dependent beam parameter)
  B_s = 0.4 #[m] beta* -  value of the beta-function at the collision point
  E_s = 2.2e-6 #[m] epsilon* - RMS normalized transverse emittance
  T_hc = 1.5e-4 #[rad] theta_hc - half crossing angle 
  k_b = 2556  #number of colliding bunches
  N_i = k_b*n_i #Intensity of the beam (only time dependent beam parameter)
  T_ph = 3360*3600 #[s] total physics time (140 days)

   #Definition of sigma* - Tranverse RMS dimension
  S_s = np.sqrt((B_s*E_s)/(B_r*G_r))

  #Definition of F(theta_c, sigma_z, sigma*) -  that accounts for the reduction in volume overlap
  # between the colliding bunches due to the presence of a crossing angle
  Fe = 1/(np.sqrt(1+(((T_hc*S_z)/(S_s))**2)))


  #Definition of the revolution frequency
  f_rev = c/(C_LHC) #Hz

  #Definition of Xi
  Xi = ((G_r*f_rev)/(4*np.pi*E_s*B_s*k_b))*Fe
  
  #Epsilon Definition
  Eps = (S_int*n_c*Xi)/f_rev
  
  return n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps

def Model_L17(t_fill, t_a):
  """Evaluates the Instantaneous luminosity, the integrated luminosity and the Total Luminosity.
    
  Parameters
  ----------
  t_fill, t_a: any 
    
  Returns
  -------
  L_inst, L_int, L_tot: any"""
  #Parameters
  n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = Parameters2017()
  
  #Instantaneous Luminosity
  L_inst=(Xi*(N_i**2))
  
  #Integrated Luminosity
  L_int = ((N_i*Xi)/(f_rev*Eps))*((N_i*Eps*f_rev*t_fill)/(1+(Eps*N_i*f_rev*t_fill)))
  
  
  #Total Luminosity  
  L_tot = (T_ph/(t_a+t_fill))*L_int
  
  return   L_inst, L_int, L_tot

def L_optimal_17(t_a):
  """Evaluates the Total Luminosity in the simplest case: the t_fill is the one that maximize the total luminosity, 
  and the t_a is given.
    
  Parameters
  ----------
  t_a: any 
    
  Returns
  -------
  t_opt, L_tot_opt, L_int_opt, t_a17, t_opt17, L_int17, L_tot17: any"""
  n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = Parameters2017()
  
  #Definition of t_opt
  t_opt = np.sqrt(t_a)/(np.sqrt(N_i*n_c*Xi*S_int)) 
  
  L_inst_opt, L_int_opt, L_tot_opt = Model_L17(t_opt, t_a)
  
  t_a17 = 8.883538480869058*3600 #expectation value 2017
  t_opt17 = np.sqrt(t_a17)/(np.sqrt(N_i*n_c*Xi*S_int))
  L_inst17, L_int17, L_tot17 = Model_L17(t_opt17, t_a17)

  return t_opt, L_tot_opt, L_int_opt, t_a17, t_opt17, L_int17, L_tot17

def Parameters2018():
  """Evaluates some parameters of the Luminosity model.
   
  Returns
  -------
  Xi, Eps, S_s, Fe: any
  n_c, k_b: int
  n_i, f_rev, B_s, E_s, B_r, G_r, S_int, N_i, T_hc, T_ph, S_z : float"""
  #Defining Constants and Variables
  n_i = 1.25e11 #Intensity of the beam (only time dependent beam parameter)
  B_s = 0.3 #[m] beta* -  value of the beta-function at the collision point
  E_s = 1.9e-6 #[m] epsilon* - RMS normalized transverse emittance
  T_hc = 1.6e-4 #[rad] theta_hc - half crossing angle 
  k_b = 2556  #number of colliding bunches
  N_i = k_b*n_i #Intensity of the beam (only time dependent beam parameter)
  T_ph = 3480*3600 #[s] total physics time (145 days)
  
  #Definition of sigma* - Tranverse RMS dimension
  S_s = np.sqrt((B_s*E_s)/(B_r*G_r))

  #Definition of F(theta_c, sigma_z, sigma*) -  that accounts for the reduction in volume overlap
  # between the colliding bunches due to the presence of a crossing angle
  Fe = 1/(np.sqrt(1+(((T_hc*S_z)/(S_s))**2)))


  #Definition of the revolution frequency
  f_rev = c/(C_LHC) #Hz

  #Definition of Xi
  Xi = ((G_r*f_rev)/(4*np.pi*E_s*B_s*k_b))*Fe
  
  #Epsilon Definition
  Eps = (S_int*n_c*Xi)/f_rev
  
  return n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps

def Model_L18(t_fill, t_a):
  """Evaluates the Instantaneous luminosity, the integrated luminosity and the Total Luminosity.
    
  Parameters
  ----------
  t_fill, t_a: any 
    
  Returns
  -------
  L_inst, L_int, L_tot: any"""
  #Parameters
  n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = Parameters2018()
  
  #Instantaneous Luminosity
  L_inst=(Xi*(N_i**2))
  
  #Integrated Luminosity
  L_int = ((N_i*Xi)/(f_rev*Eps))*((N_i*Eps*f_rev*t_fill)/(1+(Eps*N_i*f_rev*t_fill)))
  
  
  #Total Luminosity  
  L_tot = (T_ph/(t_a+t_fill))*L_int
  
  return   L_inst, L_int, L_tot

def L_optimal_18(t_a):
  """Evaluates the Total Luminosity in the simplest case: the t_fill is the one that maximize the total luminosity, 
  and the t_a is given.
    
  Parameters
  ----------
  t_a: any 
    
  Returns
  -------
  t_opt, L_opt, t_a16, t_opt16, L_opt16: any"""
  n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps = Parameters2018()
  
  #Definition of t_opt
  t_opt = np.sqrt(t_a)/(np.sqrt(N_i*n_c*Xi*S_int)) 
  
  L_inst_opt, L_int_opt, L_tot_opt = Model_L18(t_opt, t_a)
  
  t_a18 = 7.219331108579183*3600 #expectation value 2018
  t_opt18 = np.sqrt(t_a18)/(np.sqrt(N_i*n_c*Xi*S_int))
  L_inst18, L_int18, L_tot18 = Model_L18(t_opt18, t_a18)

  return t_opt, L_tot_opt, L_int_opt, t_a18, t_opt18, L_int18, L_tot18

def t_opt_eval(N_i, n_c, Xi, S_int, t_a):
  """Evaluate the optimal fill time.

  Args:
      N_i(Any):
      n_c (Any): [description]
      Xi (Any): [description]
      S_int (Any): [description]
      t_a (Any): Turn Around time

  Returns:
      t_opt: optimal fill time
  """
  return np.sqrt(t_a)/(np.sqrt(N_i*n_c*Xi*S_int)) 

def t_opt_eval_data(N_i, n_c, Xi, S_int, t_a):
  """Evaluate the optimal fill time with the data model.

  Args:
      N_i(Any):
      n_c (Any): [description]
      Xi (Any): [description]
      S_int (Any): [description]
      t_a (Any): Turn Around time

  Returns:
      t_opt: optimal fill time
  """
  return np.sqrt(t_a)/(np.sqrt(N_i*n_c*Xi*S_int)) 

def tau_opt_eval(N_i,Eps, t_a):
  """Evaluate the optimal fill time.

  Args:
      N_i(Any):
      n_c (Any): [description]
      Xi (Any): [description]
      S_int (Any): [description]
      t_a (Any): Turn Around time

  Returns:
      t_opt: optimal fill time
  """
  return np.sqrt(t_a)/(np.sqrt(N_i*Eps)) 
  