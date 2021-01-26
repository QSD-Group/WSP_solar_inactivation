# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 23:13:44 2018

@author: joy_c
An adaptation of simapex.py for WSP virus removal rate estimation
Time unit is converted from SSD to hour according to 1 SSD = 10 hrs
Only variables of interest are calculated
"""

def apex(depth, file_prefix, zen, alpha, beta, CNO3, CNO2, NPOC, CCO3, CHCO3, CBr, kP_OH, kP_CO3, kP_DOM, kP_1O2, fi_P, y_OH, y_CO3, y_1O2, y_3DOM, y_Phot, qyieldOH_CDOM, carbonateyieldCO3_CDOM, qyield1O2_CDOM, qyieldTriplet_CDOM, k_scav_OH, k_scav_CO3, k_quench_DOM):
    
    # preparations
    import pandas as pd
    import math
    import numpy as np
    
    # correct optical path length
    d = depth/(1-(math.sin(zen/180*math.pi)/1.34)**2)**0.5
    
    # define IC, CP and V
    IC = CCO3 + CHCO3 # total concentration of inorganic carbon is the sum of CO3 concentration and HCO3 concentration
    CP = 1e-15 # assumed concentration of substrate P (doesn't matter; won't influence the output)
    V = 0.1 * d # volume?
    
    # to read the input file and extract input variables
    inputf = '%s.csv' % file_prefix
    ALL = pd.read_csv(inputf) # read the input file and store it as a dataframe called ALL
    LL = list(ALL.iloc[:, 0]) # locate and return the wavelength column
    ENO3 = list(ALL.iloc[:, 1]) # locate and return the column for molar absorption coefficient of nitrate, in M^(-1)cm^(-1)
    ENO2 = list(ALL.iloc[:, 2]) # locate and return the column for molar absorption coefficient of nitrite, in M^(-1)cm^(-1)
    phi_NO2_OH = list(ALL.iloc[:, 3]) # locate and return the column for quantum yield of .OH generation by nitrite, unitless
    p0sun = list(ALL.iloc[:, 4]) # locate and return the column for incident photon flux of sunlight at the water surface, in Einstein cm^(-2)s^(-1)nm^(-1)
    phi = list(ALL.iloc[:, 5]) # locate and return the column for photolysis quantum yield of the substrate, unitless
    EP = list(ALL.iloc[:, 6]) # locate and return the column for molar absorption coefficient of the substrate, in M^(-1)cm^(-1)
    Aw = list(ALL.iloc[:, 7]) # locate and return the column for water absorbance over a 1 cm path length, in cm^(-1)
    
    # to create lists of the same length as LL to be updated
    RL = [0]*len(LL)
    paTOT = [0]*len(LL)
    paDOM = [0]*len(LL)
    paNO3 = [0]*len(LL)
    paNO2 = [0]*len(LL)
    RNO2 = [0]*len(LL)
    paP = [0]*len(LL)
    
    # to calculate absorption spectra
    phi_neg_count = 0
    
    for i in range(0, len(LL)):
        
        ANO3 = ENO3[i]*d*CNO3*100 # absorbance of nitrate in water
        ANO2 = ENO2[i]*d*CNO2*100 # absorbance of nitrite in water
        
        # absorbance of CDOM in water
        if Aw[i] > 0:
            ADOM = Aw[i]*d*100
        else:
            ADOM = NPOC*alpha*math.exp(-beta*LL[i])*d*100
            
        # photolysis quantum yield of the substrate; if it's not given in the input file, use uniform default value for all wavelength
        if phi[i] < 0:
            phi[i] = fi_P
            phi_neg_count += 1
            
        # absorbance of substrate in water
        if EP[i] > 0:
            AP = EP[i]*d*CP*100
        else:
            AP = 0.0
        
        # to calculate total absorbance of relevant chemicals/viruses in water, and update their absorbption spectra
        ATOT = ANO3 + ANO2 + ADOM + AP
        paTOT[i] = p0sun[i]*(1-10**(-ATOT))  
        paDOM[i] = paTOT[i]*ADOM/ATOT
        paNO3[i] = paTOT[i]*ANO3/ATOT
        paNO2[i] = paTOT[i]*ANO2/ATOT
        paP[i] = paTOT[i]*AP/ATOT
        RL[i] = paP[i]*phi[i] # photolysis rate of the substrate at each wavelength
        RNO2[i] = paNO2[i]*phi_NO2_OH[i] # formation rate of .OH by the nitrite
        
    
    # if photolysis quantum yields of the substrate for all wavelength are given in input file, rewrite the default value to 1
    if phi_neg_count == 0:
        fi_P = 1
            
    # to calculate the absorbed quantum flux and rate of photolysis of the substrate
    Pa_NO3 = np.trapz(paNO3, x = LL)
    Pa_NO2 = np.trapz(paNO2, x = LL)
    Pa_DOM = np.trapz(paDOM, x = LL)
    Pa_TOT = np.trapz(paTOT, x = LL)
    RPhot = np.trapz(RL, x = LL)
            
    # to calculate .OH production rate of sensitizers and scanverging rate
    ROH_NO3 = Pa_NO3*0.043*(IC + 0.0075)/(2.25*IC + 0.0075)
    ROH_NO2 = np.trapz(RNO2, x = LL)
    ROH_DOM = Pa_DOM*qyieldOH_CDOM
    ROH_TOT = ROH_NO3 + ROH_NO2 + ROH_DOM
    Scav = k_scav_OH*NPOC + (8.5e6)*CHCO3 + (3.9e8)*CCO3 + (1.1e10)*CBr + (1e10)*CNO2
    
    # half lifetime (h) and 1st-order inactivation coefficient (h^-1) of the substrate due to reaction with .OH
    if kP_OH > 0:
        t_OH = (1.897e-5)*d*Scav/(ROH_TOT*kP_OH) 
        k_OH = math.log(2)/t_OH
    else:
        k_OH = 0
        #t_OH = math.inf
    

    RCO3NO3 = ROH_NO3*11.3*IC/(IC + 0.061)
    RCO3 = RCO3NO3 + carbonateyieldCO3_CDOM*CCO3*Pa_DOM + ROH_TOT*((8.5e6)*CHCO3 + (3.9e8)*CCO3)/Scav
    
    # half lifetime (h) and 1st-order inactivation coefficient (h^-1) of the substrate due to reaction with CO3
    if kP_CO3 > 0:
        t_CO3 = (1.929e-5)*d*NPOC*k_scav_CO3/(RCO3*kP_CO3)    
        k_CO3 = math.log(2)/t_CO3
    else:
        k_CO3 = 0
        #t_CO3 = math.inf
        
    # half lifetime (h) and 1st-order inactivation coefficient (h^-1) of the substrate due to reaction with 3CDOM
    R3DOM = Pa_DOM*qyieldTriplet_CDOM
    if kP_DOM > 0:
        t_3DOM = (1.929e-5)*d*k_quench_DOM/(R3DOM*kP_DOM)  # the scavenging rate constant of 3CDOM* here is also different from default in natural water
        k_3DOM = math.log(2)/t_3DOM
    else:
        k_3DOM = 0
        #t_3DOM = math.inf

    # half lifetime (h) and 1st-order inactivation coefficient (h^-1) of the substrate due to reaction with 1O2
    R1O2 = Pa_DOM*qyield1O2_CDOM
    if kP_1O2 > 0:
        t_1O2 = (1.929e-5)*d*(2.5e5)/(R1O2*kP_1O2)
        k_1O2 = math.log(2)/t_1O2
    else:
        k_1O2 = 0
        #t_1O2 = math.inf
    
    # half lifetime (h) and 1st-order inactivation coefficient (h^-1) of the substrate due to direct photolysis
    if fi_P > 0:
        #t_Phot = math.log(2)*V*CP/(3600*RPhot)
        k_Phot = 3600*RPhot/(V*CP)
    else:
        k_Phot = 0
        #t_Phot = math.inf
    
    # total inactivation coefficient (h^-1) and half lifetime (h)
    k_tot = k_OH + k_CO3 + k_3DOM + k_1O2 + k_Phot
    #t_tot = math.log(2)/k_tot
    
    # steady-state concentration of PPRIs
    #coOH = ROH_TOT/(V*Scav)
    #coCO3 = RCO3/(V*k_scav_CO3*NPOC)   # why 3000? Should it be the scavenging rate constant of CO3-?
    #co1O2 = R1O2/(V*2.5e5)
    #co3DOM = R3DOM/(V*k_quench_DOM)   # why 2.71e6? Should it be the quenching rate constant of 3CDOM*?
    
    # intermediate formation rate constants
    #f_OH = y_OH * k_OH
    #f_CO3 = y_CO3 * k_CO3
    #f_1O2 = y_1O2 * k_1O2
    #f_3DOM = y_3DOM * k_3DOM
    #f_Phot = y_Phot * k_Phot
    #f_tot = f_OH + f_CO3 + f_1O2 + f_3DOM + f_Phot
    
    # overall formation yield of intermediate from P
    #y_tot = f_tot/k_tot

    # contributions of each mechanisms to inactivation of P
    #role_OH_P = k_OH/k_tot
    #role_CO3_P = k_CO3/k_tot
    #role_1O2_P = k_1O2/k_tot
    #role_3DOM_P = k_3DOM/k_tot
    #role_Phot_P = k_Phot/k_tot

    # contributions of each mechanisms to the formation of intermediates
    #role_OH_I = f_OH/f_tot
    #role_CO3_I = f_CO3/f_tot
    #role_1O2_I = f_1O2/f_tot
    #role_3DOM_I = f_3DOM/f_tot
    #role_Phot_I = f_Phot/f_tot
    
    # contrbution of each chemicals to the formationo of .OH
    #NO3_OH = ROH_NO3/ROH_TOT
    #NO2_OH = ROH_NO2/ROH_TOT
    #DOM_OH = ROH_DOM/ROH_TOT
    
    return k_tot