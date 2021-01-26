# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11, 2018 18:06:03

@author: Joy Cheung
"""

# derive sampling matrix
from SALib.sample.morris import sample
problem = {'num_vars': 16,
           'names': ['month', 'hour', 'altitude', 'longitude', 'latitude', 
                     'kP_OH', 'kP_CO3', 'kP_1O2', 'kP_3DOM', 'fi_P', 'alpha', 'beta', 
                     'NPOC', 'CNO3', 'CNO2', 'depth'], 
           'bounds':[[0,1], [0,1], [0,3.5], [-180,179], [-59,59]] +
                     [[0,1]]*7 + 
                     [[1,40], [1.13e-5,9.84e-4], [0,3.57e-5], [0.2,2]]}
smp = sample(problem, N = 1000, num_levels = 4) # returns a numpy.ndarray containing (16+1)*N rows and 16 columns, i.e. a 17000*16 matrix

# determine the sample values for all universal variables (not specific to virus species or water type), including epsilon factor (+-10%)
from scipy.stats import randint
smp[:,0] = randint.ppf(smp[:,0], 4, 13).astype(int) # month

from timezonefinder import TimezoneFinder
import datetime, pytz, ephem, os, re
import numpy as np

tz = []

for i in range(len(smp)):   # hour
    month = smp[i,0]
    alt = smp[i,2]  # in km
    long = smp[i,3]
    lati = smp[i,4]
    # get timezone from location
    TF = TimezoneFinder()
    try:
        tzname = TF.timezone_at(lng = long, lat = lati)
        j = 1
        while tzname == None:
            tzname = TF.closest_timezone_at(lng = long, lat = lati, delta_degree = j)
            j += 1
    except ValueError:
        print("The coordinate was out of time zone bounds.")
        break
    tzone = int(datetime.datetime.now(pytz.timezone(tzname)).utcoffset().total_seconds()/60/60) # convert to UTC offset
    tz.append([tzone, tzname])  # record time zone, and time zone name for each sample location for SMARTS, now smp's shape is 17000*18
    # calculate sunrise, sunset time
    ob = ephem.Observer()
    ob.date = '2017/%s/22' % month
    ob.lat = '%s:%s' % (int(lati), abs(round((lati - int(lati))*60)))
    ob.lon = '%s:%s' % (int(long), abs(round((long - int(long))*60)))
    ob.elev = alt*1000
    s = ephem.Sun()
    srise1 = ob.next_rising(s).datetime()
    srise2 = ob.previous_rising(s).datetime()
    sset1 = ob.next_setting(s).datetime()
    sset2 = ob.previous_setting(s).datetime()
    srise = (max((srise1.hour + srise1.minute/60), (srise2.hour + srise2.minute/60)) + 0.5 + tzone) % 24
    sset = (min((sset1.hour + sset1.minute/60), (sset2.hour + sset2.minute/60)) - 0.5 + tzone) % 24
    # hour
    smp[i,1] = round(srise + smp[i,1]*(sset - srise), 1)

mmdd = datetime.date.today().strftime("%Y%m%d")

os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations")
np.savetxt("smp_matrix_%s.csv" % mmdd, smp, delimiter=",")


#%% All sunlight irradiance samples
os.makedirs("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/OUTPUT/%s_ext" % mmdd)
os.makedirs("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/OUTPUT/%s_out" % mmdd)

# get template input file for SMARTS
os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS")
with open("smarts295.inp-template.txt") as sinpfile:
    sinp = sinpfile.read()

cmRegex = re.compile(r'comment')
sprRegex = re.compile(r'!Card 2a')
atmRegex = re.compile(r'!Card 3a')
timeRegex = re.compile(r'!Card 17a')

import shutil
n = 0

for row in smp:

    n += 1

    month = row[0]
    hour = row[1]
    alt = row[2]  # in km
    long = row[3]
    lati = row[4]
    tzone = tz[n-1][0]
    tzname = tz[n-1][1]

    # to write SMARTS input files
    inp = cmRegex.sub('sample '+ str(n), sinp) # n is the index of simulation
    inp = sprRegex.sub(str(lati)+' '+str(alt)+' '+str(0), inp) # lati is latitude, alt is altitude

    if (tzname.startswith('US/') or tzname.startswith('America/')) and not(tzname.endswith('Alaska') or tzname.endswith('Samoa')):
        inp = atmRegex.sub('USSA', inp)  # U.S. standard atmosphere
    elif lati >= 50 and (tzname.endswith('Alaska') or tzname.startswith('Canada') or tzname.startswith('Iceland')):
        if month >= 3 and month <= 8:
            inp = atmRegex.sub('SAS', inp)  # sub-Arctic summer
        else:
            inp = atmRegex.sub('SAW', inp)  # sub-Arctic winter
    elif lati >= 35:
        if month >= 3 and month <= 8:
            inp = atmRegex.sub('MLS', inp)  # mid-latitude summer
        else:
            inp = atmRegex.sub('MLW', inp)  # mid-latitude winter
    elif lati <= -35:
        if month < 3 or month > 8:
            inp = atmRegex.sub('MLS', inp)  # mid-latitude summer
        else:
            inp = atmRegex.sub('MLW', inp)  # mid-latitude winter
    elif abs(lati) < 23.5:
        inp = atmRegex.sub('TRL', inp)  # tropical
    else:
        if (lati > 0 and month >=3 and month <= 8) or (lati < 0 and (month < 3 or month > 8)):
            inp = atmRegex.sub('STS', inp) # sub-tropical summer
        else:
            inp = atmRegex.sub('STW', inp) # sub-tropical winter

    inp = timeRegex.sub('2017 '+str(month)+' 22 '+str(hour)+' '+str(lati)+' '+str(long)+' '+str(tzone), inp)

    file_inp = open("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/smarts295.inp.txt", "w")
    file_inp.write(inp)
    file_inp.close()
        
    # run SMARTS with command lines
    os.system('smarts295bat.exe')

    # To rename and remove output file containing photon flux density
    for filename in os.listdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS"):
        if filename.endswith('ext.txt'):
            shutil.move("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/%s" % filename, "C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/OUTPUT/%s_ext/%s" % (mmdd, n))
        elif filename.endswith('out.txt'):
            shutil.move("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/%s" % filename, "C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/OUTPUT/%s_out/%s" % (mmdd, n))
        elif filename.endswith('inp.txt'):
            os.unlink("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/%s" % filename)

#%% simulator

def solar_virus_Morris(mmdd, scenario, problem, sample_matrix):
    
    import sys, os, re
    import pandas as pd
    import numpy as np

    sys.path.append("C:/Users/joy_c/Dropbox/MS/EPA_solar/APEX_virus_python")
    import simapex
  
    for virus, water in scenario:
        
        os.makedirs("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_Morris/%s_%s" % (mmdd, virus, water))

        # random sampling and model running for Morris screening
        
        ssize = len(sample_matrix)
        asmp = np.append(np.arange(ssize).reshape([ssize,1])+1, sample_matrix, axis = 1) # actual input matrix + sample index  = 17000*17
        asmp = np.append(asmp, np.zeros([ssize, 2]), axis = 1) # two empty columns for zenith angles and temperatures = 17000*19; column 0 is sample index, 1~16 parameters, 17: zenith; 18: temperature
        
        if water == 'WSP':
            asmp[:,11] = 0.14 + 0.08*asmp[:,11]          # alpha
            asmp[:,12] = 0.011 + 0.004*asmp[:,12]        # beta
        elif water == 'natural':
            asmp[:,11] = 0.41 + 0.08*asmp[:,11]
            asmp[:,12] = 0.013 + 0.004*asmp[:,12]
    
        if virus == 'MS2':
            asmp[:,6] = 5.5e9 + (3.8e9)*asmp[:,6]        # kP_OH
            asmp[:,7] = 1e8 + (6e7)*asmp[:,7]            # kP_CO3
            asmp[:,8] = 1.2e8 + (18.5e8)*asmp[:,8]       # kP_1O2
            asmp[:,9] = 4.7e8 + (3.6e8)*asmp[:,9]        # kP_3DOM
            asmp[:,10] = 2.56e-4 + (2.644e-3)*asmp[:,10]   # QY: 2.56e-4~2.90e-3
        elif virus == 'PhiX174':
            asmp[:,6] = 1e9 + (1.4e9)*asmp[:,6]
            asmp[:,7] = 3e7 + (6e7)* asmp[:,7]
            asmp[:,8] = 4.7e7 + (2.2e7)*asmp[:,8]
            asmp[:,9] = 1.2e7 + (1e7)*asmp[:,9]
            asmp[:,10] = 1.4e-3*(1 + 9*asmp[:,10])         # 1.4e-3~1.4e-4
        else:
            asmp[:,6] = 2.8e9 + (2.8e9)*asmp[:,6]		
            asmp[:,7] = 5e7 + (6e7)*asmp[:,7]
            asmp[:,8] = 1.8e8 + (8e7)*asmp[:,8]
            asmp[:,9] = 4e7 + (6e7)*asmp[:,9]
            asmp[:,10] = 2.5e-5*(1 + 9*asmp[:,10])         # 2.5e-5~2.5e-4
    
        Y = [] # an empty list to store output matrix

        
        for row in asmp:
        
            n, month, hour, alt, long, lati, kP_OH, kP_CO3, kP_1O2, kP_3DOM, fi_P, alpha, beta, NPOC, CNO3, CNO2, depth = row[:17]
            
            with open("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/OUTPUT/%s_out/%s" % (mmdd, int(n))) as outfile:
                out = outfile.read()
            
            # find patterns of text with regular expression and extract intermediate inputs from SMARTS output file
            zenithRegex = re.compile(r'Zenith Angle \(apparent\) =\s+(\d+\.?\d*)')
            zenithstr = zenithRegex.search(out)
            zen = float(zenithstr.group(1))
            row[-2] = zen
            temRegex = re.compile(r'Instantaneous at site\'s altitude =\s+(\d+\.?\d*) K')
            temstr = temRegex.search(out)
            temp = float(temstr.group(1))
            row[-1] = temp
        
            # to read the SMARTS output file as a dataframe
            ext = pd.read_csv("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS/OUTPUT/%s_ext/%s" % (mmdd, int(n)), sep = " ", index_col = 0, na_values = '*********')
            ext.columns = ["p0sun, Einstein/(cm2)/s/nm"]
            ext.iloc[:,[0]] = ext.iloc[:,[0]]/6.022e23        # unit conversion of photon flux density
        
            # to write photon flux density and wavelength into APEX input file
            definp = pd.read_csv("C:/Users/joy_c/Dropbox/MS/EPA_solar/APEX_virus_python/%s.csv" % virus, index_col = 0)
            definp = definp.drop("p0sun, Einstein/(12cm2)/s/nm", axis = 1).reindex(ext.index.values)
            apexinp = definp.join(ext)
            apexinp = apexinp[["ENO3, M-1 cm-1",	"ENO2, M-1 cm-1",	"phi(NO2/OH)",	"p0sun, Einstein/(cm2)/s/nm",	 "phi", "EP, M-1 cm-1",	"Aw, cm-1"]]
            apexinp.loc[:,["ENO3, M-1 cm-1", "ENO2, M-1 cm-1", "EP, M-1 cm-1"]] = apexinp.loc[:,["ENO3, M-1 cm-1", "ENO2, M-1 cm-1", "EP, M-1 cm-1"]].interpolate(method = 'akima')
            apexinp.loc[:,["phi(NO2/OH)", "phi"]] = apexinp.loc[:,["phi(NO2/OH)", "phi"]].interpolate(method = 'values')
            apexinp.loc[:,["Aw, cm-1"]] = -1
            apexinp = apexinp.fillna(0)
            with open ("C:/Users/joy_c/Dropbox/MS/EPA_solar/APEX_virus_python/inputf.csv", 'w') as f:
                apexinp.to_csv(f)    
    
            # run APEX
            os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/APEX_virus_python")
        
            # =====================================================================================================
            # *** DEFAULT DATA ***
            CBr = 1e-11     
            y_OH = 0.1      
            y_CO3 = 1e-3    
            y_1O2 = 1e-4    
            y_3DOM = 1e-5   
            y_Phot = 1e-6 
            
            if water == 'WSP':
                CCO3 = 2e-5     # from savetable.m (probably for WSP)
                CHCO3 = 4.3e-3  
                qyieldOH_CDOM = 8.09e-5
                carbonateyieldCO3_CDOM = 1e-2 	# 10 in savetable.m
                qyield1O2_CDOM = 6.63e-3
                qyieldTriplet_CDOM = 1.1e-2
                k_scav_OH = 31760
                k_scav_CO3 = 132639 
                k_quench_DOM = 3.32e6   # pseudo 1st order quenching rate constant of 3CDOM* 
            
            elif water == 'natural':
                CCO3 = 4.94e-6  # plotgraph.m
                CHCO3 = 1.04e-3 
                qyieldOH_CDOM = 3e-5
                carbonateyieldCO3_CDOM = 6.5e-3 	# The relevant equation is: R_CO3_CDOM = carbonateyieldCO3_CDOM * CCO3 * PaCDOM
                qyield1O2_CDOM = 1.25e-3
                qyieldTriplet_CDOM = 1.28e-3
                k_scav_OH = 5e4
                k_scav_CO3 = 1e2
                k_quench_DOM = 5e5
                
            # *** END OF DEFAULT DATA ***
            # ====================================================================================================
    
    
            y = simapex.apex('inputf', depth, zen, alpha, beta, CNO3, CNO2, NPOC, CCO3, CHCO3, CBr, kP_OH, kP_CO3, kP_3DOM, kP_1O2, fi_P, y_OH, y_CO3, y_1O2, y_3DOM, y_Phot, qyieldOH_CDOM, carbonateyieldCO3_CDOM, qyield1O2_CDOM, qyieldTriplet_CDOM, k_scav_OH, k_scav_CO3, k_quench_DOM)
            Y.append([n] + y)
            
            os.unlink("C:/Users/joy_c/Dropbox/MS/EPA_solar/APEX_virus_python/inputf.csv")
            
            print("Simulation %s for %s in %s water finished" % (int(n), virus, water))    
             
        os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_Morris/%s_%s" % (mmdd, virus, water))

        np.savetxt('Morrisinp_data.csv', asmp, delimiter = ",", 
                   header = "sample,month,hour,altitude,longitude,latitude,kP_OH,kP_CO3,kP_1O2,kP_3DOM,fi_P,alpha,beta,NPOC,CNO3,CNO2,depth,zenith,temperature")
        np.savetxt('Morrisout_data.csv', np.array(Y), delimiter = ',',
                   header = "sample,k_OH,k_CO3,k_1O2,k_3DOM,k_Phot,k_tot")

        
        # computation of Morris sensitivity indices
        X = asmp[:, 1:17]
        k_total = np.array(Y)[:,-1]
        #Y = list(map(list, zip(*Y)))
    
        #k_total = np.array(Y[5])
        #k_Phot = np.array(Y[4])
        #k_OH = np.array(Y[0])
        #k_1O2 = np.array(Y[2])
        from SALib.analyze.morris import analyze

        # k_total
        Si_Morris = analyze(problem, X, k_total, num_levels = 4)
        Si = pd.DataFrame.from_dict(Si_Morris)
        Si.to_csv("Si_ktot.csv") 
    
        # k of direct photolysis
        #Si_Morris = amorris.analyze(problem, X, k_Phot, num_levels = 12, grid_jump = 6)
        #Si = pd.DataFrame.from_dict(Si_Morris)
        #Si.to_csv("Si_kPhot.csv")
    
        # k of reaction with OH
        #Si_Morris = amorris.analyze(problem, X, k_OH, num_levels = 12, grid_jump = 6)
        #Si = pd.DataFrame.from_dict(Si_Morris)
        #Si.to_csv("Si_kOH.csv")
    
        # k of reaction with 1O2
        #Si_Morris = amorris.analyze(problem, X, k_1O2, num_levels = 12, grid_jump = 6)
        #Si = pd.DataFrame.from_dict(Si_Morris)
        #Si.to_csv("Si_k1O2.csv")
    

#%%
scenario = [['MS2', 'WSP'],
            ['PhiX174', 'WSP'],
            ['Adeno', 'WSP'],
            ['MS2', 'natural'],
            ['PhiX174', 'natural'],
            ['Adeno', 'natural']]

solar_virus_Morris(mmdd, scenario, problem, smp)
