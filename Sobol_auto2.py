# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 17:13:28 2018

@author: joy_c
"""

scenario = [['MS2', 'WSP'],
            ['PhiX174', 'WSP'],
            ['Adeno', 'WSP'],
            ['MS2', 'natural'],
            ['PhiX174', 'natural'],
            ['Adeno', 'natural']]

#%% derive sampling matrix
from SALib.sample.saltelli import sample
problem = {'num_vars': 14,
           'names': ['month', 'hour', 'location', 'kP_OH', 'kP_CO3', 'kP_1O2', 'kP_3DOM', 'fi_P', 'alpha', 'beta', 
                     'NPOC', 'CNO3', 'CNO2', 'depth'], 
           'bounds':[[0,1]]*10 + [[1,40], [1.13e-5,9.84e-4], [0,3.57e-5], [0.2,2]]}
smp = sample(problem, N = 10000) # returns a numpy.ndarray containing (14+1)*2N rows and 14 columns, i.e. a 300000*14 matrix

# determine the sample values for all universal variables (not specific to virus species or water type), including epsilon factor (+-10%)
from scipy.stats import randint
smp[:,0] = randint.ppf(smp[:,0], 1, 13).astype(int) # month

from timezonefinder import TimezoneFinder
import datetime, pytz, ephem, os, re
import numpy as np

tz = []
smp = np.insert(smp, 3, 0, axis = 1) # add a column for longitude
smp = np.insert(smp, 4, 0, axis = 1) # add a column for latitude

gisdata = np.load("C:/Users/joy_c/Dropbox/MS/EPA_solar/GISdata/gisdata.npy")
gisdata = gisdata[abs(60.008333333768-0.008333333333*gisdata[:, 0]) <= 55]

for i in range(len(smp)):   

    geo = int(round(len(gisdata)*smp[i,2]))
    lati = round(- gisdata[geo,0]*0.008333333333 + 55, 1)
    long = round(gisdata[geo,1]*0.008333333333 - 180.00833333377, 1)
    alt = gisdata[geo,2]/1000    # in km
    
    month = smp[i,0]
    smp[i,2] = alt
    smp[i,3] = long
    smp[i,4] = lati

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
    tz.append([tzone, tzname])  # record location, time zone, and time zone name of each sample for SMARTS    
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
np.savetxt("smp_matrix_%s.csv" % mmdd, smp, delimiter=",")   # 300,000*16 matrix, with sample values of all universal variables
smp = np.genfromtxt("smp_matrix_%s.csv" % mmdd, delimiter = ",")

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


#%% function that export SI results

def Sobol_write(problem, y, name):
    
    # problem is the same as the one for saltelli.sample()
    # y should be a 1-dimensional numpy array of model outputs
    # name (str) is the name of the output variable y
    
    import pandas as pd
    from SALib.analyze.sobol import analyze

    Si = analyze(problem, y, print_to_console = True)
    
    first_order = pd.DataFrame({'Parameter': problem["names"], 
                                'S1': Si["S1"],
                                'S1_conf' : Si["S1_conf"],
                                'ST': Si["ST"],
                                'ST_conf': Si["ST_conf"]})
    first_order.to_csv("S1&T_%s.csv" % name)

    scnd_order = pd.DataFrame(Si["S2"], columns = problem["names"], index = problem["names"])
    scnd_order.to_csv("S2_%s.csv" % name)
    
    scnd_conf = pd.DataFrame(Si["S2_conf"], columns = problem["names"], index = problem["names"])
    scnd_conf.to_csv("S2_conf_%s.csv" % name)

#%% simulator

def solar_virus_Sobol(mmdd, scenario, problem, sample_matrix):
    
    import sys, os, re
    import pandas as pd
    import numpy as np

    sys.path.append("C:/Users/joy_c/Dropbox/MS/EPA_solar/APEX_virus_python")
    import simapex

    for virus, water in scenario:
        
        os.makedirs("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_Sobol/%s_%s" % (mmdd, virus, water))

        ssize = len(sample_matrix)
        asmp = np.append(np.arange(ssize).reshape([ssize,1])+1, sample_matrix, axis = 1) # actual input matrix + sample index  = 300,000*17
        asmp = np.append(asmp, np.zeros([ssize, 2]), axis = 1) # two empty columns for zenith angles and temperatures = 300,000*19; column 0 is sample index, 1~16 parameters, 17: zenith; 18: temperature

        # random sampling and model running for Sobol' analysis
        
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
                        
            # to read the SMARTS output file "smarts295.out.txt" as a string value
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
             
        os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_Sobol/%s_%s" % (mmdd, virus, water))

        np.savetxt('Sobolinp_data.csv', asmp, delimiter = ",", 
                   header = "sample,month,hour,altitude,longitude,latitude,kP_OH,kP_CO3,kP_1O2,kP_3DOM,fi_P,alpha,beta,NPOC,CNO3,CNO2,depth,zenith,temperature")
        np.savetxt('Sobolout_data.csv', np.array(Y), delimiter = ',',
                   header = "sample,k_OH,k_CO3,k_1O2,k_3DOM,k_Phot,k_tot")

        k_total = np.array(Y)[:,-1]

        Sobol_write(problem, k_total, "ktot")

#%%
solar_virus_Sobol(mmdd, scenario, problem, smp)

