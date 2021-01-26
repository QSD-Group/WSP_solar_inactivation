# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 10:03:35 2019

@author: joy_c

This script is dedicated to approximately estimate the virus removal rate of 
designed maturation pond with sunlight-mediated inactivation, accounting for 
uncertainties in photoreactivity parameters and water parameters.

It differs from the model in 'WSP virus removal.py':
    
    1) it takes into account the azimuth angle in estimation of khat.
    2) a different range of HRT is considered?
    3) maybe consider different L/W ratios as well?

"""
#%% Given apparent inactivaton rate constant, calculate relative effluent concentration profile

def N_hourly(N0, HRT, khat):
    
    """
    *************INPUT***************
    HRT: hydraulic retention time, in [day]
    khat: apparent total inactivation rate constant, in [h^-1]
    N0: system virus concentration at the beginning of the hour relative to influent concentration, i.e. N0/Ni, unitless
    
    *************OUTPUT**************
    N: effluent concentration at the end of the hour relative to influent concentration, i.e. N/Ni, unitless
    """
    from math import exp
    
    HRT = HRT*24 # convert day to hour
    
    a = 1/HRT + khat
    b = 1/HRT
    
    N = N0*exp(-a*1)+b/a*(1-exp(-a*1))
    
    return N

#%% Given a series of k in water columns of different depth, approximate estimation of khat of the maturation pond

def khat_approx(L, W, D, phi, z, k):
    
    """ 
    *************INPUT***************
    L: the bottom length of the maturation pond, in [m]
    W: the bottom width of the maturation pond, in [m]
    D: a numpy array of water depth from 0.01 to d with constant step, where d is the maximal depth of the pond, in [m]
    phi: azimuth angle of incident sunlight, in [°]
    z: zenith angle of incident sunlight, in [°]
    k: a numpy array of ktot returned by simapex_v2.apex(h), estimated across D, in [h^-1]
    
    *************OUTPUT**************
    khat: apparent k of the maturation pond, in [h^-1]
    """
    from math import pi, asin, sin, tan, cos
    import numpy as np
    
    d = max(D)                              # depth of the pond
    V = L*W*d                               # total volume of the pond
    phi = (phi % 180)/180*pi                # reposition phi to [0, 180] and convert it to radians
    z = z/180*pi                            # convert z to radians
    eta = asin(sin(z)/1.34)                 # refraction angle in water, in radians
    
    if L < W*abs(tan(phi)):
        l_arr = np.arange(0.5, L, 1)
        l_arr = np.append(l_arr, L)
        
        kA_L = []
        for l in l_arr:
            if l > d*tan(eta)*sin(phi):
                kA = (l/sin(phi) - d*tan(eta))*d*k[-1] + tan(eta)*np.trapz(y = k*D, x = D)
            else:
                kA = tan(eta)*np.trapz(y = k[D <= l/sin(phi)/tan(eta)]*D[D <= l/sin(phi)/tan(eta)], x= D[D <= l/sin(phi)/tan(eta)])
            kA_L.append(kA)
        
        
        khat = (2*abs(cos(phi))*np.trapz(y = np.array(kA_L), x = l_arr) + kA_L[-1]*(W*sin(phi)-L*abs(cos(phi))))/V
    
    else:
        w_arr = np.arange(0.5, W, 1)
        w_arr = np.append(w_arr, W)
        
        kA_W = []
        for w in w_arr:
            if w > d*tan(eta)*abs(cos(phi)):
                kA = (w/abs(cos(phi))- d*tan(eta))*d*k[-1] + tan(eta)*np.trapz(y = k*D, x = D)
            else:
                kA = tan(eta)*np.trapz(y = k[D <= w/tan(eta)/abs(cos(phi))]*D[D <= w/tan(eta)/abs(cos(phi))], x = D[D <= w/tan(eta)/abs(cos(phi))])
            kA_W.append(kA)
        
        khat = (2*abs(sin(phi))*np.trapz(y = np.array(kA_W), x = w_arr) + kA_W[-1]*(L*abs(cos(phi))-W*sin(phi)))/V

    return khat


#%% determine WSP design based on HRT and pond depth

def dsgn_WSP(Q, HRT, d, lw):
    
    """
    *************INPUT***************
    Q: flow rate of WSP, in [m^3/day]
    HRT: theoretical hydraulic retention time, in [day]
    d: pond depth, in [m]
    lw: the ratio of length to width of pond bottom, unitless

    *************OUTPUT**************
    L: the bottom length of the maturation pond, in [m]
    W: the bottom width of the maturation pond, in [m]
    """
    
    from math import sqrt
    
    V = Q*HRT
    W = sqrt(V/lw/d)
    L = W*lw
    
    return L, W


#%% determine when to stop simulation

def check(cps):
    
    """
    INPUT: a list of checkpoints, containing only two elements, usually the effluent concentration at the same time of two consecutive days
    OUTPUT: the percent difference between the checkpoints
    """
    
    from math import inf
    
    if 0 in cps:
        pdiff = inf        
    else: 
        x = min(cps)
        y = max(cps)
        pdiff = 1 - x/y
    
    return pdiff


#%% Simulations for k_tot profile across depth in each hour during daytime 
#   with uncertainty in photoreactivity and water quality

import numpy as np
import pandas as pd
import ephem, datetime, os, re, sys

gisdata = np.load("C:/Users/joy_c/Dropbox/MS/EPA_solar/GISdata/gisdata.npy")

lati = 42  # Western Illinois
long = -90
tzone = -6 # UTC time zone
month = 12  # June or Dec.

latid = int((60.008333333768 - lati)/0.008333333333)
lonid = int((long + 180.00833333377)/0.008333333333)

alt = gisdata[np.logical_and(gisdata[:, 0] == latid, gisdata[:, 1] == lonid)][0, 2]/1000

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

os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS")
with open("smarts295.inp-template.txt") as sinpfile:
    sinp = sinpfile.read()
    
cmRegex = re.compile(r'comment')
sprRegex = re.compile(r'!Card 2a')
atmRegex = re.compile(r'!Card 3a')
timeRegex = re.compile(r'!Card 17a')

# generate sunlight spectra in each hour from SMARTS

t = srise
zenith = []
azimuth = []
hour = []

while t <= sset - 1:
    
    hour.append(t)
    
    os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/SMARTS")
        
    for file in ["smarts295.ext.txt",
                 "smarts295.inp.txt",
                 "smarts295.out.txt"]:
        try:
            os.unlink(file)
        except FileNotFoundError:
            continue
    
    inp = cmRegex.sub('hour: '+ str(int(t)), sinp) # n is the index of simulation
    inp = sprRegex.sub(str(lati)+' '+str(alt)+' '+str(0), inp) # lati is latitude, alt is altitude
    inp = atmRegex.sub('USSA', inp)  # U.S. standard atmosphere
    inp = timeRegex.sub('2017 '+str(month)+' 22 '+str(t)+' '+str(lati)+' '+str(long)+' '+str(tzone), inp)

    file_inp = open("smarts295.inp.txt", "w")
    file_inp.write(inp)
    file_inp.close()
    
    os.system('smarts295bat.exe')
        
    # to read the SMARTS output file "smarts295.out.txt" as a string value
    with open("smarts295.out.txt") as outfile:
        out = outfile.read()
    
    # find patterns of text with regular expression and extract intermediate inputs from SMARTS output file
    zenithRegex = re.compile(r'Zenith Angle \(apparent\) =\s+(\d+\.?\d*)')
    zenithstr = zenithRegex.search(out)
    zen = float(zenithstr.group(1))
    zenith.append(zen)
    
    azimuthRegex = re.compile(r'Azimuth \(from North\) =\s+(\d+\.?\d*)')
    azimuthstr = azimuthRegex.search(out)
    azm = float(azimuthstr.group(1))
    azimuth.append(azm)
    
    # to read the SMARTS output file "smarts295.ext.txt" as a dataframe
    ext = pd.read_csv("smarts295.ext.txt", sep = " ", index_col = 0, na_values = '*********')
    ext.columns = ["p0sun, Einstein/(cm2)/s/nm"]
    ext.iloc[:,[0]] = ext.iloc[:,[0]]/6.022e23        # unit conversion of photon flux density
    
    os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/APEX_virus_python")
    
    definp = pd.read_csv("MS2.csv", index_col = 0)    
    definp = definp.drop("p0sun, Einstein/(12cm2)/s/nm", axis = 1).reindex(ext.index.values)
    apexinp = definp.join(ext)
    apexinp = apexinp[["ENO3, M-1 cm-1",	"ENO2, M-1 cm-1",	"phi(NO2/OH)",	"p0sun, Einstein/(cm2)/s/nm",	 "phi", "EP, M-1 cm-1",	"Aw, cm-1"]]
    apexinp.loc[:,["ENO3, M-1 cm-1", "ENO2, M-1 cm-1", "EP, M-1 cm-1"]] = apexinp.loc[:,["ENO3, M-1 cm-1", "ENO2, M-1 cm-1", "EP, M-1 cm-1"]].interpolate(method = 'akima')
    apexinp.loc[:,["phi(NO2/OH)", "phi"]] = apexinp.loc[:,["phi(NO2/OH)", "phi"]].interpolate(method = 'values')
    apexinp.loc[:,["Aw, cm-1"]] = -1
    apexinp = apexinp.fillna(0)
    with open ("MS2_%s_%s.csv" % (month, round(t)), 'w') as f:
        apexinp.to_csv(f)
    
    t += 1

#%%
sys.path.append("C:/Users/joy_c/Dropbox/MS/EPA_solar/codes")
import simapex_v2

smp = np.random.rand(2000, 10)
mmdd = datetime.date.today().strftime("%Y%m%d")
os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations")
np.savetxt("smp_matrix_%s.csv" % mmdd, smp, delimiter=",")

depth = np.arange(0.01, 1.51, 0.01)
l = len(depth)

os.makedirs("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_WSPrmv/samples" % mmdd)

# =====================================================================================================
# *** DEFAULT DATA FOR APEX ***
CBr = 1e-11     
y_OH = 0.1      
y_CO3 = 1e-3    
y_1O2 = 1e-4    
y_3DOM = 1e-5   
y_Phot = 1e-6 
        
CCO3 = 2e-5     # from savetable.m (for WSP)
CHCO3 = 4.3e-3  
qyieldOH_CDOM = 8.09e-5
carbonateyieldCO3_CDOM = 1e-2 	# 10 in savetable.m
qyield1O2_CDOM = 6.63e-3
qyieldTriplet_CDOM = 1.1e-2
k_scav_OH = 31760
k_scav_CO3 = 132639 
k_quench_DOM = 3.32e6   # pseudo 1st order quenching rate constant of 3CDOM* 
        
# *** END OF DEFAULT DATA ***
# ====================================================================================================

from scipy.stats import triang

n = 0
for sample in smp:
    
    n += 1
    
    # MS2
    kP_OH = 5.5e9 + (3.8e9)*sample[0]   
    kP_CO3 = 1e8 + (6e7)*sample[1]
    kP_1O2 = 1.2e8 + (18.5e8)*sample[2]	
    kP_DOM = 4.7e8 + (3.6e8)*sample[3]
    fi_P = 2.56e-4 + (2.644e-3)*sample[4]
    
    # WSP water
    alpha = 0.14 + 0.08*sample[5]        
    beta = 0.011 + 0.004*sample[6]

    NPOC = 1 + 39*sample[7]
    CNO3 = 1.13e-5 + 9.73e-4*sample[8]
    CNO2 = triang(c = 0.2, loc = 0, scale = 3.57e-5).ppf(sample[9])
    
    
    rsl = pd.DataFrame(depth, columns = ['depth'])
    
    for j in range(len(hour)):
        
        file_prefix = "C:/Users/joy_c/Dropbox/MS/EPA_solar/APEX_virus_python/MS2_%s_%s" % (month, round(hour[j]))
        z = zenith[j]
        kt = list(map(simapex_v2.apex, depth, [file_prefix]*l, [zen]*l, [alpha]*l, [beta]*l, 
                      [CNO3]*l, [CNO2]*l, [NPOC]*l, [CCO3]*l, [CHCO3]*l, [CBr]*l, [kP_OH]*l, [kP_CO3]*l, [kP_DOM]*l, 
                      [kP_1O2]*l, [fi_P]*l, [y_OH]*l, [y_CO3]*l, [y_1O2]*l, [y_3DOM]*l, [y_Phot]*l, 
                      [qyieldOH_CDOM]*l, [carbonateyieldCO3_CDOM]*l, [qyield1O2_CDOM]*l, [qyieldTriplet_CDOM]*l, 
                      [k_scav_OH]*l, [k_scav_CO3]*l, [k_quench_DOM]*l))
        kt = pd.DataFrame(kt, columns = ['%s' % round(hour[j])])
        rsl = rsl.join(kt)
    
    os.chdir("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_WSPrmv/samples" % mmdd)
    with open ("smp%s_m%s.csv" % (n, month), 'x') as f:
        rsl.to_csv(f)
            
#%% simulation for daily profile effluent virus concentration until it reaches equilibrium
import pandas as pd
import numpy as np
from math import log10

Q = 1e4  # treatment capacity, in [m3/day]
hour = list(map(round, hour))

HRT_arr = range(1, 6)    # in [days]
d_arr = np.arange(0.4, 1.51, 0.1)    # in [m]
lw_arr = [1, 2, 4, 8, 14, 20]

for HRT in HRT_arr:
    for d in d_arr:
        for lw in lw_arr:
            #lw = 1/lw
            L, W = dsgn_WSP(Q, HRT, d, lw)
            export_data = pd.DataFrame(columns = ['sample', 'eqlbrm_day', 'max_log_removal', 'min_log_removal'])
            # begin simulation of daily profile for each sample of k_tot
            for n in range(1, len(smp)+1):
                
                dt = pd.read_csv("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_WSPrmv/samples/smp%s_m%s.csv" % (mmdd, n, month), index_col = 0)
                dt = dt[dt['depth'] <= d]
                D = np.array(dt['depth'])
                
                khat_daily = {}
                
                for i in range(len(hour)):
                    k = np.array(dt['%s' % hour[i]])
                    khat_hour = khat_approx(L, W, D, azimuth[i], zenith[i], k)
                    khat_daily[hour[i]] = khat_hour
                    
                N = 1
                t = 0
                day = 1
                N_profile = [N]
                checknoon = [0, 0]
                checknight = [0, 0]
            
                while check(checknoon) > 0.001 or check(checknight) > 0.001:
                
                    h = t % 24
                    day = 1 + t//24
                    
                    if h not in hour:
                        khat = 0
                    else:
                        khat = khat_daily[h]
                
                    N = N_hourly(N, HRT, khat)
                    N_profile.append(N)
                
                    if h == 12:
                        checknoon.append(N)
                        checknoon.pop(0)
                    elif h == 0:
                        checknight.append(N)
                        checknight.pop(0)
                        
                    t += 1
            
                Ne_min = min(N_profile[-24:])
                Ne_max = max(N_profile[-24:])
                
                #Ne_midnight = checknight[1]
                #minh = h - 23 + N_profile[-24:].index(Ne_min)
                #maxh = h - 23 + N_profile[-24:].index(Ne_max)
                #if minh < 0:
                 #   minh += 24
                #if maxh < 0:
                 #   maxh += 24
            
                lgrmv_max = -log10(Ne_min)
                lgrmv_min = -log10(Ne_max)
            
            
                export_data = export_data.append(pd.Series([n, day, lgrmv_max, lgrmv_min], index = export_data.columns), ignore_index = True)
                print(n)
                
            with open("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_WSPrmv/results/winter/HRT%s_d%s_lw%s.csv" % (mmdd, HRT, round(d*100), lw), 'x') as f:
                export_data.to_csv(f, index = False)

#%% plot the profile of N
import matplotlib.pyplot as plt

fig = plt.figure(1, figsize = (5, 4))
T = np.arange(len(N_profile))
plt.plot(T, N_profile)
plt.xlabel('Hour')
plt.ylabel('$N_{t}/N_{in}$')
plt.ylim(bottom = 0, top = 1)

fig.savefig('C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_WSPrmv/results/N_profile_winter.pdf' % mmdd)

#%% explore the relationship between surface area and depth and HRT
import pandas as pd
import numpy as np

Q = 1e4  # treatment capacity, in [m3/day]
lw = 1.6

HRT = np.arange(5, 16, 2.5)
depth = np.arange(1.5, 0.3, -0.1)

area = pd.DataFrame(columns = ['5d', '7.5d', '10d', '12.5d', '15d'])

for d in depth:
    area_list = []
    for hrt in HRT:
        L, W = dsgn_WSP(Q, hrt, d, lw)
        S = L*W
        area_list.append(S)
    area = area.append(pd.Series(area_list, index = area.columns), ignore_index = True)

area.set_index(depth, inplace = True)

with open("C:/Users/joy_c/Dropbox/MS/EPA_solar/Simulations/%s_WSPrmv/results/WSP_area.csv" % mmdd, "x") as f:
    area.to_csv(f)

