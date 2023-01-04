""" 
The current state of passive sampling leads to a variety of methods and tools being used to estimate sampling rates
Considering the number of assumptions already needing to be made (constant water concentration, consistent flow effects etc.) it would be benficial for the field to have a consistent approach
This package attempts to define functions for estimating sampling rates in a simple and concise manner

These equations are for passive sampler designs which do not expect to have significant effects from flow rate (i.e. the water boundary layer [WBL] is negated by the diffusive layer)
"""
# Import necessary libraries/packages/modules

from numpy import exp
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns
import math


# These functions deal with the curvilinear period of the uptake curve, this requires curve fitting a model to estimate sorbent water coefficient (ksw) and sampling rate (rs)
def nonlin(rs, time, ksw, sorbent_mass): 
    # create the full equation for the non/curvilinear phase for a passive sampler with limited WBL intereference
    return ksw * sorbent_mass * (1 - exp(-(rs * time)/(ksw * sorbent_mass)))


def plateau(ksw, sorbent_mass):
    # CHECK THIS!!!
    return ksw * sorbent_mass


# estimate sampling rate for the curvilinear phase for a passive sampler with limited WBL intereference
# this uses the nonlin function and a curve fitting method from scipy.optimize
def two_phase_nonlin_fit(df, time_column, compound_column, time_unit = 'day', water_unit = 'mL', plot = False): 
    
    params, covs = curve_fit(nonlin, df.loc[:, time_column], df.loc[:, compound_column])
    ksw, sampling_rate = params[0], params[1]
    
    if plot == True:
        print(ksw, str(sampling_rate + water_unit + "/" + time_unit))
        plot_range = np.arange(min(time), max(time))
        fig, ax = plt.subplots()
        ax.plot(time_column, compound_column, 'ko', label="y-original")
        ax.plot(plot_range, nonlin(plot_range, *params), label="a*M*(1-exp(-(r*x)/(a*M)))", color = 'black')
        plt.xlabel('Day')
        plt.ylabel('ng per sampler / ng per mL water')
        plt.legend(loc='best', fancybox=True, shadow=True)
        plt.grid(True)
        plt.show()
    return ksw, sampling_rate


# Calculate the t1/2 i.e., you can deploy the sampler for this compound for 2 times this value to reach equilibrium
# unit agnostic
def half_time_equi(ksw, sampling_rate, sorbent_mass):
    return (math.log(2)*sorbent_mass*ksw)/sampling_rate


# These functions deal with the kinetic period of the uptake curve, here a linear regression is an appropriate approximation of the rs
def two_phase_kinetic(df, time, compound, time_unit = 'day', water_unit = 'mL', plot = False): 
    # estimate sampling rate for a passive sampler in the kinetic phase (with limited WBL interference)
    # This is essentially a simple wrapper around scipy's linear regression, slope is an approximation of sampling rate in a linear system
    if len(compound) == 0 or len(time) == 0:
        raise ValueError("Inputs must not be empty.")
    if len(compound) != len(time):
        raise ValueError(f"The two arrays need to match in length, your time array is {len(time)} long while your compound array is {len(compound)}")
    #rs = ns/(cw*t)
    results = sp.stats.linregress(time, compound)
    print(f"Sampling rate is {results[0]:.3f} ± {results[4]:.3f}{water_unit}/{time_unit}\np-value is {results[3]}\nR\u00B2 is {results[2]:.3f}")

    if plot == True:
        sns.lmplot(data = df, x = time, y = compound).set(title = compound)
        plt.show()
    

def kinetic_plot(time, compound, time_unit = 'day', water_unit = 'mL', plot = False):
    if len(compound) == 0 or len(time) == 0:
        raise ValueError("Inputs must not be empty.")
    if len(compound) != len(time):
        raise ValueError(f"The two arrays need to match in length, your time array is {len(time)} long while your compound array is {len(compound)}")
    #rs = ns/(cw*t)
    results = sp.stats.linregress(time, compound)
    print(f"Sampling rate is {results[0]:.3f} ± {results[4]:.3f}{water_unit}/{time_unit}\np-value is {results[3]}\nR\u00B2 is {results[2]:.3f}")
    if plot == True:



# This section deals with comparing the fits of models i.e. it attempts to quantify whether the dataset is in the kinetic or non-linear phase
def compare_fits(): 
    # compare a linear to a non-linear fit and provide estimators of best fit
    # this will utilise the two_phase_kinetic and two_phase_nonlin_fit functions and compare using AIC, R-squared and p-value, giving an estimation of best sampling rate
    pass

"""
One of the benefits of using passive samplers in the kinetic or non-linear (as opposed to equilibrium) phases is that they represent the concentration over the period of deployment
This allows us to back calculate a time weighted average concentration of an analyte for the period of deployment

"""

"""
General notes
Cs = Cw*Ksw(1-e^[-Ke*t])
where Cs is the concentration (μg g−1) of the analyte in the sorbent, 
Cw the TWA concentration (μg l−1) of the analyte in water, 
Ksw the POCIS-water partition constant (l g−1) and ke the elimination rate constant (d−1).


If the elimination rate ke is negligible compared to the uptake rate Ku (l g−1 d−1 or ml g−1 d−1), then the sampler acts as an infinite sink and we can redue to:
Cs = Cw*Ku*t

If we introduce the mass of the sorbent (Ms) we can change out Ku for a relationship with Rs (sampling rate)
Cs = (Cw*Rs*t)/Ms

Thus we can solve for Cw which is the time weighted average by rearranging to:
Cw = (Cs*Ms)/(Rs*t)
"""

def TWA(sampling_rate, time, sorbent_mass, conc_sorbent):
    return (conc_sorbent * sorbent_mass)/(sampling_rate * time)
    

def Ksw():
    pass

"""
Semi-permeable membrane devices (SPMD)
SPMDs are polyethylene membranes containing triolein (lipid), they are most useful for organic compounds with a LogP > 3, this section concerns their surface water uses
Examples of compounds SPMDs are typically used for are: Polycyclic Aromatic Hydrocarbons (PAHs), Polychlorinated Biphenyls (PCBs), Polybrominated Diphenyl Ethers (PBDEs), Organochlorine Pesticides, Fragrances, Dioxins and Furans
"""

# 7.1
def SPMD_flux(ki,Ci):
    return ki * Ci

# 7.4
def SPMD_mass_transfer_coef(WBL_mtc, biofilm_mtc, LDPE_mtc, biofilm_water_partition, membrane_water_partition):
    return (1/WBL_mtc) + (1/(biofilm_mtc*biofilm_water_partition) + (1/(LDPE_mtc*membrane_water_partition)))
    
# 7.8
def SPMD_sampling_rate(mass_transfer_resistance, surface_area):
    return mass_transfer_resistance * surface_area

def SPMD_equilibrium(sampler_conc, Ksw):
    return sampler_conc / Ksw