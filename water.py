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
    TODO
    return ksw * sorbent_mass


# estimate sampling rate for the curvilinear phase for a passive sampler with limited WBL intereference
# this uses the nonlin function and a curve fitting method from scipy.optimize
def two_phase_nonlin_fit(df, time_column, compound_column, time_unit = 'day', water_unit = 'mL', plot = False): 
    # use the scipy optimise curve fit function to fit the data to the nonlin equation, estimating unknown parameters
    params, covs = curve_fit(nonlin, df.loc[:, time_column], df.loc[:, compound_column])
    
    # pulls out the estimated sorbent water partitioning coefficient (ksw) and the estimated sampling rate
    ksw, sampling_rate = params[0], params[1]
    
    # if plot is called in the args this will create a non-linear uptake graph
    if plot == True:
        print(ksw, str(sampling_rate + water_unit + "/" + time_unit))
        # plot_range creates a range from the lowest to highest time point, allowing the curve fit to plot smoothly by interpolating unknowns between data points
        plot_range = np.arange(min(time), max(time))
        # substantiate the graph space
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
    TODO
    # there are two possibilities here, the second is the red book, the first is from Sarit's prepub
    v1 = (math.log(2)*sorbent_mass*ksw)/sampling_rate
    v2 = (math.log(2)/(sampling_rate/(ksw*sorbent_mass))) # Vs undefined by Alvarez, great fucking work there. Rory confirmed sorbent mass and there are multiple ways to calculate


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
        pass



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
    # this is the base calculation for determining the time weighted average based on a calibrated passive sampler
    return (conc_sorbent * sorbent_mass)/(sampling_rate * time)
    

def Ksw(conc_sorbent, conc_water):
    # base calculation for determining sorbent/water coefficient
    return (conc_sorbent/conc_water)

"""
Semi-permeable membrane devices (SPMD)
SPMDs are polyethylene membranes containing triolein (lipid), they are most useful for organic compounds with a LogP > 3, this section concerns their surface water uses
Examples of compounds SPMDs are typically used for are: Polycyclic Aromatic Hydrocarbons (PAHs), Polychlorinated Biphenyls (PCBs), Polybrominated Diphenyl Ethers (PBDEs), Organochlorine Pesticides, Fragrances, Dioxins and Furans
"""

# 7.1
def SPMD_flux(ki,Ci):
    TODO
    return ki * Ci

# 7.2
#def 

def SPMD_ksw(Vm, Kmw, KLW, VL):

    """
    Determination of Uptake Kinetics (Sampling Rates) by Lipid-Containing Semipermeable Membrane Devices (SPMDs) for Polycyclic Aromatic Hydrocarbons (PAHs) in Water
    James N. Huckins, Jimmie D. Petty, Carl E. Orazio, Jon A. Lebo, Randal C. Clark, Virginia L. Gibson, William R. Gala, and Kathy R. Echols
    Environmental Science & Technology 1999 33 (21), 3918-3923
    DOI: 10.1021/es990440u  """

    # where Vm = volume of the membrane, Kmw = volume-averaged partition coefficient for the membrane, KLW = volume-averaged partition coefficient for the triolein, VL = volume of lipid
    return (Vm*Kmw)+(VL*KLW)/(Vm+VL)

# 7.8
def SPMD_sampling_rate(mass_transfer_resistance, surface_area):
    # calculate the sampling rate for an analyte with the SPMD passive sampler
    return mass_transfer_resistance * surface_area

def SPMD_equilibrium(sampler_conc, Ksw):
    # Calculate the concentration in the water at the time of retrieval assuming SPMD is in equilibrium with surrounding environment for the relevant analyte
    return sampler_conc / Ksw