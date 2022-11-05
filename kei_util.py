'''
!     ==================================================================
!     KPP-Ecosystem-Ice (KEI, pronounced 'key') Model
!     ==================================================================
!
!     Version History (please add a version modification line below
!     if code is edited)
!     ------------------------------------------------------------------
!     Version: 0.9 (2022-01-05, Benjamin Saenz, blsaenz@gmail.com)
!
!     This model derives from Large et al. [1994], Doney et al. [1996](KPP mixing),
!     Ukita and Martinson [2001](Mixed layer - ice interactions), Saenz and Arrigo
!     [2012 and 2014] (SIESTA sea ice model), Hunke and Lipscomb [2008] (CICE v4),
!     Moore et al. [2002,2004] (CESM ecosystem model).
!
!     Python version: upcoming

Notes on fortran building (gfortran):
gfortran -c -m64 -02 mag_kinds_mod.f90
gfortran -c -m64 -02 mag_parameters_mod.f90
f2py -c -m64 -02 -fopenmp -fdec-math mag_kinds_mod.f90 mag_parameters_mod.f90 mag_calc.f90 -m mag_calc


'''
#import os,sys,shutil,csv,pickle,math,time,datetime
#from calendar import isleap
import math
import numpy as np
#import h5py
#from matplotlib.dates import date2num,num2date
from numba import jit
from netCDF4 import date2num # use these with date2num(dt,'days since %i-01-01'%year) to convert!
#from numpy import asfortranarray,ascontiguousarray
#import pandas as pd
import xarray as xr

# some fixed information regarding kei tracers and data
ecosys_tracers = {'PO4':     {'units':'',  'idx':0},
                  'NO3':     {'units':'',  'idx':1},
                  'SiO3':    {'units':'',  'idx':2},
                  'NH4':     {'units':'',  'idx':3},
                  'Fe':      {'units':'',  'idx':4},
                  'O2':      {'units':'',  'idx':5},
                  'DIC':     {'units':'',  'idx':6},
                  'ALK':     {'units':'',  'idx':7},
                  'DOC':     {'units':'',  'idx':8},
                  'spC':     {'units':'',  'idx':9},
                  'spChl':   {'units':'',  'idx':10},
                  'spCaCO3': {'units':'',  'idx':11},
                  'diatC':   {'units':'',  'idx':12},
                  'diatChl': {'units':'',  'idx':13},
                  'zooC':    {'units':'',  'idx':14},
                  'spFe':    {'units':'',  'idx':15},
                  'diatSi':  {'units':'',  'idx':16},
                  'diatFe':  {'units':'',  'idx':17},
                  'diazC':   {'units':'',  'idx':18},
                  'diazChl': {'units':'',  'idx':19},
                  'diazFe':  {'units':'',  'idx':20},
                  'DON':     {'units':'',  'idx':21},
                  'DOFe':    {'units':'',  'idx':22},
                  'DOP':     {'units':'',  'idx':23},
                  }

ma_tracers = {'maC':   {'units':'',  'idx':24},
              'maChl': {'units':'',  'idx':25},
              'maN':   {'units':'',  'idx':26},
              'maP':   {'units':'',  'idx':27},
              }

# the order of forcing vars in a forcing array - the must be the same as in fortran code
forcing_idx =  {'date': 0,
                'tau_x': 1,
                'tau_y': 2,
                'qswins': 3,
                'qlwdwn': 4,
                'tz': 5,
                'qz': 6,
                'prain': 7,
                'psnow': 8,
                'msl': 9,
                'h': 10,
                'dustf': 11,
                'divu': 12,
                'ic': 13,
                'ain': 14,
                'aout': 15
}
grid_vars = ['dm','hm','zm','f_time'] # midpoint depth of cells, at least needed for xarray
init_vars_ocn = ['t','s','u','v']
init_vars_eco = list(ecosys_tracers.keys())
forcing_vars = list(forcing_idx.keys())



@jit(nopython=True)
def TfrzC(S, Db):
    """Freezing point of water in degrees C at salinity S in PSU and pressure Db in decibars
     -- older relationshp:
    TfrzC = (-0.0575 +1.710523e-3 *sqrt(S) -2.154996e-4 *S) *S - 7.53e-4 *Db
    --- newer relationship below is compatible with sea ice model integrations
    """
    return -0.054 * S - 7.53e-4 * Db

@jit(nopython=True)
def CPSW(S, T1, P0):
    """
    # UNITS:
    #       PRESSURE        P0       DECIBARS
    #       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
    #       SALINITY        S        (IPSS-78)
    #       SPECIFIC HEAT   CPSW     J/(KG DEG C)
    # ***
    # REF: MILLERO ET AL,1973,JGR,78,4499-4507
    #       MILLERO ET AL, UNESCO REPORT NO. 38 1981 PP. 99-188.
    # PRESSURE VARIATION FROM LEAST SQUARES POLYNOMIAL
    # DEVELOPED BY FOFONOFF 1980.
    # ***
    # CHECK VALUE: CPSW = 3849.500 J/(KG DEG. C) FOR S = 40 (IPSS-78),
    # T = 40 DEG C, P0= 10000 DECIBARS
    """

    #   check that temperature is above -2
    T = T1
    if T < -2.:
        T = -2.

    #   SCALE PRESSURE TO BARS
    P = P0 / 10.

    # SQRT SALINITY FOR FRACTIONAL TERMS
    SR = math.sqrt(math.abs(S))
    # SPECIFIC HEAT CP0 FOR P=0 (MILLERO ET AL ,UNESCO 1981)
    A = (-1.38385E-3 * T + 0.1072763) * T - 7.643575
    B = (5.148E-5 * T - 4.07718E-3) * T + 0.1770383
    C = (((2.093236E-5 * T - 2.654387E-3) * T + 0.1412855) * T -3.720283) * T + 4217.4
    CP0 = (B * SR + A) * S + C

    # CP1 PRESSURE AND TEMPERATURE TERMS FOR S = 0
    A = (((1.7168E-8 * T + 2.0357E-6) * T - 3.13885E-4) * T + 1.45747E-2) * T -0.49592
    B = (((2.2956E-11 * T - 4.0027E-9) * T + 2.87533E-7) * T - 1.08645E-5) * T +2.4931E-4
    C = ((6.136E-13 * T - 6.5637E-11) * T + 2.6380E-9) * T - 5.422E-8
    CP1 = ((C * P + B) * P + A) * P

    # CP2 PRESSURE AND TEMPERATURE TERMS FOR S > 0
    A = (((-2.9179E-10 * T + 2.5941E-8) * T + 9.802E-7) * T - 1.28315E-4) * T +4.9247E-3
    B = (3.122E-8 * T - 1.517E-6) * T - 1.2331E-4
    A = (A + B * SR) * S
    B = ((1.8448E-11 * T - 2.3905E-9) * T + 1.17054E-7) * T - 2.9558E-6
    B = (B + 9.971E-8 * SR) * S
    C = (3.513E-13 * T - 1.7682E-11) * T + 5.540E-10
    C = (C - 1.4300E-12 * T * SR) * S
    CP2 = ((C * P + B) * P + A) * P

    # SPECIFIC HEAT RETURN
    return CP0 + CP1 + CP2

def write_grid_f90_code(params):
  '''write out params fortran 90 file before compilation, which defines the model domain, grid, etc.'''
  bounds_params = {
    'DZ':400

  }

def reindex_forcing(ds_in,f_time,zm):
    ds = xr.Dataset()
    nz = len(zm)
    ds['f_time'] = ('f_time'),f_time
    ds['zm'] = ('zm'),zm
    ds['zm'].attrs['units'] = 'm'
    for v in ['dm','hm'] + init_vars_ocn + init_vars_eco:
        if v in ds_in.data_vars:
            ds[v] = ('zm'),ds_in[v][0:nz]
    for v in forcing_idx.keys():
        if v in ds_in.data_vars:
            ds[v] = ('f_time'),ds_in[v]
    return ds

def salinity_correction_rate(total_precip_forcing,init_s_profile,hm,dt_seconds,simulation_days):
    #     !  FORTRAN code:
    #     ! find salinity addition from total freshwater input, hourly time step:
    #     freshwater = 3600 * &
    #         (sum(kforce%f_data(:,prain_f_ind)) &
    #          + sum(kforce%f_data(:,psnow_f_ind)))
    #     ! find mean initial salinity
    #     mean_sal = 0.
    #     weight = 0.
    #     do i=1,nzp1
    #         mean_sal = mean_sal + (X(i,2)+Sref) * hm(i)
    #         weight = weight + hm(i)
    #     enddo
    #     mean_sal = mean_sal/weight
    #     ! total salt deficit (kg) ~= fresh_mass (kg) * mean ptt (kg/kg/1000)
    #     ! salinity correction rate = total salt deficit (kg) / total days (day) / (s/day)
    #     sal_correction_rate = &
    #         freshwater * mean_sal / 1000 &
    #         / maxval(kforce%f_data(:,date_f_ind)) / 86400.  ! (kg/s)
    #
    #     ! too high b/c of evap, potential rain problem in ice model
    #     sal_correction_rate = sal_correction_rate / 2.

    # total freshwater input
    freshwater = np.sum(total_precip_forcing)*dt_seconds

    # mean initial salinity
    mean_sal = np.sum(init_s_profile*hm)/np.sum(hm)

    # total salinity correction rate, with no consideration of evap, potential ignorance of rain in ice model
    sal_correction_rate = freshwater * mean_sal / 1000. / simulation_days / 86400.  # (kg/s)

    # account for evap, ice losses by dividing by 2.  wow.
    sal_correction_rate = sal_correction_rate / 2.

    return sal_correction_rate


def write_get_set_f90_code(params):
  '''there could be a lot of get/set hardcoded variables that we want to change dynamically. We could list them
  here and have python generate the subroutines automagically.'''
  pass
