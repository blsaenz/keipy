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
import math, os
import numpy as np
#import h5py
#from matplotlib.dates import date2num,num2date
from numba import jit
from netCDF4 import date2num # use these with date2num(dt,'days since %i-01-01'%year) to convert!
#from numpy import asfortranarray,ascontiguousarray
#import pandas as pd
import xarray as xr
import netCDF4

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

ice_dim = 'nni'
snow_dim = 'nns'
nz_dim = 'zm'
t_dim = 'f_time'
flx_dim = 'nflx'

ocn_output_vars = { # single dimension output vars
               'time': {'units':'day from time_start',
                        'dim':None},
               'hmx': {'units':'m',
                       'dim':None},
               'zml': {'units':'m',
                       'dim':None},
               'atm_flux_to_ocn_surface': {'units': 'W/m^2',
                       'dim': None},
                # nz-dimension vars
               'wU': {'units':'m',
                       'dim':nz_dim},
               'wV': {'units':'m',
                       'dim':nz_dim},
               'wW': {'units':'m',
                       'dim':nz_dim},
               'wT': {'units':'m',
                       'dim':nz_dim},
               'wS': {'units':'m',
                       'dim':nz_dim},
               'wB': {'units':'m',
                       'dim':nz_dim},
               'Tprev': {'units':'C',
                       'dim':nz_dim},
               'Sprev': {'units':'psu',
                       'dim':nz_dim},
               'km': {'units':'',
                       'dim':nz_dim},
               'ks': {'units':'',
                       'dim':nz_dim},
               'kt': {'units':'',
                       'dim':nz_dim},
               'ghat': {'units':'',
                       'dim':nz_dim},
              }

eco_output_vars = {
    'tot_prod': {'units': '',
                 'dim': nz_dim},
    'sp_Fe_lim': {'units': '',
                  'dim': nz_dim},
    'sp_N_lim': {'units': '',
                 'dim': nz_dim},
    'sp_P_lim': {'units': '',
                 'dim': nz_dim},
    'sp_light_lim': {'units': '',
                     'dim': nz_dim},
    'diat_Fe_lim': {'units': '',
                    'dim': nz_dim},
    'diat_N_lim': {'units': '',
                   'dim': nz_dim},
    'diat_P_lim': {'units': '',
                   'dim': nz_dim},
    'diat_Si_lim': {'units': '',
                    'dim': nz_dim},
    'diat_light_lim': {'units': '',
                       'dim': nz_dim},
    'graze_sp': {'units': '',
                 'dim': nz_dim},
    'graze_diat': {'units': '',
                   'dim': nz_dim},
    'graze_tot': {'units': '',
                  'dim': nz_dim},
}

ice_output_vars = {
    'hi': {'units': 'm',
            'dim': None},
    'hs': {'units': 'm',
            'dim': None},
    'ni': {'units': 'valid levels',
           'dim': 'int'},
    'ns': {'units': 'valid levels',
           'dim': 'int'},
    'fice': {'units': 'fraction',
            'dim': None},
    'dzi': {'units': 'm',
           'dim': ice_dim},
    'Ti': {'units': 'C',
           'dim': ice_dim},
    'Si': {'units': 'psu',
           'dim': ice_dim},
    'dzs': {'units': 'm',
           'dim': snow_dim},
    'Ts': {'units': 'C',
           'dim': snow_dim},
    'atm_flux_to_ice_surface': {'units': 'W/m^2',
           'dim': None},
    'ice_ocean_bottom_flux': {'units': 'W/m^2',
           'dim': None},
    'ice_ocean_bottom_flux_potential': {'units': 'W/m^2',
           'dim': None},
    'total_ice_melt': {'units': 'J/m^2',
           'dim': None},
    'total_ice_freeze': {'units': 'J/m^2',
           'dim': None},
    'frazil_ice_volume': {'units': 'm^3/m^2',
           'dim': None},
    'congelation_ice_volume': {'units': 'm^3/m^2',
           'dim': None},
    'snow_ice_volume': {'units': 'm^3/m^2',
           'dim': None},
    'snow_precip_mass': {'units': 'kg/m^2',
           'dim': None},

}

output_var_data = {**ocn_output_vars, **eco_output_vars, **ice_output_vars}


output_var_data_meta = {
    'time':{'units':'days','long_name':'decimal days since start'},
    'hmx':{'units':'m','long_name':'KPP mixing depth'},
    'zml':{'units':'m','long_name':' gradient calc mixed layer thickness'},
    'atm_flux_to_ocn_surface':{'units':'W m-2','long_name':'energy flux from atm to ocean surface'},
    'wU':{'units':'m2 s-2','long_name':''},
    'wV':{'units':'m2 s-2','long_name':''},
    'wW':{'units':'m s-1','long_name':''},
    'wT':{'units':'celsius m s-1','long_name':''},
    'wS':{'units':'psu m s-1','long_name':''},
    'wB':{'units':'m s-1','long_name':''},
    'Tprev':{'units':'C','long_name':'water temperature, previous time step'},
    'Sprev':{'units':'psu','long_name':'salinty, previous time step'},
    'km':{'units':'m2 s-2','long_name':'momentum diffusivity coefficient'},
    'ks':{'units':'m2 s-1','long_name':'scalar diffusivity coefficient'},
    'kt':{'units':'m2 s-1','long_name':'temperature diffusivity coefficient'},
    'ghat':{'units':'','long_name':'gradient ghat'},
    'tot_prod':{'units':'mgC m-3','long_name':'total production'},
    'sp_Fe_lim':{'units':'fractional','long_name':'small phytoplankton Fe limtation term'},
    'sp_N_lim':{'units':'fractional','long_name':'small phytoplankton N limtation term'},
    'sp_P_lim':{'units':'fractional','long_name':'small phytoplankton P limtation term'},
    'sp_light_lim':{'units':'fractional','long_name':'small phytoplankton light limtation term'},
    'diat_Fe_lim':{'units':'fractional','long_name':'diatom phytoplankton Fe limtation term'},
    'diat_N_lim':{'units':'fractional','long_name':'diatom phytoplankton N limtation term'},
    'diat_P_lim':{'units':'fractional','long_name':'diatom phytoplankton P limtation term'},
    'diat_Si_lim':{'units':'fractional','long_name':'diatom phytoplankton Si limtation term'},
    'diat_light_lim':{'units':'fractional','long_name':'diatom phytoplankton light limtation term'},
    'graze_sp':{'units':'','long_name':'grazing of small phytos'},
    'graze_diat':{'units':'','long_name':'grazing of diatoms'},
    'graze_tot':{'units':'','long_name':'total grazing'},
    'hi':{'units':'m','long_name':'sea-ice thickness'},
    'hs':{'units':'m','long_name':'snow over sea-ice thickness'},
    'ni':{'units':'#','long_name':'number active ice layers'},
    'ns':{'units':'#','long_name':'number active snow laters'},
    'fice':{'units':'fractional','long_name':'fractional ice coverage'},
    'dzi':{'units':'m','long_name':'ice layer thicknesses'},
    'Ti':{'units':'C','long_name':'ice layer temperatures'},
    'Si':{'units':'psu','long_name':'ice layer salinities'},
    'dzs':{'units':'m','long_name':'snow layer thicknesses'},
    'Ts':{'units':'C','long_name':'snow/ice surface temperature'},
    'atm_flux_to_ice_surface':{'units':'W m-2','long_name':'energy flux from atmosphere to sea-ice surface'},
    'ice_ocean_bottom_flux':{'units':'W m-2','long_name':'energy flux to the sea-ice bottom from the PBL'},
    'ice_ocean_bottom_flux_potential':{'units':'W m-2','long_name':'ocean heat flux to ice potential'},
    'total_ice_melt':{'units':'J m-2','long_name':'total ice melted'},
    'total_ice_freeze':{'units':'J m-2','long_name':'total ice frozen'},
    'frazil_ice_volume':{'units':'m3 m-2','long_name':'total frazil ice production volume'},
    'congelation_ice_volume':{'units':'m3 m-2','long_name':'total congelation ice production volume'},
    'snow_ice_volume':{'units':'m3 m-2','long_name':'total snow ice production volume'},
    'snow_precip_mass':{'units':'kg m-2','long_name':'total snow fall over sea ice'},
    'fatm':{'units':'','long_name':''},
    'fao':{'units':'','long_name':''},
    'fai':{'units':'','long_name':''},
    'fio':{'units':'','long_name':''},
    'focn':{'units':'','long_name':''},
    'T':{'units':'C','long_name':'ocean later temperature'},
    'S':{'units':'psu','long_name':'ocean layer salinity'},
    'PO4':{'units':'mmol PO4 m-3','long_name':'ocean layer PO4 conentration'},
    'NO3':{'units':'mmol NO3 m-3','long_name':'ocean layer NO3 conentration'},
    'SiO3':{'units':'mmol SiO3 m-3','long_name':'ocean layer SiO3 conentration'},
    'NH4':{'units':'mmol NH4 m-3','long_name':'ocean layer NH4 conentration'},
    'Fe':{'units':'nmol Fe m-3','long_name':'ocean layer Fe conentration'},
    'O2':{'units':'nmol cm-3','long_name':'ocean layer O2 conentration'},
    'DIC':{'units':'mmol m-3','long_name':'ocean layer DIC conentration'},
    'ALK':{'units':'','long_name':'ocean layer Alkalinity'},
    'DOC':{'units':'mmol m-3','long_name':'ocean layer DOC conentration'},
    'spC':{'units':'mmol m-3','long_name':'small phytoplankton carbon'},
    'spChl':{'units':'mgChl m-3','long_name':'small phytoplankton Chloropyll'},
    'spCaCO3':{'units':'mmol CaCO3 m-3','long_name':'small phytoplankton CaCO3'},
    'diatC':{'units':'','long_name':'diatom phytoplankton carbon'},
    'diatChl':{'units':'mgChl m-3','long_name':'diatom phytoplankton Chloropyll'},
    'zooC':{'units':'','long_name':'heterotrophic zooplankton phytoplankton carbon'},
    'spFe':{'units':'nmol Fe m-3','long_name':'small phytoplankton Fe'},
    'diatSi':{'units':'','long_name':'diatom phytoplankton Si'},
    'diatFe':{'units':'nmol Fe m-3','long_name':'diatom phytoplankton Fe'},
    'diazC':{'units':'','long_name':'diazotroph phytoplankton carbon'},
    'diazChl':{'units':'mgChl m-3','long_name':'diazotroph phytoplankton Chloropyll'},
    'diazFe':{'units':'nmol Fe m-3','long_name':'diazotroph phytoplankton Fe'},
    'DON':{'units':'mmol N m-3','long_name':'dissolved organic N'},
    'DOFe':{'units':'nmol Fe m-3','long_name':'dissolved organic Fe'},
    'DOP':{'units':'mmol P m-3','long_name':'dissolved organic P'},
    'hour':{'units':'fractional day','long_name':'this variable is poorly named'}, # TODO: rename me
    'date':{'units':'days','long_name':'days since simulation start'},
    'tau_x':{'units':'m s-2','long_name':'x-direction 10m wind speed'},  # should rename var; tau indicates stress...
    'tau_y':{'units':'m s-2','long_name':'y-direction 10m wind speed'},  # should rename var; tau indicates stress...
    'qswins':{'units':'W m-2','long_name':'surface downward shortwave irradiance'},
    'qlwdwn':{'units':'W m-2','long_name':'surface downward longwave irradiance'},
    'tz':{'units':'C','long_name':'2m [surface] air temperature'},
    'qz':{'units':'kg m-3','long_name':'2m [surface] humidity'},
    'prain':{'units':'kg m-2 s-1','long_name':'rain rate'},
    'psnow':{'units':'kg m-2 s-1','long_name':'snow rate'},
    'msl':{'units':'mbar','long_name':'mean sea level pressure'},
    'h':{'units':'kg/kg','long_name':'specific humdity'},
    'dustf':{'units':'g Fe m-2 s-2','long_name':'dust flux atmipshere to ocean surface - re-check units'},
    'divu':{'units':'fractional','long_name':'sea-ice coverage divergence'},
    'ic':{'units':'fractional','long_name':'fractional sea-ice coverage'},
    'ain':{'units':'fractional','long_name':'sea-ice coverage influx'},
    'aout':{'units':'fractional','long_name':'sea-ice coverage outflux'},
    'zm':{'units':'m','long_name':'layer midpoint depth [negative downward]'},
    'nni':{'units':'level','long_name':'sea-ice layers'},
    'nns':{'units':'level','long_name':'snow layers'},
    'nflx':{'units':'#','long_name':'KEI flux structure count'},
 }


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
            ds[v] = ('zm'),ds_in[v][0:nz].data
    for v in forcing_idx.keys():
        if v in ds_in.data_vars:
            ds[v] = ('f_time'),ds_in[v].data
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


def compile_ERA5_met():
    #from metpy.calc import relative_humidity_from_dewpoint, relative_humidity_from_specific_humidity
    #from metpy.calc import specific_humidity_from_dewpoint, mixing_ratio_from_relative_humidity, specific_humidity_from_mixing_ratio
    #from metpy.units import units

    # lat cell 26 (0-based indexing)
    # lon cell 1 (0-based indexing)
    xi = 1
    yi = 26

    ds_irad = xr.open_dataset(r'~/Downloads/ERA5_WAP_stuff_I_forgot.nc',engine='netcdf4')

    ecmwf_vars = ['d2m','msl','msr','t2m','tp','u10','v10']
    kei_vars = ['qz','msl','psnow','tz','prain','tau_x','tau_y']

    data = {k:None for k in ecmwf_vars}
    for y in range(2007,2012):
        dsy = xr.open_dataset(os.path.join(r'/Users/blsaenz/Downloads','ERA5_WAP_%i.nc'%y),engine='netcdf4')
        for k in ecmwf_vars:
            if data[k] is None:
                data[k] = dsy[k][:,yi,xi].values
            else:
                data[k] = np.append(data[k],dsy[k][:,yi,xi].values)

    data['msr'][data['msr']<0.0] = 0.0
    data['tp'][data['tp']<0.0] = 0.0

    kdata = {}
    kdata['qlwdwn'] = ds_irad['strd'][:,yi,xi].values/3600.0
    kdata['qswins'] = ds_irad['ssrd'][:,yi,xi].values/3600.0
    kdata['qlwdwn'][kdata['qlwdwn']<0.0] = 0.0
    kdata['qswins'][kdata['qswins']<0.0] = 0.0

    for i,k in enumerate(ecmwf_vars):
        kk = kei_vars[i]
        if k == 'tp':
            # prain == tp - snow
            kdata[kk] = np.maximum(0.0,data['tp']/3600.0*1000.0 - data['msr'])
        elif k == 'd2m':
            # convert dewpoint to specififc humidity, I think
            #kdata[kk] = specific_humidity_from_dewpoint(data['msl'] * units.Pa, data[k] * units.degC).to('kg/kg')

            Tice = 250.16
            alpha = np.ones(len(data[k]))
            t2m = data['t2m']
            for i in range(len(data[k])):
                if t2m[i] <= Tice:
                    alpha[i] = 0.0
                elif t2m[i] < 273.16:
                    alpha[i] = ((t2m[i] - Tice) / (273.16 - Tice)) ** 2

            esat_w = 611.2 * np.exp(17.502 * ((data[k] - 273.16) / (data[k] - 32.19))) # over water
            esat_i = 611.2 * np.exp(22.587 * ((data[k] - 273.16) / (data[k] + 20.7))) # over ice
            Rair = 287.058 # J / kg / K
            Rvap = 461.495 # J / kg / K
            r_o_r = Rair / Rvap
            h_w = (r_o_r * esat_w) / (data['msl'] - (1 - r_o_r) * esat_w)
            h_i = (r_o_r * esat_i) / (data['msl'] - (1 - r_o_r) * esat_i)
            kdata['h'] = alpha * h_w + (1.0 - alpha) * h_i
            rhoair = data['msl'] / (data['t2m'] * (Rair - kdata['h'] * Rair + kdata['h'] * Rvap)) # kg / m ^ 3
            kdata[kk] = kdata['h'] * rhoair # kg / m ^ 3

        elif k == 't2m':
            kdata[kk] = data[k] - 273.16
        else:
            kdata[kk] = data[k]

    # add to current dataset and write out netcdf3
    for year in range(7,12):
        if year != 9:
            print('Saving over forcing:',year+2000)
            rg = netCDF4.Dataset('/Users/blsaenz/KEI_run/DATA/kf_300_%02i_ERA5_div.nc'%year,'a')
            istart = ((year+2000)-2007)*8760
            if year > 10:
                istart += 24
            flen = len(kdata['tz'][istart:])
            for fv in kei_vars + ['h']:
                rg.variables[fv][0:flen] = kdata[fv][istart:]
                rg.sync()
            rg.close()

    dude=1


    # specific_humidity_from_dewpoint(988 * units.hPa, 15 * units.degC).to('g/kg')




def write_get_set_f90_code(params):
  '''there could be a lot of get/set hardcoded variables that we want to change dynamically. We could list them
  here and have python generate the subroutines automagically.'''
  pass

#compile_ERA5_met()
#exit()
