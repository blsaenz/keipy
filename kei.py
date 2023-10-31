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


'''
import os,sys,shutil,csv,pickle,math,time,datetime
from calendar import isleap
import numpy as np
from numpy import asfortranarray,ascontiguousarray
#from numba import jit
#from netCDF4 import date2num # use these with date2num(dt,'days since %i-01-01'%year) to convert!
import xarray as xr

# import local utils and fortran modules
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import kei_util as util
from kei_util import ecosys_tracers,ma_tracers,forcing_idx,init_vars_ocn,init_vars_eco
from kei_util import ocn_output_vars,eco_output_vars,ice_output_vars,output_var_data
from kei_util import output_var_data_meta

#try:
from f90 import kei
#except:
#    print('kei not imported/available!')
#    sys.exit()


nf = len(forcing_idx)

grid_vars = ['dm','hm','zm','f_time'] # midpoint depth of cells, at least needed for xarray


# fantastic thing to have around ...
month_doy = [1,32,60,91,121,152,182,213,244,274,305,335]


class kei_parameters(object):
    #default params
    p = {
      'dt'  :           3600.0,  # [seconds] 1 Hour default time step
      'lice':           1,     # ice model enabled
      'leco':           1,     # ecosystem model enabled
      #'lma':           0,     # seaweed model enabled
    }

    def __init__(self, params={}):
        self.update(params)

    def defaults(self):
        return self.p

    def update(self,params):
        self.p.update(params)
        self.build()

    def build(self):
        '''Perform any supplementary calculations for derived parameters'''
        pass


def kei_forcing(nc_file = None, f_dict = {}, start_date=None, freq=None, legacy_nc=False):
    ''' Reads and updates forcing data into XArray dataset, which can be fed to KEI model simulation for easy
    interpolation or whatever

    All variables fed into this class should be 1D numpy arrays, except those that come from reading a netcdf file.

    f_time must be a pandas/xarray time series, I think

    Currenty, all forcing vars must have similar timing
    '''

    if nc_file is not None:
        ds_in = xr.open_dataset(nc_file)
        if legacy_nc:
            # we are going to make a new dataset with the correct time-index dimensions
            if start_date is None or freq is None:
                raise ValueError('f_time is not provided, and start_date and/or time_delta are not provided; need one of them')
            else:
                f_time_len = len(ds_in['tz'][...])
                f_time = xr.cftime_range(start=start_date,periods=f_time_len,freq=freq)
                zm = ds_in['zm'][0:-1].data
            ds = util.reindex_forcing(ds_in,f_time,zm)
            ds['msl'][...] = ds['msl'][...] * 0.01  # many old forcing netCDF files are in Pa, meed mbar
        else:
            ds = ds_in
    else:
        ds = xr.Dataset()

    for k,v in f_dict.items():
        if k in grid_vars:
            ds[k] = (k),v
        elif k in forcing_idx.keys() + init_vars_ocn + init_vars_eco:
            ds[k] = ('f_time'),v
        elif k in forcing_vars_init:
            ds[k] = ('zm'),v

    return ds


class kei_output(object):

    flx_vars = ['fatm','fao','fai','fio','focn']

    def __init__(self,out_vars,f_time,zm,nflx,nni,nns):

        self.out_vars = out_vars

        self.out_ds = xr.Dataset()
        self.out_ds['f_time'] = ('f_time'),f_time
        self.out_ds['zm'] = ('zm'),zm
        self.out_ds['zm'].attrs['units'] = 'm'
        self.out_ds['zm'].attrs['long_name'] = 'layer midpoint depth'
        self.out_ds['nni'] = ('nni'),np.arange(nni,dtype=np.int)
        self.out_ds['nns'] = ('nns'),np.arange(nns,dtype=np.int)
        self.out_ds['nflx'] = ('nflx'),np.arange(nflx,dtype=np.int)
        # add dimension var meta
        for v in ['zm','nni','nns','nflx']:
            self.out_ds[v].attrs['units'] = output_var_data_meta[v]['units']
            self.out_ds[v].attrs['long_name'] = output_var_data_meta[v]['long_name']
        for d in ['nni','nns']:
            self.out_ds[d].attrs['positive'] = 'down'
        self.out_ds['zm'].attrs['positive'] = 'up'

        len_f_time = len(f_time)
        len_zm = len(zm)



        # setup lists of the exposed numpy arrays for potential usage in numba-optimized
        # methods

        # self.vars_1D = {}
        # self.vars_int = {}
        # self.vars_2D = {}
        # self.vars_ice = {}
        # self.vars_snow = {}
        # for v in self.out_vars:
        #     if output_var_data[v]['dim'] is None:
        #         self.out_ds[v] = ('f_time'),np.full(len_f_time,np.nan,np.float32)
        #         self.vars_1D[v] = self.out_ds[v].__array__()  # is this the best way?
        #     elif output_vars[v]['dim'] == 'int':
        #         self.out_ds[v] = ('f_time'),np.full(len_f_time,0,np.int32)
        #         self.vars_int[v] = self.out_ds[v].__array__()  # is this the best way?
        #     elif output_vars[v]['dim'] == 'zm':
        #         self.out_ds[v] = ('zm','f_time'),np.full((len_zm,len_f_time),np.nan,np.float32)
        #         self.vars_2D[v] = self.out_ds[v].__array__()  # is this the best way?
        #     elif output_vars[v]['dim'] == 'nni':
        #         self.out_ds[v] = ('nni','f_time'),np.full((nni,len_f_time),np.nan,np.float32)
        #         self.vars_ice[v] = self.out_ds[v].__array__()  # is this the best way?
        #     elif output_vars[v]['dim'] == 'nns':
        #         self.out_ds[v] = ('nns','f_time'),np.full((nns,len_f_time),np.nan,np.float32)
        #         self.vars_snow[v] = self.out_ds[v].__array__()  # is this the best way?

        # classify output variables so we can make the right fortran calls to record them
        self.vars_1D = []
        self.vars_int = []
        self.vars_2D = []
        self.vars_ice = []
        self.vars_snow = []
        for v in self.out_vars:
            if output_var_data[v]['dim'] is None:
                self.out_ds[v] = ('f_time'),np.full(len_f_time,np.nan,np.float32)
                self.vars_1D.append(v)
            elif output_var_data[v]['dim'] == 'int':
                self.out_ds[v] = ('f_time'),np.full(len_f_time,0,np.int32)
                self.vars_int.append(v)
            else:
                if output_var_data[v]['dim'] == 'zm':
                    self.vars_2D.append(v)
                    self.out_ds[v] = ('zm', 'f_time'), np.full((len_zm, len_f_time),np.nan,np.float32)
                elif output_var_data[v]['dim'] == 'nni':
                    self.vars_ice.append(v)
                    self.out_ds[v] = ('nni', 'f_time'), np.full((nni, len_f_time),np.nan,np.float32)
                elif output_var_data[v]['dim'] == 'nns':
                    self.vars_snow.append(v)
                    self.out_ds[v] = ('nns', 'f_time'), np.full((nns, len_f_time), np.nan, np.float32)
            self.out_ds[v].attrs['units'] = output_var_data_meta[v]['units']
            self.out_ds[v].attrs['long_name'] = output_var_data_meta[v]['long_name']

    def store_tracers(self,Vsave,Tsave,Flxsave,leco):
        self.leco=leco

        print('Storing Fluxes...')
        for i,fv in enumerate(self.flx_vars):
            self.out_ds[fv] = ('nflx','f_time'),Flxsave[:,i,:]
            self.out_ds[fv].attrs['units'] = output_var_data_meta[fv]['units']
            self.out_ds[fv].attrs['long_name'] = output_var_data_meta[fv]['long_name']
        print('Storing Tracers...')
        self.out_ds['T'] = ('zm', 'f_time'), Tsave[:, 0, :]
        self.out_ds['T']['units'] = 'C'
        self.out_ds['S'] = ('zm', 'f_time'), Tsave[:, 1, :]
        self.out_ds['S']['units'] = 'psu'
        for i,fv in enumerate(['T','S']):
            self.out_ds[fv].attrs['units'] = output_var_data_meta[fv]['units']
            self.out_ds[fv].attrs['long_name'] = output_var_data_meta[fv]['long_name']

        if leco:
            for fv,fvd in ecosys_tracers.items():
                self.out_ds[fv] = ('zm','f_time'),Tsave[:,fvd['idx']+2,:]
                self.out_ds[fv].attrs['units'] = output_var_data_meta[fv]['units']
                self.out_ds[fv].attrs['long_name'] = output_var_data_meta[fv]['long_name']

    def store_step_outvars(self,kei,nt):

        # for v, arr in output.vars_1D.items():
        #     arr[nt] = kei.link.get_data_real(v)
        # for v, arr in output.vars_2D.items():
        #     arr[0:nz, nt] = kei.link.get_nz_data(v)
        # for v, arr in output.vars_ice.items():
        #     arr[0:nni, nt] = kei.link.get_data_ice(v)
        # for v, arr in output.vars_snow.items():
        #     arr[0:nns, nt] = kei.link.get_data_snow(v)
        # for v, arr in output.vars_int.items():
        #     arr[nt] = kei.link.get_data_int(v)

        for v in self.vars_1D:
            self.out_ds[v][nt] = kei.link.get_data_real(v)
        for v in self.vars_2D:
            self.out_ds[v][:,nt] = kei.link.get_nz_data(v)
        for v in self.vars_ice:
            self.out_ds[v][:,nt] = kei.link.get_ice_data(v)
        for v in self.vars_snow:
            self.out_ds[v][:,nt] = kei.link.get_snow_data(v)
        for v in self.vars_int:
            self.out_ds[v][nt] = kei.link.get_data_int(v)


    def write(self,out_filepath,Finterp,params):

        # create output directory


        # write KEI parameters, run parameters, any config we can think of


        # add compatibility vars for matlab plotting routines
        self.out_ds['hour'] = self.out_ds['f_time'].dt.hour / 24.0
        day = self.out_ds['f_time'].dt.day + self.out_ds['hour']
        #days must be increasing always, above 365
        day_add = 0
        nt = len(self.out_ds['hour'])
        seq_days = np.zeros(nt)
        seq_days[0] = day[0]
        for i in range(1,nt):
            if (day[i] - day[i-1]) < 0:
                day_add = int(seq_days[i-1])
            seq_days[i] = day[i] + day_add
        self.out_ds['day'] = seq_days

        # add forcing variables to output dataset
        for v in forcing_idx.keys():
            self.out_ds[v] = Finterp[v]

        # create encodings
        #compress_vars = list(self.vars_1D.keys()) + list(self.vars_2D.keys()) + list(forcing_idx.keys()) + \
        #                list(self.vars_ice.keys()) + list(self.vars_snow.keys()) + list(self.vars_flx.keys())
        compress_vars = self.out_vars + list(forcing_idx.keys()) + ['T','S'] + self.flx_vars
        if self.leco:
            compress_vars += list(ecosys_tracers.keys())
        encoding = {}
        for v in compress_vars:
            encoding[v] = {"dtype":np.float32,"zlib": True, "complevel": 4}

        # write
        self.out_ds.to_netcdf(out_filepath,mode='w',format='netcdf4',encoding=encoding)



class kei_simulation(object):

    # Compile time parameters defaults - these will be updated in fortran code for JIT compilation,
    # if requested
    f90_params = {
            'NZ':       400,  # number of vertical water layers
    }
    recompile = False

    def __init__(self,forcing_ds,t_start=None,t_end=None,f90_params=None,recompile=False):

        self.update_f90(f90_params,recompile)
        self.update_forcing(forcing_ds,t_start,t_end)


    def update_f90(self,f90_params=None,recompile=True):
        if f90_params is not None:
          self.f90_params.update(f90_params)
        self.recompile = recompile


    def update_forcing(self,forcing_ds,t_start=None,t_end=None):
        '''Save and interpolate forcing'''
        self.F0 = forcing_ds

        if t_start is None:
          self.t_start = forcing_ds.f_time[0]
        else:
          self.t_start = t_start
        if t_end is None:
          self.t_end = forcing_ds.f_time[-1]
        else:
          self.t_end = t_end


    def compute(self,params,output_path,out_vars=None,run_name='keipy'):

        # recompile fortran module if requested, then re-import kei
        # .....  to do

        # get instance of kei, perform parameter init, and query dimensions
        kei.link.kei_param_init()
        nvel = kei.kei_parameters.nvel  # total number of water velocities (2)
        nsclr = kei.kei_parameters.nsclr  # total number of tracers/scalers
        nsflxs = kei.kei_parameters.nsflxs  # total number of tracers/scalers
        nni = kei.kei_icecommon.nni  # total number of ice layers
        nns = kei.kei_icecommon.nns  # total number of snow layers
        nz = kei.kei_parameters.nz

        # prepare & interpolate forcing
        dt_str = '%iS'%params.p['dt']
        Finterp = self.F0[list(forcing_idx.keys())+['f_time']]
        Finterp = Finterp.sel(f_time=slice(self.t_start, self.t_end))
        Finterp = Finterp.resample({'f_time':dt_str}).interpolate()
        nt = Finterp.dims['f_time']

        # copy interpolated forcing into an array for fast import into kei
        Fcomp = np.zeros((nf,nt),order='F',dtype=np.float32)
        for k,idx in forcing_idx.items():
          Fcomp[idx,:] = Finterp[k][...]

        # write out link_test data
        #np.savetxt(r'/Users/blsaenz/Projects/git/keipy/test_data/kf_200_100_2000_savetxt.txt',Fcomp[:,0:1000])

        # update changeable parameters
        for k,v in params.p.items():
          if isinstance(v,int):
            kei.link.set_param_int(k,v)
          elif isinstance(v,float):
            kei.link.set_param_real(k,v)
          else:
            raise ValueError('Unknown kei parameters type:',k,type(v))
        kei.link.set_param_int('nend', nt)

        # init local storage
        Velocity = np.zeros((nz,nvel),order='F')
        Tracers = np.zeros((nz,nsclr),order='F')
        Fluxes = np.zeros((nsflxs,5),order='F')
        Vsave = np.full((nz,nvel,nt),np.nan,np.float32)
        Tsave = np.full((nz,nsclr,nt),np.nan,np.float32)
        Flxsave = np.full((nsflxs,5,nt),np.nan,np.float32)

        # init/copy fortran storage
        kei.link.set_grid(self.F0['dm'].values,
                          self.F0['hm'].values,
                          self.F0['zm'].values)
        Velocity[:,0] = self.F0['u'][0:nz]
        Velocity[:,1] = self.F0['v'][0:nz]
        Tracers[:,0] = self.F0['t'][0:nz]
        Tracers[:,1] = self.F0['s'][0:nz]
        for t,v in ecosys_tracers.items():
          Tracers[:,v['idx']+2] = self.F0[t][0:nz]
        # macroalgae tracers will go here as well, likely in separate array
        kei.link.set_tracers(Velocity,Tracers)

        # save init data for fortran_test
        # np.savetxt(r'/Users/blsaenz/Projects/git/keipy/test_data/dm_savetxt.txt',self.F0['dm'][...])
        # np.savetxt(r'/Users/blsaenz/Projects/git/keipy/test_data/hm_savetxt.txt',self.F0['hm'][...])
        # np.savetxt(r'/Users/blsaenz/Projects/git/keipy/test_data/zm_savetxt.txt',self.F0['zm'][...])
        #np.savetxt(r'/Users/blsaenz/Projects/git/keipy/test_data/X_savetxt.txt',np.transpose(Tracers))
        #np.savetxt(r'/Users/blsaenz/Projects/git/keipy/test_data/U_savetxt.txt',np.transpose(Velocity))

        # initialize
        kei.link.kei_compute_init()

        # init output
        now_str = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
        self.output_path = os.path.join(output_path,run_name + '_' + now_str)
        self.out_vars = list(output_var_data.keys()) if out_vars is None else out_vars
        if not os.path.exists(self.output_path):
            os.mkdir(self.output_path)
        else:
            raise ValueError('KEI output path exists!',self.output_path)
        output = kei_output(self.out_vars,Finterp['f_time'].data,self.F0['zm'].data,nsflxs,
                            nni=nni,nns=nns)

        # copy code
        shutil.copytree(os.path.dirname(__file__),os.path.join(self.output_path,'code'))

        # write option to txt and pickle
        params_csv = os.path.join(self.output_path,run_name+'_keipy_params.csv')
        with open(params_csv,"w",newline='', encoding='utf-8') as csvfile:
            w = csv.writer(csvfile)
            for key, val in params.p.items():
                w.writerow([key, val])

        # main loop
        for nt,time in enumerate(Finterp['f_time']):

            # load atm data to fortran
            kei.link.set_forcing(Fcomp[:,nt])

            # calc step
            kei.link.kei_compute_step(nt)

            # store step data
            #Fluxes = kei.link.get_fluxes()
            Flxsave[...,nt] = kei.link.get_fluxes()
            #Velocity,Tracers = kei.link.get_tracers()
            Vsave[...,nt],Tsave[...,nt] = kei.link.get_tracers()
            output.store_step_outvars(kei,nt)

        # write output netCDF file
        output.store_tracers(Vsave,Tsave,Flxsave,params.p['leco'])
        outfile = os.path.join(self.output_path,run_name+'.nc')
        print('Saving KEI output:',outfile)
        output.write(outfile,Finterp,params)


if __name__ == '__main__':

    params = kei_parameters()
    kf_ds = kei_forcing(r'/Users/blsaenz/KEI_run/DATA1/kf_200_100_2000.nc',start_date='2000-01-01', freq='H',legacy_nc=True)
    k = kei_simulation(kf_ds,t_start='2000-01-15',t_end='2000-08-15')
    k.compute(params,r'/Users/blsaenz/temp/keipy_output',run_name='keipytest5')

