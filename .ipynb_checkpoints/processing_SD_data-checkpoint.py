####import packages 
import pandas as pd
import SD_Project as SD
import numpy as np
import gsw
import seawater as sw
import xarray as xr
import datetime

import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
print('loaded packages')

#data cleaning
## load data in
# import key variables from data sets 
names_20=['latitude','longitude','BARO_PRES_MEAN','TEMP_AIR_MEAN','RH_MEAN','UWND_MEAN','VWND_MEAN','TEMP_CTD_MEAN','SAL_MEAN']
names_22=['latitude','longitude','BARO_PRES_MEAN','TEMP_AIR_MEAN','RH_MEAN','UWND_MEAN','VWND_MEAN','TEMP_CTD_RBR_MEAN','SAL_RBR_MEAN']
names_new=['latitude','longitude','pressure','airtemp','humidity','uwind','vwind','sst','salinity']


filepath_20= 'C:/Users/Hanna/Python/saildrone-gen_5-antarctica_circumnavigation_2019-sd1020-20190119T040000-20190803T043000-1_minutes-v1.1564884498845.nc'
DS_20= xr.open_dataset(filepath_20)
DS_20= DS_20.rename(dict(zip(names_20,names_new)))[names_new].squeeze('trajectory')

filepath_22= 'D:/Studium_der_Physik/M.Sc/Oceanography_Master/Masterthesis/SD_data/1022_new_merged.nc'
DS_22= xr.open_dataset(filepath_22)
DS_22= DS_22.rename(dict(zip(names_22,names_new)))[names_new]

filepath_23= 'D:/Studium_der_Physik/M.Sc/Oceanography_Master/Masterthesis/SD_data/1023_new_merged.nc'
DS_23= xr.open_dataset(filepath_23)
DS_23= DS_23.rename(dict(zip(names_22,names_new)))[names_new]

print('uploaded data')
# make time a coordinate and drop trajectory(except for 1020, there it is droped later) in each dataset
# explaination: the time vector provided in the Sd datasets were not contious, this might differ in other SD datasets
t_20 = DS_20['time'].values
df_20 = DS_20.to_dataframe()
df_20.index = t_20
df_20.index.names=['time']
df_20 = df_20.drop(columns=['time'])

df_22 = DS_22.to_dataframe().drop(columns=['trajectory'])
df_23 = DS_23.to_dataframe().drop(columns=['trajectory'])

def make_cont_timevector(df_index):
    '''df_index: dataframe index [datetime64[ns]], df.index 
    output: contious timevector with same start and end date as dataframe '''
    import numpy as np
    td_full= (df_index[-1]- df_index[0])
    t_min= td_full.days*24*60
    td_rest= (df_index[-1]-np.timedelta64(td_full.days,'D'))-df_index[0]
    t_min= t_min + td_rest.seconds//60
    
    date = np.array(df_index[0], dtype=np.datetime64)
    vec= date + pd.to_timedelta(np.arange(t_min+1), 'm')
    return vec

# make contious time vectors for all three datasets, reproject onto this vector, fill with nans 
df_cont={}
df= {'1020':df_20, '1022':df_22, '1023':df_23}
for frame in df:
    t_new= make_cont_timevector(df[frame].index)
    df_new= pd.DataFrame(index=t_new)
    df_new.index.names=['time']
    for var in df[frame].columns: 
        df_new[var] = pd.Series(df[frame][var], index=df[frame].index)
    df[frame]=df_new

# drop trajectory column for 1020  
df['1020'] = df['1020'].drop(columns=['trajectory'])

print('fixed time vectors')
###calculate 5 min rolling mean of measurments 
for d in df:
    df[d]['wind'] = SD.wind_to_ref_height(SD.wind(df[d]['uwind'],df[d]['vwind']), 3.6, 10) # correct wind speed to ref hight of 10m
    df[d]['wind'] = df[d]['wind'].rolling(5, center=True).mean()
    df[d]['airtemp'] =  df[d]['airtemp'].rolling(5, center=True).mean()
    df[d]['sst'] =  df[d]['sst'].rolling(5, center=True).mean()
    df[d]['humidity'] = df[d]['humidity'].rolling(5, center=True).mean()
    df[d]['salinity'] = df[d]['salinity'].rolling(5, center=True).mean()
    df[d]['humidity'] = df[d]['humidity'].rolling(5, center=True).mean()
    df[d]['pressure'] = df[d]['pressure'].rolling(5, center=True).mean()
    
    df[d]['Q_sens'] = SD.Q_sensible(df[d]['wind'], df[d]['airtemp'], df[d]['sst'])
    df[d]['Q_lat'] = -SD.Q_latent(df[d]['wind'], SD.humidity_spec_sat3(df[d]['humidity'], df[d]['airtemp'], df[d]['pressure']), SD.humidity_specific3(df[d]['humidity'], df[d]['airtemp'], df[d]['pressure']))
    df[d]['salinity']= SD.despike_sal(df[d]['salinity'], 0.1) # remove sal spikes: this takes lonnnng
# mask the complete row (corresponding to time) if one sensore mesurment is nan 
# explaination: there were missing long, lat, and intermitted missing values in measurments, 
# to calculate 
for d in df: 
    mask= np.ones(df[d].shape[0])
    for i in df[d].columns:
        mask= np.zeros(df[d].shape[0])*df[d][i]+mask

    for i in df[d].columns:
        df[d][i]=df[d][i]*mask
        
# calculate variables from measurments 
for d in df:
    df[d]['sensor_p']=np.zeros(df[d]['wind'].shape[0]) #create wind vector
    df[d]['sensor_p'][:]=0.5 #sensor depth
    df[d]['density']= gsw.rho(gsw.SA_from_SP(df[d]['salinity'],  df[d]['sensor_p'], df[d]['longitude'], df[d]['latitude']), gsw.CT_from_t(df[d]['salinity'], df[d]['sst'],  df[d]['sensor_p']),  df[d]['sensor_p'])
    #df[d]['dist_covered']= np.zeros(len(df[d]['latitude']))
    df[d]['dist_cov'] = pd.Series(sw.dist(df[d]['latitude'], df[d]['longitude'], units='km')[0],index=df[d].index[1:]) #distance covered between measuemnts 
    df[d]['dist_cov'] = df[d]['dist_cov'].replace(0.0, np.nan)
    df[d]['dist_NZ'] = df[d]['dist_cov'].cumsum()  #calculate track distance 
    df[d]['dist_cov'] = df[d]['dist_cov'].where(df[d]['dist_cov']>0.015)
    df[d]['sensor_p']=np.zeros((df[d]['wind'].shape[0]))# 
    df[d]['sensor_p'][:]=0.5 #sensor depth
    df[d]['alpha'] = gsw.alpha(gsw.SA_from_SP(df[d]['salinity'], df[d]['sensor_p'], df[d]['longitude'], df[d]['latitude']), gsw.CT_from_t(df[d]['salinity'], df[d]['sst'], df[d]['sensor_p']), df[d]['sensor_p'])
    df[d]['beta']= gsw.beta(gsw.SA_from_SP(df[d]['salinity'], df[d]['sensor_p'], df[d]['longitude'], df[d]['latitude']), gsw.CT_from_t(df[d]['salinity'], df[d]['sst'], df[d]['sensor_p']), df[d]['sensor_p'])

    R=np.abs((df[d]['alpha'][:-1]*np.diff(df[d]['sst']))/(df[d]['beta'][:-1]*np.diff(df[d]['salinity'])))
    df[d]['R'] = pd.Series(R, index=df[d].index[1:])
    
    df[d]['density_grad'] = (np.abs(np.diff(df[d]['density'])))/(df[d]['dist_cov'][1:])

#save data as csv

df['1020'].to_csv('df1020.csv')
df['1022'].to_csv('df1022.csv')
df['1023'].to_csv('df1023.csv')


print('csv files saved')