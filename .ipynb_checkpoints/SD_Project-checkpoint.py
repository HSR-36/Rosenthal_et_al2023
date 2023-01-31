# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 13:03:05 2019

@author: Hanna
"""
    
def load_data(filepath):
    """ open single data file in nc format 
    
    uses squeeze to get rid of trajectory dimension from Saildronedata
        
    filepath example C:/Users/Hanna/Python/SDtest.nc
    """
    # single file
    import xarray as xr
    DS = xr.open_dataset(filepath)
    #print(DS.var)# DS.var, DS.dims, DS.coords, and DS.attrs
    DS = DS.squeeze(dim='trajectory')
    # get rid of trajectory 
    return DS

def wind(u,v): 
    """ calculates wind from eastward and westward wind components 
        u [m/s] array
        v [m/s] array
        wind=np.sqrt((u**2)+(v**2))
    """
    import numpy as np
    w= np.sqrt((u**2)+(v**2))
    return w 

def wind_to_ref_height(wind, m_height, ref_height):
    '''Correcting wind to desired reference height (usually 10m) from Schmidt et. al. 2017
    output: wind_new=wind*(np.log(Z/zo)/np.log(zm/zo))
    input. wind [m/s]
    m_height: measurment height [m]
    ref_height: (desired) refrence height [m]'''
    import numpy as np
    Z= ref_height # refrence height in m 
    zo= 1.52e-4 # roughness length [m]
    zm= m_height # measuremnt height [m] of sensor on Saildrone
    wind_new=wind*(np.log(Z/zo)/np.log(zm/zo))
    return wind_new
    

    
#%%       
def plot_timeseries_grid(x,y,ylabel):
    """ plots timeseries of variable x with grid on 
        x must be variable
        y must be time variable 
        ylabel must be string
        
    """
    import matplotlib.pyplot as plt

    #fig = Figure(figsize=(16,6), dpi=120)
    fig, host = plt.subplots()

    fig.subplots_adjust(right=1.6)

    plt.plot(x,y,linewidth=2, label='ylabel')
   
    fig.autofmt_xdate()

    plt.title("SD 1020")
    #plt.tick_params(bottom=False)
    plt.xlabel('') 
    plt.ylabel('ylabel')

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.grid(color='grey', linestyle='-', linewidth=.5)
    #plt.savefig("SD1020_tempertures.pdf", bbox_inches='tight')
    plt.show()
#%% 
def Q_sensible(wind,temp_air, SST):
    """ calculates sensible heatflux [W/m^2]
        with  Q_sens= ((rho_air*cpa*Cd)*wind*((SST) - temp_air))
        with rho_air = 1.275, density of air [kg/m3]
        cpa =1005, specific heat capacity of air [J/kgK] 
        Cd = 0.001, drag coefficient [dimensionless]
        from ocean perspective: if Q<0: ocean losing heat 
        
    """
    rho_air = 1.275 #density of air [kg/m3]
    cpa =1005 #specific heat capacity of air [J/kgK] 
    Cd = 0.001 #drag coefficient []
    Q_s= ((rho_air*cpa*Cd)*wind*((temp_air+273.16)-(SST+273.16)));
    return Q_s
#%%
def humidity_specific(Rel_humidity,temp_air):
    """
    calculates specific humidity from rel. humidity [%] and air temp[°C]
    spec_hum= (RH*a = eso*np.exp((Lv/Rv)*((1/T0)-(1/T))))
    input: 
        Rel. humidity [% kg/kg]
        temp_air [K]
    returns rel. humidity in % [kg/kg]   
    https://earthscience.stackexchange.com/questions/
    5076/how-to-calculate-specific-humidity-with-relative-humidity-temperature-and-pres
    """
    import numpy as np
    RH=Rel_humidity*0.01 # relative humidity in percent, convert to 1 = 100% 
    Lv= 2.501e6 #specific enthalpy of vaporization (J kg−1)
    Rv= 461.52 #[J/kgK] spec. gas constant water vapor 
    T0 = 273.16 #[K] refrence temperature in Kelvin 
    T =temp_air+ 273.16 #[K] Air temperature in Kelvin
    eso= 6.113 #[hPa] or 611.3 Pa saturated water vapour pressure at refrence temp T0 
    k= (Lv/Rv)*((1/T0)-(1/T))
    a = eso*np.exp(k)
    spec_hum= (RH*a)
    return spec_hum
#%%
def humidity_specific2(Rel_humidity,temp_air, pressure_air):
    """Same approach as engeniering toolbox 
    rel.humidity =(0.622 *Rel_humidity*pws)/(p-pws)
    with pws= (np.exp(77.3450 + 0.0057*T-(7235/T)))/T**8.2
         p = ((pressure / Ra*T) *(1 + y)) / (1 + (y*Rw/Ra))
    input: 
        rel. humidity [%]
        temp_air [°C]
        pressure_air [hPa]
        
    pws [kg/m^3] density of water vapor 
    p [kg / m^3] denisty of moist air 
    
   
    results in spec. humidity as [kg/kg] as ratio not % , add *100 to get % 
    """
    import numpy as np
    T= 273.16+ temp_air
    pressure= pressure_air*100 #convert hPa to Pa
    Ra = 286.9 #- individual gas constant air (J/kg K)
    Rw = 461.5 #- individual gas constant water vapor (J/kg K)
    y= Rel_humidity/100 # humidity ratio [kg/kg]
    p = ((pressure / Ra*T) *(1 + y)) / (1 + (y*Rw/Ra))
    pws= (np.exp(77.3450 + 0.0057*T-(7235/T)))/T**8.2
    x =(0.622 *Rel_humidity*pws)/(p-pws)
    return x
#%%
def humidity_specific_saturated(atm_pressure, air_temp):
    '''satured specific humidty from Louise via Slack
    pressure in hPa
    WMO, 2008 ew=6.112*exp(17,62t/(243,12+t))
    test
    '''
    import numpy as np
    #These equations used for saturated humidity
    #Using WMO (2008) definition for ew
    ep = 0.62198 # ratio of molecular weight of water and dry air
    
    ew=6.112*np.exp((17.62*air_temp)/(243.12+air_temp)) #MWO, 2008 [ew]=hPa, [t]=°C=airtemp,
    qs= (ep*ew)/atm_pressure
    
    #ew1 = 6.112*np.exp((17.62*(SStemperature))/(243.12+SStemperature)) #in hPa
    #ew = ew1*100
    #Using AOMIP definition of specific humidity. Can't find description of where 0.378 comes from
    #shum_sat = (ep*ew)/((pressure) - (0.378*ew))
    return qs
#%%
def humidity_specific3(Rel_humidity,temp_air, pressure_air):
    """
    calculate spec. humidity from rel. humidity, air temperature and air pressure
    q=(0.622*e)/(p-e)/(1+(0.622*Rel_humidity*611*np.exp(Lv/Rv*((1/T0)-(1/T))))/(p-Rel_humidity*611*np.exp(Lv/Rv*((1/T0)-(1/T)))))
    
    input: 
        Rel. humidity [%]
        temp_air [°C]
        pressure _ air [hPa]
    source: ametsoc , american meterological society and Rebecca :) 
    """
    import numpy as np
    p=pressure_air*100 #convert from [hPa] to [Pa]
    T=273.16+temp_air #convert from [°C] to [K]
    T0= 273.16#refrence temperature in [K]
    Lv= 2.501e6 #specific enthalpy of vaporization (J kg−1)
    Rv= 461.52 #[J/kgK] spec. gas constant water vapor 
    e=(Rel_humidity/100)*611*np.exp(Lv/Rv*((1/T0)-(1/T)))
    p=pressure_air*100# pressure_air from hPa in Pa
    rv=(0.622*e)/(p-e) # rv=mixing ratio
    q=rv/(1+rv)
    return q
#%%
def humidity_spec_sat3(Rel_humidity,temp_air,pressure_air):
    import numpy as np
    p=pressure_air*100 #convert from [hPa] to [Pa]
    T=273.16+temp_air #convert from [°C] to [K]
    T0= 273.16#refrence temperature in [K]
    Lv= 2.501e6 #specific enthalpy of vaporization (J kg−1)
    Rv= 461.52 #[J/kgK] spec. gas constant water vapor 
    e=(Rel_humidity/100)*611*np.exp(Lv/Rv*((1/T0)-(1/T)))
    p=pressure_air*100# pressure_air from hPa in Pa
    rv=(0.622*e)/(p-e) # rv=mixing ratio
    rs=rv/(Rel_humidity/100)
    q_sat=rs/(1+rs)
    return q_sat
#%%
#def humidity_sepc_sat_altern(temp_air, pressure_air):
#return q_sat
def Q_latent(wind,shum_sat,shum):
    """
    returns latent heat flux: Q = rho_air*Lv*Cd*wind*(shum_sat-shum) [W/m^2]
    input:
        wind [m/s]
        shum= specific humidity 
    """
    rho_air=1.275 #density of air [kg/m3]
    Cd= 0.001 #drag coefficient 
    Lv= 2.501e6 #[J/kg] latent heat of evaporation of water 
    Q = rho_air*Lv*Cd*wind*(shum_sat-shum)               #CHECKING
    return Q
#%%
def plot_2y(x,y1,y2,label1,label2,c1,c2,xlim,fig_title,pdf_name,save=False): 
    """
    plot_2y(x,y1,y2,label1,label2,c1,c2,xlim,fig_title,pdf_name,save=False)
    plot for two y on the same x-axis 
    label='name'
    c1= 'color'
    x= time vektor 
    xlim = [datetime.date(2019, 1, 19), datetime.date(2019, 3, 27)]
    [start, end] in datetime
    pdf_name= 'name.pdf'
    
    save=False , saving figure as pdf is optional
    """
    
    import matplotlib.pyplot as plt
    fig, host = plt.subplots()

    fig.subplots_adjust(right=1.6)

    par1 = host.twinx()
    p1, = host.plot(x,y1, label=label1,color=c1)
    p2, = par1.plot(x,y2, label=label2,color= c2)

    host.set_ylabel(label1,color=c1)
    par1.set_ylabel(label2,color=c2)

    tkw = dict(size=4, width=1.5)

    fig.suptitle(fig_title, fontsize=16)

    host.set_xlim(xlim)

    host.tick_params(axis='y',colors=c1,**tkw)
    par1.tick_params(axis='y', colors=c2, **tkw)

    from matplotlib.dates import DateFormatter
    formatter = DateFormatter('%d-%m-%Y')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)

    fig.autofmt_xdate()
    if save == True:    
        fig.savefig(pdf_name, bbox_inches='tight')
    
    plt.show()
#%%
def plot_3y(x,y1,y2,y3,label1,label2,label3,c1,c2,c3,xlim,fig_title,pdf_name,save=False):
    """
    plot_3y(x,y1,y2,y3,label1,label2,label3,c1,c2,c3,xlim,fig_title,pdf_name)
    ploting three diffent variables over the same time 
    
    
    x = time array, datetime
    y1, y2, y3 must have same dimensions and same length as x
    label='name' string
    c1= 'color' ,string 
    xlim = [datetime.date(2019, 1, 19), datetime.date(2019, 3, 27)]
    [start, end] in datetime
    pdf_name= 'name.pdf', string 
    
    save=False , saving figure as pdf is optional
    """
    #import datetime
    import matplotlib.pyplot as plt
    fig, host = plt.subplots()

    fig.subplots_adjust(right=1.6)

    par1 = host.twinx()
    par2 = host.twinx()
    par2.spines["right"].set_position(("axes", 1.1))
    par2.spines["right"].set_visible(True)

    p1, = host.plot(x,y1, label=label1,color=c1)
    p2, = par1.plot(x,y2, label=label2,color=c2)
    p3, = par2.plot(x,y3, label=label3,color=c3)


    host.set_ylabel(label1),#color=c1)
    par1.set_ylabel(label2)
    par2.set_ylabel(label3)

    host.yaxis.label.set_color(c1)
    par1.yaxis.label.set_color(c2)
    par2.yaxis.label.set_color(c3)

    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=c1,**tkw)
    par1.tick_params(axis='y', colors=c2,**tkw)
    par2.tick_params(axis='y', colors=c3, **tkw)


    host.set_xlabel("Date[YYYY-MM-DD]")
    #lines = [p1, p2, p3]
    
    fig.suptitle(fig_title, fontsize=16,x=.9,y=1)
    host.set_xlim(xlim)

    fig.autofmt_xdate()

    from matplotlib.dates import DateFormatter
    formatter = DateFormatter('%d-%m-%Y')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    if save == True:    
        fig.savefig(pdf_name, bbox_inches='tight')
    plt.show()
#%%
def make_circle(): 
    import matplotlib.path as mpath
    import numpy as np
    theta=np.linspace(0,2*np.pi,100)
    center, radius= [0.5,0.5],0.5
    verts=np.vstack([np.sin(theta),np.cos(theta)]).T
    circle=mpath.Path(verts*radius+center)
    return circle
#%%
def plot_scatter_errorbars_color_correlation(x1, y1, z, xerro, yerro, xlabel, ylabel, cbarlabel, title, capsize1, m, l):
    """
    plot_scatter_errorbars_color_correlation(x1, y1, z, xerro, yerro, xlabel, ylabel, cbarlabel, title)
    
    x1, y1 , z :array like and need to have the same size
    
    xerro, yerro: errorbars if applicable 
    cbarlabel: colorbar label 
    title: figure title 
    capsize1: marker size in px 
    m: markersize
    l: linewidth
    
    """
    import matplotlib
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import stats
    from sklearn.metrics import r2_score
    from scipy.stats import ttest_ind
    #import scipy
    
    #make colormap for scatter plot
    sc = plt.scatter(x1,y1,s=0,c=z)
    #create colorbar according to the scatter plot
    clb = plt.colorbar(sc,label=cbarlabel)
    norm = matplotlib.colors.Normalize(vmin=min(z), vmax=max(z), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in z])

    #loop over each data point to plot
    for x, y, e, color in zip(x1, y1, xerro, time_color):
        plt.plot(x, y, 'o', color=color, lw=l, ms=m)
        ax= plt.errorbar(x, y,  xerr=e, yerr=yerro, capsize=capsize1, color=color, lw=l, ms=m)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    

    # plt.plot([0,25],[0,25], "r--",label='1:1')
    plt.plot([0,1500],[0,1500], "r--",label='1:1')
    #define function
    #creating OLS regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x1,y1)
    def linefitline(b):
        return intercept + slope * b
    line1=linefitline(x1)
    
    def rmse(y1,y2):
        from sklearn.metrics import mean_squared_error
        from math import sqrt
        rms = sqrt(mean_squared_error(y1, y2))
        return rms
    rmse_1=rmse(x1, y1)

    #calculate p-value (probability the null hypothesis is true): 
    t, p = ttest_ind(x1,y1)
    r2 = r2_score(y1, linefitline(x1))
    plt.plot(x1,line1, c = 'darkblue',label='Linear fit')
    plt.plot([], [], ' ', label= '$R^2=$'+ str(np.round(r2,decimals=3)))
    
    if p_value<0.0001:
        plt.plot([], [], ' ', label= 'p $\ll$ 0.001')
    else: 
        plt.plot([], [], ' ', label= 'p =' + str(np.round(p_value,decimals=3)))
    plt.plot([], [], ' ', label= 'RMSE =' + str(np.round(rmse_1,decimals=3)))
    plt.legend(loc='best')
    return clb, ax
#%%
def plot_scatter_errorbars_color(x1, y1, z, xerro, yerro, xlabel, ylabel, cbarlabel, title, capsize1, l, m):
    """
    plot_scatter_errorbars_color_correlation(x1, y1, z, xerro, yerro, xlabel, ylabel, cbarlabel, title)
    
    x1, y1 , z :array like and need to have the same size
    
    xerro, yerro: errorbars if applicable 
    cbarlabel: colorbar label 
    title: figure title 
    capsize1, size of markers in px
    m: markersize
    l: linewidth 
    
    """
    import matplotlib
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import stats
    from sklearn.metrics import r2_score
    from scipy.stats import ttest_ind
    #import scipy
    
    #make colormap for scatter plot
    sc = plt.scatter(x1,y1,s=0,c=z)
    #create colorbar according to the scatter plot
    clb = plt.colorbar(sc,label=cbarlabel)
    norm = matplotlib.colors.Normalize(vmin=min(z), vmax=max(z), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in z])

    #loop over each data point to plot
    for x, y, e, color in zip(x1, y1, xerro, time_color):
        plt.plot(x, y, 'o', color=color, lw=l, ms=m)
        plt.errorbar(x, y,  xerr=e, yerr=yerro, lw=l, capsize=capsize1, color=color, ms=m)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)

    plt.legend(loc='best')
#%%
def rmse(y1,y2):
    from sklearn.metrics import mean_squared_error
    from math import sqrt
    rms = sqrt(mean_squared_error(y1, y2))
    return rms
#%%
def nan_helper(y):
    import numpy as np
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]
#%%
def despike_sal(sal, th):
    """
    despike_sal(data, theshold)
    sal/data: data to despike (array like)
    th: threshold , maximum diff(data) to tolerate
    
    spikes larger than threshold are substituted with np.nan
    """
    import numpy as np
    sal_bar=np.zeros(len(sal))
    sal_bar[0]=sal[0]
    for i in range(0,len(sal)-1):
        if np.isfinite(sal_bar[i]):
            sal_current= sal[i]
            ds=sal[i+1]-sal_current
        else:
            sal_current= sal[i]
            ds=sal[i+1]-sal_current
        if np.abs(ds)> th: 
            sal_bar[i+1]=np.nan
        else: 
            sal_bar[i+1]=sal[i+1]
    return sal_bar
#%% NOOOOO
def find_fronts_grad2(x, t, heatflux, T_air, sal, w, frontsize_mean, front_grad, temp_max, wind_max, distance): 
    '''
    find_fronts_grad(x, t, heatflux, T_air, sal, w,  frontsize_mean, front_grad, temp_max, wind_max, distance)
    x: density gradient [kg/m^3*1/km]
    t: time 
    heatflux: Heatflux [W/m^2]
    T_air: air temperature [°C]
    sal: salinity [PSU]
    w: wind speed [m/s]
    frontsize_mean: minimal front size for detection
    front_grad: minimal front gradient for detection
    temp_max: max air temp gradient for detection
    wind_max: max wind speed grad for front detection
    distance: commulative distance from start [km]

    '''
    import numpy as np
    T_air_diff=np.diff(T_air)
    w_diff=np.diff(w)
    x_diff=np.diff(x)
    
    start=[]
    stop=[]
    front_indx=[]
    front_size= []
    front_mag=[]
    Q_change=[]
    
    
    def condition(x, x_diff, k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max):
        if np.abs(x_diff[k])>front_grad and (x[k]+x[k+1])/2 > frontsize_mean and np.abs(T_air_diff[k]) < temp_max and np.abs(w_diff[k]) < wind_max:
            return True 
        else:
            return False
    k=0                
    while k < len(x)-1: 
        if condition(x, x_diff,k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max):
            start.append(k)
            front_indx.append(k)
            f=1
            while x[k+f]> x[k]+(x[k]*0.05):
                front_indx.append(k+f)
                f=f+1
            stop.append(k+f-1)
           
            front_size.append(distance[stop[-1]]-distance[start[-1]])
            front_mag.append(np.nanmax(x[k:k+f]))
            Q_change.append(np.abs(heatflux[start[-1]]-heatflux[stop[-1]]))
            k=k+f
        else: 
            k=k+1
    return start, stop, front_indx, front_size, front_mag, Q_change

#%%DIESE FUNKTION
def find_fronts_grad_mean(x, t, heatflux, T_air_diff, sal, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, toleranz, t_ctd, R, min_dist, temp_front= False, min_length=False, min_points=4): 
    import numpy as np
    '''
    find_fronts_grad(x, t, heatflux, T_air, sal, w,  frontsize_mean, front_grad, temp_max, wind_max, distance, t_ctd)
    x: density gradient [kg/m^3*1/km]
    t: time 
    heatflux: Heatflux [W/m^2]
    T_air: air temperature [°C]
    sal: salinity [PSU]
    w: wind speed [m/s]
    frontsize_mean: minimal front size for detection
    front_grad: minimal front gradient for detection
    temp_max: max air temp gradient for detection
    wind_max: max wind speed grad for front detection
    distance: commulative distance from start [km]
    R: absolute density ratio
    Temp_front: True or False 
    min_dist: minimum distance of two consecutive points to detect front start [km]
    min_length: minimum lenght of front is 5 data points
    
    output: 0: start
            1: stop 
            2: front_index
            3: front size/length
            4: front magnitude
            5: Q_change over front
            6: sal change over front 
            7: t_ctd change over front
            8: abs Q_change over front
    '''
    start=[]
    stop=[]
    front_indx=[]
    front_size= []
    front_mag=[]
    Q_change=[]
    Q_abs_change=[]
    t_ctd_change=[]
    sal_change=[]
    true_p=[]
    x_diff=np.diff(x)
    def condition(x, x_diff, k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, min_dist):
        if np.isnan(x[k]):
            return False
        elif temp_front==False:
            if x_diff[k]>front_grad and (x[k]+x[k+1])/2 > frontsize_mean and np.abs(T_air_diff[k]) < temp_max and np.abs(w_diff[k]) < wind_max and (distance[k+1]-distance[k])>min_dist:
                return True 
            else:
                return False
        else:
            if x_diff[k]>front_grad and (x[k]+x[k+1])/2 > frontsize_mean and np.abs(T_air_diff[k]) < temp_max and np.abs(w_diff[k]) < wind_max and R[k]<1 and (distance[k+1]-distance[k])>min_dist:
                return True 
            else:
                return False
    k=0                
    while k < len(x)-2: 
        if condition(x, x_diff, k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, min_dist):
            true_points=[]
            true_points.append(k)
            
            f=1
            g=0
            front=True
            while front and k+f < len(x)-1:
                if (distance[k+f]-distance[k+g])<min_dist: 
                    f=f+1
                elif (distance[k+f]-distance[k+g])>=min_dist and (x[k+f]+x[k+g])/(f-g+1) > frontsize_mean:#x[k+f]> x[k]+(x[k]*toleranz):
                    g=f
                    true_points.append(k+g)
                    f=f+1
                else:
                    front=False
                    f=f-1
            """while x[k+f]> x[k]+(x[k]*toleranz):
                if (distance[k+f]-distance[k+g-1])>min_dist:
                    indx.append(k+g)
                    g=f
                    f=f+1
                else:
                    f=f+1
                    if k+f > len(x)-2: 
                    front_indx.append(k+f)
                    break"""
                 
            if min_length: 
                if len(true_points)>=min_points:
                    start.append(k)
                    for i in range(g+1):
                        front_indx.append(k+i)
                    stop.append(k+g)
                    for i in true_points:
                        true_p.append(i)
                    front_size.append(distance[stop[-1]]-distance[start[-1]])
                    front_mag.append(np.nanmax(x[k:k+g]))
                    Q_abs_change.append(np.abs(heatflux[start[-1]]-heatflux[stop[-1]]))
                    Q_change.append(heatflux[start[-1]]-heatflux[stop[-1]])
                    sal_change.append(sal[start[-1]]-sal[stop[-1]])
                    t_ctd_change.append(t_ctd[start[-1]]-t_ctd[stop[-1]])
                    
            else: 
                start.append(k)
                for i in range(g):
                    front_indx.append(k+i)
                stop.append(k+g)
                for i in true_points:
                    true_p.append(i)
               
                front_size.append(distance[stop[-1]]-distance[start[-1]])
                front_mag.append(np.nanmax(x[k:k+g]))
                Q_abs_change.append(np.abs(heatflux[start[-1]]-heatflux[stop[-1]]))
                Q_change.append(heatflux[start[-1]]-heatflux[stop[-1]])
                sal_change.append(sal[start[-1]]-sal[stop[-1]])
                t_ctd_change.append(t_ctd[start[-1]]-t_ctd[stop[-1]])
            k=k+f
        else: 
            k=k+1
    return start, stop, front_indx, front_size, front_mag, Q_change, sal_change, t_ctd_change, Q_abs_change,  true_p
#%%#TESTFKT
def find_fronts_grad_test(x, t, heatflux, T_air_diff, sal, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, toleranz, t_ctd, R, min_dist, temp_front= False, min_length=False, min_points=4): 
    import numpy as np
    '''
    find_fronts_grad(x, t, heatflux, T_air, sal, w,  frontsize_mean, front_grad, temp_max, wind_max, distance, t_ctd)
    x: density gradient [kg/m^3*1/km]
    t: time 
    heatflux: Heatflux [W/m^2]
    T_air: air temperature [°C]
    sal: salinity [PSU]
    w: wind speed [m/s]
    frontsize_mean: minimal front size for detection
    front_grad: minimal front gradient for detection
    temp_max: max air temp gradient for detection
    wind_max: max wind speed grad for front detection
    distance: commulative distance from start [km]
    R: absolute density ratio
    Temp_front: True or False 
    min_dist: minimum distance of two consecutive points to detect front start [km]
    min_length: minimum lenght of front is 5 data points
    
    output: 0: start
            1: stop 
            2: front_index
            3: front size/length
            4: front magnitude
            5: Q_change over front
               6: sal change over front 
            7: t_ctd change over front
            8: abs Q_change over front
    '''
    start=[]
    stop=[]
    front_indx=[]
    front_size= []
    front_mag=[]
    Q_change=[]
    Q_abs_change=[]
    t_ctd_change=[]
    sal_change=[]
    true_p=[]
    x_diff=np.diff(x)
    def condition(x, x_diff, k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, min_dist):
        if np.isnan(x[k]):
            return False
        elif temp_front==False:
            if x_diff[k]>front_grad and x[k] > frontsize_mean and np.abs(T_air_diff[k]) < temp_max and np.abs(w_diff[k]) < wind_max and (distance[k+1]-distance[k])>min_dist:
                return True 
            else:
                return False
        else:
            if x_diff[k]>front_grad and x[k] > frontsize_mean and np.abs(T_air_diff[k]) < temp_max and np.abs(w_diff[k]) < wind_max and R[k]<1 and (distance[k+1]-distance[k])>min_dist:
                return True 
            else:
                return False
    k=0                
    while k < len(x)-2: 
        if condition(x, x_diff, k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, min_dist):
            true_points=[]
            true_points.append(k)
            
            f=1
            g=0
            front=True
            while front and k+f < len(x)-1:
                if (distance[k+f]-distance[k+g])<min_dist: 
                    f=f+1
                elif (distance[k+f]-distance[k+g])>=min_dist and x[k+f] > frontsize_mean:#x[k+f]> x[k]+(x[k]*toleranz):
                    g=f
                    true_points.append(k+g)
                    f=f+1
                else:
                    front=False
                    f=f-1
            """while x[k+f]> x[k]+(x[k]*toleranz):
                if (distance[k+f]-distance[k+g-1])>min_dist:
                    indx.append(k+g)
                    g=f
                    f=f+1
                else:
                    f=f+1
                    if k+f > len(x)-2: 
                    front_indx.append(k+f)
                    break"""
                 
            if min_length: 
                if len(true_points)>=min_points:
                    start.append(k)
                    for i in range(g+1):
                        front_indx.append(k+i)
                    stop.append(k+g)
                    for i in true_points:
                        true_p.append(i)
                    front_size.append(distance[stop[-1]]-distance[start[-1]])
                    front_mag.append(np.nanmax(x[k:k+g]))
                    Q_abs_change.append(np.abs(heatflux[start[-1]]-heatflux[stop[-1]]))
                    Q_change.append(heatflux[start[-1]]-heatflux[stop[-1]])
                    sal_change.append(sal[start[-1]]-sal[stop[-1]])
                    t_ctd_change.append(t_ctd[start[-1]]-t_ctd[stop[-1]])
                    
            else: 
                start.append(k)
                for i in range(g):
                    front_indx.append(k+i)
                stop.append(k+g)
                for i in true_points:
                    true_p.append(i)
               
                front_size.append(distance[stop[-1]]-distance[start[-1]])
                front_mag.append(np.nanmax(x[k:k+g]))
                Q_abs_change.append(np.abs(heatflux[start[-1]]-heatflux[stop[-1]]))
                Q_change.append(heatflux[start[-1]]-heatflux[stop[-1]])
                sal_change.append(sal[start[-1]]-sal[stop[-1]])
                t_ctd_change.append(t_ctd[start[-1]]-t_ctd[stop[-1]])
            k=k+f
        else: 
            k=k+1
    return start, stop, front_indx, front_size, front_mag, Q_change, sal_change, t_ctd_change, Q_abs_change,  true_p
#%%
def find_fronts_grad(x, t, heatflux, T_air_diff, sal, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, toleranz, t_ctd, R, min_dist, temp_front= False, min_length=False, min_points=4): 
    import numpy as np
    '''
    find_fronts_grad(x, t, heatflux, T_air, sal, w,  frontsize_mean, front_grad, temp_max, wind_max, distance, t_ctd)
    x: density gradient [kg/m^3*1/km]
    t: time 
    heatflux: Heatflux [W/m^2]
    T_air: air temperature [°C]
    sal: salinity [PSU]
    w: wind speed [m/s]
    frontsize_mean: minimal front size for detection
    front_grad: minimal front gradient for detection
    temp_max: max air temp gradient for detection
    wind_max: max wind speed grad for front detection
    distance: commulative distance from start [km]
    R: absolute density ratio
    Temp_front: True or False 
    min_dist: minimum distance of two consecutive points to detect front start [km]
    min_length: minimum lenght of front is 5 data points
    
    output: 0: start
            1: stop 
            2: front_index
            3: front size/length
            4: front magnitude
            5: Q_change over front
            6: sal change over front 
            7: t_ctd change over front
            8: abs Q_change over front
    '''
    start=[]
    stop=[]
    front_indx=[]
    front_size= []
    front_mag=[]
    Q_change=[]
    Q_abs_change=[]
    t_ctd_change=[]
    sal_change=[]
    true_p=[]
    x_diff=np.diff(x)
    def condition(x, x_diff, k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, min_dist):
        if np.isnan(x[k]):
            return False
        elif temp_front==False:
            if x_diff[k]>front_grad and x[k] > frontsize_mean and np.abs(T_air_diff[k]) < temp_max and np.abs(w_diff[k]) < wind_max and (distance[k+1]-distance[k])>min_dist:
                return True 
            else:
                return False
        else:
            if  x[k] > frontsize_mean and np.abs(T_air_diff[k]) < temp_max and np.abs(w_diff[k]) < wind_max and R[k]<1 and (distance[k+1]-distance[k])>min_dist:
                return True 
            else:
                return False
    k=0                
    while k < len(x)-2: 
        if condition(x, x_diff, k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, min_dist):
            true_points=[]
            true_points.append(k)
            
            f=1
            g=0
            front=True
            while front and k+f < len(x)-1:
                if (distance[k+f]-distance[k+g])<min_dist: 
                    f=f+1
                elif (distance[k+f]-distance[k+g])>=min_dist and x[k+f]> frontsize_mean:#x[k+f]> x[k]+(x[k]*toleranz):
                    g=f
                    true_points.append(k+g)
                    f=f+1
                else:
                    front=False
                    f=f-1
            """while x[k+f]> x[k]+(x[k]*toleranz):
                if (distance[k+f]-distance[k+g-1])>min_dist:
                    indx.append(k+g)
                    g=f
                    f=f+1
                else:
                    f=f+1
                    if k+f > len(x)-2: 
                    front_indx.append(k+f)
                    break"""
                 
            if min_length: 
                if len(true_points)>=min_points:
                    start.append(k)
                    for i in range(g+1):
                        front_indx.append(k+i)
                    stop.append(k+g)
                    for i in true_points:
                        true_p.append(i)
                    front_size.append(distance[stop[-1]]-distance[start[-1]])
                    front_mag.append(np.nanmax(x[k:k+g]))
                    Q_abs_change.append(np.abs(heatflux[start[-1]]-heatflux[stop[-1]]))
                    Q_change.append(heatflux[start[-1]]-heatflux[stop[-1]])
                    sal_change.append(sal[start[-1]]-sal[stop[-1]])
                    t_ctd_change.append(t_ctd[start[-1]]-t_ctd[stop[-1]])
                    
            else: 
                start.append(k)
                for i in range(g):
                    front_indx.append(k+i)
                stop.append(k+g)
                for i in true_points:
                    true_p.append(i)
               
                front_size.append(distance[stop[-1]]-distance[start[-1]])
                front_mag.append(np.nanmax(x[k:k+g]))
                Q_abs_change.append(np.abs(heatflux[start[-1]]-heatflux[stop[-1]]))
                Q_change.append(heatflux[start[-1]]-heatflux[stop[-1]])
                sal_change.append(sal[start[-1]]-sal[stop[-1]])
                t_ctd_change.append(t_ctd[start[-1]]-t_ctd[stop[-1]])
            k=k+f
        else: 
            k=k+1
    return start, stop, front_indx, front_size, front_mag, Q_change, sal_change, t_ctd_change, Q_abs_change,  true_p
    
#%%#%%#FIND NOT A FRONT (to calculate)
def find_background(x, t, heatflux, T_air_diff, sal, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, toleranz, t_ctd, R, min_dist, temp_front= False, min_length=False, min_points=4): 
    import numpy as np
    '''
    find_fronts_grad(x, t, heatflux, T_air, sal, w,  frontsize_mean, front_grad, temp_max, wind_max, distance, t_ctd)
    x: density gradient [kg/m^3*1/km]
    t: time 
    heatflux: Heatflux [W/m^2]
    T_air: air temperature [°C]
    sal: salinity [PSU]
    w: wind speed [m/s]
    frontsize_mean: minimal front size for detection
    front_grad: minimal front gradient for detection
    temp_max: max air temp gradient for detection
    wind_max: max wind speed grad for front detection
    distance: commulative distance from start [km]
    R: absolute density ratio
    Temp_front: True or False 
    min_dist: minimum distance of two consecutive points to detect front start [km]
    min_length: minimum lenght of front is 5 data points
    
    output: 0: start
            1: stop 
            2: front_index
            3: front size/length
            4: front magnitude
            5: Q_change over front
               6: sal change over front 
            7: t_ctd change over front
            8: abs Q_change over front
    '''
    start=[]
    stop=[]
    b_indx=[]
    b_Qgrad=[]
    
    def condition(x, k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, min_dist):
        if np.isnan(x[k]):
            return False
        elif temp_front==False:
            if  x[k] < frontsize_mean and np.abs(T_air_diff[k]) < temp_max and np.abs(w_diff[k]) < wind_max and (distance[k+1]-distance[k])>=min_dist:
                return True 
            else:
                return False
        else:
            if  x[k] < frontsize_mean and np.abs(T_air_diff[k]) < temp_max and np.abs(w_diff[k]) < wind_max and R[k]<1 and (distance[k+1]-distance[k])>=min_dist:
                return True 
            else:
                return False
    k=0                
    while k < len(x)-2: 
        if condition(x, k, T_air_diff, w_diff, frontsize_mean, front_grad, temp_max, wind_max, distance, min_dist):
            true_points=[]
            true_points.append(k)
            
            f=1
            g=0
            b=True
            while b and k+f < len(x)-1:
                if (distance[k+f]-distance[k+g])<min_dist: 
                    f=f+1
                elif x[k+f] < frontsize_mean and np.abs(T_air_diff[k+f]) < temp_max and np.abs(w_diff[k+f]) < wind_max:#x[k+f]> x[k]+(x[k]*toleranz):
                    g=f
                    true_points.append(k+g)
                    f=f+1
                else:
                    b=False
                    f=f-1
            
            if min_length: 
                if len(true_points)>=min_points:
                    start.append(k)
                    for i in range(g+1):
                        b_indx.append(k+i)
                    stop.append(k+g)
                    
            else: 
                start.append(k)
                for i in range(g):
                    b_indx.append(k+i)
                stop.append(k+g)
            
            q_back=[]
            d_back=[]
            for i in true_points:
                q_back.append(heatflux[i])
                d_back.append(distance[i])
            if q_back:
                b_Qgrad.append(np.nanmean(np.abs(np.diff(q_back))/np.diff(d_back)))
            
            k=k+f
        else: 
            k=k+1
    return start, stop, b_indx, b_Qgrad
#%%
def decorrelation_length(data, interval_length):
    import numpy as np
    '''"""data: one variable [array] to autocorrelate
    interval_length: on this length the data should be cut and then analysed on this scale for autocorr.
    returns: decorrelation length for each interval. (should be list of length: len(data)/interval_length) """'''
    def autocorr_diy(x, N): 
        b=[]
        for i in range(len(x)-1):
            m = np.nanmean(x)
            k= (x[i]-m)**2
            b.append(k)
        s= np.nansum(b)
        f=[]  
        for n in range(N): 
            l=[]
            for i in range(len(x)-n): 
                l.append((x[i]-m)*(x[i+n]-m))
            f.append(np.nansum(l))
        return np.array(f/s)

    intv = np.arange(0, len(data), interval_length)
    small_x= []
    for i in intv: 
        small_x.append(data[i:i+interval_length-1])
    
    # devide data in equally sized arrays: 
    small_corr= []
    for i in range(len(small_x)):
        small_corr.append(autocorr_diy(small_x[i], interval_length))
    
    # remove 0. from nanas and turn into nan
    for i in range(len(small_corr)):
        small_corr[i][small_corr[i]==0] = np.nan
        
    d_length=[]
    for i in range(len(small_corr)): 
        if np.all(np.isnan(small_corr[i]))==False and len(np.where(small_corr[i]<0)[0])>0:   
            d_length.append(np.where(small_corr[i] < 0)[0][0])
        else: 
            d_length.append(np.nan)
    return d_length
