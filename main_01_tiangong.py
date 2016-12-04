"""

main_01_tiangong.py

Satellite ground track plotting.

Demonstration of satellite ground track plotting. 

Reads in the position vectors in the Earth Centered Inertial or the Earth Centered Earth Fixed Frame, 
converts these to latitude and longitude and plots the location on an equirectangular-projected map.

Author: Ashiv Dhondea, RRSG, UCT.
Date: 05 December 2016
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib import colors
import matplotlib.patches as patches
import matplotlib as mpl

import datetime as dt
import pytz
import aniso8601
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

import AstroFunctions as AstFn

# ------------------------------------------------------------------------------------------------ #
## TIANGONG 1              
tle_line1 = '1 37820U 11053A   16339.16302693  .00021051  00000-0  17357-3 0  9997';
tle_line2 = '2 37820  42.7632 305.9790 0015673 147.2225 332.7571 15.69763185297236';

line1 = (tle_line1);
line2 = (tle_line2);

satellite_obj = twoline2rv(line1, line2, wgs72);
print 'satellite number'
print satellite_obj.satnum
print 'epochyr'
print satellite_obj.epochyr
print 'epochdays'
print satellite_obj.epochdays
print 'jdsatepoch'
print satellite_obj.jdsatepoch
print 'epoch'
print satellite_obj.epoch
print 'inclination'
print math.degrees(satellite_obj.inclo)
print 'RAAN'
print math.degrees(satellite_obj.nodeo)
print 'eccentricity'
print satellite_obj.ecco
print 'argument of perigee'
print math.degrees(satellite_obj.argpo)
print 'mean anomaly'
print math.degrees(satellite_obj.mo)

delta_t = 1; #[s]
simulation_period = 720*60 ;#[s], 720 min of simulation
timevec = np.arange(0,simulation_period+delta_t,delta_t,dtype=np.float64);
x_state = np.zeros([6,len(timevec)],dtype=np.float64);
xecef = np.zeros([3,len(timevec)],dtype=np.float64);
lat = np.zeros([len(timevec)],dtype=np.float64);
lon = np.zeros([len(timevec)],dtype=np.float64);

index = 0;
current_time = timevec[index];
hrs,mins,secs = AstFn.fnSeconds_To_Hours(current_time + (satellite_obj.epoch.hour*60*60) + (satellite_obj.epoch.minute*60)+ satellite_obj.epoch.second);
dys = satellite_obj.epoch.day + int(math.ceil(hrs/24));     
if hrs >= 24:
    hrs = hrs - 24*int(math.ceil(hrs/24)) ;
        
satpos,satvel = satellite_obj.propagate(satellite_obj.epoch.year,satellite_obj.epoch.month,dys,hrs,mins,secs+(1e-6)*satellite_obj.epoch.microsecond);
x_state[0:3,index] = np.asarray(satpos);
x_state[3:6,index] = np.asarray(satvel);

tle_epoch_test = dt.datetime(year=satellite_obj.epoch.year,month=satellite_obj.epoch.month,day=int(dys),hour=int(hrs),minute=int(mins),second=int(secs),microsecond=0,tzinfo= pytz.utc);
theta_GMST =  math.radians(AstFn.fn_Convert_Datetime_to_GMST(tle_epoch_test));        
## Rotate ECI position vector by GMST angle to get ECEF position
theta_GMST = AstFn.fnZeroTo2Pi(theta_GMST);
xecef[:,index] = AstFn.fnECItoECEF(x_state[0:3,index],theta_GMST);
lat[index],lon[index] = AstFn.fnCarts_to_LatLon(xecef[:,index]);

for index in range(1,len(timevec)):
    current_time = timevec[index];
    hrs,mins,secs = AstFn.fnSeconds_To_Hours(current_time + (satellite_obj.epoch.hour*60*60) + (satellite_obj.epoch.minute*60)+ satellite_obj.epoch.second);
    dys = satellite_obj.epoch.day + int(math.ceil(hrs/24)); 
    
    if hrs >= 24:
        hrs = hrs - 24*int(math.ceil(hrs/24)) ;
        
    satpos,satvel = satellite_obj.propagate(satellite_obj.epoch.year,satellite_obj.epoch.month,dys,hrs,mins,secs+(1e-6)*satellite_obj.epoch.microsecond);
    x_state[0:3,index] = np.asarray(satpos);
    x_state[3:6,index] = np.asarray(satvel);

    tle_epoch_test = dt.datetime(year=satellite_obj.epoch.year,month=satellite_obj.epoch.month,day=int(dys),hour=int(hrs),minute=int(mins),second=int(secs),microsecond=0,tzinfo= pytz.utc);
    theta_GMST =  math.radians(AstFn.fn_Convert_Datetime_to_GMST(tle_epoch_test));        
    ## Rotate ECI position vector by GMST angle to get ECEF position
    theta_GMST = AstFn.fnZeroTo2Pi(theta_GMST);
    xecef[:,index] = AstFn.fnECItoECEF(x_state[0:3,index],theta_GMST);
    lat[index],lon[index] = AstFn.fnCarts_to_LatLon(xecef[:,index]);

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
## plot results     
coastline_data= np.loadtxt('Coastline.txt',skiprows=1)
# Coastline map data obtained from
# https://www.mathworks.com/matlabcentral/fileexchange/13439-orbital-mechanics-library/content/Groundtrack.m

w, h = plt.figaspect(0.5)
fig = plt.figure(figsize=(w,h))
ax = fig.gca()
plt.rc('text', usetex=True)
plt.rc('font', family='serif');
plt.rc('font',family='helvetica');
params = {'legend.fontsize': 8,
    'legend.handlelength': 2}
plt.rcParams.update(params)

groundtrack_title = satellite_obj.epoch.strftime('%d %B %Y')
fig.suptitle(r"\textbf{Tiangong-1 Ground Track on %s}" %groundtrack_title,fontsize=16)
plt.plot(coastline_data[:,0],coastline_data[:,1],'g');
ax.set_xlabel(r'Longitude $[\mathrm{^\circ}]$',fontsize=14)
ax.set_ylabel(r'Latitude $[\mathrm{^\circ}]$',fontsize=14)
plt.xlim(-180,180);
plt.ylim(-90,90);
plt.yticks([-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90]);
plt.xticks([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180]);

for index in range(0,len(timevec)):
    plt.plot(math.degrees(lon[index]),math.degrees(lat[index]),'b.',markersize=1);

ax.grid(True);
at = AnchoredText("AshivD",prop=dict(size=5), frameon=True,loc=4)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)
fig.savefig('main_01_tiangong_ground_track.pdf',format='pdf',bbox_inches='tight',pad_inches=0.01,dpi=1200);

