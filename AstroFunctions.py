"""
AstroFunctions.py

Astrodynamics Functions 

(Only a few functions from my original AstroFunctions.py are included here. )

Description:
Various Python functions useful for astrodynamics applications.
Most of these functions are based on Fundamentals of Astrodynamics, Vallado. 4th ed.

Author: Ashiv Dhondea, RRSG, UCT.
Date: 05 December 2016
"""
# ------------------------------------------------------------------------------------------ #
import numpy as np
import math
# ------------------------------------------------------------------------------------------ #
def fnSeconds_To_Hours(time_period):
    """
    Convert from seconds to hours, minutes and seconds.
    
    Date: 16 October 2016
    """
    num_hrs = int(time_period/(60.*60.));
    time_period =time_period - num_hrs*60.*60.;
    num_mins = int(time_period/60.);
    num_secs = time_period - num_mins*60.;
    return num_hrs,num_mins,num_secs # edit: 1/12/16: float division and multiplication

def fn_Convert_Datetime_to_GMST(datetime_object):
    """
    Converts a date and time in the datetime object form
    to GMST.
    
    Date: 05 October 2016
    """
    obj = datetime_object;
    julianday =  fnJulianDate(obj.year,obj.month,obj.day,obj.hour,obj.minute,obj.second);
    theta_GMST =  fn_Calculate_GMST(julianday);
    return theta_GMST # validated with example 3-5 in vallado.
    
def fnZeroTo2Pi(rotangle):
    """
    Wraps angle to fit in [0,2pi).
    Works in [rad] not [deg]
    Date: 7 October 2016
    """
    wrappedangle = rotangle % (2*math.pi);
    return wrappedangle
    
def fnECItoECEF(ECI,theta):
    ECEF = np.zeros([3],dtype=np.float64);
    # Rotating the ECI vector into the ECEF frame via the GST angle about the Z-axis
    ECEF = np.dot(fnRotate3(theta),ECI);
    return ECEF

def fnJulianDate(yr, mo, d, h, m, s):
    """
    Implements Algo 14 in Vallado book: JulianDate
    Date: 05 October 2016
    """
    JD = 367.0*yr - int((7*(yr+ int((mo+9)/12)))/4.0) + int((275.0*mo)/9.0) + d+ 1721013.5 + ((((s/60.0)+m)/60+h)/24.0);
    return JD # validated with example 3-4 in vallado.
    
def fn_Calculate_GMST(JD):
    """
    Calculates the Greenwich Mean Sidereal Time according to eqn 3-47 on page 188 in Vallado.
    Date: 05 October 2016
    Edit: 06 October 2016: CAUTION: theta_GMST is output in [degrees] rather than in [radians], 
    unlike most of the angles in this file.
    """
    T_UT1 = (JD - 2451545.0)/36525.0
    theta_GMST = 67310.54841 + (876600.0*60*60 + 8640184.812866)*T_UT1 + 0.093104 * T_UT1**2 - 6.2e-6 * T_UT1**3;
    
    while theta_GMST > 86400.0:
        theta_GMST = theta_GMST - 86400;

    theta_GMST = theta_GMST/240.0;
    theta_GMST = theta_GMST - 360; # in [deg] not [rad] !!!!!!!!!!!!!!!!!!!!!!!!!!!
    return theta_GMST # validated with example 3-5 in vallado.

def fnRotate3(alpha_rad):
    T = np.array([[ math.cos(alpha_rad),math.sin(alpha_rad),0], 
                  [-math.sin(alpha_rad),math.cos(alpha_rad),0],
                  [                 0,                0,1]],dtype=np.float64);
    return T # Validated against Vallado's example 2-3. 20/06/16

def fnCarts_to_LatLon(R):
    """
    function which converts ECEF position vectors to latitude and longitude
    Based on rvtolatlong.m in Richard Rieber's orbital library on mathwork.com
    
    Note that this is only suitable for groundtrack visualization, not rigorous 
    calculations.
    Date: 18 September 2016
    """
    r_delta = np.linalg.norm(R[0:1]);
    sinA = R[1]/r_delta;
    cosA = R[0]/r_delta;

    Lon = math.atan2(sinA,cosA);

    if Lon < -math.pi:
        Lon = Lon + 2*math.pi;

    Lat = math.asin(R[2]/np.linalg.norm(R));
    return Lat,Lon
