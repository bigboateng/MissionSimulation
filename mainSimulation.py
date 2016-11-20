from visual import *
from math import cos,pi, atan, tan
t = 0
dt = 1000
twoPi= 2*pi

def computeMeanAnomaly(t, T, P):
    """
    @param t: current time
    @param T: last time the body was at periphelion
    @param P: period of orbit
    @returns: Mean anomaly - (Angle between body and periphillion if orbot was circular)
    """
    return (twoPi*(t-T))/(P)


def eccentricAnomaly(e, E):
    """
    @param e: orbit eccentricity
    @param E: mean anomally
    @returns: Mean anomaly 
    """
    v = 2*atan(sqrt((1+e)/(1-e)*tan(E/2)))
    return v

print("Eccentric(0.2,pi/10) = {}".format(eccentricAnomaly(0.2,pi/10)))
