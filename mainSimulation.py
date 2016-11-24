from math import pi, cos, sin, sqrt, atan, tan
from visual import *

def mean_anomaly(t, T, P):
    return (2 * pi * (t - T))/ P
    
def eccentric_anomaly(m, e):
    """
    M = mean anomaly
    E = eccentric anomaly
    e = eccentricity

    M = E - esinE
    """
    def M(E): return E - (e * sin(E))
    E = 0
    while m > M(E):
        E += 1
    while M(E) > m:
        E -= 0.00001
    return E
    
def true_anomaly(e, E):
    s = sqrt((1+e)/(1-e))
    return 2 * atan(s*tan(E/2))

def polar_position(a, e, v):
    numerator = a * (1-e*e)
    denom = 1 + e * cos(v)
    return numerator/denom

P = 88
t = 0
T = 0
e = 0.21
a = 57.91e6
earth = sphere(make_trail=True,radius=2.4397e6,material=materials.earth)
M = mean_anomaly(t, T, P)
E = eccentric_anomaly(M, e)
v = true_anomaly(e, E)
r = polar_position(a, e, v)
x, y = r * cos(v), r*sin(v)
#print("{}, {}".format(x,y))
earth.pos = vector(r*cos(v), r*sin(v),0)
sun = sphere(pos=vector(0,0,0),radius=2.4397e6, color=color.orange)
# satellite orbit info
S_a = 10139.6e3
S_P = 15
S_T = 0
S_e = 0.2
sat = sphere(make_trail=True,radius=2.4397e6)
M = mean_anomaly(t,S_T, S_P)
E = eccentric_anomaly(M, S_e)
v = true_anomaly(S_e, E)
r = polar_position(S_a, S_e, v)
x, y = r * cos(v), r*sin(v)
sat.pos=vector(x,y,0)

while True:
    rate(50)
    M = mean_anomaly(t, T, P)
    E = eccentric_anomaly(M, e)
    v = true_anomaly(e, E)
    r = polar_position(a, e, v)
    x, y = r * cos(v), r*sin(v)
    #print("{}, {}".format(x,y))
    earth.pos = vector(x, y,0)
    #satellite case
    M = mean_anomaly(t, S_T, S_P)
    E = eccentric_anomaly(M, S_e)
    v = true_anomaly(S_e, E)
    r = polar_position(a, S_e, v)
    xx, yy = r * cos(v), r*sin(v)
    #print("{}, {}".format(x,y))
    sat.pos = vector(x + xx, y + yy,0)
    
    if x == 0:
        T = t

    if xx == 0:
        S_T = t
    t += 0.5
    
        
