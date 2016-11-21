from visual import *
import numpy
import math


class OrbitTools:
    def __init__(self):
        pass

    def mean_anomaly(self, t, p):
        """
        @param t: elapse time
        @param p: period of orbit
        returns:  mean anomaly
        M = 2pi * t
        -------
           P 
        """
        return (2 * math.pi * t) / p

    def eccentric_anomaly(self, m, e):
        """
        @param m: mean anomaly
        @param e: eccentricity
        returns E: eccentric anomaly
        
        M = E - esinE
        """
        def M(E): return E - (e * math.sin(E))
        E = 0
        while m > M(E):
            E += 1
        while M(E) > m:
            E -= 0.00001
        return E

    def v(self, e, E):
        """
        (1 - e)tan^2(theta/2) = (1 + e)tan^2(E/2)
        e = eccentricity
        theta = true anomaly
        E = eccentric anomaly
        """
        def l(theta): return (1-e)*(math.tan(theta/2))**2
        r = (1+e)*(math.tan(E/2))**2
        theta = 0
        while l(theta) < r:
            theta += 0.1
        while r < l(theta):
            theta -= 0.00001
        return [theta, 2*(math.pi - theta) + theta]
    
    def heliocentric_distance(self, a, e, E):
        """
        a = semi-major axis
        e = eccentricity
        E = eccentric anomaly

        r = a(1 - ecosE)
        """
        return a * (1 - (e * math.cos(E)))


    def calculate_position(self, e, t, p, a):
        M = self.mean_anomaly(t, p)
        E = self.eccentric_anomaly(M, e)
        if list(math.modf(float(t) / p))[0] > 0.5:
            theta = self.v(e, E)[1]
        if list(math.modf(float(t) / p))[0] < 0.5:
            theta = self.v(e, E)[0]
        r = self.heliocentric_distance(a, e, E)
        return [theta, r] 



t = 0
dt = 30
orbitTool = OrbitTools()
sun = sphere()
sun.radius = 6.95e8
sun = vector(0 , 0, 0)
sun.color = color.red
earth = sphere(make_trail=True)
earth.radius = 6.37e6
earth.material = materials.earth
r, theta = orbitTool.calculate_position(0.0167, 0, 365.25, 1.496e8)
earth.pos = vector(r*math.cos(theta), r*math.sin(theta), 0)
while True:
    rate(100)
    x = orbitTool.calculate_position(0.0167, t, 365.25, 1.496e8)
    r, theta = x[0], x[1]
    earth.pos = vector(r*math.cos(theta), r*math.sin(theta), 0)
    t += dt
    print(t)
    print(earth.pos)
    
