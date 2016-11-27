from math import pi, cos, sin, sqrt, atan, tan

from visual import *

from visual.graph import *



class Body():

    def __init__(self,e,a,P,m):

        self.e = e # eccentricity of orbot

        self.a = a # semi -major axis of orbit

        self.P = P # period of orbit

        self.T = 0 # last time it past periapsis

        self.inclination = 0 # orbit inclination

        self.r = 0 # current radius 

        self.mass = m #mass of body

        self.G = 6.674e-11 #newton gravitational constant

#        self.long_periheleon = long_periheleon #Longitude of perihelion

 #       self.long_accension = long_accension #Longitude of accension

        

    def mean_anomaly(self, t, T, P):

        return (2 * pi * (t - T))/ P

    

    def eccentric_anomaly(self, m, e):

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

    

    def true_anomaly(self, e, E):

        s = sqrt((1+e)/(1-e))

        return 2 * atan(s*tan(E/2))



    def polar_position(self, a, e, v):

        numerator = a * (1-e*e)

        denom = 1 + e * cos(v)

        return numerator/denom



    def transform_to_helioecentric(self, x, y):

        return x,y



    def get_position(self, time, parentBody=None):

        M = self.mean_anomaly(time, self.T, self.P)

        E = self.eccentric_anomaly(M, self.e)

        v = self.true_anomaly(self.e, E)

        r = self.polar_position(self.a, self.e, v)

        x, y = r * cos(v), r*sin(v)

        if parentBody:

            x += parentBody.pos.x

            y += parentBody.pos.y

        self.pos = vector(x,y,0)

        if x == 0:

            self.T = time

        return x,y

    

    def kinetic_energy(self, parentBody):

        dist = mag(parentBody.pos - self.pos)

        KE = (0.5 * self.G * parentBody.mass * self.mass) / dist

        return KE

        

def hrsTodays(hrs):

    return hrs/24



P = 2122 # 88 days = 2122 hrs

t = 0 #hrs

T = 0

e = 0.21

a = 57.91e6

#e,a,P,m

mercury_ = Body(0.2056,57.91e6,2122,3.285e23) # object to save the data

mercury = sphere(radius=2440, make_trail=True)

mercury.pos = vector(mercury_.get_position(t))

mercury_label = label(pos=mercury.pos,text='Mercury', xoffset=12,yoffset=12,font='sans')

#print("{}, {}".format(x,y))

sun = sphere(pos=vector(0,0,0),radius=695700, color=color.orange)

sun.mass = 1.989e30 #kg

# satellite orbit info

sat_ = Body(0.21,10139.6e3,12,1200) # object to save the data

sat = sphere(radius=2.4397e6/2, make_trail=True,retain=80)

# earth params



earth_ = Body(0.0167,149.60e6,8760,5.9723e24)# earth data

earth = sphere(radius=6371, make_trail=true)

earth.material = materials.earth

earth.pos = earth_.get_position(t)

earth_label = label(pos=earth.pos,text='Earth', font='sans',xoffset=-12,yoffset=10)

#f1 = gcurve(color=color.green,width=600)

#venus_ = Body(0.0067,108.21e6,5400,5.9723e24)# earth data

#venus = sphere(radius=2.4397e6, make_trail=true)

#venus.pos = venus_.get_position(t)

#time text

time_label = label(pos=(0,0.25,0), text='T +0',align='upper_left')


# ellipse temp
b = 10139.6e3 * math.sqrt(1-0.21**2)
el = shapes.ellipse(width=10139.6e3, height=b)
el.pos = sat.pos
while True:

    rate(100)

    #print("{}, {}".format(x,y))

    mercury.pos = mercury_.get_position(t)

    mercury_label.pos=mercury.pos

    sat.pos = sat_.get_position(t, mercury)
    el.pos = mercury.pos 

    earth.pos = earth_.get_position(t)

    earth_label.pos = earth.pos
    

 #   venus.pos = venus_.get_position(t)

##    f1.plot(pos=[t,mercury_.kinetic_energy(sun)])

    time_label.text = "T +{} days".format(hrsTodays(t))

    t += 0.2

    
