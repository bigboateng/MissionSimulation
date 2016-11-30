from math import pi, cos, sin, sqrt, atan, tan
from visual import *
from visual.graph import *
from visual.controls import *
import wx

window1Size = 600
Window = window(menu=false, title="Mission Simulation - By Boateng",
                x = 0, y = 0, width=window1Size, height=window1Size,
                style=wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX)
L = 320
Hgraph = 10
d = 20
w = window(width=2*(L+window.dwidth), height=L+window.dheight+window.menuheight+Hgraph,
           menus=True, title='Orbit Info',
           style=wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX)

display(window=Window, x=0, y=0,width=window1Size,height=window1Size)
controlsDisplay = w.panel # Refers to the full region of the window in which to place widgets


class Body():
    def __init__(self,e,a,P,mass):
        self.e = e # Eccentricity of orbit.
        self.a = a # Semi-Major axis of orbit.
        self.P = P # Period of orbit in hours.
        self.T = 0 # Last time body pass periapsis.
        self.inclination = 0 # Orbit inclination.
        self.mass = mass # Mass of body
        self.v = 0 #True anomaly
        self.G = 6.674e-11 #newton gravitational constant
        self.numOfOrbits = 0
        #self.long_periheleon = long_periheleon #Longitude of perihelion
        #self.long_accension = long_accension #Longitude of accension
        
    def mean_anomaly(self, t, T, P):
        """
        :param t: current time
        :param T: time it last went past periapsus
        :param P: orbit period
        :return: M
        M = (2 * pi * (T - t)) / P
        """
        return (2 * pi * (t - T))/ P
    
    def eccentric_anomaly(self, m, e):
        """
        :param m: mean anomaly
        :param e: eccentricity
        :return: eccentric anomaly (E)
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
        """
        :param e: eccentricity
        :param E: eccentric anomaly
        :return: true anomaly

        """
        s = sqrt((1+e)/(1-e))
        return 2 * atan(s*tan(E/2))

    def polar_position(self, a, e, v):
        """
        :param a: Semi-major axis
        :param e: Eccentricity
        :param v: True Anomaly
        :return: Position in polar coordinates
        """
        numerator = a * (1-e*e)
        denom = 1 + e * cos(v)
        return numerator/denom

    def transform_to_helioecentric(self, x, y):
        """

        :param x:
        :param y:
        :return:
        """
        return x,y

    def get_position(self, time, parentBody=None):
        """
        :param time: Current time
        :param parentBody: Body this one is orbiting around
        :return: Position in x,y plane
        """
        M = self.mean_anomaly(time, self.T, self.P)
        E = self.eccentric_anomaly(M, self.e)
        v = self.true_anomaly(self.e, E)
        self.v= v
        r = self.polar_position(self.a, self.e, v)
        x, y = r * cos(v), r*sin(v)
        if parentBody:
            x += parentBody.pos.x
            y += parentBody.pos.y
        self.pos = vector(x,y,0)
        if x == 0:
##            self.T = time
            self.numOfOrbits += 1
        return x,y

    def get_speed(self, parentBody):
        mu = self.G * parentBody.mass
        r = mag(parentBody.pos - self.pos)
        v = math.sqrt(mu* ((2/r)-(1/self.a)))
        return v
    
    def kinetic_energy(self, parentBody):
        """
        :param parentBody:
        :return: Kinetic energy of planet
        """
        dist = mag(parentBody.pos - self.pos)
        KE = (0.5 * self.G * parentBody.mass * self.mass) / dist
        return KE
    
    def get_anomaly(self):
        # convert to degres
        deg = (self.v * 180) / pi
        return deg

    def get_num_orbits(self):
        return self.numOfOrbits
        
def hrsTodays(hrs):
    return hrs/24

t= 0 # Time elapsed.
# Mercury orbit info/parameters.
mercury_ = Body(e=0.2056,a=57.91e6,P=2122,mass=3.285e23) # object to save the data
mercury = sphere(radius=2440,make_trail=True)
mercury.pos = vector(mercury_.get_position(t))
mercury_label = label(pos=mercury.pos,text='Mercury', xoffset=17,yoffset=12,font='sans')
# Sun orbit info/parameters.
sun = sphere(pos=vector(0,0,0),radius=695700, color=color.orange)
sun.mass = 1.989e30 #kg
# Satellite orbit info/parameters.
sat_ = Body(0.21,10139.6e3,12,1200)
sat = sphere(radius=2440,make_trail=True)
sat.pos=vector(sat_.get_position(t))
sat_label = label(pos=sat.pos,text='Satllite', font='sans',xoffset=-12,yoffset=25)
# Earth orbit info/parameters.
earth_ = Body(0.0167,149.60e6,8760,5.9723e24)
earth = sphere(radius=6371, make_trail=true)
earth.material = materials.earth
earth.pos = earth_.get_position(t)
earth_label = label(pos=earth.pos,text='Earth', font='sans',xoffset=-12,yoffset=15)
#f1 = gcurve(color=color.green,width=600)
# Venus orbit parameters.
venus_ = Body(e=0.0067,a=108.21e6,P=5400,mass=5.9723e24)# earth data
venus = sphere(radius=6052, make_trail=true)
venus.pos = venus_.get_position(t)
#time text
#time_label = label(pos=(0,0.25,0), text='T +0')
p = w.panel
timeText = wx.StaticText(p, pos=(d,4),size=(500,40), label='A 3D canvas',
              style=wx.ALIGN_CENTRE | wx.ST_NO_AUTORESIZE)
meanAnomaly = wx.StaticText(p, pos=(d,30),size=(500,40), label='A 3D canvas',
              style=wx.ALIGN_CENTRE | wx.ST_NO_AUTORESIZE)
eccenText = wx.StaticText(p, pos=(d,30*2),size=(500,40), label="e = 0.21",
              style=wx.ALIGN_CENTRE | wx.ST_NO_AUTORESIZE)
semiMajorText = wx.StaticText(p, pos=(d,30*3),size=(500,40), label="Semi-Major Axis=10139 km",
              style=wx.ALIGN_CENTRE | wx.ST_NO_AUTORESIZE)
speedText = wx.StaticText(p, pos=(d,30*4),size=(500,40), label="Speed = 100 km/s",
              style=wx.ALIGN_CENTRE | wx.ST_NO_AUTORESIZE)
numOrbitsText = wx.StaticText(p, pos=(d,30*5),size=(500,40), label="# Orbits = 34",
              style=wx.ALIGN_CENTRE | wx.ST_NO_AUTORESIZE)
font = wx.Font(30, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
timeText.SetFont(font)
meanAnomaly.SetFont(font)
eccenText.SetFont(font)
semiMajorText.SetFont(font)
speedText.SetFont(font)
numOrbitsText.SetFont(font)

# calculating radio shadow
diff = earth.pos-sat.pos
pointer = arrow(pos=sat.pos,axis=diff, shaftwidth=695700)
def setTimeText(t):
    if t <8760:
        timeText.SetLabel("T+ {:.0f} days".format(t / 24))
    else:
        yrs = t / 8760
        days = (t*8760 - yrs)/24
        timeText.SetLabel("T+ {} yrs, {} days".format(yrs, days))

def showPlanets():
    mercury.pos = mercury_.get_position(t)
    sat.pos = sat_.get_position(t, mercury)
    earth.pos = earth_.get_position(t)
    #venus.pos = venus_.get_position(t)

def showLabels():
    mercury_label.pos=mercury.pos
    earth_label.pos = earth.pos
    sat_label.pos = sat.pos

def showControls():
    setTimeText(t)
    meanAnomaly.SetLabel("True anomaly = {:0.2f} deg".format(mercury_.get_anomaly()))
    speedText.SetLabel("Speed = {:.2f} km/s".format(sat_.get_speed(mercury_)))
    numOrbitsText.SetLabel("# Orbits = {}".format(sat_.get_num_orbits()))

def calculateRadioShadow(mercury, earth):
    # work out straight line between earth and sun y = mx + b
    a = mercury
    b = earth
    m = (b.y - a.y) / (b.x - a.x)
    c = a.y - m*a.x
    # we now know the eqn of line y = mx + c
    # we know eqn of circle x^2+y^2=R^2 where R = radius of sun
    R = 695700
    B = 2 * m * c
    A = m**2 + 1
    C = -R**2 + c**2
    if (B**2 - 4*A*C > 0):
        print("SHADOW")
    else:
        print("NOT SHADOW")
    
    
    
while True:
    rate(60)
    showPlanets()
    showLabels()
    showControls()
    diff = earth.pos-sat.pos
    pointer.pos = sat.pos
    pointer.axis = diff
    calculateRadioShadow(mercury, earth)
    #f1.plot(pos=[t,mercury_.kinetic_energy(sun)])
   # time_label.text = "T +{} days".format(hrsTodays(t))
    t += 1
    
        
