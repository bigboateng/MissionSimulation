from  visual import *

G = 6.7e-11 #N.m^2/kg2
Re = 6.37e6 # Radius of earth
Rm = 1.74e6 # Radius of moon
Rs = 6.95e8 # Radius of sun
t =0
dt = 1000 #time step
Se = 1.5e11 # Distance from sun to eath
Sm = 3.48e8
Ms = 1.989e30
Me = 5.97e24
Mm = 7.347e22


shouldRun= True

earth=sphere(pos=vector(0,0,0), radius=Re, material=materials.earth, make_trail=True)
moon = sphere(pos=earth.pos + vector(-Sm,0,0), radius=Rm, color=(.7,0.7,.7),make_trail=True)
sun = vector(-Se,0,0)
earth.m = Me
moon.m = Mm

## Initial Angular velocities
we = sqrt(G*Ms/Se**3)
wm = sqrt(G*Me/Sm**3)
tmonth=28*24*60*20

earth.p = vector(0, Se*Me*we, 0)
moon.p = moon.m *(earth.p/earth.m+vector(0,-Sm*wm,0))

while t < tmonth:
    rate(100)
    re=earth.pos-sun
    rem=moon.pos-earth.pos
    res=moon.pos-sun
    Fe = G*Ms*earth.m*norm(re)/mag(re)**2+G*earth.m*moon.m*norm(rem)/mag(rem)**2
    Fm = -G*Ms*moon.m*norm(res)/mag(res)**2-G*earth.m*moon.m*norm(rem)/mag(rem)**2
    earth.p = earth.p+Fe*dt
    moon.p = moon.p+Fm*dt
    earth.pos += earth.p*dt/earth.m
    moon.pos+=moon.p*dt/moon.m
    t+=dt
