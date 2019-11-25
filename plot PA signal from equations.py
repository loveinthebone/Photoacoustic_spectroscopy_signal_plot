
#According to paper"Wavelength modulation spectroscopy: combined frequency and intensity laser modulation"
from numpy import sqrt,cos,sin,abs
import numpy as np
import matplotlib.pyplot as plt

Fi=-np.pi/2 #phase shift between the intensity and wavelength modulations
p_row=-1
p_w=-2
delta_v=0.07327 #absorption line bandwidth
m=2 #modulation index

#eq A10 A3
L=12e-3 #path length
N=100e-6 #density in ppm
S=1e-20 #line strength
I0=24e-3 #ex laser power in mW
a0=L*N*S/(3.14*delta_v)

x=np.arange(-6,6,0.1)

I_O=p_row*delta_v*x+1 #eq(10)

X=1-x**2+m**2 #eq A12
r=(X**2+4*x**2)**0.5

s0x=I0*(1-a0*((2*(r+X))**0.5)/(2*r))
s1x=I0*a0*(sqrt(2)*(-x*sqrt(r+X)+np.sign(x)*sqrt(r-X))/(m*r))

s2x=I0*a0*(-4/m**2+sqrt(2)/m**2*((r+1-x**2)*sqrt(r+X)+2*abs(x)*sqrt(r-X))/r)

s3x=-I0*a0/m**3*(16*x+sqrt(2)/r*(x**3-3*x*(r+1))*sqrt(r+X)+sqrt(2)/r*np.sign(x)*(1-3*x**2-3*r)*sqrt(r-X))


s1p=I_O*cos(Fi)*s1x-p_w*delta_v*m/2*(2*s0x+cos(2*Fi)*s2x)
s1q=I_O*sin(Fi)*s1x-p_w*delta_v*m/2*sin(2*Fi)*s2x



s2p=I_O*cos(2*Fi)*s2x-p_w*delta_v*m/2*(cos(Fi)*s1x+cos(3*Fi)*s3x)
s2q=I_O*sin(2*Fi)*s2x-p_w*delta_v*m/2*(sin(Fi)*s1x+sin(3*Fi)*s3x)


plt.plot(x,s2p)
plt.show()

plt.plot(x,s2q)
plt.show()
