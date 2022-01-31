# projet-oscillant

from scipy.integrate import ode

import numpy as np

import matplotlib.pylab as plt

 

 

### Variables

 

f0 = 1.


w0 = 2*np.pi*f0


Q0 = 10.


k0 = 20.

fs = 40. * f0 

Tw = 120.

Ts = 1/fs

N = int(Tw/Ts)

t=0

tab_t = np.linspace(0,Tw,N)

tab_f = np.linspace(0,2*f0,len(tab_t))

z0 = 0.

zp0 = 10^-9

B=10

### Energie

 

z=0.

zp=1e-9

tab_z=[z]

tab_zp=[zp]

 
Fexc= B*z**3
 

def E(t,z,zp):
    return ((w0**2/k0)*Fexc-(w0/Q0)*zp-w0**2*z)

 

def RK4(t,z,zp):

    k1 = Ts*zp
    
    l1 = Ts*E(t,z,zp)
    
    k2 = Ts*(zp + l1/2)
    
    l2 = Ts*E(t + Ts/2, z+k1/2, zp + l1/2)
    
    k3 = Ts*(zp + l2/2)
    
    l3 = Ts*E(t+Ts/2, z+k2/2, zp+l2/2)
    
    k4 = Ts*(zp + l3)
    
    l4 = Ts*E(t+Ts, z+k3, zp+l3)
    
    q1 = (k1+2*k2+2*k3+k4)/6
    
    q2 = (l1+2*l2+2*l3+l4)/6
    
    return([z+q1,zp+q2])

 

 

for i in range(int(N)-1):

    tab_z.append(RK4(t,tab_z[i],tab_zp[i])[0])
    
    tab_zp.append(RK4(t,tab_z[i],tab_zp[i])[1])
    
    np.append(t,Ts)

 
plt.title("avec force d'exitation")
plt.plot(tab_t,tab_z)
plt.title("sans force d'exitation")
plt.show()

#force d'exitation

z=0.

zp=1e-9

tab_z=[z]

tab_zp=[zp]


Aexc = 100

B = 100

Fexc = k0*Aexc*np.cos(2*np.pi*f0*t) + B*z**3

for i in range(int(N)-1):

    tab_z.append(RK4(t,tab_z[i],tab_zp[i])[0])
    
    tab_zp.append(RK4(t,tab_z[i],tab_zp[i])[1])
    
    np.append(t,Ts)
    
plt.plot(tab_t,tab_z)
plt.title("avec force d'exitation")
plt.show()
    
    
