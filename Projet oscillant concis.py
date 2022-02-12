from scipy.integrate import ode

import numpy as np

import matplotlib.pylab as plt


###############################################################################

### Variables

f0 = 10

w0 = 2.*np.pi*f0

Q0 = 10.

k0 = 20.

A0 = 1e-9

fs = 40. * f0

Tw = 120

Ts = 1./fs

N = int(Tw/Ts)

tab_f = np.linspace(0,2.*f0,N)

tab_t = np.linspace(0,Tw,N)
 
B=-10**19
 
###############################################################################

### Defintion des fonctions utilisees pour l'integration de l'eq. diff et cacluls des osccilations

def Force(t,z,Aexc,fexc,B): # Fonction de la force totale agissant sur l'oscillateur qui peut dependre de plusieurs parametres sans inté

    Fexc = k0*Aexc*np.cos(2.*np.pi*fexc*t) # force d'excitation harmonique
    
    Fnl = B*z**3 # force d'interaction non-lineaire ici nulle, mais à adapter
    
    F_tot= Fexc+Fnl; # force totale (nulle aussi donc)
    
    return F_tot
    
    

def E(t,z,zp,Aexc,fexc,B): # Fonction E de RK4 qui peut dependre de plusieurs parametres
    return ((w0**2/k0)*Force(t,z,Aexc,fexc,B)-(w0/Q0)*zp-w0**2*z)

 
 

def RK4(t,z,zp,Aexc,fexc,B):
    k1 = Ts*zp
    l1 = Ts*E(t,z,zp,Aexc,fexc,B)
    
    k2 = Ts*(zp + l1/2)
    l2 = Ts*E(t + Ts/2, z+k1/2, zp + l1/2, Aexc, fexc,B)
     
    k3 = Ts*(zp + l2/2)
    l3 = Ts*E(t+Ts/2, z+k2/2, zp+l2/2, Aexc ,fexc,B)
    
    k4 = Ts*(zp + l3)
    l4 = Ts*E(t+Ts, z+k3, zp+l3, Aexc, fexc,B)
    
    q1 = (k1+2*k2+2*k3+k4)/6
    q2 = (l1+2*l2+2*l3+l4)/6
    
    return([z+q1,zp+q2])


def Oscillations(Aexc,fexc,B):
    t = 0.
    z = 0.    
    zp= 1e-11
    tab_t = [t]
    tab_z = [z]
    tab_zp= [zp]
    if type(fexc)!=np.ndarray:
        for i in range(N-1):
            tab_z.append(RK4(t,tab_z[i],tab_zp[i], Aexc,fexc,B)[0])
            tab_zp.append(RK4(t,tab_z[i],tab_zp[i], Aexc,fexc,B)[1])
            t = t+Ts
            tab_t.append(t)
    else:
        for i in range(N-1):
            tab_z.append(RK4(t,tab_z[i],tab_zp[i], Aexc, fexc[i]/2,B)[0])
            tab_zp.append(RK4(t,tab_z[i],tab_zp[i], Aexc, fexc[i]/2,B)[1])
            t = t+Ts
            tab_t.append(t)
    return(tab_z,tab_zp)


def amp(B,f):
    z=Oscillations(A0/Q0,f,B)[0]
    zp=Oscillations(A0/Q0,f,B)[1]
    z_sommets=[]
    f_sommets=[]
    f_amp=[]
    for i in range(len(zp)):
        if zp[i]*zp[i-1]<=0:
            z_sommets.append(z[i])
            f_sommets.append(f[i])
    amp=[z_sommets[2*i+1] for i in range(int(len(z_sommets)/2))]
    f_amp=[f_sommets[2*i+1] for i in range(int(len(f_sommets)/2))]
    return(amp,f_amp,z_sommets)

########## Plots

### Oscillations sans excitation

# z_1=Oscillations(0,0,0)[0] 

# plt.plot(tab_t,z_1)
# plt.show()

# plt.plot(tab_t[0:1000],z_1[0:1000])
# plt.show()

### Oscillations à f0

# z_exc = Oscillations(A0/Q0,f0,0)[0]

# plt.plot(tab_t,z_exc)
# plt.show()

# plt.plot(tab_t[0:1000],z_exc[0:1000])
# plt.show()

### Oscillations de 0 à 2*f0

# z_freq=Oscillations(A0/Q0,tab_f,0)[0]

# plt.plot(tab_f,z_freq)
# plt.show()

# z_freq_nl = Oscillations(A0/Q0,tab_f,B)[0]

# plt.plot(tab_f,z_freq_nl)
# plt.show()

# ### Amplitude de 0 à 2*f0

# plt.plot(amp(0,tab_f)[1],amp(0,tab_f)[0])
# plt.plot(amp(B,tab_f)[1],amp(B,tab_f)[0])
# plt.show()


f_down=np.linspace(2*f0,0,N)


plt.plot(tab_t,tab_f)
plt.plot(tab_t,f_down)
plt.show()

plt.plot(tab_t,Oscillations(A0/Q0,f0,B)[0])
plt.show()






















