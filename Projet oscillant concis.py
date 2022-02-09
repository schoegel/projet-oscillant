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

Tw = 120.

 

 

Ts = 1./fs

N = int(Tw/Ts)

 

tab_f = np.linspace(0,2.*f0,N)

 

z0 = 0.

zp0 = 10^-9

 

B=-10**19

 

### Conditions initiales

t = 0.

z = 0.

zp= 1e-11

 

tab_t = [t]

tab_z = [z]

tab_zp= [zp]

 

 

###############################################################################

### Defintion des fonctions utilisees pour l'integration de l'eq. diff.

def Force(t,z,Aexc,fexc,A=B): # Fonction de la force totale agissant sur l'oscillateur qui peut dependre de plusieurs parametres sans inté

    Fexc = k0*Aexc*np.cos(2.*np.pi*fexc*t) # force d'excitation harmonique
    
    Fnl = A*z**3 # force d'interaction non-lineaire ici nulle, mais à adapter
    
    F_tot= Fexc+Fnl; # force totale (nulle aussi donc)
    
    return F_tot

 
# Fexc =  0 # force d'excitation harmonique

# Fnl = 0. # force d'interaction non-lineaire ici nulle, mais à adapter

# F_tot= Fexc+Fnl # force totale (nulle aussi donc)

# def Force(t,z,Aexc,fexc):
    
    
    


def E(t,z,zp,Aexc,fexc,A=B): # Fonction E de RK4 qui peut dependre de plusieurs parametres
    return ((w0**2/k0)*Force(t,z,Aexc,fexc,A)-(w0/Q0)*zp-w0**2*z)

 

 

def RK4(t,z,zp,Aexc,fexc,A=B): # Algorithme de RK4 en tant que tel

    k1 = Ts*zp
    
    l1 = Ts*E(t,z,zp,Aexc,fexc,A)
    
     
    
    k2 = Ts*(zp + l1/2)
    
    l2 = Ts*E(t + Ts/2, z+k1/2, zp + l1/2, Aexc, fexc,A)
    
     
    
    k3 = Ts*(zp + l2/2)
    
    l3 = Ts*E(t+Ts/2, z+k2/2, zp+l2/2, Aexc ,fexc,A)
    
     
    
    k4 = Ts*(zp + l3)
    
    l4 = Ts*E(t+Ts, z+k3, zp+l3, Aexc, fexc,A)
    
     
    
    q1 = (k1+2*k2+2*k3+k4)/6
    
    q2 = (l1+2*l2+2*l3+l4)/6
    
     
    
    return([z+q1,zp+q2])

###############################################################################

 

 

###############################################################################

### Simulation de l'oscillation sur N iterations representant une duree de simulation t=Tw

### Cas de la réponse impulsionnelle de l'oscillateur: pas d'excitation et pas de force non-lineaire non plus!!!

Aexc = 0. # Amplitude d'excitation, nulle pour avoir la reponse impulsionnelle de l'oscillateur

fexc = f0 # frequence d'excitation

for i in range(N-1): # Attention int(N) inutile, N, deja entier, suffit
    tab_z.append(RK4(t,tab_z[i],tab_zp[i], Aexc,fexc)[0])
    
    tab_zp.append(RK4(t,tab_z[i],tab_zp[i], Aexc,fexc)[1])
    
    t = t+Ts
    
    tab_t.append(t)

 

 

### Sorties graphiques

plt.figure()

plt.plot(tab_t,tab_z)

plt.xlabel("temps (s)")

plt.ylabel("oscillation z (m)")

plt.title("Reponse impulsionnelle de l'oscillateur en fonction du temps")

plt.grid()

plt.show()

 

plt.figure()

plt.plot(tab_t[1:1000],tab_z[1:1000])

plt.xlabel("temps (s)")

plt.ylabel("oscillation z (m)")

plt.title(" Même graphe que précédemment, zoom temporel")

plt.grid()

plt.show()

###############################################################################

 

### Cas de la réponse harmonique de l'oscillateur: excitation harmonique,

### mais toujours pas de force non-lineaire!!!

# Redefinition des conditions initiales, il faut effacer les var. en memoire

del t, z, zp

del Aexc, fexc

del tab_t, tab_z, tab_zp

 

t = 0.

z = 0.

zp= 1e-9

 

tab_t = [t]

tab_z = [z]

tab_zp= [zp]

 

Aexc = A0/Q0 # Amplitude d'excitation pour avoir une oscillation A0 à la résonance

fexc = f0 # frequence d'excitation egale à la resonance

for i in range(N-1): # Attention int(N) inutile, N, deja entier, suffit
    tab_z.append(RK4(t,tab_z[i],tab_zp[i], Aexc, fexc)[0])
    
    tab_zp.append(RK4(t,tab_z[i],tab_zp[i], Aexc, fexc)[1])

    t = t+Ts

    tab_t.append(t)

 

 

### Sorties graphiques

plt.figure()

plt.plot(tab_t,tab_z)

plt.xlabel("temps (s)")

plt.ylabel("oscillation z (m)")

plt.title("Reponse de l'oscillateur en fonction du temps à une excitation harmonique à f=f0")

plt.grid()

plt.show()

 

plt.figure()

plt.plot(tab_t[1:1000],tab_z[1:1000])

plt.xlabel("temps (s)")

plt.ylabel("oscillation z (m)")

plt.title(" Même graphe que précédemment, zoom temporel")

plt.grid()

 

plt.show()

###############################################################################

 

### Il faut maintenant reflechir à la question du balayage en frequence,

### puis au cas de l'interaction non-lineaire!

### A VOUS DE JOUER

 

###############################################################################




 

del t, z, zp

del Aexc, fexc

del tab_t, tab_z, tab_zp

 

t = 0.

z = 0.

zp= 1e-9

 

tab_t1 = [t]

tab_z1 = [z]

tab_zp1= [zp]



tab_t2 = [t]

tab_z2 = [z]

tab_zp2= [zp]

 

Aexc = A0/Q0 # Amplitude d'excitation pour avoir une oscillation A0 à la résonance

fexc = np.linspace(0,2*f0,N) # frequence d'excitation allant de 0 a 2f0 en N pts

 
for i in range(N-1):
    
    tab_z1.append(RK4(t,tab_z1[i],tab_zp1[i], Aexc, fexc[i]/2,0)[0])
    
    tab_zp1.append(RK4(t,tab_z1[i],tab_zp1[i], Aexc, fexc[i]/2,0)[1])

    t = t+Ts
    
    tab_t1.append(t)


for i in range(N-1):
    
    tab_z2.append(RK4(t,tab_z2[i],tab_zp2[i], Aexc, fexc[i]/2)[0])
    
    tab_zp2.append(RK4(t,tab_z2[i],tab_zp2[i], Aexc, fexc[i]/2)[1])

    t = t+Ts
    
    tab_t2.append(t)

 

tab_sommets=[]

freq_sommets=[]

 

for i in range(1,len(tab_zp1)): # On utilise l'abcisse pour laquelle la derivee zp s'annule,

    if tab_zp1[i]*tab_zp1[i-1]<=0: # ce qui nous donne la coordonnee de chaque sommet de z,

        tab_sommets.append(tab_z1[i]/A0) # on peut donc en déduire son amplitude
        
        freq_sommets.append(fexc[i]/f0)

 

tab_amp=[]


 

for i in range(1,len(tab_sommets)):
    
    tab_amp.append((abs(tab_sommets[i])+abs(tab_sommets[i-1]))/2)

 
z1_normalise = [z/A0 for z in tab_z1]

f_normalise = [f/f0 for f in fexc]

z2_normalise = [z/A0 for z in tab_z2]

 

plt.figure()

plt.title('Oscillations et amplitude en fonction de f')

#plt.plot(freq_sommets[0:-1],tab_amp, '-o', color='red')

plt.plot(f_normalise,tab_z1)

plt.plot(f_normalise,tab_z2)

plt.xlabel('Frequence d excitation fexc')

plt.ylabel('z(t) et amplitude d oscillation')

plt.grid()

plt.show()


