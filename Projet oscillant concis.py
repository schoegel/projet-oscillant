from scipy.integrate import ode

import numpy as np

import matplotlib.pylab as plt


###############################################################################

### Variables

f0 = 10

w0 = 2.*np.pi*f0

Q0 = 5

k0 = 20.

A0 = 1e-9

fs = 40. * f0

Tw = 120

Ts = 1./fs

N = int(Tw/Ts)

tab_f = np.linspace(0,2*f0,N)

f_down = np.linspace(2*f0,0,N)

tab_f_normalise = [f/f0 for f in np.linspace(0,2*f0,N)]

f_down_normalise = [f/f0 for f in np.linspace(2*f0,0,N)]

tab_t = np.linspace(0,Tw,N)

Beta = -10**19         # Coefficient Beta de l'interaction non lineaire de l'effet Duffing

d = 3*10**-10       # Distance caractéristique de l'intéraction de Van der Waals
 
Alpha = 3*10**-28       # Coefficient Alpha de l'interaction de Van der Waals
 
###############################################################################

### Definition des fonctions utilisees pour l'integration de l'eq. diff, calcul des oscillations, amplitude, et phase

def Force(t, z, Aexc, fexc, Beta, sens='croissant', Alpha=0):   # Fonction de la force totale agissant sur l'oscillateur qui peut dependre de plusieurs 
                                                                # parametres
    if sens == 'croissant':                                     
        Fexc = k0*Aexc*np.cos(2.*np.pi*fexc*t) # force d'interaction harmonique
        
    if sens == 'decroissant':
        Fexc = k0*Aexc*np.cos(2.*np.pi*fexc*(Tw-t))
        
    Fnl = Beta*z**3 ##= force d'interaction non lineaire correspondant à l'effet Duffing
    
    Fnl2 = -Alpha/((d+A0+z)**2) # force d'interaction non lineaire de Van der Waals
    
    F_tot = Fexc + Fnl + Fnl2
    
    return F_tot
    
    

def E(t, z, zp, Aexc, fexc, Beta, sens='croissant',Alpha=0): # E est l'expression de z" en fonction des autres variables du système
    return ((w0**2/k0)*Force(t,z,Aexc,fexc,Beta,sens,Alpha)-(w0/Q0)*zp-w0**2*z)

 
 

def RK4(t,z,zp,Aexc,fexc,Beta, sens='croissant',Alpha=0): # Algorithme de RK4
    k1 = Ts*zp
    l1 = Ts*E(t,z,zp,Aexc,fexc,Beta,sens,Alpha)
    
    k2 = Ts*(zp + l1/2)
    l2 = Ts*E(t + Ts/2, z+k1/2, zp + l1/2, Aexc, fexc,Beta,sens,Alpha)
     
    k3 = Ts*(zp + l2/2)
    l3 = Ts*E(t+Ts/2, z+k2/2, zp+l2/2, Aexc ,fexc,Beta,sens,Alpha)
    
    k4 = Ts*(zp + l3)
    l4 = Ts*E(t+Ts, z+k3, zp+l3, Aexc, fexc,Beta,sens,Alpha)
    
    q1 = (k1+2*k2+2*k3+k4)/6
    q2 = (l1+2*l2+2*l3+l4)/6
    
    return([z+q1,zp+q2])


def Oscillations(Aexc,fexc,Beta, sens='croissant',Alpha=0):     # Cette fonction nous permet d'automatiser l'intégration des oscillations 
    t = 0.                                                      #  à partir d'une excitation, d'une force d'interaction non lineaire
    z = 0.                                                      #  et du sens de balayage
    zp= 1e-11
    tab_t = [t]
    tab_z = [z]
    tab_zp= [zp]
    if type(fexc)!=np.ndarray:              # On teste si la frequence est une seule valeur (frequence d'excitation constante)
        for i in range(N-1):                #  ou si la frequence est sous forme de tableau et donc variable (wobulation)
            tab_z.append(RK4(t,tab_z[i],tab_zp[i], Aexc,fexc,Beta,sens,Alpha)[0])
            tab_zp.append(RK4(t,tab_z[i],tab_zp[i], Aexc,fexc,Beta,sens,Alpha)[1])
            t = t+Ts
            tab_t.append(t)
    else:
        for i in range(N-1):
            tab_z.append(RK4(t,tab_z[i],tab_zp[i], Aexc, fexc[i]/2,Beta,sens,Alpha)[0])
            tab_zp.append(RK4(t,tab_z[i],tab_zp[i], Aexc, fexc[i]/2,Beta,sens,Alpha)[1])
            t = t+Ts
            tab_t.append(t)
    return(tab_z,tab_zp)


def amp(Beta,f, sens='croissant',Alpha=0):
    z=Oscillations(A0/Q0,f,Beta,sens,Alpha)[0]
    zp=Oscillations(A0/Q0,f,Beta,sens,Alpha)[1]
    f_amp=[]
    amp=[]
    for i in range(len(zp)):
        if zp[i]*zp[i-1]<=0 and zp[i]<=zp[i-1]:         # On identifie l'abcisse de chaque zeros (ou changement de signe) de la dérivée de z, 
            amp.append(max(z[i],z[i-1]))                # ce qui permet d'identifier l'abcisse de chaque sommet et donc d'obtenir l'amplitude
            f_amp.append(f[i])
    return(amp,f_amp)
    
    
### Calculs de différentes oscillations et amplitudes
    
z_sans_excitation = Oscillations(0,0,0)[0]             
z_avec_excitation = Oscillations(A0/Q0,f0,0)[0]     
    

z_lineaire = Oscillations(A0/Q0,tab_f,0)[0]          
zp_lineaire = Oscillations(A0/Q0,tab_f,0)[1]
z_amp_lineaire = amp(0,tab_f)[0]                       
f_amp_lineaire = amp(0,tab_f)[1]                      


z_non_lineaire_1 = Oscillations(A0/Q0,tab_f,B)[0]
z_amp_non_lineaire_1 = amp(Beta,tab_f)[0]
f_amp_non_lineaire_1 = amp(Beta,tab_f)[1]


z_non_lineaire_2 = Oscillations(A0/Q0,tab_f,-B)[0]
z_amp_non_lineaire_2 = amp(-Beta,tab_f)[0]
f_amp_non_lineaire_2 = amp(-Beta,tab_f)[1]


z_amp_non_lineaire_1_down = amp(Beta,f_down,sens='decroissant')[0]
f_amp_non_lineaire_1_down = amp(Beta,f_down,sens='decroissant')[1]


z_non_lineaire_3 = Oscillations(A0/Q0,tab_f,0,Alpha=A)[0]
z_amp_non_lineaire_3 = amp(0,tab_f,Alpha=A)[0]
f_amp_non_lineaire_3 = amp(0,tab_f,Alpha=A)[1]


z_amp_non_lineaire_3_down = amp(0,f_down,sens='decroissant',Alpha=A)[0]
f_amp_non_lineaire_3_down = amp(0,f_down,sens='decroissant',Alpha=A)[1]
    
########## Plots


### Oscillations linéaires sans excitation et avec excitation constante

plt.figure(1)

plt.subplot(2,2,1)
plt.plot(tab_t,z_sans_excitation)
plt.title('Oscillations sans excitation')
plt.grid()

plt.subplot(2,2,2)
plt.plot(tab_t[0:1000],z_sans_excitation[0:1000])
plt.title('Oscillations sans excitation zoomé')
plt.grid()

plt.subplot(2,2,3)
plt.plot(tab_t,z_avec_excitation)
plt.title('Oscillations avec excitation constante f0')
plt.grid()

plt.subplot(2,2,4)
plt.plot(tab_t[0:1000],z_avec_excitation[0:1000])
plt.title('Oscillations avec excitation constante f0 zoomée')
plt.grid()

plt.show()


### Oscillations, amplitude et hystérésis avec balayage en fréquence pour force linéaire et non linéaire en -B*z^3

plt.figure(2)

plt.subplot(2,2,1)
plt.plot(tab_f,z_lineaire)
plt.plot(f_amp_lineaire,z_amp_lineaire)
plt.title('Reponse impulsionnelle Beta=0')
plt.xlabel('frequence')
plt.ylabel('oscillation')
plt.grid()
plt.show()

plt.subplot(2,2,2)
plt.plot(tab_f,z_non_lineaire_1)
plt.plot(f_amp_non_lineaire_1, z_amp_non_lineaire_1)
plt.title('Reponse impulsionnelle B=-10e19')
plt.xlabel('frequence')
plt.ylabel('oscillation')
plt.grid()

plt.subplot(2,2,3)
plt.plot(tab_f,z_non_lineaire_2)
plt.plot(f_amp_non_lineaire_2, z_amp_non_lineaire_2)
plt.title('Reponse impulsionnelle B=10e19')
plt.xlabel('frequence')
plt.ylabel('oscillation')
plt.grid()

plt.subplot(2,2,4)
plt.plot(f_amp_non_lineaire_1_down, z_amp_non_lineaire_1_down,label='Amplitude avec f decroissante')
plt.plot(f_amp_non_lineaire_1,z_amp_non_lineaire_1, label='Amplitude avec f croissante')
plt.title('Hysteresis avec B=-10e19')
plt.xlabel('frequence')
plt.ylabel('amplitude')
plt.grid()

plt.show()


### Oscillations, amplitude et hystérésis avec balayage en fréquence pour force non linéaire en alpha/(d+A0+z)^2

plt.figure(3)

plt.subplot(2,1,1)
plt.plot(tab_f,z_non_lineaire_3)
plt.plot(f_amp_non_lineaire_3,z_amp_non_lineaire_3)
plt.title('Réponse impulsionnelle A=3e-28')
plt.xlabel('frequence')
plt.ylabel('amplitude')
plt.grid()

plt.subplot(2,1,2)
plt.plot(f_amp_non_lineaire_3_down,z_amp_non_lineaire_3_down)
plt.plot(f_amp_non_lineaire_3,z_amp_non_lineaire_3)
plt.title('Hysteresis avec A=3e-28')
plt.xlabel('frequence')
plt.ylabel('amplitude')
plt.grid()

plt.show()


### Amplitude et phase par la méthode de la détection synchrone

p=600

def Lockin(z, precision=p):
    X=[]
    Y=[]
    X_moy=[]
    Y_moy=[]
    for i in range(int(N/precision)):   # On met bout a bout des sous-listes de largeur p correspondant au découpage de z
        X.append([z[k]*np.cos(2*np.pi*tab_f[k]/2*tab_t[k]) for k in range(precision*i,precision*(i+1))])
        Y.append([z[k]*np.sin(2*np.pi*tab_f[k]/2*tab_t[k]) for k in range(precision*i,precision*(i+1))])
    for k in range(int(N/precision)):
        X_moy.append(np.average(X[k]))
        Y_moy.append(np.average(Y[k]))
    Amplitude = [2*np.sqrt(X_moy[i]**2 + Y_moy[i]**2) for i in range(int(N/precision))]
    Phase = [np.arctan(Y_moy[i]/X_moy[i]) for i in range(int(N/precision))]
    return([Amplitude,Phase])


Ampli_lockin = Lockin(z_lineaire)[0]
Phase_lockin = Lockin(z_lineaire)[1]
f_lockin = np.linspace(0,2*f0,int(N/p))


### Amplitude et phase théorique

def Ampli_phase_theorique(f):
    Ampli = [A0/(np.sqrt(Q0**2*(1-(i/f0)**2)**2 + (i/f0)**2)) for i in f]
    Phase = [np.arctan(i/(f0*Q0*(1-(i/f0)**2))) for i in f]
    return(Ampli,Phase)

Ampli_theorique = Ampli_phase_theorique(tab_f)[0]
Phase_theorique = Ampli_phase_theorique(tab_f)[1]


plt.figure(4)

plt.subplot(2,2,1)
plt.plot(f_lockin, Ampli_lockin)
plt.title('Amplitude et phase numerique')
plt.grid()

plt.subplot(2,2,2)
plt.plot(tab_f,Ampli_theorique)
plt.title('Amplitude et phase theorique')
plt.grid()

plt.subplot(2,2,3)
plt.plot(f_lockin,Phase_lockin)
plt.grid()

plt.subplot(2,2,4)
plt.plot(tab_f,Phase_theorique)
plt.grid()

plt.show()


### Marge d'erreur et precision du modele

# Comparaison de la détection synchrone et modèle théorique

pourcentage = 1
erreur_ampli = []
C = 0

for i in range(len(Ampli_lockin)):
    erreur_ampli.append(abs(Ampli_lockin[i] - A0/(np.sqrt(Q0**2*(1-(f_lockin[i]/f0)**2)**2 + (f_lockin[i]/f0)**2))))
    
for j in range(len(erreur_ampli)):
    if erreur_ampli[j] <= (pourcentage/100)*A0/(np.sqrt(Q0**2*(1-(f_lockin[i]/f0)**2)**2 + (f_lockin[i]/f0)**2)):
        C+=1

ratio = C/len(erreur_ampli)*100

print(ratio, '% des points par détection synchrone sont à', pourcentage, '% ou moins de la valeur théorique' )

enveloppe_1 = [i*(1+pourcentage/100) for i in Ampli_theorique]
enveloppe_2 = [i*(1-pourcentage/100) for i in Ampli_theorique]


var = np.sqrt(np.average([i**2 for i in erreur_ampli]))

# Comparaison de la méthode des dérivées et modèle théorique

pourcentage_2 = 1

erreur_ampli_2 = []

C_2 = 0

for i in range(len(z_amp_lineaire)):
    erreur_ampli_2.append(abs(z_amp_lineaire[i] - A0/(np.sqrt(Q0**2*(1-(f_amp_lineaire[i]/f0)**2)**2 + (f_amp_lineaire[i]/f0)**2))))
    
for j in range(len(erreur_ampli_2)):
    if erreur_ampli_2[j] <= (pourcentage_2/100)*A0/(np.sqrt(Q0**2*(1-(f_amp_lineaire[i]/f0)**2)**2 + (f_amp_lineaire[i]/f0)**2)):
        C_2+=1

ratio_2 = C_2/len(erreur_ampli_2)*100

print(ratio_2, '% des points par la méthode des dérivées sont à', pourcentage_2, '% ou moins de la valeur théorique' )

var_2 = np.sqrt(np.average([i**2 for i in erreur_ampli_2]))


# plots des différents modèles et comparaison

plt.figure(5)


plt.subplot(1,2,1)
plt.plot(tab_f,enveloppe_1,color = 'k', label = 'Enveloppe à 1% autour de la valeur théorique')
plt.plot(tab_f,enveloppe_2, color = 'k')
plt.plot(f_lockin,Ampli_lockin, color = 'orangered', label = 'Amplitude par détection synchrone')
plt.plot(f_amp_lineaire,z_amp_lineaire, color = 'steelblue', label = 'Amplitude par méthode des dérivées')
plt.xlabel('Fréquence en Hz')
plt.ylabel('Amplitude en m')
plt.legend()
plt.grid()
plt.title('Comparaison des différents modèles à la valeur théorique')


plt.subplot(3,2,2)
plt.xlim(1, 2)
plt.ylim(1.92e-10, 2.15e-10)
plt.plot(tab_f,enveloppe_1,color = 'k')
plt.plot(tab_f,enveloppe_2, color = 'k')
plt.plot(f_lockin,Ampli_lockin, color = 'orangered')
plt.plot(f_amp_lineaire,z_amp_lineaire, color = 'steelblue')
plt.grid()


plt.subplot(3,2,4)
plt.xlim(9.85, 10)
plt.ylim(0.99e-9, 1.02e-9)
plt.plot(tab_f,enveloppe_1,color = 'k')
plt.plot(tab_f,enveloppe_2, color = 'k')
plt.plot(f_lockin,Ampli_lockin, color = 'orangered')
plt.plot(f_amp_lineaire,z_amp_lineaire, color = 'steelblue')
plt.grid()


plt.subplot(3,2,6)
plt.xlim(19.5, 20)
plt.ylim(0.65e-10, 0.725e-10)
plt.plot(tab_f,enveloppe_1,color = 'k')
plt.plot(tab_f,enveloppe_2, color = 'k')
plt.plot(f_lockin,Ampli_lockin, color = 'orangered')
plt.plot(f_amp_lineaire,z_amp_lineaire, color = 'steelblue')
plt.grid()

plt.show()






































