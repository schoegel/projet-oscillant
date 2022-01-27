from scipy.integrate import ode
import numpy as np
import matplotlib.pylab as plt



#variables

f0 = 1000 #s-1

q0 = 10 

k0 = 20 #N/m

#w0 = 2pif0 = sqrt(k0/m0)

#q0 = m0*w0/gam0

fs = 40*f0 

Ts = 1/fs #s

Tw = 120 #s

#N = Tw/Ts = 400000

df0 = 30 #s-1


#trace temporel de fréquence

tab_t = np.arange(0,10,Ts) 

print(tab_t)

tab_f = np.linspace(0,2000,len(tab_t))

print(tab_f)

plt.plot(tab_t,tab_f)      
plt.legend()
plt.xlabel("t")
plt.ylabel("f(t)") 
plt.title("trace temporelle de la fréquence")
plt.savefig("trace temporelle de la fréquence.pdf")       
plt.show()
