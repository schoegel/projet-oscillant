#!/usr/bin/env python3
#variables

f0 = 1000 #s-1

q0 = 100 

k0 = 20 #N/m

#w0 = 2pif0 = sqrt(k0/m0)

#q0 = m0*w0/gam0

fs = 40*f0 

Ts = 1/fs #s

Tw = 10 #s

#N = Tw/Ts = 400000

df0 = 30 #s-1


#trace temporel de fréquence

tab_t = np.arange(0,10,Ts) 

print(tab_t)

tab_f = np.arange(0,2000,Ts)

print(tab_f)

plt.plot(tab_x_rungekuttascipy,tab_y_rungekuttascipy, label = "rungekutta scipy" )      
plt.legend()
plt.xlabel("x")
plt.ylabel("y(x)") 
plt.title("rungekuttascipy")
plt.savefig("rungekuttascipy.pdf")       
plt.show()

