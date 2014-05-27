import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint

"""

general model

max lifespan = ml
ml is upper levelof 95% C.I.
rate of change of lifespan curve = sl
sl is different in time at 90% mortality and 50% mortality

delg = -R*T(ln([S]*((max lifespan/rate of change of lifespan curve)-1)))
"""
"""
N2_text = open('3.7_gibbs-N2.txt','r')
eri6_text = open('3.7_gibbs-eri6(mg379).txt','r')

N2 = N2_text.read()
eri6 = eri6_text.read()
N2_text.close()
eri6_text.close()
"""
# eri6 gibbs constants

# 25ul/50ml tBuOOH

e50ml=13.98
e50sl=23-17

# 25ul/25ml tBuOOH

e25ml=17.84
e25sl=23-17

# 25ul/10ml tBuOOH

e10ml=21.25
e10sl=33-17

# n2 gibbs constants

# 25ul/50ml tBuOOH

n50ml=16.13
n50sl=26-17

# 25ul/25ml tBuOOH

n25ml=13.63
n25sl=26-17

# 25ul/10ml tBuOOH

n10ml=11.97
n10sl=37-6

# Constants in model
R= 8.314
#T= 298.15
T=np.linspace(0,10000,2000) 
ten = 25e-6/10e-3
twofive = 25e-6/25e-3
fifty = 25e-6/50e-3

#dG = -R*T(ln(S*((max lifespan/rate of change of lifespan curve)-1)))

e50=-R*T*np.log((fifty*e50ml)/(e50sl-1))
e25=-R*T*np.log((twofive*e25ml)/(e50sl-1))
e10=-R*T*np.log((ten*e10ml)/(e50sl-1))

n50=-R*T*np.log((fifty*n50ml)/(n50sl-1))
n25=-R*T*np.log((twofive*n25ml)/(n25sl-1))
n10=-R*T*np.log((ten*n10ml)/(n10sl-1))



# Plot the solutions

plt.figure()
p0,=plt.plot(e50) 
p1,=plt.plot(e25) 
p2,=plt.plot(e10) 
p3,=plt.plot(n50) 
p4,=plt.plot(n25) 
p5,=plt.plot(n10)
plt.legend([p0,p1,p2,p3,p4,p5],["eri-6 (25/50)","eri-6 (25/25)","eri-6 (25/10)","N2 (25/50)","N2 (25/25)","N2 (25/10)"])
plt.xlabel('Temperature (Kelvin)') 
plt.ylabel('Delta G (Joules/mol)')
plt.show()
