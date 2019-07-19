import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

N = np.arange(1,1000)

i = 3

def f(x):
    return x**2 + x*4

exact_value = 18

def fi(t):
    return i*np.tanh((np.pi/2)*np.sinh(t))


def fi_diff(t):
    return i*(((np.pi/2)*np.cosh(t))/(np.cosh((np.pi/2)*np.sinh(t)))**2)
"""
def fi(t, h):
    b = 1/4
    M = np.pi/h
    a = b/np.sqrt(1+M*(np.log(1+M))/(1/(4*np.pi)))
    a2 = np.e**(-2*t-a*(1-np.e**(-t))-b*(np.e**(t)-1))

    return t/(1-a2)


def fi_diff(t,h):
    eps = 1e-100
    diff = (fi(t+eps,h)-fi(t-eps,h))/(2*eps)
    return diff
"""

def DE_formula(t, h):
    return h*(f(fi(t))*fi_diff(t)).sum()

def normal_integral(x, h):
    sum_ = np.array([])
    for i in range(len(x)-1):
        tmpx = (x[i+1]+x[i])/2
        sum_ = np.append(sum_,(h*f(tmpx)))
    return sum_.sum()

error_DEformul_n = np.array([])
error_norlmalmethod_n = np.array([])

for n in N:
    #print(h)
    #x = np.arange(-i,i,h)
    h = 1/100
    t = np.arange(-n,n,h)
    error_DEformul_n = np.append(error_DEformul_n, (DE_formula(t, h)-exact_value)**2)
    #if abs(DE_formula(t, h)-exact_value)<= 1e-10:
    #    break

    #error_norlmalmethod = np.append(error_norlmalmethod, (normal_integral(x, h)-exact_value)**2)

n = 3

error_DEformul = np.array([])
error_norlmalmethod = np.array([])

for h in 1/N:
    #print(h)
    x = np.arange(-i,i,h)
    t = np.arange(-n,n,h)
    M = np.pi/h
    error_DEformul = np.append(error_DEformul, (DE_formula(t, h)-exact_value)**2)
    error_norlmalmethod = np.append(error_norlmalmethod, (normal_integral(x, h)-exact_value)**2)

import matplotlib.pyplot as plt

#plt.plot(N, error_DEformul, label='DE_formula')
#plt.plot(N, error_norlmalmethod, label='normal_method')
#plt.plot(list([i+1 for i in range(100)]), error_norlmalmethod, label='normal_method')

plt.plot(N, error_DEformul_n)
plt.xscale('log')
plt.yscale('log')
plt.title('x^2+4x    [{0}:{1}]'.format(-n,n))
#plt.xlabel('N')
plt.xlabel('Threshold  N')
plt.ylabel('Error')
plt.legend()
plt.show()

plt.plot(N, error_DEformul, label='DE_formula')
plt.plot(N, error_norlmalmethod, label='normal_method')
plt.xscale('log')
plt.yscale('log')
plt.title('x^2+4x    [-3:3]   n={}'.format(n))
plt.xlabel('h')
#plt.xlabel('Threshold')
plt.ylabel('Error')
plt.legend()
plt.show()
