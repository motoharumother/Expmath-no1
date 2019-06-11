import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
x = np.arange(-pi, pi, pi/1000)

def bk(k):
    tmpb = (4/k**2)*np.cos(k*pi)
    return tmpb

def Sn(n,b,x):
    y = np.zeros((n))
    sin_list = np.array([np.cos((k+1)*x) for k in range(n)])
    sn = sin_list*b
    return sn.sum()


nlist = np.arange(1,100,32)

maxlist = np.zeros((len(nlist)))

y = np.zeros((len(nlist), len(x)))

plt.figure(figsize=(10,7))

for m, n in enumerate(nlist):
    b = np.zeros((n))
    for k in range(n):
        b[k] = bk(k+1)

    for i, _x in enumerate(x):
        y[m][i] = Sn(n,b,_x)
    plt.plot(x, y[m], label='n={}'.format(n), alpha=0.4)
    
er = [((x-y[i])**2).sum() for i in range(len(nlist))]
superror = [(abs(x-y[i])).max() for i in range(len(nlist))]
    
plt.plot(x, x**2, label='orignal')
plt.legend()
plt.show()

plt.plot(nlist, er)
plt.scatter(nlist, er, color='r')
plt.show()

plt.plot(nlist, superror)
plt.show()
