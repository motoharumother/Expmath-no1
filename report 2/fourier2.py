import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

pi = np.pi
x = np.arange(-pi, pi, pi/1000)

#関数定義
def originalfunc(x):
    return np.cos(4*x) + x
#係数
def ak(k):
    tmpa = integrate.quad(lambda x, c : (1/pi)*originalfunc(x)*np.cos(x*c), -pi, pi, args=k)[0]
    return tmpa

def bk(k):
    tmpb = integrate.quad(lambda x, c : (1/pi)*originalfunc(x)*np.sin(x*c), -pi, pi, args=k)[0]
    return tmpb


#Sn(x)の計算
def Sn(n,ak,bk,x):
    cos_list = np.array([np.cos((k+1)*x) for k in range(n)])
    sin_list = np.array([np.sin((k+1)*x) for k in range(n)])
    sn = cos_list*ak + sin_list*bk
    return sn.sum()


nlist = np.arange(1,16,5)

maxlist = np.zeros((len(nlist)))

y = np.zeros((len(nlist), len(x)))

a_zero = integrate.quad(lambda x: (1/pi)*originalfunc(x), -pi, pi)[0]/2

plt.figure(figsize=(10,7))
plt.plot(x, originalfunc(x), label='cos(x) + x', c='r')

for m, n in enumerate(nlist):
    a = np.zeros((n))
    b = np.zeros((n))
    for k in range(n):
        a[k] = ak(k+1)
        b[k] = bk(k+1)

    for i, _x in enumerate(x):
        y[m][i] = Sn(n=n,ak=a,bk=b,x=_x) + a_zero
    plt.plot(x, y[m], label='n={}'.format(n), alpha=0.4)
#２乗差
er = [((x-y[i])**2).sum() for i in range(len(nlist))]
#sup error
superror = [(abs(x-y[i])).max() for i in range(len(nlist))]
    

plt.legend()
plt.show()

"""
plt.plot(nlist, er)
plt.scatter(nlist, er, color='r')
plt.show()

plt.plot(nlist, superror)
plt.show()
"""
