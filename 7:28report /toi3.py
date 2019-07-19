import numpy as np
import matplotlib.pyplot as plt


T = 1
n = 50
dt = T/n
t = np.arange(0,10,dt)

y = np.zeros((2,len(t)))

y[0][0] = 2
y[1][0] = 1

def f(x):
    X = x.copy()
    X[0] = x[0]*(x[1]-2)
    X[1] = (1-x[0])*x[1]
    return X


def ExEuler(y, dt, t):
    for i in range(len(t)-1):
        y[:,i+1] = y[:,i] + dt*f(y[:,i])
    return y

def Runge_kutta(y, dt, t):
    for i in range(len(t)-1):
        y1 = y[:,i]
        y2 = y1 + (dt/2)*f(y1)
        y3 = y1 + (dt/2)*f(y2)
        y4 = y1 + dt*f(y3)
        y[:,i+1] = y1 + (dt/6)*(f(y1)+2*f(y2)+2*f(y3)+f(y4))
    return y

def SympleticEuler(y, dt, t):
    for i in range(len(t)-1):
        A = np.linalg.inv(np.array([[1,-dt*y[0,i]],[dt*y[1,i], 1]]))
        y[:,i+1] = np.dot(A,np.array([(1-2*dt)*y[0,i],(1+dt)*y[1,i]]))
    return y



#y = ExEuler(y, dt, t)
y = Runge_kutta(y, dt, t)
#y = SympleticEuler(y, dt, t)
plt.figure(figsize=(8,8))
plt.plot(y[0], y[1])
plt.xlabel('u')
plt.ylabel('v')
plt.show()
#
