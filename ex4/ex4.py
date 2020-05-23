import matplotlib.pyplot as plt
import numpy as np


def Fourth(N,a,b,initial):
    y = [0]
    step  = abs(b-a)/N
    for i in range(N):
        t   = a+step*i
        c1 = step* generation(t  ,initial)
        c2 = step* generation(t+step/2 ,initial + c1/2)
        c3 = step* generation(t+step/2,initial + c2/2)
        c4 = step* generation(t+step,initial + c3)
        
        sum_c = np.add(np.add(c1,2*c2),np.add(c3*2,c4))/6
        
        initial = np.add(initial,sum_c)
        y.append(initial[0])
    return y,step

def generation(t,arr):
    fff = (np.exp(2*t)*np.sin(t) + 2*arr[1]-2*arr[0])
    return np.array([arr[1],fff])

def main():
    N = 100
    a = 0
    b = 2
    initial = [0,0]
    y,h = Fourth(N,a,b,initial)
    plt.plot(np.linspace(a,b,N+1),y)
    plt.show()
    print(y[-1])
main()
