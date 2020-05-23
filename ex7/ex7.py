from numpy import *
from matplotlib.pyplot import *
import time

def DFT(array,delta):
    out = []  # 輸出矩陣
    omega = [] # 頻率空間
    N = len(array)
    for j in range(0,N):
        this = 0
        omega.append((j*2*pi)/(delta*N))
        for k in range(0,N):
            this += array[k]*exp(-1j*2*pi*(j-N/2)*k/N) # sum
        out.append(this/(sqrt(N)))
    return out,omega

def Fourth(N,a,b,initial):
    y = []
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
    fff = 1.15*cos(t*2/3) - 0.5*arr[1] - sin(arr[0])
    return np.array([arr[1],fff])

def main():
    N = 2**12
    print(N)
    a = 0
    b = 2**5
    delta = (b-a)/(N-1)
    y,step = Fourth(N,a,b,[0,2])
    # out, ww = DFT(y,delta)
    x_axis = linspace(a,b,N)
    # print(len(y),len(x_axis))
    # subplot(211)
    title("y")
    plot(x_axis,y)
    # subplot(212)
    # title("DFT")
    # plot(ww,log10(abs(array(out))**2))
    show()


if __name__=="__main__":
    now = time.time()
    main()
    print(time.time()-now,"s")
