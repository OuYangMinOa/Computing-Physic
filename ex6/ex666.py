from matplotlib.pyplot import *
from numpy import *
func = lambda x:exp(-x**2)*exp(2j*x)
func_ft = lambda x:1/sqrt(2)*exp(-((2-x)**2)/4)

def DFT(array):
    out = []
    N = len(array)
    for j in range(0,N-1):
        this = 0
        for k in range(0,N-1):
            this += array[k]*exp(-1j*2*pi*(j-N/2)*k/N)
        out.append(this/(sqrt(N)))
    return out
def draw_fourier(N,A):
    Delta = 2*A/(N-1)
    fc = 1/(2*Delta)
    w_max = pi/Delta
    arr = linspace(-A,A,N)
    in_array = func(arr)
    out_array = DFT(in_array)
    new, ww = [], []
    for num,i in enumerate(out_array):
        wk =((num-N/2)*2*pi/N/Delta)
        ww.append(wk)
        new.append(i*exp(1j*wk*A))
    new = array(new) / (sqrt(2*pi)/sqrt(N)/(2*A/(N-1)))

    
    subplots()
    plot(ww,(new.real),'x')
    plot(ww,func_ft(array(ww)),'+')
    yscale("log")
    grid(True)
    show()


    
if __name__=="__main__":
    A = 7
    N = 2048
    draw_fourier(N,A)
