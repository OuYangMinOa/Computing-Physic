from matplotlib.pyplot import *
from numpy import  *

class exersice6:
    def __init__(self,N,A):
        self.N = N
        self.A = A
        self.omega = 2*A/(N-1)
        self.fc = 1/(2*self.omega)
        self.im = complex(0,1)
        self.x = linspace(-A,A,N)
        self.oringinal_data = self.oringinal_func(self.x)
        plot(self.x,real(self.DFT(self.oringinal_data)))
        show()
    def DFT(self,array):
        out = []
        for j in range(self.N):
            this = 0
            for k in range(self.N):
                ww = 2*pi*(j-self.N/2)/self.N
                this += (array[k]*exp(-self.im*ww*k))
            out.append(this)
        return out
    def oringinal_func(self,x):
        return exp(-x**2)*exp(2*self.im*x)
    def exact_function(self,x):
        return 1/sqrt(2)*exp(-(2-x)**2/4)
exersice6(2048,7
          )
