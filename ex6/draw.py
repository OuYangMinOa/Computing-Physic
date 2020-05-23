from matplotlib.pyplot import *
from numpy import  *
out = []
with open(r"D:\計算物理\ex6\output.txt","r") as f:
    while True:
        line = f.readline()
        if (line==""):break
        out.append(line.strip().split())
        
new = rot90(array(out)).astype(np.float)
plot(*new[:-1],'x')
show()
