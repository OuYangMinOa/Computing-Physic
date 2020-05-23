# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy.linalg import eig
#from scipy import linalg as SA
import time
# the posibility to transition
tx = 1
ty = 1
# V(x)
U = lambda x, y:  0#0.5 if (x>3 and x<6 and y>3 and y<6) else 0
size = 30
x_cod = []
y_cod = [] # ψ''(x) = t0(ψi+1,i  -2ψi,i + ψi-1,i)
z_cod = [] # ψ''(y) = t0(ψi,i+1  -2ψi,i + ψi,i-1)
# condition (importance)
def persi(fi,fj,i,j):
    if (fi==i and fj == j):
        return 2*tx + 2*ty + U(i,j)
    elif (fi ==i and abs(fj-j)==1):
        return -1*(tx ) + U(i,j)
    elif (fj ==j and abs(fi-i)==1):
        return -1*(ty ) + U(i,j)
    else:
        return 0
# show Hamitonian (you can ignore this func)
def post(dont_ignore_column = False,dont_ignore_row=False):
    time_2 = time_1 = 0
    # if Hamitonian to big (when size = 6,7,8,9...  use '...' to ignore)
    if (size >5):
        for m,i in enumerate(re_z):
            if ( not dont_ignore_column):
                if (m>5 and m < len(i)-5):
                        if (time_2<4):
                            print('.')
                            time_2 += 1
                        continue
            time_1 = 0
            for num,j in enumerate(i):
                if (not dont_ignore_row):
                    if (num>5 and num < len(i)-5):
                        if (time_1<4):
                            print('.',end='')
                            time_1 += 1
                        continue
                print('{:4}'.format(j),end=' ')
            print(' {}'.format('ψ'),
                  int(np.array(x_cod).reshape(size**2,1)[m]),
                  ',',
                  int(np.array(y_cod).reshape(size**2,1)[m]),sep='')
        return
     # if Hamitonian not to big
    for m,i in enumerate(re_z):
        for j in i:
            print('{:4}'.format(j),end=' ')
        print(' {}'.format('ψ'),int(np.array(x_cod).reshape(size**2,1)[m]),',',int(np.array(y_cod).reshape(size**2,1)[m]),sep='')
# make x_cod y_cod (you can use mershgrid)
for i in range(size):
    for j in range(size):
        x_cod.append(j)
        y_cod.append(i)
# make Hamitonian
for i in range(size**2):
    this = []
    for j in range(size**2):
        z_cod.append(persi(x_cod[i],y_cod[i],x_cod[j],y_cod[j]))
# reshaping
(x_cod,y_cod,re_z)  =  (np.array(x_cod).reshape(size,size),
                        np.array(y_cod).reshape(size,size),
                        np.array(z_cod).reshape(size**2,size**2))
# get value and vector
value, vector       =  eig(re_z)
# sorted index
idx                 =  np.argsort(value)[::-1]
# sort value
value               =  np.array(value[idx]).reshape(size,size)
# rot vector and rotation 90 degree
vector              =  np.rot90(np.array(vector[:,idx]))[::-1]
#print vector
print(vector)
# show Hamitonian
post()
# show x_cod y_cod
# create fig
fig = plt.figure()
# make fig to 3D
ax = Axes3D(fig)
# plot value or vector (False or True)
do_or_not = 1
if (do_or_not):
    # plot moving eigenvector animation fig
    while 1:
        for i,j in enumerate(vector):
            fig.suptitle('N : '+str(i+1))
            ax.plot_surface(x_cod,y_cod,abs((np.array((j)).reshape(size,size))),cmap=plt.cm.winter)
            # wait for 0.1s
            plt.pause(0.1)
            # clear image
            plt.cla()
            if (i>30): break
    # plot stoppable image 
    ax.plot_surface(x_cod,y_cod,abs((vector[3]).reshape(size,size)),cmap=plt.cm.winter)
else:
    # plot eigenvalue
    surf = ax.plot_surface(y_cod,x_cod,abs(value),cmap=plt.cm.winter)
    fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

