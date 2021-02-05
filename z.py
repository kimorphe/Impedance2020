import numpy as np
import matplotlib.pyplot as plt

x=np.logspace(-2,3,200)
L=8.0
p0=0.5
w=(1j*x)**p0
y1=np.tanh(w*L)/w

alpha=0.4
p=p0+alpha
w=(1j*x)**p
#y2=1/(1+w)
Tf=2
y2=1/(w*Tf)      # CPE impedance 

y=y1*y2/(y1+y2)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(np.real(y), -np.imag(y))
ax.grid(True)
ax.set_aspect(1.0)
plt.show()
