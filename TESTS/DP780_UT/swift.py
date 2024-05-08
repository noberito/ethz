import matplotlib.pyplot as plt
import numpy as np

a = 500 
b = 140
c = 50

def voce(epbar):
    return(a - b * np.exp(- c * epbar))

voce = np.vectorize(voce)

x = np.linspace(0, 0.2, 100)
y = voce(x)

print(x,y)

plt.plot(x, y)
plt.show()