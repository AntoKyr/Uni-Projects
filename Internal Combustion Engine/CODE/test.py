import numpy as np
from matplotlib import pyplot as plt

print(np.polyval(np.polyfit(x=[1500, 2000, 2500, 3000, 4000, 4500, 5000], y=3*np.array([1.7, 1.9, 2.2, 2.55, 3.45, 3.9, 4.4])/4.65, deg=2), 2600))