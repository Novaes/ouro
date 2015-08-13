import time
import numpy as np
from scipy import signal as sg


x = np.mgrid[:2048, :2048][0]
y = np.random.rand(2048,2048)
start_time = time.time()
np.fft.ifft2(np.fft.fft2(x)*np.fft.fft2(y))
print(" %s " % (time.time() - start_time))
# print sg.convolve2d(x,y)