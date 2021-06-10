"""
 Module to plot wav file data
 Assumes python3 and scipy  and matplotlib are installed
 Invoke it as  python3 WavPlot.py example_input1.wav 
"""
import sys

import matplotlib.pyplot as plt
import numpy as np

from scipy.io import wavfile
import scipy.io


for arg in sys.argv[1:]:
        print("Processing:", arg)
        wav_fname=arg
        samplerate, data = wavfile.read(wav_fname)
        try:
            length = data.shape[0] / samplerate
            print(f"num Audio channels = {data.shape[1]}", f"length = {length}s")
            time = np.linspace(0., length, data.shape[0])
            for i in range(data.shape[1]):
                plt.plot(time, data[:, i], label=wav_fname.split(".")[0]+"Channel"+str(i))
        except Exception as e:
            print (str(e))       

plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")
plt.show()
	


