import matplotlib.pyplot as plt
import numpy as np

n = [100, 1000, 10000, 100000]
index = 0

plt.figure(1, figsize=(14,7))
for i in n:
    index += 1
    uniform = open("uniform_sample_" + str(i) + ".txt")
    uniform_data = []
    for line in uniform:
        uniform_data.append(np.float64(line.rstrip("\n")))

    plt.subplot(2, len(n), index)
    plt.xlabel("[ a, b ]")
    plt.ylabel("frequency")
    plt.title("Uniform Distribution(n="+str(i)+")")
    plt.hist(uniform_data, 100, density=True, color='red')
plt.figure(1, figsize=(14,7))
for i in n:
    index += 1
    gauss = open("guass_sample_" + str(i) + ".txt")
    gauss_data = []
    for line in gauss:
        gauss_data.append(np.float64(line.rstrip("\n")))

    plt.subplot(2, len(n), index)
    plt.xlabel("[ m = 0.5, s = 1.5 ]")
    plt.ylabel("frequency")
    plt.title("Gauss Distribution(n="+str(i)+")")
    plt.hist(gauss_data, 100, density=True, color='blue')
plt.tight_layout()
plt.show()