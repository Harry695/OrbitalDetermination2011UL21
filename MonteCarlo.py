from math import*
from vpython import*
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt

N = 100000

canvas(center=vector(50,50,0), xmin = 0, xmax=100, ymin=0, ymax=100)
particle = sphere(pos=vector(50, 50, 0), radius=1, make_trail=True)

def moveX(obj, dist):
    print("moved x by", dist)
    obj.pos.x += dist

def moveY(obj, dist):
    print("moved y by", dist)
    obj.pos.y += dist

for i in range(N):
    rate(10)
    # moveY(particle, 1)
    movement = random.choice([-1, 1])
    # print(movement)
    random.choice([moveX, moveY])(particle, movement)

# hypersphere
# volumeList = []

# for j in range(1000):
#     count = 0
#     for i in range(N):
#         coords = random.uniform(-1, 1, 10)
#         sum = np.sum(coords ** 2)
#         if (sqrt(sum) < 1):
#             count += 1
#             # print(sqrt(sum))
#     # NCircle / N = VCircle / VSquare
#     volume = count / N * 2**10
#     print("Volume", volume)
#     volumeList.append(volume)
# # print(volumeList)
# plt.plot(volumeList)
# print("Mean", np.mean(volumeList))
# print("StDev", np.std(volumeList))
# plt.show()