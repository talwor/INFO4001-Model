import random

print(str(random.random()))


for i in range(100):

    if random.random() < 0.01:
        print("skib" + str(i))