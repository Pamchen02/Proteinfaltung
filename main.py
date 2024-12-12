import random


class Gitter:
    def __init__(self, laenge):
        self.laenge = laenge


class Gitter_punkt:
    def __init__(self, pos):
        self.position = pos
        self.bound =[]
        self.nachbarn = []



x = [random.randint(0, 9) for i in range(5)]
y = [random.randint(0, 9) for a in range(5)]

print(x, y)
for index, i in enumerate(x):
    c_1 = Gitter_punkt((i, y[index]))
    print(c_1.position)
