import random
from alive_progress import alive_bar
import numpy as np

class Protein:
    def __init__(self, laenge):
        self.laenge = laenge
        self.start = (0, 0)
        walk = False
        while not walk:
            walk = randomwalk(self.laenge, self.start)
        self.amino_class = walk[0]
        self.positions = walk[1]
        surround_search(self.amino_class, self.positions)

    def __str__(self):
        return str(self.positions)


class Aminosaeure:
    def __init__(self, pos):
        self.position = pos
        self.amino = random.randint(1, 20)
        self.neighbour = []
        self.connected = []
        self.possible_jumps = []


def randomwalk(laenge, start):
    position = start
    position_list = []
    amino_list = []
    with alive_bar(laenge, force_tty=True) as bar:
        for i in range(laenge):
            amino_list.append(Aminosaeure(position))
            position_list.append(position)
            last_position = position
            dir_check = True
            trys = 0
            while dir_check:
                position = last_position
                direction = random.randint(1, 4)
                match direction:
                    case 1:
                        position = (position[0], position[1] + 1)
                    case 2:
                        position = (position[0], position[1] - 1)
                    case 3:
                        position = (position[0] + 1, position[1])
                    case 4:
                        position = (position[0] - 1, position[1])
                if position not in position_list:
                    dir_check = False
                trys += 1
                if trys > 4:
                    return False
            bar()
    return (amino_list, position_list)


"""
Schaut welche Aminosäure zu welcher, per kovalenter Bindung, connected ist,
sucht die nicht connectedten Nachbarn 
und welche Spots frei sind für Positions swaps bei der die kovalenten Bindungen erhalten bleiben.
"""
def surround_search(amino_list, pos_array):
    laenge = len(pos_array)
    for index, amin in enumerate(amino_list):
        direc_list = np.array(((0, 1), (0, -1), (1, 0), (-1, 0)))
        position = np.array(amin.position)
        if index > 0:
            prev_connected = amino_list[index-1]
            amin.connected.append(prev_connected)
        if index < laenge - 1:
            next_connected = amino_list[index+1]
            amin.connected.append(next_connected)
        connected_step = [x.position - position for x in amin.connected]
        direc_set = set([tuple(x) for x in direc_list])
        neigh_set = set([tuple(x) for x in connected_step])
        free_direc = np.array([x for x in direc_set - neigh_set])

        for direc in free_direc:
            destination_tuple = tuple(amin.position + direc)
            pos_list = [tuple(x) for x in pos_array]
            if destination_tuple in pos_list:
                amin.neighbour.append(amino_list[pos_list.index(destination_tuple)])
            else:
                if np.linalg.norm(np.sum(free_direc, axis=0)) > 1:
                    amin.possible_jumps.append(destination_tuple)


np.array(Protein(laenge=30).positions)
