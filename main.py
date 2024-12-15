import random
# from alive_progress import alive_bar
import numpy as np
import matplotlib.pyplot as plt

boltzmann = 1.380*10**(-23)


Diagramme = False
laenge = 30
amino_auswahl = 20


class Protein:
    def __init__(self, laenge, amino_auswahl):
        self.laenge = laenge
        self.position = (0, 0)
        self.amino_auswahl = amino_auswahl

        walk = False
        while not walk:
            walk = avoiding_randomwalk(self.laenge, self.position)
        self.amino_class_list = walk[0]
        self.amino_positions = walk[1]

        first_surround_search(self.amino_class_list, self.amino_positions)
        self.interaction_matrix = Wechselmatrix(self.amino_auswahl)
        self.energie = Energiecalc(self.amino_class_list, self.interaction_matrix)

    def __str__(self):
        return str(self.amino_positions)

    def Position_swap(self, temp):
        potentialle_swaps = [(amino_class, index) for index, amino_class in enumerate(self.amino_class_list) if amino_class.possible_jumps]
        random_amino, random_amino_index = potentialle_swaps[random.randint(0, len(potentialle_swaps)-1)]
        old_pos = random_amino.position
        random_amino.position = random_amino.possible_jumps[0]
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                check_pos = (old_pos[0]+dx, old_pos[1]+dy)
                if check_pos in self.amino_positions:
                    print("moin")

        # warum ist das nicht vernümpftig gemacht

    def Neighour_check(self, current_amino):
        direc_list = np.array(((0, 1), (0, -1), (1, 0), (-1, 0)))
        position = np.array(current_amino.position)
        connected_step = [x.position - position for x in current_amino.connected]
        direc_set = set([tuple(x) for x in direc_list])
        neigh_set = set([tuple(x) for x in connected_step])
        free_direc = np.array([x for x in direc_set - neigh_set])
        current_amino.neighbour, current_amino.possible_jumps = [], []
        for direc in free_direc:
            destination_tuple = tuple(current_amino.position + direc)
            pos_list = [tuple(x) for x in self.amino_positions]
            if destination_tuple in pos_list:
                current_amino.neighbour.append(self.amino_class_list[pos_list.index(destination_tuple)])
            else:
                if np.linalg.norm(np.sum(free_direc, axis=0)) > 1:
                    current_amino.possible_jumps.append(destination_tuple)

    def Energie_check(self, current_amino):
        self.energie -= current_amino.energie
        current_amino.energie = 0
        x = current_amino.amino_type
        if current_amino.neighbour:
            for neighbour in current_amino.neighbour:
                y = neighbour.amino_type
                current_amino.energie += self.interaction_matrix[x, y]
                self.energie += self.interaction_matrix[x, y]

    def full_check_amino(self, index):
        current_amino = self.amino_class_list[index]
        Protein.Neighour_check(self, current_amino)
        Protein.Energie_check(self, current_amino)

class Aminosaeure:
    def __init__(self, pos, amino_auswahl):
        self.position = pos
        self.amino_type = random.randint(0, amino_auswahl - 1)
        self.neighbour = []
        self.connected = []
        self.possible_jumps = []
        self.energie = 0


def avoiding_randomwalk(laenge, start):
    position = start
    position_list = []
    amino_list = []
    for i in range(laenge):
        amino_list.append(Aminosaeure(position, amino_auswahl))
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
    return amino_list, position_list


"""
Schaut welche Aminosäure zu welcher, per kovalenter Bindung, connected ist,
sucht die nicht connectedten Nachbarn 
und welche Spots frei sind für Positions swaps bei der die kovalenten Bindungen erhalten bleiben.
"""


def first_surround_search(amino_list, pos_array):
    laenge = len(pos_array)
    for index, amin in enumerate(amino_list):
        direc_list = np.array(((0, 1), (0, -1), (1, 0), (-1, 0)))
        position = np.array(amin.position)
        if index > 0:
            prev_connected = amino_list[index - 1]
            amin.connected.append(prev_connected)
        if index < laenge - 1:
            next_connected = amino_list[index + 1]
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


def Wechselmatrix(laenge):
    matrix = np.zeros((laenge, laenge))
    sigma = 1 / np.sqrt(2)
    mu = -3
    for x in range(laenge):
        for y in range(x + 1):
            matrix[x, y] = np.random.normal(loc=mu, scale=sigma)
            matrix[y, x] = matrix[x, y]

    if Diagramme:
        plt.imshow(matrix)
        plt.show()

    eigenwerte = np.linalg.eigvalsh(matrix)
    sort_eigen = np.abs(np.append(eigenwerte[eigenwerte > 0], eigenwerte[eigenwerte < 0]))

    if Diagramme:
        plt.bar(list(range(laenge)), sort_eigen)
        plt.show()

    return matrix


def Energiecalc(amino_class_list, Matrix):
    energie = 0
    for amino in amino_class_list:
        x = amino.amino_type
        if amino.neighbour:
            for neighbour in amino.neighbour:
                y = neighbour.amino_type
                amino.energie += Matrix[x, y]
                energie += Matrix[x, y]
    print(energie)
    return energie



np.array(Protein(laenge=30, amino_auswahl=20).amino_positions)
