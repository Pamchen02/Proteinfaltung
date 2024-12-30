import random
from alive_progress import alive_bar
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

boltzmann = 1 # 1.380 * 10 ** (-23)

Diagramme = False
laenge = 30
amino_auswahl = 20
Faltungs_schritte = 1000


class Protein:
    def __init__(self, laenge, amino_auswahl):
        self.laenge = laenge
        self.position = np.zeros(2)
        self.amino_auswahl = amino_auswahl

        walk = False
        while not walk:
            walk = avoiding_randomwalk(self.laenge, self.position)
        self.amino_class_list = walk[0]
        self.amino_positions = walk[1]

        self.surround_search(first=True)
        self.interaction_matrix = Wechselmatrix(self.amino_auswahl)
        self.energie = Energiecalc(self.amino_class_list, self.interaction_matrix)

    def Energie_update(self):
        self.energie = Energiecalc(self.amino_class_list, self.interaction_matrix)

    def __str__(self):
        return str(self.amino_positions)

    def clone_Amino(self, amino_class, index):
        """Erstellt nen Clone an der Position wo der Possible_Jump ist, interagiert noch nicht mit dem Protein"""
        clone = Aminosaeure(amino_class.possible_jumps[0], self.amino_auswahl, index)
        clone.connected = amino_class.connected
        clone.amino_type = amino_class.amino_type
        self.Neighbour_amino_update(clone)
        self.Energie_amino_update(clone)
        return clone

    def Position_swap(self, temp):
        potentialle_swaps = [(amino_class, index) for index, amino_class in enumerate(self.amino_class_list) if amino_class.possible_jumps]
        random_amino, random_amino_index = potentialle_swaps[random.randint(0, len(potentialle_swaps)-1)]
        random_clone = self.clone_Amino(random_amino, random_amino.index)
        energie_diff = random_clone.energie - random_amino.energie
        swap = swap_check(temp, energie_diff)  # Gibt True oder False zurück
        global debug_index
        debug_index = 1
        if swap:
            debug_index = random_amino_index
            self.energie += energie_diff
            self.amino_class_list[random_amino_index] = random_clone
            self.amino_positions[random_amino_index] = random_clone.position
            self.surround_search(first=False)
            self.energie = Energiecalc(self.amino_class_list, self.interaction_matrix)
            """position_check_list = []    # hier sind Indices drin
            old_pos, new_pos = random_amino.position, random_clone.position
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    old_check = (old_pos[0] + dx, old_pos[1] + dy)
                    new_check = (new_pos[0] + dx, new_pos[1] + dy)
                    if old_check in self.amino_positions:
                        if self.amino_positions.index(old_check) not in position_check_list:
                            position_check_list.append(self.amino_positions.index(old_check))
                    if new_check in self.amino_positions:
                        if self.amino_positions.index(new_check) not in position_check_list:
                            position_check_list.append(self.amino_positions.index(new_check))
            for index in position_check_list:
                amino_to_check = self.amino_class_list[index]
                self.full_correction_amino(amino_to_check)"""

    def Neighbour_amino_update(self, current_amino):
        # Ändert die Protein_eigenschaften nicht
        direc_list = np.array(((0, 1), (0, -1), (1, 0), (-1, 0)))
        position = np.array(current_amino.position)
        connected_step = [x.position - position for x in current_amino.connected]
        direc_set = set([tuple(x) for x in direc_list])
        neigh_set = set([tuple(x) for x in connected_step])
        free_direc = np.array([x for x in direc_set - neigh_set])
        current_amino.neighbour, current_amino.possible_jumps = [], []
        pos_list = [tuple(x) for x in self.amino_positions]
        for direc in free_direc:
            destination_tuple = tuple(current_amino.position + direc)
            if destination_tuple in pos_list:
                current_amino.neighbour.append(self.amino_class_list[pos_list.index(destination_tuple)])
        if np.linalg.norm(np.sum(connected_step, axis=0)) > 1:
            jump_location = tuple(current_amino.position + connected_step[0] + connected_step[1])
            if jump_location not in pos_list:
                current_amino.possible_jumps.append(jump_location)
        elif np.linalg.norm(np.sum(connected_step, axis=0)) == 1:
            dings = tuple(np.sum(connected_step, axis=0))
            jump_1, jump_2 = (0, 0), (0, 0)
            if dings[0] == 0:
                jump_1 = (position[0] + 1, position[1] + dings[1])
                jump_2 = (position[0] - 1, position[1] + dings[1])
            if dings[1] == 0:
                jump_1 = (position[0] + dings[0], position[1] + 1)
                jump_2 = (position[0] + dings[0], position[1] - 1)
            before = random.random()
            if before > 0.5:
                if jump_1 not in pos_list:
                    current_amino.possible_jumps.append(jump_1)
                if jump_2 not in pos_list:
                    current_amino.possible_jumps.append(jump_2)
            else:
                if jump_2 not in pos_list:
                    current_amino.possible_jumps.append(jump_2)
                if jump_1 not in pos_list:
                    current_amino.possible_jumps.append(jump_1)

    def Energie_amino_update(self, current_amino):
        # Ändert die Protein_eigenschaften nicht
        current_amino.energie = 0
        x = current_amino.amino_type
        if current_amino.neighbour:
            for neighbour in current_amino.neighbour:
                y = neighbour.amino_type
                current_amino.energie += self.interaction_matrix[x, y]

    def full_correction_amino(self, amino_class):
        # Ändert die Protein_eigenschafte Energie auf den korrigierten Wert
        self.Neighbour_amino_update(amino_class)

        self.energie -= amino_class.energie
        self.Energie_amino_update(amino_class)
        self.energie += amino_class.energie

    def surround_search(self, first):
        """ Schaut welche Aminosäure zu welcher, per kovalenter Bindung, connected ist,
        sucht die nicht connectedten Nachbarn
        und welche Spots frei sind für Positions swaps bei der die kovalenten Bindungen erhalten bleiben."""
        pos_list = self.amino_positions
        direc_list = np.array((np.asarray((0, 1)), np.asarray((0, -1)), np.asarray((-1, 0)), np.asarray((1, 0))))
        if first:
            laenge = len(self.amino_positions)
            for index, amin in enumerate(self.amino_class_list):
                if index > 0:
                    prev_connected = self.amino_class_list[index - 1]
                    amin.connected.append(prev_connected)
                if index < laenge - 1:
                    next_connected = self.amino_class_list[index + 1]
                    amin.connected.append(next_connected)


        if not first:
            for index, amin in enumerate(self.amino_class_list):
                amin.neighbour = []
                amin.possible_jumps = []

        for index, amin in enumerate(self.amino_class_list):
            position = np.array(amin.position)
            connected_step = np.asarray([amino.position - position for amino in amin.connected])
            aset = set([tuple(x) for x in connected_step])
            bset = set([tuple(x) for x in direc_list])
            free_direc = np.array([x for x in bset - aset])
            for direc in free_direc:
                destination_tuple = tuple(amin.position + direc)
                if destination_tuple in pos_list:
                    amin.neighbour.append(self.amino_class_list[pos_list.index(destination_tuple)])
            if np.linalg.norm(np.sum(connected_step, axis=0)) > 1:
                jump_location = tuple(amin.position + connected_step[0] + connected_step[1])
                if jump_location not in pos_list:
                    amin.possible_jumps.append(jump_location)
            elif np.linalg.norm(np.sum(connected_step, axis=0)) == 1:
                dings = tuple(np.sum(connected_step, axis=0))
                jump_1, jump_2 = (0, 0), (0, 0)
                if dings[0] == 0:
                    jump_1 = (position[0] + 1, position[1] + dings[1])
                    jump_2 = (position[0] - 1, position[1] + dings[1])
                if dings[1] == 0:
                    jump_1 = (position[0] + dings[0], position[1] + 1)
                    jump_2 = (position[0] + dings[0], position[1] - 1)
                before = random.random()
                if before > 0.5:
                    if jump_1 not in pos_list:
                        amin.possible_jumps.append(jump_1)
                    if jump_2 not in pos_list:
                        amin.possible_jumps.append(jump_2)
                else:
                    if jump_2 not in pos_list:
                        amin.possible_jumps.append(jump_2)
                    if jump_1 not in pos_list:
                        amin.possible_jumps.append(jump_1)


class Aminosaeure:
    def __init__(self, pos, amino_auswahl, index):
        self.position = pos
        self.index = index
        self.amino_type = random.randint(0, amino_auswahl - 1)
        self.neighbour = []
        self.connected = []
        self.possible_jumps = []
        self.Abstand_next = 1
        self.Abstand_prev = -1
        self.energie = 0


def avoiding_randomwalk(laenge, start):
    position = start
    position_list = []
    amino_list = []
    for i in range(laenge):
        amino_list.append(Aminosaeure(position, amino_auswahl, i))
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
            if not np.isin(position, position_list).all():
                dir_check = False
            trys += 1
            if trys > 4:
                return False
    return amino_list, position_list


def RandoMandoDangoLaengo(steps):
    def RandoMando(steps, start=(0, 0)):
        position = start
        position_set = {position}
        position_list = [position]
        msd_cumulative_list = [0]
        cumulative_sum = 0

        for step in range(1, steps + 1):
            dir_check = True
            tries = 0
            while dir_check:
                new_position = position
                direction = random.randint(1, 4)
                match direction:
                    case 1:
                        new_position = (position[0], position[1] + 1)
                    case 2:
                        new_position = (position[0], position[1] - 1)
                    case 3:
                        new_position = (position[0] + 1, position[1])
                    case 4:
                        new_position = (position[0] - 1, position[1])

                if selfavoiding_random_walk:
                    if new_position not in position_set:
                        position = new_position
                        position_set.add(position)
                        position_list.append(position)
                        dir_check = False
                else:
                    position = new_position
                    position_list.append(position)
                    dir_check = False

                tries += 1
                if tries > 4:
                    return None, None

            squared_distance = position[0] ** 2 + position[1] ** 2
            cumulative_sum += squared_distance
            msd_cumulative = cumulative_sum / step
            msd_cumulative_list.append(msd_cumulative)

        return position_list, msd_cumulative_list

    while True:
        positions, msd_cumulative_values = RandoMando(steps)
        if positions and msd_cumulative_values:
            return msd_cumulative_values


def nen_fit(tau, a, alpha):
    return a * tau ** alpha


def multiple_iterations(steps, iterations):
    all_msd_values = []

    for i in range(iterations):
        msd_values = RandoMandoDangoLaengo(steps)
        all_msd_values.append(msd_values)

    avg_msd = [sum(step_values) / len(step_values) for step_values in zip(*all_msd_values)]

    tau = np.arange(1, len(avg_msd) + 1)
    msd = np.array(avg_msd)

    popt, pcov = curve_fit(nen_fit, tau, msd, p0=(1, 0.5))

    msd_fit = nen_fit(tau, *popt)

    plt.figure(figsize=(12, 6))

    for msd_values in all_msd_values:
        plt.plot(range(1, len(msd_values) + 1), msd_values, alpha=0.3, linewidth=1, color='gray')

    plt.plot(tau, msd_fit, color='blue', linestyle='--', linewidth=2,
             label=f'Fit: a={popt[0]:.4f}, α={popt[1]:.4f}')

    plt.xlabel("Schritte")
    plt.ylabel("Mittlere quadratische Entfernung (MSD)")
    plt.title(f"MSD für {iterations} Iterationen und {steps} Schritte")
    plt.legend()
    plt.grid()
    plt.show()

    return avg_msd, all_msd_values, popt


selfavoiding_random_walk = False

steps = 250
iterations = 20

average_msd, all_msd, fit_params = multiple_iterations(steps, iterations)

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
        amino.energie = 0
        x = amino.amino_type
        if amino.neighbour:
            for neighbour in amino.neighbour:
                y = neighbour.amino_type
                amino.energie += Matrix[x, y]
                energie += Matrix[x, y]
    return energie


def swap_check(temp, energie_diff):
    if energie_diff > 0:
        return True
    else:
        chance = random.random()
        barrier = np.exp(energie_diff/(temp*boltzmann))
        if chance < barrier:
            return True
        else:
            return False

def check_connections(amin):
    for connection in amin.connected:
        Abstand_tuple = (amin.position[0]-connection.position[0],amin.position[1]-connection.position[1])
        direc_list = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        if Abstand_tuple not in direc_list:
            print("Problem")
    return True


Protein1 = Protein(laenge=30, amino_auswahl=20)
with alive_bar(Faltungs_schritte) as bar:
    for i in range(Faltungs_schritte):
        """plt.plot([point[0] for point in Protein1.amino_positions], [point[1] for point in Protein1.amino_positions])
        plt.scatter([point[0] for point in Protein1.amino_positions], [point[1] for point in Protein1.amino_positions], s=100)"""
        Protein1.Position_swap(1)
        """plt.plot([point[0] for point in Protein1.amino_positions], [point[1] for point in Protein1.amino_positions])
        plt.scatter([point[0] for point in Protein1.amino_positions], [point[1] for point in Protein1.amino_positions], s=100)
        plt.grid()
        plt.scatter(Protein1.amino_positions[debug_index][0], Protein1.amino_positions[debug_index][1], s=100)
        plt.show()"""
        for amino in Protein1.amino_class_list:
            check_connections(amino)
            for connect in amino.connected:
                if amino.index - connect.index != 1 and amino.index - connect.index != -1:
                    print("ALARM!", amino.index, connect.index)
        bar()

    plt.plot([point[0] for point in Protein1.amino_positions], [point[1] for point in Protein1.amino_positions])
    plt.scatter([point[0] for point in Protein1.amino_positions], [point[1] for point in Protein1.amino_positions], s=100)
    plt.grid()
    plt.show()
