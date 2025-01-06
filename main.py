import random
from alive_progress import alive_bar
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

boltzmann = 1 # 1.380 * 10 ** (-23)

Diagramme = False
laenge = 30
amino_auswahl = 20
Faltungs_schritte = 100000
Start_Temperatur = 1
Temperatur_Schritte = 10
# Im besten Fall ist Faltungs_schritte durch Temperatur_Schritte teilbar, sonst doff


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
        self.amino_positions_tuple = [tuple(position) for position in self.amino_positions]

        self.interaction_matrix = Wechselmatrix(self.amino_auswahl)
        self.surround_search()
        self.energie = Energiecalc(self.amino_class_list, self.interaction_matrix)

    def Protein_update(self):
        self.surround_search()
        self.energie = Energiecalc(self.amino_class_list, self.interaction_matrix)

    def __str__(self):
        return str(self.amino_positions)

    def Energie_amino_update(self, current_amino):
        # Ändert die Protein_eigenschaften nicht
        current_amino.energie = 0
        x = current_amino.amino_type
        if current_amino.neighbour:
            for neighbour in current_amino.neighbour:
                y = neighbour.amino_type
                current_amino.energie += self.interaction_matrix[x, y]

    def surround_search(self):
        # defines connections and neighbours for each amino-acid in the Protein
        for index, amin in enumerate(self.amino_class_list):
            self.single_search(amin)

    def single_search(self, amin):
        direc_list = np.array((np.asarray((0, 1)), np.asarray((0, -1)), np.asarray((-1, 0)), np.asarray((1, 0))))
        connected_step = []
        connected_classes = []
        if amin.index > 0:
            prev_connected = self.amino_class_list[amin.index - 1]
            connected_classes.append(prev_connected)
            connected_step.append(prev_connected.position - amin.position)
        if amin.index < laenge - 1:
            next_connected = self.amino_class_list[amin.index + 1]
            connected_classes.append(next_connected)
            connected_step.append(next_connected.position - amin.position)
        amin.connected = connected_classes
        aset = set([tuple(x) for x in connected_step])
        bset = set([tuple(x) for x in direc_list])
        free_direc = np.array([x for x in bset - aset])
        for direc in free_direc:
            neighbours = []
            if tuple(amin.position + direc) in self.amino_positions_tuple:
                neighbour_index = self.amino_positions_tuple.index(tuple(amin.position + direc))
                neighbours.append(self.amino_class_list[neighbour_index])
        amin.neighbour = neighbours

    def Position_swap(self, temp):
        random_index = np.random.randint(0, high=self.laenge)
        random_amino = self.amino_class_list[random_index]
        vector = [0, 0]

        for connection in random_amino.connected:
            vector += connection.position-random_amino.position
            if len(random_amino.connected) == 1:
                for index, cord in enumerate(vector):
                    if cord == 0:
                        chance = np.random.random()
                        if chance < 0.5:
                            cord = 1
                        else:
                            cord = -1
                        vector[index] = cord
        if tuple(random_amino.position+vector) != tuple(random_amino.position):
            if tuple(random_amino.position+vector) not in self.amino_positions_tuple:
                jump_pos = random_amino.position+vector
                clone = self.clone_Amino(random_amino, jump_pos)
                energie_diff = random_amino.energie - clone.energie
                swap = swap_check(temp, energie_diff)
                if swap:
                    self.amino_class_list[random_index] = clone
                    self.amino_positions[random_index] = jump_pos
                    self.amino_positions_tuple[random_index] = tuple(jump_pos)
                    self.Protein_update()

    def clone_Amino(self, amino_class, jump_location):
        # Erstellt nen Clone an der Position wo der Possible_Jump ist, interagiert noch nicht mit dem Protein
        clone = Aminosaeure(jump_location, self.amino_auswahl, amino_class.index)
        clone.connected = amino_class.connected
        clone.amino_type = amino_class.amino_type
        self.single_search(clone)
        self.Energie_amino_update(clone)
        return clone

class Aminosaeure:
    def __init__(self, pos, amino_auswahl, index):
        self.position = pos
        self.index = index
        self.amino_type = random.randint(0, amino_auswahl - 1)
        self.neighbour = []
        self.connected = []
        self.energie = 0

def avoiding_randomwalk(laenge, start):
    position = start
    position_array = np.zeros(shape=(laenge, 2))
    amino_list = []
    for i in range(laenge):
        amino_list.append(Aminosaeure(position, amino_auswahl, i))
        position_array[i] = position
        last_position = position
        dir_check = True
        trys = 0
        while dir_check:
            position = last_position
            direction = np.random.randint(1, 5)
            match direction:
                case 1:
                    position = np.asarray((position[0], position[1] + 1))
                case 2:
                    position = np.asarray((position[0], position[1] - 1))
                case 3:
                    position = np.asarray((position[0] + 1, position[1]))
                case 4:
                    position = np.asarray((position[0] - 1, position[1]))
            if tuple(position) not in list([tuple(position) for position in position_array]):
                dir_check = False
            trys += 1
            if trys > 4:
                return False
    return amino_list, position_array

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
        chance = np.random.random()
        barrier = np.exp(energie_diff/(temp*boltzmann))
        if chance < barrier:
            return True
        else:
            return False

def Abstand_A_O(Protein):
    A_pos, O_pos = Protein.amino_positions[0], Protein.amino_positions[-1]
    Abstand = np.sqrt((A_pos[0]-O_pos[0])**2 + (A_pos[1]-O_pos[1])**2)
    return Abstand


def Aufgabe_3():
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

        with alive_bar(iterations, force_tty=True) as bar:
            for i in range(iterations):
                msd_values = RandoMandoDangoLaengo(steps)
                all_msd_values.append(msd_values)
                bar()

        avg_msd = [sum(step_values) / len(step_values) for step_values in zip(*all_msd_values)]

        tau = np.arange(1, len(avg_msd) + 1)
        msd = np.array(avg_msd)

        # noinspection PyTupleAssignmentBalance
        popt, pcov = curve_fit(nen_fit, tau, msd, p0=(1, 0.5))

        msd_fit = nen_fit(tau, *popt)

        plt.figure(figsize=(12, 6))

        for index, msd_values in enumerate(all_msd_values):
            if index < 1000:
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

    selfavoiding_random_walk = True
    steps = 100
    iterations = 100000
    average_msd, all_msd, fit_params = multiple_iterations(steps, iterations)

def Aufgabe_4(laenge, amino_auswahl, temp, schritte):

    Protein_4 = Protein(laenge=laenge, amino_auswahl=amino_auswahl)

    plot1 = plt.subplot2grid((2, 2), (0, 0), colspan=2, rowspan=1)
    plot2 = plt.subplot2grid((2, 2), (1, 0), colspan=1, rowspan=1)
    plot3 = plt.subplot2grid((2, 2), (1, 1), colspan=1, rowspan=1)
    # Vorher Graph
    plot2.plot([point[0] for point in Protein_4.amino_positions], [point[1] for point in Protein_4.amino_positions])
    plot2.scatter([point[0] for point in Protein_4.amino_positions], [point[1] for point in Protein_4.amino_positions], s=100)
    plot2.set_xlabel(str("Abstand", Abstand_A_O(Protein_4)))
    plot2.grid()

    # Der wichtige Teil
    Energie_array = np.zeros(schritte)
    with alive_bar(schritte) as bar:
        for i in range(schritte):
            Protein_4.Position_swap(temp)
            Energie_array[i] = Protein_4.energie
            bar()

    # Nachher Graph
    plot3.plot([point[0] for point in Protein_4.amino_positions], [point[1] for point in Protein_4.amino_positions])
    plot3.scatter([point[0] for point in Protein_4.amino_positions], [point[1] for point in Protein_4.amino_positions], s=100)
    plot3.set_xlabel("Abstand", str(Abstand_A_O(Protein_4)))
    plot3.grid()

    # Energie Graph
    # plot1.plot(list(range(schritte)), Energie_array)
    plot1.scatter(list(range(schritte)), Energie_array)
    plot1.scatter([0, schritte], [Energie_array[0], Energie_array[-1]], c="red")
    plot1.plot([0, schritte], [Energie_array[0], Energie_array[-1]], c="red")
    plot1.grid()
    plt.show()

def Aufgabe_5(laenge, amino_auswahl, temp, schritte, Temp_schritte):
    # Durch np.linspace immer gleich Temperatur schritte, immer gleich viele Faltungen pro Temperatur
    temp_list = np.linspace(temp, Start_Temperatur/Temp_schritte, Temp_schritte)

    plot1 = plt.subplot2grid((2, 2), (0, 0), colspan=1, rowspan=1)
    plot2 = plt.subplot2grid((2, 2), (1, 0), colspan=1, rowspan=1)
    plot3 = plt.subplot2grid((2, 2), (0, 1), colspan=1, rowspan=1)
    plot4 = plt.subplot2grid((2, 2), (1, 1), colspan=1, rowspan=1)



    Protein_5 = Protein(laenge=laenge, amino_auswahl=amino_auswahl)
    Energie_array, Abstands_array = np.zeros(schritte), np.zeros(schritte)

    plot3.plot([point[0] for point in Protein_5.amino_positions], [point[1] for point in Protein_5.amino_positions])
    plot3.scatter([point[0] for point in Protein_5.amino_positions], [point[1] for point in Protein_5.amino_positions],
                  s=100)
    plot3.set_xlabel("Start")
    plot3.grid()

    with alive_bar(schritte) as bar:
        for index, temp in enumerate(temp_list):
            print(temp)
            for i in range(Faltungs_schritte//Temp_schritte):
                Protein_5.Position_swap(temp)
                Energie_array[i+index*Faltungs_schritte//Temp_schritte] = Protein_5.energie
                Abstands_array[i+index*Faltungs_schritte//Temp_schritte] = Abstand_A_O(Protein_5)
                bar()

    plot1.scatter(list(range(schritte)), Energie_array)
    plot1.scatter([0, schritte], [Energie_array[0], Energie_array[-1]], c="red")
    plot1.plot([0, schritte], [Energie_array[0], Energie_array[-1]], c="red")
    plot1.vlines(list(range(0,schritte, schritte//Temp_schritte)), np.min(Energie_array), np.max(Energie_array), colors="k", zorder=0)
    plot1.set_ylabel("Energie")
    plot1.grid()

    plot2.scatter(list(range(schritte)), Abstands_array)
    plot2.scatter([0, schritte], [Abstands_array[0], Abstands_array[-1]], c="red")
    plot2.plot([0, schritte], [Abstands_array[0], Abstands_array[-1]], c="red")
    plot2.vlines(list(range(0, schritte, schritte // Temp_schritte)), np.min(Abstands_array), np.max(Abstands_array), colors="k", zorder=0)
    plot2.set_ylabel("Abstand")
    plot2.grid()

    plot4.plot([point[0] for point in Protein_5.amino_positions], [point[1] for point in Protein_5.amino_positions])
    plot4.scatter([point[0] for point in Protein_5.amino_positions], [point[1] for point in Protein_5.amino_positions],
                  s=100)
    plot4.set_xlabel("Ende")
    plot4.grid()

    plt.show()

def main():
    print("YI STILL THE MAIN")
    # Aufgabe_3()
    # Aufgabe_4(laenge, amino_auswahl, temp=Start_Temperatur, schritte=Faltungs_schritte)
    Aufgabe_5(laenge, amino_auswahl, temp=Start_Temperatur, schritte=Faltungs_schritte, Temp_schritte=Temperatur_Schritte)
    print(np.linspace(Start_Temperatur, Start_Temperatur/Temperatur_Schritte, Temperatur_Schritte))

main()
