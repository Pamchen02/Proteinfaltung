"""
 def clone_Amino(self, amino_class, index):
        # Erstellt nen Clone an der Position wo der Possible_Jump ist, interagiert noch nicht mit dem Protein
        clone = Aminosaeure(amino_class.possible_jumps[0], self.amino_auswahl, index)
        clone.connected = amino_class.connected
        clone.amino_type = amino_class.amino_type
        self.Neighbour_amino_update(clone)
        self.Energie_amino_update(clone)
        return clone
"""
"""
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

            Das gehört nur so halb dazu
            position_check_list = []    # hier sind Indices drin
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
                self.full_correction_amino(amino_to_check)

"""
"""
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
        sum_conn = np.sum(connected_step, axis=0)
        if np.sqrt(sum_conn[0]**2 + sum_conn[1]**2) > 1:
            jump_location = tuple(current_amino.position + connected_step[0] + connected_step[1])
            if jump_location not in pos_list:
                current_amino.possible_jumps.append(jump_location)
        elif np.sqrt(sum_conn[0]**2 + sum_conn[1]**2) == 1:
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
"""
"""
    def full_correction_amino(self, amino_class):
        # Ändert die Protein_eigenschafte Energie auf den korrigierten Wert
        self.Neighbour_amino_update(amino_class)
        self.energie -= amino_class.energie
        self.Energie_amino_update(amino_class)
        self.energie += amino_class.energie
"""
"""    
def surround_search(self, first):
        Schaut welche Aminosäure zu welcher, per kovalenter Bindung, connected ist,
        sucht die nicht connectedten Nachbarn
        und welche Spots frei sind für Positions swaps bei der die kovalenten Bindungen erhalten bleiben.
        pos_array = self.amino_positions
        print(pos_array)
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
            print(index)
            position = np.array(amin.position)
            connected_step = np.asarray([amino.position - position for amino in amin.connected])
            print("step", connected_step)
            aset = set([tuple(x) for x in connected_step])
            bset = set([tuple(x) for x in direc_list])
            free_direc = np.array([x for x in bset - aset])
            for direc in free_direc:
                destination_array = amin.position + direc
                destination_tuple = tuple(destination_array)
                pos_tuple_list = [tuple(item) for item in pos_array]
                if destination_tuple in pos_tuple_list:
                    amin.neighbour.append(self.amino_class_list[pos_tuple_list.index(destination_tuple)])
            sum_conn = np.sum(connected_step, axis=0)

            print("sum", np.sqrt(sum_conn[0]**2 + sum_conn[1]**2))
            print("c1", amin.position + connected_step[0], connected_step[0], amin.position)
            if index > 0 and index < 29:
                print("c2", amin.position + connected_step[1], connected_step[1], amin.position)
            if np.sqrt(sum_conn[0]**2 + sum_conn[1]**2) > 1:
                jump_location = tuple(amin.position + connected_step[0] + connected_step[1])
                print("jl", jump_location)
                if jump_location not in pos_tuple_list:
                    amin.possible_jumps.append(jump_location)
            elif np.sqrt(sum_conn[0]**2 + sum_conn[1]**2) == 1:
                dings = tuple(np.sum(connected_step, axis=0))
                jump_1, jump_2 = (0, 0), (0, 0)
                if dings[0] == 0:
                    jump_1 = (position[0] + 1, position[1] + dings[1])
                    jump_2 = (position[0] - 1, position[1] + dings[1])
                if dings[1] == 0:
                    jump_1 = (position[0] + dings[0], position[1] + 1)
                    jump_2 = (position[0] + dings[0], position[1] - 1)
                print("j1", jump_1, "j2", jump_2, index)
                before = random.random()
                if before > 0.5:
                    if jump_1 not in pos_tuple_list:
                        amin.possible_jumps.append(jump_1)
                    if jump_2 not in pos_tuple_list:
                        amin.possible_jumps.append(jump_2)
                else:
                    if jump_2 not in pos_tuple_list:
                        amin.possible_jumps.append(jump_2)
                    if jump_1 not in pos_tuple_list:
                        amin.possible_jumps.append(jump_1)
            print("jumps", amin.possible_jumps)
"""
"""
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
"""
"""
def check_connections(amin):
    for connection in amin.connected:
        Abstand_tuple = (amin.position[0]-connection.position[0],amin.position[1]-connection.position[1])
        direc_list = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        if Abstand_tuple not in direc_list:
            return False
    return True
"""