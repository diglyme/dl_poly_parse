#!/usr/bin/env python3

from numpy import array, empty, fromiter, where, square, sqrt, concatenate, average
import pdb

# input and output files
HISTORY = "HISTORY"
OUT = "CoM.txt"

# Number of atoms per cage
CAGE_ATOMS = 168
# Cages per cell
CAGES = 8

WINDOW = [(2, 100, 142), (16, 72, 128), (30, 58, 156), (44, 86, 114)]

VDW = {"h": 1.2, "c": 1.7, "n": 1.55}

GUEST_ATOMS = 1
GUESTS = 1

class Cage:
    def __init__(self, n, x, y, z):
        self.n = n
        self.x = x
        self.y = y
        self.z = z
        self.atoms = len(x)

        # create temporary coordinates that can exceed unit cell
        self.per_x = x
        self.per_y = y
        self.per_z = z
        for i in range(self.atoms-1):
            xdiff = self.x[i+1] - self.x[0]
            ydiff = self.y[i+1] - self.y[0]
            zdiff = self.z[i+1] - self.z[0]

            # if distances between atom 0 and atom i are physically unreasonable
            # change coordinates of atom i to bring it close to atom 0
            if abs(xdiff) > box[self.n]/2 and xdiff < 0:
                self.per_x[i+1] += box[self.n]
            if abs(xdiff) > box[self.n]/2 and xdiff > 0:
                self.per_x[i+1] -= box[self.n]

            if abs(ydiff) > box[self.n]/2 and ydiff < 0:
                self.per_y[i+1] += box[self.n]
            if abs(ydiff) > box[self.n]/2 and ydiff > 0:
                self.per_y[i+1] -= box[self.n]

            if abs(zdiff) > box[self.n]/2 and zdiff < 0:
                self.per_z[i+1] += box[self.n]
            if abs(zdiff) > box[self.n]/2 and zdiff > 0:
                self.per_z[i+1] -= box[self.n]

    def per_position(self, index):
        return (self.per_x[index], self.per_y[index], self.per_z[index])

    def periodic_coords(self):
        pass
        
    def centre_of_mass(self, periodic=False):
        # use moved coordinates to calculate centre of mass
        com_x = sum(cage_mass*self.per_x)/sum(cage_mass)
        com_y = sum(cage_mass*self.per_y)/sum(cage_mass)
        com_z = sum(cage_mass*self.per_z)/sum(cage_mass)

        if periodic:
            return (com_x, com_y, com_z)

        # some CoMs will now be outside unit cell
        # hopefully this will return them to correct centre
        if com_x > box[self.n]/2:
            com_x -= box[self.n]
        elif com_x < -box[self.n]/2:
            com_x += box[self.n]

        if com_y > box[self.n]/2:
            com_y -= box[self.n]
        elif com_y < -box[self.n]/2:
            com_y += box[self.n]

        if com_z > box[self.n]/2:
            com_z -= box[self.n]
        elif com_z < -box[self.n]/2:
            com_z += box[self.n]

        return (com_x, com_y, com_z)

    def pore_radius(self):
        com = self.centre_of_mass(periodic=True)
        distances = [ distance(self.per_position(i), com) - VDW[cage_type[i][0]] for i in range(self.atoms) ]
        return min(distances)

    def window_radii(self):
        windows = []
        for w in WINDOW:
            A, B, C = w
            #pdb.set_trace()
            a = distance(self.per_position(B), self.per_position(A))
            b = distance(self.per_position(C), self.per_position(B))
            c = distance(self.per_position(A), self.per_position(C))

            radius = (a*b*c) / sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c))
            windows.append(radius)
        return windows


    def largest_window(self):
        return max(self.window_radii())

    def smallest_window(self):
        return min(self.window_radii())

class Guest:
    def __init__(self, n, x, y, z, types, masses):
        self.n = n
        self.x = x
        self.y = y
        self.z = z
        self.type = types
        self.mass = masses
        self.atoms = len(x)

    def centre_of_mass(self):
        if self.atoms == 1:
            com_x = sum(self.mass*self.x)/sum(self.mass)
            com_y = sum(self.mass*self.y)/sum(self.mass)
            com_z = sum(self.mass*self.z)/sum(self.mass)
            return (com_x, com_y, com_z)

        else:
            # create temporary coordinates that can exceed unit cell for CoM calculation
            tmp_x = self.x
            tmp_y = self.y
            tmp_z = self.z
            for i in range(self.atoms-1):
                xdiff = self.x[i+1] - self.x[0]
                ydiff = self.y[i+1] - self.y[0]
                zdiff = self.z[i+1] - self.z[0]

                # if distances between atom 0 and atom i are physically unreasonable
                # change coordinates of atom i to bring it close to atom 0
                if abs(xdiff) > box[self.n]/2 and xdiff < 0:
                    tmp_x[i+1] += box[self.n]
                if abs(xdiff) > box[self.n]/2 and xdiff > 0:
                    tmp_x[i+1] -= box[self.n]

                if abs(ydiff) > box[self.n]/2 and ydiff < 0:
                    tmp_y[i+1] += box[self.n]
                if abs(ydiff) > box[self.n]/2 and ydiff > 0:
                    tmp_y[i+1] -= box[self.n]

                if abs(zdiff) > box[self.n]/2 and zdiff < 0:
                    tmp_z[i+1] += box[self.n]
                if abs(zdiff) > box[self.n]/2 and zdiff > 0:
                    tmp_z[i+1] -= box[self.n]

            # use moved coordinates to calculate centre of mass
            com_x = sum(self.mass*tmp_x)/sum(self.mass)
            com_y = sum(self.mass*tmp_y)/sum(self.mass)
            com_z = sum(self.mass*tmp_z)/sum(self.mass)

            # some CoMs will now be outside unit cell
            # hopefully this will return them to correct centre
            if com_x > box[self.n]/2:
                com_x -= box[self.n]
            elif com_x < -box[self.n]/2:
                com_x += box[self.n]

            if com_y > box[self.n]/2:
                com_y -= box[self.n]
            elif com_y < -box[self.n]/2:
                com_y += box[self.n]

            if com_z > box[self.n]/2:
                com_z -= box[self.n]
            elif com_z < -box[self.n]/2:
                com_z += box[self.n]

            return (com_x, com_y, com_z)

def get_lines():
    print("Reading HISTORY file... ", end="", flush=True)
    with open(HISTORY, "r") as histfile:
        lines = histfile.readlines()
    print("done.")
    return lines

def pull_data():
    global box, step, cage_type, cage_mass
    timestep = 0
    init = True
    atom = False
    atom_num = 0
    step = []
    atom_type = []
    atom_mass = []
    box = []
    atoms_per = CAGE_ATOMS * CAGES + GUEST_ATOMS * GUESTS

    # print("Reading HISTORY file... ", end="", flush=True)
    # with open(HISTORY, "r") as histfile:
    #     lines = histfile.readlines()
    # print("done!")
    # num_lines = len(lines)
    with open(HISTORY, "r") as histfile:
        num_lines = sum(1 for l in histfile)
        histfile.seek(0)
        tot_atoms = sum(1 for l in histfile if len(l) == 43)
        histfile.seek(0)

        x = empty(tot_atoms, "float")
        y = empty(tot_atoms, "float")
        z = empty(tot_atoms, "float")

        print("Extracting data from file...")
        for i, line in enumerate(histfile):
            #print(line)
            l = line.split()

            # if len(step) > 0:
            #     if step[-1] == 86000: debug = True

            if "timestep" in l:
                timestep += 1
                step.append(int(l[1]))
                if timestep > 1: init = False

            elif len(l) == 4:
                atom = True
                if init:
                    atom_type.append(l[0])
                    atom_mass.append(float(l[2]))

            elif len(l) == 3 and atom:
                x[atom_num] = float(l[0])
                y[atom_num] = float(l[1])
                z[atom_num] = float(l[2])
                atom_num += 1
                atom = False

            elif len(l) == 3 and len(box) < timestep:
                box.append(float(l[0]))

            else: atom = False

            debug = False

            done = str(round(((i+1)/num_lines)*100, 1))
            print(done+"% done\r", end="", flush=True)


    # print("Cleaning up raw file data.")
    # lines = None

    print("Creating cage and guest objects...")
    cage_type = atom_type[:CAGE_ATOMS]
    cage_mass = fromiter(atom_mass[:CAGE_ATOMS], "float", CAGE_ATOMS)
    guest_type = atom_type[CAGE_ATOMS*CAGES:CAGE_ATOMS*CAGES+GUEST_ATOMS]
    guest_mass = fromiter(atom_mass[CAGE_ATOMS*CAGES:CAGE_ATOMS*CAGES+GUEST_ATOMS], "float", GUEST_ATOMS)

    frame = []
    steps = len(step)

    for i in range(steps):
        start = i*atoms_per
        end = (i+1)*atoms_per

        cages = []
        guests = []

        for j in range(CAGES):
            cage_start = start + j * CAGE_ATOMS
            cage_end = start + (j+1) * CAGE_ATOMS
            cages.append(Cage(i, x[cage_start:cage_end], y[cage_start:cage_end], z[cage_start:cage_end]))

        for j in range(GUESTS):
            guest_start = start + (CAGE_ATOMS+CAGES) + j * GUEST_ATOMS
            guest_end = start + (CAGE_ATOMS+CAGES) + (j+1) * GUEST_ATOMS
            guests.append(Guest(i, x[guest_start:guest_end], y[guest_start:guest_end], z[guest_start:guest_end], guest_type, guest_mass))

        frame.append({"cage": cages, "guest": guests})

        done = str(round(((i+1)/steps)*100, 1))
        print(done+"% done\r", end="", flush=True)

    print("All data extracted!")
    return frame

def distance(coords1, coords2):
    # expand tuples to variables
    x1, y1, z1 = coords1
    x2, y2, z2 = coords2
    return sqrt(square(x2-x1) + square(y2-y1) + square(z2-z1))

def visualise(frame, begin, stop):
    for i, f in enumerate(frame[begin:stop]):
        cages_com = [ c.centre_of_mass() for c in f["cage"] ]
        towrite = str(CAGES) + "\nCentres of mass " + str(i) + "\n"

        for com in cages_com:
            x, y, z = com
            towrite += " ".join(["Fe", str(x), str(y), str(z)]) + "\n"

        with open("centres.xyz", "a") as centresfile:
            centresfile.write(towrite)

def msd(frame, g=0):
    disps = []
    for i, f in enumerate(frame[1:]):
        curr_x, curr_y, curr_z = f["guest"][g].centre_of_mass()
        prev_x, prev_y, prev_z = frame[i-1]["guest"][g].centre_of_mass()

        xdiff = curr_x - prev_x
        ydiff = curr_y - prev_y
        zdiff = curr_z - prev_z

        if abs(xdiff) > box[i-1]/2 and xdiff < 0:
            curr_x += box[i-1]
        if abs(xdiff) > box[i-1]/2 and xdiff > 0:
            curr_x -= box[i-1]

        if abs(ydiff) > box[i-1]/2 and ydiff < 0:
            curr_y += box[i]
        if abs(ydiff) > box[i-1]/2 and ydiff > 0:
            curr_y -= box[i-1]

        if abs(zdiff) > box[i-1]/2 and zdiff < 0:
            curr_z += box[i-1]
        if abs(zdiff) > box[i-1]/2 and zdiff > 0:
            curr_z -= box[i-1]

        disps.append(square(curr_x-prev_x) + square(curr_y-prev_y) + square(curr_z-prev_z))

    tot_disps = [sum(disps[:i+1]) for i in range(len(disps))]
    return [average(tot_disps[:i+1]) for i in range(len(tot_disps))]

def in_cage(frame, g=0):
    in_cage = []
    for i, f in enumerate(frame):
        # calculate com for each cage in frame
        cages_com = [ c.centre_of_mass(periodic=True) for c in f["cage"] ]
        guest_com = f["guest"][g].centre_of_mass()

        # square_displacement.append(square(distance(guest_start, guest_com)))
        # msd = sum(square_displacement) / len(square_displacement)

        distances = [ distance(guest_com, com) for com in cages_com ]
        # guest is in cage with minimum distance to centre of mass
        # this gives cage number, counting from 1 rather than 0
        in_cage.append(distances.index(min(distances)) + 1)

    return in_cage

def guest_centres(frame, g=0):
    return [f["guest"][g].centre_of_mass() for f in frame]

def main():
    print("Number of guest molecules: ", end="", flush=True)
    GUESTS = int(input())
    print("Number of atoms per guest: ", end="", flush=True)
    GUEST_ATOMS = int(input())

    if GUESTS != 1 or GUEST_ATOMS < 1:
        print("Main function of this program only prints centre of mass of a single guest molecule.")
        return 1

    frame = pull_data()
    steps = len(frame)
    guest_start = frame[0]["guest"][0].centre_of_mass()
    # square_displacement = []

    with open(OUT, "w") as output:
        output.write("") # clear existing file

    print("\nCalculating centres of mass:")
    
    with open(OUT, "a") as output:
        for i, f in enumerate(frame):
            # calculate com for each cage in frame
            cages_com = [ c.centre_of_mass(periodic=True) for c in f["cage"] ]
            guest_com = f["guest"][0].centre_of_mass()

            # square_displacement.append(square(distance(guest_start, guest_com)))
            # msd = sum(square_displacement) / len(square_displacement)

            distances = [ distance(guest_com, com) for com in cages_com ]
            # guest is in cage with minimum distance to centre of mass
            # this gives cage number, counting from 1 rather than 0
            in_cage = distances.index(min(distances)) + 1

            x, y, z = guest_com

            towrite = " ".join([str(step[i]), str(x), str(y), str(z), str(in_cage)]) + "\n"

            #with open(OUT, "a") as output:
            output.write(towrite)

            done = str(round(((i+1)/steps)*100, 1))
            print(done+"% done\r", end="", flush=True)

    return 0

if __name__ == "__main__":
    main()