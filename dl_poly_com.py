#!/usr/bin/env python3

from numpy import array, where, square, sqrt, concatenate

# input and output files
HISTORY = "HISTORY"
OUT = "CoM.txt"

# Number of atoms per cage
CAGE_SIZE = 168
# Cages per cell
CAGES = 8

cage_atoms = CAGE_SIZE * CAGES

# Number of guest molecules
GUESTS = 1

class Cage:
    def __init__(self, n, x, y, z):
        self.n = n
        self.x = x
        self.y = y
        self.z = z
        self.atoms = len(x)

    def position(self, index):
        return (self.x[index], self.y[index], self.z[index])

    def centre_of_mass(self):
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
        com_x = sum(cage_mass*tmp_x)/sum(cage_mass)
        com_y = sum(cage_mass*tmp_y)/sum(cage_mass)
        com_z = sum(cage_mass*tmp_z)/sum(cage_mass)

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
        pass

    def window_radii(self):
        pass

    def largest_window(self):
        pass

    def smallest_window(self):
        pass

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
    print("done!")
    return lines

def pull_data(lines):
    global box, step, cage_type, cage_mass

    print("Pulling data from file:")

    # list which all coordinate data will go into, each list item is a timestep
    frame = []

    # list of unit cell dimension at each timestep
    box = []

    init = True # for checking if we haven't set up cage_type and cage_mass variables yet
                # for getting cage atom masses on first cage encountered
                # this assumes the atom ordering is identical for each cage!

    # list of timestep values
    step = array([ int(l.split()[1]) for l in lines if "timestep" in l])

    # indices of each time in step list when step is less than previous, discarding the first
    # this gives the index of each restart timestep
    restarts = [ i for i, s in enumerate(step) if s < step[i-1] ][1:]

    for s in restarts:
        # at each restart step, we want all future steps added to final of previous run
        step = concatenate((step[:s], step[s:]+step[s-1]))
    
    start = 1 # counter of how many starts have passed
    tot_steps = len(step)

    tot_atoms = int(lines[1].split()[-1])
    guest_atoms = tot_atoms - cage_atoms
    guest_size = guest_atoms // GUESTS

    for i in range(tot_steps):
        # each frame will be a dict containing list of cage objects and list of guest objects
        cage = []
        guest = []

        if i in restarts:
            start += 1

        start = (i+1)*4 + tot_atoms*i*2 + start*2 # two header lines for each start, 4 for each timstep
        end = start + tot_atoms*2

        # unit cell coordinate data found before atom start
        box.append(float(lines[start-1].split()[2]))

        for j in range(CAGES):
            cage_start = start + j*CAGE_SIZE*2
            cage_end = start + (j+1)*CAGE_SIZE*2

            # every second line contains coordinates
            xs = array([ float(l.split()[0]) for l in lines[cage_start+1:cage_end:2] ])
            ys = array([ float(l.split()[1]) for l in lines[cage_start+1:cage_end:2] ])
            zs = array([ float(l.split()[2]) for l in lines[cage_start+1:cage_end:2] ])

            cage.append(Cage(i, xs, ys, zs)) # create Cage object, add to list of cages

            if init: # get cage atom types and masses if they have not yet been initialised
                cage_type = [ l.split()[0] for l in lines [cage_start:cage_end:2]]
                cage_mass = array([ float(l.split()[2]) for l in lines[cage_start:cage_end:2] ])
                init = False

        for j in range(GUESTS):
            guest_start = start + cage_atoms*2 + j*guest_size*2
            guest_end = start + cage_atoms*2 + (j+1)*guest_size*2

            xs = array([ float(l.split()[0]) for l in lines[guest_start+1:guest_end:2] ])
            ys = array([ float(l.split()[1]) for l in lines[guest_start+1:guest_end:2] ])
            zs = array([ float(l.split()[2]) for l in lines[guest_start+1:guest_end:2] ])

            atom_types = [ l.split()[0] for l in lines[guest_start:guest_end:2] ]
            atom_masses = array([ float(l.split()[2]) for l in lines[guest_start:guest_end:2] ])

            guest.append(Guest(i, xs, ys, zs, atom_types, atom_masses))

        # add lists of cages and guests to list of frames
        # as a dict so they can be called by frame[n]["guest"][m]
        frame.append({"cage": cage, "guest": guest})

        done = str(round((i/tot_steps)*100, 1)) # percentage completed
        print(done+"% done\r", end="", flush=True)

    return frame

def distance(coords1, coords2):
    # expand tuples to variables
    x1, y1, z1 = coords1
    x2, y2, z2 = coords2
    return sqrt(square(x2-x1) + square(y2-y1) + square(z2-z1))

def main(frame):

    if len(frame[0]["guest"]) > 1:
        print("Main function of this program only prints centre of mass of a single guest molecule.")
        return 1

    output = open(OUT, "a")

    for i, f in enumerate(frame):
        # calculate com for each cage in frame
        cages_com = [ c.centre_of_mass() for c in f["cage"] ]
        guest_com = f["guest"][0].centre_of_mass()

        distances = [ distance(guest_com, com) for com in cages_com ]
        # guest is in cage with minimum distance to centre of mass
        # this gives cage number, counting from 1 rather than 0
        in_cage = distances.index(min(distances)) + 1

        x, y, z = guest_com

        towrite = " ".join([str(step[i]), str(x), str(y), str(z), str(in_cage)]) + "\n"

        #with open(OUT, "a") as output:
        output.write(towrite)

    output.close()

    return 0

if __name__ == "__main__":
    main(pull_data(get_lines()))