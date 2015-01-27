#!/usr/bin/env python3

from numpy import array, empty, fromiter, where, square, sqrt, concatenate, average
from progressbar import ProgressBar
import pdb
import argparse

# input and output files
HISTORY = "HISTORY"
OUT = "CoM.txt"

# Number of atoms per cage
CAGE_ATOMS = 168
# Cages per cell
CAGES = 8

WINDOW = [(2, 100, 142), (16, 72, 128), (30, 58, 156), (44, 86, 114)]

VDW = {"h": 1.2, "c": 1.7, "n": 1.55}

guest_atoms = 1
guest_num = 1

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
            return (float(self.x), float(self.y), float(self.z))

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

def pull_data(guest_num, guest_atoms, is_guest):
    global steps, box, step, cage_type, cage_mass
    timestep = 0
    init = True
    atom = False
    atom_num = 0
    step = []
    atom_type = []
    atom_mass = []
    box = []
    atoms_per = CAGE_ATOMS * CAGES + guest_atoms * guest_num

    # print("Reading HISTORY file... ", end="", flush=True)
    # with open(HISTORY, "r") as histfile:
    #     lines = histfile.readlines()
    # print("done!")
    # num_lines = len(lines)
    with open(HISTORY, "r") as histfile:
        print("Counting size of file... ", end="", flush=True)
        num_lines = sum(1 for l in histfile)
        histfile.seek(0)
        print("done! (%i lines)" % (num_lines))

        print("Counting number of atoms... ", end="", flush=True)
        tot_atoms = sum(1 for l in histfile if len(l) == 43)
        histfile.seek(0)
        print("done! (%i atoms)" % (tot_atoms))

        x = empty(tot_atoms, "float")
        y = empty(tot_atoms, "float")
        z = empty(tot_atoms, "float")

        print("\nExtracting data from file...")
        filebar = ProgressBar(maxval=num_lines).start()
        for i, line in enumerate(histfile):
            l = line.split()

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

            filebar.update(i+1)

    filebar.finish()
    steps = len(step)


    print("\nCreating cage %sobjects..." % ("and guest "*is_guest))
    objbar = ProgressBar(maxval=steps).start()
    cage_type = atom_type[:CAGE_ATOMS]
    cage_mass = fromiter(atom_mass[:CAGE_ATOMS], "float", CAGE_ATOMS)
    if is_guest:
        guest_type = atom_type[CAGE_ATOMS*CAGES:CAGE_ATOMS*CAGES+guest_atoms]
        guest_mass = fromiter(atom_mass[CAGE_ATOMS*CAGES:CAGE_ATOMS*CAGES+guest_atoms], "float", guest_atoms)

    frame = []

    for i in range(steps):
        start = i*atoms_per
        end = (i+1)*atoms_per

        cages = []

        for j in range(CAGES):
            cage_start = start + j * CAGE_ATOMS
            cage_end = start + (j+1) * CAGE_ATOMS
            cages.append(Cage(i, x[cage_start:cage_end], y[cage_start:cage_end], z[cage_start:cage_end]))

        if is_guest:
            guests = []
            for j in range(guest_num):
                guest_start = start + (CAGE_ATOMS*CAGES) + j * guest_atoms
                guest_end = start + (CAGE_ATOMS*CAGES) + (j+1) * guest_atoms
                guests.append(Guest(i, x[guest_start:guest_end], y[guest_start:guest_end], z[guest_start:guest_end], guest_type, guest_mass))

            frame.append({"cage": cages, "guest": guests})
        else:
            frame.append({"cage": cages})

        objbar.update(i+1)

    objbar.finish()
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

def get_windows(frame):
    print("\nFinding all window radii...")
    pbar = ProgressBar(maxval=len(frame)).start()
    windows = []
    for i, f in enumerate(frame):
        for c in f["cage"]:
            for w in c.window_radii():
                windows.append(w)
        pbar.update(i+1)
    pbar.finish()
    return windows

def get_pores(frame):
    print("\nFinding all pore radii...")
    pbar = ProgressBar(maxval=len(frame)).start()
    pores = []
    for i, f in enumerate(frame):
        for c in f["cage"]:
            pores.append(c.pore_radius())
        pbar.update(i+1)
    pbar.finish()
    return pores


def msd(frame, begin_at, guest=0): # this doesn't actually give MSD, needs fixing!
    print("\nFinding mean square displacement of guest %i..." % (guest+1))
    pbar = ProgressBar(maxval=len(frame)).start()
    disps = []
    for i, f in enumerate(frame[1:]):
        curr_x, curr_y, curr_z = f["guest"][guest].centre_of_mass()
        prev_x, prev_y, prev_z = frame[i-1]["guest"][guest].centre_of_mass()

        xdiff = curr_x - prev_x
        ydiff = curr_y - prev_y
        zdiff = curr_z - prev_z

        if abs(xdiff) > box[i+begin_at-1]/2 and xdiff < 0:
            curr_x += box[i+begin_at-1]
        if abs(xdiff) > box[i+begin_at-1]/2 and xdiff > 0:
            curr_x -= box[i+begin_at-1]

        if abs(ydiff) > box[i+begin_at-1]/2 and ydiff < 0:
            curr_y += box[i+begin_at-1]
        if abs(ydiff) > box[i+begin_at-1]/2 and ydiff > 0:
            curr_y -= box[i+begin_at-1]

        if abs(zdiff) > box[i+begin_at-1]/2 and zdiff < 0:
            curr_z += box[i+begin_at-1]
        if abs(zdiff) > box[i+begin_at-1]/2 and zdiff > 0:
            curr_z -= box[i+begin_at-1]

        disps.append(square(curr_x-prev_x) + square(curr_y-prev_y) + square(curr_z-prev_z))
        pbar.update(i+1)

    tot_disps = [sum(disps[:i+1]) for i in range(len(disps))]
    msd = [0.0] + [average(tot_disps[:i+1]) for i in range(len(tot_disps))]
    pbar.finish()
    return msd

def in_cage(frame, guest=0):
    print("\nCalculating cages guest %i has travelled through..." % (guest+1))
    pbar = ProgressBar(maxval=len(frame)).start()
    in_cage = []
    for i, f in enumerate(frame):
        # calculate com for each cage in frame
        cages_com = [ c.centre_of_mass(periodic=True) for c in f["cage"] ]
        guest_com = f["guest"][guest].centre_of_mass()

        # square_displacement.append(square(distance(guest_start, guest_com)))
        # msd = sum(square_displacement) / len(square_displacement)

        distances = [ distance(guest_com, com) for com in cages_com ]
        # guest is in cage with minimum distance to centre of mass
        # this gives cage number, counting from 1 rather than 0
        in_cage.append(distances.index(min(distances)) + 1)

        pbar.update(i+1)

    pbar.finish()
    return in_cage

def guest_centres(frame, guest=0):
    print("\nFinding guest %i centres of mass... " % (guest+1))
    coms = []
    pbar = ProgressBar(maxval=len(frame)).start()
    for i, f in enumerate(frame):
        coms.append(f["guest"][guest].centre_of_mass())
        pbar.update(i+1)
    pbar.finish()
    return coms


def main():
    parser = argparse.ArgumentParser(description="Calculates centres-of-mass and other positional data from a DL_POLY HISTORY file.", epilog="WARNING: currently may do strange things for multiple guests.")
    parser.add_argument("-n", "--guests", type=int, default=0, help="Number of guest molecules (default = 0)")
    parser.add_argument("-g", "--guest_atoms", type=int, default=1, help="Number of atoms in each guest (default = 1)")
    parser.add_argument("-e", "--equilib", type=int, default=0, help="Number of equilibriation steps to skip before calculating (default = 0)")
    parser.add_argument("-o", "--output", default="guest_motion.txt", help="File name for guest motion output (default = 'guest_motion.txt'")
    parser.add_argument("task", nargs="+", help="""Task(s) to run, choose from:
        com (print guest centres of mass);
        cage (print which cage the guest is in);
        msd (print mean square displacement);
        pores (print all pore radii in own file);
        windows (print all window radii in own file).""")
    args = parser.parse_args()

    tasks = [task for task in ["com", "cage", "msd", "pores", "windows"] if task in args.task]
    if len(tasks) == 0:
        print("Error: no valid task selected.\nPlease choose from: com, cage, msd, pores, windows")
        return 1

    # THIS VARIABLE CHANGE FOR TESTING PURPOSES ONLY!
    # tasks = ["com", "cage", "msd", "pores", "windows"]

    if "com" in tasks or "cage" in tasks or "msd" in tasks:
        is_guest = True
        if args.guests == 0  or args.guest_atoms == 0:
            print("Error: com, cage and msd tasks can only be performed when guest molecules are present (did you forget to set the number of guests?)")
            return 1
    else:
        is_guest = False

    frame = pull_data(args.guests, args.guest_atoms, is_guest)

    # find first recorded step after equilibriation
    begin_at = step.index(next(s for s in step if s > args.equilib))

    if "windows" in tasks:
        windows = get_windows(frame[begin_at:])
        print("Writing window radii... ", end="", flush=True)
        with open("windows.txt", "w") as window_file:
            window_file.write("\n".join([str(w) for w in windows]))
        print("done!")

    if "pores" in tasks:
        pores = get_pores(frame[begin_at:])
        print("Writing pore radii... ", end="", flush=True)
        with open("pores.txt", "w") as pore_file:
            pore_file.write("\n".join([str(p) for p in pores]))
        print("done!")

    for g in range(args.guests):
        output_data = [[str(s) for s in step[begin_at:]]]

        if "com" in tasks:
            guest_coms = guest_centres(frame[begin_at:], guest=g)
            str_coms = []
            for com in guest_coms:
                x, y, z = com
                str_coms.append(" ".join([str(x),str(y),str(z)]))
            output_data.append(str_coms)

        if "cage" in tasks:
            output_data.append([str(c) for c in in_cage(frame[begin_at:], guest=g)])

        if "msd" in tasks:
            output_data.append([str(m) for m in msd(frame[begin_at:], begin_at)])

        if len(output_data) > 1:
            if args.guests == 1:
                output = args.output
            else:
                output = "guest_" + str(g+1) + "_motion.txt"

            print("\nWriting output file... ", end="", flush=True)
            towrite = "Step " + "X Y Z "*("com" in tasks) + "In_cage "*("cage" in tasks) + "Mean_Square_Displacement"*("msd" in tasks) + "\n"

            for i in range(len(output_data[0])):
                towrite += " ".join([output_data[j][i] for j in range(len(output_data))]) + "\n"

            with open(output, "w") as output_file:
                output_file.write(towrite)
            print("done!")

    return 0

if __name__ == "__main__":
    main()