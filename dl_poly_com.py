#!/usr/bin/env python3

# dl_poly_com
# Follows the guest molecule by computing CoM of guest and cages
# For each timestep, will print CoM of guest and cage that it is in
# This can then be used to count frequency of visited cages

# To do:
#  * Calculate cage CoMs at each step, compare time difference
#  * Use file.readline(), file.tell() and file.seek() instead of reading all lines to conserve memory

from numpy import array, where, square, sqrt

# Input and output files
HISTORY = "HISTORY"
OUT = "CoM_approx.txt"

# Number of atoms per cage
CC3 = 168
# Cages per unit cell
CAGES = 8

def getLines():
    with open(HISTORY, "r") as f:
        lines = f.readlines()
    return lines

def centre_of_mass(lines):
	mass = []
	x = []
	y = []
	z = []
	
	for i, line in enumerate(lines):
		if len(line.split()) == 4:
			mass.append(float(line.split()[2]))
			x.append(float(lines[i+1].split()[0]))
			y.append(float(lines[i+1].split()[1]))
			z.append(float(lines[i+1].split()[2]))

	com_x = sum((array(mass)*array(x)))/sum(mass)
	com_y = sum((array(mass)*array(y)))/sum(mass)
	com_z = sum((array(mass)*array(z)))/sum(mass)

	return com_x, com_y, com_z

def cage_centres(lines):
	cages_x = []
	cages_y = []
	cages_z = []

	for i in range(CAGES):
		x, y, z = centre_of_mass(lines[i*CC3:(i*CC3)+CC3])
		cages_x.append(x)
		cages_y.append(y)
		cages_z.append(z)

	return array(cages_x), array(cages_y), array(cages_z)

def in_cage(guest_x, guest_y, guest_z, cages_x, cages_y, cages_z):
	distances = sqrt(square(cages_x-guest_x) + square(cages_y-guest_y) + square(cages_z-guest_z))
	index = where(distances==min(distances))
	return index[0][0] + 1

def main():
	lines = getLines()
	output = []

	cage_atoms = CC3 * CAGES
	guest_atoms = int(lines[1].split()[2]) - cage_atoms

	# each atom takes 2 lines, plus 4 header lines per timestep
	step_lines = (cage_atoms + guest_atoms) * 2 + 4

	cages_x, cages_y, cages_z = cage_centres(lines[6:6+cage_atoms*2])

	for i, line in enumerate(lines):
		if (i-2) % step_lines == 0: # first timstep comes after two header lines
			step = int(line.split()[1])
			guest_x, guest_y, guest_z = centre_of_mass(lines[i+4+cage_atoms*2:i+step_lines]) # gets guest CoM
			current_cage = in_cage(guest_x, guest_y, guest_z, cages_x, cages_y, cages_z)
			output.append(str(step)+" "+str(round(guest_x, 4))+" "+str(round(guest_y, 4))+" "+str(round(guest_z, 4))+" "+str(current_cage)+"\n")

	with open(OUT, "w") as f:
		for line in output:
			f.write(line)

if __name__ == "__main__":
	main()
