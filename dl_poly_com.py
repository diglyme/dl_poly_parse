#!/usr/bin/env python3

# dl_poly_com
# Follows the guest molecule by computing CoM of guest and cages
# For each timestep, will print CoM of guest and cage that it is in
# This can then be used to count frequency of visited cages

from numpy import array

# Input and output files
HISTORY = "HISTORY"
OUT = "CoM.txt"

# Define atomic masses
H = 1
C = 12
N = 14
O = 16

# Number of atoms per cage
CC3 = 168

def getLines():
    with open(HISTORY, "r") as f:
        lines = f.readlines()
    return lines

# def cageCoMs(lines):
# 	for i, l in enumerate(lines):

# 	return centres_list

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

def main():
	lines = getLines()
	output = []

	cage_atoms = CC3 * 8
	guest_atoms = int(lines[1].split()[2]) - cage_atoms

	# each atom takes 2 lines, plus 4 header lines per timestep
	step_lines = (cage_atoms + guest_atoms) * 2 + 4

	for i, line in enumerate(lines):
		if (i-2) % step_lines == 0: # first timstep comes after two header lines
			step = int(line.split()[1])

			x, y, z = centre_of_mass(lines[i+4+cage_atoms*2:i+step_lines]) # gets guest CoM
			output.append(str(step)+" "+str(x)+" "+str(y)+" "+str(z)+"\n")

	with open(OUT, "w") as f:
		for line in output:
			f.write(line)

if __name__ == "__main__":
	main()
