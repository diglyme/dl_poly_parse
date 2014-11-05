#!/usr/bin/env python3

# dl_poly_conv
# Converts DL_POLY REVCON file to pdb or xyz

import os

REVCON = "REVCON"
DIR = os.getcwd().split("/")[-1] # by default, use directory name as file name
# TYPES

def getLines():
	with open(REVCON, "r") as f:
		lines = f.readlines()
	return lines

def getAtoms(lines):
	# returns atoms as a list of dicts with id, atom type, x, y z
	atoms = []
	for i, l in enumerate(lines):
		if len(l.split()) == 2 and len(lines[i+1].split()) == 3:
			l1 = l.split()
			l2 = lines[i+1].split()
			# Adds a dict for each atom to list of atoms
			atoms.append(dict(zip(["index", "atom", "x", "y", "z"], [l1[1], l1[0], l2[0], l2[1], l2[2]])))

	return atoms

def xyz():
	atoms = getAtoms(getLines())

	# xyz file begins with number of atoms and a comment line
	xyzoutput = len(atoms) + "\n" + DIR + "\n"
	
	for atom in atoms:
		atom["atom"] = atom["atom"][0].capitalize()
		xyzoutput += " ".join([atom["atom"], atom["x"], atom["y"], atom["z"]) + "\n"

	with open(DIR+".xyz", "w") as f:
		f.write(xyzoutput)



def pdb():
	atoms = getAtoms(getLines())

	pdboutput = ""

	for atom in atoms:
		atom["atom"] = atom["atom"][0].capitalize()
		# round coords to 4dp so they will be 8 chars including -xx.
		# atom["x"] = str(round(float(atom["x"]), 5))[8:]
		# atom["y"] = str(round(float(atom["y"]), 5))[8:]
		# atom["z"] = str(round(float(atom["z"]), 5))[8:]
		pdboutput += "HETATM%5s %-4s UNK          %8s%8s%8s 1.000 0.000          %-2s\n" % (atom["index"], atom["atom"], atom["x"][:8], atom["y"][:8], atom["z"][:8], atom["atom"])

	with open(DIR+".pdb", "w") as f:
		f.write(pdboutput)


if __name__ == "__main__":
	pdb()
