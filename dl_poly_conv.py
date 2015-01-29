#!/usr/bin/env python3

# dl_poly_conv
# Converts DL_POLY REVCON file to pdb or xyz

import os, sys

REVCON = "REVCON"
DIR = os.getcwd().split("/")[-1] # by default, use directory name as file name
# TYPES

def get_atoms():
    # returns atoms as a list of dicts with id, atom type, x, y z
    atoms = []
    with open(REVCON, "r") as _f:
        lines = _f.readlines()
        for i, line in enumerate(lines):
            if len(line) == 43:
                l1 = line.split()
                l2 = lines[i+1].split()
                # Adds a dict for each atom to list of atoms
                atoms.append(dict(zip(["index", "atom", "x", "y", "z"],
                                      [l1[1], l1[0], l2[0], l2[1], l2[2]])))

    return atoms

def xyz():
    atoms = get_atoms()

    # xyz file begins with number of atoms and a comment line
    xyz_output = "%i\n%s\n" % (len(atoms), DIR)

    for atom in atoms:
        atom["atom"] = atom["atom"][0].capitalize()
        xyz_output += " ".join([atom["atom"], atom["x"], atom["y"], atom["z"]]) + "\n"

    with open(DIR+".xyz", "w") as _f:
        _f.write(xyz_output)

def pdb():
    atoms = get_atoms()

    pdb_output = ""

    for atom in atoms:
        atom["atom"] = atom["atom"][0].capitalize()
        # round coords to 4dp so they will be 8 chars including -xx.
        # atom["x"] = str(round(float(atom["x"]), 5))[8:]
        # atom["y"] = str(round(float(atom["y"]), 5))[8:]
        # atom["z"] = str(round(float(atom["z"]), 5))[8:]
        pdb_output += "HETATM%5s %-4s UNK          %8s%8s%8s 1.000 0.000          %-2s\n" % (atom["index"], atom["atom"], atom["x"][:8], atom["y"][:8], atom["z"][:8], atom["atom"])

    with open(DIR+".pdb", "w") as _f:
        _f.write(pdb_output)

def main():
    args = sys.argv

    if "xyz" in args:
        xyz()
        return 0

    if "pdb" in args:
        pdb()
        return 0

    print("Please use 'xyz' or 'pdb' options to choose output format.")
    return 1

if __name__ == "__main__":
    main()
