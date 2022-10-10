# About
* This code can be used to calculate dipole-dipole reorientational correlation time for interfacial water(IW) dipole.
* IW are considered as those water molecule which are within 6 Angstorm from the protein backbone atom.
* Dipole is defined as vector joining the central oxygen atom of water molecule and specific hydrogen atom attached with oxygen
atom.


# Required input file
* A sample input file of a single frame is given as input-orient.pdb
* Arrange atom according to this pdb file using gmx trjconv or any other trajectory processing tools.

# Compilation
To run the code:
* gfortan orient.f95 -o a.out
* ./a.out
