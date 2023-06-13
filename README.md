This set of utilites were created back in 2021 to create a matrix representation of the values of the laplacian of the atomic graph (The critical point of the electronic laplacian ($\nabla^2 \rho (r)$). The atomic graph is an Euler topoligical object that is made up of vertices (Yellow), Edges (Green) and Faces(Pink) (See picture below) and therefore fulfill Euler formula V+F-E= 2.
Below is an example of the atomic graph from the Hexaaquacopper(II):
![copper_aqua](figures/0_M.jpg?raw=true)
![atomic_graph](figures/cu_m.png?raw=true)

This graph is specif for a certain atom in a defined oxidation state. The script AG.py contains two classes: Ag and CP, the first is to fully represent the atomic graph of an atom and the second represents the critical points within an atom. The main() function of this script will create an Ag object that will have a square matrix of # of vertices with the values of $\nabla^2 \rho (r)$ of the connected edges.
This script is called from terminal by:

python AG.py PATH_OF_AGPVIZ (Individual file) 

An example of the result of this script with Metals/cu22.agpviz will be:

CP criteria: [0.5, 0.65] #Criteria to determine if the critical point is in the atomic graph
Euler: True # Satisfies the Euler relationship
[4(4.0),8,6] # [V,E,F]
Atom: cu
[[ 6.58939  0.      12.64917 12.65041]
 [ 0.       6.58941 12.6504  12.64917]
 [12.64917 12.6504   6.58388  0.     ]
 [12.65041 12.64917  0.       6.58386]]
Determinant: -25886.391723965295
AGD: -141.763086998

The script can also get the information of folders with a set of many .agpviz files. To do so, one must import the multiple_files function that takes two parameters:

multiple_files(PATH, name) # Here name represents the name that the user wants to give to an existing file, one can also include the path where this file should be written.

The outpout of this action is file with four columns: # of file, path of the file, metal, "[F,E,V]", "AGD"

The user of this script can get more information about the atomic graph and its components in the book chapter: https://doi.org/10.1016/B978-0-323-90891-7.00018-9
from the book, QTAIM and beyond.
