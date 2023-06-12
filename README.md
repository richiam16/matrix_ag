This set of utilites were created back in 2021 to create a matrix representation of the atomic graph (The critical point of the electronic laplacian $\nabla^2 \rho (r)$). The atomic graph is an Euler topoligical object that is made up of vertices (Yellow), Edges (Green) and Faces(Pink) (See picture below) and therefore fulfill Euler formula V+F-E= 2.
Below is an example of the atomic graph from the Hexaaquacopper(II):
![alt text](figures/0_M.jpg?raw=true)
![alt text](figures/cu_m.png?raw=true)

This graph is specif for a certain atom in a defined oxidation state. The script AG.py contains two classes: Ag and CP, the first is to fully represent the atomic graph of an atom and the second represents the critical points within an atom. The main function of this script will create an Ag object that will have a square matrix of # of vertices with the values of $\nabla^2 \rho (r)$ of the connected edges.
This script is called from terminal by:

python AG.py PATH_OF_AGPVIZ 

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
