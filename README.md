# Finite Element Method
## Description
This Project is a C++ implementation of the Finite Element Method used to simulate 2D heat conduction in a solid body. The program calculates how temperature changes over time across a mesh of nodes and elements, considering material properties and environmental conditions.

The solver handles transient time-dependent heat flow and accounts for convection on the boundaries.

## Key Features
*   **2D FEM Solver**: Uses 4-node quadrilateral elements.
*   **Numerical Integration**: Support for Gaussian integration (2-point, 3-point, and 4-point rules).
*   **Boundary Conditions**: Implementation of the 3rd type boundary condition (convection).
*   **Transient Analysis**: Simulates temperature changes step-by-step using an implicit integration scheme.
*   **Custom Equation Solver**: Includes a built-in Gaussian elimination algorithm to solve the global system of equations.


## How It Works
1.  **Data Input**: The program reads simulation data from a file named `dane.txt`. This includes:
    *   Material properties (conductivity, density, specific heat).
    *   Environment data (ambient temperature, convection coefficient alpha).
    *   Mesh data (coordinates of nodes, element connections, and boundary flags).
2.  **Matrix Generation**: 
    *   For each element, the program calculates the H Matrix (conductivity) and C Matrix (heat capacity).
    *   It uses a Jacobian matrix to transform local element coordinates to global real-world coordinates.
    *   Boundary conditions are applied to the edges marked as "BC".
3.  **Global Assembly**: Local matrices are gathered into a large Global Matrix and a Global Vector.
4.  **Time Iteration**: In every time step, the solver calculates the new temperature distribution based on the previous state and material properties.
5.  **Output**: The program prints the simulation time and the minimum/maximum temperatures in the mesh for each step.

## Input File Format (`dane.txt`)
The program expects a specific text file format containing:
*   `SimulationTime`, `TimeStep`
*   `Conductivity`, `Alfa`, `Tot` (Ambient Temp), `InitialTemp`
*   `Density`, `SpecificHeat`
*   `Nodes number`, `Elements number`
*   Lists of node coordinates and element definitions.


## Example input file &  output

SimulationTime 500

SimulationStepTime 50

Conductivity 25

Alfa 300

Tot 1200

InitialTemp 100

Density 7800

SpecificHeat 700

Nodes number 16

Elements number 9


*Node

   1,  0.100000001, 0.00499999989
   
   2, 0.0666666701, 0.00499999989
   
   3, 0.0333333351, 0.00499999989
   
   4,           0., 0.00499999989
   
   5,  0.100000001, -0.0283333343
   
   6, 0.0666666701, -0.0283333343
   
   7, 0.0333333351, -0.0283333343
   
   8,           0., -0.0283333343
   
   9,  0.100000001, -0.0616666675
   
   10, 0.0666666701, -0.0616666675
   
   11, 0.0333333351, -0.0616666675
   
   12,           0., -0.0616666675
   
   13,  0.100000001, -0.0949999988  
   
   14, 0.0666666701, -0.0949999988
   
   15, 0.0333333351, -0.0949999988
   
   16,           0., -0.0949999988
   
     
*Element, type=DC2D4

 1,  1,  2,  6,  5
 
 2,  2,  3,  7,  6
 
 3,  3,  4,  8,  7
 
 4,  5,  6, 10,  9
 
 5,  6,  7, 11, 10
 
 6,  7,  8, 12, 11
 
 7,  9, 10, 14, 13
 
 8, 10, 11, 15, 14
 
 9, 11, 12, 16, 15
 
 
*BC

1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16

Output:
<img width="428" height="186" alt="image" src="https://github.com/user-attachments/assets/4f36546c-2a5d-4ed4-bbb8-44aac4c817e5" />

## Future Improvements

I plan to refactor the code to improve its structure, safety, and maintainability:
* **Code Refactoring**: Transform the current procedural approach into a more object-oriented structure for better readability and scalability.
*  **Modern C++ Containers**: Replace raw pointers and manual memory management (`new`/`delete`) with `std::vector`. This will improve memory safety and reduce the risk of leaks.
