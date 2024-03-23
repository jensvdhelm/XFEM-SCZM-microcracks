# XFEM-SCZM-microcracks

This repository is the official implementation of the following paper.

Paper title: [Limiting microcracks and hydrogen permeability in thermoplastic composites for LH2 storage](DOI will follow once thesis is uploaded to TU Delft repository)




## Description


The model is based on a combined extended finite element method (XFEM) for microcrack initiation and propagation and a surface cohesive zone model (SCZM) for delamination. Conventional FEM requires a mesh that conforms to the geometry of the crack, which means that a crack is modelled along the edge or surface of the element. Besides, the mesh needs to be constantly updated to adapt to crack growth. A Weibull distribution is used to simulate random microcracks.

As a delamination is an inter-ply failure mode, XFEM cannot be used solely here? this is also because XFEM is limited to one crack surface per element. The solution that is used is SCZM methodology, this method uses (additional) predefined delamination surfaces and in that way makes interaction between adjacent plies and cracks possible. The way this is implemented in the model is by using zero-thickness cohesive elements in between the existing mesh elements.

Abaqus is chosen as the software to create the model because Abaqus provides the option to use Python scripting? python is therefore used to create the input file for the analysis. Within the mesoscale model, the section that is used for the finite element analysis is 10 x 10 mm with ply thicknesses dependent on the different testcases, but maximized at a total laminate thickness of 2.38 mm (23 plies). The model is specific for TC1225

Caveats and assumptions: Model is not suitable for cryogenic cycling, the enrichment mesh generation distance is assumed, ‘beta’ for the in-situ shear strength is assumed, permeability is only implemented by means of Darcy’s law, DSMC method for flow mechanism is not implemented yet. 

After the Abaqus model is generated, the following steps should still be implemented to use the Abaqus results in a second part of the Python code:
1. Read nodal connectivity from mesh. 
2. Read x, y, z nodal displacements from output database. 
3. Cross reference connectivity with cracked XFEM elements. 
4. DCOD calculation based on relative x, y, z displacements of adjacent nodes in crack elements. 
5. Calculation of crack-overlap area for individual crack networks. 
6. Sum over the entire laminate and calculate permeability 



## History

A changelog is not applicable



## Authors or Maintainers

    Jens van der Helm ([@jensvdhelm]( https://github.com/jensvdhelm), jensvanderhelm@live.nl, TU Delft (student)
    
	Arshdeep Singh Brar ([@arshdeep-brar]( https://github.com/arshdeep-brar), TU Delft (student)   



## Table of Contents

The Python code is provided in one .py file which scripts the Abaqus model. The code encompasses the following parts:
1. Importing Abaqus modules
2. Mesh parameters and material properties
3. Calculation of shear strength of the matrix per ply
4. Implementing permeability variables
5. Defining the Knudsen number
6. Initiating the sketch of the frame in Abaqus incl. boundary conditions and delamination zones
7. Assigning to each element: material, conductivity, CTE, specific heat, permeability
8. Assembly of the model incl. loading (coupled temperature displacement step)
9. XFEM assignment
10. Implementing desired temperature, heat flux and convection boundary conditions
11. Write Job input for Abaqus



## Requirements  

[![Spyder (Python 3.8)](https://img.shields.io/badge/Python-3.8-3776AB)](https://www.python.org/downloads/release/python-380/)

[![Abaqus CAE]( https://img.shields.io/badge/Abaqus-V.2020-blue?logo=Dassault%20Syst%C3%A8mes)]( https://www.3ds.com/edu/education/students/solutions/abaqus-le )




## Structure

1. Clean_integrated model.py
2. README_code.md
3. LICENSE



## License

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)  


    The contents of this repository are licensed under a **Apache License 2.0** license (see LICENSE file).



## References

N/A



## Citation


        If you want to cite this repository in your research paper, please use the following information:   
        
Reference: [Connecting 4TU.ResearchData with Git](https://data.4tu.nl/info/about-your-data/getting-started)   

Reference: [DOI 4TU.ResearchData](10.4121/4303d493-d687-4ded-981b-b1714619097a)   




## Would you like to contribute?

If you would like to contribute to this project, you can fork the project.

List of contributors:
1. Arshdeep Singh Brar (Input parameters, Abaqus model basis, delamination and XFEM assignment)
2. Jens van der Helm (Modification for TC1225, Knudsen number, in-situ shear strength and permeability)  
