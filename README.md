# Projects
##### This repository aims to present the different projects and numerical skills that I have been able to develop since 2017. 
##### Many of the projects presented here have been edited or taken out of their most practical context for this overview.

## `erosion_solver`
### User interface with model selection for erosive wear calculation 
#### 2018 - UFRJ 
- `Python` : cgi, numpy
- `HTML`
- `CSS`
- `JavaScript`

Work that integrates mathematical models for the description of erosive wear (an important problem in the oil and gas industry) to a user interface through CGI (Common Gateaway Interface).
The models were implemented with Python, in addition to the CGI script. HTML, CSS and jQuery were used for basic front-end characterization and operation.

The simulator can be seen working on the link below: 
http://www.escoamentosmultifasicos.coppe.ufrj.br/cgi-bin/erosion_models

## `crow_instability`
### Simulation and stability analysis of Crow Instability 
#### 2021 - M2 Fluid Mechanics at X.
- `Python` : FEniCS, SciPy, Matplotlib
- `gmsh`

This work integrates meshes generated with the gmsh tool with the FEniCS library, used to solve PDEs with Python through finite elements.
At first, the dynamics of a vortex dipole is simulated by means of a linearized Navier-Stokes direct resolution (baseflow.py).
Then, the tools learned in class for spectral analysis of flow stability and optimal perturbations are used (EVP.py).

## `chaos_windows`
### Optimized calculations for the logistic map 
#### (2018-2020) - UFRJ

- `Python` : Numba, numpy, matplotlib, scipy.
- `Mathematica`

In these notebooks I do relevant calculations for a research project focused on discrete chaotic dynamics, with which I was involved for 2 and a half years at the Physics Institute of the Federal University of Rio de Janeiro. 

In summary, the research was focused on a complexified version of the logistic map, in order to study the appearance of periodicity windows in the transition to chaos regime. In these years I learned a lot about code optimization, as the results depended on exhaustive iterative calculations (check average_chaos.ipynb).

All my contribution to this research was in the numerical and analytical scope. In one of these notebooks in particular, I use the Numba library to translate Python functions into optimized machine code, using an LLVM compilation library (lces.ipynb).

Most of my involvement in this project took place in Mathematica notebooks, as it is a high-level language that allows for a more fluid integration between analytical and numerical studies. I had weekly meetings with Professor Sergio Joras, who gave me numerical challenges in Mathematica related to research every week.  The various Mathematica notebooks I've programmed over the years can be found in the `mathematica` folder under `chaos`. 

## `climate_canada`
### Canada climate analysis 
#### 2021 - 2A ENSTA 

- `Python` : Matplotlib, cartopy, xarray.

Work done throughout the MF204 course, going through the main python tools used for climate change research. The xarray library was used to read _NetCDF4_ files containing global temperature information over the years. Additionally, the __cartopy__ library was used for a proper visualization of these data. Several statistical tools were used to obtain predictive climate information.

## `boiling` and `bubble`
### Simulations and tools developed for Basilisk C.
#### 2021 - PRe Interniship and M2 Fluid Mechanics

- `Basilisk C`
- `awk`
- `gnuplot`

Basilisk is a C-based library used to solve PDEs. In a practical way, it is a great tool for fluid dynamics simulations with dynamically adapted mesh.

I had the opportunity to carry out my ENSTA 2A PRe working with the simulation of bubbles in variable gravitational fields. In it, in addition to the simulations, I had to develop several scripts in C to extract dynamic data from the bubble interface.

Subsequently, in the M2 of fluid mechanics, I was in charge of implementing a modification in the program's source code in order to simulate the nucleate boiling phenomenon (proposed by the professor as a challenge, as I had already worked with the code). 

## `grades_pdf`
### Extracting information from a pdf and processing the data 
#### 2019 - Personal project

- `Python` : tabula, pandas, numpy

The test results for one of the most important universities in Rio de Janeiro had been released, but the official ranking would be released only 1 month later, due to possible revisions of the tests.
However, a friend of mine who had taken the exam was very anxious to find out what position she had been in for medical school.

With that, I made a notebook that extracts a dataframe from the files containing the name of all candidates and their respective grades, through the __tabula__ library.

Using __pandas__, I filtered out all the candidates who had only taken the biology and chemistry exams (specific for medical candidates), and calculated the grade.

With that, I was able to confirm that my friend was among the top 20 candidates; which was confirmed with its approval in the following month. 


## `heat_equation`
### Solving the heat equation using Fourier modes 
#### 2020 - 2A ENSTA 

- `Python` : Matplotlib, SciPy

Simple script made for a TD of the MF202 course. The main objective is to numerically solve the heat equation. Finite differences were used in the spectral space, through Fast Fourier Transforms.

## `jet_signal`
### Signal analysis of the velocity of a turbulent jet 
#### 2021 - M2 Fluid Mechanics at X.

- `Python` : SciPy, numpy, Matplotlib

Use of basic SciPy tools for signal analysis in the case of turbulent flow. Results on frequency spectral power density are obtained, in addition to signal autocorrelation analyses.
The results are compared with the existing turbulence theory, and the decay factor for energy is verified as in accordance with the studied power laws.

## `fourier_art`
### Drawings made with Fourier series harmonic circles 
#### 2020 - Personal project

- `Python` : SciPy, numpy, Matplotlib

Generation of an animation that reproduces the outline drawing of a selected image, using Fourier analysis. 

At first, the image is treated in order to more easily find a parametric function that continuously traverses the contour of the selected image.
Next, we find the representation with Fourier series for the parametric function in x and y. A threshold number of harmonics is selected as an approximation of the parametric equation. Finally, using animation functions, a video is obtained as seen in "julia.mp4", similar to the art that a friend had made for me ("art.jpg")


## `compressible_solver`
### Solving compressible problems numerically with the Saint-Venant equations 
#### 2021 - 2A ENSTA

- `Matlab`
- `Python` :  numpy, Matplotlib

Project that served as the final exam for the MF206 course (Introduction to Computation Fluid Mechanics). In it, different methods were used and compared to solve the Saint-Venant equations, such as the Roe scheme and the Rusanov scheme.
 
The code was used to solve the classical problem of dam failure (Riemann problem), in addition to the problem of an oscillating lake.

Finally, the code was implemented using Matlab and post-processing was performed in Python. 

## `modex`
### Post-processing of data from shock tube experiments
#### 2021 - 2A ENSTA

- `Python` : pandas, numpy, Matplotlib, SciPy

Quick use of Python for post-processing of shock wave data, performed in ENSTA's 2A MODEX course.

Equations for shock waves were computed and solved with SciPy. Finally, the experimental results were compared with the theory. 

## Other projects

Many other projects were put aside for confidentiality reasons (as in the case of the internship I did at Wikki Brasil), or because of my carelessness in having lost them. 





