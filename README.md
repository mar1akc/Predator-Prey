# Predator-Prey
A predator-prey model with an age-structured role reversal developed in
**L. Suarez,M. Cameron, W. Fagan, D. Levy, arXiv:<TBA>**

Software developers: L. Suarez and M. Cameron.

This package is designed to compute phase diagrams in the **(tau,g)**-plane for the predator-prey model with age-structured predator population and role reversal.

The parameter **tau** is the maturation age of the predator.
The parameter **g** is the consumption rate of the juvenile predator by the prey.


## Main codes:
**driver4AgreStructured.m**, **PlotDiagram.m**, **driver4ODEmodel.m**, **PlotDiagram4ODE.m**, **driver4DDEmodel.m**, **PlotDiagram4DDE.m**

You need to create folders **Data**, **DataODE_DDE**, and **ODEparams**.

The model parameters other than **tau** and **g** are assigned in **assign_parameters.m**.
The birth- and death-rate functions and the smoothed indicator functions are defined in **assign_parameter_functions.m**.
The model equations are defined in **time_step_helper.m**.

To compute data for the phase diagram, open the program **driver4AgeStructured.m**. Select the range of the range and the step for the tau values and run it. The run will take several hours to several days, depending on the range and the resolution. For each value of tau, all data generated by this code are saved in .mat files in the folder **Data** and can be used later for plotting phase diagrams and bifurcation diagrams.

To plot bifurcation diagrams and phase diagrams, open the program PlotDiagram.m. 
To plot the phase diagram, set **flag = 1** in line 9. The regions’ boundary data are saved in the file **ZonesPDE.mat** and can be used later to plot a combined phase diagram.

To plot a bifurcation diagram, set **flag = 2** in line 9 and assign a tau value in line 10.


To compute the data for the phase diagram for the ODE model derived from the age-structured model, run the file **driver4ODEmodel.m**. It first calls RunModel and Find_Equilib_driver.m to obtain the parameter D and the birth and death rates for the ODE model evaluated at the equilibrium of the age-structured model. These ODE parameters are saved in the folder ODEparams. Then, **EquilibriumSolutionsODE.m** is called, which computes the ODE model’s equilibria and limit cycles. The data are saved in the folder **DataODE_DDE**.

To compute the data for the phase diagram of the DDE model derived from the age-structured model, run driver4DDEmodel.m. It reads parameter values from the folder ODEparams and saves data for the phase and bifurcation diagrams in the folder **DataDDE_DDE**.

To plot the phase and bifurcation diagrams for the ODE and DDE models, use **PlotDiagram4ODE.m** and **PlotDiagram4DDE.m** respectively. The instructions are the same as for the age-structured model.

Our data for plotting phase and bifurcation diagrams for all models are available upon request.

## Folder JointPhaseDiagram
contains data for regions’ boundaries for all models, including the model with saturated birth rates, and a code **plot_zones.m** that reads these data and plots the joint bifurcation diagram.

