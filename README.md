# ImpedanceEstimationWithCarson

## Summary
This repository hosts a (unregistered) Julia package with code and data used to perform the estimation of up-to-four-wire impedance matrices for electric distribution networks, starting from measurement data (from smart meters, specifically).

From a mathematical standpoint, this is achieved through a joint (power system) state and impedance matrix estimation process. Such process maximizes the likelihood of time-variant electrical states (complex voltages/powers/currents) over a time series, by "fitting" time-invariant impedance matrix entries. Such entries are cast as functions of line and cable construction properties, i.e., distance between conductors, conductor thickness, etc., through Carson's equations.

The joint state and impedance estimation is implemented and solved as a nonlinear, non-convex mathematical optimization problem, using Julia's JuMP and the Ipopt wrapper.

Keywords: impedance matrix estimation, parameter estimation, power systems state estimation, unbalanced distribution system.

## Data and licenses

The source code is available under a BSD-3 license (see [LICENSE](LICENSE) file)

The impedance and profile data used in this work can be found in the `data` folder.

- Yearly fifteen minute resolution active power profiles are a subset of the Open Energy Data Initiative’s (OEDI) dataset [“End-Use Load Profiles for the U.S. Building Stock”](https://data.openei.org/submissions/4520), and released with a CC BY 4.0 license, like the original OEDI active power profiles. See `data/Readme_profiles.md` to check how the subset has been chosen, and how reactive power profiles are added to the original data.
- The linecode data for the impedances, curated by GridQube and found in `data/linecode_library`, are released under a CC-BY 4.0 license. 
- The network data we used consist of the **four-wire** extension of Feeder 1, Network 1 and Feeder 5, Network 4 of the Electricity North West Limited (ENWL) [dataset](https://ieeexplore.ieee.org/iel7/59/4374138/07051294.pdf) from the [Low Voltage Network Solutions project](https://www.enwl.co.uk/go-net-zero/innovation/smaller-projects/low-carbon-networks-fund/low-voltage-network-solutions/). The ENWL data are **three-wire** in their original form. ENWL granted permission to CSIRO to upload the extended four-wire data set on the [CSIRO Data Access Portal](https://doi.org/10.25919/jaae-vc35), with CC-BY 4.0. The data of the two feeders used for this work is also featured in `data/network_data/30load-feeder` (i.e.: Feeder 5, Network 4) and `data/network_data/eulvtf` (i.e.: Feeder 1, Network 1), in OpenDSS format. Minor changes w.r.t. the CSIRO portal are reported at the top of the `Lines_ ... .txt` files.

Once available, we will put here links to our relevant paper/preprints.

## Citing this work
We have submitted a paper relative to our work, currently under review:

@Misc{ID, <br />
author = {M. Vanin, F. Geth, R. Heidari, D. Van Hertem}, <br />
title = {Distribution System State and Impedance Estimation Augmented with Carson’s Equations}, <br />
howpublished = {}, <br />
year = {2024} <br />
}

## Running state and impedance estimation calculations

The general code is provided in the form of a Julia package, and can as such be installed using the [julia package manager](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

The functionalities of the package are defined in the functions included in `src/ImpedanceEstimationWithCarson.jl`. Any file which is present in `src/` but not listed here, or commented (e.g., `length_utils.jl`), is not used but kept for archiving purposes.

Notably, the variables, objective and constraints that build the optimization problem can be found in `src/core`, and are then "gathered" in the optimization problem in `src/prob/imp_est_carson.jl`. 

`src/io` contains function to clean / adjust / parse the input data from the sources above into a useable julia dictionary. All the information to build the optimization problem ultimately stems from such a dictionary. Except `src/io/solution.jl`, which contains functions to create CSV files from the optimization results, which store the solutions for successive plotting, archiving, etc.

The structure of the code is taken from [PowerModelsDistribution.jl](https://github.com/lanl-ansi/PowerModelsDistribution.jl), which is an essential dependency of this package.

The folder `validation_history` contains some old code used to verify that the mathematical model is bug-free. They are kept for author reference, but are not relevant for the package functionalities.

The folder `assets` simply contains the logos used in this Readme file (below).

The folder `paper_case_studies` contains the scripts that were used to generate the results as presented in our papers.

To reproduce the results of our first journal paper, you should run:
- `30loads_ug.jl`
- `eulvtf_ug.jl`
for the 30 load case, and the European Low Voltage Test feeder, respectively.
The power flow validations, in all papers, are run through the functions in `powerflow_validation.jl`, where you need to provide as input the folder path to where the impedance estimation results are stored (these are the csv files generated by the other relevant scripts in this folder).

To be able to run the paper script, you need to have the `ImpedanceEstimationWithCarson` package installed and set as a dependency (this is why they are in a separate folder with its own `Project.toml` file).


## Funding

This project results from a collaboration between KU Leuven/EnergyVille (Leuven and Genk, Belgium), CSIRO Energy (Newcastle, NSW, Australia) and GridQube (Springfield Central, QLD, Australia).

This work received funding from the Agency for Innovation and Entrepreneurship of the Flemish Government (VLAIO) and Flux50, through the strategic research project IMPRO-CAP (Grant N°. HBC 2022 0733).

This work was also supported by the Australian Department of Climate Change, Energy, the Environment and Water under the International Clean Innovation Researcher Networks (ICIRN) program grant number ICIRN000072. 

The collaboration involved a long overseas stay that received a grant from the Science Foundation: Flanders (FWO) (Grant N° V420224N).

<img src="./assets/readme/ku_leuven_logo.png" alt="KULeuven" height="120" width="300"/>
<img src="./assets/readme/ENERGYVILLE-ICOON.png" alt="EnergyVille" height="120" width="120"/>
<img src="./assets/readme/CSIRO-logo.png" alt="CSIRO" height="120" width="150"/>
<img src="./assets/readme/improcap_logo.png" alt="Improcap" height="120" width="450"/>
