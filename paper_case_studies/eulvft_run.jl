using Pkg
Pkg.activate("./")
import DataDrivenImpedanceEstimationWithCarson as _IMP
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
import Random as _RAN
import Ipopt

##
include("utils.jl")

include("eulvtf_oh_ground.jl")
include("eulvtf_ug.jl")

result_folder = "/Users/hei06j/Documents/repositories/remote/DataDrivenImpedanceEstimationWithCarson/results"

ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 600., "max_iter" => 3000) # maybe increase time/iterations??
profiles = CSV.read(_IMP.DATA_DIR*"/nrel_profiles.csv", _DF.DataFrame, ntasks = 1)
pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "print_level"=>0 )

## cases with series impedance only

for (t_start, t_end) in zip([50, 111, 161], [110, 160, 220])
     run_impedance_estimation_ug_noshunt_eulvtf(result_folder, ie_solver, pf_solver, profiles, t_start, t_end, add_meas_noise = true, length_bounds_percval=0.3, power_mult=1., exploit_squaredness = true, exploit_equal_crossection = true)
end

for (t_start, t_end) in zip([50, 111, 161], [110, 160, 220])
     run_impedance_estimation_ug_noshunt_eulvtf(result_folder, ie_solver, pf_solver, profiles, t_start, t_end, add_meas_noise = true, length_bounds_percval=0.3, power_mult=1., exploit_equal_crossection = true)
end
 
for (t_start, t_end) in zip([50, 111, 161], [110, 160, 220])
     run_impedance_estimation_ug_noshunt_eulvtf(result_folder, ie_solver, pf_solver, profiles, t_start, t_end, add_meas_noise = true, length_bounds_percval=0.3, power_mult=1., exploit_squaredness = true)
end
 
for (t_start, t_end) in zip([50, 111, 161], [110, 160, 220])
     run_impedance_estimation_ug_noshunt_eulvtf(result_folder, ie_solver, pf_solver, profiles, t_start, t_end, add_meas_noise = true, length_bounds_percval=0.3, power_mult=1.)
end
 
 
 