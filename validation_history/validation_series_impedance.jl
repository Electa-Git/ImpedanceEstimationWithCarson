import DataDrivenImpedanceEstimationWithCarson as _IMP
import PowerModelsDistribution as _PMD
import Ipopt

include("validation_utils.jl")

ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "max_iter" => 1500)

#################################################################################################
### Validation steps: ###
# BASIC CASE: P,Q,|U| from load and generator buses are used as measurements, no noise is added
# ALL-VOLTS CASE: noiseless voltage measurement is added at node `1` (zero-injection bus)
# BOUNDED CASE: near-exact bounds are added on A_p of the single-phase conductors
# OTHER CASES could play on assumptions on the branch lengths (GIS data not accurate but also not kilometres off!)

# in all cases, bus shunts are zero and are not defined as variables (only carson series impedance is sought for)
#################################################################################################

# BASIC CASE:

mn_data = set_up_validation_series(;power_mult = 15., t_start = 1, t_end = 10)
sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
sol = build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported

# note that the cumulative impedance metric is imperfect because it does not account for mutual impedance (self-only)
imp_true = get_cumulative_impedance_of_loads_from_data(mn_data, true)
imp_est = get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)


# ALL-VOLTS CASE:

add_zib_voltage!(mn_data)

sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
sol = build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported

imp_est = get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)

# BOUNDED CASE (the real value of `A_p` for linecode `2` is 30 😄)

### NB!!! sometimes this case takes the largest amount of iterations, yet it has a lower minimum then the previous cases.
###       i.e., previously stuck in a local optimum??
###       can probably pick better starting values!! 

mn_data["nw"]["1"]["linecode_map"][2]["A_p_max"] = [33, 33]
mn_data["nw"]["1"]["linecode_map"][2]["A_p_min"] = [26, 26]

sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
sol = build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported
imp_est = get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)

#################################################################################################
### DEGREE of SUCCESS of the IMPEDANCE ESTIMATION: ###
#
# EQUIVALENT MODEL: cumulative impedance path user-transformer is accurate but we cannot retrieve any linecode information. Also, voltages are close to ground-truth PF result. → this model is ok to use for (main frequency) (O)PF & co.
# ACCURATE LINECODE IMPEDANCE: all linecode impedances are accurate, but the geometry is unrealistic (multiple solutions due to multiplication of different variables → degenerate impedance model) → this model is ok to use for (O)PF & co., and maybe with some human in the loop, can reverse-engineer the linecodes
# EXACT LINECODE: correct linecode geometry is retrieved for all unique linecode variables → the best we could wish for, closest to ground truth 
# MIXED: combinations of the above (e.g., some linecodes are exactly reconstructed, some just accurate...)
#################################################################################################

#### HOW TO CHECK THE DEGREE of SUCCESS

# 1- the cumulative impedance retrieval ("EQUIVALENT MODEL") is shown above! (`imp_est` vs `imp_true`)
# 2- the "r's" and "x's" of the linecode should be the same as those in `small_validation_case/LineCode.txt`
# 3- the linecode_map geometry in the solution matches the original one (i.e., `A_p` and distances are correct. Coordinates can be all over the place, it's ok.) 
sol["solution"]["nw"]["1"]["linecode_map"]["2"]


################ FAIRLY CONFIDENT ABOUT THE VALIDATION BECAUSE:
# 1) regardless of whether branches are trustworthy or not, residuals are 0 if no measurement noise is added
#    - i.e., result is the same as 100% PMD-generated power flow
# 2) good quality solutions are found regardless, but adding bounds, etc., leads to the right solution for all linecodes (again, in noiseless scenario), hurray!
