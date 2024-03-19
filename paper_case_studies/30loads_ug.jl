import DataDrivenImpedanceEstimationWithCarson as _IMP
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
import Ipopt

include("utils.jl")

ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 500., "max_iter" => 6000)
profiles = CSV.read(_IMP.DATA_DIR*"/nrel_profiles.csv", _DF.DataFrame, ntasks = 1)
pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 200., "print_level"=>0 )

timestep_set = find_most_loaded_timesteps(profiles, 300)

run_impedance_estimation_ug_noshunt_30_load_case(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_most_restricted/", ie_solver, pf_solver, profiles, timestep_set, add_meas_noise = true, length_bounds_percval=0.3, power_mult=3., exploit_squaredness = true, exploit_equal_crossection = true)
run_impedance_estimation_ug_noshunt_30_load_case(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_cross_only/", ie_solver, pf_solver, profiles, timestep_set, add_meas_noise = true, length_bounds_percval=0.3, power_mult=1., exploit_equal_crossection = true)
run_impedance_estimation_ug_noshunt_30_load_case(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_squared_only/", ie_solver, pf_solver, profiles, timestep_set, add_meas_noise = true, length_bounds_percval=0.3, power_mult=1., exploit_squaredness = true)
run_impedance_estimation_ug_noshunt_30_load_case(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_no_restriction/", ie_solver, pf_solver, profiles, timestep_set, add_meas_noise = true, length_bounds_percval=0.3, power_mult=1.)

function run_impedance_estimation_ug_noshunt_30_load_case(result_path::String, ie_solver, pf_solver, profiles::_DF.DataFrame, timestep_set; scenario_id::Int = 1, add_meas_noise::Bool=true, power_mult::Float64=1., use_length_bounds::Bool=true, length_bounds_percval::Float64=0.10, exploit_equal_crossection::Bool=false, exploit_squaredness::Bool=false, exploit_horizontality::Bool=false)    

    data, eng, z_pu = prepare_math_eng_data(profiles)

    data, eng = build_linecode_for_ug_noshunt_30l!(data, eng, z_pu)

    mn_data, real_volts = _IMP.build_multinetwork_dsse_data(data, profiles, pf_solver; timestep_set = timestep_set, add_noise=add_meas_noise, seed = scenario_id, power_mult = power_mult)

    add_material_properties_for_ug_noshunt_30l!(mn_data, eng)

    make_all_branches_untrustworthy!(mn_data, eng)

    # materials and other carsons inputs
    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["exploit_horizontality"] = exploit_horizontality
    mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"] = exploit_equal_crossection
    mn_data["nw"]["1"]["settings"]["exploit_squaredness"] = exploit_squaredness
    mn_data["nw"]["1"]["settings"]["oh_or_ug"] = "ug"
    mn_data["nw"]["1"]["settings"]["rescaler"] = 100.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = 1.
    mn_data["nw"]["1"]["temperature"] = Dict()
    mn_data["nw"]["1"]["rho"] = Dict()
    mn_data["nw"]["1"]["alpha"] = Dict()

    if use_length_bounds add_length_bounds!(mn_data, length_bounds_percval) end

    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported by default
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)

    case = "case30loads_series_"

    _IMP.drop_results(case, result_path, "", [], sol, mn_data, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, est_volts, exploit_equal_crossection, exploit_squaredness, exploit_horizontality)

end