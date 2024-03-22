import DataDrivenImpedanceEstimationWithCarson as _IMP
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
import Ipopt

include("utils.jl")

ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 3600., "max_iter" => 10000)
profiles = CSV.read(_IMP.DATA_DIR*"/nrel_profiles.csv", _DF.DataFrame, ntasks = 1)
pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 200., "print_level"=>0 )

timestep_set = find_most_loaded_timesteps(profiles, 200)
for power_mult in [1., 2.]
    run_impedance_estimation_oh_ground_eulvtf(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results/", ie_solver, pf_solver, profiles, timestep_set, add_meas_noise = true, length_bounds_percval=0.3, power_mult=power_mult, exploit_horizontality = true, exploit_equal_crossection = true)
end

function run_impedance_estimation_oh_ground_eulvtf(result_path::String, ie_solver, pf_solver, profiles::_DF.DataFrame, timestep_set; scenario_id::Int = 1, add_meas_noise::Bool=true, power_mult::Float64=1., use_length_bounds::Bool=true, length_bounds_percval::Float64=0.10 , shunt_resistive::Bool=true, exploit_equal_crossection::Bool=false, exploit_squaredness::Bool=false, exploit_horizontality::Bool=false)    

    data, eng, z_pu = prepare_math_eng_data(profiles, feeder_name = "eulvtf", oh_or_ug = "oh")

    data, eng = build_linecode_for_oh_ground_eulvtf!(data, eng, z_pu) # assigns the set of linecodes we elected for this case and builds R,X

    mn_data, real_volts, real_vas = _IMP.build_multinetwork_dsse_data(data, profiles, pf_solver; timestep_set = timestep_set, add_noise=add_meas_noise, seed = scenario_id, power_mult = power_mult)

    loads_with_shunts = [l for (l,load) in data["load"]][1:25]
    buses_with_shunts = [data["load"][l]["load_bus"] for l in loads_with_shunts]
    gs = _RAN.rand([10., 20., 30., 40., 50., 60., 70., 80.]./z_pu, 25)
    bs = shunt_resistive ? gs.*0. : gs./(3 * z_pu)
    mn_data, real_volts, real_vas = _IMP.build_multinetwork_dsse_data_with_shunts(data, profiles, pf_solver; timestep_set = timestep_set, loads_with_shunts=loads_with_shunts, gs = gs, bs= bs, add_noise=add_meas_noise, seed = scenario_id, power_mult = power_mult)

    add_material_properties_for_oh_ground_eulvtf!(mn_data, eng, buses_with_shunts)

    make_all_branches_untrustworthy!(mn_data, eng)

    if use_length_bounds add_length_bounds!(mn_data, length_bounds_percval) end

    case = "eulvtf_oh_ground_"

    ##########################################################
    #### CASE #1: constrain both layout and cross-section ####
    ##########################################################

    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["exploit_horizontality"] = true
    mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"] = true
    mn_data["nw"]["1"]["settings"]["exploit_squaredness"] = false
    mn_data["nw"]["1"]["settings"]["oh_or_ug"] = "ug"
    mn_data["nw"]["1"]["settings"]["rescaler"] = 100.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = 1.
    mn_data["nw"]["1"]["temperature"] = Dict()
    mn_data["nw"]["1"]["rho"] = Dict()
    mn_data["nw"]["1"]["alpha"] = Dict()

    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported by default
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)
    est_vas = _IMP.build_estimated_vas_dataframe(sol, mn_data, scenario_id)

    _IMP.drop_results(case, result_path*"/eulvtf_oh_most_restricted/", "_power_mult_$(power_mult)_", [], sol, mn_data, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data["nw"]["1"]["settings"]["exploit_squaredness"], mn_data["nw"]["1"]["settings"]["exploit_horizontality"])
    _IMP.drop_shunt_results(case, result_path*"/eulvtf_oh_most_restricted/", "_power_mult_$(power_mult)_", sol, mn_data, scenario_id, shunt_resistive)

    ##########################################################
         #### CASE #2: constrain only cross-section ####
    ##########################################################

    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["exploit_horizontality"] = false 
    mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"] = true
    mn_data["nw"]["1"]["settings"]["exploit_squaredness"] = false
    mn_data["nw"]["1"]["settings"]["oh_or_ug"] = "ug"
    mn_data["nw"]["1"]["settings"]["rescaler"] = 100.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = 1.
    mn_data["nw"]["1"]["temperature"] = Dict()
    mn_data["nw"]["1"]["rho"] = Dict()
    mn_data["nw"]["1"]["alpha"] = Dict()

    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported by default
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)
    est_vas = _IMP.build_estimated_vas_dataframe(sol, mn_data, scenario_id)

    _IMP.drop_results(case, result_path*"/eulvtf_oh_cross_only/", "_power_mult_$(power_mult)_", [], sol, mn_data, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data["nw"]["1"]["settings"]["exploit_squaredness"], mn_data["nw"]["1"]["settings"]["exploit_horizontality"])
    _IMP.drop_shunt_results(case, result_path*"/eulvtf_oh_cross_only/", "_power_mult_$(power_mult)_", sol, mn_data, scenario_id, shunt_resistive)

    ##########################################################
             #### CASE #3: constrain only layout ####
    ##########################################################

    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["exploit_horizontality"] = true
    mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"] = false
    mn_data["nw"]["1"]["settings"]["exploit_squaredness"] = false
    mn_data["nw"]["1"]["settings"]["oh_or_ug"] = "ug"
    mn_data["nw"]["1"]["settings"]["rescaler"] = 100.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = 1.
    mn_data["nw"]["1"]["temperature"] = Dict()
    mn_data["nw"]["1"]["rho"] = Dict()
    mn_data["nw"]["1"]["alpha"] = Dict()

    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported by default
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)
    est_vas = _IMP.build_estimated_vas_dataframe(sol, mn_data, scenario_id)

    _IMP.drop_results(case, result_path*"/eulvtf_oh_horizontal_only/", "_power_mult_$(power_mult)_", [], sol, mn_data, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data["nw"]["1"]["settings"]["exploit_squaredness"], mn_data["nw"]["1"]["settings"]["exploit_horizontality"])
    _IMP.drop_shunt_results(case, result_path*"/eulvtf_oh_horizontal_only/", "_power_mult_$(power_mult)_", sol, mn_data, scenario_id, shunt_resistive)

    ##########################################################
                 #### CASE #4: no constraints ####
    ##########################################################

    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["exploit_horizontality"] = false 
    mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"] = false
    mn_data["nw"]["1"]["settings"]["exploit_squaredness"] = false
    mn_data["nw"]["1"]["settings"]["oh_or_ug"] = "ug"
    mn_data["nw"]["1"]["settings"]["rescaler"] = 100.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = 1.
    mn_data["nw"]["1"]["temperature"] = Dict()
    mn_data["nw"]["1"]["rho"] = Dict()
    mn_data["nw"]["1"]["alpha"] = Dict()

    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported by default
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)
    est_vas = _IMP.build_estimated_vas_dataframe(sol, mn_data, scenario_id)

    _IMP.drop_results(case, result_path*"/eulvtf_oh_no_restriction/", "_power_mult_$(power_mult)_", [], sol, mn_data, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data["nw"]["1"]["settings"]["exploit_squaredness"], mn_data["nw"]["1"]["settings"]["exploit_horizontality"])
    _IMP.drop_shunt_results(case, result_path*"/eulvtf_oh_no_restriction/", "_power_mult_$(power_mult)_", sol, mn_data, scenario_id, shunt_resistive)

end