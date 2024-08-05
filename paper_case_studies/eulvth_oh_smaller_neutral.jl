import DataDrivenImpedanceEstimationWithCarson as _IMP
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
using MKL
import Ipopt
import HSL_jll

include("utils.jl")

ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 3600., "max_iter" => 8000,  "hsllib"=> HSL_jll.libhsl_path, "linear_solver" => "ma27")
profiles = CSV.read(_IMP.DATA_DIR*"/profiles.csv", _DF.DataFrame, ntasks = 1)
pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 200., "print_level" => 0, "hsllib"=> HSL_jll.libhsl_path, "linear_solver" => "ma27" )

timestep_set = find_most_loaded_timesteps(profiles, 100) # UNTIL THERE IS NO EDD_MEAS_NOISE!
for power_mult in [2.]
    run_impedance_estimation_oh_smaller_neutral(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\results_ma27/", ie_solver, pf_solver, profiles, vcat(timestep_set[1:70], timestep_set[72:end]), add_meas_noise = true, length_bounds_percval=0.3, power_mult=power_mult)
end

function run_impedance_estimation_oh_smaller_neutral(result_path::String, ie_solver, pf_solver, profiles::_DF.DataFrame, timestep_set; scenario_id::Int = 1, add_meas_noise::Bool=true, power_mult::Float64=1., use_length_bounds::Bool=true, length_bounds_percval::Float64=0.10 , shunt_resistive::Bool=true, exploit_equal_crossection::Bool=false, exploit_squaredness::Bool=false, exploit_horizontality::Bool=false)    

    data, eng, z_pu = prepare_math_eng_data(profiles, feeder_name = "eulvtf", oh_or_ug = "oh")

    data, eng = build_linecode_for_oh_ground_eulvtf(data, eng, z_pu) # assigns the set of linecodes we elected for this case and builds R,X

    loads_with_shunts = [l for (l,load) in data["load"]][1:30]
    buses_with_shunts = [data["load"][l]["load_bus"] for l in loads_with_shunts]
    gs = _RAN.rand(z_pu./[10., 20., 40., 50., 100., 120., 150., 200., 250., 300.], 30)
    bs = shunt_resistive ? gs.*0. : gs./(3 * 1/z_pu)

    # for (b, branch) in data["branch"]
    #     if size(branch["br_r"]) == 2
    #         branch["br_r"][2,2].*= 3.
    #     end
    # end

    mn_data, real_volts, real_vas = _IMP.build_multinetwork_dsse_data_with_shunts(data, profiles, pf_solver; timestep_set = timestep_set, loads_with_shunts=loads_with_shunts, gs = gs, bs= bs, add_noise=add_meas_noise, seed = scenario_id, power_mult = power_mult)

    add_material_properties_for_oh_ground_eulvtf!(mn_data, eng, buses_with_shunts)

    make_all_branches_untrustworthy!(mn_data, eng)

    if use_length_bounds add_length_bounds!(mn_data, length_bounds_percval) end

    ### -- ONLY DIFFERENCE WITH THE CASE THE DOES HAVE THE NEUTRAL SHUNT VARIABLE -- ###
    case1 = "eulvtf_oh_"
    case2 = "eulvtf_oh_nogroundvar_"

    mn_data_ng = deepcopy(mn_data)

    for (n,nw) in mn_data_ng["nw"]
        empty!(nw["shunt"])
        for (b,bus) in nw["bus"]
            bus["imp_grounded"] = fill(false, length(bus["terminals"]))
        end
    end
 
    ####################################################################################

    ##########################################################
    #### CASE #1: constrain both layout and cross-section ####
    ##########################################################

    mn_data["nw"]["1"]["settings"]["shunt_resistive"] = shunt_resistive
    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["exploit_horizontality"] = true
    mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"] = true
    mn_data["nw"]["1"]["settings"]["exploit_squaredness"] = false
    mn_data["nw"]["1"]["settings"]["oh_or_ug"] = "oh" 
    mn_data["nw"]["1"]["settings"]["rescaler"] = 100.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = 1.
    mn_data["nw"]["1"]["temperature"] = Dict()
    mn_data["nw"]["1"]["rho"] = Dict()
    mn_data["nw"]["1"]["alpha"] = Dict()

    mn_data_ng["nw"]["1"]["settings"] = deepcopy(mn_data["nw"]["1"]["settings"])
    mn_data_ng["nw"]["1"]["temperature"] = Dict()
    mn_data_ng["nw"]["1"]["rho"] = Dict()
    mn_data_ng["nw"]["1"]["alpha"] = Dict()

    sol_ng = _IMP.solve_imp_est_carson(mn_data_ng, ie_solver)
    sol_ng = _IMP.build_rx_sol_dict(mn_data_ng, sol_ng) # completes solution information getting together things that are not reported by default
    imp_est_ng  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data_ng, sol_ng, false)
    imp_true_ng = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data_ng, true)

    est_volts_ng = _IMP.build_estimated_volts_dataframe(sol_ng, mn_data_ng, scenario_id)
    est_vas_ng = _IMP.build_estimated_vas_dataframe(sol_ng, mn_data_ng, scenario_id)

    _IMP.drop_results(case2, result_path*"/eulvtf_oh_smaller_neutral_noground_most_restricted/", "_noshunt_power_mult_$(power_mult)_", [], sol_ng, mn_data_ng, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est_ng, imp_true_ng, real_volts, real_vas, est_volts_ng, est_vas_ng, mn_data_ng["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data_ng["nw"]["1"]["settings"]["exploit_squaredness"], mn_data_ng["nw"]["1"]["settings"]["exploit_horizontality"])
    
    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported by default
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)
    est_vas = _IMP.build_estimated_vas_dataframe(sol, mn_data, scenario_id)

    _IMP.drop_results(case1, result_path*"/eulvtf_oh_smaller_neutral_most_restricted/", "_power_mult_$(power_mult)_", [], sol, mn_data, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data["nw"]["1"]["settings"]["exploit_squaredness"], mn_data["nw"]["1"]["settings"]["exploit_horizontality"])
    _IMP.drop_shunt_results(case1, result_path*"/eulvtf_oh_smaller_neutral_most_restricted/", "_power_mult_$(power_mult)_", sol, mn_data, scenario_id, shunt_resistive)

    ##########################################################
    #### CASE #2: constrain only layout ####
    ##########################################################

    mn_data["nw"]["1"]["settings"]["shunt_resistive"] = shunt_resistive
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

    mn_data_ng["nw"]["1"]["settings"] = deepcopy(mn_data["nw"]["1"]["settings"])
    mn_data_ng["nw"]["1"]["temperature"] = Dict()
    mn_data_ng["nw"]["1"]["rho"] = Dict()
    mn_data_ng["nw"]["1"]["alpha"] = Dict()

    sol_ng = _IMP.solve_imp_est_carson(mn_data_ng, ie_solver)
    sol_ng = _IMP.build_rx_sol_dict(mn_data_ng, sol_ng) # completes solution information getting together things that are not reported by default
    imp_est_ng  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data_ng, sol_ng, false)
    imp_true_ng = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data_ng, true)

    est_volts_ng = _IMP.build_estimated_volts_dataframe(sol_ng, mn_data_ng, scenario_id)
    est_vas_ng = _IMP.build_estimated_vas_dataframe(sol_ng, mn_data_ng, scenario_id)

    _IMP.drop_results(case2, result_path*"/eulvtf_oh_smaller_neutral_noground_horizontal_only/", "_noshunt_power_mult_$(power_mult)_", [], sol_ng, mn_data_ng, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est_ng, imp_true_ng, real_volts, real_vas, est_volts_ng, est_vas_ng, mn_data_ng["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data_ng["nw"]["1"]["settings"]["exploit_squaredness"], mn_data_ng["nw"]["1"]["settings"]["exploit_horizontality"])
    
    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported by default
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)
    est_vas = _IMP.build_estimated_vas_dataframe(sol, mn_data, scenario_id)

    _IMP.drop_results(case1, result_path*"/eulvtf_oh_smaller_neutral_horizontal_only/", "_power_mult_$(power_mult)_", [], sol, mn_data, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data["nw"]["1"]["settings"]["exploit_squaredness"], mn_data["nw"]["1"]["settings"]["exploit_horizontality"])
    _IMP.drop_shunt_results(case1, result_path*"/eulvtf_oh_smaller_neutral_horizontal_only/", "_power_mult_$(power_mult)_", sol, mn_data, scenario_id, shunt_resistive)

    ##########################################################
    #### CASE #3: constrain only cross-section ####
    ##########################################################

    mn_data["nw"]["1"]["settings"]["shunt_resistive"] = shunt_resistive
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

    mn_data_ng["nw"]["1"]["settings"] = deepcopy(mn_data["nw"]["1"]["settings"])
    mn_data_ng["nw"]["1"]["temperature"] = Dict()
    mn_data_ng["nw"]["1"]["rho"] = Dict()
    mn_data_ng["nw"]["1"]["alpha"] = Dict()

    sol_ng = _IMP.solve_imp_est_carson(mn_data_ng, ie_solver)
    sol_ng = _IMP.build_rx_sol_dict(mn_data_ng, sol_ng) # completes solution information getting together things that are not reported by default
    imp_est_ng  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data_ng, sol_ng, false)
    imp_true_ng = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data_ng, true)

    est_volts_ng = _IMP.build_estimated_volts_dataframe(sol_ng, mn_data_ng, scenario_id)
    est_vas_ng = _IMP.build_estimated_vas_dataframe(sol_ng, mn_data_ng, scenario_id)

    _IMP.drop_results(case2, result_path*"/eulvtf_oh_smaller_neutral_noground_cross_only/", "_noshunt_power_mult_$(power_mult)_", [], sol_ng, mn_data_ng, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est_ng, imp_true_ng, real_volts, real_vas, est_volts_ng, est_vas_ng, mn_data_ng["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data_ng["nw"]["1"]["settings"]["exploit_squaredness"], mn_data_ng["nw"]["1"]["settings"]["exploit_horizontality"])
    
    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported by default
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)
    est_vas = _IMP.build_estimated_vas_dataframe(sol, mn_data, scenario_id)

    _IMP.drop_results(case1, result_path*"/eulvtf_oh_smaller_neutral_cross_only/", "_power_mult_$(power_mult)_", [], sol, mn_data, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data["nw"]["1"]["settings"]["exploit_squaredness"], mn_data["nw"]["1"]["settings"]["exploit_horizontality"])
    _IMP.drop_shunt_results(case1, result_path*"/eulvtf_oh_smaller_neutral_cross_only/", "_power_mult_$(power_mult)_", sol, mn_data, scenario_id, shunt_resistive)

    ##########################################################
    #### CASE #4: no constraints / restrictions ####
    ##########################################################

    mn_data["nw"]["1"]["settings"]["shunt_resistive"] = shunt_resistive
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

    mn_data_ng["nw"]["1"]["settings"] = deepcopy(mn_data["nw"]["1"]["settings"])
    mn_data_ng["nw"]["1"]["temperature"] = Dict()
    mn_data_ng["nw"]["1"]["rho"] = Dict()
    mn_data_ng["nw"]["1"]["alpha"] = Dict()

    sol_ng = _IMP.solve_imp_est_carson(mn_data_ng, ie_solver)
    sol_ng = _IMP.build_rx_sol_dict(mn_data_ng, sol_ng) # completes solution information getting together things that are not reported by default
    imp_est_ng  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data_ng, sol_ng, false)
    imp_true_ng = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data_ng, true)

    est_volts_ng = _IMP.build_estimated_volts_dataframe(sol_ng, mn_data_ng, scenario_id)
    est_vas_ng = _IMP.build_estimated_vas_dataframe(sol_ng, mn_data_ng, scenario_id)

    _IMP.drop_results(case2, result_path*"/eulvtf_oh_smaller_neutral_noground_no_restriction/", "_noshunt_power_mult_$(power_mult)_", [], sol_ng, mn_data_ng, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est_ng, imp_true_ng, real_volts, real_vas, est_volts_ng, est_vas_ng, mn_data_ng["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data_ng["nw"]["1"]["settings"]["exploit_squaredness"], mn_data_ng["nw"]["1"]["settings"]["exploit_horizontality"])
    
    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported by default
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)
    est_vas = _IMP.build_estimated_vas_dataframe(sol, mn_data, scenario_id)

    _IMP.drop_results(case1, result_path*"/eulvtf_oh_smaller_neutral_no_restriction/", "_power_mult_$(power_mult)_", [], sol, mn_data, timestep_set, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"], mn_data["nw"]["1"]["settings"]["exploit_squaredness"], mn_data["nw"]["1"]["settings"]["exploit_horizontality"])
    _IMP.drop_shunt_results(case1, result_path*"/eulvtf_oh_smaller_neutral_no_restriction/", "_power_mult_$(power_mult)_", sol, mn_data, scenario_id, shunt_resistive)

end