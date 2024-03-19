import DataDrivenImpedanceEstimationWithCarson as _IMP
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
import Ipopt

include("utils.jl")


function run_impedance_estimation_ug_noshunt_eulvtf(result_path::String, ie_solver, pf_solver, profiles::_DF.DataFrame, t_start::Int, t_end::Int; scenario_id::Int = 1, add_meas_noise::Bool=true, power_mult::Float64=1., use_length_bounds::Bool=true, length_bounds_percval::Float64=0.10, exploit_equal_crossection::Bool=false, exploit_squaredness::Bool=false, exploit_horizontality::Bool=false)    

    data, eng, z_pu = prepare_math_eng_data(profiles, feeder_name = "eulvtf")

    build_linecode_for_ug_noshunt_eulvtf!(data, eng) # assigns the set of linecodes we elected for this case and builds R,X

    ######## SET MEASUREMENT SPECS
        
    max_volt_error = 2.3 # in Volts\
    max_power_error = 0.1 # kW
    original_sourcebus_id = collect(keys(data["settings"]["vbases_default"]))[1]

    σ_v = 1/3*max_volt_error/(data["settings"]["vbases_default"][original_sourcebus_id]*data["settings"]["voltage_scale_factor"])
    σ_d = 1/3*max_power_error/(data["settings"]["sbase_default"])
    σ_g = 1/3*max_power_error/(data["settings"]["sbase_default"])

    ############### CREATE MULTINETWORK DATA WITH MEASUREMENT TIMESERIES ###############
    # it runs a powerflow for each time step first, so it takes some time...

    mn_data, real_volts = _IMP.build_multinetwork_dsse_data(data, profiles, pf_solver, σ_v, σ_d, σ_g; t_start=t_start, t_end=t_end, add_noise=add_meas_noise, seed = scenario_id, power_mult = power_mult)

    add_material_properties_for_ug_noshunt_eulvtf!(mn_data)

    make_all_branches_untrustworthy!(mn_data, eng)

    for (b,bus) in mn_data["nw"]["1"]["bus"]
        bus["imp_grounded"] = fill(false, length(bus["terminals"]))
    end  

    # materials and other carsons inputs
    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["exploit_horizontality"] = exploit_horizontality
    mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"] = exploit_equal_crossection
    mn_data["nw"]["1"]["settings"]["exploit_squaredness"] = exploit_squaredness
    mn_data["nw"]["1"]["settings"]["oh_or_ug"] = "ug"
    mn_data["nw"]["1"]["settings"]["rescaler"] = 1.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = 1.
    mn_data["nw"]["1"]["temperature"] = Dict()
    mn_data["nw"]["1"]["rho"] = Dict()
    mn_data["nw"]["1"]["alpha"] = Dict()

    if use_length_bounds add_length_bounds!(mn_data, length_bounds_percval) end

    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported
    imp_est  = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)
    est_vas = _IMP.build_estimated_vas_dataframe(sol, mn_data, scenario_id)

    case = "eulvtf_series_"

    _IMP.drop_results(case, result_path, "", [], sol, mn_data, t_start, t_end, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, exploit_equal_crossection, exploit_squaredness, exploit_horizontality)

end