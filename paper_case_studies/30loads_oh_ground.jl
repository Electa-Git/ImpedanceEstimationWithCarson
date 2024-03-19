import DataDrivenImpedanceEstimationWithCarson as _IMP
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
import Random as _RAN
import Ipopt

include("utils.jl")

##### other things to customize (conf file???)
# 1) solver settings
# 2) bounds on the conductor distances?

ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 500., "max_iter" => 3000)
profiles = CSV.read(_IMP.DATA_DIR*"/nrel_profiles.csv", _DF.DataFrame, ntasks = 1)
pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "print_level"=>0 )

result_path = raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_oh_shunt_most_restricted/"

for (t_start, t_end) in zip([50, 111, 161], [110, 160, 220])
        run_impedance_estimation_oh_ground_30_load_case(result_path, ie_solver, pf_solver, profiles, t_start, t_end, add_meas_noise = true, scenario_id =
scenario_id, length_bounds_percval=0.3, power_mult=1., exploit_horizontality = true, exploit_equal_crossection = true)
end

for (t_start, t_end) in zip([50, 111, 161], [110, 160, 220])
        run_impedance_estimation_oh_ground_30_load_case(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_oh_shunt_cross_only/", ie_solver, pf_solver, profiles, t_start, t_end, add_meas_noise = true, scenario_id =
scenario_id, length_bounds_percval=0.3, power_mult=1., exploit_equal_crossection = true)
end

for (t_start, t_end) in zip([50, 111, 161], [110, 160, 220])
        run_impedance_estimation_oh_ground_30_load_case(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_oh_shunt_horiz_only/", ie_solver, pf_solver, profiles, t_start, t_end, add_meas_noise = true, scenario_id =
scenario_id, length_bounds_percval=0.3, power_mult=1., exploit_horizontality = true)
end

for (t_start, t_end) in zip([50, 111, 161], [110, 160, 220])
    run_impedance_estimation_oh_ground_30_load_case(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_oh_shunt_no_restriction/", ie_solver, pf_solver, profiles, t_start, t_end, add_meas_noise = true, scenario_id =
scenario_id, length_bounds_percval=0.3, power_mult=1.)
end

function run_impedance_estimation_oh_ground_30_load_case(result_path::String, ie_solver, pf_solver, profiles::_DF.DataFrame, t_start::Int, t_end::Int; scenario_id::Int = 1, add_meas_noise::Bool=true, power_mult::Float64=1., use_length_bounds::Bool=true, length_bounds_percval::Float64=0.10, shunt_resistive::Bool=true, exploit_equal_crossection::Bool=false, exploit_horizontality::Bool=false)    

    data, eng, z_pu = prepare_math_eng_data(profiles ;feeder_name = "30load-feeder", oh_or_ug = "oh")

    build_linecode_for_ug_oh_ground_30l!(data, eng) # assigns the set of linecodes we elected for this case and builds R,X

    ######## SET MEASUREMENT SPECS
        
    max_volt_error = 2.3 # in Volts
    max_power_error = 0.1 # kW
    original_sourcebus_id = collect(keys(data["settings"]["vbases_default"]))[1]

    σ_v = 1/3*max_volt_error/(data["settings"]["vbases_default"][original_sourcebus_id]*data["settings"]["voltage_scale_factor"])
    σ_d = 1/3*max_power_error/(data["settings"]["sbase_default"])
    σ_g = 1/3*max_power_error/(data["settings"]["sbase_default"])

    ############### CREATE MULTINETWORK DATA WITH MEASUREMENT TIMESERIES ###############
    # it runs a powerflow for each time step first, so it takes some time...

    loads_with_shunts = [l for (l,load) in data["load"]][1:25]
    buses_with_shunts = [data["load"][l]["load_bus"] for l in loads_with_shunts]
    gs = _RAN.rand([10., 20., 30., 40., 50., 60., 70., 80.]./z_pu, 25)
    bs = shunt_resistive ? gs.*0. : gs./(3 * z_pu)
    mn_data, real_volts = _IMP.build_multinetwork_dsse_data_with_shunts(data, profiles, pf_solver, σ_v, σ_d, σ_g; loads_with_shunts=loads_with_shunts, gs = gs, bs= bs, t_start=t_start, t_end=t_end, add_noise=add_meas_noise, seed = scenario_id, power_mult = power_mult)

    material_resist_dict = Dict(
        "pluto" => 46.63530931436062,
        "hydrogen" => 24.047320966903072,
        "abc2x16_lv_oh_2w_bundled" => 44.72977629032573,
        "tw2x16_lv_oh_2w_bundled" => 36.23111879516384
    )

    mn_data["nw"]["1"]["linecode_map"] = Dict{Int, Any}() 
    for (id, code) in enumerate(keys(eng["linecode"]))
        if code ∈ ["tw2x16_lv_oh_2w_bundled", "abc2x16_lv_oh_2w_bundled"]
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 2,
                "r_material" => fill(material_resist_dict[code], 2)
            )
        elseif code ∈ ["pluto", "hydrogen"]
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 4,
                "r_material" => fill(material_resist_dict[code], 4)
            )
        end
    end

    make_all_branches_untrustworthy!(mn_data, eng)

    for (b,bus) in mn_data["nw"]["1"]["bus"]
        if parse(Int, b) ∉ buses_with_shunts
            bus["imp_grounded"] = fill(false, length(bus["terminals"]))
        else
            bus["imp_grounded"] = [false, true]
        end
    end  

    # materials and other carsons inputs
    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["shunt_resistive"] = shunt_resistive
    mn_data["nw"]["1"]["settings"]["exploit_equal_crossection"] = exploit_equal_crossection
    mn_data["nw"]["1"]["settings"]["exploit_horizontality"] = exploit_horizontality
    mn_data["nw"]["1"]["settings"]["oh_or_ug"] = "oh"
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

    case = "case30loads_oh_ground_"

    _IMP.drop_results(case, result_path, "", [], sol, mn_data, t_start, t_end, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, exploit_equal_crossection, exploit_squaredness, exploit_horizontality)
    _IMP.drop_shunt_results(case, result_path, "", sol, mn_data, scenario_id, t_start, t_end, shunt_resistive)
end