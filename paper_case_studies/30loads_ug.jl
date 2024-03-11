import DataDrivenImpedanceEstimationWithCarson as _IMP
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
import Ipopt

include("utils.jl")

##### other things to customize (conf file???)
# 1) solver settings
# 2) bounds on the conductor distances?

function run_impedance_estimation_ug_noshunt_30_load_case(result_path::String, ie_solver, pf_solver, profiles::_DF.DataFrame, t_start::Int, t_end::Int; scenario_id::Int = 1, add_meas_noise::Bool=true, power_mult::Float64=1., use_length_bounds::Bool=true, length_bounds_percval::Float64=0.10)

    # ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "max_iter" => 3000)
    # profiles = CSV.read(_IMP.DATA_DIR*"/nrel_profiles.csv", _DF.DataFrame, ntasks = 1)
    # pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "print_level"=>0 )

    data, eng, z_pu = prepare_math_eng_data();

    ###################################
    ### CHANGE LINECODES OF SERVICE CABLES TO 2-WIRE (EVERYHING IS 4-WIRE IN THE BEGINNING BY CONSTRUCTION)
    ###################################

    for (b, branch) in data["branch"]
        if b ∈ ["32", "1", "2", "51", "27", "33", "28", "25", "49", "5", "43", "34", "44", "55", "37", "12", "20", "6", "7", "57", "4", "22"]
            # two-wire bits of one linecode
            r = eng["linecode"]["uglv_185al_xlpe/nyl/pvc_ug_2w_bundled"]["rs"]
            x = eng["linecode"]["uglv_185al_xlpe/nyl/pvc_ug_2w_bundled"]["xs"]
            eng["line"][branch["name"]]["linecode"] = "uglv_185al_xlpe/nyl/pvc_ug_2w_bundled"
        else
            if length(branch["f_connections"]) == 2 # two-wire bits of other linecode
                r = eng["linecode"]["ugsc_16al_xlpe/pvc_ug_2w_bundled"]["rs"]
                x = eng["linecode"]["ugsc_16al_xlpe/pvc_ug_2w_bundled"]["xs"]    
                eng["line"][branch["name"]]["linecode"] = "ugsc_16al_xlpe/pvc_ug_2w_bundled"
            else
                r = eng["linecode"]["uglv_240al_xlpe/nyl/pvc_ug_4w_bundled"]["rs"]
                x = eng["linecode"]["uglv_240al_xlpe/nyl/pvc_ug_4w_bundled"]["xs"]    
                # linecode name in `eng` for these ones is the default one        
            end
        end
        l = eng["line"][branch["name"]]["length"]
        branch["br_r"] = r.*l./z_pu
        branch["br_x"] = x.*l./z_pu
    end

    ######## SET MEASUREMENT SPECS
        
    max_volt_error = 2.3 # in Volts\
    max_power_error = 0.1 # kW
    original_sourcebus_id = collect(keys(data["settings"]["vbases_default"]))[1]

    σ_v = 1/3*max_volt_error/(data["settings"]["vbases_default"][original_sourcebus_id]*data["settings"]["voltage_scale_factor"])
    σ_d = 1/3*max_power_error/(data["settings"]["sbase_default"])
    σ_g = 1/3*max_power_error/(data["settings"]["sbase_default"])

    ############### CREATE MULTINETWORK DATA WITH MEASUREMENT TIMESERIES ###############
    # it runs a powerflow for each time step first, so it takes some time...

    # pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "print_level"=>0 )
    mn_data, real_volts = _IMP.build_multinetwork_dsse_data(data, profiles, pf_solver, σ_v, σ_d, σ_g; t_start=t_start, t_end=t_end, add_noise=add_meas_noise, seed = scenario_id, power_mult = power_mult)

    r_ac_dict = Dict(
        "uglv_185al_xlpe/nyl/pvc_ug_2w_bundled" => 0.211,
        "ugsc_16al_xlpe/pvc_ug_2w_bundled" => 2.44,
        "uglv_240al_xlpe/nyl/pvc_ug_4w_bundled" => 0.162
    )

    mn_data["nw"]["1"]["linecode_map"] = Dict{Int, Any}() 
    for (id, code) in enumerate(keys(eng["linecode"]))
        if code ∈ ["uglv_185al_xlpe/nyl/pvc_ug_2w_bundled", "ugsc_16al_xlpe/pvc_ug_2w_bundled"]
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 2,
                "r_ac" => fill(r_ac_dict[code], 2)
            )
        elseif code == "uglv_240al_xlpe/nyl/pvc_ug_4w_bundled"
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 4,
                "r_ac" => fill(r_ac_dict[code], 4)
            )
        end
    end

    make_all_branches_untrustworthy!(mn_data, eng)

    for (b,bus) in mn_data["nw"]["1"]["bus"]
        bus["imp_grounded"] = fill(false, length(bus["terminals"]))
    end  

    # materials and other carsons inputs
    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["rescaler"] = 1.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = 1.
    mn_data["nw"]["1"]["temperature"] = Dict()
    mn_data["nw"]["1"]["rho"] = Dict()
    mn_data["nw"]["1"]["alpha"] = Dict()

    # mn_data["nw"]["1"]["linecode_map"][7]["A_p_max"] = [300, 300]
    # mn_data["nw"]["1"]["linecode_map"][7]["A_p_min"] = [10, 10]

    # mn_data["nw"]["1"]["linecode_map"][9]["A_p_max"] = [300, 300]
    # mn_data["nw"]["1"]["linecode_map"][9]["A_p_min"] = [10, 10]

    # mn_data["nw"]["1"]["linecode_map"][11]["A_p_max"] = [300, 300, 300, 300]
    # mn_data["nw"]["1"]["linecode_map"][11]["A_p_min"] = [10, 10, 10, 10]

    if use_length_bounds
        for (b, branch) in mn_data["nw"]["1"]["branch"]
            branch["l_min"] = branch["orig_length"]*(1-length_bounds_percval)
            branch["l_max"] = branch["orig_length"]*(1+length_bounds_percval)
        end
    end

    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)
    sol = _IMP.build_rx_sol_dict(mn_data, sol) # completes solution information getting together things that are not reported
    imp_est = _IMP.get_cumulative_impedance_of_loads_from_sol(mn_data, sol, false)
    imp_true = _IMP.get_cumulative_impedance_of_loads_from_data(mn_data, true)

    est_volts = _IMP.build_estimated_volts_dataframe(sol, mn_data, scenario_id)

    case = "case30loads_series_"

    _IMP.drop_results(case, result_path, "", [], sol, mn_data, t_start, t_end, scenario_id, add_meas_noise, power_mult, false, false, false, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, est_volts)

end