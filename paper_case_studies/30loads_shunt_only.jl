import DataDrivenImpedanceEstimationWithCarson as _IMP
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
using MKL
import Ipopt
import Distributions as _DST
import HSL_jll

include("utils.jl")

ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 1800., "max_iter" => 8000,  "hsllib"=> HSL_jll.libhsl_path, "linear_solver" => "ma27")
profiles = CSV.read(_IMP.DATA_DIR*"/profiles.csv", _DF.DataFrame, ntasks = 1)
pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 200., "print_level" => 0, "hsllib"=> HSL_jll.libhsl_path, "linear_solver" => "ma27" )

timestep_set = find_most_loaded_timesteps(profiles, 200)

for power_mult in [1.]
    run_impedance_estimation_shunt_only_30_load_case(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\results_ma27/", ie_solver, pf_solver, profiles, timestep_set, add_meas_noise = true, power_mult=power_mult)
end

function run_impedance_estimation_shunt_only_30_load_case(result_path::String, ie_solver, pf_solver, profiles::_DF.DataFrame, timestep_set; scenario_id::Int = 1, add_meas_noise::Bool=true, power_mult::Float64=1.)    

    shunt_resistive = true

    data, eng, z_pu = prepare_math_eng_data(profiles ;feeder_name = "30load-feeder", oh_or_ug = "oh")    
    data, eng = build_linecode_for_oh_ground_30l!(data, eng, z_pu) # assigns the set of linecodes we elected for this case and builds R,X

    loads_with_shunts = [l for (l,load) in data["load"]][1:25]
    buses_with_shunts = [data["load"][l]["load_bus"] for l in loads_with_shunts]
    gs = _RAN.rand(z_pu./[10., 20., 30., 40., 50., 60., 70., 80.], 25)
    bs = shunt_resistive ? gs.*0. : gs./(3 * 1/z_pu)
    mn_data, real_volts, real_vas = _IMP.build_multinetwork_dsse_data_with_shunts(data, profiles, pf_solver; timestep_set = timestep_set, loads_with_shunts=loads_with_shunts, gs = gs, bs= bs, add_noise=add_meas_noise, seed = scenario_id, power_mult = power_mult)

    add_material_properties_for_oh_ground_30l!(mn_data, eng, buses_with_shunts)

    make_all_branches_trustworthy!(mn_data, eng) # i.e., fix the series impedances to their actual value

    case1 = "30loads_oh_"

    mn_data["nw"]["1"]["settings"]["shunt_resistive"] = shunt_resistive
    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["rescaler"] = 100.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = 1.
    mn_data["nw"]["1"]["temperature"] = Dict()
    mn_data["nw"]["1"]["rho"] = Dict()
    mn_data["nw"]["1"]["alpha"] = Dict()

    sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)

    _IMP.drop_shunt_results(case1, result_path*"/30l_shunt_only/", "_power_mult_$(power_mult)", sol, mn_data, scenario_id, shunt_resistive)

    # for (n,nw) in mn_data["nw"]
    #     for (m,meas) in nw["meas"]
    #         if meas["var"] == :vm
    #             meas["crit"] = "rwls"
    #         end
    #     end
    # end

    # sol = _IMP.solve_shunt_only_vm_in_obj_only(mn_data, ie_solver)

    # _IMP.drop_shunt_results(case1, result_path*"/30l_shunt_only/", "_power_mult_$(power_mult)_vm_wls", sol, mn_data, scenario_id, shunt_resistive)

end