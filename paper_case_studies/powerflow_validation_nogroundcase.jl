########################################################################
#### RUNS POWER FLOW VALIDATION
#### FOR THE CASE WHERE WE DO HAVE SHUNTS BUT DO NOT ESTIMATE THEM
#### (i.e., our second paper)
########################################################################

import Ipopt
import PowerModelsDistribution as _PMD
import CSV
import DataDrivenImpedanceEstimationWithCarson as _IMP
import JSON
import StatsPlots as _SP
import HSL_jll

include("utils.jl")

general_result_path = raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\results_ma27"
pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "print_level"=>0, "hsllib"=> HSL_jll.libhsl_path, "linear_solver" => "ma27")
profiles = CSV.read(_IMP.DATA_DIR*"/profiles.csv", _DF.DataFrame, ntasks = 1)
validation_timesteps = find_most_loaded_timesteps(profiles, 400)[201:end]

########################################################################
#### RUNS POWER FLOW VALIDATION
########################################################################

for folder in ["30l_oh_noground_cross_only", "30l_oh_noground_horizontal_only", "30l_oh_noground_most_restricted", "30l_oh_noground_no_restriction"]#["eulvtf_oh_noground_cross_only", "eulvtf_noground_oh_horizontal_only", "eulvtf_oh_noground_most_restricted", "eulvtf_oh_noground_no_restriction"]#["30l_oh_noground_cross_only", "30l_oh_noground_horizontal_only", "30l_oh_noground_most_restricted", "30l_oh_noground_no_restriction"]
    if isdir(joinpath(general_result_path, folder)) && occursin("oh", folder)
        feeder_name = occursin("30", folder) ? "30load-feeder" : "eulvtf"
        result_path = joinpath(general_result_path, folder)
        if feeder_name == "30load-feeder"
            true_ground_shunts = raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\results_ma27\30l_oh_most_restricted\_30loads_oh__shunts_scenario_1__power_mult_1.0__resistive_true.json"
            for power_mult in [2.0]#, 2.0, 3.0]
                extra_id = "_power_mult_$(power_mult)"
                gen_res_file = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_general_summary_scenario_1__power_mult_$(power_mult)_", f)][1]
                general_results = CSV.read(joinpath(general_result_path, folder, gen_res_file), _DF.DataFrame, ntasks = 1)
                if string(general_results.solve_status) != "TIME_LIMIT"
                    estimated_linecode = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_linecode_results_scenario_1__power_mult_$(power_mult)_", f)][1]
                    estimated_branch_length = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_length_dict_scenario_1__power_mult_$(power_mult)_", f)][1]
                    powerflow_validation_nogroundcase(feeder_name, result_path, true_ground_shunts, joinpath(general_result_path, folder,estimated_linecode), joinpath(general_result_path, folder,estimated_branch_length), profiles, extra_id, validation_timesteps, pf_solver; power_mult=power_mult)
                end
            end
        else
            true_ground_shunts = raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\results_ma27\eulvtf_oh_no_restriction\_eulvtf_oh__shunts_scenario_1__power_mult_1.0__resistive_true.json"
            for power_mult in [2.0]
                extra_id = "_power_mult_$(power_mult)"
                estimated_linecode = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_linecode_results_scenario_1__power_mult_$(power_mult)_", f)][1]
                estimated_branch_length = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_length_dict_scenario_1__power_mult_$(power_mult)_", f)][1]
                powerflow_validation_nogroundcase(feeder_name, result_path, true_ground_shunts, joinpath(general_result_path, folder,estimated_linecode), joinpath(general_result_path, folder,estimated_branch_length), profiles, extra_id, validation_timesteps, pf_solver; power_mult=power_mult)
            end
        end
    end
end
