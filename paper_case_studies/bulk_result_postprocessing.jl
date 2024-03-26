general_result_path = raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results"
validation_timesteps = find_most_loaded_timesteps(profiles, 400)[201:end]
pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "print_level"=>0 )
profiles = CSV.read(_IMP.DATA_DIR*"/nrel_profiles.csv", _DF.DataFrame, ntasks = 1)

########################################################################
#### RUNS POWER FLOW VALIDATION
########################################################################

for folder in readdir(general_result_path)
    if isdir(joinpath(general_result_path, folder)) && occursin("ug", folder)
        feeder_name = occursin("30", folder) ? "30load-feeder" : "eulvtf"
        result_path = joinpath(general_result_path, folder)
        if feeder_name == "30load-feeder"
            for power_mult in [1.0, 2.0, 3.0]
                extra_id = "_power_mult_$(power_mult)"
                gen_res_file = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_general_summary_scenario_1__power_mult_$(power_mult)_", f)][1]
                general_results = CSV.read(joinpath(general_result_path, folder, gen_res_file), _DF.DataFrame, ntasks = 1)
                if string(general_results.solve_status) != "TIME_LIMIT"
                    estimated_linecode = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_linecode_results_scenario_1__power_mult_$(power_mult)_", f)][1]
                    estimated_branch_length = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_length_dict_scenario_1__power_mult_$(power_mult)_", f)][1]
                    powerflow_validation(feeder_name, "ug", result_path, joinpath(general_result_path, folder,estimated_linecode), joinpath(general_result_path, folder,estimated_branch_length), profiles, extra_id, validation_timesteps, pf_solver; power_mult=power_mult)
                end
            end
        else
            for power_mult in [1.0, 2.0]
                extra_id = "_power_mult_$(power_mult)"
                estimated_linecode = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_linecode_results_scenario_1__power_mult_$(power_mult)_", f)][1]
                estimated_branch_length = [f for f in readdir(joinpath(general_result_path, folder)) if occursin("_length_dict_scenario_1__power_mult_$(power_mult)_", f)][1]
                powerflow_validation(feeder_name, "ug", result_path, joinpath(general_result_path, folder,estimated_linecode), joinpath(general_result_path, folder,estimated_branch_length), profiles, extra_id, validation_timesteps, pf_solver; power_mult=power_mult)
            end
        end
    end
end

########################################################################
#### PLOTS POWER FLOW VALIDATION - case by case
########################################################################

for folder in readdir(general_result_path)
    for pm in [1.0, 2.0, 3.0]
        if "pf_validation_est_vm_power_mult_$pm.csv" ∈ readdir(joinpath(general_result_path, folder))
            real_pf = CSV.read(joinpath(general_result_path, folder, "pf_validation_real_vm_power_mult_$pm.csv"),_DF.DataFrame, ntasks = 1)
            est_pf = CSV.read(joinpath(general_result_path, folder, "pf_validation_est_vm_power_mult_$pm.csv"),_DF.DataFrame, ntasks = 1)
            p = powerflow_validation_boxplot(est_pf, real_pf)
            _SP.plot!(title = "$folder - power mult $pm")
            _SP.savefig((joinpath(general_result_path, folder))*"/pf_validation_power_mult_$pm.png")
        end
    end
end

########################################################################
#### PLOTS POWER FLOW VALIDATION - crossconstraint
########################################################################

# TODO

########################################################################
#### PLOTS CUMULATIVE IMPEDANCES PER USER
########################################################################

for folder in readdir(general_result_path)
    if isdir(joinpath(general_result_path, folder))
        if any([occursin("imp_est_scenario_1__power_mult_", i) for i in readdir(joinpath(general_result_path, folder))])
            case = [i for i in readdir(joinpath(general_result_path, folder)) if occursin("imp_true_scenario_1__power_mult_", i)][1][1:7]
            true_impedance_dict = JSON.parsefile(joinpath(general_result_path, folder, [i for i in readdir(joinpath(general_result_path, folder)) if occursin("imp_true_scenario_1__power_mult_", i)][1]))
            est_impedance_path = joinpath(general_result_path, folder)
            for whatt ∈ ["Xc", "Rc", "Zc"]
                p =  cumulative_zrx_per_user(true_impedance_dict, est_impedance_path, case; whatt=whatt)
                _SP.savefig((joinpath(general_result_path, folder))*"/cumulative_$(whatt)_diff_mult.png")
            end
        end
    end
end

########################################################################
#### PLOTS CUMULATIVE IMPEDANCES BOXPLOT - across power multipliers
########################################################################

for folder in readdir(general_result_path)
    if isdir(joinpath(general_result_path, folder))
        if any([occursin("imp_est_scenario_1__power_mult_", i) for i in readdir(joinpath(general_result_path, folder))])
            case = [i for i in readdir(joinpath(general_result_path, folder)) if occursin("imp_true_scenario_1__power_mult_", i)][1][1:7]
            true_impedance_dict = JSON.parsefile(joinpath(general_result_path, folder, [i for i in readdir(joinpath(general_result_path, folder)) if occursin("imp_true_scenario_1__power_mult_", i)][1]))
            est_impedance_path = joinpath(general_result_path, folder)
            for whatt ∈ ["Xc", "Rc", "Zc"]
                p =  cumulative_zrx_boxplot(true_impedance_dict, est_impedance_path, case; whatt=whatt)
                _SP.savefig((joinpath(general_result_path, folder))*"/cumulative_$(whatt)_diff_mult_boxplot.png")
            end
        end
    end
end

########################################################################
#### PLOTS CUMULATIVE IMPEDANCES BOXPLOT - across constraint sets
########################################################################

case = "30l_ug"
for whatt in ["Zc", "Rc", "Xc"]
    for power_mult in [1.0, 2.0, 3.0]
        for folder in readdir(general_result_path)
            if isdir(joinpath(general_result_path, folder))
                if occursin(case, folder) && any([occursin("imp_est_scenario_1__power_mult_$(power_mult)", i) for i in readdir(joinpath(general_result_path, folder))])
                    true_impedance_dict = JSON.parsefile(joinpath(general_result_path, folder, [i for i in readdir(joinpath(general_result_path, folder)) if occursin("imp_true_scenario_1__power_mult_", i)][1]))
                    # est_impedance_path = joinpath(general_result_path, folder)
                    for whatt ∈ ["Xc", "Rc", "Zc"]
                        p =  cumulative_zrx_boxplot_crossconstraint(general_result_path, case, power_mult; whatt=whatt)
                    end
                end
            end
        end
    _SP.savefig(general_result_path*"/cumulative_$(whatt)_diff_boxplot_constr_pm_$power_mult.png")
    end
end

########################################################################
#### PLOT LINECODE DIFFERENCES
########################################################################

case = "30l_ug"
feeder_id = "30load-feeder"
power_mult = 3.0
rx_linecode_matrix_estimated_crossconstraint(general_result_path, case, feeder_id, power_mult)

rx_linecode_matrix_estimated_crossconstraint_perc(general_result_path, case, feeder_id, power_mult)