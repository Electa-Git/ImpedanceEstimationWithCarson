import StatsPlots as _SP
import Format

"""
plots a boxplot for every user. different sets of dots for different options / estimated paths
"""
function cumulative_zrx_per_user_boxplot(true_impedance_dict::Dict, est_impedance_path::String, case::String; whatt::String="Zc", ytickz::Vector=[])
    
    est_dicts = [joinpath(est_impedance_path, file) for file in readdir(est_impedance_path) if occursin("imp_est", file) && occursin(case, file)]

    nr_scenarios = unique([split(ed, "scenario_")[2][1] for ed in est_dicts])
    other_identifiers = unique([split(ed, "__")[end][1:end-5] for ed in est_dicts])

    idx = findfirst(x->occursin("scenario_$(minimum(nr_scenarios))", x), est_dicts)
    starter_est = est_dicts[idx]
    starter_label = split(starter_est, "__")[end][1:end-5]

    xticks = (1:length(true_impedance_dict), [parse(Int, k) for (k,v) in true_impedance_dict])

    p = _SP.scatter([(v["$(whatt)_true"]-JSON.parsefile(joinpath(est_impedance_path, starter_est))[k]["$(whatt)_est"])/v["$(whatt)_true"]*100 for (k,v) in true_impedance_dict],
                    ylabel = "($(whatt)ᵗʳᵘᵉ - $(whatt)ᵉˢᵗ)/$(whatt)ᵗʳᵘᵉ × 100 [%]", labels = "sc. $(minimum(nr_scenarios)), id. $(starter_label)",
                    xlabel = "User id. [-]", xticks = xticks)

    for ed in est_dicts
        if ed != starter_est
            sc = split(ed, "scenario_")[2][1]
            id = split(ed, "__")[end][1:end-5]
           _SP.scatter!([(v["$(whatt)_true"]-JSON.parsefile(joinpath(est_impedance_path, ed))[k]["$(whatt)_est"])/v["$(whatt)_true"]*100 for (k,v) in true_impedance_dict],
                        labels = "sc. $sc, id. $id")
        end
    end

    return p
end

"""
plots a dot for every user
"""
function cumulative_zrx_per_user(true_impedance_dict::Dict, est_impedance_dict::Dict, case::String; whatt::String="Zc", ytickz::Vector=[])

    xticks = (1:length(true_impedance_dict), [parse(Int, k) for (k,v) in true_impedance_dict])

    p = _SP.scatter([(v["$(whatt)_true"]-est_impedance_dict[k]["$(whatt)_est"])/v["$(whatt)_true"]*100 for (k,v) in true_impedance_dict],
                    ylabel = "($(whatt)ᵗʳᵘᵉ - $(whatt)ᵉˢᵗ)/$(whatt)ᵗʳᵘᵉ × 100 [%]", 
                    xlabel = "User id. [-]", xticks = xticks)

    return p
end

function rx_linecode_matrix_estimated()
    
    p = _SP.heatmap(zeros(4,4),
               c = _SP.cgrad([:blue,:white,:red]),
               clims=(-1, 1), legend = :none,
               xticks = ([], []), yticks = ([], []), xaxis = false, yaxis=false, 
               xlims=(-0.01, 4.01), ylims=(-0.01, 4.01),
               aspect_ratio = 1/1.5 )
    _SP.vline!([1,2,3], color=:lightgrey, linestyle = :dashdotdot)
    _SP.hline!([1,2,3], color=:lightgrey, linestyle = :dashdotdot)
    _SP.hline!([0.,4.], color=:black)
    _SP.vline!([0.,4], color=:black)

    _SP.annotate!(1, 2, _SP.text("mytext", :black, :center, 10))

    return p
end

function rx_linecode_matrix_exact(linecode_name::String; whatt::String="rs")

    eng = _PMD.parse_file(_IMP.NTW_DATA_DIR*"/eulvtf/Master_oh.dss", data_model = _PMD.ENGINEERING, transformations=[_PMD.transform_loops!,_PMD.remove_all_bounds!])
    mat = eng["linecode"][linecode][whatt]

    p = _SP.heatmap(zeros( size(mat)[1], size(mat)[1]),
               c = _SP.cgrad([:blue,:white,:red]),
               clims=(-1, 1), legend = :none,
               xticks = ([], []), yticks = ([], []), xaxis = false, yaxis=false, 
               xlims=(-0.01, 4.01), ylims=(-0.01, 4.01),
               aspect_ratio = 1/1.5 )
    if size(mat)[1] == 4
        _SP.vline!([1,2,3], color=:lightgrey, linestyle = :dashdotdot)
        _SP.hline!([1,2,3], color=:lightgrey, linestyle = :dashdotdot)
        _SP.hline!([0.,4.], color=:black)
        _SP.vline!([0.,4], color=:black)
        count_row = 1
        for r in 1:size(mat)[1]
            count_col = 1
            for c in 1:size(mat)[1]
                str =  Format.cfmt( "%.4e", mat[r,c]) 
                _SP.annotate!(count_row-0.5, 4.5-count_col, _SP.text(str, :black, :center, 10))
                count_col+=1
            end
            count_row+=1
        end
    else
        _SP.vline!([1], color=:lightgrey, linestyle = :dashdotdot)
        _SP.hline!([1], color=:lightgrey, linestyle = :dashdotdot)
        _SP.hline!([0.,2], color=:black)
        _SP.vline!([0.,2], color=:black)
        count_row = 1
        for r in 1:size(mat)[1]
            count_col = 1
            for c in 1:size(mat)[1]
                str =  Format.cfmt( "%.4e", mat[r,c]) 
                _SP.annotate!(count_row-0.5, 2.5-count_col, _SP.text(str, :black, :center, 10))
                count_col+=1
            end
            count_row+=1
        end
    end    

    return p
end


# folders = [raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_cross_only", raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_most_restricted", raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_squared_only", raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_no_restriction"]

function generate_summary(folders::Vector{String})
    super_summary = _DF.DataFrame()
    for f in folders
        content = readdir(f)
        for file in content 
            if occursin("general_summary", file)
                summary = CSV.read(joinpath(f, file), _DF.DataFrame, stringtype=String)
                if isempty(super_summary)
                    super_summary = summary
                else
                    for r in eachrow(summary)
                        push!(super_summary, r)
                    end
                end
            end
        end
    end
    return super_summary
end

function powerflow_validation(feeder_name, oh_or_ug, result_path::String, estimated_linecode::String, estimated_branch_length::String, profiles, extra_id::String, validation_timesteps, pf_solver; power_mult::Float64=2.)
        
    data, eng, z_pu = prepare_math_eng_data(profiles ;feeder_name = feeder_name, oh_or_ug = oh_or_ug)

    real_volts = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["time_step", "termination_status"]))
    est_volts  = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["time_step", "termination_status"]))

    real_va = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["time_step", "termination_status"]))
    est_va  = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["time_step", "termination_status"]))

    if feeder_name == "30load-feeder"
        if oh_or_ug == "ug"
            build_linecode_for_ug_noshunt_30l!(data,eng,z_pu)
        end
    end

    estimated_branch_length = JSON.parsefile(estimated_branch_length)
    estim_lc = CSV.read(estimated_linecode, _DF.DataFrame, ntasks = 1)
    estimated_linecode = CSV.File(estimated_linecode)
    estim_lc.r_est = eval.(Meta.parse.(estimated_linecode.r_est))
    estim_lc.x_est = eval.(Meta.parse.(estimated_linecode.x_est))
    
    linecode_dict = Dict(row["linecode_name"] => Dict(
        "xs"   => row["x_est"],
        "rs"   => row["r_est"])
        for row in eachrow(estim_lc)
    )

    estimated_data = deepcopy(data)

    for (b,branch) in estimated_data["branch"]
        linecode = linecode_dict[eng["line"][branch["name"]]["linecode"]]
        branch["br_r"] = linecode["rs"].*estimated_branch_length[b]["length_est"]./(z_pu*1000)
        branch["br_x"] = linecode["xs"].*estimated_branch_length[b]["length_est"]./(z_pu*1000)
    end

    for (ts_id, ts) in enumerate(validation_timesteps)
        
        _IMP.insert_profiles!(data, profiles, ts, power_mult=power_mult)
        _IMP.insert_profiles!(estimated_data, profiles, ts, power_mult=power_mult)

        real_pf_results = _PMD.solve_mc_opf(data, _PMD.IVRENPowerModel, pf_solver)
        est_pf_results = _PMD.solve_mc_opf(estimated_data, _PMD.IVRENPowerModel, pf_solver)
        
        # converts vr and vi to vm (phase to neutral)
        _IMP.pf_solution_to_voltage_magnitudes!(real_pf_results) 
        _IMP.pf_solution_to_voltage_magnitudes!(est_pf_results) 
        _IMP.pf_solution_to_voltage_angles!(real_pf_results)
        _IMP.pf_solution_to_voltage_angles!(est_pf_results)

        push!(real_volts, vcat([real_pf_results["solution"]["bus"]["$(load["load_bus"])"]["vm"][1] for (l,load) in data["load"]], [ts, real_pf_results["termination_status"]]))
        push!(est_volts, vcat([est_pf_results["solution"]["bus"]["$(load["load_bus"])"]["vm"][1] for (l,load) in data["load"]], [ts, est_pf_results["termination_status"]]))
        push!(real_va, vcat([real_pf_results["solution"]["bus"]["$(load["load_bus"])"]["va"][1] for (l,load) in data["load"]], [ts, real_pf_results["termination_status"]]))
        push!(est_va, vcat([est_pf_results["solution"]["bus"]["$(load["load_bus"])"]["va"][1] for (l,load) in data["load"]], [ts, est_pf_results["termination_status"]]))

    end

    CSV.write(result_path*"/pf_validation_real_vm"*extra_id*".csv", real_volts)
    CSV.write(result_path*"/pf_validation_est_vm"*extra_id*".csv" ,  est_volts)
    CSV.write(result_path*"/pf_validation_real_va"*extra_id*".csv", real_va)
    CSV.write(result_path*"/pf_validation_est_va"*extra_id*".csv" ,  est_va)
end

function powerflow_validation_boxplot(est_pf::Union{String, _DF.DataFrame}, real_pf::Union{String, _DF.DataFrame}; per_user::Bool=false, yticks=[], show_sm_info::Bool=false)
    est_pf = est_pf isa String ? CSV.read(est_pf, _DF.DataFrame, ntasks=1) : est_pf
    real_pf = real_pf isa String ? CSV.read(real_pf, _DF.DataFrame, ntasks=1) : real_pf

    diff_df = deepcopy(est_pf)
    for row in 1:size(diff_df)[1]
        diff_df[row:row,1:end-2] .= abs.(est_pf[row:row,1:end-2] .- real_pf[row:row,1:end-2])
    end

    if !per_user # we do a single boxplot, not a boxplot per each user
        all_diffs = [i for i in vec(Matrix(diff_df)) if (i isa Float64)]
        p = _SP.boxplot(all_diffs*240, ylabel="| |U|ᵗʳᵘᵉ-|U|ᵉˢᵗ | [V]", xticks = ((),()), label = "", legend=false)
        if !isempty(yticks)
            _SP.plot!(yticks=yticks)
        end        
        if show_sm_info
            _SP.hline!([240*0.005], label = "3σ |U| error", legend=true)
            _SP.hline!([240*2*0.005/3], label = "2σ |U| error", legend=true)
            _SP.hline!([240*0.005/3], label = "σ |U| error", legend=true)
        end
    end
    return p
end

# cumulative_zrx_per_user(JSON.parsefile(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_most_restricted\_case30loads_series__imp_true_scenario_1_.json"), JSON.parsefile(raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\paper_results\30l_ug_most_restricted\_case30loads_series__imp_est_scenario_1__ts_50_te_110.json"), " "; whatt="Xc")