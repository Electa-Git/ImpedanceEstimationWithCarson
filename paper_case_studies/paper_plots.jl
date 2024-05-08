import StatsPlots as _SP
using LaTeXStrings

general_result_path = raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\results_ma27"

########################################################
#####                                              #####
#####         POWER FLOW VALIDATION PLOTS          #####
#####                                              #####
########################################################

function plot_30l_ug_pf_validation(general_result_path)
    case = "30l_ug"
    for folder in [i for i in readdir(general_result_path) ] 
        if isdir(joinpath(general_result_path, folder))
            if "pf_validation_est_vm_power_mult_1.0.csv" ∈ readdir(joinpath(general_result_path, folder))
                powerflow_validation_paperplot(general_result_path, case, ytickz = [0., 0.5, 1., 1.5])
                _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper.png")
                _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper.pdf")
            end
        end
    end
end

function plot_30l_oh_pf_validation(general_result_path)
    case = "30l_oh"
    for folder in [i for i in readdir(general_result_path) ] 
        if isdir(joinpath(general_result_path, folder))
            if "pf_validation_est_vm_power_mult_1.0.csv" ∈ readdir(joinpath(general_result_path, folder))
                powerflow_validation_paperplot(general_result_path, case, ytickz = [0., 0.5, 1., 1.5, 2., 2.5])
                _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper.png")
                _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper.pdf")
            end
        end
    end
end

function plot_eu_ug_pf_validation(general_result_path)
    case = "eulvtf_ug"
    for folder in [i for i in readdir(general_result_path) ] 
        if isdir(joinpath(general_result_path, folder))
            if "pf_validation_est_vm_power_mult_1.0.csv" ∈ readdir(joinpath(general_result_path, folder))
                powerflow_validation_paperplot(general_result_path, case, ytickz = [0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.])
                _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper.png")
                _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper.pdf")
            end
        end
    end
end

function plot_eu_oh_pf_validation(general_result_path)
    case = "eulvtf_oh"
    for folder in [i for i in readdir(general_result_path) ] 
        if isdir(joinpath(general_result_path, folder))
            if "pf_validation_est_vm_power_mult_1.0.csv" ∈ readdir(joinpath(general_result_path, folder))
                powerflow_validation_paperplot(general_result_path, case, ytickz = [0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.])
                _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper.png")
                _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper.pdf")
            end
        end
    end
end

function powerflow_validation_paperplot(general_result_path::String, case::String; ytickz::Vector=[])

    folders = [f for f in readdir(general_result_path) if occursin(case, f) && !occursin(".png", f)]
    est_dfs = []
    real_dfs = []
    legend_filenames = []
    for folder in folders 
        if isdir(joinpath(general_result_path, folder))
            for file in readdir(joinpath(general_result_path, folder))
                if occursin("pf_validation_real", file) && occursin("vm_power_mult_1.0", file)
                    real_pf = CSV.read(joinpath(general_result_path, folder, "pf_validation_real_vm_power_mult_1.0.csv"),_DF.DataFrame, ntasks = 1)
                    est_pf  = CSV.read(joinpath(general_result_path, folder, "pf_validation_est_vm_power_mult_1.0.csv"),_DF.DataFrame, ntasks = 1)
                    if occursin("no_restriction", folder) push!(legend_filenames, L"\textrm{No  \, \, rest.}") end
                    if occursin("most_restricted", folder) push!(legend_filenames, L"A_p+\mathcal{G}  \, \, \textrm{rest.}") end
                    if occursin("cross", folder) push!(legend_filenames, L"A_p  \, \, \textrm{rest.}") end
                    if (occursin("horizontal", folder) || occursin("square", folder)) push!(legend_filenames, L"\mathcal{G}  \, \, \textrm{rest.}") end
                    push!(est_dfs, est_pf)
                    push!(real_dfs, real_pf)
                end
            end
        end
    end

    starter_est = est_dfs[1]
    starter_real = real_dfs[1]

    diff_df = deepcopy(starter_est)
    for row in 1:size(diff_df)[1]
        diff_df[row:row,1:end-2] .= abs.(starter_est[row:row,1:end-2] .- starter_real[row:row,1:end-2])
    end
    all_diffs = [i for i in vec(Matrix(diff_df)) if (i isa Float64)]

    p = _SP.boxplot(all_diffs*240, ylabel= L"\Delta \textrm{U}^{\textrm{mag}} \, \, [\textrm{V}]",#L"| \textrm{U}^{\textrm{mag, true}}-\textrm{U}^{\textrm{mag, est}}| \, \, [\textrm{V}]", 
                                   color = :lightgrey, legend = false, xtickfontsize=11,ytickfontsize=11, ylabelfontsize=14)

    # adds the other cases
    for i in 2:length(est_dfs)
        diff_df = deepcopy(est_dfs[i])
        for row in 1:size(diff_df)[1]
            diff_df[row:row,1:end-2] .= abs.(est_dfs[i][row:row,1:end-2] .- real_dfs[i][row:row,1:end-2])
        end
        all_diffs = [i for i in vec(Matrix(diff_df)) if (i isa Float64)]
        p = _SP.boxplot!(all_diffs*240, color = :lightgrey, xticks = ((1:4),legend_filenames))
    end
       
    if !isempty(ytickz)
        _SP.plot!(yticks=(ytickz, [L"%$i" for i in ytickz]))
    end 

    return p
end

########################################################
#####                                              #####
#####         CUMULATIVE IMPEDANCE PLOTS           #####
#####                                              #####
########################################################

function plot_30l_ug_cumulative_diffs_rxz_mult_1(general_result_path)
    case = "30l_ug"
    for folder in [i for i in readdir(general_result_path) ]
        if isdir(joinpath(general_result_path, folder))
            if occursin(case, folder) && any([occursin("imp_est_scenario_1__power_mult_1.0", i) for i in readdir(joinpath(general_result_path, folder))])
                for whatt in ["Xc", "Rc", "Zc"]
                    ylims = whatt == "Xc" ? (-0.4, 30.4) : ()
                    if whatt == "Xc" yticks = [i for i in 0:5:30] end 
                    if whatt == "Rc" yticks = [i for i in 0:100:300] end  
                    if whatt == "Zc" yticks = [i for i in 0:10:40] end  
                    p =  cumulative_zrx_boxplot_crossconstraint_paper(general_result_path, case, 1.0; whatt=whatt, ytickz = yticks, ylims = ylims)
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm1.png")
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm1.pdf")
                end
            end
        end
    end
end

function plot_30l_ug_cumulative_diffs_rxz_mult_3(general_result_path)
    case = "30l_ug"
    for folder in [i for i in readdir(general_result_path) ]
        if isdir(joinpath(general_result_path, folder))
            if occursin(case, folder) && any([occursin("imp_est_scenario_1__power_mult_3.0", i) for i in readdir(joinpath(general_result_path, folder))])
                for whatt in ["Xc", "Rc", "Zc"]
                    ylims = whatt == "Xc" ? (-0.4, 30.4) : ()  
                    yticks = whatt == "Xc" ? [i for i in 0:5:30] : []  
                    p =  cumulative_zrx_boxplot_crossconstraint_paper(general_result_path, case, 3.0; whatt=whatt, ytickz = yticks, ylims = ylims)
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm3.png")
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm3.pdf")
                end
            end
        end
    end
end

function plot_30l_oh_cumulative_diffs_rxz_mult_1(general_result_path)
    case = "30l_oh"
    for folder in [i for i in readdir(general_result_path) ]
        if isdir(joinpath(general_result_path, folder))
            if occursin(case, folder) && any([occursin("imp_est_scenario_1__power_mult_1.0", i) for i in readdir(joinpath(general_result_path, folder))])
                for whatt in ["Xc", "Rc", "Zc"]
                    ylims = whatt == "Xc" ? (-0.4, 30.4) : ()
                    if whatt == "Xc" yticks = [i for i in 0:5:30] end 
                    if whatt == "Rc" yticks = [i for i in 0:20:100] end  
                    if whatt == "Zc" yticks = [i for i in 0:10:60] end  
                    p =  cumulative_zrx_boxplot_crossconstraint_paper(general_result_path, case, 1.0; whatt=whatt, ytickz = yticks, ylims = ylims)
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm1.png")
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm1.pdf")
                end
            end
        end
    end
end

function plot_eu_ug_cumulative_diffs_rxz_mult_1(general_result_path)
    case = "eulvtf_ug"
    for folder in [i for i in readdir(general_result_path) ]
        if isdir(joinpath(general_result_path, folder))
            if occursin(case, folder) && any([occursin("imp_est_scenario_1__power_mult_1.0", i) for i in readdir(joinpath(general_result_path, folder))])
                for whatt in ["Xc", "Rc", "Zc"]
                    ylims = whatt == "Xc" ? (-0.4, 30.4) : ()
                    if whatt == "Xc" yticks = [i for i in 0:5:30] end 
                    if whatt == "Rc" yticks = [i for i in 0:20:100] end  
                    if whatt == "Zc" yticks = [i for i in 0:10:60] end  
                    p =  cumulative_zrx_boxplot_crossconstraint_paper(general_result_path, case, 1.0; whatt=whatt, ytickz = yticks, ylims = ylims)
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm1.png")
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm1.pdf")
                end
            end
        end
    end
end

function plot_eu_oh_cumulative_diffs_rxz_mult_1(general_result_path)
    case = "eulvtf_oh"
    for folder in [i for i in readdir(general_result_path) ]
        if isdir(joinpath(general_result_path, folder))
            if occursin(case, folder) && any([occursin("imp_est_scenario_1__power_mult_1.0", i) for i in readdir(joinpath(general_result_path, folder))])
                for whatt in ["Xc", "Rc", "Zc"]
                    ylims = whatt == "Xc" ? (-0.4, 30.4) : ()
                    if whatt == "Xc" yticks = [i for i in 0:5:30] end 
                    if whatt == "Rc" yticks = [i for i in 0:20:100] end  
                    if whatt == "Zc" yticks = [i for i in 0:10:60] end  
                    p =  cumulative_zrx_boxplot_crossconstraint_paper(general_result_path, case, 1.0; whatt=whatt, ytickz = yticks, ylims = ylims)
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm1.png")
                    _SP.savefig(general_result_path*"/$(case)_cumulative_$(whatt)_diff_boxplot_constr_paper_pm1.pdf")
                end
            end
        end
    end
end

function cumulative_zrx_boxplot_crossconstraint_paper(general_result_path::String, case::String, power_mult::Float64; whatt::String="Zc", ytickz::Vector=[], ylims=())
    
    folders = [f for f in readdir(general_result_path) if occursin(case, f) && !occursin(".png", f) && !occursin(".pdf", f)]
    est_dicts = []
    legend_filenames = []
    for folder in folders 
        if isdir(joinpath(general_result_path, folder))
            for file in readdir(joinpath(general_result_path, folder))
                if occursin("imp_est", file) && occursin("power_mult_$(power_mult)", file)
                    if occursin("no_restriction", folder) push!(legend_filenames, L"\textrm{No  \, \, rest.}") end
                    if occursin("most_restricted", folder) push!(legend_filenames, L"A_p+\mathcal{G}  \, \, \textrm{rest.}") end
                    if occursin("cross", folder) push!(legend_filenames, L"A_p  \, \, \textrm{rest.}") end
                    if (occursin("horizontal", folder) || occursin("square", folder)) push!(legend_filenames, L"\mathcal{G}  \, \, \textrm{rest.}") end
                    push!(est_dicts, joinpath(general_result_path, folder, file))
                end
            end
        end
    end

    starter_est = est_dicts[1]

    true_impedance_dict = JSON.parsefile(replace(starter_est, "est" => "true"))
    ylab = L"(|%$whatt^{\textrm{true}} - %$whatt^{\textrm{est}}|)/%$whatt^{\textrm{true}} \times \textrm{100} \, \, [\%]"
    p = _SP.boxplot([abs(v["$(whatt)_true"]-JSON.parsefile(starter_est)[k]["$(whatt)_est"])/v["$(whatt)_true"]*100 for (k,v) in true_impedance_dict],
                    ylabel = ylab, 
                    legend = false, color = :lightgrey, xtickfontsize=11,ytickfontsize=11, ylabelfontsize=14)

    for ed in est_dicts
        if ed != starter_est
           _SP.boxplot!([abs(v["$(whatt)_true"]-JSON.parsefile(ed)[k]["$(whatt)_est"])/v["$(whatt)_true"]*100 for (k,v) in true_impedance_dict],
            xticks = ((1:4),legend_filenames), color = :lightgrey
           )
        end
    end

    if !isempty(ytickz)
        _SP.plot!(yticks=(ytickz, [L"%$i" for i in ytickz]))
    end 

    if !isempty(ylims) _SP.plot!(ylims=ylims) end

    return p
end

########################################################
#########                                      #########
#########             LINECODE PLOTS           #########
#########                                      #########
########################################################

function linecode_plots(general_result_path::String, case::String, feeder_id::String, power_mult::Float64)
    master = occursin("oh", case) ? "Master_oh.dss" : "Master_ug.dss"
    eng = _PMD.parse_file(_IMP.NTW_DATA_DIR*"/"*feeder_id*"/"*master, data_model = _PMD.ENGINEERING, transformations=[_PMD.transform_loops!,_PMD.remove_all_bounds!])

    for linecode in keys(eng["linecode"])

        mat_rs = eng["linecode"][linecode]["rs"].*1000
        mat_xs = eng["linecode"][linecode]["xs"].*1000

        wires = size(mat_rs)[1]

        dict = Dict(
            "30l_ug_cross_only" => L"A_p \textrm{rest.}",
            "30l_ug_most_restricted" => L"A_p+\mathcal{G} \textrm{rest.}",
            "30l_ug_no_restriction"  => L"\textrm{No rest.}",
            "30l_ug_squared_only"   => L"\mathcal{G} \textrm{rest.}",
            "30l_oh_horizontal_only" => L"A_p \textrm{rest.}",
            "30l_oh_most_restricted" => L"A_p+\mathcal{G} \textrm{rest.}",
            "30l_oh_no_restriction"  => L"\textrm{No rest.}",
            "30l_oh_squared_only"   => L"\mathcal{G} \textrm{rest.}",
            "eulvtf_ug_cross_only" => L"A_p \textrm{rest.}",
            "eulvtf_ug_most_restricted" => L"A_p+\mathcal{G} \textrm{rest.}",
            "eulvtf_ug_no_restriction"  => L"\textrm{No rest.}",
            "eulvtf_ug_squared_only"   => L"\mathcal{G} \textrm{rest.}",
            "eulvtf_oh_horizontal_only" => L"A_p \textrm{rest.}",
            "eulvtf_oh_most_restricted" => L"A_p+\mathcal{G} \textrm{rest.}",
            "eulvtf_oh_no_restriction"  => L"\textrm{No rest.}",
            "eulvtf_oh_squared_only"   => L"\mathcal{G} \textrm{rest.}",
        )

        if wires == 4
            _SP.scatter([1:14], [mat_rs[1,1], mat_rs[2,2], mat_rs[3,3], mat_rs[4,4], mat_xs[1,1], mat_xs[2,2], mat_xs[3,3], mat_xs[4,4],
                                    mat_xs[1,2], mat_xs[1,3], mat_xs[1,4], mat_xs[2,3], mat_xs[2,4], mat_xs[3,4]], 
                                    xticks=((1:14), [L"R\_{aa}", L"R\_{bb}", L"R\_{cc}", L"R\_{nn}", L"X\_{aa}", L"X\_{bb}", L"X\_{cc}", L"X\_{nn}", L"X\_{ab}", L"X\_{ac}", L"X\_{an}", L"X\_{bc}", L"X\_{bn}", L"X\_{cn}"]), 
                                    legend=:bottomright, ylabel=L"[\Omega/km]", label = L"\textrm{True}", xtickfontsize=11,ytickfontsize=11, ylabelfontsize=14)
        else
            _SP.scatter([1:5], [mat_rs[1,1], mat_rs[2,2], mat_xs[1,1], mat_xs[2,2], mat_xs[1,2]], xticks=((1:5), [L"R\_{pp}", L"R\_{nn}", L"X\_{pp}", L"X\_{nn}", L"R\_{pn}"]), legend=:topright, label = L"\textrm{True}",
            xtickfontsize=11,ytickfontsize=11, ylabelfontsize=14)
        end
        folders = [f for f in readdir(general_result_path) if occursin(case, f)]

        for folder in folders 
            if isdir(joinpath(general_result_path, folder))
                for file in readdir(joinpath(general_result_path, folder))
                    if occursin("linecode_results", file) && occursin("power_mult_$(power_mult)", file)
                        df = CSV.read(joinpath(general_result_path, folder, file), _DF.DataFrame, ntasks = 1)
                        if linecode ∈ df.linecode_name
                            x_est = filter(x->x.linecode_name.==linecode, df).x_est[1][2:end-2]
                            r_est = filter(x->x.linecode_name.==linecode, df).r_est[1][2:end-2]
                            x_est = replace(x_est, ";" => "")
                            r_est = replace(r_est, ";" => "")
                            xs = parse.(Float64, split(x_est, " "))
                            rs = parse.(Float64, split(r_est, " "))
                            xs = reshape(xs, (wires,wires))
                            rs = reshape(rs, (wires,wires))
                            if wires == 4
                                _SP.scatter!(Vector(1:14), [rs[1,1], rs[2,2], rs[3,3], rs[4,4], xs[1,1], xs[2,2], xs[3,3], xs[4,4], xs[1,2], xs[1,3], xs[1,4], xs[2,3], xs[2,4], xs[3,4]], label = "$folder")
                            else
                                _SP.scatter!(Vector(1:5), [rs[1], rs[4], xs[1], xs[4], xs[2]], label = "$folder")
                            end
                        end
                    end
                end
            end
        end
        _SP.savefig(joinpath(general_result_path, "linecode_$(split(linecode, "/")[1])_results_case_$(case)_power_mult_$(power_mult).png"))
    end
end