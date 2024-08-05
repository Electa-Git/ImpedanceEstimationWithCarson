##############################################################################
### This file contains plotting functions that would                         
### recreate the papers' figures.
### The functions in this file are not particularly efficient / clean, sorry.
##############################################################################

import StatsPlots as _SP
using LaTeXStrings
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
import DataDrivenImpedanceEstimationWithCarson as _IMP
import JSON

general_result_path = raw"C:\Users\mvanin\OneDrive - KU Leuven\Desktop\repos\DataDrivenImpedanceEstimationWithCarson\results_ma27/"

########################################################
#####                                              #####
#####         POWER FLOW VALIDATION PLOTS          #####
#####                                              #####
########################################################

#########################
### PLOTS FOR PAPER 1 ###
#########################

# powerflow validation
plot_30l_ug_pf_validation(general_result_path)
plot_eu_ug_pf_validation(general_result_path)

# cumulative impedance
plot_30l_ug_cumulative_diffs_rxz_mult_1(general_result_path)
plot_30l_ug_cumulative_diffs_rxz_mult_3(general_result_path)

# linecode plots
linecode_plots(general_result_path, "30l_ug")
linecode_plots(general_result_path, "eulvtf_ug")

# carson input plots
carson_input_plots(general_result_path, "30l_ug")
carson_input_plots(general_result_path, "eulvtf_ug")

# branch length plots
branch_length_plots(general_result_path, "30l_ug")
branch_length_plots(general_result_path, "eulvtf_ug")

#########################
### PLOTS FOR PAPER 2 ###
#########################

# powerflow validation
# plot_30l_oh_pf_validation(general_result_path, power_mult=2.0)
plot_eu_oh_pf_validation(general_result_path, power_mult=2.0, ylims = (-0.1, 2.2))
# plot_30l_oh_noground_pf_validation(general_result_path, power_mult=2.0)
plot_eu_oh_noground_pf_validation(general_result_path, power_mult=2.0, ylims = (-0.1, 2.2))

# cumulative impedance
# plot_30l_oh_cumulative_diffs_rxz_mult_1(general_result_path, power_mult=2.0, ylims=(-0.5, 50))
plot_eu_oh_cumulative_diffs_rxz_mult_1(general_result_path, power_mult=2.0, ylims=(-0.5, 25))
# plot_30l_oh_noground_cumulative_diffs_rxz_mult_1(general_result_path, power_mult=2.0,ylims=(-0.5, 50))
plot_eu_oh_cumulative_diffs_rxz_noground(general_result_path, power_mult=2.0, ylims=(-0.5, 25))

# linecode plots
# linecode_plots(general_result_path, "30l_oh", power_mult=2.0)
# linecode_plots(general_result_path, "30l_oh_noground", power_mult=2.0)
linecode_plots(general_result_path, "eulvtf_oh"; feeder_id = "eulvtf", power_mult=2.0)
linecode_plots(general_result_path, "eulvtf_oh_noground"; feeder_id = "eulvtf", power_mult=2.0)

# carson input plots
# carson_input_plots(general_result_path, "30l_oh", power_mult=2.0)
# carson_input_plots(general_result_path, "30l_oh_noground", power_mult=2.0)
carson_input_plots(general_result_path, "eulvtf_oh", power_mult=2.0)
carson_input_plots(general_result_path, "eulvtf_oh_noground", power_mult=2.0)

# branch length plots
# branch_length_plots(general_result_path, "30l_oh", power_mult=2.0)
branch_length_plots(general_result_path, "eulvtf_oh", power_mult=2.0)
# branch_length_plots(general_result_path, "30l_oh_noground", power_mult=2.0)
branch_length_plots(general_result_path, "eulvtf_oh_noground", power_mult=2.0)

# shunt plots
# shunt_plots_scatter(general_result_path, "30l_oh", power_mult=2.0)
shunt_plots_scatter(general_result_path, "eulvtf_oh",power_mult=2.0)

#################### PLOTS with no NOISE ###################



#######################################################################################################
#######################################################################################################
#######################################################################################################

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

function plot_30l_oh_pf_validation(general_result_path; power_mult::Float64=1.0)
    case = "30l_oh"
    folders = ["30l_oh_cross_only", "30l_oh_horizontal_only", "30l_oh_most_restricted", "30l_oh_no_restriction"]
    powerflow_validation_paperplot(general_result_path, folders, power_mult=power_mult, ytickz = [0., 0.5, 1., 1.5, 2., 2.5])
    _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper_pm_$power_mult.png")
    _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper_pm_$power_mult.pdf")
end

function plot_30l_oh_noground_pf_validation(general_result_path; power_mult::Float64=1.0)
    case = "30l_oh_noground"
    folders = ["30l_oh_noground_cross_only", "30l_oh_noground_horizontal_only", "30l_oh_noground_most_restricted", "30l_oh_noground_no_restriction"]
    powerflow_validation_paperplot(general_result_path, folders, power_mult = power_mult, ytickz = [0., 0.5, 1., 1.5, 2., 2.5])
    _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper_pm_$power_mult.png")
    _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper_pm_$power_mult.pdf")
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

function plot_eu_oh_pf_validation(general_result_path; power_mult::Float64=1.0, ylims=(-0.1, 2.0))
    case = "eulvtf_oh"
    folders = ["eulvtf_oh_cross_only", "eulvtf_oh_horizontal_only", "eulvtf_oh_most_restricted", "eulvtf_oh_no_restriction"]
    p = powerflow_validation_paperplot(general_result_path, folders, power_mult=power_mult, ytickz = [0., 0.5, 1., 1.5, 2., 2.5])
    _SP.ylims!(ylims)
    _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper_power_mult_$(power_mult).png")
    _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper_power_mult_$(power_mult).pdf")
end

function plot_eu_oh_noground_pf_validation(general_result_path; power_mult::Float64=1.0, ylims = (-0.1, 2.0))
    case = "eulvtf_oh_noground"
    folders = ["eulvtf_oh_noground_cross_only", "eulvtf_oh_noground_horizontal_only", "eulvtf_oh_noground_most_restricted", "eulvtf_oh_noground_no_restriction"]
    p = powerflow_validation_paperplot(general_result_path, folders, power_mult = power_mult, ytickz = [0., 0.5, 1., 1.5, 2., 2.5])
    _SP.ylims!(ylims)
    _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper_$(power_mult).png")
    _SP.savefig(general_result_path*"/$(case)_pf_validation_constr_paper_$(power_mult).pdf")
end


function powerflow_validation_paperplot(general_result_path::String, case::String; ytickz::Vector=[], power_mult::Float64=1.0)

    folders = [f for f in readdir(general_result_path) if occursin(case, f) && !(occursin(".png", f) || occursin(".pdf", f))]
    est_dfs = []
    real_dfs = []
    legend_filenames = []
    for folder in folders 
        if isdir(joinpath(general_result_path, folder))
            for file in readdir(joinpath(general_result_path, folder))
                if occursin("pf_validation_real", file) && occursin("vm_power_mult_$(power_mult)", file)
                    real_pf = CSV.read(joinpath(general_result_path, folder, "pf_validation_real_vm_power_mult_$(power_mult).csv"),_DF.DataFrame, ntasks = 1)
                    est_pf  = CSV.read(joinpath(general_result_path, folder, "pf_validation_est_vm_power_mult_$(power_mult).csv"),_DF.DataFrame, ntasks = 1)
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

    # adds the other cases > 2
    display(est_dfs)
    if length(est_dfs) > 2
        for i in 2:length(est_dfs)
            diff_df = deepcopy(est_dfs[i])
            for row in 1:size(diff_df)[1]
                diff_df[row:row,1:end-2] .= abs.(est_dfs[i][row:row,1:end-2] .- real_dfs[i][row:row,1:end-2])
            end
            all_diffs = [i for i in vec(Matrix(diff_df)) if (i isa Float64)]
            p = _SP.boxplot!(all_diffs*240, color = :lightgrey, xticks = ((1:4),legend_filenames))
        end
    end
       
    if !isempty(ytickz)
        _SP.plot!(yticks=(ytickz, [L"%$i" for i in ytickz]))
    end 

    return p
end

function powerflow_validation_paperplot(general_result_path::String, folders::Vector{String}; power_mult::Float64=1.0, ytickz::Vector=[])

    est_dfs = []
    real_dfs = []
    legend_filenames = []
    for folder in folders 
        for file in readdir(joinpath(general_result_path, folder))
            if occursin("pf_validation_real_vm", file) # && occursin("vm_noground_power_mult_1.0", file)
                real_pf = CSV.read(joinpath(general_result_path, folder, file),_DF.DataFrame, ntasks = 1)
                est_file = [f for f in readdir(joinpath(general_result_path, folder)) if (occursin("pf_validation_est_vm", f) && occursin("$power_mult", f))][1]
                est_pf  = CSV.read(joinpath(general_result_path, folder, est_file),_DF.DataFrame, ntasks = 1)
                if occursin("no_restriction", folder) push!(legend_filenames, L"\textrm{No  \, \, rest.}") end
                if occursin("most_restricted", folder) push!(legend_filenames, L"A_p+\mathcal{G}  \, \, \textrm{rest.}") end
                if occursin("cross", folder) push!(legend_filenames, L"A_p  \, \, \textrm{rest.}") end
                if (occursin("horizontal", folder) || occursin("square", folder)) push!(legend_filenames, L"\mathcal{G}  \, \, \textrm{rest.}") end
                push!(est_dfs, est_pf)
                push!(real_dfs, real_pf)
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

function powerflow_validation_paperplot_nogroundcase(general_result_path::String, case::String; ytickz::Vector=[])

    est_dfs = []
    real_dfs = []
    legend_filenames = []
    folders = case == "eulvtf_oh_noground" ? ["eulvtf_oh_noground_cross_only", "eulvtf_oh_noground_horizontal_only", "eulvtf_oh_noground_most_restricted", "eulvtf_oh_noground_no_restriction"] : ["30l_oh_noground_cross_only", "30l_oh_noground_horizontal_only", "30l_oh_noground_most_restricted", "30l_oh_noground_no_restriction"] 
    for folder in folders
        if isdir(joinpath(general_result_path, folder))
            for file in readdir(joinpath(general_result_path, folder))
                if occursin("pf_validation_real_vm_noground", file) && occursin("_power_mult_1.0", file)
                    real_pf = CSV.read(joinpath(general_result_path, folder, "pf_validation_real_vm_noground_power_mult_1.0.csv"),_DF.DataFrame, ntasks = 1)
                    est_pf  = CSV.read(joinpath(general_result_path, folder, "pf_validation_est_vm_noground_power_mult_1.0.csv"),_DF.DataFrame, ntasks = 1)
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

function plot_30l_oh_cumulative_diffs_rxz_mult_1(general_result_path; power_mult=1.0, ylims = ())
    folders = ["30l_oh_cross_only", "30l_oh_horizontal_only", "30l_oh_most_restricted", "30l_oh_no_restriction"] 
    for whatt in ["Xc", "Rc", "Zc"]
        # ylims = whatt == "Xc" ? (-0.4, 30.4) : ()
        if whatt == "Xc" yticks = [i for i in 0:5:30] end 
        if whatt == "Rc" yticks = [i for i in 0:20:100] end  
        if whatt == "Zc" yticks = [i for i in 0:10:60] end  
        p =  cumulative_zrx_boxplot_crossconstraint_paper(general_result_path, folders, power_mult; whatt=whatt, ytickz = yticks, ylims = ylims)
        _SP.savefig(general_result_path*"/30l_oh_cumulative_$(whatt)_diff_boxplot_constr_paper_pm_$power_mult.png")
        _SP.savefig(general_result_path*"/30l_oh_cumulative_$(whatt)_diff_boxplot_constr_paper_pm_$power_mult.pdf")
    end
end

function plot_30l_oh_noground_cumulative_diffs_rxz_mult_1(general_result_path;  power_mult=1.0, ylims=())
    folders = ["30l_oh_noground_cross_only", "30l_oh_noground_horizontal_only", "30l_oh_noground_most_restricted", "30l_oh_noground_no_restriction"] 
    for folder in folders
        if  any([occursin("imp_est_scenario_1__power_mult_1.0", i) for i in readdir(joinpath(general_result_path, folder))])
            for whatt in ["Xc", "Rc", "Zc"]
                # ylims = whatt == "Xc" ? (-0.4, 30.4) : ()
                if whatt == "Xc" yticks = [i for i in 0:5:30] end 
                if whatt == "Rc" yticks = [i for i in 0:20:100] end  
                if whatt == "Zc" yticks = [i for i in 0:10:60] end  
                p =  cumulative_zrx_boxplot_crossconstraint_paper(general_result_path, folders, power_mult; whatt=whatt, ytickz = yticks, ylims = ylims)
                _SP.savefig(general_result_path*"/30l_oh_noground_cumulative_$(whatt)_diff_boxplot_constr_paper_pm_$power_mult.png")
                _SP.savefig(general_result_path*"/30l_oh_noground_cumulative_$(whatt)_diff_boxplot_constr_paper_pm_$power_mult.pdf")
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

function plot_eu_oh_cumulative_diffs_rxz_mult_1(general_result_path; power_mult::Float64=1., ylims=())
    folders = ["eulvtf_oh_cross_only", "eulvtf_oh_horizontal_only", "eulvtf_oh_most_restricted", "eulvtf_oh_no_restriction"] 
    for folder in folders
        if  any([occursin("imp_est_scenario_1__power_mult_1.0", i) for i in readdir(joinpath(general_result_path, folder))])
            for whatt in ["Xc", "Rc", "Zc"]
                # if whatt == "Xc" ylims = (-0.4, 30.4) end
                if whatt == "Xc" yticks = [i for i in 0:5:30] end 
                if whatt == "Rc" yticks = [i for i in 0:20:100] end  
                if whatt == "Zc" yticks = [i for i in 0:10:60] end  
                p =  cumulative_zrx_boxplot_crossconstraint_paper(general_result_path, folders, power_mult; whatt=whatt, ytickz = yticks, ylims = ylims)
                _SP.savefig(general_result_path*"/eulvtf_oh_cumulative_$(whatt)_diff_boxplot_constr_paper_pm_$power_mult.png")
                _SP.savefig(general_result_path*"/eulvtf_oh_cumulative_$(whatt)_diff_boxplot_constr_paper_pm_$power_mult.pdf")
            end
        end
    end
end

function plot_eu_oh_cumulative_diffs_rxz_noground(general_result_path; power_mult::Float64=1.0, ylims=())
    folders = ["eulvtf_oh_noground_cross_only", "eulvtf_oh_noground_horizontal_only", "eulvtf_oh_noground_most_restricted", "eulvtf_oh_noground_no_restriction"] 
    for folder in folders
        if  any([occursin("imp_est_scenario_1__power_mult_$power_mult", i) for i in readdir(joinpath(general_result_path, folder))])
            for whatt in ["Xc", "Rc", "Zc"]
                # if whatt == "Xc" ylims = (-0.4, 30.4) end
                if whatt == "Xc" yticks = [i for i in 0:5:30] end 
                if whatt == "Rc" yticks = [i for i in 0:20:100] end  
                if whatt == "Zc" yticks = [i for i in 0:10:60] end  
                p =  cumulative_zrx_boxplot_crossconstraint_paper(general_result_path, folders, power_mult; whatt=whatt, ytickz = yticks, ylims = ylims)
                _SP.savefig(general_result_path*"/eulvtf_oh_noground_cumulative_$(whatt)_diff_boxplot_constr_paper_pm_$power_mult.png")
                _SP.savefig(general_result_path*"/eulvtf_oh_noground_cumulative_$(whatt)_diff_boxplot_constr_paper_pm_$power_mult.pdf")
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

function cumulative_zrx_boxplot_crossconstraint_paper(general_result_path::String, folders::Vector{String}, power_mult::Float64; whatt::String="Zc", ytickz::Vector=[], ylims=())
    
    est_dicts = []
    legend_filenames = []
    for folder in folders
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

    starter_est = est_dicts[1]

    true_impedance_dict = JSON.parsefile(replace(starter_est, "_est" => "_true"))
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

function linecode_plots(general_result_path::String, case::String; feeder_id::String="30load-feeder", power_mult::Float64=1.0)
    master = occursin("oh", case) ? "Master_oh.dss" : "Master_ug.dss"
    eng = _PMD.parse_file(_IMP.NTW_DATA_DIR*"/"*feeder_id*"/"*master, data_model = _PMD.ENGINEERING, transformations=[_PMD.transform_loops!,_PMD.remove_all_bounds!])

    yticks_dict = Dict(
        "pluto"                                    => [i for i in 0.2:0.1:0.8],
        "uglv_120cu_xlpe/nyl/pvc_ug_4w_bundled"    => [i for i in 0.25:0.1:0.75],
        "ugsc_16cu_xlpe/nyl/pvc_ug_2w_bundled"     => vcat([0.14], [i for i in 0.75:0.1:0.9], [1.5, 2.2]),
        "hydrogen"                                 => vcat([i for i in 0.13:0.02:0.25], [i for i in 0.35:0.05:0.5], [i for i in 0.6:0.1:0.8]),
        "ugsc_25cu_xlpe/nyl/pvc_ug_2w_bundled"     => vcat([i for i in 1.:0.5:4.1], [i for i in 0.6:0.1:1]),
        "uglv_185al_xlpe/nyl/pvc_ug_4w_bundled"    => [i for i in 0.2:0.4:3.5],
        "abc2x16_lv_oh_2w_bundled"                 => vcat([i for i in 0.6:0.1:1.0],[i for i in 2.2:0.4:4.8]),
        "ugsc_16al_xlpe/pvc_ug_2w_bundled"         => vcat([i for i in 0.2:0.4:3.5], [i for i in 0.7:0.1:0.9]),
        "uglv_185al_xlpe/nyl/pvc_ug_2w_bundled"    => vcat([0.25], [i for i in 0.55:0.1:0.85], [2.7, 3.1, 3.4]),
        "tw2x16_lv_oh_2w_bundled"                  => vcat([i for i in 0.55:0.1:0.95],[i for i in 2.4:0.3:4.0]),
        "uglv_240al_xlpe/nyl/pvc_ug_4w_bundled"    => vcat([i for i in 0.2:0.025:0.25], [i for i in 0.3:0.1:0.8])
    )

    for linecode in keys(eng["linecode"])

        mat_rs = eng["linecode"][linecode]["rs"].*1000
        mat_xs = eng["linecode"][linecode]["xs"].*1000

        wires = size(mat_rs)[1]

        legend_dict = Dict(
            "30l_ug_cross_only" => L"A_p \, \, \textrm{rest.}",
            "30l_ug_most_restricted" => L"A_p+\mathcal{G} \, \, \textrm{rest.}",
            "30l_ug_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "30l_ug_squared_only"   => L"\mathcal{G} \, \, \textrm{rest.}",
            "30l_oh_cross_only" => L"A_p \, \, \textrm{rest.}",
            "30l_oh_most_restricted" => L"A_p+\mathcal{G} \, \, \textrm{rest.}",
            "30l_oh_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "30l_oh_horizontal_only"   => L"\mathcal{G}  \, \, \textrm{rest.}",
            "30l_oh_noground_cross_only" => L"A_p \, \, \textrm{rest.}",
            "30l_oh_noground_most_restricted" => L"A_p+\mathcal{G} \, \, \textrm{rest.}",
            "30l_oh_noground_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "30l_oh_noground_horizontal_only"   => L"\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_ug_cross_only" => L"A_p  \, \, \textrm{rest.}",
            "eulvtf_ug_most_restricted" => L"A_p+\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_ug_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "eulvtf_ug_squared_only"   => L"\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_oh_cross_only" => L"A_p  \, \, \textrm{rest.}",
            "eulvtf_oh_most_restricted" => L"A_p+\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_oh_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "eulvtf_oh_horizontal_only"   => L"\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_oh_noground_cross_only" => L"A_p  \, \, \textrm{rest.}",
            "eulvtf_oh_noground_most_restricted" => L"A_p+\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_oh_noground_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "eulvtf_oh_noground_horizontal_only"   => L"\mathcal{G}  \, \, \textrm{rest.}"
        )

        markershape_dict = Dict(
            "30l_ug_cross_only" => :circle,
            "30l_ug_most_restricted" => :diamond,
            "30l_ug_no_restriction"  => :rect,
            "30l_ug_squared_only"   =>  :utriangle,
            "30l_oh_cross_only" => :circle,
            "30l_oh_most_restricted" => :diamond,
            "30l_oh_no_restriction"  => :rect,
            "30l_oh_horizontal_only"   =>  :utriangle,
            "30l_oh_noground_cross_only" => :circle,
            "30l_oh_noground_most_restricted" => :diamond,
            "30l_oh_noground_no_restriction"  => :rect,
            "30l_oh_noground_horizontal_only"   =>  :utriangle,
            "eulvtf_ug_cross_only" => :circle,
            "eulvtf_ug_most_restricted" => :diamond,
            "eulvtf_ug_no_restriction"  => :rect,
            "eulvtf_ug_squared_only"   =>  :utriangle,
            "eulvtf_oh_cross_only" => :circle,
            "eulvtf_oh_most_restricted" => :diamond,
            "eulvtf_oh_no_restriction"  => :rect,
            "eulvtf_oh_horizontal_only"   =>  :utriangle,
            "eulvtf_oh_noground_cross_only" => :circle,
            "eulvtf_oh_noground_most_restricted" => :diamond,
            "eulvtf_oh_noground_no_restriction"  => :rect,
            "eulvtf_oh_noground_horizontal_only"   =>  :utriangle
        )

        markersize_dict = Dict(
            "30l_ug_cross_only" =>8,
            "30l_ug_most_restricted" =>4,
            "30l_ug_no_restriction"  =>5,
            "30l_ug_squared_only"   =>5,
            "30l_oh_cross_only" =>8,
            "30l_oh_most_restricted" =>4,
            "30l_oh_no_restriction"  =>5,
            "30l_oh_horizontal_only"   =>5,
            "30l_oh_noground_cross_only" =>8,
            "30l_oh_noground_most_restricted" =>4,
            "30l_oh_noground_no_restriction"  =>5,
            "30l_oh_noground_horizontal_only"   =>5,
            "eulvtf_ug_cross_only" =>8,
            "eulvtf_ug_most_restricted" =>4,
            "eulvtf_ug_no_restriction"  =>5,
            "eulvtf_ug_squared_only"   =>4,
            "eulvtf_oh_cross_only" =>8,
            "eulvtf_oh_most_restricted" =>4,
            "eulvtf_oh_no_restriction"  =>5,
            "eulvtf_oh_horizontal_only"   =>4,
            "eulvtf_oh_noground_cross_only" =>8,
            "eulvtf_oh_noground_most_restricted" =>4,
            "eulvtf_oh_noground_no_restriction"  =>5,
            "eulvtf_oh_noground_horizontal_only"   =>4
        )

        if wires == 4
            p = _SP.scatter([1:14], [mat_rs[1,1], mat_rs[2,2], mat_rs[3,3], mat_rs[4,4], mat_xs[1,1], mat_xs[2,2], mat_xs[3,3], mat_xs[4,4],
                                    mat_xs[1,2], mat_xs[1,3], mat_xs[1,4], mat_xs[2,3], mat_xs[2,4], mat_xs[3,4]], 
                                    xticks=((1:14), [L"R_{aa}", L"R_{bb}", L"R_{cc}", L"R_{nn}", L"X_{aa}", L"X_{bb}", L"X_{cc}", L"X_{nn}", L"X_{ab}", L"X_{ac}", L"X_{an}", L"X_{bc}", L"X_{bn}", L"X_{cn}"]), 
                                    legend=:bottomright, ylabel=L"[\Omega/km]", label = L"\textrm{True}", xtickfontsize=13,ytickfontsize=12, ylabelfontsize=14, legendfontsize=12, markershape = :xcross, ms=6, color="black",  yscale=:log10)
        else
            legendloc = linecode == "ugsc_16cu_xlpe/nyl/pvc_ug_2w_bundled" ? :bottomright : :topright
            p = _SP.scatter([1:5], [mat_rs[1,1], mat_rs[2,2], mat_xs[1,1], mat_xs[2,2], mat_xs[1,2]], xticks=((1:5), [L"R_{pp}", L"R_{nn}", L"X_{pp}", L"X_{nn}", L"X_{pn}"]), legend=legendloc, 
                                    ylabel=L"[\Omega/km]", label = L"\textrm{True}", xtickfontsize=13,ytickfontsize=12, ylabelfontsize=14, legendfontsize=12, markershape = :xcross, ms=6, color="black",  yscale=:log10)
        end

        folders = [f for f in readdir(general_result_path) if occursin(case, f) && isdir(joinpath(general_result_path, f))]

        if length(folders) > 4 
            folders = filter(x->!occursin("noground", x), folders)
        end

        return_p = false
        for folder in folders 
            if isdir(joinpath(general_result_path, folder))
                for file in readdir(joinpath(general_result_path, folder))
                    if occursin("linecode_results", file) && occursin("power_mult_$(power_mult)", file) #&& occursin("ground", file))
                        df = CSV.read(joinpath(general_result_path, folder, file), _DF.DataFrame, ntasks = 1)
                        if linecode ∈ df.linecode_name
                            return_p = true
                            x_est = filter(x->x.linecode_name.==linecode, df).x_est[1][2:end-2]
                            r_est = filter(x->x.linecode_name.==linecode, df).r_est[1][2:end-2]
                            x_est = replace(x_est, ";" => "")
                            r_est = replace(r_est, ";" => "")
                            xs = parse.(Float64, split(x_est, " "))
                            rs = parse.(Float64, split(r_est, " "))
                            xs = reshape(xs, (wires,wires))
                            rs = reshape(rs, (wires,wires))
                            if wires == 4
                                _SP.scatter!(Vector(1:14), [rs[1,1], rs[2,2], rs[3,3], rs[4,4], xs[1,1], xs[2,2], xs[3,3], xs[4,4], xs[1,2], xs[1,3], xs[1,4], xs[2,3], xs[2,4], xs[3,4]], label = legend_dict["$folder"],
                                markershape = markershape_dict["$folder"],  color="black", ms = markersize_dict["$folder"], mc=:white)
                            else
                                _SP.scatter!(Vector(1:5), [rs[1], rs[4], xs[1], xs[4], xs[2]], label = legend_dict["$folder"], markershape = markershape_dict["$folder"],  color="black", ms = markersize_dict["$folder"], mc=:white)
                            end
                        end
                    end
                end
            end
        end
        _SP.yticks!(yticks_dict["$linecode"], [L"%$i" for i in yticks_dict["$linecode"]])
        if return_p
            _SP.savefig(joinpath(general_result_path, "linecode_$(split(linecode, "/")[1])_results_case_$(case)_power_mult_$(power_mult).png"))
            _SP.savefig(joinpath(general_result_path, "linecode_$(split(linecode, "/")[1])_results_case_$(case)_power_mult_$(power_mult).pdf"))
        end
    end
end

########################################################
#########                                      #########
#########         CARSON INPUT PLOTS           #########
#########                                      #########
########################################################

function carson_input_plots(general_result_path::String, case::String; power_mult::Float64=1.0)

    ground_truth = CSV.read(_IMP.DATA_DIR*"/linecode_library/linecode_geometries.csv", _DF.DataFrame, ntasks=1)

    yticks_dict = Dict(
        "pluto"                                    => [i for i in -100:20:100],
        "uglv_120cu_xlpe/nyl/pvc_ug_4w_bundled"    => [i for i in -20:5:15],
        "ugsc_16cu_xlpe/nyl/pvc_ug_2w_bundled"     => [-1200, 25],
        "hydrogen"                                 => [i for i in -100:20:100],
        "ugsc_25cu_xlpe/nyl/pvc_ug_2w_bundled"     => [-700, 40, 70],
        "uglv_185al_xlpe/nyl/pvc_ug_4w_bundled"    => [i for i in -20:20:100],
        "abc2x16_lv_oh_2w_bundled"                 => occursin("eulvtf", case) ? vcat([-190, -150], [-5, 30, 50, 70]) : vcat([-1200], [25]),
        "ugsc_16al_xlpe/pvc_ug_2w_bundled"         => [-950, 25],
        "uglv_185al_xlpe/nyl/pvc_ug_2w_bundled"    => [-300, 100],
        "tw2x16_lv_oh_2w_bundled"                  => vcat([-1200], [25]),
        "uglv_240al_xlpe/nyl/pvc_ug_4w_bundled"    => [i for i in -30:15:90]
    )

    for linecode in ground_truth.linecode_name

        lc = filter(x->x.linecode_name == linecode, ground_truth)
        wires =  count(==(','),lc.A_p_true[1]) == 3 ? 4 : 2

        legend_dict = Dict(
            "30l_ug_cross_only" => L"A_p \, \, \textrm{rest.}",
            "30l_ug_most_restricted" => L"A_p+\mathcal{G} \, \, \textrm{rest.}",
            "30l_ug_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "30l_ug_squared_only"   => L"\mathcal{G} \, \, \textrm{rest.}",
            "30l_oh_cross_only" => L"A_p \, \, \textrm{rest.}",
            "30l_oh_most_restricted" => L"A_p+\mathcal{G} \, \, \textrm{rest.}",
            "30l_oh_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "30l_oh_horizontal_only"   => L"\mathcal{G}  \, \, \textrm{rest.}",
            "30l_oh_noground_cross_only" => L"A_p \, \, \textrm{rest.}",
            "30l_oh_noground_most_restricted" => L"A_p+\mathcal{G} \, \, \textrm{rest.}",
            "30l_oh_noground_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "30l_oh_noground_horizontal_only"   => L"\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_ug_cross_only" => L"A_p  \, \, \textrm{rest.}",
            "eulvtf_ug_most_restricted" => L"A_p+\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_ug_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "eulvtf_ug_squared_only"   => L"\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_oh_cross_only" => L"A_p  \, \, \textrm{rest.}",
            "eulvtf_oh_most_restricted" => L"A_p+\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_oh_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "eulvtf_oh_horizontal_only"   => L"\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_oh_noground_cross_only" => L"A_p  \, \, \textrm{rest.}",
            "eulvtf_oh_noground_most_restricted" => L"A_p+\mathcal{G}  \, \, \textrm{rest.}",
            "eulvtf_oh_noground_no_restriction"  => L"\textrm{No} \, \, \textrm{rest.}",
            "eulvtf_oh_noground_horizontal_only"   => L"\mathcal{G}  \, \, \textrm{rest.}"
        )

        markershape_dict = Dict(
            "30l_ug_cross_only" => :circle,
            "30l_ug_most_restricted" => :diamond,
            "30l_ug_no_restriction"  => :rect,
            "30l_ug_squared_only"   =>  :utriangle,
            "30l_oh_cross_only" => :circle,
            "30l_oh_most_restricted" => :diamond,
            "30l_oh_no_restriction"  => :rect,
            "30l_oh_horizontal_only"   =>  :utriangle,
            "30l_oh_noground_cross_only" => :circle,
            "30l_oh_noground_most_restricted" => :diamond,
            "30l_oh_noground_no_restriction"  => :rect,
            "30l_oh_noground_horizontal_only"   =>  :utriangle,
            "eulvtf_ug_cross_only" => :circle,
            "eulvtf_ug_most_restricted" => :diamond,
            "eulvtf_ug_no_restriction"  => :rect,
            "eulvtf_ug_squared_only"   =>  :utriangle,
            "eulvtf_oh_cross_only" => :circle,
            "eulvtf_oh_most_restricted" => :diamond,
            "eulvtf_oh_no_restriction"  => :rect,
            "eulvtf_oh_horizontal_only"   =>  :utriangle,
            "eulvtf_oh_noground_cross_only" => :circle,
            "eulvtf_oh_noground_most_restricted" => :diamond,
            "eulvtf_oh_noground_no_restriction"  => :rect,
            "eulvtf_oh_noground_horizontal_only"   =>  :utriangle
        )

        markersize_dict = Dict(
            "30l_ug_cross_only" =>8,
            "30l_ug_most_restricted" =>4,
            "30l_ug_no_restriction"  =>5,
            "30l_ug_squared_only"   =>5,
            "30l_oh_cross_only" =>8,
            "30l_oh_most_restricted" =>4,
            "30l_oh_no_restriction"  =>5,
            "30l_oh_horizontal_only"   =>5,
            "30l_oh_noground_cross_only" =>8,
            "30l_oh_noground_most_restricted" =>4,
            "30l_oh_noground_no_restriction"  =>5,
            "30l_oh_noground_horizontal_only"   =>5,
            "eulvtf_ug_cross_only" =>8,
            "eulvtf_ug_most_restricted" =>4,
            "eulvtf_ug_no_restriction"  =>5,
            "eulvtf_ug_squared_only"   =>4,
            "eulvtf_oh_cross_only" =>8,
            "eulvtf_oh_most_restricted" =>4,
            "eulvtf_oh_no_restriction"  =>5,
            "eulvtf_oh_horizontal_only"   =>4,
            "eulvtf_oh_noground_cross_only" =>8,
            "eulvtf_oh_noground_most_restricted" =>4,
            "eulvtf_oh_noground_no_restriction"  =>5,
            "eulvtf_oh_noground_horizontal_only"   =>4
        )

        d_true = parse.(Float64, split(chop(lc.dist_true[1], head=1),','))
        A_true = parse.(Float64, split(chop(lc.A_p_true[1], head=1),','))

        p = _SP.plot(legend=:bottomleft, ylabel=L"\textrm{Relative} \, \, \textrm{error} \, \, [\%]", label = L"\textrm{True}", xtickfontsize=13,ytickfontsize=12, ylabelfontsize=14, legendfontsize=12)

        folders = [f for f in readdir(general_result_path) if occursin(case, f) && isdir(joinpath(general_result_path, f))]

        if length(folders) > 4 
            folders = filter(x->!occursin("noground", x), folders)
        end

        return_p = false
        for folder in folders 
            if isdir(joinpath(general_result_path, folder))
                for file in readdir(joinpath(general_result_path, folder))
                    if occursin("linecode_results", file) && occursin("power_mult_$(power_mult)", file)
                        df = CSV.read(joinpath(general_result_path, folder, file), _DF.DataFrame, ntasks = 1)
                        if linecode ∈ df.linecode_name
                            df = filter(x->x.linecode_name.==linecode, df)
                            return_p = true
                            d_est = parse.(Float64, split(chop(string(df.dist_est[1]), head=1),','))
                            A_est = parse.(Float64, split(chop(string(df.A_p_est[1]) , head=1),','))
                            if wires == 4
                                d_est = expand_distance_results(case, d_est)
                                p = _SP.scatter!([1:10], (vcat(A_true, d_true).-vcat(A_est, d_est))./vcat(A_true, d_true)*100, 
                                    xticks=([i for i in 1:10], [L"A_{a}", L"A_{b}", L"A_{c}", L"A_{n}", L"D_{ab}", L"D_{ac}", L"D_{bc}", L"D_{an}", L"D_{bn}", L"D_{cn}"]), 
                                    legend=:bottomleft, label = legend_dict["$folder"], markershape = markershape_dict["$folder"],  color="black", ms = markersize_dict["$folder"], mc=:white)
                            else
                                p = _SP.scatter!([0.5:0.5:1.5], (vcat(A_true, d_true).-vcat(A_est, d_est))./vcat(A_true, d_true)*100,  xticks=([i for i in 0.5:0.5:1.5], [L"A_{p}", L"A_{n}", L"D_{pn}"]), 
                                    legend=:topright, label = legend_dict["$folder"], markershape = markershape_dict["$folder"],  color="black", ms = markersize_dict["$folder"], mc=:white)
                            end
                        end
                    end
                end
            end
        end
        _SP.yticks!(yticks_dict["$linecode"], [L"%$i" for i in yticks_dict["$linecode"]])
        if occursin("pluto", linecode) || occursin("hydrogen", linecode)
            _SP.plot!(ylims=(-100, 100), yscale=:log10)
        end
        if return_p
            _SP.savefig(joinpath(general_result_path, "carsoninput_$(split(linecode, "/")[1])_case_$(case)_power_mult_$(power_mult).png"))
            _SP.savefig(joinpath(general_result_path, "carsoninput_$(split(linecode, "/")[1])_case_$(case)_power_mult_$(power_mult).pdf"))
        end
    end
end

function expand_distance_results(case::String, d_est::Vector{Float64})
    if length(d_est) < 6 # this is a reduced distance variable case, needs to be mapped back to 6x vector
        if occursin("ug", case)
            d_exp = [d_est[1], d_est[1], d_est[1]*sqrt(2), d_est[1]/sqrt(2)+d_est[2], sqrt(d_est[2]^2+d_est[1]^2/2), sqrt(d_est[2]^2+d_est[1]^2/2)]
        else # overhead line
            d_exp = [d_est[1], d_est[1]+d_est[2], d_est[2], sum(d_est), d_est[2]+d_est[3], d_est[3]]
        end
        return d_exp
    else
        return d_est
    end
end

########################################################
#########                                      #########
#########      BRANCH LENGTH PLOTS             #########
#########                                      #########
########################################################

function branch_length_plots(general_result_path::String, case::String; power_mult::Float64=1.0)

    profiles = CSV.read(_IMP.DATA_DIR*"/profiles.csv", _DF.DataFrame, ntasks = 1)[1:2,:]
    oh_or_ug = occursin("oh", case) ? "oh" : "ug"
    feeder_name = occursin("30", case) ? "30load-feeder" : "eulvtf"

    data, eng, z_pu = prepare_math_eng_data(profiles; feeder_name = feeder_name, oh_or_ug = oh_or_ug);
    
    service_branches = [b for (b,br) in data["branch"] if length(br["t_connections"])==2]
    main_branches = [b for (b,br) in data["branch"] if b ∉ service_branches]

    ps = _SP.plot(legend=false, ylabel=L"\textrm{Relative} \, \, \textrm{error} \, \, [\%] - \textrm{Service} \, \, \textrm{cable}", label = L"\textrm{True}", xtickfontsize=13,ytickfontsize=12, ylabelfontsize=14, legendfontsize=12)
    pm = _SP.plot(legend=false, ylabel=L"\textrm{Relative} \, \, \textrm{error} \, \, [\%] - \textrm{Main} \, \, \textrm{cable}" , label = L"\textrm{True}", xtickfontsize=13,ytickfontsize=12, ylabelfontsize=14, legendfontsize=12)

    folders = [f for f in readdir(general_result_path) if occursin(case, f) && isdir(joinpath(general_result_path, f))]

    if length(folders) > 4 
        folders = filter(x->!occursin("noground", x), folders)
    end

    legend_filenames = []
    xposition = 1
    for folder in folders
        if occursin("no_restriction", folder) push!(legend_filenames, L"\textrm{No  \, \, rest.}") end
        if occursin("most_restricted", folder) push!(legend_filenames, L"A_p+\mathcal{G}  \, \, \textrm{rest.}") end
        if occursin("cross", folder) push!(legend_filenames, L"A_p  \, \, \textrm{rest.}") end
        if (occursin("horizontal", folder) || occursin("square", folder)) push!(legend_filenames, L"\mathcal{G}  \, \, \textrm{rest.}") end
        for file in readdir(joinpath(general_result_path, folder))
            if occursin("length_dict", file) && occursin("power_mult_$(power_mult)", file)
                dict = JSON.parsefile(joinpath(general_result_path, folder, file))
                service_est = []
                main_est = []
                service_true = []
                main_true = []
                for (b,br) in dict
                    if b ∈ service_branches
                        push!(service_est, br["length_est"])
                        push!(service_true, br["length_true"])
                    else
                        push!(main_est, br["length_est"])
                        push!(main_true, br["length_true"])
                    end
                end    
                pm = _SP.boxplot!(pm, [xposition],(main_true.-main_est)./main_true*100, yticks = ([i for i in -30:10:30], [L"%$i" for i in -30:10:30]),
                    legend=:bottomright, color="lightgrey", xticks = ((1:4),legend_filenames))
                pm = _SP.dotplot!(pm,  [xposition],(main_true.-main_est)./main_true*100, color="grey", legend=false)
                ps = _SP.boxplot!(ps,  [xposition],(service_true.-service_est)./service_true*100, yticks = ([i for i in -30:10:30], [L"%$i" for i in -30:10:30]),
                    legend=:bottomright,  color="lightgrey",xticks = ((1:4),legend_filenames))     
                ps = _SP.dotplot!(ps,  [xposition],(service_true.-service_est)./service_true*100, color="grey", legend=false)
                xposition+=1
            end
        end
    end
    _SP.savefig(pm, joinpath(general_result_path, "branch_length_main_case_$(case)_power_mult_$(power_mult).png"))
    _SP.savefig(pm, joinpath(general_result_path, "branch_length_main_case_$(case)_power_mult_$(power_mult).pdf"))
    _SP.savefig(ps, joinpath(general_result_path, "branch_length_service_case_$(case)_power_mult_$(power_mult).png"))
    _SP.savefig(ps, joinpath(general_result_path, "branch_length_service_case_$(case)_power_mult_$(power_mult).pdf"))
end

########################################################
#########                                      #########
#########              SHUNT PLOTS             #########
#########                                      #########
########################################################

function shunt_plots_scatter(general_result_path::String, case::String; power_mult::Float64=1.0)

    markershape_dict = Dict(
        "30l_ug_cross_only" => :circle,
        "30l_ug_most_restricted" => :diamond,
        "30l_ug_no_restriction"  => :rect,
        "30l_ug_squared_only"   =>  :utriangle,
        "30l_oh_cross_only" => :circle,
        "30l_oh_most_restricted" => :diamond,
        "30l_oh_no_restriction"  => :rect,
        "30l_oh_horizontal_only"   =>  :utriangle,
        "eulvtf_ug_cross_only" => :circle,
        "eulvtf_ug_most_restricted" => :diamond,
        "eulvtf_ug_no_restriction"  => :rect,
        "eulvtf_ug_squared_only"   =>  :utriangle,
        "eulvtf_oh_cross_only" => :circle,
        "eulvtf_oh_most_restricted" => :diamond,
        "eulvtf_oh_no_restriction"  => :rect,
        "eulvtf_oh_horizontal_only"   =>  :utriangle
    )

    markersize_dict = Dict(
        "30l_ug_cross_only" =>8,
        "30l_ug_most_restricted" =>4,
        "30l_ug_no_restriction"  =>5,
        "30l_ug_squared_only"   =>5,
        "30l_oh_cross_only" =>8,
        "30l_oh_most_restricted" =>4,
        "30l_oh_no_restriction"  =>5,
        "30l_oh_horizontal_only"   =>5,
        "eulvtf_ug_cross_only" =>8,
        "eulvtf_ug_most_restricted" =>4,
        "eulvtf_ug_no_restriction"  =>5,
        "eulvtf_ug_squared_only"   =>4,
        "eulvtf_oh_cross_only" =>8,
        "eulvtf_oh_most_restricted" =>4,
        "eulvtf_oh_no_restriction"  =>5,
        "eulvtf_oh_horizontal_only"   =>4
    )

    folders = [f for f in readdir(general_result_path) if occursin(case, f)]

    filename = case == "30l_oh" ? "30l_oh_most_restricted/_30loads_oh__shunts_scenario_1__power_mult_2.0__resistive_true.json" : "eulvtf_oh_most_restricted/_eulvtf_oh__shunts_scenario_1__power_mult_2.0__resistive_true.json"
    xtkz = case == "30l_oh" ? (5:5:25, [L"%$i" for i in 5:5:25]) : (5:5:57, [L"%$i" for i in 5:5:57])

    dict_truth = JSON.parsefile(joinpath(general_result_path, filename))

    shunt_true = []
    for (_, sh) in dict_truth
        push!(shunt_true, 1/sh["shunt_true"]["gs"]*57.68533333333333)
    end

    _SP.scatter(shunt_true, ylabel=L"1/y^{sh} \, \, [\Omega]", xlabel = L"\textrm{Shunt} \, \, \textrm{id.} \, \, \textrm{[-]}", label = L"\textrm{True}", 
                xticks = xtkz, xtickfontsize=13,ytickfontsize=12, ylabelfontsize=14, legendfontsize=12, markershape = :xcross, ms=6, color="black",
                yticks = ([20, 40, 60, 80, 100, 150, 250, 300], [L"%$i" for i in [20, 40, 60, 80, 100, 150, 250, 300]]))

    legend_filenames = []
    for folder in folders
        if isdir(joinpath(general_result_path, folder))
            for file in readdir(joinpath(general_result_path, folder))
                if occursin("__shunts", file) && occursin("power_mult_$(power_mult)", file)
                    if occursin("no_restriction", folder) push!(legend_filenames, L"\textrm{No  \, \, rest.}") end
                    if occursin("most_restricted", folder) push!(legend_filenames, L"A_p+\mathcal{G}  \, \, \textrm{rest.}") end
                    if occursin("cross", folder) push!(legend_filenames, L"A_p  \, \, \textrm{rest.}") end
                    if (occursin("horizontal", folder) || occursin("square", folder)) push!(legend_filenames, L"\mathcal{G}  \, \, \textrm{rest.}") end            
                    dict = JSON.parsefile(joinpath(general_result_path, folder, file))
                    shunt_est  = []
                    for (_, sh) in dict
                        push!(shunt_est, 1/sh["shunt_est"]["gs"]*57.68533333333333)
                    end
                    _SP.scatter!(shunt_est , color="grey", label=legend_filenames[end], markershape = markershape_dict["$folder"], ms = markersize_dict["$folder"], mc=:white)
                end
            end
        end
    end
    
    _SP.savefig(joinpath(general_result_path, "shunt_case_scatter_$(case)_power_mult_$(power_mult).png"))
    _SP.savefig(joinpath(general_result_path, "shunt_case_scatter_$(case)_power_mult_$(power_mult).pdf"))
end