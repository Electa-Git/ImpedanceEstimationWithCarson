import StatsPlots as _SP
import Format

"""
plots a boxplot for every 
"""
function cumulative_zrx_per_user(true_impedance_dict::Dict, est_impedance_path::String, case::String; whatt::String="Zc", ytickz::Vector=[])
    
    est_dicts = [joinpath(est_impedance_path, file) for file in readdir(est_impedance_path) if occursin("imp_est", file) && occursin(case, file)]

    nr_scenarios = unique([split(ed, "scenario_")[2][1] for ed in est_dicts])
    other_identifiers = unique([split(ed, "__")[end][1:end-5] for ed in est_dicts])

    idx = findfirst(x->occursin("scenario_$(minimum(nr_scenarios))", x), est_dicts)
    starter_est = est_dicts[idx]
    starter_label = split(starter_est, "__")[end][1:end-5]

    xticks = (1:length(true_impedance_dict), [parse(Int, k) for (k,v) in true_impedance_dict])

    p = _SP.scatter([(v["$(whatt)_true"]-JSON.parsefile(joinpath(est_impedance_path, starter_est))[k]["Zc_est"])/v["$(whatt)_true"]*100 for (k,v) in true_impedance_dict],
                    ylabel = "($(whatt)ᵗʳᵘᵉ - $(whatt)ᵉˢᵗ)/$(whatt)ᵗʳᵘᵉ × 100 [%]", labels = "sc. $(minimum(nr_scenarios)), id. $(starter_label)",
                    xlabel = "User id. [-]", xticks = xticks)

    for ed in est_dicts
        if ed != starter_est
            sc = split(ed, "scenario_")[2][1]
            id = split(ed, "__")[end][1:end-5]
           _SP.scatter!([(v["$(whatt)_true"]-JSON.parsefile(joinpath(est_impedance_path, ed))[k]["Zc_est"])/v["$(whatt)_true"]*100 for (k,v) in true_impedance_dict],
                        labels = "sc. $sc, id. $id")
        end
    end

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
