using Graphs, SimpleWeightedGraphs

function carson_constants()
    c₁ = 3.28084e-3
    c₂ = 8.0252
    r_pq = 0.049348 # not a variable at all --> TODO explain in the paper and make a table of the number of variables
    return c₁, c₂, r_pq
end

function build_rx_sol_dict(mn_data::Dict, sol::Dict)

    if haskey(mn_data["nw"]["1"]["rho"], "rho_value")
        ρ = mn_data["nw"]["1"]["rho"]["rho_value"]
        α = mn_data["nw"]["1"]["alpha"]["alpha_value"]
        T = mn_data["nw"]["1"]["temperature"]["temperature_value"]
    end

    c₁, c₂, r_pq = carson_constants()

    for (lc_id, lc) in sol["solution"]["nw"]["1"]["linecode_map"]
        A_p = lc["A_p"] 
        gmr = lc["gmr"]
        dist = [k for k in keys(lc) if occursin("dij", k)]
        Dij = lc[dist[1]]

        r_init = fill(r_pq, size(A_p))
        if haskey(mn_data["nw"]["1"]["rho"], "rho_value")
            lc["r"] = r_init.+_LA.diagm([ρ/A*(1+(α*(T-20))) for A in A_p])
        elseif haskey(mn_data["nw"]["1"]["linecode_map"][parse(Int, lc_id)], "r_ac")
            lc["r"] = r_init.+_LA.diagm([r_ac for r_ac in mn_data["nw"]["1"]["linecode_map"][parse(Int, lc_id)]["r_ac"]])
        elseif haskey(mn_data["nw"]["1"]["linecode_map"][parse(Int, lc_id)], "r_material")
            lc["r"] = r_init.+_LA.diagm([mn_data["nw"]["1"]["linecode_map"][parse(Int, lc_id)]["r_material"][i]/A_p[i] for i in 1:length(A_p)])
        end
    
        x = zeros(size(A_p)[1], size(A_p)[1])
        exploit_horizontality = haskey(mn_data["nw"]["1"]["settings"], "exploit_horizontality") && mn_data["nw"]["1"]["settings"]["exploit_horizontality"]
        exploit_squaredness   = haskey(mn_data["nw"]["1"]["settings"], "exploit_squaredness")   && mn_data["nw"]["1"]["settings"]["exploit_squaredness"]

        for c in CartesianIndices(x)
            if c[1] == c[2]
                x[c] = 0.062832* ( c₂ + log(1) - log( c₁*gmr[c[1]] ) ) 
            elseif c[1] > c[2]
                # do nothing, symmetry exploited below
            else
                distance_indices = [(1,2), (1,3), (2,3), (1,4), (2,4), (3,4)]
                idx = findfirst(x->x == (c[1], c[2]), distance_indices)
                if !exploit_horizontality && !exploit_squaredness
                    x[c[1], c[2]] = x[c[2], c[1]] = 0.062832* ( c₂ + log(1) - log( c₁*Dij[idx] ) ) 
                elseif exploit_horizontality
                    if idx ∈ [1, 3, 6] 
                        dist_expr = Dij[Int(floor(idx/3))+1] # Int(floor(idx/3))+1 maps from distance indices 1,3,6 to the length three vector [(1,2), (2,3), (3,4)] which are the only distances in te horizantal case
                    elseif idx == 2 
                        dist_expr = Dij[1]+Dij[2]
                    elseif idx == 4
                        dist_expr = Dij[1]+Dij[2]+Dij[3]
                    else
                        dist_expr = Dij[2]+Dij[3]
                    end                        
                    x[c[1], c[2]] = x[c[2], c[1]] = 0.062832 * (c₂ + log(1) - log( c₁*dist_expr ))
                elseif exploit_squaredness
                    if idx ∈ [1, 2] 
                        dist_expr = Dij[1]
                    elseif idx ∈ [5, 6] 
                        dist_expr = Dij[4]
                    elseif idx == 3
                        dist_expr = Dij[3]
                    elseif idx == 4
                        dist_expr = Dij[2]
                    end
                    x[c[1], c[2]] = x[c[2], c[1]] = 0.062832 * (c₂ + log(1) - log( c₁*dist_expr ))
                end
            end
        end
        lc["x"] = x
    end

    for (b, branch) in sol["solution"]["nw"]["1"]["branch"]
        if length(split(b, ",")) == 1
            l = branch["l"]
            lc_name = mn_data["nw"]["1"]["branch"][b]["orig_linecode"]
            lc_id = [l for (l,lc) in mn_data["nw"]["1"]["linecode_map"] if lc["name"] == lc_name][1]
            branch["R"] = l*sol["solution"]["nw"]["1"]["linecode_map"]["$lc_id"]["r"]
            branch["X"] = l*sol["solution"]["nw"]["1"]["linecode_map"]["$lc_id"]["x"]
        end
    end

    for (_, nw) in sol["solution"]["nw"]
        for (b, bus) in nw["bus"]
            bus["vmn"] = sqrt.((bus["vr"][1:(end-1)].-bus["vr"][end]).^2+(bus["vi"][1:(end-1)].-bus["vi"][end]).^2)
        end
    end

    return sol
end

function reconstruct_powers(mn_data::Dict, sol::Dict) # TODO: superfluous function probably, remove once 100% sure about validation
    for (n, nw) in sol["solution"]["nw"]
        for (l,load) in nw["load"]
            load_bus = mn_data["nw"]["1"]["load"][l]["load_bus"]
            load_branch = [br for (br, branch) in mn_data["nw"][n]["branch"] if branch["t_bus"] == load_bus][1]
            arc = "($load_branch, $(mn_data["nw"][n]["branch"]["$load_branch"]["f_bus"]), $load_bus)"
            cr = nw["branch"][arc]["cr"]
            ci = nw["branch"][arc]["ci"]
            vr_p = nw["bus"]["$load_bus"]["vr"][1:(end-1)]
            vi_p = nw["bus"]["$load_bus"]["vi"][1:(end-1)]
            vr_n = nw["bus"]["$load_bus"]["vr"][end]
            vi_n = nw["bus"]["$load_bus"]["vi"][end]
            load["p_ime"] = (vr_p .- vr_n).*cr.+(vi_p.-vi_n).*ci 
            load["q_ime"] = -(vr_p.-vr_n).*ci.+(vi_p.-vi_n).*cr
        end
    end
    for (n, nw) in sol["solution"]["nw"]
        for (b,branch) in nw["branch"]
            if occursin("(", b)
                f_bus = parse(Int, split(b, ",")[2])
                t_bus = parse(Int, split(b, ",")[3][end-1])
                cr = branch["cr"]
                ci = branch["ci"]
                vr_p = nw["bus"]["$f_bus"]["vr"][1:(end-1)]
                vi_p = nw["bus"]["$f_bus"]["vi"][1:(end-1)]    
                vr_n = nw["bus"]["$f_bus"]["vr"][end]
                vi_n = nw["bus"]["$f_bus"]["vi"][end]    
                if length(vr_p) == 3 && length(cr) == 2 # this is a load arc
                    load = [load["index"] for (l, load) in mn_data["nw"]["1"]["load"] if load["load_bus"] ∈ [f_bus, t_bus]][1]
                    loadphase = mn_data["nw"]["1"]["load"]["$load"]["connections"][1]
                    branch["p"] =  (vr_p[loadphase] .- vr_n).*cr[1].+(vi_p.-vi_n).*ci[1]
                    branch["q"] = -(vr_p[loadphase] .- vr_n).*ci[1].+(vi_p.-vi_n).*cr[1]
                else
                    branch["p"] =  (vr_p .- vr_n).*cr[1:(end-1)].+(vi_p.-vi_n).*ci[1:(end-1)]
                    branch["q"] = -(vr_p .- vr_n).*ci[1:(end-1)].+(vi_p.-vi_n).*cr[1:(end-1)]
                end
            end
        end
    end
    return sol
end

function get_cumulative_impedance_of_loads_from_data(mn_data::Dict, multiply_z_pu::Bool)

    data = haskey(mn_data, "nw") ? mn_data["nw"]["1"] : mn_data
    loads = [load["index"] for (l,load) in data["load"]]

    map_dict = Dict{String, Any}()
    cnt = 1 
    for (b,bus) in data["bus"]
        map_dict[b] = cnt
        cnt+=1
    end

    transfo_bus = data["gen"]["1"]["gen_bus"]#[bus["index"] for (b, bus) in data["bus"] if bus["bus_type"]==3][1]
    graph_trfo_bus = map_dict["$transfo_bus"]

    load_dist_dict = Dict{String, Any}()

    for load in loads
        g_r = SimpleWeightedGraph(length(data["bus"]))
        g_x = SimpleWeightedGraph(length(data["bus"]))    
        load_dist_dict["$load"] = Dict{String, Any}()
        load_bus = data["load"]["$load"]["load_bus"]
        conn = data["load"]["$load"]["connections"][1] # ignore neutral at this point
        for (b, branch) in data["branch"]
            if conn ∈ branch["f_connections"]
                bfb = map_dict["$(branch["f_bus"])"]
                tfb = map_dict["$(branch["t_bus"])"]
                Graphs.add_edge!(g_r, bfb, tfb, branch["br_r"][1])                
                Graphs.add_edge!(g_x, bfb, tfb, branch["br_x"][1])
            end
        end
        dijkstra_r = Graphs.dijkstra_shortest_paths(g_r, graph_trfo_bus)
        dijkstra_x = Graphs.dijkstra_shortest_paths(g_x, graph_trfo_bus)
        pu = multiply_z_pu ? data["settings"]["z_pu"] : 1.
        load_dist_dict["$load"]["Rc_true"] = dijkstra_r.dists[map_dict["$load_bus"]]*pu
        load_dist_dict["$load"]["Xc_true"] = dijkstra_x.dists[map_dict["$load_bus"]]*pu
        load_dist_dict["$load"]["Zc_true"] = sqrt(load_dist_dict["$load"]["Rc_true"]^2+load_dist_dict["$load"]["Xc_true"]^2)
    end
    return load_dist_dict
end

function get_cumulative_impedance_of_loads_from_sol(mn_data::Dict, sol::Dict, divide_z_pu::Bool)

    data = haskey(mn_data, "nw") ? mn_data["nw"]["1"] : mn_data
    loads = [load["index"] for (l,load) in data["load"]]

    map_dict = Dict{String, Any}()
    cnt = 1 
    for (b,bus) in data["bus"]
        map_dict[b] = cnt
        cnt+=1
    end

    transfo_bus = data["gen"]["1"]["gen_bus"]#[bus["index"] for (b, bus) in data["bus"] if bus["bus_type"]==3][1]
    graph_trfo_bus = map_dict["$transfo_bus"]

    load_dist_dict = Dict{String, Any}()

    for load in loads
        g_r = SimpleWeightedGraph(length(data["bus"]))
        g_x = SimpleWeightedGraph(length(data["bus"]))    
        load_dist_dict["$load"] = Dict{String, Any}()
        load_bus = data["load"]["$load"]["load_bus"]
        conn = data["load"]["$load"]["connections"][1] # ignore neutral at this point
        for (b, branch) in data["branch"]
            if conn[1] ∈ branch["f_connections"]
                bfb = map_dict["$(branch["f_bus"])"]
                tfb = map_dict["$(branch["t_bus"])"]
                Graphs.add_edge!(g_r, bfb, tfb, sol["solution"]["nw"]["1"]["branch"][b]["R"][1])                
                Graphs.add_edge!(g_x, bfb, tfb, sol["solution"]["nw"]["1"]["branch"][b]["X"][1])
            end
        end
        dijkstra_r = Graphs.dijkstra_shortest_paths(g_r, graph_trfo_bus)
        dijkstra_x = Graphs.dijkstra_shortest_paths(g_x, graph_trfo_bus)
        pu = divide_z_pu ? 1/data["settings"]["z_pu"] : 1.
        load_dist_dict["$load"]["Rc_est"] = dijkstra_r.dists[map_dict["$load_bus"]]*pu
        load_dist_dict["$load"]["Xc_est"] = dijkstra_x.dists[map_dict["$load_bus"]]*pu
        load_dist_dict["$load"]["Zc_est"] = sqrt(load_dist_dict["$load"]["Rc_est"]^2+load_dist_dict["$load"]["Xc_est"]^2)
    end
    return load_dist_dict
end

function drop_results(case, result_path, other_string, summary_df, sol, mn_data, timestep_set, seed, add_meas_noise, power_mult, A_p_bounds, dist_bounds, r_ac_error, use_length_bounds, length_bounds_percval, imp_est, imp_true, real_volts, real_vas, est_volts, est_vas, exploit_equal_crossection, exploit_squareness, exploit_horizontality; save_summary::Bool=true)
    
    # drop impedances as JSONs
    imp_est  = JSON.json(imp_est)
    imp_true = JSON.json(imp_true)
    open("$(result_path)_$(case)_imp_est_scenario_$(seed)_$(other_string).json","w") do f 
        write(f, imp_est) 
    end
    open("$(result_path)_$(case)_imp_true_scenario_$(seed)_$(other_string).json","w") do f 
        write(f, imp_true) 
    end

    # drop voltages as CSVs
    real_volts |> CSV.write("$(result_path)_$(case)_real_volts_scenario_$(seed)_$(other_string).csv")
    est_volts |> CSV.write("$(result_path)_$(case)_est_volts_scenario_$(seed)_$(other_string).csv")

    real_vas |> CSV.write("$(result_path)_$(case)_real_vas_scenario_$(seed)_$(other_string).csv")
    est_vas |> CSV.write("$(result_path)_$(case)_est_vas_scenario_$(seed)_$(other_string).csv")

    # drop general result summary
    summary_df = build_summary_results_csv(case, result_path, other_string, summary_df, sol, timestep_set, seed, add_meas_noise, power_mult, A_p_bounds, dist_bounds, r_ac_error, use_length_bounds, length_bounds_percval, save_summary, exploit_equal_crossection, exploit_squareness, exploit_horizontality)
    if save_summary  summary_df |> CSV.write("$(result_path)_$(case)_general_summary_scenario_$(seed)_$(other_string).csv") end

    # linecode estimation result
    df_linecode, length_dict = build_linecode_results(sol, mn_data, seed)
    df_linecode |> CSV.write("$(result_path)_$(case)_linecode_results_scenario_$(seed)_$(other_string).csv")

    open("$(result_path)_$(case)_length_dict_scenario_$(seed)_$(other_string).json","w") do f 
        write(f, length_dict) 
    end

end

function drop_shunt_results(case, result_path, other_string, sol, mn_data, seed, shunt_resistive)

    shunt_dict = Dict{String, Any}()
    for (s, shunt) in mn_data["nw"]["1"]["shunt"]
        shunt_sol = sol["solution"]["nw"]["1"]["bus"]["$(shunt["shunt_bus"])"]

        if shunt_resistive 
            shunt_dict[s] = Dict("shunt_est" =>  Dict("gs" => shunt_sol["g_sh"]), 
                                 "shunt_true" => Dict("gs"=> shunt["gs"][end, end]))
        else
            shunt_dict[s] = Dict("shunt_est" => Dict("gs" => shunt_sol["g_sh"], 
                                                    "bs" => shunt_sol["b_sh"]),
                                "shunt_true" => Dict("gs"=> shunt["gs"][end, end], 
                                                    "bs"=> shunt["bs"][end, end]))
        end
    end
    shunt_dict = JSON.json(shunt_dict)
    open("$(result_path)_$(case)_shunts_scenario_$(seed)_$(other_string)_resistive_$(shunt_resistive).json","w") do f 
        write(f, shunt_dict) 
    end

end

function build_linecode_results(sol, mn_data, seed)
    df_linecode = _DF.DataFrame(fill([], 7), ["linecode_id", "linecode_name", "A_p_est", "dist_est", "x_est", "r_est", "scenario_id"])
    for (l, linecode) in sol["solution"]["nw"]["1"]["linecode_map"]
        dij = haskey(linecode, "dij_2w") ? "dij_2w" : "dij_4w" 
        push!(df_linecode, [parse(Int, l), mn_data["nw"]["1"]["linecode_map"][parse(Int, l)]["name"], linecode["A_p"], linecode[dij], linecode["x"], linecode["r"], seed])
    end
    length_dict = Dict{String, Any}()
    for (b, branch) in sol["solution"]["nw"]["1"]["branch"]
        length_dict[b] = Dict("length_est" => branch["l"]*1000, # ×1000 so they are both in meters
                              "length_true" => mn_data["nw"]["1"]["branch"][b]["orig_length"])
    end
    return df_linecode, JSON.json(length_dict)
end

function build_estimated_volts_dataframe(sol::Dict, mn_data::Dict, seed::Int)
    est_volts = _DF.DataFrame(fill([], length(mn_data["nw"]["1"]["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in mn_data["nw"]["1"]["load"]], ["scenario_id", "time_step"]))
    for (n, nw) in sol["solution"]["nw"]
        for (i, bus) in nw["bus"]
            bus["vm"] = sqrt.( (bus["vr"][1:(end-1)] .- bus["vr"][end]).^2 + (bus["vi"][1:(end-1)] .- bus["vi"][end]).^2 )
        end
        push!(est_volts, vcat([nw["bus"]["$(load["load_bus"])"]["vm"][1] for (l,load) in mn_data["nw"]["1"]["load"]], [seed, parse(Int, n)]))
    end
    return est_volts
end

function build_estimated_vas_dataframe(sol::Dict, mn_data::Dict, seed::Int)
    est_vas = _DF.DataFrame(fill([], length(mn_data["nw"]["1"]["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in mn_data["nw"]["1"]["load"]], ["scenario_id", "time_step"]))
    for (n, nw) in sol["solution"]["nw"]
        for (i, bus) in nw["bus"]
            bus["va"] = [atan(bus["vi"][i] / bus["vr"][i] ) for i in 1:(length(bus["vr"])-1)]
        end
        push!(est_vas, vcat([nw["bus"]["$(load["load_bus"])"]["va"][1] for (l,load) in mn_data["nw"]["1"]["load"]], [seed, parse(Int, n)]))
    end
    return est_vas
end

function build_summary_results_csv(case, result_path, other_string, summary_df, sol, timestep_set ,seed, add_meas_noise, power_mult, A_p_bounds, dist_bounds, r_ac_error, use_length_bounds, length_bounds_percval, save_as_csv, exploit_equal_crossection, exploit_squareness, exploit_horizontality)
    if isempty(summary_df)
        summary_df = _DF.DataFrame(fill([], 14), ["scenario_id", "length_bounds", "meas_noise", "power_mult", "A_p_bounds", "dist_bounds", "r_ac_error", "exploit_equal_crossection", "exploit_squareness", "exploit_horizontality", "solve_time", "solve_status", "objective", "timestep_set"])
    end
    lb = use_length_bounds ? length_bounds_percval : false
    push!(summary_df, [seed, lb, add_meas_noise, power_mult, A_p_bounds, dist_bounds, r_ac_error, exploit_equal_crossection, exploit_squareness, exploit_horizontality, sol["solve_time"], sol["termination_status"], sol["objective"], timestep_set])
    if save_as_csv
        summary_df |> CSV.write("$(result_path)_$(case)_general_summary_scenario_$(seed)_$(other_string).csv")
    end
    return summary_df
end