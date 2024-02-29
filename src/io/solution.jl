using Graphs, SimpleWeightedGraphs

function carson_constants()
    c₁ = 3.28084e-3
    c₂ = 8.0252
    r_pq = 0.049348 # not a variable at all --> TODO explain in the paper and make a table of the number of variables
    return c₁, c₂, r_pq
end

function build_rx_sol_dict(mn_data::Dict, sol::Dict)
    ρ = mn_data["nw"]["1"]["rho"]["rho_value"]
    α = mn_data["nw"]["1"]["alpha"]["alpha_value"]
    T = mn_data["nw"]["1"]["temperature"]["temperature_value"]

    c₁, c₂, r_pq = carson_constants()

    for (_, lc) in sol["solution"]["nw"]["1"]["linecode_map"]
        A_p = lc["A_p"] 
        gmr = lc["gmr"]
        dist = [k for k in keys(lc) if occursin("dij", k)]
        Dij = lc[dist[1]]

        r_init = fill(r_pq, size(A_p))
        lc["r"] = r_init.+diagm([ρ/A*(1+(α*(T-20))) for A in A_p])
        lc["r"]
    
        x = zeros(size(A_p)[1], size(A_p)[1])
        for c in CartesianIndices(x)
            if c[1] == c[2]
                x[c] = 0.062832* ( c₂ + log(1) - log( c₁*gmr[c[1]] ) ) 
            elseif c[1] > c[2]
                # do nothing, symmetry exploited below
            else
                distance_indices = [(1,2), (1,3), (2,3), (1,4), (2,4), (3,4)]
                idx = findfirst(x->x == (c[1], c[2]), distance_indices)
                x[c[1], c[2]] = x[c[2], c[1]] = 0.062832* ( c₂ + log(1) - log( c₁*Dij[idx] ) ) 
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