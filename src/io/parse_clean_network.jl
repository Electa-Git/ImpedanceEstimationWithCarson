"""
Cleans a MATHEMATICAL data dictionary
"""
function clean_4w_data!(ntw::Dict, profiles_df::_DF.DataFrame; eng::Dict=Dict{String, Any}(), merge_buses_diff_linecodes::Bool = false)
    if merge_buses_diff_linecodes 
        remove_all_superfluous_buses!(ntw) 
        add_length!(ntw, eng)
        #TODO add function to add length even when we remove all superfluous buses
    end
    assign_load_to_parquet_id!(ntw, profiles_df)
    #### below removes virtual voltage source (but the transformer I removed in the eng model, thus beforehand)
    vsource_branch, vsource_bus, new_slackbus = find_voltage_source_branch_bus(ntw) # ← imported from src/core/benchmarking.jl
    ntw["gen"]["1"]["gen_bus"] = new_slackbus
    ntw["bus"]["$new_slackbus"] = deepcopy(ntw["bus"]["$vsource_bus"])
    ntw["bus"]["$new_slackbus"]["bus_i"] = new_slackbus
    ntw["bus"]["$new_slackbus"]["index"] = new_slackbus
    delete!(ntw["branch"], vsource_branch)
    delete!(ntw["bus"], "$vsource_bus")
    return ntw
end
"""
Goes get the length of the branches, which are stored in the engineering model
and adds them to the mathematical model.
"""
function add_length!(ntw::Dict, eng::Dict)
    for (_, br) in ntw["branch"]
        br["orig_length"] = eng["line"][br["name"]]["length"]
    end
    return ntw
end
"""
simply adds the degree (as defined in graph theory) of all nodes/buses in the feeder
"""
function add_degree_to_bus!(data)
    for (b, bus) in data["bus"]
        bus["degree"] = 0
        for (br, branch) in data["branch"]
            if branch["t_bus"] == bus["index"] || branch["f_bus"] == bus["index"]
                bus["degree"] += 1
            end
        end
    end
end
"""
removes all the buses that have degree 2 and are not connected to anything (they are electrically meaningless but the GIS does not know).
the impedance of the resulting segment is the sum of the impedances of all removed segments
"""
function remove_all_superfluous_buses!(data::Dict)
    @assert !haskey(data, "nw") "Please use `remove_all_intermediate_buses_mn` for multinetwork data dicts like this one"
    load_buses = ["$(load["load_bus"])" for (_, load) in data["load"]]
    gen_buses = ["$(gen["gen_bus"])" for (_, gen) in data["gen"]]
    add_degree_to_bus!(data)
    for lb in load_buses @assert data["bus"][lb]["degree"] == 1 "Load $lb is on the main cable, add a small connection cable with the appropriate util function!" end
    to_delete = [b for (b, bus) in data["bus"] if (b ∉ union!(gen_buses, load_buses) && bus["degree"] <= 2)]
    for db in to_delete
        data["bus"][db]["adjacent_buses"] = []
        data["bus"][db]["inout_branches"] = []
        for (br, branch) in data["branch"]
            if branch["f_bus"] == parse(Int, db) || branch["t_bus"] == parse(Int, db)
                push!(data["bus"][db]["inout_branches"], br)
                if branch["f_bus"] != parse(Int, db)
                    push!(data["bus"][db]["adjacent_buses"], "$(branch["f_bus"])")
                else
                    push!(data["bus"][db]["adjacent_buses"], "$(branch["t_bus"])")
                end
            end
        end
    end
    while !isempty(to_delete)
        for db in to_delete
            if any([b ∈ to_delete for b in data["bus"][db]["adjacent_buses"]])
                deletable_adj_bus = [b for b in data["bus"][db]["adjacent_buses"] if b ∈ to_delete][1]
                other_adj_bus = [b for b in data["bus"][db]["adjacent_buses"] if b != deletable_adj_bus][1]
                deletable_adj_bus_branches = data["bus"][deletable_adj_bus]["inout_branches"]
                delete_branch = first(intersect(Set(data["bus"][db]["inout_branches"]), Set(deletable_adj_bus_branches)))
                preserve_branch = [br for br in data["bus"][db]["inout_branches"] if br != delete_branch][1]
                
                Req = (data["branch"][preserve_branch]["br_r"] .+ data["branch"][delete_branch]["br_r"])
                Xeq = (data["branch"][preserve_branch]["br_x"] .+ data["branch"][delete_branch]["br_x"]) 
                data["branch"][preserve_branch]["br_r"] = Req
                data["branch"][preserve_branch]["br_x"] = Xeq
                
                data["branch"][preserve_branch]["f_bus"] = parse(Int64, other_adj_bus)
                data["branch"][preserve_branch]["t_bus"] = parse(Int64, deletable_adj_bus)
                data["bus"][deletable_adj_bus]["adjacent_buses"] = filter(x->x!=db, data["bus"][deletable_adj_bus]["adjacent_buses"])
                push!(data["bus"][deletable_adj_bus]["adjacent_buses"], other_adj_bus)
                data["bus"][deletable_adj_bus]["inout_branches"] = filter(x->x!=delete_branch, data["bus"][deletable_adj_bus]["inout_branches"])
                push!(data["bus"][deletable_adj_bus]["inout_branches"], preserve_branch)
                delete!(data["branch"], delete_branch)
            else
                delete_branch = data["bus"][db]["inout_branches"][1]
                preserve_branch = [br for br in data["bus"][db]["inout_branches"] if br != delete_branch][1]
                Req = (data["branch"][preserve_branch]["br_r"] .+ data["branch"][delete_branch]["br_r"])
                Xeq = (data["branch"][preserve_branch]["br_x"] .+ data["branch"][delete_branch]["br_x"]) 
                data["branch"][preserve_branch]["br_r"] = Req
                data["branch"][preserve_branch]["br_x"] = Xeq
                
                delete!(data["branch"], delete_branch)
                data["branch"][data["bus"][db]["inout_branches"][2]]["f_bus"] = parse(Int64, data["bus"][db]["adjacent_buses"][1])
                data["branch"][data["bus"][db]["inout_branches"][2]]["t_bus"] = parse(Int64, data["bus"][db]["adjacent_buses"][2])
            end
            delete!(data["bus"], db)
            to_delete = filter(x->x!=db, to_delete)
        end
    end
    # the lines below make sure that the orientation of the branch at the slack bus is from slack_bus to --> rest of feeder
    ref_bus = [bus["index"] for (_,bus) in data["bus"] if bus["bus_type"] == 3][1]
    ref_branch_fr = [b for (b, br) in data["branch"] if br["f_bus"] == ref_bus]
    if isempty(ref_branch_fr) 
        ref_branch_to = [b for (b, br) in data["branch"] if br["t_bus"] == ref_bus][1]
        f_bus = data["branch"][ref_branch_to]["f_bus"]
        data["branch"][ref_branch_to]["f_bus"] = ref_bus
        data["branch"][ref_branch_to]["t_bus"] = f_bus
    end
    return data
end
"""
Assign laods in the data dict to some nrel profile
"""
function assign_load_to_parquet_id!(data::Dict, df::_DF.DataFrame)
    parquet_ids = names(df)
    count = 1
    for (_, load) in data["load"]
        load["parquet_id"] = split(parquet_ids[count], "_")[end]
        count+=1
    end
    return data
end
""" 
timestep by timestep, reads the PQ profiles and adds the demand at each user
"""
function insert_profiles!(data, df, timestep)
    power_unit = data["settings"]["sbase"]
    @assert power_unit == 1e5 "The profiles are in kW, but the power_unit seems different. Please fix."
    for (_, load) in data["load"]
        p_id = load["parquet_id"]
        load["pd"] = [df[timestep, "P_kW_"*p_id]/power_unit]
        load["qd"] = [df[timestep, "Q_kVAr_"*p_id]/power_unit]
    end
end

function make_loadbuses_loadbranches_singlephase!(data::Dict)
    @assert !haskey(data, "nw") "Please use `make_loadbuses_loadbranches_singlephase_mn` for multinetwork data dicts like this one"
    load_buses = [load["load_bus"] for (l, load) in data["load"]]
    load_connections = Dict("$(load["load_bus"])" => load["connections"] for (l, load) in data["load"])

    load_branches_fr = Dict("$b" => "$(branch["f_bus"])" for (b, branch) in data["branch"] if branch["f_bus"] ∈ load_buses)
    load_branches_to = Dict("$b" => "$(branch["t_bus"])" for (b, branch) in data["branch"] if branch["t_bus"] ∈ load_buses)    

    load_branches = merge(load_branches_fr, load_branches_to)

    for (b, bus) in data["bus"]
        if bus["index"] ∈ load_buses 
            bus["terminals"] = load_connections[b]
            bus["vmin"] = bus["vmin"][load_connections[b]]
            bus["vmax"] = bus["vmax"][load_connections[b]]
            if haskey(bus, "pf_va") bus["pf_va"] = bus["pf_va"][load_connections[b]] end
            if haskey(bus, "pf_vm") bus["pf_vm"] = bus["pf_vm"][load_connections[b]] end
            bus["grounded"] = bus["grounded"][load_connections[b]]
        end
    end

    for (branch, loadbus) in load_branches
        terms = data["bus"][loadbus]["terminals"]
        data["branch"][branch]["f_connections"] = terms
        data["branch"][branch]["t_connections"] = terms
        if length(terms) == 1
            data["branch"][branch]["g_fr"] = data["branch"][branch]["g_to"] = data["branch"][branch]["b_fr"] = data["branch"][branch]["b_to"] = zeros(1,1)
            data["branch"][branch]["br_r"] = ones(1,1)*data["branch"][branch]["br_r"][terms[1], terms[1]]
            data["branch"][branch]["br_x"] = ones(1,1)*data["branch"][branch]["br_x"][terms[1], terms[1]]
        end
    end
end

function make_loadbuses_loadbranches_singlephase_validationcase!(data::Dict)
    @assert !haskey(data, "nw") "Please use `make_loadbuses_loadbranches_singlephase_mn` for multinetwork data dicts like this one"
    load_buses = [load["load_bus"] for (l, load) in data["load"]]
    load_connections = Dict("$(load["load_bus"])" => load["connections"] for (l, load) in data["load"])

    load_branches_fr = Dict("$b" => "$(branch["f_bus"])" for (b, branch) in data["branch"] if branch["f_bus"] ∈ load_buses)
    load_branches_to = Dict("$b" => "$(branch["t_bus"])" for (b, branch) in data["branch"] if branch["t_bus"] ∈ load_buses)    

    load_branches = merge(load_branches_fr, load_branches_to)

    for (b, bus) in data["bus"]
        if bus["index"] ∈ load_buses 
            bus["terminals"] = load_connections[b]
            if haskey(bus, "pf_va") bus["pf_va"] = bus["pf_va"][load_connections[b]] end
            if haskey(bus, "pf_vm") bus["pf_vm"] = bus["pf_vm"][load_connections[b]] end
        end
    end

    for (branch, loadbus) in load_branches
        terms = data["bus"][loadbus]["terminals"]
        data["branch"][branch]["f_connections"] = terms
        data["branch"][branch]["t_connections"] = terms
        if length(terms) == 1
            data["branch"][branch]["g_fr"] = data["branch"][branch]["g_to"] = data["branch"][branch]["b_fr"] = data["branch"][branch]["b_to"] = zeros(1,1)
            data["branch"][branch]["br_r"] = ones(1,1)*data["branch"][branch]["br_r"][terms[1], terms[1]]
            data["branch"][branch]["br_x"] = ones(1,1)*data["branch"][branch]["br_x"][terms[1], terms[1]]
        end
    end
end