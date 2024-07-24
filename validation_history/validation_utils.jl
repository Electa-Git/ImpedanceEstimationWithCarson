import CSV
import DataFrames as _DF
import LinearAlgebra: diagm
import Distributions as _DST

function add_zib_voltage!(mn_data::Dict)
    for (n,nw) in mn_data["nw"]
        m = maximum(parse.(Int, keys(nw["meas"])))
        nw["meas"]["$(m+1)"] = Dict{String, Any}(
            "dst" => [_DST.Normal(μ, 0.00319208) for μ in nw["bus"]["1"]["pf_vm"]],
            "var" => :vm,
            "cmp" => :bus,
            "cmp_id" => 1
        )
    end
    return mn_data    
end

function get_clean_validation_case(eng, profiles_df)

    data = _PMD.transform_data_model(eng, kron_reduce=false, phase_project=false)
    z_pu = (data["settings"]["voltage_scale_factor"]*data["settings"]["vbases_default"]["3"])^2/(data["settings"]["power_scale_factor"])

    load_buses = [load["load_bus"] for (_, load) in data["load"]]

    _IMP.clean_4w_data!(data, profiles_df, merge_buses_diff_linecodes = false, eng = eng)
    _PMD.add_start_vrvi!(data)
    
    _IMP.make_loadbuses_loadbranches_singlephase!(data) # TODO: or change .dss file?
    
    for (_,bus) in data["bus"]
        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer") && bus["index"] ∉ load_buses 
            bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
            bus["grounded"] .=  0
            bus["imp_grounded"] = fill(false, length(bus["terminals"]))
        elseif bus["index"] ∈ load_buses 
            bus["imp_grounded"] = [false, true]  # IMPORTANT! add an imp ground to all users' neutrals
            bus["grounded"] = [false, false]
            bus["vm_pair_lb"] = [(bus["terminals"][1], bus["terminals"][2], 0.9)]
            bus["vm_pair_ub"] = [(bus["terminals"][1], bus["terminals"][2], 1.1)]
            bus["vr_start"] = bus["vr_start"][[bus["terminals"][1], bus["terminals"][2]]]
            bus["vi_start"] = bus["vi_start"][[bus["terminals"][1], bus["terminals"][2]]]      
        elseif bus["bus_type"] == 3
            bus["imp_grounded"] = fill(false, length(bus["terminals"]))
        end
    end
    
    for (_, branch) in data["branch"]
        if length(branch["t_connections"]) == 2
            i = branch["t_connections"][1]
            j = branch["t_connections"][2]
            for key in ["rate_a", "rate_b", "rate_c", "c_rating_a", "c_rating_b", "c_rating_c", "angmin", "angmax"]
                branch[key] = branch[key][branch["t_connections"]]
            end
            for key in ["g_fr", "g_to", "b_fr", "b_to"]
                branch[key] = [branch[key][i,i] branch[key][i,j]; branch[key][j,i] branch[key][j,j]]
            end
            r = eng["linecode"]["lc6"]["rs"]
            x = eng["linecode"]["lc6"]["xs"]
            l = eng["line"][branch["name"]]["length"]
            branch["br_r"] = r.*l./z_pu
            branch["br_x"] = x.*l./z_pu
        end
    end
    return data
end

function set_up_validation_series(;power_mult::Float64 = 25., t_start::Int=1, t_end::Int=16)

    profiles_df = CSV.read(joinpath(_IMP.DATA_DIR, "profiles.csv"), _DF.DataFrame,  ntasks = 1)

    ntw_path = joinpath(_IMP.DATA_DIR, "opendss/small_validation_case/Master.dss")
    eng = _PMD.parse_file(ntw_path, transformations=[_PMD.transform_loops!,_PMD.remove_all_bounds!]) # note that after this the impedances in linecodes are in Ω/m
    _IMP.rm_enwl_transformer!(eng)
    _IMP.reduce_enwl_lines_eng!(eng)
    eng["settings"]["sbase_default"] = 1
    
    data = get_clean_validation_case(eng, profiles_df)
    _IMP.add_length!(data, eng)
    
    # I need to do the below because could not add a linecode with 2 wires directly in .dss, the PMD .dss parser didn't like it
    for (l, line) in eng["line"]
        if l != "line8"
            line["linecode"] = "lc6"
        end
    end
    
    ######## SET MEASUREMENT SPECS
    
    max_volt_error = 2.3 # in Volts\
    max_power_error = 0.1 # kW
    original_sourcebus_id = collect(keys(data["settings"]["vbases_default"]))[1]
    
    σ_v = 1/3*max_volt_error/(data["settings"]["vbases_default"][original_sourcebus_id]*data["settings"]["voltage_scale_factor"])
    σ_d = 1/3*max_power_error/(data["settings"]["sbase_default"])
    σ_g = 1/3*max_power_error/(data["settings"]["sbase_default"])
    
    ############### CREATE MULTINETWORK DATA WITH MEASUREMENT TIMESERIES ###############
    # it runs a powerflow for each time step first, so it takes some time...
    
    pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "print_level"=>0 )
    mn_data = _IMP.build_multinetwork_dsse_data(data, profiles_df, pf_solver, σ_v, σ_d, σ_g; t_start=t_start, t_end=t_end, add_noise=false, power_mult = power_mult)
    
    mn_data["nw"]["1"]["linecode_map"] = Dict{Int, Any}() 
    for (id, code) in enumerate(keys(eng["linecode"]))
        if code == "lc6"
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 2
            )
        else
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 4
            )
        end
    end
    
    make_all_branches_untrustworthy!(mn_data, eng)
    
    z_pu = (data["settings"]["voltage_scale_factor"]*data["settings"]["vbases_default"]["3"])^2/(data["settings"]["power_scale_factor"])
    T = 55
    α = .00429 
    ρ = 10 
    μ = 1
    
    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["rescaler"] = 1.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = μ
    mn_data["nw"]["1"]["temperature"] = Dict("temperature_value" => T)
    mn_data["nw"]["1"]["rho"] = Dict("rho_value" => ρ)
    mn_data["nw"]["1"]["alpha"] = Dict("alpha_value" => α)
    
    ############### SOLVE IMPEDANCE ESTIMATION ###############
    
    # ignore shunt impedances (they are zero in the power flow anyway)
    for (b,bus) in mn_data["nw"]["1"]["bus"]
        bus["imp_grounded"] = fill(false, length(bus["terminals"]))
    end    
    return mn_data
end

function set_up_validation_shunt(;power_mult::Float64 = 25.)

    profiles_df = CSV.read(joinpath(_IMP.DATA_DIR, "profiles.csv"), _DF.DataFrame,  ntasks = 1)

    ntw_path = raw"data/Master.dss"
    eng = _PMD.parse_file(ntw_path, transformations=[_PMD.transform_loops!,_PMD.remove_all_bounds!]) # note that after this the impedances in linecodes are in Ω/m
    _IMP.rm_enwl_transformer!(eng)
    _IMP.reduce_enwl_lines_eng!(eng)
    eng["settings"]["sbase_default"] = 1
    
    data = get_clean_validation_case(eng, profiles_df)
    _IMP.add_length!(data, eng)
    
    # I need to do the below because could not add a linecode with 2 wires directly in .dss, the PMD .dss parser didn't like it
    for (l, line) in eng["line"]
        if l != "line8"
            line["linecode"] = "lc6"
        end
    end
    
    ######## SET MEASUREMENT SPECS
    
    max_volt_error = 2.3 # in Volts\
    max_power_error = 0.1 # kW
    original_sourcebus_id = collect(keys(data["settings"]["vbases_default"]))[1]
    
    σ_v = 1/3*max_volt_error/(data["settings"]["vbases_default"][original_sourcebus_id]*data["settings"]["voltage_scale_factor"])
    σ_d = 1/3*max_power_error/(data["settings"]["sbase_default"])
    σ_g = 1/3*max_power_error/(data["settings"]["sbase_default"])
    
    ############### CREATE MULTINETWORK DATA WITH MEASUREMENT TIMESERIES ###############
    # it runs a powerflow for each time step first, so it takes some time...
    
    pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "print_level"=>0 )
    mn_data, real_volts, real_vas = _IMP.build_multinetwork_dsse_data_with_shunts(data, profiles_df, pf_solver, add_noise=false, power_mult = power_mult)
        
    z_pu = (data["settings"]["voltage_scale_factor"]*data["settings"]["vbases_default"]["3"])^2/(data["settings"]["power_scale_factor"])
    T = 55
    α = .00429 
    ρ = 10 
    μ = 1
    
    mn_data["nw"]["1"]["settings"]["z_pu"] = z_pu
    mn_data["nw"]["1"]["settings"]["rescaler"] = 1.
    mn_data["nw"]["1"]["settings"]["mu_rel"] = μ
    mn_data["nw"]["1"]["temperature"] = Dict("temperature_value" => T)
    mn_data["nw"]["1"]["rho"] = Dict("rho_value" => ρ)
    mn_data["nw"]["1"]["alpha"] = Dict("alpha_value" => α)

    return mn_data
end

function make_all_branches_untrustworthy!(mn_data, eng)
    for (k, val) in mn_data["nw"]["1"]["linecode_map"]
        for (b, branch) in mn_data["nw"]["1"]["branch"] # all branch info useful for the model is stored in nw 1. Better than copying it to all nws, but there might be better ways still..
            branch["untrustworthy_branch"] = true # if not untrustworthy, we could fix the impedance. Default assumption is all are untrustworthy
            if eng["line"][branch["name"]]["linecode"] == val["name"]
                branch["orig_linecode"] = val["name"]
                branch["linecode_id"] = k
            end
        end
    end
    return mn_data
end

function make_all_branches_trustworthy!(mn_data, eng)
    for (k, val) in mn_data["nw"]["1"]["linecode_map"]
        for (b, branch) in mn_data["nw"]["1"]["branch"] 
            branch["untrustworthy_branch"] = false 
            if eng["line"][branch["name"]]["linecode"] == val["name"]
                branch["orig_linecode"] = val["name"]
                branch["linecode_id"] = k
            end
        end
    end
    return mn_data
end

function carson_constants()
    c₁ = 3.28084e-3
    c₂ = 8.0252
    r_pq = 0.049348 
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


# the function below has been used to generate the impedance matrices!
# A_p must be in mm²
# x, y must be in mm
# ρ must be in 
# α must be in
# T must be in degrees centigrade
# the results are in Ω/km
function forward_carson_linecode_generation(A_p, x, y, ρ, α, T, μ)
    # constants
    c₁ = 3.28084e-3
    c₂ = 8.0252
    r_pq = 0.049348
    distance_indices = [(1,2), (1,3), (2,3), (1,4), (2,4), (3,4)]

    # dependent variables
    gmr = exp(-μ/4)*[sqrt(A/π) for A in A_p] # is in mm

    if length(A_p) == 2
        Dij = sqrt((x[2]-x[1])^2+(y[2]-y[1])^2) #all lengths will be in mm
    elseif length(A_p) == 3
        Dij = [sqrt((x[d[2]]-x[d[1]])^2+(y[d[2]]-y[d[1]])^2) for d in distance_indices[1:3]]
    elseif length(A_p) == 4
        Dij = [sqrt((x[d[2]]-x[d[1]])^2+(y[d[2]]-y[d[1]])^2) for d in distance_indices]
    end

    # resistance
    r = fill(r_pq, length(A_p), length(A_p))+diagm([ρ/A*(1+(α*(T-20))) for A in A_p]) # this is i

    # reactance
    x = zeros(length(A_p), length(A_p))
    for c in CartesianIndices(x)
        if c[1] == c[2]
            x[c] = 0.062832 * ( c₂ + log(1/ ( c₁*gmr[c[1]] ) ) ) 
        elseif c[1] > c[2]
            # do nothing, symmetry exploited below
        else
            
            idx = findfirst(x->x == (c[1], c[2]), distance_indices)
            x[c[1], c[2]] = x[c[2], c[1]] = 0.062832 * ( c₂ + log(1 / ( c₁*Dij[idx] ) ) ) 
        end
    end
    return r, x
end
