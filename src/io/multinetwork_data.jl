function build_multinetwork_dsse_data(data::Dict, df::_DF.DataFrame, pf_solver, σ_v::Float64, σ_d::Float64, σ_g::Float64; t_start::Int=1, t_end::Int=20, add_noise::Bool=false, seed::Int64=2, power_mult::Float64=1.0)

    mn_data = Dict{String, Any}()
    mn_data["multinetwork"] = true
    mn_data["nw"] = Dict{String, Any}()


    for key in ["name", "map", "se_settings", "data_model", "bus_lookup", "per_unit"]
        if haskey(data, key) mn_data[key] = data[key] end
    end

    real_volts = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["scenario_id", "time_step"]))

    power_unit = data["settings"]["sbase"]
    @info "The power unit in use is $power_unit and the load multiplication factor is $power_mult"

    for (ts_id, ts) in enumerate(t_start:t_end)
        
        # add info / initialize multinetwork dict
        mn_data["nw"]["$ts_id"] = Dict{String, Any}()
        mn_data["nw"]["$ts_id"]["original_timestep"] = ts
        mn_data["nw"]["$ts_id"]["per_unit"] = true

        # assign unbalance voltage magnitude at ref_bus # TODO do we want sth like this??
        # ref_bus = [b for (b,bus) in data["bus"] if bus["bus_type"] ==3][1]
        # data["bus"][ref_bus]["vm"] = vcat([randn(_RAN.MersenneTwister(ts)), randn(_RAN.MersenneTwister(ts+1000)), randn(_RAN.MersenneTwister(ts+2000))]./800 .+1 , 0)

        _build_dictionary_entries(data, mn_data, ts_id)
        insert_profiles!(data, df, ts, power_mult=power_mult) # inserts NREL load profiles for powerflow

        # Solves the opf, i.e., which gives noiseless input
        pf_results = _PMD.solve_mc_opf(data, _PMD.IVRENPowerModel, pf_solver)
        
        if string(pf_results["termination_status"]) ∉ ["LOCALLY_SOLVED"]
            @warn "power flow result for time step $ts_id converged to: $(pf_results["termination_status"])"
        end
        # converts vr and vi to vm (phase to neutral)
        pf_solution_to_voltage_magnitudes!(pf_results) 

        push!(real_volts, vcat([pf_results["solution"]["bus"]["$(load["load_bus"])"]["vm"][1] for (l,load) in data["load"]], [seed, ts]))

        # adds the powerflow results (P, Q, |U|) of this timestep on the mn dict, for future reference/comparison
        add_pf_result_to_mn_data!(mn_data["nw"]["$ts_id"], pf_results)

        # converts the powerflow results into (noisy or not) measurements
        add_measurements!(data, pf_results, σ_v, σ_d, σ_g, add_noise = add_noise, seed = seed, include_transfo_meas = false)

        # store this timestep in multinetwork dict
        mn_data["nw"]["$ts_id"]["meas"] = deepcopy(data["meas"]);

    end
    return mn_data, real_volts
end

function build_multinetwork_dsse_data_with_shunts(data::Dict, df::_DF.DataFrame, pf_solver, σ_v::Float64, σ_d::Float64, σ_g::Float64; t_start::Int=1, t_end::Int=20, add_noise::Bool=false, loads_with_shunts::Vector{String} = ["1"], gs::Vector{Float64} = [50.], bs::Vector{Float64} = [15.], seed::Int64=2, power_mult::Float64=1.0)

    mn_data = Dict{String, Any}()
    mn_data["multinetwork"] = true
    mn_data["nw"] = Dict{String, Any}()

    for key in ["name", "map", "se_settings", "data_model", "bus_lookup", "per_unit"]
        if haskey(data, key) mn_data[key] = data[key] end
    end

    real_volts = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["scenario_id", "time_step"]))

    count_shunts = 1
    for (l, load) in data["load"]
        if l ∈ loads_with_shunts
            data["shunt"][l] = Dict{String, Any}(
                "shunt_bus" => load["load_bus"],
                "connections" => load["connections"],
                "status" => 1,
                "dispatchable" => 0,
                "gs" => [0. 0.; 0. gs[count_shunts]] ,
                "bs" => [ 0.  0.; 0. bs[count_shunts]]
            )
            count_shunts+=1
        end
    end

    power_unit = data["settings"]["sbase"]
    @info "The power unit in use is $power_unit and the load multiplication factor is $power_mult"

    for (ts_id, ts) in enumerate(t_start:t_end)
        
        # add info / initialize multinetwork dict
        mn_data["nw"]["$ts_id"] = Dict{String, Any}()
        mn_data["nw"]["$ts_id"]["original_timestep"] = ts
        mn_data["nw"]["$ts_id"]["per_unit"] = true

        # assign unbalance voltage magnitude at ref_bus # TODO do we want sth like this??
        # ref_bus = [b for (b,bus) in data["bus"] if bus["bus_type"] ==3][1]
        # data["bus"][ref_bus]["vm"] = vcat([randn(_RAN.MersenneTwister(ts)), randn(_RAN.MersenneTwister(ts+1000)), randn(_RAN.MersenneTwister(ts+2000))]./800 .+1 , 0)

        _build_dictionary_entries(data, mn_data, ts_id)
        insert_profiles!(data, df, ts, power_mult=power_mult) # inserts NREL load profiles for powerflow

        # Solves the opf, i.e., which gives noiseless input
        pf_results = _PMD.solve_mc_opf(data, _PMD.IVRENPowerModel, pf_solver)
        
        if string(pf_results["termination_status"]) ∉ ["LOCALLY_SOLVED"]
            @warn "power flow result for time step $ts_id converged to: $(pf_results["termination_status"])"
        end
        # converts vr and vi to vm (phase to neutral)
        pf_solution_to_voltage_magnitudes!(pf_results)
        
        push!(real_volts, vcat([pf_results["solution"]["bus"]["$(load["load_bus"])"]["vm"][1] for (l,load) in data["load"]], [seed, ts]))

        # adds the powerflow results (P, Q, |U|) of this timestep on the mn dict, for future reference/comparison
        add_pf_result_to_mn_data!(mn_data["nw"]["$ts_id"], pf_results)

        # converts the powerflow results into (noisy or not) measurements
        add_measurements!(data, pf_results, σ_v, σ_d, σ_g, add_noise = add_noise, seed = seed, include_transfo_meas = true)

        # store this timestep in multinetwork dict
        mn_data["nw"]["$ts_id"]["meas"] = deepcopy(data["meas"])

    end
    return mn_data, real_volts
end


function _build_dictionary_entries(data, mn_data, nwentry::Int)
    for key in filter( x-> x != "meas", setdiff(keys(data), keys(mn_data)))
        mn_data["nw"]["$nwentry"][key] = deepcopy(data[key])
    end
end
""" 
timestep by timestep, reads the PQ profiles and adds the demand at each user
"""
function insert_profiles!(data, df, timestep; power_mult::Float64=1.)
    power_unit = data["settings"]["sbase"]
    for (_, load) in data["load"]
        p_id = load["parquet_id"]
        load["pd"] = [df[timestep, "P_kW_"*p_id]/power_unit]*power_mult
        load["qd"] = [df[timestep, "Q_kVAr_"*p_id]/power_unit]*power_mult
    end
end

function pf_solution_to_voltage_magnitudes!(sol::Dict)
    for (i, bus) in sol["solution"]["bus"]
        bus["vm"] = sqrt.( (bus["vr"][1:(end-1)] .- bus["vr"][end]).^2 + (bus["vi"][1:(end-1)] .- bus["vi"][end]).^2 )
    end
end

function add_pf_result_to_mn_data!(mn_data::Dict, pf_results::Dict)
    for (b, bus) in pf_results["solution"]["bus"]
        mn_data["bus"][b]["pf_vm"] = deepcopy(bus["vm"])
    end
    for (l, load) in pf_results["solution"]["load"]
        mn_data["load"][l]["pf_pd"] = deepcopy(load["pd"])
        mn_data["load"][l]["pf_qd"] = deepcopy(load["qd"])
    end
    mn_data["gen"]["1"]["pf_pg"] = deepcopy(pf_results["solution"]["gen"]["1"]["pg"])
    mn_data["gen"]["1"]["pf_qg"] = deepcopy(pf_results["solution"]["gen"]["1"]["qg"])
end

function add_measurements!(data::Dict, pf_results::Dict, σ_v::Float64, σ_d::Float64, σ_g::Float64; add_noise::Bool = true, seed::Int = 1, include_transfo_meas::Bool = false)
    
    data["meas"] = Dict{String, Any}()
    m_id = 1
    
    for (l, load) in pf_results["solution"]["load"]

        load_bus = data["load"][l]["load_bus"]
        randRNG = [_RAN.seed!(seed+load_bus+100*i) for i in 1:3] # not sure this really makes sense...

        vm_dst = [_DST.Normal{Float64}(res, σ_v) for res in pf_results["solution"]["bus"]["$load_bus"]["vm"]]
        pd_dst = [_DST.Normal{Float64}(res, σ_d) for res in load["pd"]]
        qd_dst = [_DST.Normal{Float64}(res, σ_d) for res in load["qd"]]

        if add_noise 
            vm_dst = [_DST.Normal{Float64}(_RAN.rand(randRNG[i], d), σ_v) for (i,d) in enumerate(vm_dst)] 
            pd_dst = [_DST.Normal{Float64}(_RAN.rand(randRNG[i], d), σ_d) for (i,d) in enumerate(pd_dst)] 
            qd_dst = [_DST.Normal{Float64}(_RAN.rand(randRNG[i], d), σ_d) for (i,d) in enumerate(qd_dst)] 
        end

        # add voltage magnitude measurement
        data["meas"]["$m_id"] = Dict{String, Any}(
            "var"    => :vm,
            "cmp"    => :bus,
            "cmp_id" => load_bus,
            "dst"    => vm_dst
        ) 

        # add active power measurement
        data["meas"]["$(m_id+1)"] = Dict{String, Any}(
            "var"    => :pd,
            "cmp"    => :load,
            "cmp_id" => data["load"][l]["index"],
            "dst"    => pd_dst
        ) 

        # add reactive power measurement
        data["meas"]["$(m_id+2)"] = Dict{String, Any}(
            "var"    => :qd,
            "cmp"    => :load,
            "cmp_id" => data["load"][l]["index"],
            "dst"    => qd_dst
        ) 

        m_id += 3

    end

    if include_transfo_meas

        m_id = maximum(parse.(Int, keys(data["meas"])))+1

        gen_bus = data["gen"]["1"]["gen_bus"]
        randRNG = [_RAN.seed!(seed+100*i) for i in 1:3] # not sure this really makes sense...

        vm_dst = [_DST.Normal{Float64}(res, σ_v) for res in pf_results["solution"]["bus"]["$gen_bus"]["vm"]]
        pg_dst = [_DST.Normal{Float64}(res, σ_g) for res in pf_results["solution"]["gen"]["1"]["pg"]]
        qg_dst = [_DST.Normal{Float64}(res, σ_g) for res in pf_results["solution"]["gen"]["1"]["qg"]]

        if add_noise 
            vm_dst = [_DST.Normal{Float64}(_RAN.rand(randRNG[i], d), σ_v) for (i,d) in enumerate(vm_dst)] 
            pg_dst = [_DST.Normal{Float64}(_RAN.rand(randRNG[i], d), σ_g) for (i,d) in enumerate(pg_dst)] 
            qg_dst = [_DST.Normal{Float64}(_RAN.rand(randRNG[i], d), σ_g) for (i,d) in enumerate(qg_dst)] 
        end

        # add voltage magnitude measurement
        data["meas"]["$m_id"] = Dict{String, Any}(
            "var"    => :vm,
            "cmp"    => :bus,
            "cmp_id" => gen_bus,
            "dst"    => vm_dst
        ) 

        # add active power measurement
        data["meas"]["$(m_id+1)"] = Dict{String, Any}(
            "var"    => :pg,
            "cmp"    => :gen,
            "cmp_id" => 1,
            "dst"    => pg_dst
        ) 

        # add reactive power measurement
        data["meas"]["$(m_id+2)"] = Dict{String, Any}(
            "var"    => :qg,
            "cmp"    => :gen,
            "cmp_id" => 1,
            "dst"    => qg_dst
        ) 

    end

end