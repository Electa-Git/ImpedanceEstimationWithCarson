function build_multinetwork_dsse_data(data::Dict, df::_DF.DataFrame, pf_solver; timestep_set::Union{Vector{Int64}, UnitRange{Int64}} = 1:20, add_noise::Bool=false, seed::Int64=2, power_mult::Float64=1.0)
    mn_data = Dict{String, Any}()
    mn_data["multinetwork"] = true
    mn_data["nw"] = Dict{String, Any}()


    for key in ["name", "map", "se_settings", "data_model", "bus_lookup", "per_unit"]
        if haskey(data, key) mn_data[key] = data[key] end
    end

    real_volts = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["scenario_id", "time_step"]))
    real_vas = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["scenario_id", "time_step"]))

    power_unit = data["settings"]["sbase"]
    @info "The power unit in use is $power_unit and the load multiplication factor is $power_mult"

    for (ts_id, ts) in enumerate(timestep_set)
        
        # add info / initialize multinetwork dict
        mn_data["nw"]["$ts_id"] = Dict{String, Any}()
        mn_data["nw"]["$ts_id"]["original_timestep"] = ts
        mn_data["nw"]["$ts_id"]["per_unit"] = true

        _build_dictionary_entries(data, mn_data, ts_id)
        insert_profiles!(data, df, ts, power_mult=power_mult) # inserts NREL load profiles for powerflow

        # Solves the opf, i.e., which gives noiseless input
        pf_results = _PMD.solve_mc_opf(data, _PMD.IVRENPowerModel, pf_solver)
        
        if string(pf_results["termination_status"]) ∉ ["LOCALLY_SOLVED"]
            @warn "power flow result for time step $ts_id converged to: $(pf_results["termination_status"])"
        end
        # converts vr and vi to vm (phase to neutral)
        pf_solution_to_voltage_magnitudes!(pf_results) 
        pf_solution_to_voltage_angles!(pf_results) 

        push!(real_volts, vcat([pf_results["solution"]["bus"]["$(load["load_bus"])"]["vm"][1] for (l,load) in data["load"]], [seed, ts]))
        push!(real_vas, vcat([pf_results["solution"]["bus"]["$(load["load_bus"])"]["va"][1] for (l,load) in data["load"]], [seed, ts]))

        # adds the powerflow results (P, Q, |U|) of this timestep on the mn dict, for future reference/comparison
        add_pf_result_to_mn_data!(mn_data["nw"]["$ts_id"], pf_results)

        # converts the powerflow results into (noisy or not) measurements
        add_measurements!(data, pf_results; add_noise = add_noise, seed = seed, include_transfo_meas = false)

        # store this timestep in multinetwork dict
        mn_data["nw"]["$ts_id"]["meas"] = deepcopy(data["meas"]);

    end
    return mn_data, real_volts, real_vas
end

function build_multinetwork_dsse_data_with_shunts(data::Dict, df::_DF.DataFrame, pf_solver; timestep_set::Union{Vector{Int64}, UnitRange{Int64}} = 1:20, add_noise::Bool=false, loads_with_shunts::Vector{String} = ["1"], gs::Vector{Float64} = [50.], bs::Vector{Float64} = [15.], seed::Int64=2, power_mult::Float64=1.0)

    mn_data = Dict{String, Any}()
    mn_data["multinetwork"] = true
    mn_data["nw"] = Dict{String, Any}()

    for key in ["name", "map", "se_settings", "data_model", "bus_lookup", "per_unit"]
        if haskey(data, key) mn_data[key] = data[key] end
    end

    real_volts = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["scenario_id", "time_step"]))
    real_vas = _DF.DataFrame(fill([], length(data["load"])+2), vcat(["load_$(l)_ph_$(load["connections"][1])" for (l,load) in data["load"]], ["scenario_id", "time_step"]))
    
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

    for (ts_id, ts) in enumerate(timestep_set)
        
        # add info / initialize multinetwork dict
        mn_data["nw"]["$ts_id"] = Dict{String, Any}()
        mn_data["nw"]["$ts_id"]["original_timestep"] = ts
        mn_data["nw"]["$ts_id"]["per_unit"] = true

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
        push!(real_vas, vcat([pf_results["solution"]["bus"]["$(load["load_bus"])"]["va"][1] for (l,load) in data["load"]], [seed, ts]))

        # adds the powerflow results (P, Q, |U|) of this timestep on the mn dict, for future reference/comparison
        add_pf_result_to_mn_data!(mn_data["nw"]["$ts_id"], pf_results)

        # converts the powerflow results into (noisy or not) measurements
        add_measurements!(data, pf_results; add_noise = add_noise, seed = seed, include_transfo_meas = true)

        # store this timestep in multinetwork dict
        mn_data["nw"]["$ts_id"]["meas"] = deepcopy(data["meas"])

    end
    return mn_data, real_volts, real_vas
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

function pf_solution_to_voltage_angles!(sol::Dict)
    for (i, bus) in sol["solution"]["bus"]
        bus["va"] = [atan(bus["vi"][i] / bus["vr"][i] ) for i in 1:(length(bus["vr"])-1)]
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

function add_measurements!(data::Dict, pf_results::Dict; add_noise::Bool = true, seed::Int = 1, include_transfo_meas::Bool = false)
    
    data["meas"] = Dict{String, Any}()
    m_id = 1
    # original_sourcebus_id = collect(keys(data["settings"]["vbases_default"]))[1]

    max_volt_error = 0.005 # in percentual
    max_p_error = 0.01 # in percentual 
    max_q_error = 2*max_p_error # in percentual
    
    for (l, load) in pf_results["solution"]["load"]

        load_bus = data["load"][l]["load_bus"]
        randRNG = [_RAN.seed!(seed+load_bus+100*i) for i in 1:3] # not sure this seeding really makes sense...

        σ_v_mult = 1/3*max_volt_error #/(data["settings"]["vbases_default"][original_sourcebus_id]*data["settings"]["voltage_scale_factor"])
        σ_p_mult = 1/3*max_p_error    #/(data["settings"]["sbase_default"])
        σ_q_mult = 1/3*max_q_error    #/(data["settings"]["sbase_default"])

        vm_dst = [_DST.Normal{Float64}(res, σ_v_mult*res) for res in pf_results["solution"]["bus"]["$load_bus"]["vm"]]
        pd_dst = [_DST.Normal{Float64}(res, σ_p_mult*res) for res in load["pd"]]
        qd_dst = [_DST.Normal{Float64}(res, σ_q_mult*res) for res in load["qd"]]

        if add_noise 
            vm_meas = [_RAN.rand(randRNG[i], d) for (i,d) in enumerate(vm_dst)]
            pd_meas = [_RAN.rand(randRNG[i], d) for (i,d) in enumerate(pd_dst)]
            qd_meas = [_RAN.rand(randRNG[i], d) for (i,d) in enumerate(qd_dst)]
            vm_dst = [_DST.Normal{Float64}(m, σ_v_mult*m) for m in vm_meas] 
            pd_dst = [_DST.Normal{Float64}(m, maximum([σ_p_mult*m, 5e-4])) for m in pd_meas] 
            qd_dst = [_DST.Normal{Float64}(m, maximum([σ_q_mult*m, 5e-4])) for m in qd_meas] 
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
        nothing
        # NB: there used to be code but is deprecated. check older version to retrieve it if needed
    end
end