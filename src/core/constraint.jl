"""
    constraint_mc_residual

Equality constraint that describes the residual definition, which depends on the
criterion assigned to each individual measurement in data["meas"]["m"]["crit"].
"""
function constraint_mc_residual(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_IM.nw_id_default)

    cmp_id = get_cmp_id(pm, nw, i)
    res = _PMD.var(pm, nw, :res, i)
    dst = _PMD.ref(pm, nw, :meas, i, "dst")
    rsc = _PMD.ref(pm, 1, :settings)["rescaler"]
    cmp =  _PMD.ref(pm, nw, :meas, i, "cmp")

    if _PMD.ref(pm, nw, :meas, i, "var") == :l
        conns = [1]
    else
        conns =  get_active_connections(pm, nw, cmp, cmp_id)
        bus_id = _PMD.ref(pm, nw, :meas, i, "var") == :vm ? cmp_id : get_cmp_bus(pm, nw, cmp, cmp_id)    
        # get neutral voltage of the bus
        vr_n = _PMD.var(pm, nw, :vr, bus_id)[4]
        vi_n = _PMD.var(pm, nw, :vi, bus_id)[4]
    end

    for (idx, c) in enumerate(setdiff(conns, [4])) # setdiff(...,[4]) excludes the neutral, which is never measured!

        if  _PMD.ref(pm, nw, :meas, i, "var") == :vm    
            
            μ, σ = (_DST.mean(dst[idx]), _DST.std(dst[idx]))

            vr_p = _PMD.var(pm, nw, :vr, bus_id)[c]
            vi_p = _PMD.var(pm, nw, :vi, bus_id)[c]
            
            JuMP.@constraint(pm.model,
                res[idx] * rsc * σ >=   ((vr_p-vr_n)^2+(vi_p-vi_n)^2 - μ^2) 
            )
            JuMP.@constraint(pm.model,
                res[idx] * rsc * σ >= - ((vr_p-vr_n)^2+(vi_p-vi_n)^2 - μ^2)
            )
        
        elseif _PMD.ref(pm, nw, :meas, i, "var") ∈ [:pg, :pd, :qg, :qd]
            
            μ, σ = (_DST.mean(dst[idx]), _DST.std(dst[idx]))

            cr, ci = get_currents(pm, nw, cmp, cmp_id, conns, c) # gets load or generator current variables depending on which of the two it is
            
            vr_p = _PMD.var(pm, nw, :vr, bus_id)[c]
            vi_p = _PMD.var(pm, nw, :vi, bus_id)[c]

            if _PMD.ref(pm, nw, :meas, i, "var") ∈ [:pg, :pd] 
                JuMP.@constraint(pm.model,
                    res[idx] * rsc * σ >=   ( (vr_p - vr_n)*cr+(vi_p-vi_n)*ci - μ) 
                )
                JuMP.@constraint(pm.model,
                    res[idx] * rsc * σ >= - ( (vr_p - vr_n)*cr+(vi_p-vi_n)*ci - μ)
                )
            else
                JuMP.@constraint(pm.model,
                    res[idx] * rsc * σ >=   ( -(vr_p-vr_n)*ci+(vi_p-vi_n)*cr - μ) 
                )
                JuMP.@constraint(pm.model,
                    res[idx] * rsc * σ >= - ( -(vr_p-vr_n)*ci+(vi_p-vi_n)*cr - μ)
                )
            end
        elseif _PMD.ref(pm, nw, :meas, i, "var") == :l

                μ, σ = (_DST.mean(dst), _DST.std(dst))

                JuMP.@constraint(pm.model,
                    res[idx] * σ >=   (_PMD.var(pm, nw, :l, _PMD.ref(pm, nw, :meas, i, "cmp_id")) - μ) 
                )
                JuMP.@constraint(pm.model,
                    res[idx] * σ >= - (_PMD.var(pm, nw, :l, _PMD.ref(pm, nw, :meas, i, "cmp_id"))- μ)
                )
        else
            @error "Sorry, measurement $(_PMD.ref(pm, nw, :meas, i, "var")) not supported (yet)."
        end
    end
end
"""
    get_cmp_id(pm, nw, i)
Retrieves the id of component i. This is a tuple if the component is a branch. Otherwise, it is an Int.
"""
function get_cmp_id(pm, nw, i)
    if  _PMD.ref(pm, nw, :meas, i, "cmp") == :branch
        branch_id = _PMD.ref(pm, nw, :meas, i, "cmp_id")
        cmp_id = (branch_id, _PMD.ref(pm,nw,:branch, branch_id)["f_bus"], _PMD.ref(pm,nw,:branch,branch_id)["t_bus"])
    else
        cmp_id = _PMD.ref(pm, nw, :meas, i, "cmp_id")
    end
    return cmp_id
end
"""
    function get_active_connections(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, cmp_type::Symbol, cmp_id::Int)
Returns the list of terminals, connections or t_ and f_connections, depending on the type of the component.
"""
function get_active_connections(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, cmp_type::Symbol, cmp_id::Int)
    if cmp_type == :bus
       active_conn = _PMD.ref(pm, nw, :bus, cmp_id)["terminals"]
   elseif cmp_type ∈ [:gen, :load]
       active_conn = _PMD.ref(pm, nw, cmp_type, cmp_id)["connections"]
   elseif cmp_type == :branch
       active_conn = intersect(_PMD.ref(pm, nw, :branch, cmp_id)["f_connections"], _PMD.ref(pm, nw, :branch, cmp_id)["t_connections"])
   else
       error("Measurements for component of type $cmp_type are not supported")
   end
   return active_conn
end

function constraint_mc_current_balance(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_IM.nw_id_default)

    bus = _PMD.ref(pm, 1, :bus, i)
    bus_arcs = _PMD.ref(pm, 1, :bus_arcs_conns_branch, i)
    bus_gens = _PMD.ref(pm, 1, :bus_conns_gen, i)
    bus_loads = _PMD.ref(pm, 1, :bus_conns_load, i)
    bus_shunts = _PMD.ref(pm, nw, :bus_conns_shunt, i)

    terminals = bus["terminals"]
    imp_grounded = bus["imp_grounded"]

    vr = _PMD.var(pm, nw, :vr, i)
    vi = _PMD.var(pm, nw, :vi, i)

    cr    = get(_PMD.var(pm, nw),   :cr_bus,  Dict()); _PMD._check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(_PMD.var(pm, nw),   :ci_bus,  Dict()); _PMD._check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crd   = get(_PMD.var(pm, nw),   :crd_bus, Dict()); _PMD._check_var_keys(crd, bus_loads, "real current", "load")
    cid   = get(_PMD.var(pm, nw),   :cid_bus, Dict()); _PMD._check_var_keys(cid, bus_loads, "imaginary current", "load")
    crg   = get(_PMD.var(pm, nw),   :crg_bus, Dict()); _PMD._check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(_PMD.var(pm, nw),   :cig_bus, Dict()); _PMD._check_var_keys(cig, bus_gens, "imaginary current", "generator")

    imperfectly_grounded_terminals = [(idx, t) for (idx,t) in enumerate(terminals) if imp_grounded[idx]]
    terminals_not_grounded_at_all = [(idx, t) for (idx,t) in enumerate(terminals) if !bus["grounded"][idx]] # the PMD way
    total_ungrounded_terminals = unique(vcat(imperfectly_grounded_terminals, terminals_not_grounded_at_all))

    Gt, Bt = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    for (idx, t) in total_ungrounded_terminals
        if isempty(imperfectly_grounded_terminals)           
            # do nothing
        else
            n = size(Gt)[1]
            Gt = Matrix(undef, n, n)
            Bt = Matrix(undef, n, n)
            Gt .= 0.
            Bt .= 0.
            Gt[n,n] = _PMD.var(pm, 1, :g_sh, i) # this brutal overwriting allows you to leave the bs,gs, in the shunt entries in the dictionary ;)
            Bt[n,n] = _PMD.var(pm, 1, :b_sh, i) # this brutal overwriting allows you to leave the bs,gs, in the shunt entries in the dictionary ;)
        end

        JuMP.@constraint(pm.model,    sum(cr[a][t] for (a, conns) in bus_arcs if t in conns)
                                    ==
                                      sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(crd[d][t]         for (d, conns) in bus_loads if t in conns) 
                                    - sum(Gt[idx, jdx]*vr[u] - Bt[idx, jdx]*vi[u] for (jdx,u) in total_ungrounded_terminals)
                                    )

        JuMP.@constraint(pm.model,    sum(ci[a][t] for (a, conns) in bus_arcs if t in conns)
                                    ==
                                      sum(cig[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(cid[d][t]         for (d, conns) in bus_loads if t in conns) 
                                    - sum(Gt[idx, jdx]*vi[u] + Bt[idx, jdx]*vr[u] for (jdx,u) in total_ungrounded_terminals)
                                    )
    end
end
"""
Length is excluded, that's added sperately for each branch
"""
function carson_impedance_expressions(pm::_PMD.AbstractExplicitNeutralIVRModel)

    c₁, c₂, z_pu, r_pq = carson_constants(pm)

    # get materials properties
    T = haskey(_PMD.ref(pm, 1, :temperature), "temperature_value") ? _PMD.ref(pm, 1, :temperature)["temperature_value"] : @info "Temperature value not given" 
    if haskey(_PMD.ref(pm, 1, :temperature), "temperature_value") @info "Temperature is being fixed at value $(_PMD.ref(pm, 1, :temperature)["temperature_value"]) °C" end
    ρ = haskey(_PMD.ref(pm, 1, :rho), "rho_value") ? _PMD.ref(pm, 1, :rho)["rho_value"] : @info "rho value not given"
    if haskey(_PMD.ref(pm, 1, :rho), "rho_value") @info "ρ is being fixed at value $(_PMD.ref(pm, 1, :rho)["rho_value"])" end
    α = haskey(_PMD.ref(pm, 1, :alpha), "alpha_value") ? _PMD.ref(pm, 1, :alpha)["alpha_value"] : @info "alpha value not given"
    if haskey(_PMD.ref(pm, 1, :alpha), "alpha_value") @info "α is being fixed at value $(_PMD.ref(pm, 1, :alpha)["alpha_value"])" end

    # initialize expression containers
    _PMD.var(pm, 1)[:r_ac] = Dict{Int, Any}()
    _PMD.var(pm, 1)[:x] = Dict{Int, Any}()

    for i in _PMD.ids(pm, 1, :linecode_map)

        n_wires = _PMD.ref(pm, 1, :linecode_map, i)["n_wires"]

        # get carson variables
        A_p = _PMD.var(pm, 1, :A_p, i) 
        gmr = _PMD.var(pm, 1, :gmr, i)
        Dij = _PMD.var(pm, 1, :dij, i)  
        
        if haskey(_PMD.ref(pm, 1, :linecode_map, i), "r_ac")     
            @info "r_ac's are being fixed"
            _PMD.var(pm, 1, :r_ac)[i] = _PMD.ref(pm, 1, :linecode_map, i)["r_ac"]
        else
            _PMD.var(pm, 1, :r_ac)[i] = JuMP.@expression(pm.model, [ρ/A_p[i]*(1+(α*(T-20))) for i in 1:n_wires])
        end

        # build X matrix
        x = zeros(JuMP.NonlinearExpr, (n_wires, n_wires))
        for c in CartesianIndices(x)
            if c[1] == c[2]
                x[c] = JuMP.@expression(pm.model, 
                                            0.062832 * ( c₂ + log(1) - log( c₁*gmr[c[1]] ) ) 
                                        )
            elseif c[1] > c[2]
                # do nothing, symmetry exploited below
            else
                distance_indices = [(1,2), (1,3), (2,3), (1,4), (2,4), (3,4)]
                idx = findfirst(x->x == (c[1], c[2]), distance_indices)
                x[c[1], c[2]] = x[c[2], c[1]] = JuMP.@expression(pm.model, 
                                            0.062832 * (c₂ + log(1) - log( c₁*Dij[idx] )    ) 
                                        )
            end
        end

        _PMD.var(pm, 1, :x)[i] = x

    end
end

function constraint_mc_bus_voltage_drop(pm::_PMD.AbstractExplicitNeutralIVRModel, i::Int; nw::Int=_PMD.nw_id_default)

    # get all info and cars that do not depend on impedance construction
    branch = _PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_arc = (i, f_bus, t_bus)

    vr_fr = [_PMD.var(pm, nw, :vr, f_bus)[c] for c in branch["f_connections"]]
    vi_fr = [_PMD.var(pm, nw, :vi, f_bus)[c] for c in branch["f_connections"]]

    vr_to = [_PMD.var(pm, nw, :vr, t_bus)[c] for c in branch["t_connections"]]
    vi_to = [_PMD.var(pm, nw, :vi, t_bus)[c] for c in branch["t_connections"]]

    cr_fr = [_PMD.var(pm, nw, :cr_bus, f_arc)[c] for c in branch["f_connections"]]
    ci_fr = [_PMD.var(pm, nw, :ci_bus, f_arc)[c] for c in branch["f_connections"]]

    # build impedances 
    if !_PMD.ref(pm, 1, :branch, i, "untrustworthy_branch")
        R = branch["br_r"]
        X = branch["br_x"] 
    else
        c₁, c₂, z_pu, r_pq = carson_constants(pm)
        l   = _PMD.var(pm, 1, :l, i)
        # get linecode stuff for this branch
        lc    = _PMD.ref(pm, 1, :branch, i)["orig_linecode"]
        lc_id = [key for (key, val) in _PMD.ref(pm, 1, :linecode_map) if val["name"] == lc][1]
        n_wires = [val["n_wires"] for (key, val) in _PMD.ref(pm, 1, :linecode_map) if val["name"] == lc][1]

        r = fill(r_pq, (n_wires, n_wires))
        r += _LA.diagm(_PMD.var(pm, 1, :r_ac, lc_id)) # TODO: make aux var for this one to streamline ad?
        x = _PMD.var(pm, 1, :x, lc_id) # TODO: make aux var for this one to streamline ad?

        R = JuMP.@expression(pm.model, r.*l./z_pu)
        X = JuMP.@expression(pm.model, x.*l./z_pu)

    end

    JuMP.@constraint(pm.model, vr_to .== vr_fr - R*cr_fr + X*ci_fr)
    JuMP.@constraint(pm.model, vi_to .== vi_fr - R*ci_fr - X*cr_fr)

end

function constraint_mc_current_from_to(pm::_PMD.AbstractExplicitNeutralIVRModel, i::Int; nw::Int=_IM.nw_id_default)
    
    branch = _PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_arc = (i, f_bus, t_bus)
    t_arc = (i, t_bus, f_bus)
    f_connections = branch["f_connections"]
    t_connections = branch["t_connections"]

    cr_fr =  _PMD.var(pm, nw, :cr, f_arc)
    ci_fr =  _PMD.var(pm, nw, :ci, f_arc)

    cr_to =  _PMD.var(pm, nw, :cr, t_arc)
    ci_to =  _PMD.var(pm, nw, :ci, t_arc)

    _PMD.var(pm, nw, :cr_bus)[f_arc] = _PMD._merge_bus_flows(pm, cr_fr, f_connections)
    _PMD.var(pm, nw, :ci_bus)[f_arc] = _PMD._merge_bus_flows(pm, ci_fr, f_connections)
    _PMD.var(pm, nw, :cr_bus)[t_arc] = _PMD._merge_bus_flows(pm, cr_to, t_connections)
    _PMD.var(pm, nw, :ci_bus)[t_arc] = _PMD._merge_bus_flows(pm, ci_to, t_connections)

    JuMP.@constraint(pm.model, cr_fr .== -cr_to)
    JuMP.@constraint(pm.model, ci_fr .== -ci_to)

end
