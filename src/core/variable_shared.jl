"""
Variable that defines the residuals that end up in the objective function.
    Unfortunately cannot get rid of these because we need them as auxiliary vars to relax
    the WLAV in two inequality constraints.
"""
function variable_mc_residual(  pm::_PMD.AbstractUnbalancedPowerModel;
                                nw::Int=_PMD.nw_id_default, bounded::Bool=true,
                                report::Bool=true)

    connections = Dict(i => length(meas["dst"]) for (i,meas) in _PMD.ref(pm, nw, :meas) )

    res = _PMD.var(pm, nw)[:res] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:connections[i]], base_name = "$(nw)_res_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :meas, i), "res_start", c, 0.0)
            ) for i in _PMD.ids(pm, nw, :meas)
        )
    if bounded
        for (i, meas) in _PMD.ref(pm, nw, :meas), c in 1:connections[i]
            JuMP.set_lower_bound(res[i][c], 0.0)
            res_max = haskey(_PMD.ref(pm, nw, :meas, i), "res_max") ? meas["res_max"] : 5e3
            JuMP.set_upper_bound(res[i][c], res_max)
        end
    end
    report && _IM.sol_component_value(pm,:pmd, nw, :meas, :res, _PMD.ids(pm, nw, :meas), res)
end
"""
	function variable_mc_branch_current(
		pm::AbstractExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates total current variables `:cr` and `:ci`,
series current variables `:csr` and `:csi`,
and placeholder dictionaries for the terminal current flows `:cr_bus` and `:ci_bus`
    NOTE: w.r.t. PMD we only have a series current defined for the neutral,
    in the other cases, the series is identical to the total current (no shunt current).
"""
function variable_mc_branch_current(pm::_PMD.AbstractExplicitNeutralIVRModel; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)

    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in _PMD.ref(pm, nw, :branch))

    cr = _PMD.var(pm, nw)[:cr] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_cr_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :branch, l), "cr_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch)
    )

    ci = _PMD.var(pm, nw)[:ci] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_ci_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :branch, l), "ci_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch)
    )

    _PMD.var(pm, nw)[:cr_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    _PMD.var(pm, nw)[:ci_bus] = Dict{Tuple{Int,Int,Int}, Any}()

    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :branch, :cr_fr, :cr_to, _PMD.ref(pm, nw, :arcs_branch_from), _PMD.ref(pm, nw, :arcs_branch_to), cr)
    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :branch, :ci_fr, :ci_to, _PMD.ref(pm, nw, :arcs_branch_from), _PMD.ref(pm, nw, :arcs_branch_to), ci)

end

function variable_mc_load_current(pm::_PMD.AbstractExplicitNeutralIVRModel; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_load_current_real(pm; nw=nw, kwargs...)
    variable_mc_load_current_imag(pm; nw=nw, kwargs...)
end

function variable_mc_load_current_real(pm::_PMD.AbstractExplicitNeutralIVRModel;
                                 nw::Int=_IM.nw_id_default, bounded::Bool=true, report::Bool=true)

    int_dim = Dict(i => _PMD._infer_int_dim_unit(load, false) for (i,load) in _PMD.ref(pm, nw, :load))

    crd = _PMD.var(pm, nw)[:crd] = Dict(i => JuMP.@variable(pm.model,
          [c in 1:int_dim[i]], base_name="$(nw)_crd_$(i)"
        ) for i in _PMD.ids(pm, nw, :load)
    )

    _PMD.var(pm, nw)[:crd_bus] = Dict{Int, Any}()

    for i in _PMD.ids(pm, nw, :load)
        _PMD.var(pm, nw, :crd_bus)[i] = _PMD._merge_bus_flows(pm, [crd[i]..., -sum(crd[i])], _PMD.ref(pm, nw, :load, i)["connections"])
    end

    report && _IM.sol_component_value(pm, :pmd, nw, :load, :crd, _PMD.ids(pm, nw, :load), crd)
end

function variable_mc_load_current_imag(pm::_PMD.AbstractExplicitNeutralIVRModel; nw::Int=_IM.nw_id_default, bounded::Bool=true, report::Bool=true, meas_start::Bool=false)

    int_dim = Dict(i => _PMD._infer_int_dim_unit(load, false) for (i,load) in _PMD.ref(pm, nw, :load))

    cid = _PMD.var(pm, nw)[:cid] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cid_$(i)"
        ) for i in _PMD.ids(pm, nw, :load)
    )
    
    _PMD.var(pm, nw)[:cid_bus] = Dict{Int, Any}()

    for i in _PMD.ids(pm, nw, :load)
        _PMD.var(pm, nw, :cid_bus)[i] = _PMD._merge_bus_flows(pm, [cid[i]..., -sum(cid[i])], _PMD.ref(pm, nw, :load, i)["connections"])
    end

    report && _IM.sol_component_value(pm, :pmd, nw, :load, :cid, _PMD.ids(pm, nw, :load), cid)

end

function variable_mc_generator_current(pm::_PMD.AbstractExplicitNeutralIVRModel; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)

    int_dim = Dict(i => _PMD._infer_int_dim_unit(gen, false) for (i,gen) in _PMD.ref(pm, nw, :gen))

    crg = _PMD.var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_crg_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :gen, i), "crg_start", c, 0.0)
        ) for i in _PMD.ids(pm, nw, :gen)
    )

    cig = _PMD.var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cig_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :gen, i), "cig_start", c, 0.0)
        ) for i in _PMD.ids(pm, nw, :gen)
    )

    _PMD.var(pm, nw)[:crg_bus] = Dict{Int, Any}()
    _PMD.var(pm, nw)[:cig_bus] = Dict{Int, Any}()

    for i in _PMD.ids(pm, nw, :gen)
        _PMD.var(pm, nw, :crg_bus)[i] = _PMD._merge_bus_flows(pm,[crg[i]..., -sum(crg[i])], _PMD.ref(pm, nw, :gen, i)["connections"])
        _PMD.var(pm, nw, :cig_bus)[i] = _PMD._merge_bus_flows(pm,[cig[i]..., -sum(cig[i])], _PMD.ref(pm, nw, :gen, i)["connections"])
    end

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen, :crg, _PMD.ids(pm, nw, :gen), crg)
    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen, :cig, _PMD.ids(pm, nw, :gen), cig)

end