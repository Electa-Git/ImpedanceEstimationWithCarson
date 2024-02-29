"""
Carson's equations give impedance in Ω/km. This multiplies to get the Ω.
SO expect the result to be in km!!
"""
function variable_segment_length(pm::_PMD.AbstractExplicitNeutralIVRModel; bounded::Bool = true, report::Bool=true)

    nw = 1 # variable only defined at the first timestep is what this is. not elegant, to change, but just to be sure in the meantime

    l = _PMD.var(pm, nw)[:l] = Dict(i => JuMP.@variable(pm.model,
                base_name="length_$(i)",
                lower_bound = 1e-5, # lower bound is 100 cm unless specified otherwise (see below)
                start = _PMD.comp_start_value(_PMD.ref(pm, nw, :branch, i), "l_start", 0.001) # 1 meters default start value if no `l_start` exists
                ) for i in _PMD.ids(pm, nw, :branch) 
                     if _PMD.ref(pm, nw, :branch, i, "untrustworthy_branch")
            ) 

    if bounded
        for (i,branch) in _PMD.ref(pm, nw, :branch)
            if _PMD.ref(pm, nw, :branch, i, "untrustworthy_branch")
                haskey(branch, "l_max") ?  JuMP.set_upper_bound(l[i], branch["l_max"]) : JuMP.set_upper_bound(l[i],  100.0) # upper bound is 10 km
                if haskey(branch, "l_min") JuMP.set_lower_bound(l[i], branch["l_min"]) end # otherwise takes the one stated in the definition above
            end
        end
    end

    rprt_ids = [i for i ∈ _PMD.ids(pm, nw, :branch) if _PMD.ref(pm, nw, :branch, i, "untrustworthy_branch")] 
    report && _IM.sol_component_value(pm, :pmd, nw, :branch, :l, rprt_ids, l)
end
"""
Only one variable per conductor per linecode
"""
function variable_cross_section_and_gmr(pm::_PMD.AbstractExplicitNeutralIVRModel, μ; bounded::Bool = true, report::Bool=true)

    nw = 1 # variable only defined at the first timestep is what this is. not elegant, to change, but just to be sure in the meantime

    conn = Dict( i => [wire for wire in 1:l["n_wires"]] for (i,l) in _PMD.ref(pm, 1, :linecode_map))
    
    A_p = _PMD.var(pm, nw)[:A_p] = Dict(i => JuMP.@variable(pm.model,
                [t in conn[i]],
                base_name="A_p_$(i)",
                lower_bound = 5, # 
                start = _PMD.comp_start_value(linecode, "A_p_start", t, 30.) #
                ) for (i, linecode) in _PMD.ref(pm, 1, :linecode_map) 
            ) 

    gmr = _PMD.var(pm, nw)[:gmr] = Dict(i => JuMP.@variable(pm.model,
            [t in conn[i]],
            base_name="gmr_$(i)",
            lower_bound = 1, # 
            start = _PMD.comp_start_value(linecode, "gmr_start", t, 3.) #
            ) for (i, linecode) in _PMD.ref(pm, 1, :linecode_map) 
        ) 

    if bounded
        for (i,linecode) in _PMD.ref(pm, 1, :linecode_map)
            for idx in conn[i]
                   haskey(linecode, "A_p_max") ?  JuMP.set_upper_bound(A_p[i][idx], linecode["A_p_max"][idx]) : JuMP.set_upper_bound(A_p[i][idx],  200.0) #
                 if haskey(linecode, "A_p_min")    JuMP.set_lower_bound(A_p[i][idx], linecode["A_p_min"][idx]) end # otherwise takes the one stated in the definition above
            end
        end
    end

    for (i,l) in _PMD.ref(pm, nw, :linecode_map)
        for j in 1:l["n_wires"]
            JuMP.@constraint(pm.model, gmr[i][j]^2 == exp(-μ/2)*A_p[i][j]/π)
        end
    end

    report && _IM.sol_component_value(pm, :pmd, 1, :linecode_map, :A_p, _PMD.ids(pm, 1, :linecode_map), A_p)
    report && _IM.sol_component_value(pm, :pmd, 1, :linecode_map, :gmr, _PMD.ids(pm, 1, :linecode_map), gmr)

end
"""
coordinates are in mm
"""
function variable_coordinates_and_distances(pm::_PMD.AbstractExplicitNeutralIVRModel; report::Bool=true) 

    nw = 1 # variable only defined at the first timestep is what this is. not elegant, to change, but just to be sure in the meantime

    two_w   = [ i for (i,l) in _PMD.ref(pm, nw, :linecode_map) if l["n_wires"] == 2 ]
    three_w = [ i for (i,l) in _PMD.ref(pm, nw, :linecode_map) if l["n_wires"] == 3 ]
    four_w  = [ i for (i,l) in _PMD.ref(pm, nw, :linecode_map) if l["n_wires"] == 4 ]

    xcoord = _PMD.var(pm, nw)[:xcoord] = Dict(i => JuMP.@variable(pm.model,
                [j = 2:l["n_wires"]], # starts from 2 because first cond is (0,0)
                base_name="xcoord_code_$(i)",
                lower_bound = 0, 
                upper_bound =  2e3, 
                ) for (i,l) in _PMD.ref(pm, 1, :linecode_map) if i ∈ vcat(three_w, four_w)
            )

    for (i,_) in _PMD.ref(pm, 1, :linecode_map) 
        if i ∈ vcat(three_w, four_w)
            JuMP.@constraint(pm.model, xcoord[i][3]>=xcoord[i][2])
            if i ∈ four_w
                JuMP.@constraint(pm.model, xcoord[i][4]>=xcoord[i][3])
            end
        end
    end
            
    ycoord = _PMD.var(pm, nw)[:ycoord] = Dict(i => JuMP.@variable(pm.model,
                [j = 2:l["n_wires"]], # starts from 2 because first cond is (0,0)
                base_name="ycoord_code_$(i)",
                lower_bound = -2e3, 
                upper_bound =  0, 
            ) for (i,l) in _PMD.ref(pm, 1, :linecode_map) if i ∈ vcat(three_w, four_w)
        )

    dij_2w = Dict(i => JuMP.@variable(pm.model,
            [1:1], # making it a vector for consistency with the rest (access through indices later on)
            base_name="dij_2w_code_$i",
            lower_bound = 5, 
            upper_bound = 2e3, 
            ) for (i,l) in _PMD.ref(pm, 1, :linecode_map) if i ∈ two_w
        ) 

    distance_indices = [(1,2), (1,3), (2,3), (1,4), (2,4), (3,4)]

    dij_3w = Dict(i => JuMP.@variable(pm.model,
                [j = 1:3],
                base_name="dij_3w_code_$(i)",
                lower_bound = 10,
                upper_bound =  2e3, # two meters
                ) for (i,l) in _PMD.ref(pm, 1, :linecode_map) if i ∈ three_w
            )

    dij_4w = Dict(i => JuMP.@variable(pm.model,
                [j = 1:6],
                base_name="dij_3w_code_$(i)",
                lower_bound = 10,
                upper_bound =  2e3, 
                ) for (i,l) in _PMD.ref(pm, 1, :linecode_map) if i ∈ four_w
            )

    _PMD.var(pm, nw)[:dij] = merge(dij_2w, dij_3w, dij_4w)

    for (i,l) in _PMD.ref(pm, nw, :linecode_map)
        if i ∈ three_w
            for (id, ij) in enumerate(distance_indices[1:3])
                x0 = ij[1] == 1 ? 0. : xcoord[i][ij[1]] # because the coordinates of the first conductor are fixed to zero
                y0 = ij[1] == 1 ? 0. : ycoord[i][ij[1]]
                JuMP.@constraint(pm.model, dij_3w[i][id]^2 == (x0-xcoord[i][ij[2]])^2+(y0-ycoord[i][ij[2]])^2 )
            end
        elseif i ∈ four_w
            for (id, ij) in enumerate(distance_indices)
                x0 = ij[1] == 1 ? 0. : xcoord[i][ij[1]] # because the coordinates of the first conductor are fixed to zero
                y0 = ij[1] == 1 ? 0. : ycoord[i][ij[1]]
                JuMP.@constraint(pm.model, dij_4w[i][id]^2 == (x0-xcoord[i][ij[2]])^2+(y0-ycoord[i][ij[2]])^2 )
            end
        end
    end

    report_idx_2w_dij = [i for (i,l) in _PMD.ref(pm, 1, :linecode_map) if i ∈ two_w]
    report_idx_3w_dij = [i for (i,l) in _PMD.ref(pm, 1, :linecode_map) if i ∈ three_w]
    report_idx_4w_dij = [i for (i,l) in _PMD.ref(pm, 1, :linecode_map) if i ∈ four_w]
    report_idx_coords = [i for (i,l) in _PMD.ref(pm, 1, :linecode_map) if i ∈ vcat(three_w, four_w)]

    report && _IM.sol_component_value(pm, :pmd, 1, :linecode_map, :dij_2w, report_idx_2w_dij, dij_2w)
    report && _IM.sol_component_value(pm, :pmd, 1, :linecode_map, :dij_3w, report_idx_3w_dij, dij_3w)
    report && _IM.sol_component_value(pm, :pmd, 1, :linecode_map, :dij_4w, report_idx_4w_dij, dij_4w)
    report && _IM.sol_component_value(pm, :pmd, 1, :linecode_map, :xcoord, report_idx_coords, xcoord)
    report && _IM.sol_component_value(pm, :pmd, 1, :linecode_map, :ycoord, report_idx_coords, ycoord)

end

function variable_bus_shunt_admittance(pm::_PMD.AbstractExplicitNeutralIVRModel; bounded::Bool=true, report::Bool=true)

    nw = 1 # variable only defined at the first timestep is what this is. not elegant, to change, but just to be sure in the meantime

    g_sh = _PMD.var(pm, nw)[:g_sh] = Dict(i => JuMP.@variable(pm.model,
                base_name="g_sh_$(i)",
                lower_bound = 0., 
                start = _PMD.comp_start_value(_PMD.ref(pm, nw, :bus, i), "g_sh_start", 5.) 
                ) for i in _PMD.ids(pm, nw, :bus) 
                     if any(_PMD.ref(pm, nw, :bus, i, "imp_grounded"))
            )
    
    b_sh = _PMD.var(pm, nw)[:b_sh] = Dict(i => JuMP.@variable(pm.model,
                base_name="b_sh_$(i)",
                lower_bound = 0., 
                start = _PMD.comp_start_value(_PMD.ref(pm, nw, :bus, i), "b_sh_start", 5.) 
                ) for i in _PMD.ids(pm, nw, :bus) 
                     if any(_PMD.ref(pm, nw, :bus, i, "imp_grounded"))
            )

    if bounded
        for (i,bus) in _PMD.ref(pm, nw, :bus)
            if any(_PMD.ref(pm, nw, :bus, i, "imp_grounded"))
                # g_shunt  
                haskey(bus, "g_sh_max") ?  JuMP.set_upper_bound(g_sh[i], bus["g_sh_max"]) : JuMP.set_upper_bound(g_sh[i], 100) # TODO: better default val
                if haskey(bus, "g_sh_min") JuMP.set_lower_bound(g_sh[i], bus["g_sh_min"]) end
                # b_shunt
                haskey(bus, "b_sh_max") ?  JuMP.set_upper_bound(b_sh[i], bus["b_sh_max"]) : JuMP.set_upper_bound(b_sh[i], 100) # TODO: better default val
                if haskey(bus, "b_sh_min") JuMP.set_lower_bound(b_sh[i], bus["b_sh_min"]) end
            end
        end
    end

    rprt_ids = [i for i in _PMD.ids(pm, nw, :bus) if any(_PMD.ref(pm, nw, :bus, i, "imp_grounded"))] 
    report && _IM.sol_component_value(pm, :pmd, nw, :bus, :g_sh, rprt_ids, g_sh)
    report && _IM.sol_component_value(pm, :pmd, nw, :bus, :b_sh, rprt_ids, b_sh)

end

function variable_carson_series(pm)

    variable_segment_length(pm) # branch length is always in the variable set regardless (one length per each branch)
    variable_coordinates_and_distances(pm) # do not depend on any material property

    # ↓ temperature and material properties that might be variable or not
    # TODO: build a separate one per linecode  instead of assuming all same material
    μ = haskey(_PMD.ref(pm, 1, :settings), "mu_rel") ? _PMD.ref(pm, 1, :settings)["mu_rel"] : variable_mu_rel(pm)
    if μ isa Real @info "μ is being fixed" end

    # ↓ variables that depend on material properties
    variable_cross_section_and_gmr(pm, μ)

end

function carson_constants(pm)
    z_pu = _PMD.ref(pm, 1, :settings)["z_pu"]
    c₁ = 3.28084e-3
    c₂ = 8.0252
    r_pq = 0.049348 # not a variable at all --> TODO explain in the paper and make a table of the number of variables
    return c₁, c₂, z_pu, r_pq
end
