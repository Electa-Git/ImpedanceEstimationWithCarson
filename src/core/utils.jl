"""
Simple util that finds out to which bus component `cmp` (:load or :gen) `cmp_id` is connected
"""
get_cmp_bus(pm, nw, cmp, cmp_id) = 
    cmp == :gen ? (return _PMD.ref(pm, nw, cmp, cmp_id)["gen_bus"]) : (return _PMD.ref(pm, nw, cmp, cmp_id)["load_bus"])
"""
Simple util that gets current variables of generator and load buses
"""
function get_currents(pm, nw, cmp, cmp_id, conns, c) 
    if cmp == :gen 
        return _PMD.var(pm, nw, :crg_bus, cmp_id)[c], _PMD.var(pm, nw, :cig_bus, cmp_id)[c] 
    else    
        return _PMD.var(pm, nw, :crd_bus, cmp_id)[c], _PMD.var(pm, nw, :cid_bus, cmp_id)[c]
    end
end