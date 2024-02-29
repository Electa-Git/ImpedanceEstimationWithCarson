"""
    identical to `objective_mc_se` in PMDSE, which is not a dependency because it does not support 4-wire models, yet
"""
function objective_minimize_residuals(pm::_PMD.AbstractUnbalancedPowerModel)
    return JuMP.@objective(pm.model, Min,
    sum(
        sum(
            sum(_PMD.var(pm, nw, :res, i)[idx] for idx in 1:length(_PMD.var(pm, nw, :res, i)) )
        for i in _PMD.ids(pm, nw, :meas))
    for (nw, nw_ref) in _PMD.nws(pm) )
    )
end