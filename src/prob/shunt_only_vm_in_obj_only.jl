function solve_shunt_only_vm_in_obj_only(data, solver)
    _PMD.solve_mc_model(data, _PMD.IVRQuadraticENPowerModel, solver, build_shunt_only_vm_in_obj_only, multinetwork=true)
end
"""
Builder for impedance estimation with carson's equation
"""
function build_shunt_only_vm_in_obj_only(pm::_PMD.IVRQuadraticENPowerModel)
    
    # Time varying variables
    for (n, _) in _PMD.nws(pm)
        _PMD.variable_mc_bus_voltage(pm, nw=n)
        variable_mc_branch_current(pm, nw=n)
        variable_mc_load_current(pm, nw=n)  
        variable_mc_generator_current(pm, nw=n)
        variable_mc_residual_vm(pm, nw=n)
    end

    # Time independent variables (i.e., shunt variables)
    variable_bus_shunt_admittance(pm)

    # Time varying constraints
    for (n, _) in _PMD.nws(pm)
        for i in _PMD.ids(pm, n, :bus)

            if i in _PMD.ids(pm, n, :ref_buses)
                _PMD.constraint_mc_theta_ref(pm, i, nw = n) # never fix reference voltage magnitude in DSSE, theta only
            end

            _PMD.constraint_mc_voltage_absolute(pm, i, nw = n) 
            _PMD.constraint_mc_voltage_pairwise(pm, i, nw = n) 
        end

        for id in _PMD.ids(pm, n, :meas)
            if _PMD.ref(pm, n, :meas, id)["var"] == :vm
                constraint_mc_residual_vm(pm, id, nw = n)
            end
        end

        for i in _PMD.ids(pm, n, :branch)
            constraint_mc_current_from_to(pm, i, nw = n)  
            constraint_mc_bus_voltage_drop(pm, i, nw = n) # function of series impedance (i.e., carson and all)
        end

        for i in _PMD.ids(pm, n, :bus)
            constraint_mc_current_balance(pm, i, nw = n)
        end

    end

    # Objective
    objective_minimize_residuals_vm(pm)
end