import Random as _RAN

function prepare_math_eng_data(;feeder_name::String="30load-feeder", oh_or_ug::String="ug")
    eng = _PMD.parse_file(_IMP.NTW_DATA_DIR*"/"*feeder_name*"/Master_$(oh_or_ug).dss", data_model = _PMD.ENGINEERING, transformations=[_PMD.transform_loops!,_PMD.remove_all_bounds!])
    _IMP.rm_enwl_transformer!(eng)
    _IMP.reduce_enwl_lines_eng!(eng)
    eng["settings"]["sbase_default"] = 1
    data = _PMD.transform_data_model(eng, kron_reduce=false, phase_project=false)
    z_pu = (data["settings"]["voltage_scale_factor"]*data["settings"]["vbases_default"][collect(keys(data["settings"]["vbases_default"]))[1]])^2/(data["settings"]["power_scale_factor"])
    
    load_buses = [load["load_bus"] for (_, load) in data["load"]]
    _IMP.clean_4w_data!(data, profiles, merge_buses_diff_linecodes = false, eng = eng)
    _PMD.add_start_vrvi!(data)
    
    _IMP.make_loadbuses_loadbranches_singlephase!(data) 
    
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
        end
    end

    _IMP.add_length!(data, eng)
    return data, eng, z_pu
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