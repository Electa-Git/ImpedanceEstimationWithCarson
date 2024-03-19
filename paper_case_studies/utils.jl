import Random as _RAN

function add_length_bounds!(mn_data::Dict, length_bounds_percval::Float64)
    for (_, branch) in mn_data["nw"]["1"]["branch"]
        branch["l_min"] = branch["orig_length"]*(1-length_bounds_percval)/1000 # / 1000 because length data is in m but length var is in km
        branch["l_max"] = branch["orig_length"]*(1+length_bounds_percval)/1000 # / 1000 because length data is in m but length var is in km
    end
    return mn_data
end

function add_material_properties_for_ug_noshunt_eulvtf!(mn_data, eng)
    
    material_resist_dict = Dict(
        "ugsc_25cu_xlpe/nyl/pvc_ug_2w_bundled" => 26.210307508899646,
        "ugsc_16cu_xlpe/nyl/pvc_ug_2w_bundled" => 26.600493316475497,
        "uglv_120cu_xlpe/nyl/pvc_ug_4w_bundled" => 27.782140561270225
    )

    mn_data["nw"]["1"]["linecode_map"] = Dict{Int, Any}() 
    for (id, code) in enumerate(keys(eng["linecode"]))
        if code ∈ ["ugsc_16cu_xlpe/nyl/pvc_ug_2w_bundled", "ugsc_25cu_xlpe/nyl/pvc_ug_2w_bundled"]
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 2,
                "r_material" => fill(material_resist_dict[code], 2)
            )
        elseif code == "uglv_120cu_xlpe/nyl/pvc_ug_4w_bundled"
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 4,
                "r_material" => fill(material_resist_dict[code], 4)
            )
        end
    end
    return mn_data
end


function build_linecode_for_ug_noshunt_eulvtf!(data, eng, z_pu)

    for (b, branch) in data["branch"]
        if b ∈ ["32", "2", "105", "109", "74", "41", "51", "27", "75", "42", "33", "28", "63", "93", "26", "10", "77", "59", "5", "89", "62", "43", "90", "39"]
            # two-wire bits of one linecode
            r = eng["linecode"]["ugsc_16cu_xlpe/nyl/pvc_ug_2w_bundled"]["rs"]
            x = eng["linecode"]["ugsc_16cu_xlpe/nyl/pvc_ug_2w_bundled"]["xs"]
            eng["line"][branch["name"]]["linecode"] = "ugsc_16cu_xlpe/nyl/pvc_ug_2w_bundled"
        else
            if length(branch["f_connections"]) == 2 # two-wire bits of other linecode
                r = eng["linecode"]["ugsc_25cu_xlpe/nyl/pvc_ug_2w_bundled"]["rs"]
                x = eng["linecode"]["ugsc_25cu_xlpe/nyl/pvc_ug_2w_bundled"]["xs"]    
                eng["line"][branch["name"]]["linecode"] = "ugsc_25cu_xlpe/nyl/pvc_ug_2w_bundled"
            else
                r = eng["linecode"]["uglv_120cu_xlpe/nyl/pvc_ug_4w_bundled"]["rs"]
                x = eng["linecode"]["uglv_120cu_xlpe/nyl/pvc_ug_4w_bundled"]["xs"]    
                # linecode name in `eng` for these ones is the default one        
            end
        end
        l = eng["line"][branch["name"]]["length"]
        branch["br_r"] = r.*l./z_pu
        branch["br_x"] = x.*l./z_pu
    end
    return data, eng
end

function build_linecode_for_ug_noshunt_30l!(data, eng, z_pu) # assigns the set of linecodes we elected for this case and builds R,X
    for (b, branch) in data["branch"]
        if b ∈ ["32", "1", "2", "51", "27", "33", "28", "25", "49", "5", "43", "34", "44", "55", "37", "12", "20", "6", "7", "57", "4", "22"]
            # two-wire bits of one linecode
            r = eng["linecode"]["uglv_185al_xlpe/nyl/pvc_ug_2w_bundled"]["rs"]
            x = eng["linecode"]["uglv_185al_xlpe/nyl/pvc_ug_2w_bundled"]["xs"]
            eng["line"][branch["name"]]["linecode"] = "uglv_185al_xlpe/nyl/pvc_ug_2w_bundled"
        else
            if length(branch["f_connections"]) == 2 # two-wire bits of other linecode
                r = eng["linecode"]["ugsc_16al_xlpe/pvc_ug_2w_bundled"]["rs"]
                x = eng["linecode"]["ugsc_16al_xlpe/pvc_ug_2w_bundled"]["xs"]    
                eng["line"][branch["name"]]["linecode"] = "ugsc_16al_xlpe/pvc_ug_2w_bundled"
            else
                r = eng["linecode"]["uglv_240al_xlpe/nyl/pvc_ug_4w_bundled"]["rs"]
                x = eng["linecode"]["uglv_240al_xlpe/nyl/pvc_ug_4w_bundled"]["xs"]    
                # linecode name in `eng` for these ones is the default one        
            end
        end
        l = eng["line"][branch["name"]]["length"]
        branch["br_r"] = r.*l./z_pu
        branch["br_x"] = x.*l./z_pu
    end
    return data,eng
end

function add_material_properties_for_ug_noshunt_30l!(mn_data, eng)
    material_resist_dict = Dict(
        "uglv_185al_xlpe/nyl/pvc_ug_2w_bundled" => 45.66553107812399,
        "ugsc_16al_xlpe/pvc_ug_2w_bundled" => 44.15319979061239,
        "uglv_240al_xlpe/nyl/pvc_ug_4w_bundled" => 43.076513156374084
    )

    mn_data["nw"]["1"]["linecode_map"] = Dict{Int, Any}() 
    for (id, code) in enumerate(keys(eng["linecode"]))
        if code ∈ ["uglv_185al_xlpe/nyl/pvc_ug_2w_bundled", "ugsc_16al_xlpe/pvc_ug_2w_bundled"]
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 2,
                "r_material" => fill(material_resist_dict[code], 2)
            )
        elseif code == "uglv_240al_xlpe/nyl/pvc_ug_4w_bundled"
            mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
                "name" => code,
                "n_wires" => 4,
                "r_material" => fill(material_resist_dict[code], 4)
            )
        end
    end

    for (b,bus) in mn_data["nw"]["1"]["bus"]
        bus["imp_grounded"] = fill(false, length(bus["terminals"]))
    end  

    return mn_data
end

function build_linecode_for_oh_ground_30l!(data, eng) # assigns the set of linecodes we elected for this case and builds R,X
    for (b, branch) in data["branch"]
        if length(branch["f_connections"]) == 4 && b ∈ ["29", "54", "41", "53", "42", "50", "52", "26", "10", "24", "23", "31", "39", "17", "47", "9"]
            # 4-wire branches of a linecode which is not the default one
            r = eng["linecode"]["hydrogen"]["rs"]
            x = eng["linecode"]["hydrogen"]["xs"]
            eng["line"][branch["name"]]["linecode"] = "hydrogen"
        end
        if b ∈ ["32", "1", "2", "51", "27", "33", "28", "25", "49", "5", "43", "34", "44", "55", "37", "12", "20", "6", "7", "57", "4", "22"]
            # two-wire bits of one linecode
            r = eng["linecode"]["abc2x16_lv_oh_2w_bundled"]["rs"]
            x = eng["linecode"]["abc2x16_lv_oh_2w_bundled"]["xs"]
            eng["line"][branch["name"]]["linecode"] = "abc2x16_lv_oh_2w_bundled"
        else
            if length(branch["f_connections"]) == 2 # two-wire bits of other linecode
                r = eng["linecode"]["tw2x16_lv_oh_2w_bundled"]["rs"]
                x = eng["linecode"]["tw2x16_lv_oh_2w_bundled"]["xs"]    
                eng["line"][branch["name"]]["linecode"] = "tw2x16_lv_oh_2w_bundled"
            else # default four-wire linecode        
                r = eng["linecode"]["pluto"]["rs"]
                x = eng["linecode"]["pluto"]["xs"]     
            end
        end
        l = eng["line"][branch["name"]]["length"]
        branch["br_r"] = r.*l./z_pu
        branch["br_x"] = x.*l./z_pu
    end
end

function prepare_math_eng_data(profiles;feeder_name::String="30load-feeder", oh_or_ug::String="ug")
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
            bus["vm_pair_lb"] = [(1, 4, 0.8);(2, 4, 0.8);(3, 4, 0.8)]
            bus["vm_pair_ub"] = [(1, 4, 1.2);(2, 4, 1.2);(3, 4, 1.2)]
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

function find_most_loaded_timesteps(profiles::_DF.DataFrame, nr_timesteps::Int)
    sorted_df = sort(
        _DF.DataFrame(idx = 1:size(profiles)[1], val = sum(eachcol(profiles))),
        [:val, _DF.order(:idx)], rev=true
        )
    return sorted_df.idx[1:nr_timesteps]
end