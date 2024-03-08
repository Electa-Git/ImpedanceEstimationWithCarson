import DataDrivenImpedanceEstimationWithCarson as _IMP
import CSV
import DataFrames as _DF
import PowerModelsDistribution as _PMD
import Ipopt

include("utils.jl")

##### things to customize (conf file???)
# 1) solver settings
# 2) smart meter noise??
# 3) number (and subset!) of timesteps
# 4) length Y/N and bounds
# 5) power multiplier?
# 6) different samples on the errors?

ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "max_iter" => 3000)
profiles = CSV.read(_IMP.DATA_DIR*"/nrel_profiles.csv", _DF.DataFrame, ntasks = 1)

data, eng, z_pu = prepare_math_eng_data()

###################################
### CHANGE LINECODES OF SERVICE CABLES TO 2-WIRE (EVERYHING IS 4-WIRE IN THE BEGINNING BY CONSTRUCTION)
###################################

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
        end
    end
    l = eng["line"][branch["name"]]["length"]
    branch["br_r"] = r.*l./z_pu
    branch["br_x"] = x.*l./z_pu
end

############### CREATE MULTINETWORK DATA WITH MEASUREMENT TIMESERIES ###############
# it runs a powerflow for each time step first, so it takes some time...

pf_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "print_level"=>0 )
mn_data = _IMP.build_multinetwork_dsse_data(data, profiles_df, pf_solver, σ_v, σ_d, σ_g; t_start=t_start, t_end=t_end, add_noise=false, power_mult = power_mult)

mn_data["nw"]["1"]["linecode_map"] = Dict{Int, Any}() 
for (id, code) in enumerate(keys(eng["linecode"]))
    if code == "lc6"
        mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
            "name" => code,
            "n_wires" => 2
        )
    else
        mn_data["nw"]["1"]["linecode_map"][id] = Dict{String, Any}(
            "name" => code,
            "n_wires" => 4
        )
    end
end

make_all_branches_untrustworthy!(mn_data, eng)

for (b,bus) in mn_data["nw"]["1"]["bus"]
    bus["imp_grounded"] = fill(false, length(bus["terminals"]))
end  