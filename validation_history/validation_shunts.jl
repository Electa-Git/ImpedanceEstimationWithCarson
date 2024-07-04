import DataDrivenImpedanceEstimationWithCarson as _IMP
import PowerModelsDistribution as _PMD
import Ipopt

include("validation_utils.jl")

ie_solver = _PMD.optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 100., "max_iter" => 2000)

# get data
mn_data = set_up_validation_shunt(;power_mult = 15.)

# fix series impedance
for (b, branch) in mn_data["nw"]["1"]["branch"]
    branch["untrustworthy_branch"] = false
end

load_buses = [ load["load_bus"] for (_, load) in mn_data["nw"]["1"]["load"] ]

# add imperfect grounding flag to one load buses
for (b, bus) in mn_data["nw"]["1"]["bus"]
    if bus["index"] != 4
        bus["imp_grounded"] = fill(false, length(bus["terminals"]))
    else
        bus["imp_grounded"] = vcat(fill(false, length(bus["terminals"])-1), true)
    end
end    

mn_data["nw"]["1"]["linecode_map"] = Dict{Int, Any}() # needed to initialize problem

# for (n,nw) in mn_data["nw"]
#     for (m,meas) in nw["meas"]
#         if meas["var"] == :vm
#             delete!(mn_data["nw"][n]["meas"], m)
#         end
#     end
# end

sol = _IMP.solve_imp_est_carson(mn_data, ie_solver)

mn_data_free = deepcopy(mn_data)

for (b, bus) in mn_data_free["nw"]["1"]["bus"]
    bus["imp_grounded"] = fill(false, length(bus["terminals"]))
end   

sol_free = _IMP.solve_imp_est_carson(mn_data_free, ie_solver)

for (n, nw) in sol["solution"]["nw"]
    for (b, bus) in nw["bus"]
        bus["vm"] = sqrt.(bus["vi"].^2+bus["vr"].^2)
        bus["vmn"] = bus["vm"][1:(end-1)].-bus["vm"][end]
    end
end

for (n, nw) in sol_free["solution"]["nw"]
    for (b, bus) in nw["bus"]
        bus["vm"] = sqrt.(bus["vi"].^2+bus["vr"].^2)
        bus["vmn"] = bus["vm"][1:(end-1)].-bus["vm"][end]
    end
end

sol["solution"]["nw"]["1"]["bus"]
sol_free["solution"]["nw"]["1"]["bus"]