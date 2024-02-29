"""
 Random sample of length value from uniform distribution
 between the original (exact) length and plus/minus (pm) some percent
 TODO: seed + customize distribution?
"""
sample_start_length_uniform(lgth::Float64; pm::Float64 = 0.1) = 
     _RAN.rand(_DST.Uniform(lgth*(1-pm), lgth*(1+pm)))

"""
TBD: do we add length as measurement / likelihood for every time step or just one????
     do we do a gaussian description even though the sample is from a uniform dist?
This function adds the length guesses to the measurements and puts them in the objective
    through residuals. 
"""
function add_length_likelihood!(mn_data::Dict, branch::Dict; pm::Float64 = 0.1)
    m_id = maximum(parse.(Int, keys(mn_data["nw"]["1"]["meas"])))
    mn_data["nw"]["1"]["meas"]["$(m_id+1)"] = Dict{String, Any}(
        "var" => :l,
        "dst" => _DST.Normal(branch["l_start"], pm/3),
        "cmp" => :branch,
        "cmp_id" => branch["index"]
    )
end