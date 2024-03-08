module DataDrivenImpedanceEstimationWithCarson

    import DataFrames as _DF
    import Distributions as _DST
    import Graphs, SimpleWeightedGraphs
    import InfrastructureModels as _IM
    import JuMP
    import PowerModelsDistribution as _PMD

    const DATA_DIR = joinpath(dirname(@__DIR__), "data")
    const NTW_DATA_DIR = joinpath(DATA_DIR, "network_data")
    const LINECODE_DATA_DIR = joinpath(DATA_DIR, "linecode_library")

    include("core/constraint.jl")
    include("core/objective.jl")
    include("core/utils.jl")
    include("core/variable_carson.jl")
    include("core/variable_shared.jl")

    include("prob/imp_est_carson.jl")

    include("io/length_utils.jl")
    include("io/multinetwork_data.jl")
    include("io/parse_clean_network.jl")
    include("io/reduce_eng_model.jl")
    include("io/solution.jl")

end