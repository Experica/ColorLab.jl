module ColorLab

include("color_data.jl")
include("transformation.jl")
include("color_algorithm.jl")

# export all symbols
for n in names(@__MODULE__, all=true)
    if Base.isidentifier(n) && n ∉ (nameof(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end
