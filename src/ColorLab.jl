module ColorLab

include("color_data.jl")
include("color_algorithm.jl")

# exportall
# for n in names(current_module(), all=true)
#     if Base.isidentifier(n) && n ∉ (Symbol(current_module()), :eval)
#         @eval export $n
#     end
# end

end # module
