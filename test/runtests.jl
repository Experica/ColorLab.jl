using ColorLab,LinearAlgebra,Test,YAML

@testset "ColorLab" begin
    
    include("transformtest.jl")
    include("colortest.jl")

end