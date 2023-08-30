"De-Augmenting Homogeneous Coordinate to Cartesian Coordinate"
dehomovector(x::AbstractVector) = x[1:end-1]
"De-Augmenting Homogeneous Coordinates(each column) to Cartesian Coordinates"
dehomovector(x::AbstractMatrix) = x[1:end-1,:]
"De-Augmenting Homogeneous Transformation Matrix to Linear Transformation Matrix"
dehomomatrix(x::AbstractMatrix) = x[1:end-1,1:end-1]
"Augmenting Cartesian Coordinate to Homogeneous Coordinate"
homovector(x::AbstractVector) = [x;one(eltype(x))]
"Augmenting Cartesian Coordinates(each column) to Homogeneous Coordinates"
function homovector(x::AbstractMatrix)
    hx = ones(eltype(x),size(x).+(1,0))
    hx[1:end-1,:] = x
    hx
end
"Augmenting Linear Transformation Matrix to Homogeneous Transformation Matrix"
function homomatrix(x::AbstractMatrix)
    hx = zeros(eltype(x),size(x).+1)
    hx[1:end-1,1:end-1] = x
    hx[end,end] = one(eltype(x))
    hx
end

ConvertMatrix(y::AbstractVector,x::AbstractVector;di=1:3,ishomo=true) = @views ConvertMatrix(stack(i->i[di],y),stack(i->i[di],x);ishomo)
"Solve System of Linear Equations to get converting matrices between X and Y"
function ConvertMatrix(y::AbstractMatrix,x::AbstractMatrix;ishomo=true)
    X2Y = y/x
    if ishomo
        X2Y = homomatrix(X2Y)
    end
    Y2X = inv(X2Y)
    return X2Y,Y2X
end

TranslateMatrix(t::AbstractVector) = TranslateMatrix(t[1],t[2],t[3])
TranslateMatrix(;x::Real=0,y::Real=0,z::Real=0) = TranslateMatrix(x,y,z)
function TranslateMatrix(x::Real,y::Real,z::Real)
    tm = [1.0  0.0  0.0  x;
          0.0  1.0  0.0  y;
          0.0  0.0  1.0  z;
          0.0  0.0  0.0  1.0]
end
ScaleMatrix(s::AbstractVector) = ScaleMatrix(s[1],s[2],s[3])
ScaleMatrix(;x::Real=1,y::Real=1,z::Real=1) = ScaleMatrix(x,y,z)
function ScaleMatrix(x::Real,y::Real,z::Real)
    sm = [x     0.0   0.0   0.0;
          0.0   y     0.0   0.0;
          0.0   0.0   z     0.0;
          0.0   0.0   0.0   1.0]
end
function RotateXMatrix(θ::Real)
    sinθ,cosθ = sincos(θ)
    rm = [1.0   0.0   0.0   0.0;
          0.0  cosθ  -sinθ  0.0;
          0.0  sinθ   cosθ  0.0;
          0.0   0.0   0.0   1.0]
end
function RotateYMatrix(θ::Real)
    sinθ,cosθ = sincos(θ)
    rm = [cosθ   0.0   sinθ   0.0;
          0.0    1.0   0.0    0.0;
         -sinθ   0.0   cosθ   0.0;
          0.0    0.0   0.0    1.0]
end
function RotateZMatrix(θ::Real)
    sinθ,cosθ = sincos(θ)
    rm = [cosθ   -sinθ   0.0   0.0;
          sinθ    cosθ   0.0   0.0;
          0.0     0.0    1.0   0.0;
          0.0     0.0    0.0   1.0]
end
function RotateMatrix(θ::Real,dir::AbstractVector)
    sinθ,cosθ = sincos(θ)
    udir = normalize(dir)
    x=udir[1];y=udir[2];z=udir[3]
    rm = [x*x*(1-cosθ)+cosθ    y*x*(1-cosθ)-z*sinθ   z*x*(1-cosθ)+y*sinθ   0.0;
          x*y*(1-cosθ)+z*sinθ  y*y*(1-cosθ)+cosθ     z*y*(1-cosθ)-x*sinθ   0.0;
          x*z*(1-cosθ)-y*sinθ  y*z*(1-cosθ)+x*sinθ   z*z*(1-cosθ)+cosθ     0.0;
          0.0                       0.0                     0.0            1.0]
end

"Convert Absolute Coordinates[X,Y,Z] to Relative Coordinates[x,y,z] where ``x=X/X+Y+Z, y=Y/X+Y+Z, z=Z/X+Y+Z``"
divsum(x::AbstractVecOrMat) = x./sum(x,dims=1)

"""
Intersection point of a line and a plane.
points of a line are defined as a direction(Dₗ) through a point(Pₗ): P = Pₗ + λDₗ , where λ is a scaler.
points of a plane are defined as a plane through a point(Pₚ) and with normal vector(Nₚ) : Nₚᵀ(P - Pₚ) = 0 , where Nᵀ is the transpose of N.

return the point of intersection and if it's on direction.
"""
function intersectlineplane(Pₗ,Dₗ,Pₚ,Nₚ)
    NₚᵀDₗ = Nₚ'*Dₗ
    NₚᵀDₗ == 0 && return nothing,nothing
    λ = Nₚ'*(Pₚ - Pₗ) / NₚᵀDₗ
    return Pₗ + λ*Dₗ, λ >=0
end
"""
Intersection point of a line and the six faces of the unit cube with origin as a vertex and three axies as edges[0:1,0:1,0:1].
points of the line are defined as a direction(Dₗ) through a point(Pₗ).

return the intersection point on direction.
"""
function intersectlineunitcube(Pₗ,Dₗ)
    ps = [zeros(3,3) ones(3,3)]
    ns = [Matrix{Float64}(I,3,3) Matrix{Float64}(I,3,3)]
    for i in 1:6
        p,isdir = intersectlineplane(Pₗ,Dₗ,ps[:,i],ns[:,i])
        if !isnothing(p) && isdir && all(i -> -eps(1.0)<= i <=1+eps(1.0),p)
            return clamp!(p,0.0,1.0)
        end
    end
    return nothing
end
"""
Points of a line segment defined by two points P₀ and P₁.

- d: points density of unit line length
"""
function linepoints(P₀,P₁;d=10)
    d = max(nextfloat(0.0),d)
    Dₗ = P₁-P₀
    l =norm(Dₗ)
    stack(λ -> P₀.+λ*Dₗ , 0:1/(d*l):1)
end