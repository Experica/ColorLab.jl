using Dierckx,LinearAlgebra

"Tristimulus matching values of single wavelength unit spectral"
matchlambda(λ::Real,cmf::AbstractMatrix) = @views [Spline1D(cmf[:,1],cmf[:,i+1],k=3,bc="zero")(λ) for i in 1:3]

"Tristimulus values, based on spectral measurement and matching functions"
function matchcolors(C,λ,I,cmf::AbstractMatrix;f=nothing)
    MC = similar(C)
    @views w = cmf[:,1]
    if isnothing(f)
        s = 1
    else
        Δw = w[2] - w[1]
        s = f*Δw
    end
    for i in eachindex(C)
        MI = Spline1D(λ[i],I[i],k=3,bc="zero")(w)'
        @views MC[i] = s*[MI*cmf[:,2], MI*cmf[:,3], MI*cmf[:,4]]
    end
    return MC,C
end

"Cone activations of colors, based on spectral measurement and cone fundamentals"
matchcolors_LMS(C,λ,I;observer=10,f=nothing) = matchcolors(C,λ,I,observer==10 ? sscone10le : sscone2le;f)

"""
CIE “physiologically-relevant” XYZ of colors, based on following formula:

``XYZ = 683∫S(λ)x̄ȳz̄(λ)dλ``

where S(λ) is the power spectrum, x̄ȳz̄(λ) are the matching functions and Y will be the luminance(cd/m²).
"""
function matchcolors_XYZ(C,λ,I;observer=10,f=683)
    if observer == 10
        conef = sscone10le
        M = LMSToXYZ10
    else
        conef = sscone2le
        M = LMSToXYZ2
    end
    # Only 1nm CMFs have strictly equal integrals, so that equal energy spectrum has [1/3, 1/3] CIE xy Chromaticity.
    matchcolors(C,λ,I,cmf(M,conef;wi=1:10:size(conef,1));f)
end

"new color matching functions by transform existing ones"
cmf(M::AbstractMatrix,f::AbstractMatrix;wi=:) = @views hcat(f[wi,1],(M*f[wi,2:end]')')
"Michelson Contrast, where ``michelson(Lmax, Lmin) = weber(Lmax, Lmean), Lmean = (Lmax+Lmin)/2``"
contrast_michelson(Lmax,Lmin) = (Lmax-Lmin)/(Lmax+Lmin)
"Weber Contrast, where ``weber(L, Lb) = michelson(L, 2Lb-L)``"
contrast_weber(L,Lb) = (L/Lb)-1
"Converting Matrices between RGB and LMS color spaces, based on spectral measurement and cone fundamentals"
RGBLMSMatrix(C,λ,I;observer=10,f=nothing,ishomo=true) = ConvertMatrix(matchcolors_LMS(C,λ,I;observer,f)...;ishomo)
"Converting Matrices between RGB and CIE XYZ color spaces, based on spectral measurement and CIE xyz matching functions"
RGBXYZMatrix(C,λ,I;observer=10,f=683,ishomo=true) = ConvertMatrix(matchcolors_XYZ(C,λ,I;observer,f)...;ishomo)

"Convert CIE XYZ to CIE xyY."
function XYZ2xyY(XYZ::AbstractMatrix)
    size(XYZ,1) == 4 && (XYZ=dehomovector(XYZ))
    xyz = divsum(XYZ)
    @views xyz[3,:]=XYZ[2,:]
    replace!(xyz,NaN=>0)
    return xyz
end
"Convert CIE xyY to CIE XYZ."
function xyY2XYZ(xyY::AbstractMatrix;ishomo=true)
    @views s = xyY[3:3,:]./xyY[2:2,:]
    @views XYZ = s.*vcat(xyY[1:2,:], 1 .- sum(xyY[1:2,:],dims=1))
    replace!(XYZ,NaN=>0)
    ishomo && (XYZ = homovector(XYZ))
    return XYZ
end

"""
Converting Matrices between LMS and Cone Contrast(Weber) color spaces.
(DH Brainard, Cone contrast and opponent modulation color spaces, human color vision, 1996)
"""
function LMSContrastMatrix(bg::AbstractVector)
    # translate origin to bg to get differential cone activation
    t = TranslateMatrix(-bg)
    # scale relative to bg to get cone contrast
    s = ScaleMatrix(inv.(bg))
    LMSToContrast = s*t
    ContrastToLMS = inv(LMSToContrast)
    return LMSToContrast,ContrastToLMS
end

"""
Converting Matrices between differential LMS relative to background and DKL[L+M, L-M, S-(L+M)] color spaces.
(DH Brainard, Cone contrast and opponent modulation color spaces, human color vision, 1996)
"""
function dLMSDKLMatrix(bg::AbstractVector;cone=nothing,v=nothing,isnorm=true,ishomo=true)
    wl = 1;wm = 1
    if !isnothing(cone) && !isnothing(v)
        @views wl,wm = v[:,2]'/cone[:,2:3]'
    end
    dLMSToDKL = [wl        wm                 0.0;
                 1.0  -bg[1]/bg[2]            0.0;
                -wl       -wm     (wl*bg[1] + wm*bg[2])/bg[3]]
    if isnorm
        # Each column of the inverse of dLMSToDKL is the differential LMS relative to background that isolating each DKL mechanism
        dlms_dkliso = inv(dLMSToDKL)
        # Cone contrast relative to bg
        cc = dlms_dkliso./bg[1:3]
        # Pooled cone contrast of each differential LMS relative to background that isolating each DKL mechanism
        pcc = [norm(i) for i in eachslice(cc,dims=2)]
        # Scale differential LMS by its pooled cone contrast
        udlms_dkliso = dlms_dkliso./pcc'
        # Rescale dLMSToDKL so that differential LMS which isolating DKL mechanism and having unit pooled cone contrast will result unit DKL
        dLMSToDKL = inv(dLMSToDKL*udlms_dkliso)*dLMSToDKL
    end
    if ishomo
        dLMSToDKL = homomatrix(dLMSToDKL)
    end
    DKLTodLMS = inv(dLMSToDKL)
    return dLMSToDKL,DKLTodLMS
end

"Converting Matrices between LMS and DKL[L+M, L-M, S-(L+M)] color spaces."
function LMSDKLMatrix(bg::AbstractVector;observer=10,isnorm=true)
    if observer == 10
        cone = sscone10le
        v = v10le
    else
        cone = sscone2le
        v = v2le
    end
    dLMSToDKL,DKLTodLMS = dLMSDKLMatrix(bg;cone,v,isnorm,ishomo=true)
    LMSToDKL = dLMSToDKL*TranslateMatrix(-bg)
    DKLToLMS = inv(LMSToDKL)
    return LMSToDKL,DKLToLMS
end

"Desaturate perceptible but not displayable CIE colors(each column) into the gamut of a display"
function desaturate2gamut!(x::AbstractMatrix)
    for i in eachslice(x,dims=2)
        min = minimum(i)
        if min<0
            i.-=min
        end
        max = maximum(i)
        if max>0
            i./=max
        end
    end
    return x
end


function cam16uniquehueindex(h′)
    h′ᵢ = zeros(Int,size(h′))
    for j in 1:length(h′)
        for i in 1:4
            if cam16uniquehue[1,i] <= h′[j] < cam16uniquehue[1,i+1]
                h′ᵢ[j]=i
                break
            end
        end
    end
    return h′ᵢ
end
function cam16uniquehuecomposition(H,i)
    Pl = round(Int,cam16uniquehue[3,i+1]-H)
    Pr = round(Int,H-cam16uniquehue[3,i])
    return (cam16uniquehue[4,i]=>Pl,cam16uniquehue[4,i+1]=>Pr)
end
function cam16nonlinearcompression(x,Fl;checksign=true)
    if checksign
        p = (Fl*abs(x)/100.0)^0.42
        sign(x)*400*p/(27.13 + p) + 0.1
    else
        p = (Fl*x/100.0)^0.42
        400*p/(27.13 + p) + 0.1
    end
end
"""
CAM16 Viewing Conditions

W: White in test illuminant [Xw, Yw, Zw]
Yb: Background in test conditions
La: Luminance of test adapting field (cd/m2)
Surround: Surround condition {Average, Dim, Dark}, Nc and F are functions of c, and their values can be linearly interpolated.
          Sr = Lsw / Ldw, where Lsw is the luminance of reference white in surround and Ldw in the display area.
          Sr == 0: Dark
          0 < Sr < 0.2: Dim
          Sr >= 0.2: Average
"""
function cam16view(;W=100*WP_D65,Surround=:Average,La=40,Yb=20)
    F,c,Nc = cam16surround[cam16surround[:,1].==Surround,2:end]
    k = 1.0/(5*La + 1)
    Fl = La*(k^4) + 0.1 * ((1 - k^4)^2) * ((5*La)^(1/3))

    Yw = W[2]
    n = Yb/Yw
    Nbb = Ncb = 0.725*(1/n)^0.2
    z = 1.48 + sqrt(n)
    D = clamp(F * (1 - (1/3.6)*exp(-(La+42)/92)), 0, 1)

    RGBw = cat16*W
    Da = Diagonal([D*Yw/RGBw[i] + 1 - D for i in 1:3])
    RGBwc = Da*RGBw
    RGBaw = cam16nonlinearcompression.(RGBwc,Fl,checksign=false)
    Aw = Nbb*(([2 1 1/20]*RGBaw)[1] - 0.305)

    return (c=c,Nc=Nc,Nbb=Nbb,Ncb=Ncb,Fl=Fl,n=n,z=z,Da=Da,Aw=Aw)
end
"""
Convert CIE XYZ to CIE Color Appearance Model 2016
Li, C., Li, Z., Wang, Z., Xu, Y., Luo, M.R., Cui, G., Melgosa, M., Brill, M.H., and Pointer, M. (2017). Comprehensive color solutions: CAM16, CAT16, and CAM16‐UCS.

return:
J is the lightness
C is the chroma
h is the hue angle
Q is the brightness
M is the colourfulness
s is the saturation
"""
function XYZ2CAM16(XYZ;Fl,Nc,Ncb,Nbb,n,c,z,Da,Aw)
    # Chromatic Adaptation
    RGBc = Da*cat16*XYZ
    # Non Linear Response Compression
    RGBa = cam16nonlinearcompression.(RGBc,Fl,checksign=true)
    # Perceptual Attribute Correlates
    a = ([1 -12/11 1/11]*RGBa)[:]
    b = ([1/9 1/9 -2/9]*RGBa)[:]
    h = mod.(2pi .+ atan.(b,a), 2pi)

    h′ = rad2deg.(h)
    i = h′ .< cam16uniquehue[1,1]
    h′[i] = h′[i] .+ 360

    h′ᵢ = cam16uniquehueindex(h′)
    H = cam16uniquehue[3,h′ᵢ] .+ (100*(h′.-cam16uniquehue[1,h′ᵢ])./cam16uniquehue[2,h′ᵢ]) ./ ( (h′.-cam16uniquehue[1,h′ᵢ])./cam16uniquehue[2,h′ᵢ] .+ (cam16uniquehue[1,h′ᵢ.+1] .- h′)./cam16uniquehue[2,h′ᵢ.+1] )
    Hc = cam16uniquehuecomposition.(H,h′ᵢ)

    et = 0.25 .* (cos.(deg2rad.(h′) .+ 2) .+ 3.8)
    t = (50000/13)*Nc*Ncb .* et .* sqrt.(a.^2 .+ b.^2) ./ (([1 1 21/20]*RGBa)[:])
    A = (Nbb.*([2 1 1/20]*RGBa .- 0.305))[:]

    J = 100 .* (A./Aw).^(c*z)
    Q = (4/c) .* sqrt.(J./100) .* (Aw + 4)*Fl^0.25
    C = t.^0.9 .* sqrt.(J./100) .* (1.64 - 0.29^n)^0.73
    M = (Fl^0.25).*C
    s = 100 .* sqrt.(M./Q)

    return (J=J,Q=Q,C=C,M=M,h=h,s=s,H=H,Hc=Hc)
end
"""
CAM16 Uniform Color Space in polar[J′, M′, h] or cartesian[J′, a′, b′]
Luo, M.R., Cui, G., and Li, C. (2006). Uniform colour spaces based on CIECAM02 colour appearance model. Color Research & Application 31, 320–330.
"""
function CAM16UCS(J,M,h;c1=0.007,c2=0.0228,form=:cartesian)
    J′ = (1 + 100c1).*J ./ (1 .+ c1.*J)
    M′ = log.(1 .+ c2.*M) ./ c2
    a′ = M′.*cos.(h)
    b′ = M′.*sin.(h)
    return form == :cartesian ? hcat(J′,a′,b′)' : hcat(J′,M′,h)'
end
function CAM16UCSinv(UCS;c1=0.007,c2=0.0228,form=:cartesian)
    J′ = UCS[1,:]
    if form == :cartesian
        M′ = sqrt.(UCS[2,:].^2 .+ UCS[3,:].^2)
        h = mod.(2pi .+ atan.(UCS[3,:],UCS[2,:]), 2pi)
    else
        M′ = UCS[2,:]
        h = UCS[3,:]
    end
    M = (exp.(c2.*M′) .- 1) ./ c2
    J = J′ ./ (1 .+ 100c1 .- c1.*J′)
    return (J=J,M=M,h=h)
end
function cam16nonlinearcompressioninv(x,Fl)
    p = x-0.1
    sign(p)*(100/Fl)*(27.13*abs(p) / (400-abs(p)))^(1/0.42)
end
"""
Convert CIE Color Appearance Model 2016 to CIE XYZ
"""
function CAM162XYZ(J,M,h;Fl,Nc,Ncb,Nbb,n,c,z,Da,Aw)
    C = M ./ (Fl^0.25)

    t = (C ./ (sqrt.(J./100) .* (1.64 - 0.29^n)^0.73)).^(1/0.9)
    et = 0.25 .* (cos.(h .+ 2) .+ 3.8)
    A = Aw.*(J./100).^(1/(c*z))
    p2 = A./Nbb .+ 0.305
    p3 = 21/20
    sh = sin.(h)
    ch = cos.(h)

    a=b=zeros(length(h))
    for i in 1:length(h)
        if t[i] !=0
            p1 = 50000/13 * Nc*Ncb * et[i] * (1 / t[i])
            if abs(sh[i]) >= abs(ch[i])
                p4 = p1/sh[i]
                b[i] = p2[i]*(2+p3)*(460/1403) / ( p4+(2+p3)*(220/1403)*(ch[i]/sh[i]) - (27/1403) + p3*(6300/1403) )
                a[i] = b[i]*ch[i]/sh[i]
            else
                p5 = p1/ch[i]
                a[i] = p2[i]*(2+p3)*(460/1403) / ( p5+(2+p3)*(220/1403) - (27/1403 - p3*6300/1403)*sh[i]/ch[i] )
                b[i] = a[i]*sh[i]/ch[i]
            end
        end
    end

    RGBa = [460/1403 451/1403 288/1403;
            460/1403 -891/1403 -261/1403;
            460/1403 -220/1403 -6300/1403]*([p2 a b]')
    RGBc = cam16nonlinearcompressioninv.(RGBa,Fl)
    XYZ = cat16\Da\RGBc
end
