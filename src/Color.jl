using Plots,YAML,Makie,StatsMakie
include("color_algorithm.jl")

"Get digital RGB color spectral measured from a specific display"
function RGBSpectral(measurement)
    C = map(i->parse.(Float64,split(i)),measurement["Color"])
    λ = measurement["WL"]
    I = measurement["Spectral"]
    return C,λ,I
end

# plot cone fundamentals
Plots.plot(sscone2le[:,1],sscone2le[:,2:end],linewidth=2,color=["Red" "Green" "Blue"],xlabel="Wavelength (nm)",ylabel="Sensitivity",label=["L" "M" "S"],title="Cone Fundamentals(2deg)")
foreach(i->savefig("Cone Fundamentals(2deg)$i"),[".png",".svg"])

Plots.plot(sscone10le[:,1],sscone10le[:,2:end],linewidth=2,color=["Red" "Green" "Blue"],xlabel="Wavelength (nm)",ylabel="Sensitivity",label=["L" "M" "S"],title="Cone Fundamentals(10deg)")
foreach(i->savefig("Cone Fundamentals(10deg)$i"),[".png",".svg"])

# plot rgb color matching funcionts
Plots.plot(sbrgb2[:,1],sbrgb2[:,2:end],linewidth=2,color=["Red" "Green" "Blue"],xlabel="Wavelength (nm)",ylabel="Tristimulus Value",label=["r" "g" "b"],title="Color Matching Functions_rgb(2deg)")
Plots.scatter!(sbrgb_primary',zeros(3)',markersize=8,markerstrokewidth=0,markercolor=["Red" "Green" "Blue"],legend=false)
foreach(i->savefig("Color Matching Functions_rgb(2deg)$i"),[".png",".svg"])

Plots.plot(sbrgb10[:,1],sbrgb10[:,2:end],linewidth=2,color=["Red" "Green" "Blue"],xlabel="Wavelength (nm)",ylabel="Tristimulus Value",label=["r" "g" "b"],title="Color Matching Functions_rgb(10deg)")
Plots.scatter!(sbrgb_primary',zeros(3)',markersize=8,markerstrokewidth=0,markercolor=["Red" "Green" "Blue"],legend=false)
foreach(i->savefig("Color Matching Functions_rgb(10deg)$i"),[".png",".svg"])

# plot xyz color matching funcionts
xyz2 = newcmf(sscone2le,LMSToXYZ2)
Plots.plot(xyz2[:,1],xyz2[:,2:end],linewidth=2,color=["Red" "Green" "Blue"],xlabel="Wavelength (nm)",ylabel="Tristimulus Value",label=["x" "y" "z"],title="Color Matching Functions_xyz(2deg)")
foreach(i->savefig("Color Matching Functions_xyz(2deg)$i"),[".png",".svg"])

xyz10 = newcmf(sscone2le,LMSToXYZ10)
Plots.plot(xyz10[:,1],xyz10[:,2:end],linewidth=2,color=["Red" "Green" "Blue"],xlabel="Wavelength (nm)",ylabel="Tristimulus Value",label=["x" "y" "z"],title="Color Matching Functions_xyz(10deg)")
foreach(i->savefig("Color Matching Functions_xyz(10deg)$i"),[".png",".svg"])



# ASUS ROG Swift PG279Q IPS LCD
displayname = "ROGPG279Q"
# ViewSonic VX3276mhd IPS LCD
displayname = "VX3276mhd"
# Sony Trinitron CRT
displayname = "Trinitron"

resultdir = "./$displayname"
mkpath(resultdir)

config = YAML.load_file("CommandConfig - Display Measurement.yaml")
RGBToLMS,LMSToRGB = RGBLMSMatrix(RGBSpectral(config["Display"][displayname]["SpectralMeasurement"])...)
RGBToXYZ,XYZToRGB = RGBXYZMatrix(RGBSpectral(config["Display"][displayname]["SpectralMeasurement"])...)


# unit cube digital colors
ur=0:0.01:1
uc=[[i,j,k] for i=ur,j=ur,k=ur][:]
ucm=hcat(uc...)
uqcm=quavectors(ucm)

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGB(i...) for i in uc],markersize=0.01,transparency=true)
record(s,"Unit Color Space.mp4",1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

# Transformation between specific display RGB and LMS Spaces
lms_rgb = RGBToLMS*uqcm
rgb_lms = LMSToRGB*uqcm

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(lms_rgb[1,:],lms_rgb[2,:],lms_rgb[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("R / L","G / M","B / S")
record(s,joinpath(resultdir,"Unit RGB To LMS Space.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(rgb_lms[1,:],rgb_lms[2,:],rgb_lms[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("R / L","G / M","B / S")
record(s,joinpath(resultdir,"Unit LMS To RGB Space.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

# Transformation between specific display RGB and CIE XYZ Spaces
XYZ_rgb = RGBToXYZ*uqcm
rgb_XYZ = XYZToRGB*uqcm

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(XYZ_rgb[1,:],XYZ_rgb[2,:],XYZ_rgb[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("R / X","G / Y","B / Z")
record(s,joinpath(resultdir,"Unit RGB To XYZ Space.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(rgb_XYZ[1,:],rgb_XYZ[2,:],rgb_XYZ[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("R / X","G / Y","B / Z")
record(s,joinpath(resultdir,"Unit XYZ To RGB Space.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

# Transformation between LMS and Cone Contrast
bg = [0.5,0.5,0.5]
LMSToContrast,ContrastToLMS = LMSContrastMatrix(bg)
cc_lms = LMSToContrast*uqcm
lms_cc = ContrastToLMS*uqcm

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(cc_lms[1,:],cc_lms[2,:],cc_lms[3,:],color=[RGBA(i...,0.1) for i in uc],markersize=0.01,transparency=true)

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,0.1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(lms_cc[1,:],lms_cc[2,:],lms_cc[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)

bg = RGBToLMS*[0.5,0.5,0.5,1]
LMSToContrast,ContrastToLMS = LMSContrastMatrix(bg)
cc_lms = LMSToContrast*lms_rgb

s=Makie.scatter(lms_rgb[1,:],lms_rgb[2,:],lms_rgb[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(cc_lms[1,:],cc_lms[2,:],cc_lms[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)


# Transformation between LMS and DKL Spaces
bg = [0.5,0.5,0.5]
LMSToDKL,DKLToLMS = LMSDKLMatrix(bg,isnorm=true)
dkl_lms = LMSToDKL*uqcm
lms_dkl = DKLToLMS*uqcm

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(dkl_lms[1,:],dkl_lms[2,:],dkl_lms[3,:],color=[RGBA(i...,0.1) for i in uc],markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L / L+M","M / L-M","S / S-(L+M)")
record(s,"Unit LMS To DKL$bg Space.mp4",1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(lms_dkl[1,:],lms_dkl[2,:],lms_dkl[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L / L+M","M / L-M","S / S-(L+M)")

bg = RGBToLMS*[0.5,0.5,0.5,1]
LMSToDKL,DKLToLMS = LMSDKLMatrix(bg,isnorm=true)
dkl_lms = LMSToDKL*lms_rgb

s=Makie.scatter(lms_rgb[1,:],lms_rgb[2,:],lms_rgb[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(dkl_lms[1,:],dkl_lms[2,:],dkl_lms[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L / L+M","M / L-M","S / S-(L+M)")
record(s,joinpath(resultdir,"LMS To DKL$bg Space.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end



# Cone Isolating Modulation Range
lms_gray= RGBToLMS*[0.5,0.5,0.5,1]

s=Makie.surface(StatsMakie.density,lms_rgb[2,:],lms_rgb[3,:],colormap=:reds)
lines!([0,0].+lms_gray[2],[0,0].+lms_gray[3],[0,1],linewidth=3,color=:gray)
s.center=false
s[Axis][:names,:axisnames]=("M","S","Liso Strength")

s=Makie.surface(StatsMakie.density,lms_rgb[1,:],lms_rgb[3,:],colormap=:greens)
lines!([0,0].+lms_gray[1],[0,0].+lms_gray[3],[0,1],linewidth=3,color=:gray)
s.center=false
s[Axis][:names,:axisnames]=("L","S","Miso Strength")

s=Makie.surface(StatsMakie.density,lms_rgb[1,:],lms_rgb[2,:],colormap=:blues)
lines!([0,0].+lms_gray[1],[0,0].+lms_gray[2],[0,1],linewidth=3,color=:gray)
s.center=false
s[Axis][:names,:axisnames]=("L","M","Siso Strength")

# Cone Isolating RGB through a color
th = [0.5,0.5,0.5]
# each column of LMSToRGB is the cone isolating RGB direction
ConeIsoRGBVec = trimatrix(LMSToRGB)
# Scale RGB Vector into Unit Cube
ConeIsoRGBVec./=maximum(abs.(ConeIsoRGBVec),dims=1)
# since through color is the center of RGB cube, the line intersects at two symmatric points on the faces of unit cube
minc = th.-0.5*ConeIsoRGBVec;maxc = th.+0.5*ConeIsoRGBVec
ConeIsoRGB =map((i,j)->collect(i.+j), minc,[0:0.001:1].*(maxc.-minc))

Liso = hcat(ConeIsoRGB[:,1]...)'
Lisoc = [RGB(Liso[:,i]...) for i in 1:size(Liso,2)]
Lisolms=RGBToLMS*quavectors(Liso)

Miso = hcat(ConeIsoRGB[:,2]...)'
Misoc = [RGB(Miso[:,i]...) for i in 1:size(Miso,2)]
Misolms=RGBToLMS*quavectors(Miso)

Siso = hcat(ConeIsoRGB[:,3]...)'
Sisoc = [RGB(Siso[:,i]...) for i in 1:size(Siso,2)]
Sisolms=RGBToLMS*quavectors(Siso)


s=Makie.scatter(lms_rgb[1,:],lms_rgb[2,:],lms_rgb[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(Lisolms[1,:],Lisolms[2,:],Lisolms[3,:],color=Lisoc,markersize=0.01,transparency=true)
Makie.scatter!(Misolms[1,:],Misolms[2,:],Misolms[3,:],color=Misoc,markersize=0.01,transparency=true)
Makie.scatter!(Sisolms[1,:],Sisolms[2,:],Sisolms[3,:],color=Sisoc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L","M","S")
record(s,joinpath(resultdir,"ConeIsolating_ThroughGray_LMSSpace.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(Liso[1,:],Liso[2,:],Liso[3,:],color=Lisoc,markersize=0.01,transparency=true)
Makie.scatter!(Miso[1,:],Miso[2,:],Miso[3,:],color=Misoc,markersize=0.01,transparency=true)
Makie.scatter!(Siso[1,:],Siso[2,:],Siso[3,:],color=Sisoc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("R","G","B")
record(s,joinpath(resultdir,"ConeIsolating_ThroughGray_RGBSpace.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end


# DKL Isolating RGB through a background color
bg = RGBToLMS*[th;1]
LMSToDKL,DKLToLMS = LMSDKLMatrix(bg,isnorm=true)
# convert each column of DKLToLMS which are the DKL isolating LMS to RGB direction
DKLIsoRGBVec = trimatrix(LMSToRGB*DKLToLMS)
# Scale RGB Vector into Unit Cube
DKLIsoRGBVec./=maximum(abs.(DKLIsoRGBVec),dims=1)
# since through color is the center of RGB cube, the line intersects at two symmatric points on the faces of unit cube
minc = th.-0.5*DKLIsoRGBVec;maxc = th.+0.5*DKLIsoRGBVec
DKLIsoRGB =map((i,j)->collect(i.+j), minc,[0:0.001:1].*(maxc.-minc))

Lumiso = hcat(DKLIsoRGB[:,1]...)'
Lumisoc = [RGB(Lumiso[:,i]...) for i in 1:size(Lumiso,2)]
Lumisolms=RGBToLMS*quavectors(Lumiso)
Lumisodkl = LMSToDKL*Lumisolms

LMiso = hcat(DKLIsoRGB[:,2]...)'
LMisoc = [RGB(LMiso[:,i]...) for i in 1:size(LMiso,2)]
LMisolms=RGBToLMS*quavectors(LMiso)
LMisodkl = LMSToDKL*LMisolms

SLMiso = hcat(DKLIsoRGB[:,3]...)'
SLMisoc = [RGB(SLMiso[:,i]...) for i in 1:size(SLMiso,2)]
SLMisolms=RGBToLMS*quavectors(SLMiso)
SLMisodkl = LMSToDKL*SLMisolms

s=Makie.scatter(dkl_lms[1,:],dkl_lms[2,:],dkl_lms[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(Lumisodkl[1,:],Lumisodkl[2,:],Lumisodkl[3,:],color=Lumisoc,markersize=0.01,transparency=true)
Makie.scatter!(LMisodkl[1,:],LMisodkl[2,:],LMisodkl[3,:],color=LMisoc,markersize=0.01,transparency=true)
Makie.scatter!(SLMisodkl[1,:],SLMisodkl[2,:],SLMisodkl[3,:],color=SLMisoc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L+M","L-M","S-(L+M)")
record(s,joinpath(resultdir,"DKLIsolating_ThroughGray_DKLSpace.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

s=Makie.scatter(lms_rgb[1,:],lms_rgb[2,:],lms_rgb[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(Lumisolms[1,:],Lumisolms[2,:],Lumisolms[3,:],color=Lumisoc,markersize=0.01,transparency=true)
Makie.scatter!(LMisolms[1,:],LMisolms[2,:],LMisolms[3,:],color=LMisoc,markersize=0.01,transparency=true)
Makie.scatter!(SLMisolms[1,:],SLMisolms[2,:],SLMisolms[3,:],color=SLMisoc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L","M","S")
record(s,joinpath(resultdir,"DKLIsolating_ThroughGray_LMSSpace.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(Lumiso[1,:],Lumiso[2,:],Lumiso[3,:],color=Lumisoc,markersize=0.01,transparency=true)
Makie.scatter!(LMiso[1,:],LMiso[2,:],LMiso[3,:],color=LMisoc,markersize=0.01,transparency=true)
Makie.scatter!(SLMiso[1,:],SLMiso[2,:],SLMiso[3,:],color=SLMisoc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("R","G","B")
record(s,joinpath(resultdir,"DKLIsolating_ThroughGray_RGBSpace.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end


# DKL Isoluminance plane
DKLToRGB = LMSToRGB*DKLToLMS
RGBToDKL = LMSToDKL*RGBToLMS
lum = 0
# rotate l-m direction around l+m axis within the Isoluminance plane
DKLIsoLumRGBVec = trivectors(hcat(map(i->DKLToRGB*RotateXYZMatrix(deg2rad(i),dims=1)*[0,1,0,0], 0:0.05:179.95)...))
# find Intersections of Isoluminance directions with faces of unit RGB cube
minmaxc = hcat(map(i->intersectlineunitorigincube((DKLToRGB*[lum,0,0,1])[1:3],DKLIsoLumRGBVec[:,i]),1:size(DKLIsoLumRGBVec,2))...)
minc=minmaxc[:,1:2:end];maxc=minmaxc[:,2:2:end]
DKLIsoLumRGB =map((i,j)->collect(i.+j), minc,[0:0.001:1].*(maxc.-minc))

IsoLum=hcat(map(i->hcat(DKLIsoLumRGB[:,i]...)',1:size(DKLIsoLumRGB,2))...)
IsoLumc = [RGB(IsoLum[:,i]...) for i in 1:size(IsoLum,2)]
IsoLumlms=RGBToLMS*quavectors(IsoLum)
IsoLumdkl = LMSToDKL*IsoLumlms

s=Makie.scatter(IsoLumdkl[2,:],IsoLumdkl[3,:],color=IsoLumc,markersize=0.01,scale_plot=false)
s[Axis][:names,:axisnames]=("L-M","S-(L+M)")
s=title(s,"DKL Isoluminance Plane through $lum")
s.center=false
Makie.save(joinpath(resultdir,"DKL Isoluminance Plane through $lum.png"),s)

# combine several Isoluminance planes
IsoLum = hcat(IsoLum,hcat(map(i->hcat(DKLIsoLumRGB[:,i]...)',1:size(DKLIsoLumRGB,2))...))

IsoLumc = [RGB(IsoLum[:,i]...) for i in 1:size(IsoLum,2)]
IsoLumlms=RGBToLMS*quavectors(IsoLum)
IsoLumdkl = LMSToDKL*IsoLumlms

s=Makie.scatter(dkl_lms[1,:],dkl_lms[2,:],dkl_lms[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(IsoLumdkl[1,:],IsoLumdkl[2,:],IsoLumdkl[3,:],color=IsoLumc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L+M","L-M","S-(L+M)")
record(s,joinpath(resultdir,"DKLIsoluminancePlane_DKLSpace.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

s=Makie.scatter(lms_rgb[1,:],lms_rgb[2,:],lms_rgb[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(IsoLumlms[1,:],IsoLumlms[2,:],IsoLumlms[3,:],color=IsoLumc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L","M","S")
record(s,joinpath(resultdir,"DKLIsoluminancePlane_LMSSpace.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(IsoLum[1,:],IsoLum[2,:],IsoLum[3,:],color=IsoLumc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("R","G","B")
record(s,joinpath(resultdir,"DKLIsoluminancePlane_RGBSpace.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end


# rgb color matching chromaticity
ploc = divsum(hcat([matchlambda(w,sbrgb10) for w in sbrgb_primary]...))
wloc = divsum(hcat([matchlambda(w,sbrgb10) for w in 390:830]...))

s=Makie.lines(wloc[1,:],wloc[2,:],wloc[3,:],color=:black,linewidth=3,transparency=true)
Makie.scatter!(ploc[1,:],ploc[2,:],ploc[3,:],color=[:red,:green,:blue],markersize=0.05,transparency=true)
Makie.lines!(wloc[1,:],wloc[2,:],color=:gray40,linewidth=3,transparency=true)
Makie.scatter!(ploc[1,:],ploc[2,:],color=[:red,:green,:blue],markersize=0.05,transparency=true)

aw =  [450,500,550,600]
awloc = divsum(hcat([matchlambda(w,sbrgb10) for w in aw]...))
Makie.scatter!(awloc[1,:],awloc[2,:],color=:gray20,markersize=0.03,transparency=true)
annotations!(string.(" ",aw),Point2.(awloc[1,:],awloc[2,:]),textsize=0.08)
s[Axis][:names,:axisnames]=("r","g","b")
record(s,"Color Matching Chromaticity_rgb.mp4",1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end


# xyz color matching chromaticity
d = 250
wloc = divsum(hcat([matchlambda(w,xyz10) for w in 390:0.15:830]...))
aw =  [450,500,550,600,650]
awloc = divsum(hcat([matchlambda(w,xyz10) for w in aw]...))
realxy = hcat([linepoints(wloc[1:2,1],wloc[1:2,i],d=d) for i in 2:size(wloc,2)]...)
realxyz = vcat(realxy,1 .- sum(realxy,dims=1))

xyz_rgb_primary = divsum(trimatrix(RGBToXYZ))
xyz_white = divsum(trimatrix(RGBToXYZ)*[1,1,1])
xyzToRGB = inv(xyz_rgb_primary)
xyzToRGB ./= xyzToRGB*xyz_white
rgb_realxyz = desaturate2gamut!(xyzToRGB*realxyz)
gamuttriangle = [xyz_rgb_primary xyz_rgb_primary[:,1:1]]

s = Makie.scatter(realxyz[1,:],realxyz[2,:],realxyz[3,:],color=[RGBA(rgb_realxyz[:,i]...,1) for i in 1:size(rgb_realxyz,2)],markersize=0.005,transparency=true)

s = Makie.scatter(realxy[1,:],realxy[2,:],color=[RGBA(rgb_realxyz[:,i]...,1) for i in 1:size(rgb_realxyz,2)],markersize=0.005,transparency=true,limits = FRect(0,0,0.8,0.8))
Makie.lines!(wloc[1,:],wloc[2,:],color=:black,linewidth=4)
Makie.scatter!(awloc[1,:],awloc[2,:],color=:black,markersize=0.01)
annotations!(string.(aw),[Point2(awloc[1,i]>0.3 ? awloc[1,i]+0.01 : awloc[1,i]-0.05,awloc[2,i]) for i in 1:length(aw)],textsize=0.02)
Makie.lines!(gamuttriangle[1,:],gamuttriangle[2,:],gamuttriangle[3,:],linewidth=0.5,scale_plot=false)
s=title(s,"RGB gamut in CIE xy chromaticity")
s.center=false
Makie.save(joinpath(resultdir,"RGB gamut in CIE xy chromaticity.png"),s)


# CIECAM16 Uniform Color Space
vc = cam16view(Surround=:Dark)
cam = XYZ2CAM16(100*trivectors(XYZ_rgb);vc...)
camucs = CAM16UCS(cam.J,cam.M,cam.h,form=:cartesian)

s=Makie.scatter(camucs[1,:],camucs[2,:],camucs[3,:],color=[RGBA(i...,1) for i in uc],markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("J′","a′","b′")
record(s,"CAM16 Uniform Color Space.mp4",1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end



s = Makie.scatter(camucs[1,:],camucs[2,:],camucs[3,:],color=[RGBA(rgb_realxyz[:,i]...,1) for i in 1:size(rgb_realxyz,2)],markersize=0.005,transparency=true)



cami = CAM16UCSinv(camucs,form=:cartesian)
t = CAM162XYZ(cami.J,cami.M,cami.h;vc...)

s=Makie.lines(StatsMakie.histogram(nbins=200),cam02.J,colormap=:reds)

lines(StatsMakie.histogram,camucs[3,:])


lines(StatsMakie.histogram,cami.h)

extrema(cam.h' .-cami.h)


t.-trivectors(100*XYZ_rgb)
sum(t.-trivectors(100*XYZ_rgb))


xyz_rgb = divsum(trivectors(XYZ_rgb))



s = Makie.scatter(realxyz[1,:],realxyz[2,:],realxyz[3,:],color=[RGBA(rgb_realxyz[:,i]...,1) for i in 1:size(rgb_realxyz,2)],markersize=0.005,transparency=true)

s.center=false
s[Axis][:names,:axisnames]=("R / L","G / M","B / S")
record(s,joinpath(resultdir,"Unit RGB To LMS Space.mp4"),1:360/5,framerate=12) do i
    rotate_cam!(s,5pi/180,0.0,0.0)
end
















tIsoLumdkl = TranslateXYZMatrix(x=-0.5)*IsoLumdkl
Makie.scatter!(tIsoLumdkl[1,:],tIsoLumdkl[2,:],tIsoLumdkl[3,:],color=IsoLumc,markersize=0.01,transparency=true)
tIsoLumrgb = trivectors(DKLToRGB*tIsoLumdkl)
tIsoLumrgb = clamp.(tIsoLumrgb,0,1)
tIsoLumrgb./=maximum(abs.(tIsoLumrgb),dims=1)
Makie.scatter!(tIsoLumrgb[1,:],tIsoLumrgb[2,:],tIsoLumrgb[3,:],color=IsoLumc,markersize=0.01,transparency=true)


DKLIsoLumRGBVec = trivectors(hcat([DKLToRGB*RotateXYZMatrix(deg2rad(i>90 ? -j : j),dims=3)*RotateXYZMatrix(deg2rad(i),dims=1)*[0,1,0,0] for i=0:1:180,j=[75]][:]...))

DKLIsoLumRGBVec./=maximum(abs.(DKLIsoLumRGBVec),dims=1)
minc = th.-0.5*DKLIsoLumRGBVec;maxc = th.+0.5*DKLIsoLumRGBVec
DKLIsoLumRGB =map((i,j)->collect(i.+j), minc,[0:0.01:1].*(maxc.-minc))

IsoLum=hcat(map(i->hcat(DKLIsoLumRGB[:,i]...)',1:size(DKLIsoLumRGB,2))...)
IsoLumc = [RGB(IsoLum[:,i]...) for i in 1:size(IsoLum,2)]
IsoLumlms=RGBToLMS*quavectors(IsoLum)
IsoLumdkl = LMSToDKL*IsoLumlms

s=Makie.scatter(dkl_lms[1,:],dkl_lms[2,:],dkl_lms[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(IsoLumdkl[1,:],IsoLumdkl[2,:],IsoLumdkl[3,:],color=IsoLumc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L+M","L-M","S-(L+M)")

s=Makie.scatter(lms_rgb[1,:],lms_rgb[2,:],lms_rgb[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(IsoLumlms[1,:],IsoLumlms[2,:],IsoLumlms[3,:],color=IsoLumc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("L","M","S")

s=Makie.scatter(ucm[1,:],ucm[2,:],ucm[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(IsoLum[1,:],IsoLum[2,:],IsoLum[3,:],color=IsoLumc,markersize=0.01,transparency=true)
s.center=false
s[Axis][:names,:axisnames]=("R","G","B")




t=Makie.scatter(lms_rgb[1,:],lms_rgb[2,:],lms_rgb[3,:],color=[RGBA(i...,0.01) for i in uc],markersize=0.01,transparency=true)
Makie.scatter!(IsoLumlms[1,:],IsoLumlms[2,:],IsoLumlms[3,:],color=IsoLumc,markersize=0.01,transparency=true)



s=Makie.scatter(minc[1,:],minc[2,:],minc[3,:],color=[RGBA(minc[:,i]...,1) for i in size(minc,2)],markersize=0.01,transparency=true)
Makie.scatter!(maxc[1,:],maxc[2,:],maxc[3,:],color=[RGBA(maxc[:,i]...,1) for i in size(maxc,2)],markersize=0.01,transparency=false)




s=Makie.scatter(minc[1,:],minc[2,:],minc[3,:],color=:black,markersize=0.01,transparency=true)

Makie.scatter!(s,maxc[1,:],maxc[2,:],maxc[3,:],color=:red,markersize=0.01,transparency=true)

t=mapreduce(i->hcat(minc[:,i],maxc[:,i]),hcat,1:size(minc,2))

linesegments!(s,t[1,:],t[2,:],t[3,:])

s=Makie.scatter(maxc[1,:],maxc[2,:],maxc[3,:],color=[RGB(maxc[:,i]...) for i in size(maxc,2)],markersize=10,transparency=false)

[RGB(minc[:,i]...) for i in 1:size(minc,2)]

[RGB(maxc[:,i]...) for i in 1:size(maxc,2)]
