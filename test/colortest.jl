"Get RGB color spectra measured from a specific display"
function RGBSpectra(measurement)
    C = map(i->parse.(Float64,split(i)),measurement["Color"])
    λ = measurement["WL"]
    I = measurement["Spectral"]
    return C,λ,I
end


displayname = "ROGPG279Q" # ASUS ROG Swift PG279Q IPS LCD
displaymeasure = YAML.load_file(joinpath(@__DIR__,"DisplayMeasurement.yaml"))
Cs = RGBSpectra(displaymeasure[displayname]["SpectralMeasurement"])
RGBToLMS,LMSToRGB = RGBLMSMatrix(Cs...;observer=10)
RGBToXYZ,XYZToRGB = RGBXYZMatrix(Cs...;observer=10)

RGBs = homovector(rand(3,10))
XYZs = RGBToXYZ*RGBs
xyYs = XYZ2xyY(XYZs)
@test xyY2XYZ(xyYs) ≈ XYZs

BGrgb = [0.5,0.5,0.5,1]
BGlms = RGBToLMS*BGrgb
LMSToContrast,ContrastToLMS = LMSContrastMatrix(BGlms)

CCs = LMSToContrast*RGBToLMS*RGBs
@test LMSToRGB*ContrastToLMS*CCs ≈ RGBs 

LMSToDKL,DKLToLMS = LMSDKLMatrix(BGlms;observer=10)
DKLs = LMSToDKL*RGBToLMS*RGBs
@test LMSToRGB*DKLToLMS*DKLs ≈ RGBs

x = rand(-3:0.1:3,3,6)
desaturate2gamut!(x)
