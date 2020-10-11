using ResSimAD
using DelimitedFiles

bar2psi = 14.50377
cm2bbl = 6.289811
cm2cft = 35.31466
β = cm2bbl / bar2psi
tranx = β * readdlm(joinpath(@__DIR__, "TransX.txt"), skipstart=1)[:, 1]
trany = β * readdlm(joinpath(@__DIR__, "TransY.txt"), skipstart=1)[:, 1]
##
options = Dict()

# Grid and Rock
options["nx"] = 10; options["ny"] = 10; options["nz"] = 1;
options["v"] = 2858291.712 * cm2cft
options["d"] = 3280.
options["tranx"] = tranx
options["trany"] = trany
options["tranz"] = 1.0
options["poro"] = 0.2;

# Fluid
options["fluid"] = "OW"
options["po"] = 200*bar2psi; options["sw"] = 0.1;

options["PVDO"] = joinpath(@__DIR__, "PVDO.DAT")
options["PVTW"] = joinpath(@__DIR__, "PVTW.DAT")
options["SWOF"] = joinpath(@__DIR__, "SWOF.DAT")

# Preparation of random well locations
wellLoc = [[[4, 6], [6, 5]], [[3, 1], [9, 2], [1, 9], [8, 8]]]

wellIndex = [[70.106492, 24.897387], [16.752537, 152.879541, 94.949679, 33.423720]]


# Wells

# injectors
options["injectors"] = []

for ind = 1:2
    w = Dict();
    w["name"] = "I$ind"; w["perforation"] = [(wellLoc[1][ind][1],wellLoc[1][ind][2],1)]; 
    w["mode"] = "bhp"; w["target"] = 250*bar2psi; w["wi"] = [wellIndex[1][ind] * β];
    push!(options["injectors"], w)
end

# producers
options["producers"] = []
for ind = 1:4
    w = Dict();
    w["name"] = "P$ind"; w["perforation"] = [(wellLoc[2][ind][1],wellLoc[2][ind][2],1)];
    w["mode"] = "bhp"; w["target"] = 150*bar2psi; w["wi"] = [wellIndex[2][ind] * β];
    push!(options["producers"], w)
end


# Schedule
options["dt0"] = 0.1;
options["dt_max"] = 30.;
options["t_end"] = 10000.0;
options["min_err"] = 1.0e-6


## 
sim = Sim(options)

runsim(sim)

## Compare results with ADGPRS

using Plots
using CSV
using DelimitedFiles
using DataFrames
# Load ADGPRS results
file = joinpath(@__DIR__, "ADGPRS","OUTPUT.rates.txt")

header = String.(readdlm(file)[1,:])
data = readdlm(file,skipstart=1)
df = DataFrame(data)
rename!(df, header)

plts = []

for ind = 1:4
    # plot ressimad results
    t = get_well_rates(sim, "P$ind", "time")
    qo = get_well_rates(sim, "P$ind", "orat")
    ppp = plot(t, qo,label="ResSimAD")
    # plot ADGPRS results
    t = Float64.(df[:Day][1:end-1])
    qo =  cm2bbl *Float64.(df["PRD$ind:OPR"][1:end-1])
    plot!(ppp, t, qo, label="ADGPRS")
    push!(plts, ppp)
end
plot(plts..., layout=(2,2))


plts = []
for ind = 1:4
    # plot ressimad results
    t = get_well_rates(sim, "P$ind", "time")
    qo = get_well_rates(sim, "P$ind", "wrat")
    ppp = plot(t, qo,label="ResSimAD")
    # plot ADGPRS results
    t = Float64.(df[:Day][1:end-1])
    qo =  cm2bbl *Float64.(df["PRD$ind:WPR"][1:end-1])
    plot!(ppp, t, qo, label="ADGPRS")
    push!(plts, ppp)
end
plot(plts..., layout=(2,2), legend=:left)

plts = []
for ind = 1:2
    # plot ressimad results
    t = get_well_rates(sim, "I$ind", "time")
    q = -get_well_rates(sim, "I$ind", "wrat")
    ppp = plot(t, q,label="ResSimAD")
    # plot ADGPRS results
    t = Float64.(df[:Day][1:end-1])
    qo =  -cm2bbl *Float64.(df["INJ$ind:WIR"][1:end-1])
    plot!(ppp, t, qo, label="ADGPRS")
    push!(plts, ppp)
end
plot(plts..., layout=(1,2), legend=:left)
