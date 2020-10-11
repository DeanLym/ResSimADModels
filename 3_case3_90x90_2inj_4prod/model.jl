using ResSimAD

##
options = Dict()

# Grid and Rock

options["nx"] = 90; options["ny"] = 90; options["nz"] = 1;
options["dx"] = 107.6115; options["dy"] = 107.6115; options["dz"] = 107.6115;
options["d"] = 3280.
options["perm"] = joinpath(@__DIR__, "PERMMULTIPLIER.dat")

options["poro"] = 0.2;


# Fluid
options["fluid"] = "OW"
options["po"] = 2900.75; options["sw"] = 0.1;

options["PVDO"] = joinpath(@__DIR__, "PVDO.DAT")
options["PVTW"] = joinpath(@__DIR__, "PVTW.DAT")
options["SWOF"] = joinpath(@__DIR__, "SWOF.DAT")


# Preparation of random well locations
wellLoc = [[[38, 50], [73, 32]], [[21, 23], [77, 10], [8, 71], [64, 84]]]

# Wells

# injectors
options["injectors"] = []

for ind = 1:2
    w = Dict();
    w["name"] = "I$ind"; w["perforation"] = [(wellLoc[1][ind][1],wellLoc[1][ind][2],1)]; w["radius"] = 0.5;
    w["mode"] = "bhp"; w["target"] = 250*14.5038;
    push!(options["injectors"], w)
end

# producers
options["producers"] = []
for ind = 1:4
    w = Dict();
    w["name"] = "P$ind"; w["perforation"] = [(wellLoc[2][ind][1],wellLoc[2][ind][2],1)]; w["radius"] = 0.5;
    w["mode"] = "bhp"; w["target"] = 150*14.5038;
    push!(options["producers"], w)
end


# Schedule

options["dt0"] = 1.;
options["dt_max"] = 100.;
options["t_end"] = 10000.0;
options["min_err"] = 1.0e-3


## 
sim = Sim(options)


runsim(sim)

## Compare results with ADGPRS

using Plots
using CSV
using DelimitedFiles

file = joinpath(@__DIR__, "ADGPRS","OUTPUT.rates.txt")

header = String.(readdlm(file)[1,:])
data = readdlm(file,skipstart=1)
df = DataFrame(data)
rename!(df, header)

α = 6.28981
plts = []
for ind = 1:4
    # plot ressimad results
    t = get_well_rates(sim, "P$ind", "time")
    qo = get_well_rates(sim, "P$ind", "orat")
    ppp = plot(t, qo,label="ResSimAD")
    # plot ADGPRS results
    t = Float64.(df[:Day][1:end-1])
    qo =  α *Float64.(df["PRD$ind:OPR"][1:end-1])
    plot!(ppp, t, qo, label="ADGPRS")
    push!(plts, ppp)
end
plot(plts..., layout=(2,2))

## Run step by step and retrieve information
sim = Sim(options)

λo = []
for t in [10.0, 30.0, 50.0]
    step_to(sim, t)
    push!(λo, get_data(sim, "λo"))
end


f = value(sim.λo) ./ (value(sim.λo)+ value(sim.λw))

