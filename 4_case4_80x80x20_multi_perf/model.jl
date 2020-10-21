using ResSimAD
using DelimitedFiles

## Specify input
# Grid and Rock
options = Dict()
options["nx"] = 80; options["ny"] = 80; options["nz"] = 20;
options["dx"] = 65.6168; options["dy"] = 65.6168; options["dz"] = 6.56168;
options["tops"] = 26184.38;
options["poro"] = joinpath(@__DIR__, "poro.dat");
options["perm"] = joinpath(@__DIR__, "PERMX.dat");
# Fluid
options["fluid"] = "OW"
options["equil"] = (26184.38, 4713.73);
options["sw"] = 0.1;
options["PVDO"] = joinpath(@__DIR__, "PVDO.DAT")
options["PVTW"] = get_example_data("PVTW.DAT")
options["SWOF"] = joinpath(@__DIR__, "SWOF.DAT")
# Wells
prods = Dict([
    ("P1", (9, 65, 1, 5)),
    ("P2", (23, 40, 14, 18)),
    ("P3", (44, 25, 6, 10)),
    ("P4", (55, 44, 9, 13)),
    ("P5", (70, 62, 1, 5)),
])
num_prod = length(prods)
options["producers"] = [];
for ind = 1:num_prod
    well = Dict()
    wn = "P$ind"
    x, y, z1, z2 = prods[wn]
    well["name"] = wn;
    well["perforation"] = [(x, y, z) for z = z1:z2]
    well["radius"] = 0.328084;
    well["mode"] = "bhp";
    well["target"] = 4496.17;
    push!(options["producers"], well);
end


injs = Dict([
    ("I1", (16, 15, 1, 6)),
    ("I2", (36, 58, 11, 16)),
    ("I3", (62, 16, 7, 11)),
])

num_inj = length(injs)
options["injectors"] = []
for ind = 1:num_inj
    well = Dict()
    wn = "I$ind"
    x, y, z1, z2 = injs[wn]
    well["name"] = wn;
    well["perforation"] = [(x, y, z) for z = z1:z2]
    well["radius"] = 0.328084;
    well["mode"] = "bhp";
    well["target"] = 4786.25;
    push!(options["injectors"], well);
end

# Nonlinear solver options
options["max_newton_iter"] = 50
options["min_err"] = 5.0e-3

# Schedule
options["dt0"] = 0.01
options["dt_max"] = 50.; options["t_end"] = 1000.0;

##
 
sim = Sim(options)

runsim(sim)

step_to(sim, 1000.0)

## Load ADGPRS results
using DelimitedFiles
using DataFrames
results = Dict()
results["ADGPRS"] = Dict()

file = joinpath(@__DIR__, "ADGPRS", "STATES_VARS.rates.txt")
header = String.(readdlm(file)[1,:])
data = readdlm(file,skipstart=1)
df = DataFrame(data)
rename!(df, header)

data_types = Dict([
    ("P1", [("oil rate", "OPR"), ("water rate", "WPR",)]),
    ("P2", [("oil rate", "OPR"), ("water rate", "WPR",)]),
    ("P3", [("oil rate", "OPR"), ("water rate", "WPR",)]),
    ("P4", [("oil rate", "OPR"), ("water rate", "WPR",)]),
    ("P5", [("oil rate", "OPR"), ("water rate", "WPR",)]),
    ("I1", [("water inj. rate", "WIR")]),
    ("I2", [("water inj. rate", "WIR")]),
    ("I3", [("water inj. rate", "WIR")]),
])

results["ADGPRS"]["Day"] = df[!, "Day"]

wnames = collect(keys(data_types))

for (iw, wn) in enumerate(wnames)
    for (key, col) in data_types[wn]
        results["ADGPRS"][wn * " " * key] = df[!, string(wn, ":", col)] *6.28981
    end
end


res = results["ADGPRS"]

using Plots
using Plots.PlotMeasures
# production curves
for ind = 1:num_prod
    ofs = 5
    wn = "P$ind"
    t = get_well_rates(sim, wn, "TIME")
    qo = get_well_rates(sim, wn, "ORAT");
    qw = get_well_rates(sim, wn, "WRAT")
    p1 = plot(t[ofs:end], qo[ofs:end], color=:black, marker=true, label="ResSimAD",
        xlabel="Day", ylabel="Oil Rate (STB/Day)",
        title="$wn Oil Rate")
    plot!(p1, res["Day"][ofs:end], res["$wn oil rate"][ofs:end], color=:red, label="ADGPRS")

    p2 = plot(t[ofs:end], qw[ofs:end], color=:black, marker=true,label="ResSimAD",
     xlabel="Day", ylabel="Water Rate (STB/Day)",
     title="$wn Water Rate")
    plot!(p2, res["Day"][ofs:end], res["$wn water rate"][ofs:end], color=:red, label="ADGPRS")

    display(plot(p1, p2, layout=(1,2), legend=true, size=(720,220), bottom_margin = 10px))
end

