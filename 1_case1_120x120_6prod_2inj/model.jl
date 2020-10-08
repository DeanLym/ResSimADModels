using Revise
using ResSimAD
using HDF5
using Plots
using Statistics
using ResSimAD.LinearSolver: BICGSTAB_ILU_DUNE_ISTL_Solver, setup_lsolver

## Load permeability
permx = h5read(joinpath(@__DIR__, "perm.h5"), "data")
## rescale permeability
logk = log.(permx)
logk = (logk .- minimum(logk)) / (maximum(logk) - minimum(logk))
logk = logk * 5.5 .+ 1.5
permx = exp.(logk)

## Specify options
options = Dict()

options["nx"] = 120; options["ny"] = 120; options["nz"] = 1;
options["dx"] = 80.; options["dy"] = 80.; options["dz"] = 20.;
options["tops"] = 10000.;
options["perm"] = permx;
options["poro"] = 0.2;

options["fluid"] = "OW"
options["sw"] = 0.1
options["po"] = 6000.

options["PVDO"] = get_example_data("PVDO.DAT")

options["PVTW"] = get_example_data("PVTW.DAT")

options["SWOF"] = get_example_data("SWOF.DAT")

# options["PVCDO"] = Dict([
#     ("pref", 14.7),
#     ("bref", 1.03),
#     ("c", 0.0),
#     ("μref", 3.2),
#     ("cμ", 0.0),
# ])

# options["PVTW"] = Dict([
#     ("pref", 14.7),
#     ("bref", 1.001),
#     ("c", 1.0e-6),
#     ("μref", 0.3),
#     ("cμ", 0.0),
# ])

# options["SWOF"] = joinpath(@__DIR__, "SWOF.DAT")

options["producers"] = [];
num_prod = 6
well_locs = [(20, 20, 1), (60, 20, 1), (100, 20, 1), (20, 100, 1), (60, 100, 1), (100, 100, 1)]
for ind = 1:num_prod
    well = Dict();
    well["name"] = "P$ind";
    well["perforation"] = [well_locs[ind]];
    well["radius"] = 0.5;
    well["mode"] = "orat";
    well["target"] = 800.;
    well["limits"] = [("min_bhp", 5000.0)]
    push!(options["producers"], well);
end

options["injectors"] = []
num_inj = 2
well_locs = [(44, 60, 1), (76, 60, 1)]
for ind = 1:num_inj
    well = Dict();
    well["name"] = "I$ind";
    well["perforation"] = [well_locs[ind]];
    well["radius"] = 0.5;
    well["mode"] = "wrat";
    well["target"] = -2000.;
    push!(options["injectors"], well);
end

options["dt0"] = 0.1
options["dt_max"] = 30.; options["t_end"] = 3000.0;
options["min_err"] = 1.0e-3;

options["linear_solver"] = "GMRES_ILU"


## Run simluation
sim = Sim(options)

runsim(sim)

sim.lsolver = BICGSTAB_ILU_DUNE_ISTL_Solver()

setup_lsolver(sim.lsolver, sim.reservoir.grid)

runsim(sim)

## Plot results
using Plots
plts = []
for k = 1:6
    wn = "P$k"
    t = get_well_rates(sim, wn, "TIME")
    qo = get_well_rates(sim, wn, "ORAT")
    qw = get_well_rates(sim, wn, "WRAT")
    p1 = plot(t, qo, label="ORAT")
    plot!(p1, t, qw, label="WRAT")
    xlabel!(p1, "Day")
    ylabel!(p1, "Rate (STB/Day)")
    title!(p1, wn)
    push!(plts, p1)
end
plot(plts..., layout=(2,3), legend=false, size=(700, 300))

plts = []
for k = 1:6
    wn = "P$k"
    t = get_well_rates(sim, wn, "TIME")
    bhp = get_well_rates(sim, wn, "WBHP")
    p1 = plot(t, bhp, label="WBHP")
    xlabel!(p1, "Day")
    ylabel!(p1, "BHP (psi)")
    title!(p1, wn)
    push!(plts, p1)
end
plot(plts..., layout=(2,3), legend=false, size=(700, 300))

plts = []
for k = 1:2
    wn = "I$k"
    t = get_well_rates(sim, wn, "TIME")
    qw = get_well_rates(sim, wn, "WRAT")
    bhp = get_well_rates(sim, wn, "WBHP")
    p1 = plot(t, -qw, label="WRAT")
    xlabel!(p1, "Day")
    ylabel!(p1, "Inj. rate (STB/Day)")
    title!(p1, wn)
    push!(plts, p1)
    p1 = plot(t, bhp, label="BHP")
    xlabel!(p1, "Day")
    ylabel!(p1, "BHP (psi)")
    title!(p1, wn)
    push!(plts, p1)
end

plot(plts..., layout=(2,2))
