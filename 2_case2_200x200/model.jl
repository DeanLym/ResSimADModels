using ResSimAD
using HDF5
using CSV
using Plots
using Statistics
using DelimitedFiles
using JSON


# Load logk
logk = readdlm(joinpath(@__DIR__, "logk_gaussian_200x200.txt"))[:, 1]
# scale and transform
logk = logk * 0.3 .+ 6.5
perm = exp.(logk)
maximum(perm)
minimum(perm)
## Specify options
options = Dict()

options["nx"] = 200; options["ny"] = 200; options["nz"] = 1;
options["dx"] = 200.; options["dy"] = 200.; options["dz"] = 40.;
options["tops"] = 10000.;
options["perm"] = perm;
options["poro"] = 0.2;

options["fluid"] = "OW"
options["sw"] = 0.1
options["po"] = 6000.

options["PVDO"] = get_example_data("PVDO.DAT")

options["PVTW"] = get_example_data("PVTW.DAT")

options["SWOF"] = get_example_data("SWOF.DAT")

options["producers"] = []

prod_locs = Dict()

for i = 1:10
    for j = 1:10
        ind = 10*(j-1) + i
        well_loc = (10 + 20*(i-1), 10+20*(j-1), 1)
        wn = "P$ind"
        well = Dict();
        well["name"] = wn;
        well["perforation"] = [well_loc];
        well["radius"] = 0.5;
        well["mode"] = "bhp";
        well["target"] = 5500.;
        push!(options["producers"], well);
        prod_locs[wn] = well_loc
    end
end

inj_locs = Dict()
options["injectors"] = []
for i = 1:9
    for j = 1:9
        ind = 9*(j-1) + i
        well_loc = (20 + 20*(i-1), 20+20*(j-1), 1)
        well = Dict();
        wn = "I$ind"
        well["name"] = wn;
        well["perforation"] = [well_loc];
        well["radius"] = 0.5;
        well["mode"] = "bhp";
        well["target"] = 6500.;
        push!(options["injectors"], well);
        inj_locs[wn] = well_loc
    end
end

num_inj = length(options["injectors"])

options["dt0"] = 0.1
options["dt_max"] = 30.; options["t_end"] = 300.0;
options["min_err"] = 1.0e-3;

options["linear_solver"] = "BICGSTAB_ILU_DUNE_ISTL"


## Run simluation
sim = Sim(options)

# Run sim with random controls
nt = 30
for it = 1:nt
    tend = 100.0 * it
    for ind = 1:num_inj
        bhp = rand() * 1000. + 6100.
        change_well_target(sim, "I$ind", bhp)
    end
    step_to(sim, tend)
end

## Plot results
plts = []
for k = 1:6
    wn = "P$(2*k+30)"
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
    qo = get_well_rates(sim, wn, "ORAT")
    qw = get_well_rates(sim, wn, "WRAT")
    p1 = plot(t, qw ./ (qo .+ qw), label="WWCT")
    xlabel!(p1, "Day")
    ylabel!(p1, "Water cut")
    title!(p1, wn)
    push!(plts, p1)
end
plot(plts..., layout=(2,3), legend=false, size=(700, 300))

p1 = plot(legend=false)
for k = 1:100
    wn = "P$k"
    t = get_well_rates(sim, wn, "TIME")
    qo = get_well_rates(sim, wn, "ORAT")
    qw = get_well_rates(sim, wn, "WRAT")
    plot!(p1, t, qw ./ (qo .+ qw))
    xlabel!(p1, "Day")
    ylabel!(p1, "Water cut")
    title!(p1, wn)
end
plot(p1)
plot(plts..., layout=(2,3), legend=false, size=(700, 300))



plts = []
for k = 1:6
    wn = "I$k"
    t = get_well_rates(sim, wn, "TIME")
    qw = get_well_rates(sim, wn, "WRAT")
    p1 = plot(t, -qw, label="WRAT")
    xlabel!(p1, "Day")
    ylabel!(p1, "Rate (STB/Day)")
    title!(p1, wn)
    push!(plts, p1)
end
plot(plts..., layout=(2,3), legend=false, size=(700, 300))

for it = 1:5
    t = it*100.
    qo = get_state_map(sim, "sw", t)

    p1 = heatmap(reshape(qo, 200, 200), color=cgrad(:jet))
    display(p1)
end

for it = 1:5
    t = it*100.
    qo = get_state_map(sim, "po", t)

    p1 = heatmap(reshape(qo, 200, 200), color=cgrad(:jet))
    display(p1)
end

# Save data
for well in values(sim.facility)
    wn = well.name
    CSV.write(joinpath(@__DIR__, "well_rates", "$wn.csv"), well.results)
end

# Save well_locations
open(joinpath(@__DIR__, "prod_locs.json"), "w") do f
    JSON.print(f, prod_locs, 4)
end

open(joinpath(@__DIR__, "inj_locs.json"), "w") do f
    JSON.print(f, inj_locs, 4)
end

# Save well index
well_inds = Dict()
for well in values(sim.facility)
    wn = well.name
    well_inds[wn] = well.ind[1]
end

open(joinpath(@__DIR__, "well_inds.json"), "w") do f
    JSON.print(f, well_inds, 4)
end