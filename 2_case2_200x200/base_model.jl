using ResSimAD
using DelimitedFiles
using Statistics
using Random

# Load logk
logk = readdlm(joinpath(@__DIR__, "logk_gaussian_200x200.txt"))[:, 1]
# scale and transform
logk = logk * 0.3 .+ 6.5
perm = exp.(logk)

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
well_status = Dict()
for i = 1:10
    for j = 1:10
        ind = 10*(j-1) + i
        well_loc = (10 + 20*(i-1), 10+20*(j-1), 1)
        wn = "P$ind"
        well = Dict();
        well["name"] = wn;
        well["perforation"] = [well_loc];
        well["radius"] = 0.5;
        well["mode"] = "shut";
        well["target"] = 5500.;
        push!(options["producers"], well);
        prod_locs[wn] = well_loc
        well_status[wn] = "closed"
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
        well["mode"] = "shut";
        well["target"] = 6500.;
        push!(options["injectors"], well);
        inj_locs[wn] = well_loc
        well_status[wn] = "closed"
    end
end

num_prod = length(options["producers"])
num_inj = length(options["injectors"])

options["dt0"] = 0.1
options["dt_max"] = 30.; options["t_end"] = 300.0;
options["min_err"] = 1.0e-3;
# options["linear_solver_backend"] = "Julia"
## Run simulation for 1500 days
sim = Sim(options)

rng = MersenneTwister(12345)

nt = 15
for it = 1:nt
    tend = 100.0 * it
    for ind = 1:num_prod
        wn = "P$ind"
        if well_status[wn] == "closed"
            if rand(rng) < 0.18
                bhp = rand(rng) * 500. + 5000.
                change_well_mode(sim, wn, "bhp", bhp)
                well_status[wn] = "open"
            end
        elseif well_status[wn] == "open"
            if rand(rng) < 0.8
                bhp = rand(rng) * 500. + 5000.
                change_well_mode(sim, wn, "bhp", bhp)
                well_status[wn] = "open"
            else
                change_well_mode(sim, wn, "shut", 0.0)
                well_status[wn] = "shut"
            end
        else
            if rand(rng) < 0.4
                bhp = rand(rng) * 500. + 5000.
                change_well_mode(sim, wn, "bhp", bhp)
                well_status[wn] = "open"
            else
                change_well_mode(sim, wn, "shut", 0.0)
                well_status[wn] = "shut"
            end
        end
    end
    for ind = 1:num_inj
        wn = "I$ind"
        if well_status[wn] == "closed"
            if rand(rng) < 0.18
                bhp = rand(rng) * 500. + 6500.
                change_well_mode(sim, wn, "bhp", bhp)
                well_status[wn] = "open"
            end
        elseif well_status[wn] == "open"
            if rand(rng) < 0.8
                bhp = rand(rng) * 500. + 6500.
                change_well_mode(sim, wn, "bhp", bhp)
                well_status[wn] = "open"
            else
                change_well_mode(sim, wn, "shut", 0.0)
                well_status[wn] = "shut"
            end
        else
            if rand(rng) < 0.4
                bhp = rand(rng) * 500. + 6500.
                change_well_mode(sim, wn, "bhp", bhp)
                well_status[wn] = "open"
            else
                change_well_mode(sim, wn, "shut", 0.0)
                well_status[wn] = "shut"
            end
        end
    end
    change_dt(sim, 2.0)
    step_to(sim, tend)
end