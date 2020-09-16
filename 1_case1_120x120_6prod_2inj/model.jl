using Revise
using ResSimAD

options = Dict()

options["nx"] = 120; options["ny"] = 120; options["nz"] = 1;
options["dx"] = 20.; options["dy"] = 20.; options["dz"] = 10.;
options["tops"] = 10000.;
options["perm"] = 200.;
options["poro"] = 0.2;

options["fluid"] = "OW";
options["sw"] = 0.1;
options["po"] = 6000.;

options["PVCDO"] = Dict([
    ("pref", 14.7),
    ("bref", 1.03),
    ("c", 0.0),
    ("μref", 3.2),
    ("cμ", 0.0),
])

options["PVTW"] = Dict([
    ("pref", 14.7),
    ("bref", 1.001),
    ("c", 1.0e-6),
    ("μref", 0.8),
    ("cμ", 0.0),
])

options["SWOF"] = joinpath(@__DIR__, "SWOF.DAT");

options["producers"] = [];
num_prod = 6
well_locs = [(20, 20, 1), (60, 20, 1), (100, 20, 1), (20, 100, 1), (60, 100, 1), (100, 100, 1)]
for ind = 1:num_prod
    well = Dict();
    well["name"] = "P$ind";
    well["perforation"] = [well_locs[ind]];
    well["radius"] = 0.5;
    well["mode"] = "bhp";
    well["target"] = 5900.;
    push!(options["producers"], well);
end

options["injectors"] = [];
num_inj = 2
well_locs = [(44, 60, 1), (76, 60, 1)]
for ind = 1:num_inj
    well = Dict();
    well["name"] = "I$ind";
    well["perforation"] = [well_locs[ind]];
    well["radius"] = 0.5;
    well["mode"] = "wrat";
    well["target"] = -100.;
    push!(options["injectors"], well);
end

options["dt0"] = 0.1
options["dt_max"] = 30.; options["t_end"] = 3000.0;
options["min_err"] = 1.0e-3;

options["linear_solver"] = "GMRES_CPR";

sim = Sim(options)

runsim(sim)