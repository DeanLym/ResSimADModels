### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ f4c5cc0a-0513-11eb-1630-2d527cb3810e
begin
	using ResSimAD
	using HDF5
	using CSV
	using Plots
	using Statistics
	using DelimitedFiles
	using JSON
end

# ╔═╡ 408067f2-0514-11eb-1ae4-2b7b8f026784
begin
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
end

# ╔═╡ 51c50a24-0514-11eb-3b34-a79bdc87d622
sim = Sim(options)

# ╔═╡ 66d7ceb0-0514-11eb-1805-11b892f29047
runsim(sim)

# ╔═╡ 7beca9c4-0514-11eb-2db6-cb30568bf529


# ╔═╡ Cell order:
# ╠═f4c5cc0a-0513-11eb-1630-2d527cb3810e
# ╠═408067f2-0514-11eb-1ae4-2b7b8f026784
# ╠═51c50a24-0514-11eb-3b34-a79bdc87d622
# ╠═66d7ceb0-0514-11eb-1805-11b892f29047
# ╠═7beca9c4-0514-11eb-2db6-cb30568bf529
