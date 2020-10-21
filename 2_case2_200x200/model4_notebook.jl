### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 12726064-1195-11eb-380a-0b808c9dc24f
begin
	using ResSimAD
	using DelimitedFiles
	using Statistics
	using Random
	using Plots
	md"""
	###### 1. Import Modules
	"""
end

# ╔═╡ d8c5e688-11a3-11eb-20db-372e8b5f77eb
begin
	using ResSimAD: get_grid_index
	using ResSimAD: isproducer
	## Define function to plot wells
	
	function plot_wells_2d(plt, sim; markersize=2, color=:white)
	    for w in values(sim.facility)
	        i, j, _ = get_grid_index(sim.reservoir.grid, w.ind[1])
	        if isproducer(w)
	            marker = :circle
	        else
	            marker = :dtriangle
	        end
	        scatter!(plt, [j,], [i,], m=(marker, markersize, color), legend=false)
	    end
	end
	nothing
end

# ╔═╡ 17030f0c-119a-11eb-04f8-85058aff795d
using PlutoUI

# ╔═╡ 71340626-11a8-11eb-3910-c3c56ac68c3e
using Markdown

# ╔═╡ 991d1b32-1194-11eb-0828-cd091f95c65a
begin
	# Load logk
	logk = readdlm(joinpath(@__DIR__, "logk_gaussian_200x200.txt"))[:, 1]
	# scale and transform
	logk = logk * 0.3 .+ 6.5
	perm = exp.(logk)
	md"""
	##### 2. Load permeability
	"""
end

# ╔═╡ 10a0270c-1196-11eb-19fc-4fe0c59578e6
begin
	function get_sim(p0, sw, t0)
		options = Dict()

		options["nx"] = 200; options["ny"] = 200; options["nz"] = 1;
		options["dx"] = 200.; options["dy"] = 200.; options["dz"] = 40.;
		options["tops"] = 10000.;
		options["perm"] = perm;
		options["poro"] = 0.2;

		options["fluid"] = "OW"
		options["sw"] = sw
		options["po"] = p0

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



		options["dt0"] = 0.1
		options["dt_max"] = 30.; options["t_end"] = 300.0;
		options["min_err"] = 1.0e-3;

		sim = Sim(options)
		sim.scheduler.t_current = t0
		change_dt(sim, 2.0)
		return sim, options, well_status, prod_locs, inj_locs
	end
	md"""
	##### 3. Define history simluation model
	"""
end

# ╔═╡ 15d1a4ee-1196-11eb-2793-970d5a9c6bf0
begin
	sim_base, options, well_status, prod_locs, inj_locs = get_sim(6000.0, 0.1, 0.0)

	num_prod = length(options["producers"])
	num_inj = length(options["injectors"])
	
	nt = 15
	for it = 1:nt
	    tend = 100.0 * it
	    for ind = 1:num_prod
	        wn = "P$ind"
	        if well_status[wn] == "closed"
	            if rand() < 0.18
	                bhp = rand() * 500. + 5000.
	                change_well_mode(sim_base, wn, "bhp", bhp)
	                well_status[wn] = "open"
	            end
	        elseif well_status[wn] == "open"
	            if rand() < 0.8
	                bhp = rand() * 500. + 5000.
	                change_well_mode(sim_base, wn, "bhp", bhp)
	                well_status[wn] = "open"
	            else
	                change_well_mode(sim_base, wn, "shut", 0.0)
	                well_status[wn] = "shut"
	            end
	        else
	            if rand() < 0.4
	                bhp = rand() * 500. + 5000.
	                change_well_mode(sim_base, wn, "bhp", bhp)
	                well_status[wn] = "open"
	            else
	                change_well_mode(sim_base, wn, "shut", 0.0)
	                well_status[wn] = "shut"
	            end
	        end
	    end
	    for ind = 1:num_inj
	        wn = "I$ind"
	        if well_status[wn] == "closed"
	            if rand() < 0.18
	                bhp = rand() * 500. + 6500.
	                change_well_mode(sim_base, wn, "bhp", bhp)
	                well_status[wn] = "open"
	            end
	        elseif well_status[wn] == "open"
	            if rand() < 0.8
	                bhp = rand() * 500. + 6500.
	                change_well_mode(sim_base, wn, "bhp", bhp)
	                well_status[wn] = "open"
	            else
	                change_well_mode(sim_base, wn, "shut", 0.0)
	                well_status[wn] = "shut"
	            end
	        else
	            if rand() < 0.4
	                bhp = rand() * 500. + 6500.
	                change_well_mode(sim_base, wn, "bhp", bhp)
	                well_status[wn] = "open"
	            else
	                change_well_mode(sim_base, wn, "shut", 0.0)
	                well_status[wn] = "shut"
	            end
	        end
	    end
	    change_dt(sim_base, 2.0)
	    step_to(sim_base, tend)
	end
end

# ╔═╡ 107f13da-1198-11eb-292b-7b2b6dd5170d
begin
	t1 = 1500.
	po = get_state_map(sim_base, "po", t1)
	sw = get_state_map(sim_base, "sw", t1)
	md"""
	##### 4. Get restart po and sw
	"""
end

# ╔═╡ 82be1800-1199-11eb-341b-cdefb292506f
begin
	prod_bhp = rand((num_prod)) * 1000. .+ 4500.
	inj_rate = -(rand((num_inj)) * 20000. .+ 10000.)
	sim, options_new, well_status_new, tmp1, tmp2 = get_sim(po, sw, t1)
	for ind = 1:num_prod
		wn = "P$ind"
		bhp = prod_bhp[ind]
		change_well_mode(sim, wn, "bhp", bhp)
	end
	for ind = 1:num_inj
        wn = "I$ind"
        wrat = inj_rate[ind]
        change_well_mode(sim, wn, "wrat", wrat)
    end
	change_dt(sim, 2.0)
    step_to(sim, 1800.)
	md"""
	##### 5. Run base model with random future controls
	"""
end

# ╔═╡ 7ef75844-11a1-11eb-13e4-73b2d4baba6b
begin
	# Compute additional COPT
	copt = []
	cwpt = []
	tt2 = vcat([1500.,], get_well_rates(sim, "P1", "time"))
	
	for ind=1:num_prod
		wn = "P$ind"
		qoo = get_well_rates(sim, wn, "ORAT")
		qww = get_well_rates(sim, wn, "WRAT")
		push!(copt, sum(diff(tt2) .* qoo))
		push!(cwpt, sum(diff(tt2) .* qww))
	end
	sum(copt)
	md"""
	##### 6. Compute additional COPT from base model
	###### COPT base model: $(round((sum(copt) / 1e6), digits=3)) million bbl.
	###### CWPT base model: $(round((sum(cwpt) / 1e6), digits=3)) million bbl."""
end

# ╔═╡ bd204d62-11a1-11eb-25c5-f710a738608e
md"""
##### 7. Optimize water injection
"""

# ╔═╡ e8342ac6-11a3-11eb-285e-bdfe24cb49a6
md"""
###### 7.1 Find target injectors 
"""

# ╔═╡ 12e65244-119a-11eb-2c46-6fc449069932
begin
	wnames_ = sort(String.(collect(keys(inj_locs))))
	h1 = @bind well_names PlutoUI.MultiSelect(wnames_; default=[wnames_[1]])
	md"""
	**Well names:** $(h1)"""
end

# ╔═╡ ae7acd58-11a3-11eb-3289-71aedcfc8e51
begin
	ttt = get_well_rates(sim, well_names[1], "TIME");
	hslider = @bind tstep Slider(1:length(ttt))
	md"""
	Time step $(hslider)"""
end

# ╔═╡ 4b76be6e-119c-11eb-2fd6-2544fbc1fde2
begin
	cmap = cgrad(:jet)
	# Plot state maps
	sw2 = get_state_map(sim_base, "sw", 1500.);
	p4 = heatmap(reshape(sw2, sim_base.nx, sim_base.ny), color=cmap, title="Sw at day $(round(1500., digits=3))", clim=(0., 1.0));
	plot_wells_2d(p4, sim);
	for wn in collect(keys(inj_locs))
		xxx,yyy, zzz = inj_locs[wn]
		annotate!(p4, [(yyy, xxx, text(wn, 10, :left, :top, :black))])
	end
	for wn in well_names
		if startswith(wn, "P")
			xxx,yyy, zzz = prod_locs[wn]
		else
			xxx,yyy, zzz = inj_locs[wn]
		end
		scatter!(p4, [yyy,], [xxx,], m=(:circle, 4, :red), legend=false)
	end

	plot(p4,size=(500,500))

end

# ╔═╡ b4d6655e-11a3-11eb-14b7-454a4c0a7a1e
begin
	# Plot state maps
	sw33 = get_state_map(sim_base, "sw", 1500.);
	sw_ = get_state_map(sim, "sw", ttt[tstep]);
	p5 = heatmap(reshape(sw33, sim.nx, sim.ny), color=cmap, title="Sw at day 1500.", clim=(0., 1.0));
	plot_wells_2d(p5, sim);
	p6 = heatmap(reshape(sw_, sim.nx, sim.ny), color=cmap, title="Sw at day $(round(ttt[tstep], digits=3))", clim=(0., 1.0));
	plot_wells_2d(p6, sim);
	for wn in well_names
		if startswith(wn, "P")
			xxx,yyy, zzz = prod_locs[wn]
		else
			xxx,yyy, zzz = inj_locs[wn]
		end
		scatter!(p5, [yyy,], [xxx,], m=(:circle, 4, :red), legend=false)
		scatter!(p6, [yyy,], [xxx,], m=(:circle, 4, :red), legend=false)
	end
	plot(p5, p6, layout=(1,2), size=(600,260))
end

# ╔═╡ b4a18654-11a8-11eb-190a-c5b53bb28fce
begin
	h5 = @bind dwrat PlutoUI.NumberField(1000. : 40000.)
	md"""
	##### 8. Select injection rate increasing amount

	Increase injection rate by $(h5) STB/Day for selected injectors
	
	"""
end



# ╔═╡ d17cdbe0-11a3-11eb-23d0-5f250403ded2
begin
	sim_opt, aa, bb, cc, dd = get_sim(po, sw, t1)
	for ind = 1:num_prod
		wn = "P$ind"
		bhp = prod_bhp[ind]
		change_well_mode(sim_opt, wn, "bhp", bhp)
	end
	n1 = length(well_names)
	n2 = num_inj - n1
	tune_down = n1 * dwrat / n2
	println(tune_down)
	for ind = 1:num_inj
        wn = "I$ind"
        wrat = inj_rate[ind] + tune_down
        change_well_mode(sim_opt, wn, "wrat", wrat)
    end
	# Adjust injector rates
	for wn in well_names
		ind = parse(Int64, replace(wn, "I" => ""))
		wrat = inj_rate[ind] - dwrat
		change_well_mode(sim_opt, wn, "wrat", wrat)
	end
	
	change_dt(sim_opt, 2.0)
    step_to(sim_opt, 1800.)
	md"""
	##### 8. Run optimized model 
	"""
end

# ╔═╡ c63a7eba-119f-11eb-3c6f-8d9017d0e9d2
begin
	# Compute additional COPT
	copt_opt = []
	cwpt_opt = []
	tt3 = vcat([1500.,], get_well_rates(sim_opt, "P1", "time"))
	for ind=1:num_prod
		wn = "P$ind"
		qoo = get_well_rates(sim_opt, wn, "ORAT")
		qww = get_well_rates(sim_opt, wn, "WRAT")
		push!(copt_opt, sum(diff(tt3) .* qoo))
		push!(cwpt_opt, sum(diff(tt3) .* qww))
	end
	copt_sum = round((sum(copt) / 1e6), digits=3)
	cwpt_sum = round((sum(cwpt) / 1e6), digits=3)
	copt_opt_sum = round((sum(copt_opt) / 1e6), digits=3)
	cwpt_opt_sum = round((sum(cwpt_opt) / 1e6), digits=3)
	Markdown.parse(
	"""
	
	##### 9. Improvements from optimized model
	
	  Data   | Base Model | Optimized Model | Improvement
	:---:   | :---:| :---:| :---:
	FOPT (MBBL) | $copt_sum |$copt_opt_sum | $(round(copt_opt_sum - copt_sum, digits=3))|
	FWPT (MBBL) | $cwpt_sum |$cwpt_opt_sum | $(round(cwpt_sum - cwpt_opt_sum, digits=3))|
	
	"""
	)
end

# ╔═╡ 9fa2a0c0-11a4-11eb-3173-bfc775ecf370
begin
	# Plot state maps
	sw3 = get_state_map(sim_opt, "sw", 1800.);
	sw4 = get_state_map(sim, "sw", 1800.);
	p8 = heatmap(reshape(sw4, sim.nx, sim.ny), color=cmap, title="Sw base model", clim=(0., 1.0));
	plot_wells_2d(p8, sim_opt);
	p7 = heatmap(reshape(sw3, sim_opt.nx, sim_opt.ny), color=cmap, title="Sw opt model", clim=(0., 1.0));
	plot_wells_2d(p7, sim_opt);
	for wn in well_names
		if startswith(wn, "P")
			xxx,yyy, zzz = prod_locs[wn]
		else
			xxx,yyy, zzz = inj_locs[wn]
		end
		scatter!(p7, [yyy,], [xxx,], m=(:circle, 4, :red), legend=false)
		scatter!(p8, [yyy,], [xxx,], m=(:circle, 4, :red), legend=false)
	end
	plot(p8, p7, size=(680,280))
end

# ╔═╡ Cell order:
# ╟─12726064-1195-11eb-380a-0b808c9dc24f
# ╟─d8c5e688-11a3-11eb-20db-372e8b5f77eb
# ╟─991d1b32-1194-11eb-0828-cd091f95c65a
# ╟─10a0270c-1196-11eb-19fc-4fe0c59578e6
# ╟─15d1a4ee-1196-11eb-2793-970d5a9c6bf0
# ╟─107f13da-1198-11eb-292b-7b2b6dd5170d
# ╟─82be1800-1199-11eb-341b-cdefb292506f
# ╟─7ef75844-11a1-11eb-13e4-73b2d4baba6b
# ╟─ae7acd58-11a3-11eb-3289-71aedcfc8e51
# ╟─b4d6655e-11a3-11eb-14b7-454a4c0a7a1e
# ╟─bd204d62-11a1-11eb-25c5-f710a738608e
# ╟─17030f0c-119a-11eb-04f8-85058aff795d
# ╟─e8342ac6-11a3-11eb-285e-bdfe24cb49a6
# ╟─12e65244-119a-11eb-2c46-6fc449069932
# ╟─4b76be6e-119c-11eb-2fd6-2544fbc1fde2
# ╟─b4a18654-11a8-11eb-190a-c5b53bb28fce
# ╟─d17cdbe0-11a3-11eb-23d0-5f250403ded2
# ╟─71340626-11a8-11eb-3910-c3c56ac68c3e
# ╟─c63a7eba-119f-11eb-3c6f-8d9017d0e9d2
# ╟─9fa2a0c0-11a4-11eb-3173-bfc775ecf370
