using ResSimAD
using HDF5
using CSV
using Plots
using Statistics
using DelimitedFiles
using JSON
using Distributed
using Random


## Launch subprocesses
nrun = 10
addprocs(nrun)
workers()
##
@everywhere using ResSimAD
@everywhere using CSV
@everywhere using DelimitedFiles
@everywhere using Random
@everywhere function run(irun)
    include(joinpath(@__DIR__, "base_model.jl"))
    rng = MersenneTwister(12345*irun)
    nt2 = 20
    for it = nt+1:nt2
        tend = 100.0*it
        for ind = 1:num_prod
            wn = "P$ind"
            bhp = rand(rng) * 5000. + 900.
            change_well_mode(sim, wn, "bhp", bhp)
        end

        for ind = 1:num_inj
            wn = "I$ind"
            wrat = -(rand(rng) * 20000. + 10000.)
            change_well_mode(sim, wn, "wrat", wrat)
        end
        change_dt(sim, 2.0)
        step_to(sim, tend)
    end

    save_dir = joinpath(@__DIR__, "well_rates_model4")
    subdir = joinpath(save_dir, "run$irun")
    mkdir(subdir)
    # Save data
    for well in values(sim2.facility)
        wn = well.name
        CSV.write(joinpath(subdir, "$wn.csv"), well.results)
    end
end

data_ref = []
for (irun, worker) in enumerate(workers())
    push!(data_ref, @spawnat worker run(irun))
end

results = []
for (irun, worker) in enumerate(workers())
    push!(results, fetch(data_ref[irun]))
end


a = 1

for worker in workers()
    rmprocs(worker)
end

using ResSimAD
@doc change_well_mode

workers()

