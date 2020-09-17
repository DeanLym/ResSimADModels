## Load perm
perm = Vector{Float64}()

fid = open(joinpath(@__DIR__, "Perm_1.GRDECL"), "r")
lines = readlines(fid)
close(fid)

for line in lines
    append!(perm, parse.(Float64, split(line)))
end

## Save perm to h5 file
using HDF5
h5write(joinpath(@__DIR__, "perm.h5"), "data", perm)

using Plots
cmap = cgrad(:jet)

heatmap(reshape(log.(perm), 120, 120), color=cmap)