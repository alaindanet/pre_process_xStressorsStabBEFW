using EcologicalNetworksDynamics
using SparseArrays
using DataFrames
using Arrow

dir = ARGS[1]
println("Directory: $(dir)")
#dir = /mnt/parscratch/users/bi1ahd/sim/simCSh

raw_files = readdir(dir)
#all_files = raw_files[.!occursin.(raw_files)]
all_files = raw_files
mask_ts_files = occursin.("ts", all_files)
files = all_files[.!mask_ts_files]

println("files: $(files)")

new_dir = joinpath(joinpath(dir, "without_cvsp"))

if !isdir(new_dir)
    mkdir(new_dir)
end

for i in files[isfile.([joinpath([dir, i]) for i in files])]
    println("$i")
    df = DataFrame(Arrow.Table(joinpath([dir, i])))
    select!(df, Not([:cv_sp]))
    println("$(names(df))")
    #Arrow.write(joinpath([new_dir, i]), df)
    println("$(joinpath([new_dir, i]))")

end
