

println(@__DIR__)
include("qising.jl")

using JLD2
using LaTeXStrings
using Plots
using ProgressMeter

_DATA_PATH0 = "run/data"
datapath = _DATA_PATH0
if !isdir(datapath)
    mkpath(datapath)
end
figpath = _DATA_PATH0

N = 200

dgx = pi / 400
gxs = pi/4:dgx:(pi/2-dgx)
gs = tan.(gxs)
# gs = [4, 8, 16, 32, 64]
len_gs = length(gs)

sumzizjs = Array{Float64}(undef, len_gs)
sumzizjs2 = Array{Float64}(undef, len_gs)
slopes = Array{Float64}(undef, len_gs)

pbar = Progress(len_gs; desc="run len_gs avgzizjs data..")
print(pbar)
for (j, g) in enumerate(gs)
    sites, psi0, energy, avgz, avgzizjs, ee = gen_tising(N, g; bc="obc")

    sm = sum(avgzizjs[j1, j2] for j1 in 1:N for j2 in j1+1:N)
    sm2 = sum(avgzizjs[j1, j2]^2 for j1 in 1:N for j2 in j1+1:N)
    sumzizjs[j] = sm
    sumzizjs2[j] = sm2

    s = Int(floor(N / 4))
    # e = Int(ceil(N / 4 * 3))
    # e = s + 25
    # y = log.(avgzizjs[s:e, s])
    z = avgzizjs[s:end, s]
    ixs = 1E-6 .< z .< 1E-1
    x = collect(s:N)[ixs]
    y = log.(z[ixs])
    n = length(x)
    slp = (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x .^ 2) - sum(x)^2)
    slopes[j] = slp

    next!(pbar)
end


dat_to_save = Dict(
    "params" => collect(gs),
    "sumzizjs" => sumzizjs,
    "sumzizjs2" => sumzizjs2,
    "slopes" => slopes,
)
save("$datapath/gs-avgzizjs2-slopes.jld2", dat_to_save)



ixs = 7:76
plot(
    gs[ixs], res[ixs],
    xlabel="g",
    ylabel="sumzizjs",
    label="sumzizjs",
    framestyle=:box,
    # widen=false,
    # dpi=500,
)
plot!(
    twinx(),
    gs[ixs], (gs[ixs] .- 1) .^ (-1),
    ylabel="1/(g-1)",
    color=:red,
    label="1/(g-1)"
)
savefig("$figpath/gs-sumzizjs-1.png")


