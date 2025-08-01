

println(@__DIR__)

using Distributed
addprocs(3)
nworkers()

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

@everywhere begin
    include("../src/julia/qising.jl")
end

# include("../src/julia/qising.jl")

##################################################
# fix g=1.5, p=0.4. compute r2 vs. 1/N. (|i-j|=50)

g = 1.5
p = 0.4
bc = "obc"
dN = 10
Ns = 60:dN:200
len_Ns = length(Ns)

renyi2cs = Array{Float64}(undef, len_Ns)
rho12s = Array{Float64}(undef, len_Ns)

pbar = Progress(len_Ns; desc="run 1/N data..")
print(pbar)

for (j, nsize) in enumerate(Ns)
    println(pbar)

    mid = Int(round(nsize / 2))
    s, e = mid - 25, mid + 25
    ry2c, r12 = gen_swssb_single(nsize, g, p, s, e; bc=bc)
    renyi2cs[j] = ry2c
    rho12s[j] = r12

    next!(pbar)
end

dat_to_save = Dict(
    "params" => collect(Ns),
    "renyi2cs" => renyi2cs,
    "rho12s" => rho12s
)
save("$datapath/c1.jld2", dat_to_save)

c1data = load("$datapath/c1.jld2")
Ns = c1data["params"]
renyi2cs = c1data["renyi2cs"]
rho12s = c1data["rho12s"]


plot(
    Ns, renyi2cs,
    xlabel="system size N",
    ylabel="Renyi2 correlator",
    label=L"\frac{Tr(\rho Z_i Z_j \rho Z_i Z_j)}{Tr(\rho^2)}",
    title="Renyi2 correlator\n g=$g p=$p i-j=50 $bc\n",
    linewidth=2,
    framestyle=:box,
    widen=false,
    dpi=500,
);
# plot!(twinx(),
#     Ns, rho12s,
#     label=L"\log\frac{Tr(\rho_1\rho_2)}{Tr(\rho_1^2)}",
#     ylabel="rho12 diff",
#     color=:red
# )
# savefig("$datapath/N1.png")
savefig("$figpath/N.png")
##################################################
