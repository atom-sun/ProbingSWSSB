
##################################################
# fix N, g, p. compute Renyi2 correlator c2 vs. |i-j|. 

@everywhere begin
    println(@__DIR__)
    include("qising.jl")
    println("threads:    ", Threads.nthreads())

    using Distributed
    using JLD2

    # absolute path
    _DATA_PATH0 = "/public1/home/m8s000158/swssb/julia/dump"

    println("nworkers:  ", nworkers())

    # g: 10, 2
    # p: 0.1, 0.2, 0.3, 0.4
    g = 10
    p = 0.4
    nsize = 100
    bc = "obc"
    cfg_str = "N$nsize-g$g-p$p-$bc"

    path0 = "$_DATA_PATH0/$cfg_str"
    datapath = "$path0/data"
    if !isdir(datapath)
        mkpath(datapath)
    end
    figpath = "$path0/figs"
    if !isdir(figpath)
        mkpath(figpath)
    end

    midN = Int(round(nsize / 2))
    dis = 1:midN
    len_dis = length(dis)
    js1 = [Int(midN - di + 1) for di in dis]
    js2 = [Int(midN + di) for di in dis]

    densmat0, sites = prepare_groundstate(nsize, g; bc=bc)
    densmat1, densmat2, densmat3 = qchannel_and_symmetric_operation(densmat0, sites, p)
    rho12 = inner(densmat2, densmat1) / inner(densmat1, densmat1)
    rho12 = log(rho12)
    println("rho12 diff:  ", rho12)
end


# parallel run c2.
@sync @distributed for (j, (s, e)) in collect(enumerate(zip(js1, js2)))
    println("$j:   start..")
    c2, c0 = renyi2corr_zizj(densmat1, sites, s, e)
    save("$datapath/$j.jld2", "$j", [c2, c0])
    println("$j job ($s -- $e):  c2:  $c2,   c0:  $c0")
    println("$j:   ($s, $e)  done..")
end


# dump data.
using ProgressMeter
pbar = Progress(len_dis; desc="run 1/|i-j| data..")
print(pbar)

corr2_all = Array{Float64}(undef, len_dis)
corr0_all = Array{Float64}(undef, len_dis)

for (j, (s, e)) in enumerate(zip(js1, js2))
    println("$j:   start..")
    dat = load("$datapath/$j.jld2")
    c2, c0 = dat["$j"]
    corr2_all[j] = c2
    corr0_all[j] = c0
    println("$j job ($s -- $e):  c2:  $c2,   c0:  $c0")
    println("$j:   ($s, $e)  done..")
    next!(pbar)
end

data = Dict(
    "dis" => collect(dis),
    "corr2_all" => corr2_all,
    "corr0_all" => corr0_all,
    "rho12" => rho12,
)
save("$datapath/ij-all.jld2", data)


# load data.
data = load("$datapath/ij-all.jld2")
dis = data["dis"]
corr2_all = data["corr2_all"]
corr0_all = data["corr0_all"]
rho12 = data["rho12"]


# plots.
using Plots

x = 2 * collect(dis) .- 1
plot(
    x, corr2_all,
    xlabel="site distance |i - j| ",
    ylabel="Renyi2 correlator",
    label=L"\frac{Tr(\rho Z_i Z_j \rho Z_i Z_j)}{Tr(\rho^2)}",
    title="Renyi2 correlator\n g=$g p=$p $nsize $bc\n",
    linewidth=2,
    framestyle=:box,
    widen=false,
    dpi=500,
)
savefig("$figpath/ij.png")


plot(
    x, corr2_all,
    yaxis=:log,
    xlabel="site distance |i - j| ",
    ylabel=L"\log " * "Renyi2",
    label=L"\log\frac{Tr(\rho Z_i Z_j \rho Z_i Z_j)}{Tr(\rho^2)}",
    title="Renyi2 correlator\n g=$g p=$p $nsize $bc\n",
    linewidth=2,
    framestyle=:box,
    widen=false,
    dpi=500,
)
savefig("$figpath/ij-log.png")

##################################################

