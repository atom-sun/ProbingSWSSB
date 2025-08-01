
println(@__DIR__)
include("qising.jl")


# run...
@time gen_tising(4, 0.5)
@time gen_tising(100, 0.5)
@time gen_tising(100, 1)
@time gen_tising(100, 1.5)


@time gen_tising(4, 0.5; bc="obc", plotfig=true)
@time gen_tising(100, 0.5; bc="obc", plotfig=true)
@time gen_tising(100, 1; bc="obc", plotfig=true)
@time gen_tising(100, 1.5; bc="obc", plotfig=true)


@time gen_tising(4, 0.5, bc="pbc", plotfig=true)
@time gen_tising(100, 0.5, bc="pbc", plotfig=true)
@time gen_tising(100, 1, bc="pbc", plotfig=true)
@time gen_tising(100, 1.5, bc="pbc", plotfig=true)


gs = 0.0:0.1:2.0;
nsize = 100;

for bc in ["obc", "pbc"]
    zs = []
    for hx in gs
        sites, psi0, energy, avgz, avgzizjs, ee =
            gen_tising(nsize, hx; bc=bc, plotfig=true)
        push!(zs, avgz)
    end

    plot(
        gs, zs,
        xlabel="g",
        ylabel=L"\Delta",
        label=L"\langle Z_j\rangle" * " averaged over j",
        title="Order parameters " * L"\langle Z\rangle" * "\nN=$nsize  $bc\n",
        linewidth=2,
        framestyle=:box,
        widen=false,
        dpi=500,
    )
    savefig("$_FIG_PATH0/orderparameters-N$nsize-$bc.png")
end


# apply kraus quantum channel.
nsize = 4;
g = 1.5;

@time sites, psi0, energy, avgz, avgzizjs, ee =
    gen_tising(nsize, g; bc="obc", plotfig=false)
@time densmat0 = outer(psi0', psi0);

p = 0.1;
@time densmat1 =
    kkraus_channel_zz(densmat0, sites, p; boundary_conditions="open")


# check all expectation values.
oostrs = "I", "X", "Y", "Z"
oosall = []

for j1 in oostrs, j2 in oostrs, j3 in oostrs, j4 in oostrs
    oo = OpSum()
    oo += j1, 1, j2, 2, j3, 3, j4, 4
    oo = MPO(oo, sites)
    push!(oosall, oo)
end


ovals0 = []
ovals1 = []
for j in eachindex(oosall)
    oo = oosall[j]

    v = tr(apply(oo, densmat0))
    v = round(v, digits=12)
    v = convert(Float64, v)
    push!(ovals0, v)

    v = tr(apply(oo, densmat1))
    v = round(v, digits=12)
    v = convert(Float64, v)
    push!(ovals1, v)
end


# check swssb calculation.
@time gen_swssb_single(100, 1.5, 0.4, 25, 75; bc="obc")

