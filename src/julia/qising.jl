
using FileIO
using Plots
using ITensors, ITensorMPS
using LaTeXStrings


_FIG_PATH0 = "run/figs"


# Hamiltonian
# written in terms of Pauli matrices instead of Sz/Sx.
function transising_mpo(nsize, hx; boundary_conditions="open")
    @assert (boundary_conditions in ["open", "periodic"]) "illegal boundary conditions"

    sites = siteinds("S=1/2", nsize)

    os = OpSum()
    for j = 1:nsize-1
        os += "Z", j, "Z", j + 1
        os += hx, "X", j
    end
    os += hx, "X", nsize
    if boundary_conditions == "periodic"
        os += "Z", nsize, "Z", 1
    end
    os = -1.0 * os

    H = MPO(os, sites)
    return H, sites
end


# config string. 
function config_str(cfg)

    # extract parameters from cfg. 
    nsize = cfg["nsize"]
    hx = cfg["hx"]
    nsweeps = cfg["nsweeps"]
    bonddims = cfg["bonddims"]
    ctf = cfg["cutoff"]
    bc = cfg["boundary_conditions"]

    rp = Dict("." => "_")
    cfg_short_name = replace("N$nsize-g$hx", rp...)
    # cfg_long_name = ["$k$v" for (k, v) in pairs(cfg)]
    cfg_long_name = replace(
        "N$nsize-g$hx-$bc-nsweeps$nsweeps-bonddim$(bonddims[end])-cutoff$ctf",
        rp...)
    figpth = "$_FIG_PATH0/$cfg_long_name"
    if !isdir(figpth)
        mkpath(figpth)
    end
    return figpth, cfg_short_name, cfg_long_name
end


# entanglement entropy
function seentropy(psi, b::Integer; ignore=1E-10)
    psi = orthogonalize(psi, b)
    U, S, V = svd(psi[b], (linkinds(psi, b - 1)..., siteinds(psi, b)...))
    svn = 0.0
    sry2 = 0.0
    for n = 1:dim(S, 1)
        p = S[n, n]^2
        if abs(p) < ignore
            continue
        end
        svn -= p * log(p)
        sry2 += p^2
    end
    sry2 = -log(sry2)
    return svn, sry2
end


function plot_ee(psi, cfg; plotfig=true)
    nsize = cfg["nsize"]

    ee = [[seentropy(psi, j)...] for j in 1:nsize]
    ee = [[[0, 0]]; ee]
    ee = reduce(hcat, ee)
    ee = ee'

    if !plotfig
        return ee
    end

    hx = cfg["hx"]
    bc = cfg["boundary_conditions"]
    figpth, cfg_short_name, _ = config_str(cfg)

    plot(
        0:nsize, ee,
        xlabel="subsystem size A",
        ylabel=L"S_A",
        label=[L"S_{vN}" L"S_{Renyi2}"],
        title="Entanglement Entropy\n N=$nsize g=$hx $bc\n",
        linewidth=2,
        framestyle=:box,
        widen=false,
        dpi=500,
    )
    savefig("$figpth/$cfg_short_name-ee.png")
    return ee
end


# order parameters.  <Zj>. 
function order_params_z(psi)
    avgzjs = expect(psi, "Z")
    avgzj = abs(sum(avgzjs))
    return avgzjs, avgzj
end


# correlation functions.  <ZiZj> - <Zi><Zj>. 
function plot_corr_z(psi, cfg; plotfig=true)
    avgzizjs = correlation_matrix(psi, "Z", "Z")
    avgzjs, avgz = order_params_z(psi)
    corrzizjs = avgzizjs - avgzjs * avgzjs'

    if !plotfig
        return avgzizjs, avgz
    end

    nsize = cfg["nsize"]
    hx = cfg["hx"]
    bc = cfg["boundary_conditions"]
    figpth, cfg_short_name, _ = config_str(cfg)

    s, e = Int(floor(nsize / 4)), Int(ceil(nsize / 4 * 3))
    s = max(1, s)
    e = min(length(psi), e)

    plot(
        s:e, corrzizjs[s:e, s],
        xlabel="site distance " * L"|i-j|",
        ylabel=L"G_{ij}",
        label=L"\langle Z_iZ_j\rangle - \langle Z_i\rangle\langle Z_j\rangle",
        title="Correlation functions ZiZj\n N=$nsize g=$hx $bc\n",
        linewidth=2,
        framestyle=:box,
        widen=false,
        dpi=500,
    )
    savefig("$figpth/$cfg_short_name-corr.png")

    plot(
        s:e, corrzizjs[s:e, s],
        yaxis=:log,
        xlabel="site distance " * L"|i-j|",
        ylabel=L"\log G_{ij}",
        label=L"\langle Z_iZ_j\rangle - \langle Z_i\rangle\langle Z_j\rangle",
        title="Correlation functions ZiZj\n N=$nsize g=$hx $bc\n",
        linewidth=2,
        framestyle=:box,
        widen=false,
        dpi=500,
    )
    savefig("$figpth/$cfg_short_name-corrlog.png")

    plot(
        s:e, corrzizjs[s:e, s],
        xaxis=:log,
        yaxis=:log,
        xlabel="log site distance:  " * L"\log|i-j|",
        ylabel=L"\log G_{ij}",
        label=L"\langle Z_iZ_j\rangle - \langle Z_i\rangle\langle Z_j\rangle",
        title="Correlation functions ZiZj\n N=$nsize g=$hx $bc\n",
        linewidth=2,
        framestyle=:box,
        widen=false,
        dpi=500,
    )
    savefig("$figpth/$cfg_short_name-corrloglog.png")

    return avgzizjs, avgz
end


function gen_tising(
    nsize, hx;
    bc="obc",
    nsweeps=10, bonddims=[10, 20, 100, 200, 400, 800], ctf=1E-13,
    plotfig=false,
)

    cfg = Dict(
        "nsize" => nsize,
        "hx" => hx,
        "boundary_conditions" => bc,
        "nsweeps" => nsweeps,
        "bonddims" => bonddims,
        "cutoff" => ctf,
    )

    # Hamiltonian and physical sites.
    bcdict = Dict("pbc" => "periodic", "obc" => "open")
    H, sites = transising_mpo(nsize, hx; boundary_conditions=bcdict[bc])

    # Prepare initial state MPS
    state = ["Up" for n = 1:nsize]
    # state = [isodd(n) ? "Up" : "Dn" for n=1:nsize]
    psi0_i = MPS(sites, state)

    # Do 10 sweeps of DMRG, gradually
    # increasing the maximum MPS
    # bond dimension
    sweeps = Sweeps(nsweeps)
    setmaxdim!(sweeps, bonddims...)  # bond dimension
    setcutoff!(sweeps, ctf)  # cutoff

    # Run the DMRG algorithm
    energy, psi0 = dmrg(H, psi0_i, sweeps)

    # show results.
    avgzizjs, avgz = plot_corr_z(psi0, cfg; plotfig=plotfig)
    ee = plot_ee(psi0, cfg; plotfig=plotfig)

    return sites, psi0, energy, avgz, avgzizjs, ee
end


function kkraus_channel_zz(densmat, sites, p;
    boundary_conditions="open", ctf=1E-13,
)

    """
    KKrause quantum channel ZZ is defined as below.

    KK_{i,i+1} = (1 - p) \rho + p Z_iZ_j \rho Z_iZ_j
    kraus1 = sqrt(1-p) I
    kraus2 = sqrt(p)ZiZj
    Apply on (i,i+1): kraus1 \rho kraus1' + kraus2 \rho kraus2'
    sum over all (i,i+1).

    """

    @assert (boundary_conditions in ["open", "periodic"]) "illegal boundary conditions"
    if boundary_conditions == "periodic"
        @warn "BE CAUTIOUS WITH PBC! siteinds/linkinds be wrong itentified!"
    end

    nsize = length(densmat)

    @assert nsize == length(sites) "sites not match with density matrix"
    @assert abs(tr(densmat) - 1) < 1E-3 "density matrix trace not 1"
    # TODO: test rho hermittian.
    # rho should have proper site inds for trace calculation.
    @assert siteinds(densmat) == [[s', s] for s in sites] "improper siteinds of rho"

    # define two-qubit gates operation order.
    # eg. [(1, 2), (3, 4), (5, 6), (2, 3), (4, 5)], for N=6, obc.
    gsidx = vcat(collect(zip(1:2:nsize, 2:2:nsize)),
        collect(zip(2:2:nsize, 3:2:nsize)))
    if boundary_conditions == "periodic"
        gsidx = vcat(gsidx, [(nsize, 1)])
    end

    # apply k1 \rho k1 + k2 \rho k2
    for (j1, j2) in gsidx

        # initialize for each pair of sites
        ss = [sites[j1], sites[j2]]
        densmat = orthogonalize(densmat, j1)
        dm = densmat[j1] * densmat[j2]  # l1,j1',j1,j2',j2,l2

        # define k1 = sqrt(1 - p)
        os = OpSum()
        os += sqrt(1 - p), "I", 1, "I", 2
        k1 = MPO(os, ss)
        gate1 = k1[1] * k1[2]  # j1',j1,j2',j2

        # apply k1 dm k1
        wf1 = prime(gate1, 1)  # j1'',j1',j2'',j2'
        wf1 = wf1 * dm  # l1,j1'',j1,j2'',j2,l2
        wf1 = prime(wf1, 1, ss)  # l1,j1'',j1',j2'',j2',l2
        wf1 = wf1 * gate1  # l1,j1'',j1,j2'',j2,l2
        wf1 = prime(wf1, -1, ss'')  # l1,j1',j1,j2',j2,l2

        # define k2 = sqrt(p)ZiZj
        os = OpSum()
        os += sqrt(p), "Z", 1, "Z", 2
        k2 = MPO(os, ss)
        gate2 = k2[1] * k2[2]  # j1',j1,j2',j2

        # apply k2 dm k2
        wf2 = prime(gate2, 1)  # j1'',j1',j2'',j2'
        wf2 = wf2 * dm  # l1,j1'',j1,j2'',j2,l2
        wf2 = prime(wf2, 1, ss)  # l1,j1'',j1',j2'',j2',l2
        wf2 = wf2 * gate2  # l1,j1'',j1,j2'',j2,l2
        wf2 = prime(wf2, -1, ss'')  # l1,j1',j1,j2',j2,l2

        wf = wf1 + wf2

        idx = uniqueinds(densmat[j1], densmat[j2])  # l1,j1',j1
        U, S, V = svd(wf, idx; cutoff=ctf)

        # update density matrix
        densmat[j1] = U
        densmat[j2] = S * V

    end

    # BE CAUTIOUS WITH PBC!
    # siteinds/linkinds may be wrong identified!!! (bug with ITensors)
    # TODO: deal with periodic boundary_conditions

    return densmat

end


function strong_weak_sym_operation_x(rho::MPO, sites)
    """Strong/weak symmetric operation XXXXX on density matrix rho.

      * Strong: operate XXXXX (all X) from one side.

      * Weak: operate XXXXX from each sides.

    """

    @assert length(sites) == length(rho) "sites not match with density matrix"
    @assert abs(tr(rho) - 1) < 1E-3 "density matrix trace not 1"
    # TODO: test rho hermittian.

    rho_strong = copy(rho)
    rho_weak = copy(rho)

    for (j, (sj, r1)) in enumerate(zip(sites, rho))
        ss = sites[j:j]
        os = OpSum() + ("X", 1)
        gate = MPO(os, ss)
        ox = gate[1]

        r2 = prime(ox, 1) * r1
        r3 = prime(r2, 1, ss) * ox

        r2 = prime(r2, -1, ss'')
        rho_strong[j] = r2

        r3 = prime(r3, -1, ss'')
        rho_weak[j] = r3
    end

    return rho_strong, rho_weak
end


function renyi2corr_zizj(rho::MPO, sites, i::Int64, j::Int64)
    """Calculate renyi2 correlator with charge operator Z.

        Tr(rho ZiZj rho ZiZj) / Tr(rho^2)

    """

    @assert length(sites) == length(rho) "sites not match with density matrix"
    @assert abs(tr(rho) - 1) < 1E-3 "density matrix trace not 1"
    # TODO: test rho hermittian.
    @assert i > 0 && j > 0 "i and j must be positive"
    @assert max(i, j) <= length(rho) "i, j exceed system size"

    # rho should have proper site inds for trace calculation.
    # eg. PBC channel can produce improper site inds of rho.
    @assert siteinds(rho) == [[s', s] for s in sites] "improper siteinds of rho"

    # DON'T DO NON-LOCAL OPERATIONS ON MPS/MPO !!!
    # such as Zi Zj when |i-j|->inf.
    # apply Zi Zj order-by-order.

    rz = copy(rho)
    for k_site in (j, i)
        ss = sites[k_site:k_site]
        os = OpSum() + ("Z", 1)
        gate = MPO(os, ss)
        oz = gate[1]
        rz[k_site] = prime(prime(oz, 1) * rz[k_site], -1, ss'')
    end

    # compute renyi2 correlator.
    # c2 = inner(rz, rz) / inner(rho, rho)  # wrong!!!
    c2 = tr(apply(rz, rz)) / inner(rho, rho)
    # @assert c2 <= 1.01 "renyi2corr should not be larger than 1."

    # compute conventional correlator.
    c0 = tr(rz)

    return c2, c0
end


function prepare_groundstate(nsize, hx; bc="obc", plotfig=false, kwargs...)
    sites, psi0, energy, avgz, avgzizjs, ee =
        gen_tising(nsize, hx; bc=bc, plotfig=plotfig, kwargs...)
    densmat0 = outer(psi0', psi0)
    return densmat0, sites
end


function qchannel_and_symmetric_operation(densmat0, sites, dec_strength)
    densmat1 = kkraus_channel_zz(densmat0, sites, dec_strength)
    densmat2, densmat3 = strong_weak_sym_operation_x(densmat1, sites)
    return densmat1, densmat2, densmat3
end


function gen_swssb_single(
    nsize, hx, dec_strength, i_site, j_site;
    bc="obc",
    plotfig=false,
    kwargs...
)
    """
    hx: transverse field strength g.
    dec_strength: decoherence strength p of Krause quantum channel.
    i_site, j_site: site labels to compute renyi2 correlator.

    """
    densmat0, sites =
        prepare_groundstate(nsize, hx; bc=bc, plotfig=plotfig, kwargs...)

    # quantum channel default use OBC. Problematic with PBC.
    densmat1, densmat2, densmat3 =
        qchannel_and_symmetric_operation(densmat0, sites, dec_strength)

    ry2c, c0 = renyi2corr_zizj(densmat1, sites, i_site, j_site)
    rho12 = inner(densmat2, densmat1) / inner(densmat1, densmat1)
    rho12 = log(rho12)

    return ry2c, c0, rho12
end

