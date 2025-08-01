
println(@__DIR__)

using JLD2
using LaTeXStrings
using Plots


_DATA_PATH = "data/nn_ij"

nsize = 100
bc = "obc"
dis = 1:50

# g: 10, 2
g = 10
ps = (0.1, 0.2, 0.3, 0.4)
ys = []
for p in ps
    cfg_str = "N$nsize-g$g-p$p-$bc"
    data = load("$_DATA_PATH/$cfg_str/data/ij-all.jld2")
    y = data["corr2_all"]
    push!(ys, y)
    @assert data["dis"] == dis
end


x = 2 * dis .- 1
plot(
    x, ys,
    yaxis=:log,
    # yscale=:log10,
    xlim=(-1, 101),
    ylim=(10^-10 * 0.5, 2),
    xlabel=L"|i - j| ",
    ylabel=L"\log C_2",
    label=[L"p=0.1" L"p=0.2" L"p=0.3" L"p=0.4"],
    title="g=$g",
    linewidth=2,
    framestyle=:box,
    widen=false,
    marker=:circle,
    markersize=4,
    dpi=500
)
savefig("$_DATA_PATH/g$g.png")

