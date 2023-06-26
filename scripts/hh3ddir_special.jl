using DrWatson
@quickactivate "Junctions_KC_CUT"

using BEAST, CompScienceMeshes, LinearAlgebra
using JLD2

using Junctions_KC_CUT




function runsim(;h, κ)

    HS = Helmholtz3D.hypersingular(gamma=κ*im)
    WS = Helmholtz3D.singlelayer(wavenumber=κ)
    N = BEAST.Identity()
    uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=normalize(x̂ + ŷ + ẑ))
    # e = ∂n(uⁱ)
    f = BEAST.ScalarTrace(uⁱ)

    width = 1.0
    I1 = meshrectangle(width, width/2, h)
    I2 = CompScienceMeshes.translate(I1, -width*x̂)
    I3 = CompScienceMeshes.translate(I1, -width*x̂ - 0.5*width*ŷ)
    I4 = CompScienceMeshes.translate(I1, -0.5*width*ŷ)
    I5 = CompScienceMeshes.rotate(I1, -pi/2*ŷ)

    G1 = weld(I5, I2, I3, I4)
    G2 = weld(-I5, I1, I4, I3)
    G3 = weld(-I1, -I2, -I3, -I4)

    G = weld(I1,I2,I3,I4,I5)

    @assert CompScienceMeshes.isoriented(G1)
    @assert CompScienceMeshes.isoriented(G2)
    @assert CompScienceMeshes.isoriented(G3)

    ∂G1 = boundary(G1)
    ∂G2 = boundary(G2)
    ∂G3 = boundary(G3)

    V1 = setminus(skeleton(G1,0), skeleton(∂G1,0))
    V2 = setminus(skeleton(G2,0), skeleton(∂G2,0))
    V3 = setminus(skeleton(G3,0), skeleton(∂G3,0))

    X1 = lagrangecxd0(G1)
    X2 = lagrangecxd0(G2)
    X3 = lagrangecxd0(G3)

    Y1 = BEAST.duallagrangec0d1(G1, G1, ∂G1)
    Y2 = BEAST.duallagrangec0d1(G2, G2, ∂G2)
    Y3 = BEAST.duallagrangec0d1(G3, G3, ∂G3)

    X = X1 × X2 × X3
    Y = Y1 × Y2 × Y3

    @hilbertspace j[1:3]
    @hilbertspace k[1:3]

    fx = assemble(@discretise(f[k], k ∈ X))

    Nxy = assemble(@discretise(BEAST.diag(N)[j,k], j∈X, k∈Y))
    Dyx = BEAST.GMRESSolver(Nxy, tol=2e-12, restart=250, verbose=false)
    Dxy = BEAST.GMRESSolver(transpose(Nxy), tol=2e-12, restart=250, verbose=false)

    WSxx = assemble(@discretise(WS[k,j], j ∈ X, k ∈ X))
    HSyy = assemble(@discretise(BEAST.diag(HS)[k,j], k∈Y, j∈Y))

    u1, ch1 = solve(BEAST.GMRESSolver(WSxx, tol=2e-5, restart=250), fx)

    P = Dxy * HSyy * Dyx
    u2, ch2 = solve(BEAST.GMRESSolver(P*WSxx, tol=2e-5, restart=250), P*fx)

    xs = range(-1.5*width, 1.5*width, length=50) .+ 0.000127
    zs = range(-0.5*width, 1.5*width, length=50) .+ 0.000127
    pts = [point(x,width/4,z) for x in xs, z in zs]

    # scat1 = potential(BEAST.HH3DHyperSingularNear(wavenumber=κ), pts, u1, X)
    # scat2 = potential(BEAST.HH3DHyperSingularNear(wavenumber=κ), pts, u2, X)

    # ∇uⁱ = BEAST.gradient(uⁱ)
    # inc = ∇uⁱ.(pts)
    # near1 = scat1 .+ inc
    # near2 = scat2 .+ inc
    near1 = nothing
    near2 = nothing

    @show ch1.iters
    @show ch2.iters

    return (
        u1=u1, ch1=ch1.iters, near1=near1,
        u2=u2, ch2=ch2.iters, near2=near2)
end

h = [0.1, 0.075, 0.05, 0.025]
h = [0.015]
κ = [10.0]

simname = splitext(basename(@__FILE__))[1]
configs = dict_list(@dict h κ)
for config in configs
    results = runsim(;config...)
    results = tostringdict(merge(config, ntuple2dict(results)))
    @tagsave(datadir("simulations", simname, savename(config,"jld2")), results)
end



# error()

# using Plotly
# Gplt = G3 
# plt1 = patch(Gplt, opacity=0.5)
# plt2 = CompScienceMeshes.wireframe(Gplt)
# Plotly.plot([plt1, plt2])

# Plots.heatmap(zs, xs, real.(getindex.(near1,1)), clims=(-15,15))
# Plots.heatmap(zs, xs, real.(getindex.(near1,3)), clims=(-15,15))

# plt_pts = Plotly.scatter3d(
#     x=vec(getindex.(pts,1)),
#     y=vec(getindex.(pts,2)),
#     z=vec(getindex.(pts,3)))
# Plotly.Plot([plt1, plt_pts])

# @show xs[40]
# Plots.plot(real.(getindex.(near1[40,:],3)), ylims=(-15,15))
# Plots.scatter!(real.(getindex.(near2[40,:],3)), ylims=(-15,15))

# @show ys[60]
# Plots.plot(real.(getindex.(near1[60,:],3)), ylims=(-15,15))
# Plots.scatter!(real.(getindex.(near2[60,:],3)), ylims=(-15,15))

# @show zs[50]
# Plots.plot(real.(getindex.(near1[:,50],1)), ylims=(-15,15))
# Plots.scatter!(real.(getindex.(near2[:,50],1)), ylims=(-15,15))

# Visualise the geometry
