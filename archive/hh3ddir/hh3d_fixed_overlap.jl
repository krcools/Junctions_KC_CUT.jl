# import Pkg
# Pkg.activate(".")

using BEAST, CompScienceMeshes, LinearAlgebra
using JLD2

using Plots
plotly()

using Junctions_KC_CUT

width   = 1.0
height  = 0.5
overlap = 0.2

κ = 1.0; γ = im*κ
HS = Helmholtz3D.hypersingular(gamma=γ)
WS = Helmholtz3D.singlelayer(wavenumber=κ)
N = BEAST.Identity()

uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=normalize(ŷ + ẑ))
e = ∂n(uⁱ)
f = BEAST.ScalarTrace(uⁱ)


nearsys = BEAST.DoubleNumWiltonSauterQStrat(1, 2, 2, 3, 3, 3, 3, 3)
nearpcd = BEAST.DoubleNumWiltonSauterQStrat(1, 1, 1, 2, 2, 2, 2, 2)
farstrat  = BEAST.DoubleNumQStrat(1,2)

smat(op,tfs,bfs) = BEAST.assemble(op,tfs,bfs; quadstrat=nearsys)
pmat(op,tfs,bfs) = BEAST.assemble(op,tfs,bfs; quadstrat=nearpcd)
# hmat(op,tfs,bfs) = AdaptiveCrossApproximation.h1compress(op,tfs,bfs; nearstrat=nearstrat,farstrat=farstrat)
# mat = dmat

hs = [0.1, 0.05, 0.025, 0.0125]
hs = [0.1, 0.05, 0.025]
iters_classic = Int[]
iters_precond = Int[]


# for h in hs
function runsim(;h, κ)
    G1 = meshrectangle(width, height, h)
    G21 = meshrectangle(width, overlap, h)
    G22 = meshrectangle(width, height-overlap, h)
    CompScienceMeshes.translate!(G22, point(0, overlap, 0))
    CompScienceMeshes.rotate!(G21, 1.0π * x̂)
    CompScienceMeshes.rotate!(G22, 1.0π * x̂)
    G3 = CompScienceMeshes.rotate(G1, 1.5π * x̂)

    junction = meshsegment(1.0, 1.0, 3)
    on_junction = overlap_gpredicate(junction)

    G12 = weld(G1,-G2)
    G23 = weld(G2,-G3)
    G31 = weld(G3,-G1);

    ∂G12 = boundary(G12)
    ∂G23 = boundary(G23)
    ∂G31 = boundary(G31)

    V12 = setminus(skeleton(G12,0), skeleton(∂G12, 0))
    V23 = setminus(skeleton(G23,0), skeleton(∂G23, 0))
    V31 = setminus(skeleton(G31,0), skeleton(∂G31, 0))

    # G23_edges = skeleton(G23,1)
    # G23_junction_edges = submesh(c -> on_junction(chart(G23_edges,c)), G23_edges)
    # G23_junction_nodes = skeleton(G23_junction_edges, 0)

    # on_V12 = overlap_gpredicate(V12)
    # in_G23_junction_nodes = in(G23_junction_nodes)
    # V̂23 = submesh(V23) do node
    #     ch = chart(V23, node)
    #     in_G23_junction_nodes(node) && return true
    #     on_V12(ch) && return false
    #     return true
    # end

    # in_V̂23 = in(V̂23)
    # Ĝ23 = submesh(G23) do face
    #     index(face[1]) |> in_V̂23 && return true
    #     index(face[2]) |> in_V̂23 && return true
    #     index(face[3]) |> in_V̂23 && return true
    #     return false
    # end

    X12 = lagrangecxd0(G12)
    # X23 = lagrangecxd0(G23)
    X23 = lagrangecxd0(G3)
    # X31 = lagrangecxd0(G31)

    Y12 = BEAST.duallagrangec0d1(G12, G12, ∂G12)
    # Y23 = BEAST.duallagrangec0d1(G23, G23, ∂G23)
    Y23 = BEAST.duallagrangec0d1(G3, G3, boundary(G3))
    # Y31 = BEAST.duallagrangec0d1(G31, G31, ∂G31)

    X = X12 × X23
    Y = Y12 × Y23

    @hilbertspace j[1:2]
    @hilbertspace k[1:2]

    fx = assemble(@discretise(f[k], k ∈ X))

    WSxx = assemble(@discretise(WS[k,j], j ∈ X, k ∈ X), materialize=smat)

    Nxy = assemble(@discretise(BEAST.diag(N)[j,k], j∈X, k∈Y))
    Dyx = BEAST.GMRESSolver(Nxy, tol=2e-12, restart=250, verbose=false)
    Dxy = BEAST.GMRESSolver(transpose(Nxy), tol=2e-12, restart=250, verbose=false)

    HSyy = assemble(@discretise(BEAST.diag(HS)[k,j], k∈Y, j∈Y), materialize=pmat)
    
    P = Dxy * HSyy * Dyx

    u1, ch1 = solve(BEAST.GMRESSolver(WSxx, tol=2e-5, restart=250), fx)
    u2, ch2 = solve(BEAST.GMRESSolver(P*WSxx, tol=2e-5, restart=250), P*fx)

    @show ch1.iters
    @show ch2.iters

    push!(iters_classic, ch1.iters)
    push!(iters_precond, ch2.iters)
end


# fcr1, geo1 = facecurrents(u1[j[1]], X12)
# fcr2, geo2 = facecurrents(u1[j[2]], X23)

# pt1 = (patch(geo1, (norm.(fcr1))))
# pt2 = (patch(geo2, (norm.(fcr2))))
# Plotly.plot([pt1, pt2])

jldsave("hh3d_strip_dir.jld2"; hs, iters_classic, iters_precond)
@load "hh3d_strip_dir.jld2"

error()

plotly()
plot(hs, (iters_classic))
plot!(hs, (iters_precond))#' Check correctness

phi = 0.0
thetas = range(0,pi,length=200)
farpts = [point(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)) for theta in thetas]

nfp = BEAST.HH3DFarField(wavenumber=κ)
near1 = BEAST.potential(nfp, farpts, u1, X, type=ComplexF64)
near2 = BEAST.potential(nfp, farpts, u2, X, type=ComplexF64)
inc = uⁱ.(farpts)

plot(thetas, abs.(vec(near1)))
scatter!(thetas, abs.(vec(near2)))
 
#' Visualise spectrum

mWSxx = Matrix(WSxx)
mHSyy = Matrix(HSyy)

mDxy = Matrix(Dxy)
mDyx = Matrix(Dyx)

Q = mDxy * mHSyy * mDyx * mWSxx

w = eigvals(Q)
plot(exp.(2pi*im*range(0,1,length=200)))
scatter!(w)

mSxx = BEAST.convert_to_dense(Sxx)
mSyy = BEAST.convert_to_dense(Syy)

Z = iN * mSxx;
W = iN' * mSyy * iN * mSxx;

wZ = eigvals(Matrix(Z))
wW = eigvals(Matrix(W))

plot(exp.(im*range(0,2pi,length=200)))
scatter!(wZ)
scatter!(wW)