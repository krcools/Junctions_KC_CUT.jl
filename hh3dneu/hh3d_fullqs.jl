#  Helmholtz 3D Full Quotient Space Neumann problem

# import Pkg
# Pkg.activate(".")


using BEAST, CompScienceMeshes, LinearAlgebra
using JLD2

include("utils.jl")

using AdaptiveCrossApproximation
AdaptiveCrossApproximation.blockassembler(op,Y,X;quadstrat) = BEAST.blockassembler(op,Y,X;quadstrat)
AdaptiveCrossApproximation.scalartype(op,Y,X) = BEAST.scalartype(op,Y,X)
AdaptiveCrossApproximation.positions(X) = BEAST.positions(X)
AdaptiveCrossApproximation.numfunctions(X) = BEAST.numfunctions(X)

# width, height, h = 1.0, 0.5, 0.025
width, height = 1.0, 0.5

κ = 1.0; γ = im*κ
HS = Helmholtz3D.hypersingular(gamma=γ)
WS = Helmholtz3D.singlelayer(wavenumber=κ)
N = BEAST.Identity()
uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=ẑ)
e = ∂n(uⁱ)


Top = BEAST.Helmholtz3DOp
Tsp = BEAST.LagrangeBasis
Trf = BEAST.LagrangeRefSpace

function BEAST.quaddata(op::Top, tref::Trf, bref::Trf,
    tels, bels, qs::BEAST.DoubleNumQStrat)

    qs = BEAST.DoubleNumWiltonSauterQStrat(qs.outer_rule, qs.inner_rule, 1, 1, 1, 1, 1, 1)
    BEAST.quaddata(op, tref, bref, tels, bels, qs)
end

function BEAST.quadrule(op::Top, tref::Trf, bref::Trf,
    i ,τ, j, σ, qd, qs::BEAST.DoubleNumQStrat)

    return BEAST.DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j])
end

nearstrat = BEAST.DoubleNumWiltonSauterQStrat(1, 2, 2, 3, 3, 3, 3, 3)
farstrat  = BEAST.DoubleNumQStrat(1,2)

dmat(op,tfs,bfs) = BEAST.assemble(op,tfs,bfs; quadstrat=nearstrat)
hmat(op,tfs,bfs) = AdaptiveCrossApproximation.h1compress(op,tfs,bfs; nearstrat=nearstrat,farstrat=farstrat)
mat = hmat

hs = [0.1, 0.05, 0.025, 0.0125]
# hs = [0.1, 0.05, 0.025, 0.0125, 0.00625]
# hs = [0.00625]
iters_classic = Int[]
iters_precond = Int[]

h = 0.1

for h in hs
    G1 = meshrectangle(width, height, h)
    G2 = CompScienceMeshes.rotate(G1, 0.5π * x̂)
    G3 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
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

    X12 = lagrangec0d1(G12, V12)
    X23 = lagrangec0d1(G23, V23)
    X31 = lagrangec0d1(G31, V31)

    Y12 = BEAST.duallagrangecxd0(G12, V12)
    Y23 = BEAST.duallagrangecxd0(G23, V23)
    Y31 = BEAST.duallagrangecxd0(G31, V31)

    X = X12 × X23 × X31
    Y = Y12 × Y23 × Y31

    @hilbertspace j[1:3]
    @hilbertspace k[1:3]

    ex = assemble(@discretise(e[k], k ∈ X))

    HSxx = assemble(@discretise(HS[k,j], j ∈ X, k ∈ X), materialize=mat)

    Nxy = assemble(@discretise(BEAST.diag(N)[j,k], j∈X, k∈Y))
    Dyx = BEAST.GMRESSolver(Nxy, tol=2e-5, restart=250, verbose=false)
    Dxy = BEAST.GMRESSolver(transpose(Nxy), tol=2e-5, restart=250, verbose=false)

    WSyy = assemble(@discretise(BEAST.diag(WS)[j,k], j∈Y, k∈Y), materialize=mat)

    P = Dxy * WSyy * Dyx

    u1, ch1 = solve(BEAST.GMRESSolver(HSxx,tol=2e-5, maxiter=2000, restart=250), ex)
    u2, ch2 = solve(BEAST.GMRESSolver(P*HSxx, tol=2e-5, maxiter=2000, restart=250), P*ex)

    @show ch1.iters
    @show ch2.iters

    push!(iters_classic, ch1.iters)
    push!(iters_precond, ch2.iters)
end

jldsave("hh3d_fullqs.jld2"; hs, iters_classic, iters_precond)
@load "hh3d_fullqs.jld2"

plotly()
plot(hs, log10.(iters_classic))
plot!(hs, log10.(iters_precond))

error()

phi = pi/2
thetas = range(0,pi,length=200)
farpts = [10*point(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)) for theta in thetas]

nfp = BEAST.HH3DDoubleLayerNear(wavenumber=κ)
near1 = BEAST.potential(nfp, farpts, u1, X, type=ComplexF64)
near2 = BEAST.potential(nfp, farpts, u2, X, type=ComplexF64)
inc = uⁱ.(farpts)
# tot = near - inc

using Plots
# contour(zs,ys,imag.(near-inc),levels=50)
plot(thetas, abs.(vec(near1)))
scatter!(thetas, abs.(vec(near2)))


mHSxx = BEAST.convert_to_dense(HSxx)
mWSyy = BEAST.convert_to_dense(WSyy)
mDyx = Matrix(Dyx)
mDxy = copy(mDyx')

Z = BEAST.convert_to_dense(HSxx)
W = mDxy * mWSyy * mDyx * mHSxx;

wZ = eigvals(Matrix(Z))
wW = eigvals(Matrix(W))

plotly()
plot(exp.(im*range(0,2pi,length=200)), label=nothing)
scatter!(wZ, label="classic")
scatter!(wW, label="fully reduced MT")