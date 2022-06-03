# import Pkg
# Pkg.activate(".")

using BEAST, CompScienceMeshes, LinearAlgebra
using JLD2

using Plots
plotly()

include("utils.jl")

width, height = 1.0, 0.5

κ = 1.0; γ = im*κ
HS = Helmholtz3D.hypersingular(gamma=γ)
WS = Helmholtz3D.singlelayer(wavenumber=κ)
N = BEAST.Identity()

uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=normalize(ŷ + ẑ))
e = ∂n(uⁱ)
f = BEAST.ScalarTrace(uⁱ)


nearstrat = BEAST.DoubleNumWiltonSauterQStrat(1, 2, 2, 3, 3, 3, 3, 3)
farstrat  = BEAST.DoubleNumQStrat(1,2)

# BEAST.@defaultquadstrat (WS, X, X)  BEAST.DoubleNumWiltonSauterQStrat(4, 4, -1, -1, 6, 10, 8, 3)
# BEAST.@defaultquadstrat (HS, Y, Y) BEAST.DoubleNumWiltonSauterQStrat(1, 1, -1, -1, 2, 2, 2, 2)

dmat(op,tfs,bfs) = BEAST.assemble(op,tfs,bfs; quadstrat=nearstrat)
# hmat(op,tfs,bfs) = AdaptiveCrossApproximation.h1compress(op,tfs,bfs; nearstrat=nearstrat,farstrat=farstrat)
mat = dmat

hs = [0.1, 0.05, 0.025, 0.0125]
iters_classic = Int[]
iters_precond = Int[]

h = 0.025
# for h in hs
    G1 = meshrectangle(width, height, h)
    G2 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
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

    X12 = lagrangecxd0(G12)
    X23 = lagrangecxd0(G23)
    X31 = lagrangecxd0(G31)

    Y12 = BEAST.duallagrangec0d1(G12, G12, ∂G12)
    Y23 = BEAST.duallagrangec0d1(G23, G23, ∂G23)
    Y31 = BEAST.duallagrangec0d1(G31, G31, ∂G31)

    X = X12 × X23 × X31
    Y = Y12 × Y23 × Y31

    @hilbertspace j[1:3]
    @hilbertspace k[1:3]

    fx = assemble(@discretise(f[k], k ∈ X))

    WSxx = assemble(@discretise(WS[k,j], j ∈ X, k ∈ X))

    Nxy = assemble(@discretise(BEAST.diag(N)[j,k], j∈X, k∈Y))
    Dyx = BEAST.GMRESSolver(Nxy, tol=2e-12, restart=250, verbose=false)
    Dxy = BEAST.GMRESSolver(transpose(Nxy), tol=2e-12, restart=250, verbose=false)

    HSyy = assemble(@discretise(BEAST.diag(HS)[k,j], k∈Y, j∈Y), materialize=mat)
    
    P = Dxy * HSyy * Dyx

    u1, ch1 = solve(BEAST.GMRESSolver(WSxx, tol=2e-5, restart=250), fx)
    u2, ch2 = solve(BEAST.GMRESSolver(P*WSxx, tol=2e-5, restart=250), P*fx)

    @show ch1.iters
    @show ch2.iters

    push!(iters_classic, ch1.iters)
    push!(iters_precond, ch2.iters)
# end

# error()

jldsave("hh3d_fullqs_dir.jld2"; hs, iters_classic, iters_precond)
@load "hh3d_fullqs_dir.jld2"

plotly()
plot(hs, (iters_classic))
plot!(hs, (iters_precond))

error()

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

mSxx = BEAST.convert_to_dense(Sxx)
mSyy = BEAST.convert_to_dense(Syy)

Z = Matrix(WSxx)
W = Matrix(P*WSxx)

# Z = iN * mSxx;
# W = iN' * mSyy * iN * mSxx;

wZ = eigvals(Matrix(Z))
wW = eigvals(Matrix(W))

plot(exp.(im*range(0,2pi,length=200)))
scatter!(wZ)
scatter!(wW)