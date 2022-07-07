using DrWatson
@quickactivate "Junctions_KC_CUT"

using BEAST, CompScienceMeshes, LinearAlgebra
using JLD2

using Junctions_KC_CUT

using AdaptiveCrossApproximation
AdaptiveCrossApproximation.blockassembler(op,Y,X;quadstrat) = BEAST.blockassembler(op,Y,X;quadstrat)
AdaptiveCrossApproximation.scalartype(op,Y,X) = BEAST.scalartype(op,Y,X)
AdaptiveCrossApproximation.positions(X) = BEAST.positions(X)
AdaptiveCrossApproximation.numfunctions(X) = BEAST.numfunctions(X)

width, height = 1.0, 0.5
overlap = 0.2

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
hmat(op,tfs,bfs) = AdaptiveCrossApproximation.h1compress(op,tfs,bfs;
    nearstrat=nearstrat,farstrat=farstrat)
mat = dmat

phi = pi/2
thetas = range(0,pi,length=200)
farpts = [10*point(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)) for theta in thetas]

function runsim(;h, κ)

    HS = Helmholtz3D.hypersingular(gamma=im*κ)
    WS = Helmholtz3D.singlelayer(wavenumber=κ)
    N = BEAST.Identity()
    uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=normalize(ŷ + ẑ))
    e = ∂n(uⁱ)
    f = BEAST.ScalarTrace(uⁱ)

    G1 = meshrectangle(width, height, h)
    G21 = meshrectangle(width, overlap, h)
    G22 = meshrectangle(width, height-overlap, h)
    CompScienceMeshes.translate!(G22, point(0, overlap, 0))
    CompScienceMeshes.rotate!(G21, 1.0π * x̂)
    CompScienceMeshes.rotate!(G22, 1.0π * x̂)
    G2 = weld(G21, G22)
    G3 = CompScienceMeshes.rotate(G1, 1.5π * x̂)

    # junction = meshsegment(1.0, 1.0, 3)
    # on_junction = overlap_gpredicate(junction)

    G12 = weld(G1,-G2)
    G23 = weld(G2,-G3)
    G31 = weld(G3,-G1);

    ∂G12 = boundary(G12)
    ∂G23 = boundary(G23)
    ∂G31 = boundary(G31)

    # V12 = setminus(skeleton(G12,0), skeleton(∂G12, 0))
    # V23 = setminus(skeleton(G23,0), skeleton(∂G23, 0))
    # V31 = setminus(skeleton(G31,0), skeleton(∂G31, 0))

    Ĝ23 = weld(G21, -G3)
    # V̂23 = setminus(skeleton(Ĝ23,0), skeleton(boundary(Ĝ23),0))

    X12 = lagrangecxd0(G12)
    X23 = lagrangecxd0(Ĝ23)

    Y12 = BEAST.duallagrangec0d1(G12, G12, ∂G12)
    Y23 = BEAST.duallagrangec0d1(Ĝ23, Ĝ23, boundary(Ĝ23))

    X = X12 × X23
    Y = Y12 × Y23

    @hilbertspace j[1:2]
    @hilbertspace k[1:2]

    fx = assemble(@discretise(f[k], k ∈ X))
    WSxx = assemble(@discretise(WS[k,j], j ∈ X, k ∈ X), materialize=mat)
    Nxy = assemble(@discretise(BEAST.diag(N)[j,k], j∈X, k∈Y))
    Dyx = BEAST.GMRESSolver(Nxy, tol=2e-12, restart=250, verbose=false)
    Dxy = BEAST.GMRESSolver(transpose(Nxy), tol=2e-12, restart=250, verbose=false)
    HSyy = assemble(@discretise(BEAST.diag(HS)[k,j], k∈Y, j∈Y), materialize=mat)
    
    P = Dxy * HSyy * Dyx

    u1, ch1 = solve(BEAST.GMRESSolver(WSxx, tol=2e-5, restart=250), fx)
    u2, ch2 = solve(BEAST.GMRESSolver(P*WSxx, tol=2e-5, restart=250), P*fx)

    @show ch1.iters
    @show ch2.iters

    nfp = BEAST.HH3DFarField(wavenumber=κ)
    near1 = BEAST.potential(nfp, farpts, u1, X, type=ComplexF64)
    near2 = BEAST.potential(nfp, farpts, u2, X, type=ComplexF64)

    u1, ch1, u2, ch2, near1, near2, X
end


function makesim(d::Dict)
    @unpack h, κ = d
    u1, ch1, u2, ch2, near1, near2, X = runsim(;h, κ)
    fulld = merge(d, Dict(
        "u1" => u1,
        "u2" => u2,
        "ch1" => ch1,
        "ch2" => ch2,
        "near1" => near1,
        "near2" => near2
    ))
end

method = splitext(basename(@__FILE__))[1]
h = [0.1, 0.05, 0.025, 0.0125]
κ = [1.0, 10.0]

params = @strdict h κ
dicts = dict_list(params)
for (i,d) in enumerate(dicts)
    @show d
    f = makesim(d)
    @tagsave(datadir("simulations", method, savename(d,"jld2")), f)
end

#' Visualise spectrum

# mWSxx = Matrix(WSxx)
# mHSyy = Matrix(HSyy)

# mDxy = Matrix(Dxy)
# mDyx = Matrix(Dyx)

# Q = mDxy * mHSyy * mDyx * mWSxx

# w = eigvals(Q)
# plot(exp.(2pi*im*range(0,1,length=200)))
# scatter!(w)

# mSxx = BEAST.convert_to_dense(Sxx)
# mSyy = BEAST.convert_to_dense(Syy)

# Z = iN * mSxx;
# W = iN' * mSyy * iN * mSxx;

# wZ = eigvals(Matrix(Z))
# wW = eigvals(Matrix(W))
