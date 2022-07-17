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

Top = BEAST.MWSingleLayer3D
Tsp = BEAST.RTRefSpace
Trf = BEAST.RTRefSpace

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

    SL = Maxwell3D.singlelayer(wavenumber=κ)
    N = NCross()
    E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
    e = (n × E) × n;

    G1 = meshrectangle(width, height, h)
    G2 = CompScienceMeshes.rotate(G1, 0.5π * x̂)
    G3 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
    junction = meshsegment(1.0, 1.0, 3)
    on_junction = overlap_gpredicate(junction)

    G1_edges = skeleton(G1,1)
    junction_edges = submesh(G1_edges) do edge
        ch = chart(G1_edges, edge)
        on_junction(ch) && return true
        return false
    end
    junction_nodes = skeleton(junction_edges,0)
    on_junction_nodes = overlap_gpredicate(junction_nodes)

    G12 = weld(G1,-G2)
    G23 = weld(G2,-G3)
    G31 = weld(G3,-G1)
    G123 = weld(G1,G2,G3)

    edges12_all = skeleton(G12,1)
    edges23_all = skeleton(G23,1)
    on_bnd23 = in(boundary(G23))

    on_G12 = overlap_gpredicate(edges12_all)
    edges23_act = submesh(edges23_all) do edge
        ch = chart(edges23_all,edge)
        v1 = simplex((ch.vertices[1],), Val{0})
        v2 = simplex((ch.vertices[2],), Val{0})
        on_bnd23(edge) && return false
        on_junction_nodes(v1) && return true
        on_junction_nodes(v2) && return true
        on_G12(ch) && return false
        return true
    end

    in_edges23_act = in(edges23_act)
    Ĝ23 = submesh(G23) do face
        index(face[1],face[2]) |> in_edges23_act && return true
        index(face[2],face[3]) |> in_edges23_act && return true
        index(face[3],face[1]) |> in_edges23_act && return true
        return false
    end

    X12 = raviartthomas(G12)
    X23 = raviartthomas(Ĝ23, edges23_act)

    Y12 = buffachristiansen(G12)
    Y23 = buffachristiansen(Ĝ23, edges=edges23_act)


    X = X12 × X23
    Y = Y12 × Y23

    @show numfunctions(X)
    @show numfunctions(Y)


    @hilbertspace p[1:2]
    @hilbertspace q[1:2]
    Skj = @discretise SL[p,q] p∈X q∈X
    Sxx = assemble(Skj.bilform, Skj.test_space_dict, Skj.trial_space_dict; materialize=dmat)
    ex = BEAST.assemble(@discretise e[p] p∈X)
    
    @hilbertspace j[1:2]
    @hilbertspace k[1:2]

    Nxy = assemble(@discretise(BEAST.diag(N)[j,k], j∈X, k∈Y))
    Dyx = BEAST.GMRESSolver(Nxy, tol=2e-12, restart=250, verbose=false)
    Dxy = BEAST.GMRESSolver(transpose(Nxy), tol=2e-12, restart=250, verbose=false)

    diagSkj = @discretise BEAST.diag(SL)[p,q] p∈Y q∈Y
    Syy = BEAST.assemble(diagSkj; materialize=dmat);

    # P = iN' * Syy * iN
    P = Dxy * Syy * Dyx

    u1, ch1 = solve(BEAST.GMRESSolver(Sxx,tol=2e-5, restart=250), ex)
    u2, ch2 = solve(BEAST.GMRESSolver(P*Sxx, tol=2e-5, restart=250), P*ex)

    @show ch1.iters
    @show ch2.iters

    Φ, Θ = [0.0], range(0,stop=π,length=100)
    pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]

    near1 = potential(MWFarField3D(wavenumber=κ), pts, u1, X)
    near2 = ffd_strip_pc = potential(MWFarField3D(wavenumber=κ), pts, u2, X)

    u1, ch1.iters, u2, ch2.iters, near1, near2, X
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

#' Visualise the spectrum

# mSxx = BEAST.convert_to_dense(Sxx)
# mSyy = BEAST.convert_to_dense(Syy)

# Z = mSxx
# W = iN' * mSyy * iN * mSxx;

# wZ = eigvals(Matrix(Z))
# wW = eigvals(Matrix(W))

# plot(exp.(im*range(0,2pi,length=200)))
# scatter!(wZ)
# scatter!(wW)


# Study the various kernels
# HS = Maxwell3D.singlelayer(gamma=0.0, alpha=0.0, beta=1.0)
# Id = BEAST.Identity()

# Z12 = BEAST.lagrangecxd0(G12)
# Z23 = BEAST.lagrangecxd0(Ĝ23)
# Z = Z12 × Z23

# W12 = BEAST.duallagrangecxd0(G12)
# W23 = BEAST.duallagrangecxd0(Ĝ23)
# W = W12 × W23

# DX = assemble(Id, Z, divergence(X))
# HX = assemble(HS, X, X)

# DY = assemble(Id, W, divergence(Y))
# HY = assemble(HS, Y, Y)

# Nx = BEAST.NCross()
# NYX = assemble(Nx, Y, X)

# Q = HY * iN * HX