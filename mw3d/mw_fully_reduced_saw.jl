import Pkg
Pkg.activate(".")

using BEAST, CompScienceMeshes, LinearAlgebra
using JLD2

using Plots
plotly()

include("utils.jl")

using AdaptiveCrossApproximation
AdaptiveCrossApproximation.blockassembler(op,Y,X;quadstrat) = BEAST.blockassembler(op,Y,X;quadstrat)
AdaptiveCrossApproximation.scalartype(op,Y,X) = BEAST.scalartype(op,Y,X)
AdaptiveCrossApproximation.positions(X) = BEAST.positions(X)
AdaptiveCrossApproximation.numfunctions(X) = BEAST.numfunctions(X)

width, height = 1.0, 0.5

κ = 1.0
SL = Maxwell3D.singlelayer(wavenumber=κ)
N = NCross()

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n;

Tp = typeof(SL)
X0 = raviartthomas(meshsphere(radius=1.0, h=0.35))
x0 = BEAST.refspace(X0)
function BEAST.quaddata(op::Tp, tref::typeof(x0), bref::typeof(x0),
    tels, bels, qs::BEAST.DoubleNumQStrat)

    qs = BEAST.DoubleNumWiltonSauterQStrat(qs.outer_rule, qs.inner_rule, 1, 1, 1, 1, 1, 1)
    BEAST.quaddata(op, tref, bref, tels, bels, qs)
end
function BEAST.quadrule(op::Tp, tref::typeof(x0), bref::typeof(x0),
    i ,τ, j, σ, qd, qs::BEAST.DoubleNumQStrat)

    return BEAST.DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j])
end

nearstrat = BEAST.DoubleNumWiltonSauterQStrat(1, 2, 3, 4, 3, 3, 3, 3)
farstrat  = BEAST.DoubleNumQStrat(1,2)
hmat(op,tfs,bfs) = AdaptiveCrossApproximation.h1compress(op,tfs,bfs;nearstrat=nearstrat,farstrat=farstrat)
dmat(op,tfs,bfs) = BEAST.assemble(op,tfs,bfs; quadstrat=nearstrat)

# hs = [0.2, 0.1, 0.05, 0.025, 0.0125]
hs = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625]
# hs = [0.2, 0.1, 0.05]
iters_classic = Int[]
iters_sawtooth = Int[]

h = last(hs)
# for h in hs
    G1 = meshrectangle(width, height, h)
    G2 = CompScienceMeshes.rotate(G1, 0.5π * x̂)
    G3 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
    junction = meshsegment(1.0, 1.0, 3)

    G12 = weld(G1,-G2)
    G23 = weld(G2,-G3)
    G31 = weld(G3,-G1);

    edges12_all = skeleton(G12,1)
    edges23_all = skeleton(G23,1)
    on_bnd23 = in(boundary(G23))
    on_junction = overlap_gpredicate(junction)
    on_G12 = overlap_gpredicate(edges12_all)
    edges23_act = submesh(edges23_all) do edge
        ch = chart(edges23_all,edge)
        on_bnd23(edge) && return false
        on_junction(ch) && return true
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
    Y = Y12 × Y23;


    # @hilbertspace j
    # @hilbertspace k
    # mtefie = @discretise(
    #     SL[k,j] == e[k],
    #     j ∈ X, k ∈ X);

    @hilbertspace p[1:2]
    @hilbertspace q[1:2]
    Skj = @discretise SL[p,q] p∈X q∈X
    Sxx = assemble(Skj.bilform, Skj.test_space_dict, Skj.trial_space_dict; materialize=hmat)        

    # BEAST.@defaultquadstrat (SL,X12,X12) BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(2, 3, 6, 7, 5, 5, 4, 3)
    # Sxx = BEAST.sysmatrix(mtefie)
    # ex = BEAST.rhs(mtefie)
    ex = BEAST.assemble(@discretise e[p] p∈X)

    # N1 = assemble(N, X12, Y12)
    # N2 = assemble(N, X23, Y23)
    # iN = blkdiagm(inv(Matrix(N1)), inv(Matrix(N2)))

    @hilbertspace j[1:2]
    @hilbertspace k[1:2]
    Nxy = assemble(@discretise(BEAST.diag(N)[j,k], j∈X, k∈Y))

    Dyx = BEAST.GMRESSolver(Nxy, tol=2e-5, restart=250, verbose=false)
    Dxy = BEAST.GMRESSolver(transpose(Nxy), tol=2e-5, restart=250, verbose=false)

    diagSkj = @discretise BEAST.diag(SL)[p,q] p∈Y q∈Y
    Syy = BEAST.assemble(diagSkj; materialize=hmat);

    # @hilbertspace j[1:2]
    # @hilbertspace k[1:2]
    # precond = @discretise(
    #     SL[k[1],j[1]] + SL[k[2],j[2]] == e[k[1]] + e[k[2]],
    #     j[1] ∈ Y12, j[2] ∈ Y23, k[1] ∈ Y12, k[2] ∈ Y23,);

    # BEAST.@defaultquadstrat (SL,Y12,Y12) BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(1, 2, 3, 4, 3, 3, 3, 3)
    # Syy = BEAST.sysmatrix(precond);

    # P = iN' * Syy * iN
    P = Dxy * Syy * Dyx

    u1, ch1 = solve(BEAST.GMRESSolver(Sxx,tol=2e-5,restart=250), ex)
    u2, ch2 = solve(BEAST.GMRESSolver(P*Sxx, tol=2e-5,restart=250), P*ex)

    @show ch1.iters
    @show ch2.iters

    push!(iters_classic, ch1.iters)
    push!(iters_sawtooth, ch2.iters)
# end

error()

plot(log10.(ch1.data[:resnorm]), label="MT")
plot!(log10.(ch2.data[:resnorm]), label="MT diag. Cald. precond.")

#' Check for correctness

Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]

ffd_saw_np    = potential(MWFarField3D(wavenumber=κ), pts, u1, X)
ffd_saw_pc = potential(MWFarField3D(wavenumber=κ), pts, u2, X)

plot!(norm.(ffd_saw_np), label="MT")
scatter!(norm.(ffd_saw_pc), label="MT - diag. Cald. precond.")

#' Visualise the spectrum

mSxx = BEAST.convert_to_dense(Sxx)
mSyy = BEAST.convert_to_dense(Syy)

Z = BEAST.convert_to_dense(Sxx)
W = iN' * mSyy * iN * mSxx;

wZ = eigvals(Matrix(Z))
wW = eigvals(Matrix(W))

plot(exp.(im*range(0,2pi,length=200)))
scatter!(wZ)
scatter!(wW)

# Study the various kernels
HS = Maxwell3D.singlelayer(gamma=0.0, alpha=0.0, beta=1.0)
Id = BEAST.Identity()

Z12 = BEAST.lagrangecxd0(G12)
Z23 = BEAST.lagrangecxd0(Ĝ23)
Z = Z12 × Z23

W12 = BEAST.duallagrangecxd0(G12)
W23 = BEAST.duallagrangecxd0(Ĝ23)
W = W12 × W23

DX = assemble(Id, Z, divergence(X))
HX = assemble(HS, X, X)

DY = assemble(Id, W, divergence(Y))
HY = assemble(HS, Y, Y)

Q = HY * iN * HX