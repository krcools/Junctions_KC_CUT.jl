using DrWatson
@quickactivate "Junctions_KC_CUT"

using BEAST, CompScienceMeshes, LinearAlgebra
# using JLD2

using Junctions_KC_CUT

width, height = 1.0, 0.5
overlap = 0.2

nearstrat = BEAST.DoubleNumWiltonSauterQStrat(6, 7, 6, 7, 7, 7, 7, 7)
farstrat  = BEAST.DoubleNumQStrat(1,2)

dmat(op,tfs,bfs) = BEAST.assemble(op,tfs,bfs; quadstrat=nearstrat)
hmat(op,tfs,bfs) = AdaptiveCrossApproximation.h1compress(op,tfs,bfs; nearstrat=nearstrat,farstrat=farstrat)
mat = dmat

h = [0.1, 0.05, 0.025, 0.0125]
κ = [1.0, 10.0]

h = 0.1
κ = 1.0e-3
κ = 0.00001
κ = 1.0
γ = im*1
γ = im*κ

# function runsim(;h, κ)

    SL = Maxwell3D.singlelayer(wavenumber=κ)
    WS = Maxwell3D.weaklysingular(wavenumber=κ)
    HS = Maxwell3D.hypersingular(wavenumber=κ)
    N = NCross()
    E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
    e = (n × E) × n;

    G1 = meshrectangle(width, height, h)
    G21 = meshrectangle(width, overlap, h)
    G22 = meshrectangle(width, height-overlap, h)
    CompScienceMeshes.translate!(G22, point(0, overlap, 0))
    CompScienceMeshes.rotate!(G21, 1.0π * x̂)
    CompScienceMeshes.rotate!(G22, 1.0π * x̂)
    G2 = weld(G21, G22)
    G3 = CompScienceMeshes.rotate(G1, 1.5π * x̂)

    G12 = weld(G1,-G2)
    G23 = weld(G2,-G3)
    G31 = weld(G3,-G1)

    # Ĝ23 = weld(G21, -G3)

    # G1 = meshrectangle(width, height, h)
    # G2 = CompScienceMeshes.rotate(G1, 0.5π * x̂)
    # G3 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
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

    # G12 = weld(G1,-G2)
    # G23 = weld(G2,-G3)
    # G31 = weld(G3,-G1)
    # G123 = weld(G1,G2,G3)

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

    edges12 = setminus(skeleton(G12,1), boundary(G12))
    edges23 = setminus(skeleton(Ĝ23,1), boundary(Ĝ23))

    verts12 = setminus(skeleton(G12,0), skeleton(boundary(G12),0))
    verts23 = setminus(skeleton(Ĝ23,0), skeleton(boundary(Ĝ23),0))

    Σ12 = Matrix(connectivity(G12, edges12, sign))
    Σ23 = Matrix(connectivity(Ĝ23, edges23, sign))

    Λ12 = Matrix(connectivity(G12, verts12, sign))
    Λ23 = Matrix(connectivity(Ĝ23, verts23, sign))

    PΣ12 = Σ12 * pinv(Σ12'*Σ12) * Σ12'
    PΣ23 = Σ23 * pinv(Σ23'*Σ23) * Σ23'

    PΛH12 = I - PΣ12
    PΛH23 = I - PΣ23

    ℙΛ12 = Λ12 * pinv(Λ12'*Λ12) * Λ12'
    ℙΛ23 = Λ23 * pinv(Λ23'*Λ23) * Λ23'
    
    ℙHΣ12 = I - ℙΛ12
    ℙHΣ23 = I - ℙΛ23

    # γ = ι*κ, γ <- gamma
    MR12 = γ * PΣ12 + PΛH12
    MR23 = γ * PΣ23 + PΛH23

    # Mu12 = 1/γ * PΣ12 + PΛH12
    # Mu23 = 1/γ * PΣ23 + PΛH23

    ML12 = PΣ12 + 1/γ * PΛH12
    ML23 = PΣ23 + 1/γ * PΛH23

    𝕄R12 = γ * ℙΛ12 + ℙHΣ12
    𝕄R23 = γ * ℙΛ23 + ℙHΣ23

    𝕄L12 = ℙΛ12 + 1/γ * ℙHΣ12
    𝕄L23 = ℙΛ23 + 1/γ * ℙHΣ23

    X12 = raviartthomas(G12)
    X23 = raviartthomas(Ĝ23)

    Y12 = buffachristiansen(G12)
    Y23 = buffachristiansen(Ĝ23)

    X = X12 × X23
    Y = Y12 × Y23

    @hilbertspace p[1:2]
    @hilbertspace q[1:2]

    MR = assemble(@discretise MR12[q[1],p[1]]+MR23[q[2],p[2]] q∈X p∈X)
    ML = assemble(@discretise ML12[q[1],p[1]]+ML23[q[2],p[2]] q∈X p∈X)

    MRΣ = assemble(@discretise γ*PΣ12[q[1],p[1]] + γ*PΣ23[q[2],p[2]] q∈X p∈X)
    MRΛH = assemble(@discretise PΛH12[q[1],p[1]] + PΛH23[q[2],p[2]] q∈X p∈X)

    MLΣ = assemble(@discretise PΣ12[q[1],p[1]] + PΣ23[q[2],p[2]] q∈X p∈X)
    MLΛH = assemble(@discretise 1/γ*PΛH12[q[1],p[1]] + 1/γ*PΛH23[q[2],p[2]] q∈X p∈X)

    𝕄R = assemble(@discretise 𝕄R12[q[1],p[1]]+𝕄R23[q[2],p[2]] q∈X p∈X)
    𝕄L = assemble(@discretise 𝕄L12[q[1],p[1]]+𝕄L23[q[2],p[2]] q∈X p∈X)

    𝕄RΛ = assemble(@discretise γ*ℙΛ12[q[1],p[1]] + γ*ℙΛ23[q[2],p[2]] q∈X p∈X)
    𝕄RHΣ = assemble(@discretise ℙHΣ12[q[1],p[1]] + ℙHΣ23[q[2],p[2]] q∈X p∈X)

    𝕄LΛ = assemble(@discretise ℙΛ12[q[1],p[1]] + ℙΛ23[q[2],p[2]] q∈X p∈X)
    𝕄LHΣ = assemble(@discretise 1/γ*ℙHΣ12[q[1],p[1]] + 1/γ*ℙHΣ23[q[2],p[2]] q∈X p∈X)

    Sxx = assemble(@discretise(SL[p,q], p∈X, q∈X), materialize=mat)
    WSxx = assemble(@discretise(WS[p,q], p∈X, q∈X), materialize=mat)
    # HSxx = assemble(@discretise(HS[p,q], p∈X, q∈X), materialize=mat)
    ex = BEAST.assemble(@discretise e[p] p∈X)
    
    sys0 = Sxx
    sys1 = MLΣ * Sxx * MRΣ + MLΛH * WSxx * MRΣ + MLΣ * WSxx * MRΛH + MLΛH * WSxx * MRΛH

    rhs0 = ex
    rhs1 = ML * ex

    @hilbertspace j[1:2]
    @hilbertspace k[1:2]

    Nxy = assemble(@discretise(BEAST.diag(N)[j,k], j∈X, k∈Y))
    Dyx = BEAST.GMRESSolver(Nxy, tol=2e-12, restart=250, verbose=false)
    Dxy = BEAST.GMRESSolver(transpose(Nxy), tol=2e-12, restart=250, verbose=false)

    diagSkj = @discretise BEAST.diag(SL)[p,q] p∈Y q∈Y

    Syy = assemble(@discretise(BEAST.diag(SL)[p,q], p∈Y, q∈Y), materialize=mat)
    WSyy = assemble(@discretise(BEAST.diag(WS)[p,q], p∈Y, q∈Y), materialize=mat)
    # HSxx = assemble(@discretise(HS[p,q], p∈X, q∈X), materialize=mat)
    # Syy = assemble(diagSkj; materialize=mat)
    pcd = 𝕄LHΣ * Syy * MRΣ + 𝕄LHΣ * WSyy * MRΣ + 𝕄LHΣ * WSyy * MRΛH + 𝕄LHΣ * WSyy * MRΛH

    # pcd = 𝕄L * Syy * 𝕄R

    P = Dxy * pcd * Dyx
    sys2 = P * sys1
    rhs2 = P * rhs1

    u0, ch0 = solve(BEAST.GMRESSolver(sys0, tol=2e-8, restart=250), rhs0)
    v1, ch1 = solve(BEAST.GMRESSolver(sys1, tol=2e-8, restart=250), rhs1)
    v2, ch2 = solve(BEAST.GMRESSolver(sys2, tol=2e-8, restart=250), rhs2)

    u1 = MR * v1
    u2 = MR * v2

    # error()
    # @show ch1.iters
    # @show ch2.iters

    Φ, Θ = [0.0], range(0,stop=π,length=40)
    pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]

    near0 = potential(MWFarField3D(wavenumber=κ), pts, u0, X)

    farfield = MWFarField3D(wavenumber=κ)
    near1Σ = potential(farfield, pts, MRΣ*v1, X)
    near1ΛH = potential(farfield, pts, MRΛH*v1, X)
    near1 = near1Σ .+ near1ΛH
    # near1 = potential(MWFarField3D(wavenumber=κ), pts, u1, X)

    near1Σ = potential(farfield, pts, MRΣ*u1, X)
    near1ΛH = potential(farfield, pts, MRΛH*u1, X)
    near1 = near1Σ .+ near1ΛH
    near2 = potential(MWFarField3D(wavenumber=κ), pts, u2, X)

    # u1, ch1.iters, u2, ch2.iters, near1, near2, X
# end

plot();
plot!(Θ, norm.(near0));
scatter!(Θ, norm.(near1));
scatter!(Θ, norm.(near2));
display(plot!())
# scatter!(Θ, norm.(near2))

error()

using LinearAlgebra
using Plots
plotly()
w0 = eigvals(Matrix(sys0))
w1 = eigvals(Matrix(sys1))
# w2 = eigvals(Matrix(P * sys))
plot(exp.(2pi*im*range(0,1,length=200)))
scatter!(w0)
scatter!(w1)
scatter!(w2)


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

using AlgebraicMultigrid
A = poisson(100)
b = rand(100);
solve(A, b, RugeStubenAMG(), maxiter=1, abstol=1e-6)
