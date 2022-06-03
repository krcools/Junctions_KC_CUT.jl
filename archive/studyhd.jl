import Pkg
Pkg.activate(".")

using BEAST, CompScienceMeshes, LinearAlgebra

using Plots
plotly()

include("utils.jl")

width, height, h = 1.0, 1/2, 1/2
G1 = meshrectangle(width, height, h)
G2 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
G3 = CompScienceMeshes.rotate(G1, 1.5π * x̂)
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
on_bnd23 = in(boundary(G23))

G123 = weld(G1,G2,G3)

edges12_all = skeleton(G12,1)
edges23_all = skeleton(G23,1)
on_G12 = overlap_gpredicate(edges12_all)


edges23_act_strip = submesh(edges23_all) do edge
    ch = chart(edges23_all,edge)
    v1 = simplex((ch.vertices[1],), Val{0})
    v2 = simplex((ch.vertices[2],), Val{0})
    on_bnd23(edge) && return false
    on_junction_nodes(v1) && return true
    on_junction_nodes(v2) && return true
    on_G12(ch) && return false
    return true
end

edges23_act_jagged = submesh(edges23_all) do edge
    ch = chart(edges23_all,edge)
    on_bnd23(edge) && return false
    on_junction(ch) && return true
    on_G12(ch) && return false
    return true
end


function support(screen, active_edges)
    in_active_edges = in(active_edges)
    Ĝ23 = submesh(screen) do face
        index(face[1],face[2]) |> in_active_edges && return true
        index(face[2],face[3]) |> in_active_edges && return true
        index(face[3],face[1]) |> in_active_edges && return true
        return false
    end    
end

Ĝ23_strip = support(G23, edges23_act_strip)
Ĝ23_jagged = support(G23, edges23_act_jagged)

X12 = raviartthomas(G12)
X23 = raviartthomas(G23)
X31 = raviartthomas(G31)

X23_jagged = raviartthomas(Ĝ23_jagged, edges23_act_jagged)
X23_strip = raviartthomas(Ĝ23_strip, edges23_act_strip)
# X23_all = raviartthomas(G23)

Y12 = buffachristiansen(G12)
Y23 = buffachristiansen(G23)
Y31 = buffachristiansen(G31)

Y23_jagged = buffachristiansen(Ĝ23_jagged, edges=edges23_act_jagged)
Y23_strip = buffachristiansen(Ĝ23_strip, edges=edges23_act_strip)
# Y23_all = buffachristiansen(G23)

X_jagged = X12 × X23_jagged
Y_jagged = Y12 × Y23_jagged

X_strip = X12 × X23_strip
Y_strip = Y12 × Y23_strip

X_all = X12 × X23 × X31
Y_all = Y12 × Y23 × Y31

κ = 1.0
SL = Maxwell3D.singlelayer(wavenumber=κ)
HS = Maxwell3D.singlelayer(gamma=0.0, alpha=0.0, beta=1.0)
WS = Maxwell3D.singlelayer(gamma=0.0, alpha=1.0, beta=0.0)
Id = BEAST.Identity()
N = NCross()

diag = BEAST.diag

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n;

BEAST.@defaultquadstrat (SL,X12,X12) BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(2, 3, 6, 7, 5, 5, 4, 3)
BEAST.@defaultquadstrat (SL,Y12,Y12) BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(2, 2, 3, 4, 3, 3, 3, 3)

@hilbertspace j
@hilbertspace k

X = X_all
Y = Y_all

SLxx = Matrix(assemble(@discretise(
    SL[k,j], j ∈ X, k ∈ X)));

SLyy_diag = Matrix(assemble(@discretise(
    diag(SL)[k,j], k∈Y, j∈Y)));

HSxx = Matrix(assemble(@discretise(
    HS[k,j], j ∈ X, k ∈ X)));

HSyy = Matrix(assemble(@discretise(
    HS[k,j], j ∈ Y, k ∈ Y)));

HSyy_diag = Matrix(assemble(@discretise(
    diag(HS)[k,j], k∈Y, j∈Y)));

WSyy_diag = Matrix(assemble(@discretise(
        diag(WS)[k,j], k∈Y, j∈Y)));

HSxx_diag = Matrix(assemble(@discretise(
    diag(HS)[k,j], k∈X, j∈X)));

Nxy = Matrix(assemble(@discretise(
    diag(N)[j,k], j∈X, k∈Y)));

Dyx = inv(Nxy)
Dxy = copy(transpose(Dyx))

WSxx = Matrix(assemble(WS, X, X))

DoFs = size(HSxx,1)
# primal physloops + primal single trace
@show NSx = DoFs - rank(HSxx, atol=1e-8)
# dual physloops (not relevant)
@show NSy = DoFs - rank(HSyy, atol=1e-8)

# Q: dimension of (primal geoloops) ∩ (primal singletraces)
primal_geoloops = nullspace(HSxx, atol=1e-8)
nullspace(primal_geoloops'*WSxx*primal_geoloops, atol=1e-8)


@show NSx + NSy - DoFs
@show DoFs - rank(WSxx)


# dual geoloops
@show DoFs - rank(HSyy_diag, atol=1e-8)
#primal geoloops
@show DoFs - rank(HSxx_diag, atol=1e-8)
#primal physloops
@show DoFs - rank(HSxx, atol=1e-8)
#single trace space
@show DoFs - rank(WSxx)

w1 = eigvals(SLxx)
w2 = eigvals(Dxy*SLyy_diag*Dyx*SLxx)

plot(exp.(2pi*im*range(0,1,length=200)))
scatter!(w1)
scatter!(w2)

WSHS = Dxy*WSyy_diag*Dyx*HSxx
HSWS = Dxy*HSyy_diag*Dyx*WSxx