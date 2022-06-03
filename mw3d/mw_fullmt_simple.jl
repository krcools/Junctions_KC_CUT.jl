import Pkg
Pkg.activate(".")

using BEAST, CompScienceMeshes, LinearAlgebra

using Plots
plotly()

include("utils.jl")

width, height = 1.0, 0.5

h = 0.1
G1 = meshrectangle(width, height, h)
G2 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
G3 = CompScienceMeshes.rotate(G1, 1.5π * x̂)

G12 = weld(G1,-G2)
G23 = weld(G2,-G3)
G31 = weld(G3,-G1);

G21 = weld(G2,-G1)

X12 = raviartthomas(G12)
X23 = raviartthomas(G23)
X31 = raviartthomas(G31)

X21 = raviartthomas(G21)

Y12 = buffachristiansen(G12)
Y23 = buffachristiansen(G23)
Y31 = buffachristiansen(G31);

Y21 = buffachristiansen(G21)

# X = X12 × X23 × X31
# Y = Y12 × Y23 × Y31;

X = X12 × X21
Y = Y12 × Y21

κ = 1.0
SL = Maxwell3D.singlelayer(wavenumber=κ)
N = NCross()

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n;

@hilbertspace j
@hilbertspace k
mtefie = @discretise(
    SL[k,j] == e[k],
    j ∈ X, k ∈ X);

BEAST.@defaultquadstrat (SL,X12,X12) BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(2, 3, 6, 7, 5, 5, 4, 3)
Sxx = BEAST.sysmatrix(mtefie)
ex = BEAST.rhs(mtefie)

N1 = assemble(N, X12, Y12)
N2 = assemble(N, X21, Y21)
iN = blkdiagm(inv(Matrix(N1)), inv(Matrix(N2)))

@hilbertspace j[1:2]
@hilbertspace k[1:2]
precond = @discretise(
    SL[k[1],j[1]] + SL[k[2],j[2]]== e[k[1]],
    j[1] ∈ Y12, j[2] ∈ Y23, k[1] ∈ Y12, k[2] ∈ Y23);

Syy_nondiag = assemble(SL, Y, Y)

BEAST.@defaultquadstrat (SL,Y12,Y12) BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(2, 2, 3, 4, 3, 3, 3, 3)
Syy = BEAST.sysmatrix(precond);

P = iN' * Syy * iN

u1, ch1 = solve(BEAST.GMRESSolver(Sxx,tol=2e-5), ex)
u2, ch2 = solve(BEAST.GMRESSolver(P*Sxx, tol=2e-5), P*ex)

@show ch1.iters
@show ch2.iters

plot(log10.(ch1.data[:resnorm]), label="MT")
plot!(log10.(ch2.data[:resnorm]), label="MT diag. Cald. precond.")

error()
#' Plot the far field

Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]

ffd_mt    = potential(MWFarField3D(wavenumber=κ), pts, u1, X)
ffd_mt_pc = potential(MWFarField3D(wavenumber=κ), pts, u2, X)

plot(norm.(ffd_mt), label="MT")
scatter!(norm.(ffd_mt_pc), label="MT - diag. Cald. precond.")

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