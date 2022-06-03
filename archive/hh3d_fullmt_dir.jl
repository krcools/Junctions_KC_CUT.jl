# import Pkg
# Pkg.activate(".")

using BEAST, CompScienceMeshes, LinearAlgebra

using Plots
plotly()

include("utils.jl")

width, height, h = 1.0, 0.5, 0.025
G1 = meshrectangle(width, height, h)
G2 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
G3 = CompScienceMeshes.rotate(G1, 1.5π * x̂)

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

κ = 1.0; γ = im*κ
HS = Helmholtz3D.hypersingular(gamma=γ)
WS = Helmholtz3D.singlelayer(wavenumber=κ)
N = BEAST.Identity()

BEAST.@defaultquadstrat (WS, X, X)  BEAST.DoubleNumWiltonSauterQStrat(4, 4, -1, -1, 6, 10, 8, 3)
BEAST.@defaultquadstrat (HS, Y, Y) BEAST.DoubleNumWiltonSauterQStrat(2, 2, -1, -1, 2, 2, 2, 2)

uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=normalize(ŷ + ẑ))
e = ∂n(uⁱ)
f = BEAST.ScalarTrace(uⁱ)

@hilbertspace j
@hilbertspace k
mt = @discretise(
    WS[k,j] == f[k],
    j ∈ X, k ∈ X)

ex = BEAST.rhs(mt)
Sxx = BEAST.sysmatrix(mt)

@hilbertspace j[1:3]
@hilbertspace k[1:3]

duality = @discretise(
    N[k[1],j[1]] + N[k[2],j[2]] + N[k[3],j[3]] == e[k[1]],
    j[1] ∈ Y12, j[2] ∈ Y23, j[3] ∈ Y31, k[1] ∈ X12, k[2] ∈ X23, k[3] ∈ X31)

Nxy = BEAST.sysmatrix(duality)
iN = BEAST.GMRESSolver(Nxy; tol=1e-12, verbose=false, restart=15)

@hilbertspace j[1:3]
@hilbertspace k[1:3]

precond = @discretise(
    HS[k[1],j[1]] + HS[k[2],j[2]] + HS[k[3],j[3]] == e[k[1]] + e[k[2]] + e[k[3]],
    j[1] ∈ Y12, j[2] ∈ Y23, j[3] ∈ Y31, k[1] ∈ Y12, k[2] ∈ Y23, k[3] ∈ Y31);

Syy = BEAST.sysmatrix(precond);

P = iN' * Syy * iN

u1, ch1 = solve(BEAST.GMRESSolver(Sxx,tol=2e-5), ex)
u2, ch2 = solve(BEAST.GMRESSolver(P*Sxx, tol=2e-5, restart=200), P*ex)

@show ch1.iters
@show ch2.iters

plot()
plot!(log10.(ch1.data[:resnorm]))
plot!(log10.(ch2.data[:resnorm]))

#' Check correctness

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

Z = iN * mSxx;
W = iN' * mSyy * iN * mSxx;

wZ = eigvals(Matrix(Z))
wW = eigvals(Matrix(W))

plot(exp.(im*range(0,2pi,length=200)))
scatter!(wZ)
scatter!(wW)