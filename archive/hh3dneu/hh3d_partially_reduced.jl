# import Pkg
# Pkg.activate(".")

using BEAST, CompScienceMeshes, LinearAlgebra

include("utils.jl")

width, height, h = 1.0, 0.5, 0.05
G1 = meshrectangle(width, height, h)
G2 = CompScienceMeshes.rotate(G1, 1.0π * x̂)
G3 = CompScienceMeshes.rotate(G1, 1.5π * x̂)

G12 = weld(G1,-G2)
G23 = weld(G2,-G3)
G31 = weld(G3,-G1);

# import Plotly
# pt1 = patch(G12)
# pt2 = patch(G23)
# pt3 = patch(G31)
# Plotly.plot([pt1,pt2,pt3])

∂G12 = boundary(G12)
∂G23 = boundary(G23)
∂G31 = boundary(G31)

V12 = setminus(skeleton(G12,0), skeleton(∂G12, 0))
V23 = setminus(skeleton(G23,0), skeleton(∂G23, 0))
V31 = setminus(skeleton(G31,0), skeleton(∂G31, 0))

X12 = lagrangec0d1(G12, V12)
X23 = lagrangec0d1(G23, V23)
X31 = lagrangec0d1(G31, V31)

Y12 = BEAST.duallagrangecxd0(G12, V12)
Y23 = BEAST.duallagrangecxd0(G23, V23)
Y31 = BEAST.duallagrangecxd0(G31, V31)

X = X12 × X23
Y = Y12 × Y23

κ = 1.0; γ = im*κ
HS = Helmholtz3D.hypersingular(gamma=γ)

uⁱ = Helmholtz3D.planewave(wavenumber=κ, direction=normalize(ŷ + ẑ))
e = ∂n(uⁱ)

@hilbertspace j
@hilbertspace k
mt = @discretise(
    HS[k,j] == e[k],
    j ∈ X, k ∈ X)

# BEAST.@defaultquadstrat (SL,X12,X12) BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(2, 3, 6, 7, 5, 5, 4, 3)
Sxx = BEAST.sysmatrix(mt)
ex = BEAST.rhs(mt)

@hilbertspace j[1:2]
@hilbertspace k[1:2]
N = BEAST.Identity()
duality = @discretise(
    N[k[1],j[1]] + N[k[2],j[2]] == e[k[1]],
    j[1] ∈ Y12, j[2] ∈ Y23,k[1] ∈ X12, k[2] ∈ X23)

Nxy = BEAST.sysmatrix(duality)
iN = BEAST.GMRESSolver(Nxy; tol=1e-5, verbose=false)

@hilbertspace j[1:2]
@hilbertspace k[1:2]
WS = Helmholtz3D.singlelayer(wavenumber=κ)
precond = @discretise(
    WS[k[1],j[1]] + WS[k[2],j[2]] == e[k[1]],
    j[1] ∈ Y12, j[2] ∈ Y23, k[1] ∈ Y12, k[2] ∈ Y23)

# BEAST.@defaultquadstrat (SL,Y12,Y12) BEAST.DoubleNumWiltonSauterQStrat{Int64, Int64}(2, 2, 3, 4, 3, 3, 3, 3)
Syy = BEAST.sysmatrix(precond);

P = iN' * Syy * iN

u1, ch1 = solve(BEAST.GMRESSolver(Sxx,tol=2e-5), ex)
u2, ch2 = solve(BEAST.GMRESSolver(P*Sxx, tol=2e-5), P*ex)

@show ch1.iters
@show ch2.iters

phi = pi/2
thetas = range(0,pi,length=200)
farpts = [10*point(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)) for theta in thetas]

nfp = BEAST.HH3DDoubleLayerNear(wavenumber=κ)
near1 = BEAST.potential(nfp, farpts, u1, X, type=ComplexF64)
near2 = BEAST.potential(nfp, farpts, u2, X, type=ComplexF64)
inc = uⁱ.(farpts)
# tot = near - inc

#' Check if the far field matches

using Plots
plotly()
# contour(zs,ys,imag.(near-inc),levels=50)
plot(thetas, abs.(vec(near1)), label="classic")
scatter!(thetas, abs.(vec(near2)), label="partially reduced")

#' Compute and plot the eigenvalues

mSxx = BEAST.convert_to_dense(Sxx)
mSyy = BEAST.convert_to_dense(Syy)

Z = BEAST.convert_to_dense(Sxx)
W = iN' * mSyy * iN * mSxx;

wZ = eigvals(Matrix(Z))
wW = eigvals(Matrix(W))

plotly()
plot(exp.(im*range(0,2pi,length=200)))
scatter!(wZ)
scatter!(wW)