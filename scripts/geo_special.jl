
using DrWatson
@quickactivate "Junctions_KC_CUT"

using BEAST, CompScienceMeshes, LinearAlgebra
using JLD2

using Junctions_KC_CUT

h = 0.1

width = 1.0
I1 = meshrectangle(width, width/2, h)
I2 = CompScienceMeshes.translate(I1, -width*x̂)
I3 = CompScienceMeshes.translate(I1, -width*x̂ - 0.5*width*ŷ)
I4 = CompScienceMeshes.translate(I1, -0.5*width*ŷ)
I5 = CompScienceMeshes.rotate(I1, -pi/2*ŷ)

G1 = weld(I5, I2, I3, I4)
G2 = weld(-I5, I1, I4, I3)
G3 = weld(-I1, -I2, -I3, -I4)

using Plotly

G1 = CompScienceMeshes.translate(G1, h*point(-1,0,1))
G2 = CompScienceMeshes.translate(G2, h*point(1,0,1))
G3 = CompScienceMeshes.translate(G3, h*point(0,0,-1))

pt1 = CompScienceMeshes.patch(G1, color="red", opacity=1.0)
wf1 = CompScienceMeshes.wireframe(G1, width=4)

pt2= CompScienceMeshes.patch(G2, color="blue", opacity=1.0)
wf2 = CompScienceMeshes.wireframe(G2, width=4)

pt3 = CompScienceMeshes.patch(G3, color="green", opacity=1.0)
wf3 = CompScienceMeshes.wireframe(G3, width=4)
Plotly.plot([pt1, pt2, pt3, wf1, wf2, wf3],
    Plotly.Layout(;
        showlegend = false
    ))


G1 = weld(I5, I2, I3)
G2 = weld(-I5, I1, I4)
G3 = weld(-I1, -I2, -I3, -I4)

# Plots.heatmap(zs, xs, real.(getindex.(near1,1)), clims=(-15,15))
# Plots.heatmap(zs, xs, real.(getindex.(near1,3)), clims=(-15,15))

# plt_pts = Plotly.scatter3d(
#     x=vec(getindex.(pts,1)),
#     y=vec(getindex.(pts,2)),
#     z=vec(getindex.(pts,3)))
# Plotly.Plot([plt1, plt_pts])

# @show xs[40]
# Plots.plot(real.(getindex.(near1[40,:],3)), ylims=(-15,15))
# Plots.scatter!(real.(getindex.(near2[40,:],3)), ylims=(-15,15))

# @show ys[60]
# Plots.plot(real.(getindex.(near1[60,:],3)), ylims=(-15,15))
# Plots.scatter!(real.(getindex.(near2[60,:],3)), ylims=(-15,15))

# @show zs[50]
# Plots.plot(real.(getindex.(near1[:,50],1)), ylims=(-15,15))
# Plots.scatter!(real.(getindex.(near2[:,50],1)), ylims=(-15,15))

# Create meshless geos for type A, B, C multi-screens
h = 0.1

