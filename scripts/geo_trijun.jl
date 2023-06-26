using DrWatson
@quickactivate "Junctions_KC_CUT"

using BEAST, CompScienceMeshes, LinearAlgebra
using JLD2
# using Infiltrator

using Junctions_KC_CUT

width, height = 1.0, 0.5
overlap = 0.2

h = 0.1

G1 = meshrectangle(width, height, h)
G21 = meshrectangle(width, overlap, h)
G22 = meshrectangle(width, height-overlap, h)
CompScienceMeshes.translate!(G22, point(0, overlap, 0))
CompScienceMeshes.rotate!(G21, 0.5π * x̂)
CompScienceMeshes.rotate!(G22, 0.5π * x̂)
G2 = weld(G21, G22)
G3 = CompScienceMeshes.rotate(G1, π * x̂)

G12 = weld(G1,-G2)
G23 = weld(G2,-G3)
G31 = weld(G3,-G1)

import Plotly
G12 = CompScienceMeshes.translate(G12, h*point(0,1,1))
G23 = CompScienceMeshes.translate(G23, h*point(0,-1,1))
G31 = CompScienceMeshes.translate(G31, h*point(0,0,-1))

pt1 = patch(G12, color="red")
pt2 = patch(G23, color="blue")
pt3 = patch(G31, color="green")

wf1 = wireframe(G12, width=4)
wf2 = wireframe(G23, width=4)
wf3 = wireframe(G31, width=4)

layout = Layout(
    showlegend = false,
    scene_camera=attr(
        up     = attr(x=0, y=0, z=1),
        center = attr(x=0, y=0, z=0),
        eye    = attr(x=-1.25, y=-1.25, z=1.25)),
    scene_aspectmode="data")

plt1 = Plotly.plot([pt1, pt2, pt3, wf1, wf2, wf3], layout)
plt2 = Plotly.plot([pt1, pt2, wf1, wf2], layout)


G12 = weld(G1,-G2)
G23 = weld(G21,-G3)

G12 = CompScienceMeshes.translate(G12, h*point(0,1,1))
G23 = CompScienceMeshes.translate(G23, h*point(0,-1,1))

pt1 = patch(G12, color="red")
pt2 = patch(G23, color="blue")

wf1 = wireframe(G12, width=4)
wf2 = wireframe(G23, width=4)

plt3 = Plotly.plot([pt1, pt2, wf1, wf2], layout)

G12 = weld(G1,-G2)
G23 = weld(G2,-G3)

∂G12 = boundary(G12)
∂G23 = boundary(G23)

V12 = setminus(skeleton(G12,0), skeleton(∂G12, 0))
V23 = setminus(skeleton(G23,0), skeleton(∂G23, 0))

edges23 = skeleton(G23,1)
on_junction = overlap_gpredicate(meshsegment(1.0, 1.0, 3))
junction_edges = submesh((m,c) -> on_junction(chart(edges23,c)), edges23)
junction_nodes = skeleton(junction_edges, 0)

on_V12 = overlap_gpredicate(V12)
in_junction_nodes = in(junction_nodes)
V̂23 = submesh(V23) do m,node
    ch = chart(V23, node)
    in_junction_nodes(m,node) && return true
    on_V12(ch) && return false
    return true
end

C23 = connectivity(G23, V̂23)
using SparseArrays
G23 = submesh(G23) do m,f
    isempty(nzrange(C23,f)) && return false
    return true
end

G12 = CompScienceMeshes.translate(G12, h*point(0,1,1))
G23 = CompScienceMeshes.translate(G23, h*point(0,-1,1))

pt1 = patch(G12, color="red")
pt2 = patch(G23, color="blue")

wf1 = wireframe(G12, width=4)
wf2 = wireframe(G23, width=4)

plt4 = Plotly.plot([pt1, pt2, wf1, wf2], layout)

pltall = [plt1 plt2; plt3 plt4];
pltall.plot.layout[:showlegend] = false
pltall.plot.layout[:scene1_camera_eye] = attr(x=-1.25, y=-1.25, z=0.5)
pltall.plot.layout[:scene2_camera_eye] = attr(x=-1.25, y=-1.25, z=0.5) 
pltall.plot.layout[:scene3_camera_eye] = attr(x=-1.25, y=-1.25, z=0.5) 
pltall.plot.layout[:scene4_camera_eye] = attr(x=-1.25, y=-1.25, z=0.5) 
pltall.plot.layout[:scene1_zaxis_range] = [-0.1, 0.6]
pltall.plot.layout[:scene2_zaxis_range] = [-0.1, 0.6]
pltall.plot.layout[:scene3_zaxis_range] = [-0.1, 0.6]
pltall.plot.layout[:scene4_zaxis_range] = [-0.1, 0.6] 
display(pltall)

# plt1 = Plotly.plot([wf1, wf2, wf3, pt1, pt2, pt3], layout)
display(plt1)
# Plotly.savefig("geoall.png")



# Examples of multi-screens of type A, B, C
width, height = 1.0, 0.5
overlap = 0.2

h = 0.1

G1 = meshrectangle(width, height, h)
G21 = meshrectangle(width, overlap, h)
G22 = meshrectangle(width, height-overlap, h)
CompScienceMeshes.translate!(G22, point(0, overlap, 0))
CompScienceMeshes.rotate!(G21, 0.5π * x̂)
CompScienceMeshes.rotate!(G22, 0.5π * x̂)
G2 = weld(G21, G22)
G3 = CompScienceMeshes.rotate(G1, π * x̂)