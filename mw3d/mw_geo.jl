import Pkg
Pkg.activate(".")

using BEAST, CompScienceMeshes, LinearAlgebra

import Plots
Plots.plotly()
import Plotly

include("utils.jl")

function CompScienceMeshes.patch(Γ::CompScienceMeshes.AbstractMesh;
    intensity=nothing,
    caxis=nothing,
    showscale=true,
    color="red",
    opacity=1.0)

    v = vertexarray(Γ)
    c = cellarray(Γ)

    x = v[:,1];    y = v[:,2];    z = v[:,3]
    i = c[:,1].-1; j = c[:,2].-1; k = c[:,3].-1

    if intensity == nothing
        return Plotly.mesh3d(;
            x=x, y=y, z=z,
            i=i, j=j, k=k,
            color=color,
            opacity=opacity,
        )
    end

    m, M = extrema(intensity)
    if caxis != nothing
        m, M = caxis
    end

    s = Plotly.mesh3d(;
        x=x, y=y, z=z,
        i=i, j=j, k=k,
        
        intensitymode="cell",
        intensity=intensity,
        colorscale="Viridis",
        showscale=showscale,
        cmin=m,
        cmax=M,
        opacity=opacity,
    )

end

width, height, h = 1.0, 1/2, 1/8
h = 0.0125
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

δ = 0.05
G12_disp = CompScienceMeshes.translate(G12, δ*point(0,1,1))
G23_disp = CompScienceMeshes.translate(G23, δ*point(0,-1,1))
G31_disp = CompScienceMeshes.translate(G31, δ*point(0,0,-1))

Ĝ23_jagged_disp = CompScienceMeshes.translate(Ĝ23_jagged, δ*point(0,-1,1))
Ĝ23_strip_disp = CompScienceMeshes.translate(Ĝ23_strip, δ*point(0,-1,1))

G123_disp = CompScienceMeshes.weld(G12_disp, Ĝ23_jagged_disp)
I = zeros(length(G123_disp))
I[1:length(G12_disp)] .= 1
Plotly.plot(patch(G123_disp, intensity=I))

pt1 = patch(G12_disp, color="red", opacity=0.75)
pt2 = patch(Ĝ23_jagged_disp, color="blue", opacity=0.75)
pt3 = patch(Ĝ23_strip_disp, color="blue", opacity=0.75)
pt4 = patch(G23_disp, color="blue", opacity=0.75)
pt5 = patch(G31_disp, color="green", opacity=0.75)

wf1 = CompScienceMeshes.wireframe(skeleton(G12_disp,1), width=2)
wf2 = CompScienceMeshes.wireframe(skeleton(Ĝ23_jagged_disp,1), width=2)
wf3 = CompScienceMeshes.wireframe(skeleton(Ĝ23_strip_disp,1), width=2)
wf4 = CompScienceMeshes.wireframe(skeleton(G23_disp,1), width=2)
wf5 = CompScienceMeshes.wireframe(skeleton(G31_disp,1), width=2)

Plotly.plot([pt1,pt3,wf1,wf3], Plotly.Layout(scene_aspectmode="data"))
Plotly.plot([pt1,pt4,pt5,wf1,wf4,wf5], Plotly.Layout(scene_aspectmode="data"))
