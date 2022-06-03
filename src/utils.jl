function blkdiagm(blks...)
    m = sum(size(b,1) for b in blks)
    n = sum(size(b,2) for b in blks)
    A = zeros(eltype(first(blks)), m, n)
    o1 = 1
    o2 = 1
    for b in blks
        m1 = size(b,1)
        n1 = size(b,2)
        A[o1:o1+m1-1,o2:o2+n1-1] .= b
        o1 += m1
        o2 += n1
    end
    A
end

setminus(A::CompScienceMeshes.AbstractMesh, B::CompScienceMeshes.AbstractMesh) = submesh(!in(B),A)

function scatterplot(m::CompScienceMeshes.AbstractMesh)

    @assert dimension(m) == 0
    @assert universedimension(m) == 3

    V = vertices(m)
    V1 = [V[c[1]] for c in m]

    Plotly.scatter3d(x=getindex.(V1,1), y=getindex.(V1,2), z=getindex.(V1,3))
end

import Plotly
function showfn(space,i)
    geo = geometry(space)
    T = coordtype(geo)
    X = Dict{Int,T}()
    Y = Dict{Int,T}()
    Z = Dict{Int,T}()
    U = Dict{Int,T}()
    V = Dict{Int,T}()
    W = Dict{Int,T}()
    for sh in space.fns[i]
        chrt = chart(geo, cells(geo)[sh.cellid])
        nbd = CompScienceMeshes.center(chrt)
        vals = refspace(space)(nbd)
        x,y,z = cartesian(nbd)
        # @show vals[sh.refid].value
        u,v,w = vals[sh.refid].value
        # @show x, y, z
        # @show u, v, w
        X[sh.cellid] = x
        Y[sh.cellid] = y
        Z[sh.cellid] = z
        U[sh.cellid] = get(U,sh.cellid,zero(T)) + sh.coeff * u
        V[sh.cellid] = get(V,sh.cellid,zero(T)) + sh.coeff * v
        W[sh.cellid] = get(W,sh.cellid,zero(T)) + sh.coeff * w
    end
    X = collect(values(X))
    Y = collect(values(Y))
    Z = collect(values(Z))
    U = collect(values(U))
    V = collect(values(V))
    W = collect(values(W))
    Plotly.cone(x=X,y=Y,z=Z,u=U,v=V,w=W)
end

function compress!(space)
    T = scalartype(space)
    for (i,fn) in pairs(space.fns)
        shapes = Dict{Tuple{Int,Int},T}()
        for shape in fn
            v = get(shapes, (shape.cellid, shape.refid), zero(T))
            shapes[(shape.cellid, shape.refid)] = v + shape.coeff
            # set!(shapes, (shape.cellid, shape.refid), v + shape.coeff)
        end
        space.fns[i] = [BEAST.Shape(k[1], k[2], v) for (k,v) in shapes]
    end
end

# const CSM = CompScienceMeshes
# function CSM.patch(Γ::CSM.AbstractMesh, fcr=nothing; caxis=nothing, showscale=true)

#     v = vertexarray(Γ)
#     c = cellarray(Γ)

#     x = v[:,1];    y = v[:,2];    z = v[:,3]
#     i = c[:,1].-1; j = c[:,2].-1; k = c[:,3].-1

#     if fcr == nothing
#         a = [cartesian(CompScienceMeshes.center(chart(Γ,cell)))[3] for cell in cells(Γ)]
#     else
#         a = fcr
#     end

#     m, M = extrema(a)
#     if caxis != nothing
#         m, M = caxis
#     end

#     s = Plotly.mesh3d(;
#         x=x, y=y, z=z,
#         i=i, j=j, k=k,
        
#         intensitymode="cell",
#         intensity=a,
#         colorscale="Viridis",
#         showscale=showscale,
#         cmin=m,
#         cmax=M
#     )

# end