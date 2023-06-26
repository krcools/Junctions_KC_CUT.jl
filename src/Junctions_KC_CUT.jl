module Junctions_KC_CUT

export setminus

using BEAST
using CompScienceMeshes
# using AdaptiveCrossApproximation

# AdaptiveCrossApproximation.blockassembler(op,Y,X;quadstrat) = BEAST.blockassembler(op,Y,X;quadstrat)
# AdaptiveCrossApproximation.scalartype(op,Y,X) = BEAST.scalartype(op,Y,X)
# AdaptiveCrossApproximation.positions(X) = BEAST.positions(X)
# AdaptiveCrossApproximation.numfunctions(X) = BEAST.numfunctions(X)


# Top = BEAST.MWSingleLayer3D
# Tsp = BEAST.RTRefSpace
# Trf = BEAST.RTRefSpace

# function BEAST.quaddata(op::Top, tref::Trf, bref::Trf,
#     tels, bels, qs::BEAST.DoubleNumQStrat)

#     qs = BEAST.DoubleNumWiltonSauterQStrat(qs.outer_rule, qs.inner_rule, 1, 1, 1, 1, 1, 1)
#     BEAST.quaddata(op, tref, bref, tels, bels, qs)
# end

# function BEAST.quadrule(op::Top, tref::Trf, bref::Trf,
#     i ,τ, j, σ, qd, qs::BEAST.DoubleNumQStrat)

#     return BEAST.DoubleQuadRule(
#         qd.tpoints[1,i],
#         qd.bpoints[1,j])
# end


include("utils.jl")

end # module
