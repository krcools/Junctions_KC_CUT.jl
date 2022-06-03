using JLD2
using Plots

gr()

strip = jldopen("hh3d_strip.jld2")
fullqs = jldopen("hh3d_fullqs.jld2")
partial = jldopen("hh3d_partial.jld2")

hs = fullqs["hs"]

## Start plotting
plot(hs, fullqs["iters_classic"], marker=(:dot,6), width=2, label="full QS - no precond.")
plot!(hs, strip["iters_classic"][eachindex(hs)], marker=(:dot,6), width=2, label="strip red. - no precond.")
plot!(hs, fullqs["iters_precond"], marker=(:dot,6), width=2, label="full QS - preconditioned")
plot!(hs, partial["iters_precond"], marker=(:dot,6), width=2, label="partial red. - preconditioned")
plot!(hs, strip["iters_precond"][1:end-1], marker=(:dot,6), width=2, label="strip red. - preconditioned")

savefig("hh3d_neu_iters.png")
plot!()

strip = jldopen("hh3d_strip_dir.jld2")
fullqs = jldopen("hh3d_fullqs_dir.jld2")
partial = jldopen("hh3d_partial_dir.jld2")

hs = strip["hs"]

plot(hs, fullqs["iters_classic"], marker=(:dot,6), width=2, label="full QS - no precond.")
plot!(hs, strip["iters_classic"], marker=(:dot,6), width=2, label="strip red. - no precond.")
plot!(hs, fullqs["iters_precond"], marker=(:dot,6), width=2, label="full QS - preconditioned")
plot!(hs, partial["iters_precond"], marker=(:dot,6), width=2, label="partial red. - preconditioned")
plot!(hs, strip["iters_precond"], marker=(:dot,6), width=2, label="strip red. - preconditioned")
plot!(title="Dirichlet Problem")

savefig("hh3d_dir_iters.png")
