using DrWatson
@quickactivate "Junctions_KC_CUT"

using DataFrames
using Plots, StatsPlots

function convertdf(df)
    h = df.h
    ch1 = [c.iters for c in df.ch1]
    ch2 = [c.iters for c in df.ch2]
    df = DataFrame(;h, ch1, ch2)
end

methods = [
    "hh3ddir_fullmultitrace",
    "hh3ddir_partialredux",
    "hh3ddir_fixedoverlap",
    "hh3ddir_singlestrip",
    "hh3dneu_fullmultitrace",
    "hh3dneu_partialredux",
    "hh3dneu_fixedoverlap",
    "hh3dneu_singlestrip"
]

for method in methods
    df = collect_results(datadir("simulations", method))

    df1 = df[df.κ .== 1.0,[:h,:ch1,:ch2]]
    df2 = df[df.κ .== 10.0,[:h,:ch1,:ch2]]

    df1 = convertdf(df1)
    df2 = convertdf(df2)

    latexify(df1, env=:table)
    latexify(df2, env=:table)

    figdir = datadir("figs")
    title = string(method) * "[k=1]"
    @df df1 plot(:h, [:ch1,:ch2], title=title, label=["no precond." "Calderon precond."])
    savefig(joinpath(figdir, title * ".png"));

    title = string(method) * "[k=10]"
    @df df2 plot(:h, [:ch1,:ch2], title=title, label=["no precond." "Calderon precond."])
    savefig(joinpath(figdir, title * ".png"));
end