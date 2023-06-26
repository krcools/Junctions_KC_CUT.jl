using DrWatson
@quickactivate "Junctions_KC_CUT"

using DataFrames
using Plots, StatsPlots
# using Latexify
gr()

function convertdf(df)
    h = df.h
    ch1 = [c.iters for c in df.ch1]
    ch2 = [c.iters for c in df.ch2]
    # df[!,:numit1] .= ch1
    # df[:,:numit2] .= ch2
    df = DataFrame(;h, ch1, ch2, u1=df[!,:u1], u2=df[!,:u2], κ=df[!,:κ])
end

sims = [
    "hh3ddir_fullmultitrace",
    "hh3ddir_partialredux",
    "hh3ddir_fixedoverlap",
    "hh3ddir_singlestrip",
    "hh3dneu_fullmultitrace",
    "hh3dneu_partialredux",
    "hh3dneu_fixedoverlap",
    "hh3dneu_singlestrip"
]

problems = ["hh3ddir", "hh3dneu"]

reductions = ["fullmultitrace", "partialredux", "fixedoverlap", "singlestrip"]
linecolors = ["red", "blue", "green", "black"]
wavenumbers = [1.0, 10.0]

figdir = datadir("figs")
for κ in wavenumbers
    for problem in problems
        title = string(problem) * "[k=$κ]"
        plt = Plots.plot(title=title)
        for (r,reduction) in enumerate(reductions)
            method = problem * "_" * reduction
            df = collect_results(datadir("simulations", method))
            df = convertdf(df)

            numdofs = [length(df[i,:u1]) for i in axes(df,1)]
            df[!,:numdofs] .= numdofs
            df[!,:scaledits1] .= (numdofs .* df[!,:ch1])
            df[!,:scaledits2] .= (numdofs .* df[!,:ch2])
            
            df = df[df.κ .== κ,[:h,:ch1,:ch2,:numdofs]]

            linestyle = [:solid :solid]
            if reduction == "singlestrip"
                linestyle = [:dot :solid]
            end 
            plt = @df df Plots.plot!(plt, :h, [:ch1,:ch2],
                title=title,
                color=linecolors[r],
                marker=[:circle :square],
                label=["$(reduction) NP" "$(reduction) CP"],
                xlabel="h",
                ylabel="Number of iterations",
                linestyle=linestyle)
        end
        display(plt)
        Plots.savefig(joinpath(figdir, title * ".png"));
    end
end

# for method in sims
#     df = collect_results(datadir("simulations", method))

#     df1 = df[df.κ .== 1.0,[:h,:ch1,:ch2]]
#     df2 = df[df.κ .== 10.0,[:h,:ch1,:ch2]]

#     df1 = convertdf(df1)
#     df2 = convertdf(df2)

#     latexify(df1, env=:table)
#     latexify(df2, env=:table)

#     figdir = datadir("figs")
#     title = string(method) * "[k=1]"
#     plt = @df df1 Plots.plot(:h, [:ch1,:ch2],
#         title=title,
#         marker=:o,
#         label=["no precond." "Calderon precond."],
#         xlabel="h",
#         ylabel="Number of iterations")
#     display(plt)
#     Plots.savefig(joinpath(figdir, title * ".png"));

#     title = string(method) * "[k=10]"
#     plt = @df df2 Plots.plot(:h, [:ch1,:ch2],
#         title=title,
#         marker=:o,
#         label=["no precond." "Calderon precond."],
#         xlabel="h",
#         ylabel="Number of iterations")
#     display(plt)
#     Plots.savefig(joinpath(figdir, title * ".png"));
# end

method = sims[1]
df1 = collect_results(datadir("simulations", method))
numdofs = [length(df1[i,:u1]) for i in axes(df,1)]