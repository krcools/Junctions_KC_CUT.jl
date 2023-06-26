using DrWatson
@quickactivate "Junctions_KC_CUT"

using DataFrames
using Plots, StatsPlots
# using Latexify

problems = ["hh3ddir"]

reductions = ["special"]
linecolors = ["red"]
wavenumber = [10.0]

figdir = datadir("figs")
for κ in wavenumbers
    for problem in problems
        title = string(problem) * "[k=$κ]"
        plt = Plots.plot(title=title)
        for (r,reduction) in enumerate(reductions)
            method = problem * "_" * reduction
            df = collect_results(datadir("simulations", method))
            # df = convertdf(df)

            numdofs = [length(df[i,:u1]) for i in axes(df,1)]
            df[!,:numdofs] .= numdofs
            df[!,:scaledits1] .= (numdofs .* df[!,:ch1])
            df[!,:scaledits2] .= (numdofs .* df[!,:ch2])
            
            df = df[df.κ .== κ,[:h,:ch1,:ch2,:numdofs]]

            plt = @df df Plots.plot!(plt, :h, [:ch1,:ch2],
                title=title,
                color=linecolors[r],
                marker=[:circle :square],
                label=["$(reduction) NP" "$(reduction) CP"],
                xlabel="h",
                ylabel="Number of iterations",
                ylims = [0,70])
        end
        display(plt)
        Plots.savefig(joinpath(figdir, title * ".png"));
    end
end