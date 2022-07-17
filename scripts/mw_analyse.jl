using DrWatson
@quickactivate "Junctions_KC_CUT"

using DataFrames
# using Plots, StatsPlots
using Plotly

methods = [
    "mw_engineering",
    "mw_fullmultitrace",
    "mw_partialredux",
    "mw_fixedoverlap",
    "mw_singlestrip",
]

method_names = [
    "Engineering",
    "Full MT",
    "Partial MT",
    "Fixed Overlap",
    "Single Strip"    
]

line_colors = [
    "black",
    "blue",
    "red",
    "green",
    "magenta"
]

linsolvers = ["Classic", "Precond"]
iterfields = [:ch1, :ch2]
marker_symbols = ["circle", "square"]
mode = "lines+markers"
marker_size = 12

freq = 10.0
traces = []

# method = methods[1]
# method_name = method_names[1]
# line_color = line_colors[1]
for (method, method_name, line_color) in zip(methods, method_names, line_colors)

    @show method
    df = collect_results(datadir("simulations", method))
    df1 = df[df.κ .== freq, [:h, :ch1, :ch2]]
    x = df1[:,:h]

    # linsolver = linsolvers[1]
    # iterfield = iterfields[1]
    for (linsolver, iterfield, marker_symbol) in zip(linsolvers, iterfields, marker_symbols)

        y = df1[:,iterfield]
        name = method_name * " - " * linsolver
        tr = Plotly.scatter(;x, y, mode, name, marker_size, marker_symbol, line_color)
        push!(traces, tr)
    end
end

traces = [t for t in traces]
layout = Layout(;
    title="wavenumber k = $(freq)",
    xaxis_title="mesh size h",
    yaxis_title="#GMRES iterations")
plt = Plotly.plot(traces, layout)
Plotly.post(plt)
