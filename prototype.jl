using DataFrames, Dash, DashHtmlComponents, DashCoreComponents, UrlDownload, PlotlyJS, JSON3
# Packages we will use throughout this notebook
using UMAP
using Makie
using XLSX
using VegaDatasets
using DataFrames
using MultivariateStats
using RDatasets
using StatsBase
using Statistics
using LinearAlgebra
#using Plots
using ScikitLearn
using MLBase
using Distances
using CSV

df6 = DataFrame(urldownload("https://raw.githubusercontent.com/plotly/datasets/master/country_indicators.csv"))

rename!(df6, Dict(:"Year" => "year"))

dropmissing!(df6)

available_indicators = unique(df6[:, "Indicator Name"])
years = unique(df6[:, "year"])

function generate_table(dataframe, max_rows = 10)
    html_table([
        html_thead(html_tr([html_th(col) for col in names(dataframe)])),
        html_tbody([
            html_tr([html_td(dataframe[r, c]) for c in names(dataframe)]) for r = 1:min(nrow(dataframe), max_rows)
        ]),
    ])
end

features = DataFrame(CSV.File("data/Kidney_Sample_Annotations.txt"))
PCM = DataFrame(CSV.File("data/Kidney_Q3Norm_TargetCountMatrix.txt"))
M = Matrix(PCM[:,2:end])
data = M
data = (data .- mean(data,dims = 2))./ std(data,dims=2)
p = fit(PCA,data,maxoutdim=2)
Yte = MultivariateStats.transform(p, data) #notice that Yte[:,1] is the same as P'*(data[1,:]-mean(p))

#Getting Labels
segmentDisplayNames=features[!,r"SegmentDi."]
healthy=[!occursin(r"disease.",seg) for seg in segmentDisplayNames[!,1]];
glomeruli = [occursin(r".Geo.",seg) for seg in segmentDisplayNames[!,1]];
distTub = [occursin(r".Pan.",seg) for seg in segmentDisplayNames[!,1]];
proxTub = [occursin(r".neg",seg) for seg in segmentDisplayNames[!,1]];


app = dash()

app.layout = html_div() do


     html_div(
        children = [
            html_div(
                children = [
                    dcc_dropdown(
                        id = "crossfilter-xaxis-column",
                        options = [
                            (label = i, value = i)
                            for i in available_indicators
                        ],
                        value = "Fertility rate, total (births per woman)",
                    ),
                    dcc_radioitems(
                        id = "crossfilter-xaxis-type",
                        options = [
                            (label = i, value = i) for i in ["linear", "log"]
                        ],
                        value = "linear",
                    ),
                ],
                style = (width = "49%", display = "inline-block"),
            ),
            html_div(
                children = [
                    dcc_dropdown(
                        id = "crossfilter-yaxis-column",
                        options = [
                            (label = i, value = i)
                            for i in available_indicators
                        ],
                        value = "Life expectancy at birth, total (years)",
                    ),
                    dcc_radioitems(
                        id = "crossfilter-yaxis-type",
                        options = [
                            (label = i, value = i) for i in ["linear", "log"]
                        ],
                        value = "linear",
                    ),
                ],
                style = (
                    width = "49%",
                    float = "right",
                    display = "inline-block",
                ),
            ),
        ],
        style = (
            borderBottom = "thin lightgrey solid",
            backgroundColor = "rgb(250, 250, 250)",
            padding = "10px 5px",
        ),
    ),
    html_div(
        children = [
            dcc_graph(id = "crossfilter-indicator-scatter"),
            dcc_slider(
                id = "crossfilter-year-slider",
                min = minimum(years),
                max = maximum(years),
                marks = Dict([Symbol(v) => Symbol(v) for v in years]),
                value = minimum(years),
                step = nothing,
            ),
        ],
        style = (
            width = "49%",
            display = "inline-block"
        ),
    ),
    html_div(
        children = [
            html_div(
                children = [
                    dcc_markdown("
                    # DE
                    "),
                    html_pre(id = "DE-plot"),
                ]
            ),
            html_div(
                children = [
                    dcc_markdown("
                    # Spatial stuff
                    "),
                    html_pre(id = "spatialCellDecon"),
                ]
            ),
        ],
        style = (width = "49%", display = "inline-block"),
    )
end

callback!(
    app,
    Output("crossfilter-indicator-scatter", "figure"),
    Input("crossfilter-xaxis-column", "value"),
    Input("crossfilter-yaxis-column", "value"),
    Input("crossfilter-xaxis-type", "value"),
    Input("crossfilter-yaxis-type", "value"),
    Input("crossfilter-year-slider", "value"),
) do xaxis_column_name, yaxis_column_name, xaxis_type, yaxis_type, year_slider_value

    df6f = df6[df6.year .== year_slider_value, :]

    return Plot(
        df6f[df6f[!, Symbol("Indicator Name")] .== xaxis_column_name, :Value],
        df6f[df6f[!, Symbol("Indicator Name")] .== yaxis_column_name, :Value],
        Layout(
            xaxis_type = xaxis_type == "Linear" ? "linear" : "log",
            xaxis_title = xaxis_column_name,
            yaxis_title = yaxis_column_name,
            yaxis_type = yaxis_type == "Linear" ? "linear" : "log",
            hovermode = "closest",
            height = 450,
        ),
        kind = "scatter",
        text = df6f[
            df6f[!, Symbol("Indicator Name")] .== yaxis_column_name,
            Symbol("Country Name"),
        ],
        customdata = df6f[
            df6f[!, Symbol("Indicator Name")] .== yaxis_column_name,
            Symbol("Country Name"),
        ],
        mode = "markers",
        marker_size = 15,
        marker_opacity = 0.5,
        marker_line_width = 0.5,
        marker_line_color = "white",
    )
end

callback!(
    app,
    Output("hover-data", "children"),
    Input("basic-interactions", "hoverData"),
) do hover_data
    index = isnothing(hover_data) ? 1 : hover_data.points[1].pointIndex
    return generate_table(features[index+1:index+1, :], 28)
end

callback!(
    app,
    Output("selected-data", "children"),
    Input("basic-interactions", "selectedData"),
) do selected_data
    selectedpoints = 1:2
    if selected_data != nothing
        selectedpoints = [p[:pointIndex] + 1 for p in selected_data.points]
    end
    return generate_table(features[selectedpoints, :], 10)
end

run_server(app, "0.0.0.0", debug=true)