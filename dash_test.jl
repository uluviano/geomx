using DataFrames, Dash, DashHtmlComponents, DashCoreComponents, UrlDownload, PlotlyJS, JSON3
using MultivariateStats
# using RDatasets
using StatsBase
using Statistics
# using LinearAlgebra
#using Plots
# using ScikitLearn
# using MLBase
using Distances
using CSV

# df6 = DataFrame(urldownload("https://raw.githubusercontent.com/plotly/datasets/master/country_indicators.csv"))

# rename!(df6, Dict(:"Year" => "year"))

# dropmissing!(df6)

# available_indicators = unique(df6[:, "Indicator Name"])
# years = unique(df6[:, "year"])

# print(available_indicators)

features = DataFrame(CSV.File("data/Kidney_Sample_Annotations.txt"))
PCM = DataFrame(CSV.File("data/Kidney_Q3Norm_TargetCountMatrix.txt"))
M = Matrix(PCM[:,2:end])
data = M
data = (data .- mean(data,dims = 2))./ std(data,dims=2)
p = fit(PCA,data,maxoutdim=2)
Yte = MultivariateStats.transform(p, data) #notice that Yte[:,1] is the same as P'*(data[1,:]-mean(p))
print(size(Yte))


for name in names(features)
    print(":", name, ", ")
end


print('\n')
segmentDisplayNames=features[!,r"SegmentDi."]


#Getting Labels
segmentDisplayNames=features[!,r"SegmentDi."]
healthy = [!occursin(r"disease.",seg) for seg in segmentDisplayNames[!,1]];
glomeruli = [occursin(r".Geo.",seg) for seg in segmentDisplayNames[!,1]];
distTub = [occursin(r".Pan.",seg) for seg in segmentDisplayNames[!,1]];
proxTub = [occursin(r".neg",seg) for seg in segmentDisplayNames[!,1]];


available_indicators = ["sick", "glomeruli"]

# print(segmentDisplayNames)
print(size(glomeruli))
print(size(distTub))

df5 = DataFrame(
    x = [1, 2, 1, 2],
    y = [1, 2, 3, 4],
    customdata = [1, 2, 3, 4],
    fruit = ["apple", "apple", "orange", "orange"],
)


app = dash()

app.layout = html_div() do
    html_div(
        children = [
            html_div(
                children = [
                    dcc_dropdown(
                        id = "crossfilter-roi-type",
                        options = [
                            (label = i, value = i)
                            for i in available_indicators
                        ],
                        value = "sick",
                    ),
                ],
                style = (width = "49%", display = "inline-block"),
            ),

            dcc_graph(id = "crossfilter-indicator-scatter"),
        ],
        style = (
        borderBottom = "thin lightgrey solid",
        backgroundColor = "rgb(250, 250, 250)",
        padding = "10px 5px",
        ),
    )
    # html_div(
    #     children = [
    #         dcc_graph(id = "crossfilter-indicator-scatter"),
    #     ],
    #     style = (
    #         width = "49%",
    #         display = "inline-block"
    #     ),
    # )
end



# Callbacks!

callback!(
    app,
    Output("crossfilter-indicator-scatter", "figure"),
    Input("crossfilter-roi-type", "value"),
    # Input("crossfilter-yaxis-column", "value"),
    # Input("crossfilter-xaxis-type", "value"),
    # Input("crossfilter-yaxis-type", "value"),
    # Input("crossfilter-year-slider", "value"),
) do roi_type

    # df6f = df6[df6.year .== year_slider_value, :]

    if roi_type == "sick"
        regex = "disease."
        selection = [occursin(Regex(regex),seg) for seg in segmentDisplayNames[!,1]];
    end

    if roi_type == "glomeruli"
        regex = "disease."
        selection = [occursin(r".Geo.",seg) for seg in segmentDisplayNames[!,1]];
    end

    print(size(selection))
    print(size(Yte))

    return Plot(
        Yte[1,selection],
        Yte[2,selection],
        # Layout(
        #     xaxis_type = xaxis_type == "Linear" ? "linear" : "log",
        #     xaxis_title = xaxis_column_name,
        #     yaxis_title = yaxis_column_name,
        #     yaxis_type = yaxis_type == "Linear" ? "linear" : "log",
        #     hovermode = "closest",
        #     height = 450,
        # ),
        kind = "scatter",
        mode = "markers",
        marker_size = 15,
        marker_opacity = 0.5,
        marker_line_width = 0.5,
        marker_line_color = "white",
    )
end


run_server(app, "0.0.0.0", debug=true)