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
    dcc_graph(id = "basic-interactions", figure = (
        #You could use some nice list comprehensions here
        data = [
            (
                x = Yte[1,healthy],
                y = Yte[2,healthy],
                name = "(Ctrl)",
                mode = "markers",
                text = segmentDisplayNames[ healthy,:] ,  
                customdata = segmentDisplayNames[ healthy,:],
               
                marker = (size = 12,)
            ),
            (
                x = Yte[1,.!healthy],
                y = Yte[2,.!healthy],
                name = "(DKD)",
                mode = "markers",
                text = segmentDisplayNames[ .!healthy, :] ,  
                customdata = segmentDisplayNames[ .!healthy,:],
               
                marker = (size = 12,)
            )
        ],
        layout = (clickmode = "event+select",
        xaxis_title = "pca component1",
        yaxis_title = "pca component2"
        )
    )),

    html_div(
        children = [
            html_div(
                children = [
                    dcc_markdown("
                    **Hover Data**

                    Mouse over values in the graph.
                    "),
                    html_pre(id = "hover-data"),
                ],
            ),
            html_div(
                children = [
                    dcc_markdown("
                    **Click Data**

                    Click on points in the graph.
                    "),
                    html_pre(id = "click-data"),
                ],
            ),
            html_div(
                children = [
                    dcc_markdown("
                    **Selection Data**

                    Choose the lasso or rectangle tool in the graph's menu
                    bar and then select points in the graph.

                    Note that if `layout.clickmode = 'event+select'`, selection data also
                    accumulates (or un-accumulates) selected data if you hold down the shift
                    button while clicking.
                    "),
                    html_pre(id = "selected-data"),
                ],
            ),
            html_div(
                children = [
                    dcc_markdown("
                    **Zoom and Relayout Data**

                    Click and drag on the graph to zoom or click on the zoom
                    buttons in the graph's menu bar.
                    Clicking on legend items will also fire
                    this event.
                    "),
                    html_pre(id = "relayout-data"),
                ],
            ),
        ],
    )

end

callback!(
    app,
    Output("hover-data", "children"),
    Input("basic-interactions", "hoverData"),
) do hover_data
    index = isnothing(hover_data) ? 1 : hover_data.points[1].pointIndex
    return generate_table(features[index+1:index+1, :], 10)
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