using DataFrames, Dash, DashHtmlComponents, DashCoreComponents, UrlDownload, PlotlyJS, JSON3
using Statistics
using LinearAlgebra
using CSV

app = dash()

Yte=Matrix(CSV.read("PCA_matrix.txt",DataFrame))
#PCM = DataFrame(CSV.File("data/Kidney_Q3Norm_TargetCountMatrix.txt"))
features = DataFrame(CSV.File("data/Kidney_Sample_Annotations.txt"))


function generate_table(dataframe, max_rows = 10)
    html_table([
        html_thead(html_tr([html_th(col) for col in names(dataframe)])),
        html_tbody([
            html_tr([html_td(dataframe[r, c]) for c in names(dataframe)]) for r = 1:min(nrow(dataframe), max_rows)
        ]),
    ])
end


app = dash()

available_indicators = ["disease_status","region"#,"pathology"]
]

patients = unique(features.SlideName)

app.layout = html_div() do
    html_div(
        children = [
            dcc_dropdown(
                id = "crossfilter-roi-type",
                options = [
                    (label = i, value = i)
                    for i in available_indicators
                ],
                value = "region",
            ),
            dcc_graph(
                id = "graph-1",
            ),

            dcc_dropdown(
                id = "crossfilter-patient",
                options = [
                    (label = i, value = i)
                    for i in patients
                ],
                multi = true,
                value =patients[1],
            ),
            dcc_graph(
                id = "graph-2",
            ),

            dcc_dropdown(
                id = "group1",
                options = [
                    (label = i, value = i)
                    for i in unique(features.SegmentDisplayName)
                ],
                multi = true,
                value =features[1,"SegmentDisplayName"],
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
        ]
    )
end

callback!(
    app,
    Output("graph-1", "figure"),
    Input("crossfilter-roi-type", "value"),
    # Input("crossfilter-yaxis-column", "value"),
    # Input("crossfilter-xaxis-type", "value"),
    # Input("crossfilter-yaxis-type", "value"),
    # Input("crossfilter-year-slider", "value"),
) do roi_type

    # df6f = df6[df6.year .== year_slider_value, :]

    selector = roi_type
    statusList = unique(features[!,selector])
    plotData = [    ( x = Yte[1, features[!,selector] .== status], y = Yte[2, features[!,selector] .== status],  type = "scatter", name = status, mode = "markers",) for status in statusList]

    return (
        data = plotData,
        layout = (
            title = "By Structure",
            xaxis_title = "pca_1",
            yaxis_title = "pca_2",
        ),
    )

end

callback!(
    app,
    Output("graph-2", "figure"),
    Input("crossfilter-patient", "value"),
    # Input("crossfilter-yaxis-column", "value"),
    # Input("crossfilter-xaxis-type", "value"),
    # Input("crossfilter-yaxis-type", "value"),
    # Input("crossfilter-year-slider", "value"),
) do patientList

    # df6f = df6[df6.year .== year_slider_value, :]

    selector = "SlideName"
    #this list should be changed by the checkboxes
    plotData = [    ( x = Yte[1, features[!,selector] .== status], y = Yte[2, features[!,selector] .== status],  type = "scatter", name = status, mode = "markers",) for status in patientList]

    return (
        data = plotData,
        layout = (
            title = "By Patient",
            xaxis_title = "pca_1",
            yaxis_title = "pca_2",
        ),
    )
end

callback!(
    app,
    Output("group1", "value"),
    Input("graph-1", "selectedData"),
) do selected_data
    selectedpoints = 1:2
    if selected_data != nothing
        selectedpoints = [p[:pointIndex] + 1 for p in selected_data.points]
    end
    #print(features[selectedpoints, "SegmentDisplayName"])
    return features[selectedpoints, "SegmentDisplayName"]
end

run_server(app, "0.0.0.0", debug=true)