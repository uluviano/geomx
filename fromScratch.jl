using DataFrames, Dash, DashHtmlComponents, DashCoreComponents, UrlDownload, PlotlyJS, JSON3
using Statistics
using LinearAlgebra
using CSV

app = dash()

Yte=Matrix(CSV.read("PCA_matrix.txt",DataFrame))
#PCM = DataFrame(CSV.File("data/Kidney_Q3Norm_TargetCountMatrix.txt"))
features = DataFrame(CSV.File("data/Kidney_Sample_Annotations.txt"))
structuresDict = Dict("abnormal"=>"Glom (Abnormal)","healthy"=>"Glom (Healthy)"," PanCK" => "Tub. Distal", " neg" => "Tub. Proximal")
comprehensiveStates = [structuresDict[ismissing(row.pathology) ? split(row["SegmentDisplayName"],"|")[3] : row.pathology] for row in eachrow(features)]

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
                    for i in values(structuresDict)
                ],
                multi = true,
                value = first(values(structuresDict)),
            ),
            ##Add checkbox for healthy vs DKD
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
                id = "stage",
                options = [
                    (label = i, value = i)
                    for i in unique(features.SegmentDisplayName)
                ],
                multi = true,
                value =features[1,"SegmentDisplayName"],
            ),
            html_button(id = "g1", children = "Make group 1", n_clicks = 0),
            html_button(id = "g2", children = "Make group 2", n_clicks = 0),
            html_div(
                children = [
                    dcc_markdown("
                    **Group 1**
                    "),
                    html_pre(id = "g1-header"),
                ],
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
                    **Group 2**
                    "),
                    html_pre(id = "g2-header"),
                ],
            ),
            dcc_dropdown(
                id = "group2",
                options = [
                    (label = i, value = i)
                    for i in unique(features.SegmentDisplayName)
                ],
                multi = true,
                value =features[1,"SegmentDisplayName"],
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
) do roiList

    # df6f = df6[df6.year .== year_slider_value, :]

    plotData = [    ( x = Yte[1, comprehensiveStates .== status], y = Yte[2, comprehensiveStates .== status],  type = "scatter", name = status, mode = "markers", text = features[comprehensiveStates .== status,"SlideName"], customdata = features[comprehensiveStates .== status,"SegmentDisplayName"]) for status in roiList]

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
    Input("stage", "value"),
) do patientList, selected_data

    # df6f = df6[df6.year .== year_slider_value, :]
    
    sdList = [x in selected_data for x in features[!,"SegmentDisplayName"]]
    selector = "SlideName"
    #this list should be changed by the checkboxes
    plotData = [    ( x = Yte[1, (features[!,selector] .== status ) .& sdList ], y = Yte[2, (features[!,selector] .== status) .& sdList],  type = "scatter", name = status, mode = "markers", text = comprehensiveStates[(features[!,selector] .== status) .& sdList]) for status in patientList if any((features[!,selector] .== status) .& sdList)]

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
    Output("stage", "value"),
    Input("graph-1", "selectedData"),
) do selected_data
    selectedpoints = 1:2
    if selected_data != nothing
        selectedpoints = [p[:customdata] for p in selected_data.points]
        #print(selectedpoints)
    end
    #print(features[selectedpoints, "SegmentDisplayName"])
    #return features[selectedpoints, "SegmentDisplayName"]
    return selectedpoints
end

#buttons
callback!(
    app,
    Output("group1", "value"),
    Input("g1", "n_clicks"),
    State("stage", "value"),
) do clicks, input_1
    return input_1
end

callback!(
    app,
    Output("group2", "value"),
    Input("g2", "n_clicks"),
    State("stage", "value"),
) do clicks, input_1
    return input_1
end

run_server(app, "0.0.0.0", debug=true)