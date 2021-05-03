#Pkg.add("DataFrames")
 #Pkg.add("Dash")
 #Pkg.add("DashHtmlComponents")
 #Pkg.add("DashCoreComponents")
  #Pkg.add("UrlDownload")
   #Pkg.add("PlotlyJS")
    #Pkg.add("JSON3")
#Pkg.add("Statistics")
#Pkg.add("LinearAlgebra")
#Pkg.add("CSV")
#Pkg.add("RCall")
using DataFrames, Dash, DashHtmlComponents, DashCoreComponents, UrlDownload, PlotlyJS, JSON3
using Statistics
using LinearAlgebra
using CSV

app = dash()

features = DataFrame(CSV.File("data/Kidney_Sample_Annotations.txt"))
Yte=Matrix(CSV.read("PCA_matrix.txt",DataFrame))
PCM = DataFrame(CSV.File("data/Kidney_Q3Norm_TargetCountMatrix.txt"))
CDC = DataFrame(CSV.File("data/Kidney_Spatial_Decon.txt"))
rename!(CDC, names(CDC)[2:end] .=> features.Sample_ID)

structuresDict = Dict("abnormal"=>"Glom (Abnormal)","healthy"=>"Glom (Healthy)"," PanCK" => "Tub. Distal", " neg" => "Tub. Proximal")
comprehensiveStates = [structuresDict[ismissing(row.pathology) ? split(row["SegmentDisplayName"],"|")[3] : row.pathology] for row in eachrow(features)]
genes = PCM.TargetName
cells = CDC.Alias
insertcols!(features,"states"=>comprehensiveStates)

function generate_table(dataframe, max_rows = 10)
    html_table([
        html_thead(html_tr([html_th(col) for col in names(dataframe)])),
        html_tbody([
            html_tr([html_td(dataframe[r, c]) for c in names(dataframe)]) for r = 1:min(nrow(dataframe), max_rows)
        ]),
    ])
end


app = dash()

#available_indicators = ["disease_status","region","pathology"]]

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
            ##Sync to patients in crossfilter-patient
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
            #Group analysis
            html_div(
                children = [
                    dcc_markdown("
                    **BubblePlot**
                    "),
                    html_pre(id = "bubble-header"),
                ],
            ),
            #restrict to patients present in group1 intersection group2
            dcc_dropdown(
                id = "patientsBubble",
                options = [
                    (label = i, value = i)
                    for i in patients
                ],
                value =patients[1],
            ),
            dcc_dropdown(
                id = "geneBubble",
                options = [
                    (label = i, value = i)
                    for i in genes
                ],
                value = genes[1],
            ),

            dcc_graph(
                id = "bubblePlot1",
            ),
            dcc_dropdown(
                id = "cellBubble",
                options = [
                    (label = i, value = i)
                    for i in cells
                ],
                value = cells[1],
            ),

            dcc_graph(
                id = "bubblePlot2",
            ),

            dcc_graph(
                id = "bar",
            ),
        ]
    )
end

callback!(
    app,
    Output("graph-1", "figure"),
    Input("crossfilter-roi-type", "value"),
) do roiList

    # df6f = df6[df6.year .== year_slider_value, :]

    plotData = [    ( x = Yte[1, comprehensiveStates .== status], y = Yte[2, comprehensiveStates .== status],  type = "scatter", name = status, mode = "markers", text = features[comprehensiveStates .== status,"SlideName"], customdata = features[comprehensiveStates .== status,"SegmentDisplayName"]) for status in roiList]

    return (
        data = plotData,
        layout = (
            title = "By Structure",
            xaxis= (title = "pca_1", ),
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

callback!(
    app,
    Output("bubblePlot1", "figure"),
    Input("patientsBubble", "value"),
    Input("geneBubble", "value"),
    Input("group1", "value"),
    Input("group2", "value"),
) do patient, selectedGene, group1,group2
patientRegions= features[!,"SlideName"].==patient
points = features[patientRegions,["SlideName", "SegmentDisplayName", "RoiReportX", "RoiReportY", "states"]]
groupListsInPatients= [ 
    [ any(x .== group1) for x in points[!,"SegmentDisplayName"]], 
    [ any(x .== group2) for x in points[!,"SegmentDisplayName"]] 
    ]
#print(groupListsInPatients)
normalizedCounts =Any[ ]

maxCounts = mean(Vector(PCM[PCM[!,"TargetName" ] .== selectedGene,points[!,:].SegmentDisplayName][1,:]))
for group in groupListsInPatients
    tmp = Vector(PCM[PCM[!,"TargetName" ] .== selectedGene,points[group,:].SegmentDisplayName][1,:])
    #print(tmp)
    push!(normalizedCounts,50 .* tmp./maxCounts) 
end
#Plots.scatter(points.ROICoordinateX,points.ROICoordinateY,color=1,label="Tubules")
#Plots.title!(patient)
#Plots.scatter!(points[gloms,:].ROICoordinateX,points[gloms,:].ROICoordinateY,color=2,label="Gloms")
#print(normalizedCounts)

####PENDING#####
#text = structuresDict
    plotData = [    ( x =points[group,:].RoiReportX, y = points[group,:].RoiReportY,  type = "scatter", name = string(i), mode = "markers", marker = (size=normalizedCounts[i], symbol = "circle", ), text = points[group,:].states, customdata = points[group,:].SegmentDisplayName )  for(i,  group) in enumerate( groupListsInPatients)]

    return ( data = plotData, layout = ( title = "Regions in Slide", xaxis_title = "x", yaxis_title = "y",),)

end


callback!(
    app,
    Output("bubblePlot2", "figure"),
    Input("patientsBubble", "value"),
    Input("cellBubble", "value"),
    Input("group1", "value"),
    Input("group2", "value"),
) do patient, selectedCell, group1,group2
patientRegions= features[!,"SlideName"].==patient
points = features[patientRegions,["SlideName", "SegmentDisplayName", "RoiReportX", "RoiReportY", "states","Sample_ID"]]
groupListsInPatients= [ 
    [ any(x .== group1) for x in points[!,"SegmentDisplayName"]], 
    [ any(x .== group2) for x in points[!,"SegmentDisplayName"]] 
    ]
#print(groupListsInPatients)
normalizedCounts =Any[ ]

maxCounts = maximum(Vector(CDC[CDC[!,"Alias" ] .== selectedCell,points[!,:].Sample_ID][1,:]))
for group in groupListsInPatients
    tmp = Vector(CDC[CDC[!,"Alias" ] .== selectedCell,points[group,:].Sample_ID][1,:])
    #print(tmp)
    push!(normalizedCounts,50 .* tmp./maxCounts) 
end
#Plots.scatter(points.ROICoordinateX,points.ROICoordinateY,color=1,label="Tubules")
#Plots.title!(patient)
#Plots.scatter!(points[gloms,:].ROICoordinateX,points[gloms,:].ROICoordinateY,color=2,label="Gloms")
#print(normalizedCounts)

####PENDING#####
#text = structuresDict
    plotData = [    ( x =points[group,:].RoiReportX, y = points[group,:].RoiReportY,  type = "scatter", name = string(i), mode = "markers", marker = (size=normalizedCounts[i], symbol = "circle", ), text = points[group,:].states, customdata = points[group,:].SegmentDisplayName )  for(i,  group) in enumerate( groupListsInPatients)]

    return ( data = plotData, layout = ( title = "Regions in Slide", xaxis_title = "x", yaxis_title = "y",),)

end

callback!(
    app,
    Output("bar", "figure"),
    Input("patientsBubble", "value"),
    #Input("cellBubble", "value"),
    Input("group1", "value"),
    #Input("group2", "value"),
) do patient,  group1
patientRegions= features[!,"SlideName"].==patient
points = features[patientRegions,["SlideName", "SegmentDisplayName", "RoiReportX", "RoiReportY", "states","Sample_ID"]]
group=  [ any(x .== group1) for x in points[!,"SegmentDisplayName"]] 
    
#print(groupListsInPatients)


    plotData = [    ( x =points[group,:].SegmentDisplayName, y = Vector(CDC[CDC[!,"Alias" ] .== cell,points[group,:].Sample_ID][1,:]),  type = "bar", name = cell, text = cell, customdata = cell ) for cell in cells  ]
    

    return ( data = plotData, layout = ( title = "baar",barmode = "stack",),)

end
run_server(app, "0.0.0.0", debug=true)