using DataFrames, Dash, DashHtmlComponents, DashCoreComponents, UrlDownload, PlotlyJS, JSON3
using Statistics
using LinearAlgebra
using CSV
using RCall

app = dash()

#calling R libraries

    R"""
    source("RScripts/preamble.R")
    """


#dataframes for R analysis
genon = DataFrame(CSV.File("RScripts/GO_final.txt"))
cellNew = DataFrame(CSV.File("RScripts/cell_new_final.txt"))


Yte=Matrix(CSV.read("PCA_matrix_log2.txt",DataFrame))
PCM = DataFrame(CSV.File("data/Kidney_Q3Norm_TargetCountMatrix.txt"))
print(PCM[2,2])
foreach(n -> PCM[!, n] = log2.(PCM[!, n]), names(PCM)[2:end])
print(PCM[2,2])
features = DataFrame(CSV.File("data/Kidney_Sample_Annotations.txt"))
structuresDict = Dict("abnormal"=>"Glom (Abnormal)","healthy"=>"Glom (Healthy)"," PanCK" => "Tub. Distal", " neg" => "Tub. Proximal")
comprehensiveStates = [structuresDict[ismissing(row.pathology) ? split(row["SegmentDisplayName"],"|")[3] : row.pathology] for row in eachrow(features)]
genes = PCM.TargetName
diseaseStatusMarkers = replace(features[:,"disease_status"], "DKD" => "triangle-up", "normal" => "circle")
# diseaseStatusMarkers .= relpace.(diseaseStatusMarkers, "DKD", "square")
insertcols!(features,"states"=>comprehensiveStates)
insertcols!(features, "barID" =>["$(row[1]) | $(row[2])" for row in eachrow( features[!,["ROILabel","states"]]) ])

CellTypes = DataFrame(CSV.File("data/Cell_Types_for_Spatial_Decon.txt"))
GO = DataFrame(CSV.File("data/Kidney_ssGSEA.txt"))
rename!(GO, names(GO)[2:end] .=> features.Sample_ID)
#group immune together
CellTypes[end-15:end,9].="Immune"
CDC = DataFrame(CSV.File("data/Kidney_Spatial_Decon.txt"))
rename!(CDC, names(CDC)[2:end] .=> features.Sample_ID)
insertcols!(CDC,  "cellGroup" => [CellTypes[ CellTypes.Alias .== row[1],9][1] for row in eachrow(CDC)])
to_group = names(CDC[!, Not(["Alias", "cellGroup"])])
newCDC = combine(groupby(CDC,"cellGroup"), to_group .=> sum .=> to_group)
cells = CDC.Alias
cellGroups = newCDC.cellGroup


diseaseOptions = [
    Dict("label" => "Normal", "value" => "normal"),
    Dict("label" => "DKD", "value" => "DKD"),
]

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
patientsNormal = patients[ occursin.("normal", patients) ]
patientsDKD = patients[ occursin.("disease", patients) ]

app.layout = html_div() do
    html_div(
        children = [

    	    html_h1(
				"Data Visualization App",
        		style = Dict("color" => "#0066cc", "textAlign" => "left"),
    		),
		    dcc_markdown("## Data selection"),
        	html_div(
        		children = [
        		html_div(
		            dcc_dropdown(
		                id = "crossfilter-roi-type",
		                options = [
		                    (label = i, value = i)
		                    for i in values(structuresDict)
		                ],
		                multi = true,
		                # value = first(values(structuresDict)),
		                placeholder = "Select type of structure",
		                style = (width = "70%", display = "block", marginBottom = "5px") 
		            ),
	            	# style = (width = "100%") 
	            ),
	            
	            ##Sync to patients in crossfilter-patient
	            html_div(
		            dcc_dropdown(
		                id = "crossfilter-patient",
		                options = [
		                    (label = i, value = i)
		                    for i in patients
		                ],
		                multi = true,
		                # value = [patients[1]],
		                placeholder = "Select patient",
		                style = (width = "70%", display = "block", marginBottom = "5px") 
		            ),
	            ),

	            ##Add checkbox for healthy vs DKD
	            html_button(id = "buttonAll", children = "All", n_clicks = 0),
	            html_button(id = "buttonDKD", children = "DKD", n_clicks = 0),
	            html_button(id = "buttonNormal", children = "Normal", n_clicks = 0),
	            html_button(id = "buttonClear", children = "Clear", n_clicks = 0),
	            ],
	            style = (borderBottom = "thin lightgrey solid", backgroundColor = "rgb(250, 250, 250)", padding = "10px 5px",),
            ),

        	html_div(
        		children = [
		            dcc_graph(
		                id = "graph-1",
		            ),
	            ],
                style = ( width = "49%", display = "inline-block"),
            ),

        	html_div(
        		children = [
		            dcc_graph(
		                id = "graph-2",
		            ),
	            ],
                style = ( width = "49%", display = "inline-block"),
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
                value =features[1:2,"SegmentDisplayName"],
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
                value =features[3:4,"SegmentDisplayName"],
            ),

            html_div(id = "dummy", children = " e"),

            html_button(id = "RRun", children = "Run R", n_clicks = 0),

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
                    for i in cellGroups
                ],
                value = cellGroups[1],
            ),

            dcc_graph(
                id = "bubblePlot2",
            ),

        	html_div(
        		children = [
                    dcc_graph(
                        id = "bar1",
                    ),
	            ],
                style = ( width = "49%", display = "inline-block"),
            ),

        	html_div(
        		children = [
                    dcc_graph(
                        id = "bar2",
                    ),
	            ],
                style = ( width = "49%", display = "inline-block"),
            ),
        ],
        style = (
            borderBottom = "thin lightgrey solid",
            backgroundColor = "rgb(250, 250, 250)",
            padding = "10px 5px",
        ),
    )
end


callback!(
    app,
    Output("crossfilter-patient", "value"),
    Input("buttonNormal", "n_clicks"),
    Input("buttonDKD", "n_clicks"),
    Input("buttonAll", "n_clicks"),
    Input("buttonClear", "n_clicks"),
    State("crossfilter-patient", "value"),
) do clickNormal, clickDKD, clickAll, clickClear, patientList
	
    ctx = callback_context()
    if isempty(ctx.triggered)
    	return patientList
    end

    trigger = ctx.triggered[1]
    button_id = split(trigger.prop_id, ".")[1]

    if button_id == "buttonAll"
    	return patients
    end
    if button_id == "buttonClear"
    	return []
    end
    if button_id == "buttonDKD"
    	#patientList = union(patientList,patientsDKD)
    	 patientList = patientsDKD
    end
    if button_id == "buttonNormal"
    	#patientList = union(patientList,patientsNormal)
    	 patientList = patientsNormal
    end

    return patientList

end



callback!(
    app,
    Output("graph-1", "figure"),
    Input("crossfilter-roi-type", "value"),
    # Input("crossfilter-disease-status", "value"),
    Input("crossfilter-patient", "value"),
) do roiList, patientList

	# if no selection return empty plot
	if roiList == nothing
        return (data = [], layout = (title = "By Structure", xaxis_title = "", yaxis_title = "", clickmode = "event+select",))
    end

    # diseaseFilter = true
    # if !isempty(diseaseStatus)
    # 	for ds in diseaseStatus
    # 		diseaseFilter = (features[!,"disease_status"] .== ds) .| diseaseFilter
    # 	end
    # end


    patientFilter = false
    if patientList != nothing
    	patientFilter = [slide in patientList for slide in features.SlideName]
    end

    accumFilter = patientFilter # .& diseaseFilter 

    plotData = [    ( x = Yte[1, (comprehensiveStates .== status) .& accumFilter], y = Yte[2, (comprehensiveStates .== status) .& accumFilter],  type = "scatter", name = status, mode = "markers", marker = (size = 8, symbol = diseaseStatusMarkers[ (comprehensiveStates .== status) .& accumFilter ], ), text = features[(comprehensiveStates .== status) .& accumFilter,"SlideName" ], customdata = features[(comprehensiveStates .== status) .& accumFilter,"SegmentDisplayName" ]) for status in roiList]

    return (
        data = plotData,
        layout = (
            title = "By Structure",
            xaxis = (title = "pca_1",),
            yaxis = (title = "pca_2",),
            clickmode = "event+select",
        ),
    )

end

callback!(
    app,
    Output("graph-2", "figure"),
    Input("crossfilter-patient", "value"),
    Input("stage", "value"),
    Input("crossfilter-roi-type", "value"),
) do patientList, selected_data, roiList

	# if no selection return empty plot
	if patientList == nothing
        return ( data = [], layout =  "title" => "By Patient"  )
    end

    
    selectedFilter = true
    # if selected_data != nothing
    # 	selectedFilter = [x in selected_data for x in features[!,"SegmentDisplayName"]]
    # end

	roiFilter = false
	if roiList != nothing
    	roiFilter = [x in roiList for x in comprehensiveStates]
    end

    accumFilter = selectedFilter .& roiFilter
    # accumFilter = true

    selector = "SlideName"

    plotData = [    ( x = Yte[1, (features[!,selector] .== status ) .& accumFilter ], y = Yte[2, (features[!,selector] .== status) .& accumFilter],  type = "scatter", name = status, mode = "markers", text = comprehensiveStates[(features[!,selector] .== status) .& accumFilter]	) for status in patientList]

    return (
        data = plotData,
        layout = (
            title = "By Patient",
            xaxis = (title = "pca_1",),
            yaxis = (title = "pca_2",),
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
    Output("dummy", "children"),
    Input("RRun", "n_clicks"),
    State("group1", "value"),
    State("group2", "value"),
) do clicks, input_1, input_2
    clicks == 0 && return " "
    CSV.write("matrix_group2_group1.txt",rename(PCM,names(PCM)[1] => "gene")[!,vcat(["gene"], input_1, input_2)], delim='\t')
    n1 = length(input_1)
    n2 = length(input_2)
    @rput n1
    @rput n2
    CSV.write("matrix_Design_final.txt", DataFrame("sample" => vcat(input_1, input_2), "group1" => convert(Vector{UInt8},vcat(trues(n1),falses(n2))), "group2" => convert(Vector{UInt8},vcat(falses(n1),trues(n2)))), delim = '\t')
groupListsInPatients= [ 
    [ any(x .== input_1) for x in features[!,"SegmentDisplayName"]], 
    [ any(x .== input_2) for x in features[!,"SegmentDisplayName"]] 
    ]

    cnf1 = DataFrame((newCDC.cellGroup .=> [Vector(row) for row in eachrow( newCDC[!,features[groupListsInPatients[1],"Sample_ID"]])])... , "Group" => "Group1")
    cnf2 = DataFrame((newCDC.cellGroup .=> [Vector(row) for row in eachrow( newCDC[!,features[groupListsInPatients[2],"Sample_ID"]])])... , "Group" => "Group2")
    CSV.write("cell_new_final.txt", vcat(cnf1,cnf2), delim = '\t')

    gf1 = DataFrame((GO.Column1 .=> [Vector(row) for row in eachrow( GO[!,features[groupListsInPatients[1],"Sample_ID"]])])... , "Group" => "Group1")
    gf2 = DataFrame((GO.Column1 .=> [Vector(row) for row in eachrow( GO[!,features[groupListsInPatients[2],"Sample_ID"]])])... , "Group" => "Group2")
    CSV.write("GO_final.txt", vcat(gf1,gf2), delim = '\t')

    R"""
    source("RScripts/de.R")
    """


    return "DE and all"
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
#print(any(any.(groupListsInPatients)))
normalizedCounts =Any[ ]

maxCounts = mean(Vector(PCM[PCM[!,"TargetName" ] .== selectedGene,points[!,:].SegmentDisplayName][1,:]))
for group in groupListsInPatients
    if any(group)
        tmp = Vector(PCM[PCM[!,"TargetName" ] .== selectedGene,points[group,:].SegmentDisplayName][1,:])
        #print(tmp)
        push!(normalizedCounts,50 .* tmp./maxCounts) 
    end
end
#Plots.scatter(points.ROICoordinateX,points.ROICoordinateY,color=1,label="Tubules")
#Plots.title!(patient)
#Plots.scatter!(points[gloms,:].ROICoordinateX,points[gloms,:].ROICoordinateY,color=2,label="Gloms")
#print(normalizedCounts)

####PENDING#####
#text = structuresDict
    plotData = [    ( x =points[group,:].RoiReportX, y = points[group,:].RoiReportY,  type = "scatter", name = string(i), mode = "markers", marker = (size=normalizedCounts[i], symbol = "circle", ), text = points[group,:].states, customdata = points[group,:].SegmentDisplayName )  for(i,  group) in enumerate( groupListsInPatients) if any(group)]

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
#print(any(any.(groupListsInPatients)))
normalizedCounts =Any[ ]

maxCounts = maximum(Vector(newCDC[newCDC[!,"cellGroup" ] .== selectedCell,points[!,:].Sample_ID][1,:]))
for group in groupListsInPatients
    if any(group)
        tmp = Vector(newCDC[newCDC[!,"cellGroup" ] .== selectedCell,points[group,:].Sample_ID][1,:])
        #print(tmp)
        push!(normalizedCounts,50 .* tmp./maxCounts) 
    end
end

####PENDING#####
#text = structuresDict
    plotData = [    ( x =points[group,:].RoiReportX, y = points[group,:].RoiReportY,  type = "scatter", name = string(i), mode = "markers", marker = (size=normalizedCounts[i], symbol = "circle", ), text = points[group,:].states, customdata = points[group,:].SegmentDisplayName )  for(i,  group) in enumerate( groupListsInPatients) if any(group)]

    return ( data = plotData, layout = ( title = "Regions in Slide", xaxis_title = "x", yaxis_title = "y",),)

end

callback!(
    app,
    Output("bar1", "figure"),
    Input("patientsBubble", "value"),
    #Input("cellBubble", "value"),
    Input("group1", "value"),
    #Input("group2", "value"),
) do patient,  group1
#patientRegions= features[!,"SlideName"].==patient
points = features[!,["SlideName", "SegmentDisplayName", "RoiReportX", "RoiReportY", "states","barID","Sample_ID"]]
group=  [ any(x .== group1) for x in points[!,"SegmentDisplayName"]] 
    
#print(groupListsInPatients)


    plotData = [    ( x =points[group,:].SegmentDisplayName, y = Vector(newCDC[newCDC[!,"cellGroup" ] .== cell,points[group,:].Sample_ID][1,:]),  type = "bar", name = cell, text = cell, customdata = cell ) for cell in cellGroups  ]
    

    return ( data = plotData, layout = ( title = "Group Two",barmode = "stack",
            xaxis = (title = "ROIs",),
            yaxis = (title = "CD Cell Fraction",),
    )
    )

end

callback!(
    app,
    Output("bar2", "figure"),
    Input("patientsBubble", "value"),
    #Input("cellBubble", "value"),
    Input("group2", "value"),
) do patient,  group1
#patientRegions= features[!,"SlideName"].==patient
points = features[!,["SlideName", "SegmentDisplayName", "RoiReportX", "RoiReportY", "states","barID","Sample_ID"]]
group=  [ any(x .== group1) for x in points[!,"SegmentDisplayName"]] 
    
#print(groupListsInPatients)


    plotData = [    ( x =points[group,:].SegmentDisplayName, y = Vector(newCDC[newCDC[!,"cellGroup" ] .== cell,points[group,:].Sample_ID][1,:]),  type = "bar", name = cell, text = cell, customdata = cell ) for cell in cellGroups  ]
    

    return ( data = plotData, layout = ( title = "Group Two",barmode = "stack",
            xaxis = (title = "ROIs",),
            yaxis = (title = "CD Cell Fraction",),
    showlegend = false),)

end

run_server(app, "0.0.0.0", debug=true)