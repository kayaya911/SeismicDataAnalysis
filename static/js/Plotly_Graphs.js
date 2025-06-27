// stricter parsing and error handling
"use strict";

function ChannelList_UniqueID(ID) {
    let Indx = ChannelList.findIndex(x => x.Unique_ID === ID);
    return Indx;
}

function Channel_Click(CB) {

    // Determine which channel is clicked
    let Indx    = ChannelList_UniqueID(CB.id);
    let Div_ID  = "Div_ID_" + ChannelList[Indx].Unique_ID;

    // Show or Hide the related Div-Graph 
    if (CB.checked) {
        document.getElementById(Div_ID).style.display = 'flex';
        ChannelList[Indx].Selected = true;
        RePosition_Graph(Indx);
    }
    else {
        document.getElementById(Div_ID).style.display = 'none';
        ChannelList[Indx].Selected = false;
    }

    if      (PageNo == 0) { console.log("Page No: " + PageNo);}
    else if (PageNo == 1) { console.log("Page No: " + PageNo);}
    else if (PageNo == 2) { console.log("Page No: " + PageNo);}
    else if (PageNo == 3) { console.log("Page No: " + PageNo);}
    else if (PageNo == 4) { console.log("Page No: " + PageNo);}

}

function RePosition_Graph(ChNum) {
    //Decleration of variables 
    let Div_ID  = "Div_ID_" + ChannelList[ChNum].Unique_ID;
    let a1 = document.getElementById(Div_ID).offsetTop;
    let a2 = document.getElementById("Parameter_Container").offsetHeight;

    // RePossition the Scrollbar of the container-div so that the Graph is located at the top of container-div 
    document.getElementById("Graphs_LoadData_Container").scrollTop = a1 - a2 - 100;
}

function Stats_Table_Units_Update(Channel, a, opt) {

    // Check variable 
    if (opt == null) {opt=true;}
    
    // Clear the content of cell in the table 
    let Unit_Cell_ID = "Unit_Cell_ID_" + Channel.Unique_ID;
    let Unit_Plot_ID = "Unit_Plot_ID_" + Channel.Unique_ID;
    let Previous_Ind = document.getElementById(Unit_Plot_ID).selectedIndex;
    let Unit_List    = UnitList(a);
    document.getElementById(Unit_Cell_ID).innerHTML = "";
    
    // Create the select element using (Unit_List.Units)
    let select = Units_Select_Element(Channel, Unit_List.Units, Unit_Plot_ID);
    select.setAttribute("onchange", "Plot_Update("+ Channel.ChNum +")" );
    
    // Automatically select the index of the select-element if applicable 
    if (opt) { select.selectedIndex = Unit_List.UnitNum.indexOf(a);}
    else     { select.selectedIndex = Previous_Ind; }

    // Asign select to cell element 
    document.getElementById(Unit_Cell_ID).appendChild(select);
}

function Create_Plotly_Graph(Container_Id, Channel) {

    // Decleration of variables 
    let Indx_ChNum, Indx, Container_Div, Channel_Info, SubWrapper, PlotArea, PlotArea_ID, Div_ID, PlotInfoArea, Span_Title, Span_ChannelInfo, Table_Statistics;

    Indx_ChNum = ChannelList_UniqueID(Channel.Unique_ID);
    Indx = "("+ (Indx_ChNum+1).toString().padStart(3, '0') + ")";

    Container_Div = document.getElementById(Container_Id);
    Channel_Info  = "(Ch-" + Channel.ChNum + ") (" +  Channel.Orientation  + ") (" + Number(Channel.FSamp).toFixed(3).toString() + " Hz)";
    Div_ID        = "Div_ID_" + Channel.Unique_ID;
    PlotArea_ID   = "PlotArea_ID_" + Channel.Unique_ID;

    // Create a div element
    SubWrapper = document.createElement('div');
    SubWrapper.setAttribute('id', Div_ID);
    SubWrapper.setAttribute('class', 'Plotly_Main_Container');

    // Plotly Graph_Div
    PlotArea    = document.createElement('div');
    PlotArea.id =  PlotArea_ID;
    PlotArea.className = "Plotly_Div";

    // Plotly Information_Div
    PlotInfoArea = document.createElement('div');
    PlotInfoArea.className = "Plotly_Info_Div";

    // Span_Title (File Name)
    Span_Title = document.createElement('span');
    Span_Title.textContent = Channel.FileName;
    Span_Title.className = "Plotly_Span_Title";

    // Span_ChannelInfo (Ch-Num) (Azimut) (Sampling)
    Span_ChannelInfo = document.createElement('span');
    Span_ChannelInfo.textContent = Indx + " " + Channel_Info;
    Span_ChannelInfo.className = "Plotly_Span_ChannelInfo";3

    // Table_Statistics
    Table_Statistics = Table_Statisctics(Channel);

    // Append all div
    Container_Div.appendChild(SubWrapper);
    SubWrapper.appendChild(PlotArea);
    SubWrapper.appendChild(PlotInfoArea)
    PlotInfoArea.appendChild(Span_Title);
    PlotInfoArea.appendChild(Span_ChannelInfo);
    PlotInfoArea.appendChild(Table_Statistics);

    // Create the plot after the DOM is created 
    Plot_Create(PlotArea, Channel); 
    
}

function Table_Statisctics(Channel) {

    let Table_Statistics, Table_Body, row, cell, opt;
    let Statictics_Peak_ID, Statictics_Mean_ID, Statictics_RMS_ID;
    let FilterType_ID, FilterRow_ID, BaseLine_ID, BaseLineRow_ID, FilterResp_ID, GraphUnitRow_ID,FilterFFT_ID;
    let Span, input, div1, div2, label, Unit_List, j, Unit_Cell_ID, Unit_Plot_ID, select;

    Statictics_Peak_ID    = "Statictics_Peak_ID_" + Channel.Unique_ID;
    Statictics_Mean_ID    = "Statictics_Mean_ID_" + Channel.Unique_ID;
    Statictics_RMS_ID     = "Statictics_RMS_ID_" + Channel.Unique_ID;
    Unit_Cell_ID          = "Unit_Cell_ID_" + Channel.Unique_ID;
    BaseLine_ID           = "BaseLine_ID" + Channel.Unique_ID;
    BaseLineRow_ID        = "BaseLineRow_ID" + Channel.Unique_ID;
    FilterType_ID         = "FilterType_ID_" + Channel.Unique_ID;
    FilterRow_ID          = "FilterRow_ID_" + Channel.Unique_ID;
    FilterResp_ID         = "FilterResp_ID_" + Channel.Unique_ID;
    FilterFFT_ID          = "FilterFFT_ID_" + Channel.Unique_ID;
    GraphUnitRow_ID       = "GraphUnitRow_ID_" + Channel.Unique_ID;

    Table_Statistics = document.createElement('table');
    Table_Statistics.setAttribute('class', 'Plotly_Stat_Table')

    Table_Body = document.createElement('tbody');
    Table_Body.setAttribute('class', '')
    Table_Statistics.appendChild(Table_Body);

    // Add row to Table_Body
    row = Table_Statistics.insertRow(-1);
    row.setAttribute('class', '');
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Peak";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.innerHTML = Channel.Peak.toPrecision(4);
    cell.id        = Statictics_Peak_ID;

    // Assign class to the row
    row = Table_Statistics.insertRow(-1);
    row.setAttribute('class', '');
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Mean";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.innerHTML = Channel.Mean.toPrecision(4);  
    cell.id        = Statictics_Mean_ID;

    // Assign class to the row
    row = Table_Statistics.insertRow(-1);
    row.setAttribute('class', '');
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "RMS";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.innerHTML = Channel.RMS.toPrecision(4);
    cell.id        = Statictics_RMS_ID;

    // Assign class to the row
    row = Table_Statistics.insertRow(-1);
    row.setAttribute('class', '');
    row.setAttribute('id', GraphUnitRow_ID);
    row.style.display = "none";
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Graph Unit";
    
    Unit_List = UnitList(TypeAndUnit(Channel.TypeAndUnits).Unit);
    Unit_Plot_ID   = "Unit_Plot_ID_" + Channel.Unique_ID;
    select  = document.createElement('select');
    select.setAttribute("id", Unit_Plot_ID);
    select.setAttribute('class', 'form-select form-select-sm');
    for (j = 0; j < Unit_List.Units.length; j++) {
        opt = document.createElement("option");
        opt.value = Unit_List.Units[j];
        opt.text = Unit_List.Units[j];
        select.add(opt, null);
    }
    select.setAttribute("onchange", "Plot_Update("+ Channel.ChNum +")" );
    select.selectedIndex = Unit_List.UnitNum.indexOf(Channel.Unit) ;
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.id        = Unit_Cell_ID;
    cell.appendChild(select);

    // Assign class to the row
    row = Table_Statistics.insertRow(-1);
    row.setAttribute('class', '');
    row.setAttribute('id', BaseLineRow_ID);
    row.style.display = "none";
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Baseline";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.innerHTML = "";
    cell.id        = BaseLine_ID;
    
    // Assign class to the row
    Span = document.createElement('span');
    Span.setAttribute('class', 'Filter_Div_Span');
    Span.textContent = 'Filter';
    input = document.createElement('input');
    input.setAttribute('class', 'toggle');
    input.setAttribute('type', 'checkbox');
    input.setAttribute('id', FilterResp_ID);
    input.setAttribute('onclick', 'FilterResponse(this)');
    label = document.createElement('label');
    label.setAttribute('for', FilterResp_ID);
    label.setAttribute('class', "toggle-label");
    div1 = document.createElement('div');
    div1.setAttribute('class', "Filter_Div");
    div1.appendChild(Span);
    div1.appendChild(input);
    div1.appendChild(label);

    Span = document.createElement('span');
    Span.setAttribute('class', 'Filter_Div_Span');
    Span.textContent = 'FFT';
    input = document.createElement('input');
    input.setAttribute('class', 'toggle');
    input.setAttribute('type', 'checkbox');
    input.setAttribute('id', FilterFFT_ID);
    input.setAttribute('onclick', 'FilterResponse(this)');
    label = document.createElement('label');
    label.setAttribute('for', FilterFFT_ID);
    label.setAttribute('class', "toggle-label");
    div2 = document.createElement('div');
    div2.setAttribute('class', "Filter_Div");
    div2.appendChild(Span);
    div2.appendChild(input);
    div2.appendChild(label);

    row = Table_Statistics.insertRow(-1);
    row.setAttribute('class', '');
    row.setAttribute('id', FilterRow_ID);
    row.style.display = "none";
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.appendChild(div1);
    cell.appendChild(div2);
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.innerHTML = "";
    cell.id        = FilterType_ID;

    return Table_Statistics;
}

function Plot_Create(Div_ID, Channel) {

    // Declaration of variables
    let traces=[], trace1, trace2, trace3, layout, config;
    let yLabel  = "<b>" + Channel.TypeString + "   [" + Channel.UnitString + "] <b>";
    let y2Label = "<b><b>";
    let xLabel  = "<b><b>";

    // Define Traces (Raw Data and Filtered Data)
    trace1 = {
        x             : Channel.time,
        y             : Mult(Channel.data, Channel.ScaleFactor),
        mode          : 'lines',
        type          : 'scatter',
        yaxis         : "y1",
        name          : '<b>Raw Data<b>',   // Legend name
        opacity       : 1.0,
        visible       : true,  // Show this trace
        line          : {color: 'blue', width: 1.00, dash: 'solid' },
        showlegend    : false,
    };
    trace2 = {
        x             : [],
        y             : [],   // This should be   Channel.Results.Filter.data
        mode          : 'lines',
        type          : 'scatter',
        yaxis         : "y1",
        name          : '<b>Filtered Data<b>',  // Legend name
        opacity       : 1.0,
        visible       : false,   // Hide this trace
        line          : {color: 'red', width: 1.00, dash: 'solid' },
        showlegend    : false,
    };
    trace3 = {
        x             : [],
        y             : [],   // This should be   Channel.Results.Filter.data
        mode          : 'lines',
        type          : 'scatter',
        yaxis         : "y2",
        name          : '<b>Phase<b>',  // Legend name
        opacity       : 1.0,
        visible       : false,   // Hide this trace
        line          : {color: 'green', width: 1.00, dash: 'solid' },
        showlegend    : false,
    };
    if (Channel.Results.Filter.hasOwnProperty("FilteredData")) {
        trace2.x = Channel.time,
        trace2.x = Channel.Results.Filter.FilteredData;
    }
    traces.push(trace1);
    traces.push(trace2);
    traces.push(trace3);

    // Define Layout
    layout = {
        xaxis           : { zeroline: false, automargin: true, tickfont: { size: 15 },                    linecolor: 'black', linewidth: 1, mirror: true, title: {text: xLabel,  standoff: 5, font: {family: "Arial", size: 17} }, autorange: true },
        yaxis           : { zeroline: true,  automargin: true, tickfont: { size: 15 }, tickformat: '.2e', linecolor: 'black', linewidth: 1, mirror: true, title: {text: yLabel,  standoff: 5, font: {family: "Arial", size: 17} }, autorange: true },
        yaxis2          : { zeroline: false, automargin: true, tickfont: { size: 15 },                    linecolor: 'black', linewidth: 1, mirror: true, title: {text: y2Label, standoff: 5, font: {family: "Arial", size: 17} }, autorange: true, overlaying: 'y', side: 'right', showticklabels: false },
        plot_bgcolor    : 'rgb(255,255,255)',
        paper_bgcolor   : 'rgb(255,255,255)',
        legend          : { x: 1, y:0.85, xanchor: 'right', orientation: 'v', font: {size: 14}},
        autosize        : true,
        margin          : {t: 20, r:20, b:5, l:5},
        shapes          : [],
    };

    config = {
        responsive              : true,
        displayModeBar          : true,    // show-hide floating toolbar all together 
        modeBarButtonsToRemove  : [],
        displaylogo             : false,   // Romoves the Ployly logo from toolbar
        useResizeHandler        : true,    // Enables Plotly's resize event listener
        showTips                : false,
        scrollZoom              : false,   // Enable mouse wheel zooming
    }
    // Display using Plotly
    Plotly.newPlot(Div_ID, traces, layout, config).then( function(gd) {
        // Get the ModeBar element
        var modeBar = gd.querySelector('.modebar');

        // Hide ModeBar initially
        modeBar.style.display = 'none';

        // Add mouse event listeners to the graph div
        gd.addEventListener('mouseenter', function() {
            modeBar.style.display = 'block';
        });

        gd.addEventListener('mouseleave', function() {
            modeBar.style.display = 'none';
        });
    });

}

async function Plot_Update(ChNum) {
    
    // Declaration of variables
    let layout_update, traces;
    let res     = GraphUnits(ChNum);
    let yLabel  = "<b>" + ChannelList[ChNum].TypeString + "   [" + ChannelList[ChNum].UnitString + "] <b>";
    let y2Label = "<b>Phase<b>";
    let xLabel  = "<b><b>"
    let Div_ID  = "PlotArea_ID_" + ChannelList[ChNum].Unique_ID;
    let PlotArea_ID   = "Div_ID_" + ChannelList[ChNum].Unique_ID;

    let Statictics_Peak_ID = "Statictics_Peak_ID_" + ChannelList[ChNum].Unique_ID;
    let Statictics_Mean_ID = "Statictics_Mean_ID_" + ChannelList[ChNum].Unique_ID;
    let Statictics_RMS_ID  = "Statictics_RMS_ID_" + ChannelList[ChNum].Unique_ID;
    let BaseLine_ID        = "BaseLine_ID" + ChannelList[ChNum].Unique_ID;
    let BaseLineRow_ID     = "BaseLineRow_ID" + ChannelList[ChNum].Unique_ID;
    let FilterType_ID      = "FilterType_ID_" + ChannelList[ChNum].Unique_ID;
    let FilterRow_ID       = "FilterRow_ID_" + ChannelList[ChNum].Unique_ID;
    let FilterResp_ID      = "FilterResp_ID_" + ChannelList[ChNum].Unique_ID;
    let FilterFFT_ID       = "FilterFFT_ID_" + ChannelList[ChNum].Unique_ID;
    let GraphUnitRow_ID    = "GraphUnitRow_ID_" + ChannelList[ChNum].Unique_ID;
    let Unit_Cell_ID       = "Unit_Cell_ID_" + ChannelList[ChNum].Unique_ID;
    let Unit_Plot_ID       = "Unit_Plot_ID_" + ChannelList[ChNum].Unique_ID;

    // Get the existing traces in the Div_ID, then update it
    traces = document.getElementById(Div_ID).data;

    // Get the existing layout in the Div_ID, then update it
    layout_update = document.getElementById(Div_ID).layout;
    layout_update = {
        xaxis  : { zeroline: false, automargin: true, tickfont: { size: 15 },                    linecolor: 'black', linewidth: 1, mirror: true, title: {text: xLabel,  standoff: 5, font: {family: "Arial", size: 17} }, autorange: true },
        yaxis  : { zeroline: true,  automargin: true, tickfont: { size: 15 }, tickformat: '.2e', linecolor: 'black', linewidth: 1, mirror: true, title: {text: yLabel,  standoff: 5, font: {family: "Arial", size: 17} }, autorange: true },
        yaxis2 : { zeroline: false, automargin: true, tickfont: { size: 15 },                    linecolor: 'green', linewidth: 1, mirror: true, title: {text: y2Label, standoff: 5, font: {family: "Arial", size: 17} }, autorange: true, overlaying: 'y', side: 'right', showticklabels: true },
    };

    // Update Graphe for this channel
    if (PageNo == 0) {

        // Update Filter row, BaseLine row and GraphUnit row in table 
        document.getElementById(FilterType_ID).innerHTML        = "";
        document.getElementById(BaseLine_ID).innerHTML          = "";

        document.getElementById(GraphUnitRow_ID).style.display  = "none";
        document.getElementById(FilterRow_ID).style.display     = "none";
        document.getElementById(BaseLineRow_ID).style.display   = "none";

        document.getElementById(FilterType_ID).className        = "Plotly_Stat_Body_Td_Right";

        // Update trace for rawData
        traces[0].x = ChannelList[ChNum].time;
        traces[0].y = Mult( ChannelList[ChNum].data,  (ChannelList[ChNum].ScaleFactor * res.ScaleFactor) );

        // Update the Statistics of RawData in table - scaled to user-specified unit
        StatTableUpdate(ChannelList[ChNum], res);

        // Update trace infomation
        traces[0].visible = true;
        traces[0].opacity = 1.00;
        traces[0].line    = {color: 'blue', width: 1.50, dash: 'solid' };
        traces[0].name    = '<b>Raw Data<b>';

        traces[1].visible = false;
        traces[1].opacity = 0.35;
        traces[1].line    = {color: 'grey', width: 1.00, dash: 'solid' };
        traces[1].name    = '<b>Filtered Data<b>';

        traces[2].visible = false;
        traces[2].opacity = 0.35;
        traces[2].line    = {color: 'grey', width: 1.00, dash: 'solid' };
        traces[2].name    = '<b>Phase<b>';

        // Update y axis of Graph
        layout_update.yaxis.title.text      = res.yLabel;
        layout_update.yaxis2.showticklabels = false;
        layout_update.yaxis2.title.text     = "";
    }
    else if (PageNo == 1) {

        // Update Filter row, BaseLine row and GraphUnit row in table 
        document.getElementById( Statictics_Peak_ID ).innerHTML = "";
        document.getElementById( Statictics_Mean_ID ).innerHTML = "";
        document.getElementById( Statictics_RMS_ID  ).innerHTML = "";
        document.getElementById(FilterType_ID).innerHTML        = "";
        document.getElementById(BaseLine_ID).innerHTML          = "";
        document.getElementById(GraphUnitRow_ID).style.display  = "table-row";

        // Update the trace of RawData
        traces[0].x = ChannelList[ChNum].time;
        traces[0].y = Detrend(Mult(ChannelList[ChNum].data, (ChannelList[ChNum].ScaleFactor * res.ScaleFactor)), 0);
        
        // Update the trace of FilteredData as empty array
        traces[1].x = [],
        traces[1].y = [];

        // Continue only if Channel.Results.Filter has the attributes of "FilteredData"
        // This mean this channel contains filtered data
        if ( ChannelList[ChNum].Results.Filter.hasOwnProperty("FilteredData") ) { 

            // Continue if FilteredData contains data
            if  (ChannelList[ChNum].Results.Filter.FilteredData != undefined) {
                
                if ((!document.getElementById(FilterResp_ID).checked) && (!document.getElementById(FilterFFT_ID).checked)) {
                    // Show Filtered Data
                    // Update the trace of FilteredData
                    traces[1].x = ChannelList[ChNum].time,
                    traces[1].y = Mult(ChannelList[ChNum].Results.Filter.FilteredData, res.ScaleFactor);

                    // Update the trace of FilteredData
                    traces[2].x = [],
                    traces[2].y = [];
                }
                else if (document.getElementById(FilterResp_ID).checked) {
                    // Show Filter Response 
                    traces[0].x = [];
                    traces[0].y = [];

                    // Update the trace of RawData
                    traces[1].x = ChannelList[ChNum].Results.Filter.f;
                    traces[1].y = ChannelList[ChNum].Results.Filter.Mag;

                    // Update the trace of FilteredData
                    traces[2].x = ChannelList[ChNum].Results.Filter.f,
                    traces[2].y = ChannelList[ChNum].Results.Filter.Angle;
                }
                else if (document.getElementById(FilterFFT_ID).checked) {
                    // Show FFT of Filtered Data and FFT of Raw Data 
                    
                    // Calculate Fourier Amplitude Spectrum of FilteredData and RawData 
                    let Mag_FD=[], Angle_FD=[], Mag_RD=[], Angle_RD=[], f=[];
                    [Mag_FD, Angle_FD, f] = FourierSpec(Mult(ChannelList[ChNum].Results.Filter.FilteredData, res.ScaleFactor),                         ChannelList[ChNum].FSamp);
                    [Mag_RD, Angle_RD, f] = FourierSpec(Detrend(Mult(ChannelList[ChNum].data, (ChannelList[ChNum].ScaleFactor * res.ScaleFactor)), 0), ChannelList[ChNum].FSamp);

                    traces[0].x = f;
                    traces[0].y = Mag_RD;

                    traces[1].x = f;
                    traces[1].y = Mag_FD;

                    traces[2].x = [];
                    traces[2].y = [];

                }
                // Update the Statistics of FilteredData
                StatTableUpdate(ChannelList[ChNum], res);
            }
        }

        if ((!document.getElementById(FilterResp_ID).checked) && (!document.getElementById(FilterFFT_ID).checked)) {
            // Update trace infomation
            traces[0].visible     = true;
            traces[0].opacity     = 0.35;
            traces[0].line        = {color: 'grey', width: 1.00, dash: 'solid' };
            traces[0].name        = '<b>Raw Data<b>';
            traces[0].showlegend  = false,

            traces[1].visible     = true;
            traces[1].opacity     = 1.00;
            traces[1].line        = {color: 'blue', width: 1.50, dash: 'solid' };
            traces[1].name        = '<b>Filtered Data<b>';
            traces[1].showlegend  = false,

            traces[2].visible     = false;
            traces[2].opacity     = 0.35;
            traces[2].line        = {color: 'green', width: 1.00, dash: 'solid' };
            traces[2].name        = '<b>Phase<b>';
            traces[2].showlegend  = false,

            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";
            layout_update.yaxis.title.text      = res.yLabel;
        }
        else if (document.getElementById(FilterResp_ID).checked) {
            // Update trace infomation
            traces[0].visible     = false;
            traces[0].opacity     = 0.35;
            traces[0].line        = {color: 'grey', width: 1.00, dash: 'solid' };
            traces[0].name        = '<b>Raw Data<b>';
            traces[0].showlegend  = false,

            traces[1].visible     = true;
            traces[1].opacity     = 1.00;
            traces[1].line        = {color: 'blue', width: 1.50, dash: 'solid' };
            traces[1].name        = '<b>Magnitude<b>';
            traces[1].showlegend  = true,

            traces[2].visible     = true;
            traces[2].opacity     = 1.00;
            traces[2].line        = {color: 'green', width: 1.50, dash: 'solid' };
            traces[2].name        = '<b>Phase<b>';
            traces[2].showlegend  = true,

            layout_update.yaxis2.showticklabels = true;
            layout_update.yaxis2.title.text     = y2Label;
            layout_update.yaxis.title.text      = "<b>Magnitude<b>";
        }
        else if (document.getElementById(FilterFFT_ID).checked) {
            // Update trace infomation
            traces[0].visible     = true;
            traces[0].opacity     = 0.35;
            traces[0].line        = {color: 'grey', width: 1.00, dash: 'solid' };
            traces[0].name        = '<b>Raw Data<b>';
            traces[0].showlegend  = true,

            traces[1].visible     = true;
            traces[1].opacity     = 1.00;
            traces[1].line        = {color: 'blue', width: 1.50, dash: 'solid' };
            traces[1].name        = '<b>Filtered Data<b>';
            traces[1].showlegend  = true,

            traces[2].visible     = false;
            traces[2].opacity     = 0.35;
            traces[2].line        = {color: 'green', width: 1.00, dash: 'solid' };
            traces[2].name        = '<b>Phase<b>';
            traces[2].showlegend  = false,

            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";
            layout_update.yaxis.title.text      = res.yLabel_FFT;
        }
        
    }
    else if (PageNo == 2) {

        // Hide (do not show) the plot if it is not an acceleration record
        if (ChannelList[ChNum].Type == 0) { document.getElementById(PlotArea_ID).style.display = "flex"; }
        else                              { document.getElementById(PlotArea_ID).style.display = "none"; }

        // Update the Statistics of FilteredData as empty strings
        document.getElementById( Statictics_Peak_ID ).innerHTML = "";
        document.getElementById( Statictics_Mean_ID ).innerHTML = "";
        document.getElementById( Statictics_RMS_ID  ).innerHTML = "";
        document.getElementById(FilterType_ID).innerHTML        = "";
        document.getElementById(BaseLine_ID).innerHTML          = "";

        let Indx_Acc  = document.getElementById("Int_Acceleration").checked;
        let Indx_Vel  = document.getElementById("Int_Velocity").checked;
        let Indx_Disp = document.getElementById("Int_Displacement").checked;

        // Update the trace of FilteredData as empty array
        traces[0].x = [],
        traces[0].y = [];

        // Continue only if Channel.Results.Integral has the attributes of "Disp"
        // This mean this channel contains displacement data
        if ( ChannelList[ChNum].Results.Integral.hasOwnProperty("Disp") ) {

            // Continue if FilteredDate contains data
            if  (ChannelList[ChNum].Results.Integral.Disp != undefined) {

                if ((!document.getElementById(FilterResp_ID).checked) && (!document.getElementById(FilterFFT_ID).checked)) {
                    // Show Integrated Data
                    // Update the trace of RawData
                    traces[0].x = ChannelList[ChNum].time;
                    if (Indx_Acc)  { traces[0].y = Mult(ChannelList[ChNum].data,   (ChannelList[ChNum].ScaleFactor * res.ScaleFactor)); }
                    if (Indx_Vel)  { traces[0].y = Mult(ChannelList[ChNum].Results.Integral.Vel,  res.ScaleFactor);  }
                    if (Indx_Disp) { traces[0].y = Mult(ChannelList[ChNum].Results.Integral.Disp, res.ScaleFactor);  }

                    // Empty trace
                    traces[1].x = [],
                    traces[1].y = [];

                    // Empty trace
                    traces[2].x = [],
                    traces[2].y = [];
                }
                else if (document.getElementById(FilterResp_ID).checked) {
                    // Show Filter Response 
                    traces[0].x = [];
                    traces[0].y = [];

                    // Update the trace of Filter Response
                    traces[1].x = ChannelList[ChNum].Results.Integral.f;
                    traces[1].y = ChannelList[ChNum].Results.Integral.Mag;

                    // Update the trace of Filter Response
                    traces[2].x = ChannelList[ChNum].Results.Integral.f,
                    traces[2].y = ChannelList[ChNum].Results.Integral.Angle;
                }
                else if (document.getElementById(FilterFFT_ID).checked) {
                    // Show FFT of Integrated Data 
                    
                    // Calculate Fourier Amplitude Spectrum of FilteredData and RawData 
                    let Mag_FD=[], Angle_FD=[], f=[];

                    if (Indx_Acc)  { [Mag_FD, Angle_FD, f] = FourierSpec(Mult(ChannelList[ChNum].data,   (ChannelList[ChNum].ScaleFactor * res.ScaleFactor)), ChannelList[ChNum].FSamp);  }
                    if (Indx_Vel)  { [Mag_FD, Angle_FD, f] = FourierSpec(Mult(ChannelList[ChNum].Results.Integral.Vel,   res.ScaleFactor),                    ChannelList[ChNum].FSamp);  }
                    if (Indx_Disp) { [Mag_FD, Angle_FD, f] = FourierSpec(Mult(ChannelList[ChNum].Results.Integral.Disp,  res.ScaleFactor),                    ChannelList[ChNum].FSamp);  }

                    traces[0].x = f;
                    traces[0].y = Mag_FD;

                    traces[1].x = [];
                    traces[1].y = [];

                    traces[2].x = [];
                    traces[2].y = [];

                }

                // Update the Statistics of FilteredData
                StatTableUpdate(ChannelList[ChNum], res);
            }
        }

        if ((!document.getElementById(FilterResp_ID).checked) && (!document.getElementById(FilterFFT_ID).checked)) {
            // Update trace infomation
            traces[0].visible     = true;
            traces[0].opacity     = 1.00;
            traces[0].line        = {color: 'blue', width: 1.00, dash: 'solid' };
            traces[0].name        = '<b>Raw Data<b>';
            traces[0].showlegend  = false,

            traces[1].visible     = false;
            traces[1].opacity     = 0.35;
            traces[1].line        = {color: 'blue', width: 1.50, dash: 'solid' };
            traces[1].name        = '<b>Filtered Data<b>';
            traces[1].showlegend  = false,

            traces[2].visible     = false;
            traces[2].opacity     = 0.35;
            traces[2].line        = {color: 'green', width: 1.00, dash: 'solid' };
            traces[2].name        = '<b>Phase<b>';
            traces[1].showlegend  = false,

            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";
            if (Indx_Acc)    { yLabel  = res.yLabel       }
            if (Indx_Vel)    { yLabel =  res.yLabel_vel   }
            if (Indx_Disp)   { yLabel =  res.yLabel_disp  }
            layout_update.yaxis.title.text = yLabel;
        }
        else if (document.getElementById(FilterResp_ID).checked) {
            // Update trace infomation
            traces[0].visible     = false;
            traces[0].opacity     = 0.35;
            traces[0].line        = {color: 'grey', width: 1.00, dash: 'solid' };
            traces[0].name        = '<b>Raw Data<b>';
            traces[0].showlegend  = false,

            traces[1].visible     = true;
            traces[1].opacity     = 1.00;
            traces[1].line        = {color: 'blue', width: 1.50, dash: 'solid' };
            traces[1].name        = '<b>Magnitude<b>';
            traces[1].showlegend  = true,

            traces[2].visible     = true;
            traces[2].opacity     = 1.00;
            traces[2].line        = {color: 'green', width: 1.50, dash: 'solid' };
            traces[2].name        = '<b>Phase<b>';
            traces[2].showlegend  = true,

            layout_update.yaxis2.showticklabels = true;
            layout_update.yaxis2.title.text     = y2Label;
            layout_update.yaxis.title.text      = "<b>Magnitude<b>";
        }
        else if (document.getElementById(FilterFFT_ID).checked) {
            // Update trace infomation
            traces[0].visible     = true;
            traces[0].opacity     = 1.00;
            traces[0].line        = {color: 'blue', width: 1.00, dash: 'solid' };
            traces[0].name        = '<b>Raw Data<b>';
            traces[0].showlegend  = false,

            traces[1].visible     = false;
            traces[1].opacity     = 1.00;
            traces[1].line        = {color: 'blue', width: 1.50, dash: 'solid' };
            traces[1].name        = '<b>Filtered Data<b>';
            traces[1].showlegend  = false,

            traces[2].visible     = false;
            traces[2].opacity     = 0.35;
            traces[2].line        = {color: 'green', width: 1.00, dash: 'solid' };
            traces[2].name        = '<b>Phase<b>';
            traces[2].showlegend  = false,

            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";
            layout_update.yaxis.title.text      = res.yLabel_FFT;
        }
    }
    
    // Update the graph in the DIV.
    Plotly.update(Div_ID, traces, layout_update);

    


    // Updates all values in the Statictics Table
    function StatTableUpdate(Channel, res) {
        // Decleration of Variables
        let Unit_List, Indx_Acc, Indx_Vel, Indx_Disp;

        // Update Statistics
        if (PageNo == 0) {
            // Raw Data Statistics
            document.getElementById( Statictics_Peak_ID ).innerHTML = (Channel.Peak * res.ScaleFactor).toPrecision(4);
            document.getElementById( Statictics_Mean_ID ).innerHTML = (Channel.Mean * res.ScaleFactor).toPrecision(4);
            document.getElementById( Statictics_RMS_ID  ).innerHTML = (Channel.RMS  * res.ScaleFactor).toPrecision(4);

            // Get the list of Units for this channel
            Unit_List = UnitList(Channel.Unit);
        }
        else if (PageNo == 1) {
            // Filtered Data Statistics
            document.getElementById( Statictics_Peak_ID ).innerHTML = (Channel.Results.Filter.Peak * res.ScaleFactor).toPrecision(4);
            document.getElementById( Statictics_Mean_ID ).innerHTML = (Channel.Results.Filter.Mean * res.ScaleFactor).toPrecision(4);
            document.getElementById( Statictics_RMS_ID  ).innerHTML = (Channel.Results.Filter.RMS  * res.ScaleFactor).toPrecision(4);

            // Update Filter Info on table 
            let FilterInfo = ChannelList[ChNum].Results.Filter.FilterName_String;
            if (FilterInfo.toUpperCase() == ("None").toUpperCase()) {
                document.getElementById(FilterType_ID).innerHTML    = FilterInfo;
                document.getElementById(FilterRow_ID).style.display = "table-row";
                document.getElementById(FilterType_ID).className    = "Plotly_Stat_Body_Td_Right";
            }
            else {
                if (ChannelList[ChNum].Results.Filter.ErrorMessage != undefined) {
                    document.getElementById(FilterType_ID).innerHTML    =  ChannelList[ChNum].Results.Filter.ErrorMessage;
                    document.getElementById(FilterRow_ID).style.display = "table-row";
                    document.getElementById(FilterType_ID).className    = "Plotly_Stat_Body_Td_Right_Red";
                }
                else {
                    FilterInfo += " " + ChannelList[ChNum].Results.Filter.FilterType_String;
                    FilterInfo += " " + ChannelList[ChNum].Results.Filter.FilterBand;
                    document.getElementById(FilterType_ID).innerHTML    = FilterInfo;
                    document.getElementById(FilterRow_ID).style.display = "table-row";
                    document.getElementById(FilterType_ID).className    = "Plotly_Stat_Body_Td_Right";
                }
            }

            // Update Baseline Correction Info on table
            let BaseLineInfo = ChannelList[ChNum].Results.Filter.BaselineCorrection_String; 
            document.getElementById(BaseLine_ID).innerHTML        = BaseLineInfo;
            document.getElementById(BaseLineRow_ID).style.display = "table-row";

            // Get the list of Units for this channel
            Unit_List = UnitList(Channel.Unit);
        }
        else if (PageNo == 2) {
            // Integral Statistics

            Indx_Acc  = document.getElementById("Int_Acceleration").checked;
            Indx_Vel  = document.getElementById("Int_Velocity").checked;
            Indx_Disp = document.getElementById("Int_Displacement").checked;
            
            if (Indx_Acc) {
                document.getElementById( Statictics_Peak_ID ).innerHTML = (Channel.Peak * res.ScaleFactor).toPrecision(4);
                document.getElementById( Statictics_Mean_ID ).innerHTML = (Channel.Mean * res.ScaleFactor).toPrecision(4);
                document.getElementById( Statictics_RMS_ID  ).innerHTML = (Channel.RMS  * res.ScaleFactor).toPrecision(4);

                // Get the list of Units for this channel
                Unit_List = UnitList(Channel.Unit);

            }
            else if (Indx_Vel) {
                document.getElementById( Statictics_Peak_ID ).innerHTML = (Channel.Results.Integral.Peak_Vel * res.ScaleFactor).toPrecision(4);
                document.getElementById( Statictics_Mean_ID ).innerHTML = (Channel.Results.Integral.Mean_Vel * res.ScaleFactor).toPrecision(4);
                document.getElementById( Statictics_RMS_ID  ).innerHTML = (Channel.Results.Integral.RMS_Vel  * res.ScaleFactor).toPrecision(4);

                // Get the list of Units for this channel
                Unit_List = UnitList(4);
            }
            else if (Indx_Disp) {
                document.getElementById( Statictics_Peak_ID ).innerHTML = (Channel.Results.Integral.Peak_Disp * res.ScaleFactor).toPrecision(4);
                document.getElementById( Statictics_Mean_ID ).innerHTML = (Channel.Results.Integral.Mean_Disp * res.ScaleFactor).toPrecision(4);
                document.getElementById( Statictics_RMS_ID  ).innerHTML = (Channel.Results.Integral.RMS_Disp  * res.ScaleFactor).toPrecision(4);

                // Get the list of Units for this channel
                Unit_List = UnitList(8);
            }

            // Update Integral-Info on table 
            let IntegralInfo = Channel.Results.Integral.FilterName_String;
            if (IntegralInfo.toUpperCase() == ("None").toUpperCase()) {
                document.getElementById(FilterType_ID).innerHTML    = IntegralInfo;
                document.getElementById(FilterRow_ID).style.display = "table-row";
                document.getElementById(FilterType_ID).className    = "Plotly_Stat_Body_Td_Right";
                document.getElementById(GraphUnitRow_ID).style.display = "table-row";
            }
            else {
                if (Channel.Results.Integral.ErrorMessage != undefined) {
                    document.getElementById(FilterType_ID).innerHTML    =  Channel.Results.Integral.ErrorMessage;
                    document.getElementById(FilterRow_ID).style.display = "table-row";
                    document.getElementById(FilterType_ID).className    = "Plotly_Stat_Body_Td_Right_Red";
                    document.getElementById(GraphUnitRow_ID).style.display = "table-row";
                }
                else {
                    IntegralInfo += " " + Channel.Results.Integral.FilterType_String;
                    IntegralInfo += " " + Channel.Results.Integral.FilterBand;
                    document.getElementById(FilterType_ID).innerHTML    = IntegralInfo;
                    document.getElementById(FilterRow_ID).style.display = "table-row";
                    document.getElementById(FilterType_ID).className    = "Plotly_Stat_Body_Td_Right";
                    document.getElementById(GraphUnitRow_ID).style.display = "table-row";
                }
            }

            // Update Baseline Correction Info on table
            let BaseLineInfo = Channel.Results.Integral.BaselineCorrection_String;
            document.getElementById(BaseLine_ID).innerHTML        = BaseLineInfo;
            document.getElementById(BaseLineRow_ID).style.display = "table-row";

        }

        // Update Graph-Unit 
        // Empties the content of cell-element in table
        document.getElementById(Unit_Cell_ID).innerHTML = "";
        
        // Create the select element using (Unit_List.Units)
        let select = Create_Select_Element(Unit_List.Units, Unit_Plot_ID);
        select.setAttribute("onchange", "Graph_unit_Change(" + Channel.ChNum + ")" );
        if      (PageNo == 1) { select.selectedIndex = Channel.Results.Filter.Unit_Ind;   }
        else if (PageNo == 2) { select.selectedIndex = Channel.Results.Integral.Unit_Ind; }
        
        // Asign select-element to cell-element in table
        document.getElementById(Unit_Cell_ID).appendChild(select);

    }

    // Creates new Select-Element for Graph-Unit List
    function Create_Select_Element(Unit_List, ID) {
        
        // Decleration of variables
        let j, opt, select;

        // Create select element and populate it 
        select  = document.createElement('select');
        select.setAttribute("id", ID);
        select.setAttribute('class', 'form-select form-select-sm');
        
        // All options for the select element 
        for (j = 0; j < Unit_List.length; j++) {
            opt = document.createElement("option");
            opt.value = Unit_List[j];
            opt.text = Unit_List[j];
            select.add(opt, null);
        }
        return select;
    }

}

function Graph_unit_Change(ChNum) {

    let Unit_Plot_ID = "Unit_Plot_ID_" + ChannelList[ChNum].Unique_ID;
    let Indx         = document.getElementById(Unit_Plot_ID).selectedIndex;

    if      (PageNo == 0) {  }
    else if (PageNo == 1) { ChannelList[ChNum].Results.Filter.Unit_Ind   = Indx;   Plot_Update(ChannelList[ChNum].ChNum); }
    else if (PageNo == 2) { ChannelList[ChNum].Results.Integral.Unit_Ind = Indx;   Plot_Update(ChannelList[ChNum].ChNum); }
}

function FilterResponse(el) {
    // Decleartion of variables 
    let ChNum, ID, ID_2;

    if (el.id.includes("FilterResp_ID_")) { 

        ID     = el.id.replace("FilterResp_ID_", ""); 
        ChNum  = ChannelList_UniqueID(ID);
        ID_2   = "FilterFFT_ID_" + ChannelList[ChNum].Unique_ID;  // uncheck FFT
        document.getElementById(ID_2).checked = false; 

    }
    else if (el.id.includes("FilterFFT_ID_")) { 

        ID     = el.id.replace("FilterFFT_ID_",  ""); 
        ChNum  = ChannelList_UniqueID(ID);
        ID_2   = "FilterResp_ID_" + ChannelList[ChNum].Unique_ID;  // uncheck Filter Response
        document.getElementById(ID_2).checked = false; 

    }
    
    Plot_Update(ChNum);
}

function Plot_Create_Old(Div_ID, NumOfTrace, NumOfTrace_On_Y2_Axis) {

    // Declaration of variables
    let i, traces=[], layout, config, Leg_Name, temp, data

    // Check the default agrument 
    if (NumOfTrace == null)   { NumOfTrace = 1; }

    data = Array(1).fill(0);

    // Update legend
    Leg_Name = '<b> Trace <b>';

    // Define Traces (Data)
    for (i=0; i<NumOfTrace; i++) {

        temp = "y1";

        // Some trace might be at the right side of the gragh
        if ((NumOfTrace_On_Y2_Axis != null) && (i>=NumOfTrace_On_Y2_Axis) ) { 
           temp = 'y2';
        }

        traces.push({ 
            x             : [0,1,2,3,4,5,6,7,8,9],
            y             : [0,-5,4,-2,-6,8,9,5,-1,5],
            mode          : 'lines',
            type          : 'scatter',
            yaxis         : temp,
            name          : Leg_Name,
            line          : {color: 'blue', width: 0.75, dash: 'solid' },
            showlegend    : false,
        })
    };

    // Define Layout
    layout = {
        xaxis           : { zeroline: false, automargin: true, tickfont: { size: 15 }, linecolor: 'black', linewidth: 1, mirror: true, title: {text: '<b>Time [s]<b>'   , standoff: 5, font: {family: "Arial", size: 15} }, autorange: true },
        yaxis           : { zeroline: true,  automargin: true, tickfont: { size: 15 }, linecolor: 'black', linewidth: 1, mirror: true, title: {text: '<b> Amplitude <b>', standoff: 5, font: {family: "Arial", size: 15} }, autorange: true },
        plot_bgcolor    : 'rgb(255,255,255)',
        paper_bgcolor   : 'rgb(255,255,255)',
        legend          : { x: 1, y:1, xanchor: 'right', orientation: 'v', font: {size: 16}},
        autosize        : true,
        margin          : {t: 20, r:20, b:5, l:5},
        shapes          : [],
    };

    if (NumOfTrace_On_Y2_Axis != null) { 
        layout.yaxis2 = { zeroline: false, overlaying: 'y', side: 'right', automargin: true, tickfont: { size: 15 }, linecolor: 'black', linewidth: 1, mirror: true, title: {text: '<b> Amplitude <b>', standoff: 20, font: {family: "Arial", size: 15} }};
    }

    config = {
        responsive              : true,
        displayModeBar          : true,    // show-hide floating toolbar all together 
        modeBarButtonsToRemove  : [],
        displaylogo             : false,   // Romoves the Ployly logo from toolbar
        useResizeHandler        : true,    // Enables Plotly's resize event listener
        showTips                : false,
        scrollZoom              : false,   // Enable mouse wheel zooming
    }
    // Display using Plotly
    Plotly.newPlot(Div_ID, traces, layout, config);
}
