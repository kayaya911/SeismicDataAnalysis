// mathJS library is with lowercase math.multiply()

// stricter parsing and error handling
"use strict";

//-------------------------------------------------------------------------------------------------------------
async function Plotly_Graph_Update(ChNum) {
    
    // Return if no graph is created for this channel
    if (!ChannelList[ChNum].PlotGraph) { return; } 
    else { document.getElementById("Div_ID_"+ChannelList[ChNum].Unique_ID).style.display = 'flex'; }

    // Decleration of variables 
    let traces, layout_update, res, res_RawData, timeData, res_FilteredData, FilterInfo;
    let Indx, Indx_Acc, Indx_Vel, Indx_Disp, yTitle;
    let IsFilter_CheckBox_Selected, IsFFT_CheckBox_Selected, DisplayData;
    let plotCount, temp, temp1, Table, T1, T2;
    let res1=[], res2=[], res3=[], tD, LS1, LS2, LS3, LS4, LT1='', LT2='', LT3='', LT4='', tD1, tD2, tD3, tD4;

    let PlotArea_ID            = "PlotArea_ID_"        + ChannelList[ChNum].Unique_ID;
    let Statictics_Peak_ID     = "Statictics_Peak_ID_" + ChannelList[ChNum].Unique_ID;
    let Statictics_Mean_ID     = "Statictics_Mean_ID_" + ChannelList[ChNum].Unique_ID;
    let Statictics_RMS_ID      = "Statictics_RMS_ID_"  + ChannelList[ChNum].Unique_ID;
    let Statictics_Residual_ID = "Statictics_Residual_ID_"  + ChannelList[ChNum].Unique_ID;
    let GraphUnitRow_ID        = "GraphUnitRow_ID_"    + ChannelList[ChNum].Unique_ID;
    let Unit_Cell_ID           = "Unit_Cell_ID_"       + ChannelList[ChNum].Unique_ID;
    let BaseLine_ID            = "BaseLine_ID_"        + ChannelList[ChNum].Unique_ID;
    let BaseLineRow_ID         = "BaseLineRow_ID_"     + ChannelList[ChNum].Unique_ID;
    let SDOF_Row_ID            = "SDOF_Row_ID_"        + ChannelList[ChNum].Unique_ID;
    let SDOF_Plot_ID           = "SDOF_Plot_ID_"       + ChannelList[ChNum].Unique_ID;
    let SDOF_Method_ID         = "SDOF_Method_ID_"     + ChannelList[ChNum].Unique_ID;
    let SDOF_Method_Row_ID     = "SDOF_Method_Row_ID_" + ChannelList[ChNum].Unique_ID;
    let FilterType_ID          = "FilterType_ID_"      + ChannelList[ChNum].Unique_ID;
    let FilterRow_ID           = "FilterRow_ID_"       + ChannelList[ChNum].Unique_ID;
    let FilterResp_ID          = "FilterResp_ID_"      + ChannelList[ChNum].Unique_ID;
    let FilterFFT_ID           = "FilterFFT_ID_"       + ChannelList[ChNum].Unique_ID;
    let Unit_Plot_ID           = "Unit_Plot_ID_"       + ChannelList[ChNum].Unique_ID;

    // Get the existing traces in the PlotArea_ID
    traces = document.getElementById(PlotArea_ID).data;

    // Get the existing layout in the PlotArea_ID
    layout_update = document.getElementById(PlotArea_ID).layout;

    // Clear the graph by default then populate it later.
    Plotly_Clear_Graph(ChNum);

    // Update Graph for this channel
    if (PageNo == 0) {
        // Rawdata page

        // Scale the Rawdata to the user-specified unit in Plotly Graph (info table)
        res = Convert_Data_To_Graph_Unit(ChannelList[ChNum].data, ChNum);

        // Update trace for rawData
        traces[0].x             = ChannelList[ChNum].time;
        traces[0].y             = Multiply(res.Data, ChannelList[ChNum].ScaleFactor);
        traces[0].mode          = 'lines',
        traces[0].marker        = { color: 'blue', size: 5, symbol: 'circle' },
        traces[0].yaxis         = "y1",
        traces[0].visible       = true;
        traces[0].opacity       = 1.00;
        traces[0].line          = {color: 'blue', width: 1.50, dash: 'solid' };
        traces[0].name          = '<b>Raw Data<b>';  // legend title
        traces[0].showlegend    = false;             // Don't show legend

        // Update y-axis of Plotly Graph
        layout_update.yaxis.title.text      = res.yTitle;
        layout_update.yaxis2.showticklabels = false;
        layout_update.yaxis2.title.text     = "";

        // Update the Statistics of RawData in table - scaled to user-specified unit
        document.getElementById( Statictics_Peak_ID ).innerHTML = res.Peak.toPrecision(4);
        document.getElementById( Statictics_Mean_ID ).innerHTML = res.Mean.toPrecision(4);
        document.getElementById( Statictics_RMS_ID  ).innerHTML = res.RMS.toPrecision(4);

        // Asign select-element to cell-element in table
        document.getElementById(Unit_Cell_ID).innerHTML = "";
        document.getElementById(Unit_Cell_ID).appendChild(res.Or_Units);

    } 
    else if (PageNo == 1) {
        // Filtering Page 

        // Scale the Rawdata to the user-specified unit on Plotly Graph (info table)
        res_RawData = Convert_Data_To_Graph_Unit(ChannelList[ChNum].data, ChNum);

        // Make sure that the filter analysis is successfully completed
        if (ChannelList[ChNum].Results.Filter.IsAnalysisCompleted) {

            // Scale the FilteredData to the user-specified unit on Plotly Graph (info table)
            res_FilteredData = Convert_Data_To_Graph_Unit(ChannelList[ChNum].Results.Filter.FilteredData, ChNum); 

            // Get the status of two checkboxes in the Info table on Plotly Graph
            IsFilter_CheckBox_Selected = document.getElementById(FilterResp_ID).checked;
            IsFFT_CheckBox_Selected    = document.getElementById(FilterFFT_ID).checked;
            
            if (!IsFilter_CheckBox_Selected && !IsFFT_CheckBox_Selected) {
                
                // Plot Raw data in trace[0]
                traces[0].x           = ChannelList[ChNum].time;
                traces[0].y           = Detrend(Multiply(res_RawData.Data, ChannelList[ChNum].ScaleFactor), 0); // Remove mean
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'grey', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 0.35;
                traces[0].line        = {color: 'grey', width: 1.00, dash: 'solid' };
                traces[0].name        = '<b>Raw Data<b>';   // legend title
                traces[0].showlegend  = true;               // Show legend

                // Plot the Filtered data in tarace[1]
                traces[1].x           = ChannelList[ChNum].time,
                traces[1].y           = res_FilteredData.Data;
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y1",
                traces[1].visible     = true;
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[1].name        = '<b>Filtered Data<b>';  // legend title
                traces[1].showlegend  = true;                   // Show legend 

                layout_update.yaxis.title.text      = res_FilteredData.yTitle;   // This is the unit that user wants to see on the graph.
                layout_update.yaxis2.showticklabels = false;
                layout_update.yaxis2.title.text     = "";
            } 
            else if (IsFilter_CheckBox_Selected) {
                
                // Plot Filter Magnitude Response in trace[1]
                traces[0].x           = ChannelList[ChNum].Results.Filter.f;
                traces[0].y           = ChannelList[ChNum].Results.Filter.Mag;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = '<b>Magnitude<b>';      // legend title
                traces[0].showlegend  = true;                   // Show legend 

                // Plot Filter Phase Response in trace[2]
                traces[1].x           = ChannelList[ChNum].Results.Filter.f;
                traces[1].y           = ChannelList[ChNum].Results.Filter.Angle;
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'green', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y2",
                traces[1].visible     = true;
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'green', width: 1.50, dash: 'solid' };
                traces[1].name        = '<b>Phase<b>';          // legend title
                traces[1].showlegend  = true;                   // Show legend 

                layout_update.yaxis.title.text      = res_FilteredData.yTitle_FFT;   
                layout_update.yaxis2.showticklabels = true; 
                layout_update.yaxis2.title.text     = res_FilteredData.y2Title;
                
            }
            else if (IsFFT_CheckBox_Selected) {
                // Show FFT of Filtered Data and FFT of Raw Data 
                    
                // Declaration of variables 
                let Mag_FD=[], Angle_FD=[], Mag_RD=[], Angle_RD=[], f=[];
                
                // FFT of Filtered data
                [Mag_FD, Angle_FD, f] = FourierSpec(res_FilteredData.Data,  ChannelList[ChNum].FSamp);
                
                // FFT of Raw data
                [Mag_RD, Angle_RD, f] = FourierSpec(Detrend(Multiply(res_RawData.Data, ChannelList[ChNum].ScaleFactor), 0),   ChannelList[ChNum].FSamp);

                // Plot FFT Magnitude Response of Raw data
                traces[0].x           = f;
                traces[0].y           = Mag_RD;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'grey', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 0.35;
                traces[0].line        = {color: 'grey', width: 1.00, dash: 'solid' };
                traces[0].name        = '<b>Raw Data<b>';       // legend title
                traces[0].showlegend  = true;                   // Show legend 

                // Plot FFT Magnitude Response of Filtered data
                traces[1].x           = f;
                traces[1].y           = Mag_FD;
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y1",
                traces[1].visible     = true;
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[1].name        = '<b>Filtered Data<b>';      // legend title
                traces[1].showlegend  = true;                       // Show legend 

                layout_update.yaxis.title.text      = res_FilteredData.yTitle_FFT; 
                layout_update.yaxis2.showticklabels = false;
                layout_update.yaxis2.title.text     = "";

            }

            // Update the Statistics of RawData in table - scaled to user-specified unit
            document.getElementById( Statictics_Peak_ID ).innerHTML = res_FilteredData.Peak.toPrecision(4);
            document.getElementById( Statictics_Mean_ID ).innerHTML = res_FilteredData.Mean.toPrecision(4);
            document.getElementById( Statictics_RMS_ID  ).innerHTML = res_FilteredData.RMS.toPrecision(4);

            // Show Baseline-Row in Info table
            document.getElementById(BaseLine_ID).innerHTML = ChannelList[ChNum].Results.Filter.BaselineCorrection_String;

            // Show Filter_ID-Row in InfoBar
            FilterInfo  = ChannelList[ChNum].Results.Filter.FilterName_String;
            FilterInfo += "<br>" + ChannelList[ChNum].Results.Filter.FilterType_String;
            FilterInfo += " " + ChannelList[ChNum].Results.Filter.FilterBand;
            FilterInfo += "<br> Zero Phase: " + ChannelList[ChNum].Results.Filter.ZeroPhase;
            document.getElementById(FilterType_ID).innerHTML = FilterInfo;

        } 
        else {
            // 
            traces[0].x           = ChannelList[ChNum].time;
            traces[0].y           = Multiply(res_RawData.Data, ChannelList[ChNum].ScaleFactor);
            traces[0].mode        = 'lines',
            traces[0].marker      = { color: 'grey', size: 5, symbol: 'circle' },
            traces[0].visible     = true;
            traces[0].opacity     = 0.35;
            traces[0].line        = {color: 'grey', width: 1.00, dash: 'solid' };
            traces[0].name        = '';                 // legend title
            traces[0].showlegend  = false;              // don't show legend 

            layout_update.yaxis.title.text      = res_RawData.yTitle;   // This is the unit that user wants to see on the graph.
            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";

            // Empty the Statistics
            document.getElementById( Statictics_Peak_ID ).innerHTML = '';
            document.getElementById( Statictics_Mean_ID ).innerHTML = '';
            document.getElementById( Statictics_RMS_ID  ).innerHTML = '';

            // Show Baseline-Row in Info table
            document.getElementById(BaseLine_ID).innerHTML = '';

            // Show Filter_ID-Row in InfoBar
            document.getElementById(FilterType_ID).innerHTML = '';

        }

        // Asign select-element to cell-element in table
        document.getElementById(Unit_Cell_ID).innerHTML = "";
        document.getElementById(Unit_Cell_ID).appendChild(res_RawData.Or_Units);

    }
    else if (PageNo == 2) {
        // Integration Page

        // Find out which results to plot (Acc, Vel or Disp)
        Indx_Acc  = document.getElementById("Int_Acceleration").checked;
        Indx_Vel  = document.getElementById("Int_Velocity").checked;
        Indx_Disp = document.getElementById("Int_Displacement").checked;

        // Make sure that the integration analysis is successfully completed
        if (ChannelList[ChNum].Results.Integral.IsAnalysisCompleted) {

            // Scale data (Acc, Vel or Disp) to the user-specified unit on Plotly Graph (info table)
            if      (Indx_Acc)  { res = Convert_Data_To_Graph_Unit(Multiply(ChannelList[ChNum].data, ChannelList[ChNum].ScaleFactor),  ChNum);   yTitle = res.yTitle;       }
            else if (Indx_Vel)  { res = Convert_Data_To_Graph_Unit(ChannelList[ChNum].Results.Integral.Vel,                            ChNum);   yTitle = res.yTitle_Vel;   }
            else if (Indx_Disp) { res = Convert_Data_To_Graph_Unit(ChannelList[ChNum].Results.Integral.Disp,                           ChNum);   yTitle = res.yTitle_Disp;  }

            // Get checkbox status
            IsFilter_CheckBox_Selected = document.getElementById(FilterResp_ID).checked;
            IsFFT_CheckBox_Selected    = document.getElementById(FilterFFT_ID).checked;
            
            if (!IsFilter_CheckBox_Selected && !IsFFT_CheckBox_Selected) {

                traces[0].x           = ChannelList[ChNum].time;
                traces[0].y           = res.Data;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = '';       // legend title
                traces[0].showlegend  = false;    // Don't show legend

                layout_update.yaxis.title.text      = yTitle;   // This is the unit that user wants to see on the graph.
                layout_update.yaxis2.showticklabels = false;
                layout_update.yaxis2.title.text     = "";
            }
            else if (IsFilter_CheckBox_Selected) {
                
                // Plot Filter Magnitude Response in trace[1]
                traces[0].x           = ChannelList[ChNum].Results.Integral.f;
                traces[0].y           = ChannelList[ChNum].Results.Integral.Mag;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = '<b>Magnitude<b>';      // legend title 
                traces[0].showlegend  = true;                   // Show legend

                // Plot Filter Phase Response in trace[2]
                traces[1].x           = ChannelList[ChNum].Results.Integral.f;
                traces[1].y           = ChannelList[ChNum].Results.Integral.Angle;
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'green', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y2",
                traces[1].visible     = true;
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'green', width: 1.50, dash: 'solid' };
                traces[1].name        = '<b>Phase<b>';          // legend title 
                traces[1].showlegend  = true;                   // Show legend

                layout_update.yaxis.title.text      = res.yTitle_FFT;
                layout_update.yaxis2.showticklabels = true;
                layout_update.yaxis2.title.text     = res.y2Title;
            }
            else if (IsFFT_CheckBox_Selected) {
                // Show FFT of Integrated Data
                    
                // Declaration of variables 
                let Mag_FD=[], Angle_FD=[], f=[];
                
                // FFT of Filtered data
                [Mag_FD, Angle_FD, f] = FourierSpec(res.Data,  ChannelList[ChNum].FSamp);
                
                // Empty trace[1] - Not used 
                traces[0].x           = f;
                traces[0].y           = Mag_FD;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = '';                     // legend title
                traces[0].showlegend  = false;                  // Show legend

                layout_update.yaxis.title.text      = res.yTitle_FFT;
                layout_update.yaxis2.showticklabels = false;
                layout_update.yaxis2.title.text     = "";
            }

            // Update the Statistics of (RawData, Velocity, Displacement) in table - scaled to user-specified unit
            document.getElementById( Statictics_Peak_ID     ).innerHTML = res.Peak.toPrecision(4);
            document.getElementById( Statictics_Mean_ID     ).innerHTML = res.Mean.toPrecision(4);
            document.getElementById( Statictics_RMS_ID      ).innerHTML = res.RMS.toPrecision(4);
            document.getElementById( Statictics_Residual_ID ).innerHTML = res.Residual.toPrecision(4);

            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = ChannelList[ChNum].Results.Integral.BaselineCorrection_String;

            // Show Filter_ID-Row in InfoBar
            FilterInfo = ChannelList[ChNum].Results.Integral.FilterName_String;
            FilterInfo += "<br>" + ChannelList[ChNum].Results.Integral.FilterType_String;
            FilterInfo += " " + ChannelList[ChNum].Results.Integral.FilterBand;
            FilterInfo += "<br> Zero Phase: " + ChannelList[ChNum].Results.Integral.ZeroPhase;
            document.getElementById(FilterType_ID).innerHTML    = FilterInfo;

        } else {

            res      = Convert_Data_To_Graph_Unit(ChannelList[ChNum].data,  ChNum);    
            res.Data = Multiply(res.Data, ChannelList[ChNum].ScaleFactor);

            // 
            traces[0].x = ChannelList[ChNum].time;
            if      (Indx_Acc)  { traces[0].y = res.Data;   yTitle = res.yTitle;       }
            else if (Indx_Vel)  { traces[0].y = [];         yTitle = res.yTitle_Vel;   }
            else if (Indx_Disp) { traces[0].y = [];         yTitle = res.yTitle_Disp;  }

            traces[0].mode        = 'lines',
            traces[0].marker      = { color: 'grey', size: 5, symbol: 'circle' },
            traces[0].visible     = true;
            traces[0].opacity     = 0.35;
            traces[0].line        = {color: 'grey', width: 1.00, dash: 'solid' };
            traces[0].name        = '';         // legend title
            traces[0].showlegend  = false;      // Don't show legend    

            layout_update.yaxis.title.text      = yTitle;
            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";

            // Empty the Statistics
            document.getElementById( Statictics_Peak_ID      ).innerHTML = '';
            document.getElementById( Statictics_Mean_ID      ).innerHTML = '';
            document.getElementById( Statictics_RMS_ID       ).innerHTML = '';
            document.getElementById( Statictics_Residual_ID  ).innerHTML = '';

            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = '';

            // Show Filter_ID-Row in InfoBar
            document.getElementById(FilterType_ID).innerHTML    = '';

        }

        // Asign select-element to cell-element in table
        document.getElementById(Unit_Cell_ID).innerHTML = "";
        if      (Indx_Acc)  { document.getElementById(Unit_Cell_ID).appendChild(res.Or_Units);     }
        else if (Indx_Vel)  { document.getElementById(Unit_Cell_ID).appendChild(res.Vel_Select);   }
        else if (Indx_Disp) { document.getElementById(Unit_Cell_ID).appendChild(res.Disp_Select);  }

    }
    else if (PageNo == 3) {
        // SDOF Page

        // Get the index number of the SDOF_Analysis Method
        Indx = document.getElementById('SDOF_Analysis').selectedIndex;

        // Make sure that the filter analysis is successfully completed
        if (ChannelList[ChNum].Results.SDOF.IsAnalysisCompleted && (ChannelList[ChNum].Results.SDOF.AnalysisMethod == Indx)) {

            // Scale the data to the user-specified unit in Plotly Graph (info table)
            DisplayData  = ChannelList[ChNum].Results.SDOF.DisplayData;
            timeData     = ChannelList[ChNum].time;
            res1         = []; 
            res2         = [];
            res3         = [];
            tD1          = timeData;
            tD2          = [];
            tD3          = [];
            tD4          = [];
            LS1          = false;
            LS2          = false;
            LS3          = false;
            LS4          = false;

            if      (DisplayData == "Disp"    ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Disp,     ChNum );  
                                                if (ChannelList[ChNum].Results.SDOF.AnalysisMethod == 1) {
                                                    res1  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Up,       ChNum ); 
                                                    res2  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Uc,       ChNum );
                                                    res3  = [];
                                                    tD1  = timeData; tD2  = timeData; tD3  = timeData; tD4  = []; 
                                                    LS1   = true; LS2= true; LS3 = true;
                                                    LT1  = 'Displacement'; 
                                                    LT2  = 'Steady State'; 
                                                    LT3  = 'Transient';
                                                }
            }
            else if (DisplayData == "Vel"     ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Vel,      ChNum );  }
            else if (DisplayData == "acc"     ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.acc,      ChNum );  }
            else if (DisplayData == "Acc"     ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Acc,      ChNum );  }
            else if (DisplayData == "Fs"      ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Fs,       ChNum );  }
            else if (DisplayData == "Fc"      ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Fc,       ChNum );  }
            else if (DisplayData == "Fi"      ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Fi,       ChNum );  }
            else if (DisplayData == "HarFor"  ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.HarForce, ChNum );  }
            else if (DisplayData == 'E'       ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Ek,       ChNum );  
                                                  res1 = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Ed,       ChNum );
                                                  res2 = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Es,       ChNum );
                                                  res3 = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Ei,       ChNum );
                                                  tD1  = timeData; tD2  = timeData; tD3  = timeData; tD4  = timeData; 
                                                  LS1  = true;  LS2  = true;  LS3  = true;  LS4  = true;
                                                  LT1  = 'Kinetic Energy'; 
                                                  LT2  = 'Damping Energy'; 
                                                  LT3  = 'Strain Energy'; 
                                                  LT4  = 'Input Energy'; 
                                                }
            
            else if (DisplayData == "Hyst"    ) { res  = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Fs,       ChNum );
                                                  res1 = Convert_Data_To_Graph_Unit_SDOF(ChannelList[ChNum].Results.SDOF.Disp,     ChNum );
                                                  tD1  = res1.Data;

                                                  let xMin = Min(res1.Data).val;
                                                  let xMax = Max(res1.Data).val;
                                                  let yMin = Min(res.Data).val;
                                                  let yMax = Max(res.Data).val;
                                                  let cc   = 0.1;
                                                  layout_update.xaxis.autorange = false;
                                                  layout_update.yaxis.autorange = false;
                                                  layout_update.xaxis.range     = [xMin*(1+cc), xMax*(1+cc)];
                                                  layout_update.yaxis.range     = [yMin*(1+cc), yMax*(1+cc)];
                                                }

            // Get the status of two checkboxes in the Infor table on Plotly Graph
            IsFFT_CheckBox_Selected    = document.getElementById(FilterFFT_ID).checked;

            if (!IsFilter_CheckBox_Selected && !IsFFT_CheckBox_Selected) {
                
                // Plot Raw data in trace[0]
                traces[0].x           = tD1;
                traces[0].y           = res.Data;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = LT1;                // legend title
                traces[0].showlegend  = LS1;       // Show legend

                traces[1].x           = tD2;
                traces[1].y           = res1.Data;
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'red', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y1",
                traces[1].visible     = 'legendonly';
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'red', width: 1.50, dash: 'solid' };
                traces[1].name        = LT2;                // legend title
                traces[1].showlegend  = LS2;       // Show legend


                traces[2].x           = tD3;
                traces[2].y           = res2.Data;
                traces[2].mode        = 'lines',
                traces[2].marker      = { color: 'green', size: 5, symbol: 'circle' },
                traces[2].yaxis       = "y1",
                traces[2].visible     = 'legendonly';
                traces[2].opacity     = 1.00;
                traces[2].line        = {color: 'green', width: 1.50, dash: 'solid' };
                traces[2].name        = LT3;                // legend title
                traces[2].showlegend  = LS3;       // Show legend

                traces[3].x           = tD4;
                traces[3].y           = res3.Data;
                traces[3].mode        = 'lines',
                traces[3].marker      = { color: 'orange', size: 5, symbol: 'circle' },
                traces[3].yaxis       = "y1",
                traces[3].visible     = 'legendonly';
                traces[3].opacity     = 1.00;
                traces[3].line        = {color: 'orange', width: 1.50, dash: 'solid' };
                traces[3].name        = LT4;               // legend title
                traces[3].showlegend  = LS4;      // Show legend
                
                layout_update.yaxis.title.text      = res.yTitle;   // This is the unit that user wants to see on the graph.
                layout_update.yaxis2.showticklabels = false;
                layout_update.yaxis2.title.text     = "";

            } 
            else if (IsFFT_CheckBox_Selected) {
                // Show FFT of Filtered Data and FFT of Raw Data 
                    
                // Declaration of variables 
                let Mag_FD_0=[], Angle_FD_0=[], f_0=[];
                let Mag_FD_1=[], Angle_FD_1=[], f_1=[];
                let Mag_FD_2=[], Angle_FD_2=[], f_2=[];
                let Mag_FD_3=[], Angle_FD_3=[], f_3=[];
                
                // FFT of SDOF-data
                let FSamp = 1 / ChannelList[ChNum].Results.SDOF.SamplingInterval;
                [Mag_FD_0, Angle_FD_0, f_0] = FourierSpec(res.Data,   FSamp);
                if (res1.Data != null || res1.Data != undefined) { [Mag_FD_1, Angle_FD_1, f_1] = FourierSpec(res1.Data,  FSamp); }
                if (res2.Data != null || res2.Data != undefined) { [Mag_FD_2, Angle_FD_2, f_2] = FourierSpec(res2.Data,  FSamp); }
                if (res3.Data != null || res3.Data != undefined) { [Mag_FD_3, Angle_FD_3, f_3] = FourierSpec(res3.Data,  FSamp); }
                
                // Plot FFT Magnitude Response of Raw data
                traces[0].x           = f_0;
                traces[0].y           = Mag_FD_0;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = LT1;         // legend title
                traces[0].showlegend  = LS1;         // Show legend 

                // Plot FFT Magnitude Response of Raw data
                traces[1].x           = f_1;
                traces[1].y           = Mag_FD_1;
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'red', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y1",
                traces[1].visible     = 'legendonly';
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'red', width: 0.80, dash: 'solid' };
                traces[1].name        = LT2;         // legend title
                traces[1].showlegend  = LS2;         // Show legend 

                // Plot FFT Magnitude Response of Raw data
                traces[2].x           = f_2;
                traces[2].y           = Mag_FD_2;
                traces[2].mode        = 'lines',
                traces[2].marker      = { color: 'green', size: 5, symbol: 'circle' },
                traces[2].yaxis       = "y1",
                traces[2].visible     = 'legendonly';
                traces[2].opacity     = 1.00;
                traces[2].line        = {color: 'green', width: 0.80, dash: 'solid' };
                traces[2].name        = LT3;           // legend title
                traces[2].showlegend  = LS3;           // Show legend 

                traces[3].x           = f_2;
                traces[3].y           = Mag_FD_2;
                traces[3].mode        = 'lines',
                traces[3].marker      = { color: 'orange', size: 5, symbol: 'circle' },
                traces[3].yaxis       = "y1",
                traces[3].visible     = 'legendonly';
                traces[3].opacity     = 1.00;
                traces[3].line        = {color: 'orange', width: 0.80, dash: 'solid' };
                traces[3].name        = LT4;           // legend title
                traces[3].showlegend  = LS4;           // Show legend 

                layout_update.xaxis.autorange = true;
                layout_update.yaxis.autorange = true;

                layout_update.yaxis.title.text      = res.yTitle_FFT; 
                layout_update.yaxis2.showticklabels = false;
                layout_update.yaxis2.title.text     = "";

            }

            // Update the Statistics of (RawData, Velocity, Displacement) in table - scaled to user-specified unit
            document.getElementById( Statictics_Peak_ID ).innerHTML = res.Peak.toPrecision(4);
            document.getElementById( Statictics_Mean_ID ).innerHTML = res.Mean.toPrecision(4);
            document.getElementById( Statictics_RMS_ID  ).innerHTML = res.RMS.toPrecision(4);

            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = ChannelList[ChNum].Results.SDOF.FiltPar.BaselineCorrection_String;

            // Show Filter_ID-Row in InfoBar
            FilterInfo  = ChannelList[ChNum].Results.SDOF.FiltPar.FilterName_String;
            FilterInfo += "<br>" + ChannelList[ChNum].Results.SDOF.FiltPar.FilterType_String;
            FilterInfo += " " + ChannelList[ChNum].Results.SDOF.FiltPar.FilterBand;
            FilterInfo += "<br> Zero Phase: " + ChannelList[ChNum].Results.SDOF.FiltPar.ZeroPhase;
            document.getElementById(FilterType_ID).innerHTML = FilterInfo;

            // Assign the analysis method
            document.getElementById(SDOF_Method_ID).innerHTML = ChannelList[ChNum].Results.SDOF.AnalysisMethod_string;

        }
        else {

            layout_update.yaxis.title.text      = '';   // This is the unit that user wants to see on the graph.
            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";

            // Empty the Statistics
            document.getElementById( Statictics_Peak_ID ).innerHTML = '';
            document.getElementById( Statictics_Mean_ID ).innerHTML = '';
            document.getElementById( Statictics_RMS_ID  ).innerHTML = '';

            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = '';

            // Show Filter_ID-Row in InfoBar
            document.getElementById(FilterType_ID).innerHTML = '';
        }
        
    }
    else if (PageNo == 4) {
        // Response Spectrum Page 

        // Get the index number of the SDOF_Analysis Method
        Indx = document.getElementById('ResSpec_AnalysisMethod').selectedIndex;

        // Make sure that the filter analysis is successfully completed
        if (ChannelList[ChNum].Results.ResSpec.IsAnalysisCompleted && (ChannelList[ChNum].Results.ResSpec.AnalysisMethod == Indx)) {

            if      (Indx == 0) { plotCount = ChannelList[ChNum].Results.ResSpec.DampingRatioCount;  temp = 'ksi'; }
            else if (Indx == 1) { plotCount = ChannelList[ChNum].Results.ResSpec.DuctilityCount;     temp = 'mu';  }
            else if (Indx == 2) { plotCount = ChannelList[ChNum].Results.ResSpec.DuctilityCount;     temp = 'mu';  }

            // Scale the data to the user-specified unit in Plotly Graph (info table)
            DisplayData = ChannelList[ChNum].Results.ResSpec.DisplayData;
            timeData    = ChannelList[ChNum].Results.ResSpec.T;

            if      (DisplayData == "Disp"    ) { res  = Convert_Data_To_Graph_Unit_ResSpec(ChannelList[ChNum].Results.ResSpec.SD,   ChNum );  }
            else if (DisplayData == "Vel"     ) { res  = Convert_Data_To_Graph_Unit_ResSpec(ChannelList[ChNum].Results.ResSpec.SV,   ChNum );  }
            else if (DisplayData == "acc"     ) { res  = Convert_Data_To_Graph_Unit_ResSpec(ChannelList[ChNum].Results.ResSpec.Sa,   ChNum );  }
            else if (DisplayData == "Acc"     ) { res  = Convert_Data_To_Graph_Unit_ResSpec(ChannelList[ChNum].Results.ResSpec.SA,   ChNum );  }
            else if (DisplayData == "PAcc"    ) { res  = Convert_Data_To_Graph_Unit_ResSpec(ChannelList[ChNum].Results.ResSpec.SPa,  ChNum );  }
            else if (DisplayData == "PVel"    ) { res  = Convert_Data_To_Graph_Unit_ResSpec(ChannelList[ChNum].Results.ResSpec.SPv,  ChNum );  }

            // Get checkbox status
            IsFilter_CheckBox_Selected = document.getElementById(FilterResp_ID).checked;

            if (!IsFilter_CheckBox_Selected && !IsFFT_CheckBox_Selected) {
                // Plot Spectrum in trace[0]
                if (plotCount >= 1) {
                    
                    traces[0].x           = timeData;
                    traces[0].y           = res.Data[0];
                    traces[0].mode        = 'lines+markers',
                    traces[0].marker      = { color: 'red', size: 5, symbol: 'circle' },
                    traces[0].yaxis       = "y1",
                    traces[0].visible     = true;
                    traces[0].opacity     = 1.00;
                    traces[0].line        = {color: 'red', width: 1.50, dash: 'solid' };
                    traces[0].name        = temp + " : " + ChannelList[ChNum].Results.ResSpec[temp + "_1"].toString();    // legend title
                    traces[0].showlegend  = true;                    // Show legend
                }

                // Plot Spectrum in tarace[1]
                if (plotCount >= 2) {
                    traces[1].x           = timeData;
                    traces[1].y           = res.Data[1];
                    traces[1].mode        = 'lines+markers',
                    traces[1].marker      = { color: 'green', size: 5, symbol: 'circle' },
                    traces[1].yaxis       = "y1",
                    traces[1].visible     = true;
                    traces[1].opacity     = 1.00;
                    traces[1].line        = {color: 'green', width: 1.50, dash: 'solid' };
                    traces[1].name        = temp + " : " + ChannelList[ChNum].Results.ResSpec[temp + "_2"].toString();  // legend title
                    traces[1].showlegend  = true;                   // Show legend 
                }
                
                // Plot Spectrum in tarace[2]
                if (plotCount >= 3) {
                    traces[2].x           = timeData;
                    traces[2].y           = res.Data[2];
                    traces[2].mode        = 'lines+markers',
                    traces[2].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                    traces[2].yaxis       = "y1",
                    traces[2].visible     = true;
                    traces[2].opacity     = 1.00;
                    traces[2].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                    traces[2].name        = temp + " : " + ChannelList[ChNum].Results.ResSpec[temp + "_3"].toString();    // legend title
                    traces[2].showlegend  = true;                   // Show legend 
                }

                // Plot Spectrum in tarace[3]
                if (plotCount >= 4) {
                    traces[3].x           = timeData;
                    traces[3].y           = res.Data[3];
                    traces[3].mode        = 'lines+markers',
                    traces[3].marker      = { color: 'orange', size: 5, symbol: 'circle' },
                    traces[3].yaxis       = "y1",
                    traces[3].visible     = true;
                    traces[3].opacity     = 1.00;
                    traces[3].line        = {color: 'orange', width: 1.50, dash: 'solid' };
                    traces[3].name        = temp + " : " + ChannelList[ChNum].Results.ResSpec[temp + "_4"].toString();     // legend title
                    traces[3].showlegend  = true;                   // Show legend 
                }

                layout_update.yaxis.title.text      = res.yTitle;   // This is the unit that user wants to see on the graph.
                layout_update.yaxis2.showticklabels = false;
                layout_update.yaxis2.title.text     = "";

                layout_update.yaxis.rangemode  = 'tozero';
            }
            else if (IsFilter_CheckBox_Selected) {

                // Plot Filter Magnitude Response in trace[1]
                traces[0].x           = ChannelList[ChNum].Results.ResSpec.FiltPar.f;
                traces[0].y           = ChannelList[ChNum].Results.ResSpec.FiltPar.Mag;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = '<b>Magnitude<b>';      // legend title 
                traces[0].showlegend  = true;                   // Show legend

                // Plot Filter Phase Response in trace[2]
                traces[1].x           = ChannelList[ChNum].Results.ResSpec.FiltPar.f;
                traces[1].y           = ChannelList[ChNum].Results.ResSpec.FiltPar.Angle;
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'green', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y2",
                traces[1].visible     = true;
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'green', width: 1.50, dash: 'solid' };
                traces[1].name        = '<b>Phase<b>';          // legend title 
                traces[1].showlegend  = true;                   // Show legend

                layout_update.yaxis.title.text      = res.yTitle_FFT;
                layout_update.yaxis2.showticklabels = true;
                layout_update.yaxis2.title.text     = res.y2Title;

            }

            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = ChannelList[ChNum].Results.ResSpec.FiltPar.BaselineCorrection_String;

            // Show Filter_ID-Row in InfoBar
            FilterInfo  = ChannelList[ChNum].Results.ResSpec.FiltPar.FilterName_String;
            FilterInfo += "<br>" + ChannelList[ChNum].Results.ResSpec.FiltPar.FilterType_String;
            FilterInfo += " " + ChannelList[ChNum].Results.ResSpec.FiltPar.FilterBand;
            FilterInfo += "<br> Zero Phase: " + ChannelList[ChNum].Results.ResSpec.FiltPar.ZeroPhase;
            document.getElementById(FilterType_ID).innerHTML = FilterInfo;

            // Assign the analysis method
            document.getElementById(SDOF_Method_ID).innerHTML = ChannelList[ChNum].Results.ResSpec.AnalysisMethod_string;
        }
        else {

            layout_update.yaxis.title.text      = '';   // This is the unit that user wants to see on the graph.
            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";

            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = '';

            // Show Filter_ID-Row in InfoBar
            document.getElementById(FilterType_ID).innerHTML = '';
        }

    }
    else if (PageNo == 5) {
        // Spectrum Page 

        if (ChannelList[ChNum].Results.Spectrum.IsAnalysisCompleted) {

            DisplayData = ChannelList[ChNum].Results.Spectrum.DisplayData;

            if      (DisplayData == "FFT"  ) { res  = Convert_Data_To_Graph_Unit_Spectrum(ChannelList[ChNum].Results.Spectrum.FFT,             ChNum );   }
            else if (DisplayData == "POW"  ) { res  = Convert_Data_To_Graph_Unit_Spectrum(ChannelList[ChNum].Results.Spectrum.PowerSpectrum,   ChNum );   }
            else if (DisplayData == "PSD"  ) { res  = Convert_Data_To_Graph_Unit_Spectrum(ChannelList[ChNum].Results.Spectrum.PSD,             ChNum );   }
            else if (DisplayData == "HET"  ) { res  = Convert_Data_To_Graph_Unit_Spectrum(ChannelList[ChNum].Results.Spectrum.Spectrogram,     ChNum );   }
            
            if (["FFT", "POW", "PSD"].includes(DisplayData)) {
                traces[0].x           = ChannelList[ChNum].Results.Spectrum.Freq_vector;
                traces[0].y           = res.Data;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.00, dash: 'solid' };
                traces[0].name        = '<b>'+ res.leg +'<b>';   // legend title
                traces[0].showlegend  = false;                   // Show legend
            }
            else if (["HET"].includes(DisplayData)) { 
                // Heat map  ( time-frequency graph)

                const z = ChannelList[ChNum].Results.Spectrum.Freq_vector.map((_, fi) => ChannelList[ChNum].Results.Spectrum.tVec.map((_, ti) => res.Data[ti][fi]));
                
                traces[0].type        = 'heatmap';
                traces[0].colorscale  = 'Bluered';
                traces[0].x           = ChannelList[ChNum].Results.Spectrum.tVec;
                traces[0].y           = ChannelList[ChNum].Results.Spectrum.Freq_vector;
                traces[0].z           = z;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.00, dash: 'solid' };
                traces[0].name        = '<b>'+ res.leg +'<b>';   // legend title
                traces[0].showlegend  = false;                   // Show legend

            }
            

            layout_update.yaxis.title.text      = res.yTitle;   // This is the unit that user wants to see on the graph.
            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";

            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = ChannelList[ChNum].Results.Spectrum.FiltPar.BaselineCorrection_String;

            // Show Filter_ID-Row in InfoBar
            FilterInfo  = ChannelList[ChNum].Results.Spectrum.FiltPar.FilterName_String;
            FilterInfo += "<br>" + ChannelList[ChNum].Results.Spectrum.FiltPar.FilterType_String;
            FilterInfo += " " + ChannelList[ChNum].Results.Spectrum.FiltPar.FilterBand;
            FilterInfo += "<br> Zero Phase: " + ChannelList[ChNum].Results.Spectrum.FiltPar.ZeroPhase;
            document.getElementById(FilterType_ID).innerHTML = FilterInfo;
            
        }


    }
    else if (PageNo == 6) {
        // SM Parameters Page 

        if (ChannelList[ChNum].Results.SM_Parameters.IsAnalysisCompleted) {

            DisplayData = ChannelList[ChNum].Results.SM_Parameters.DisplayData;
            timeData    = ChannelList[ChNum].time;

            if      (DisplayData == "AI"   ) { res  = Convert_Data_To_Graph_Unit_SM_Par(ChannelList[ChNum].Results.SM_Parameters.AI,  ChNum );   }
            else if (DisplayData == "CAV"  ) { res  = Convert_Data_To_Graph_Unit_SM_Par(ChannelList[ChNum].Results.SM_Parameters.CAV, ChNum );   }
            
            // SM Parameter
            traces[0].x           = timeData;
            traces[0].y           = res.Data;
            traces[0].mode        = 'lines',
            traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
            traces[0].yaxis       = "y1",
            traces[0].visible     = true;
            traces[0].opacity     = 1.00;
            traces[0].line        = {color: 'blue', width: 2.00, dash: 'solid' };
            traces[0].name        = '<b>'+ res.leg +'<b>';   // legend title
            traces[0].showlegend  = false;                   // Show legend

            // Raw data 
            traces[1].x           = timeData;
            traces[1].y           = Detrend(Multiply(ChannelList[ChNum].data, ChannelList[ChNum].ScaleFactor), 0);   // SM parameters are calculated using baseline corrcted data
            traces[1].mode        = 'lines',
            traces[1].marker      = { color: 'grey', size: 5, symbol: 'circle' },
            traces[1].yaxis       = "y2",
            traces[1].visible     = true;
            traces[1].opacity     = 0.80;
            traces[1].line        = {color: 'grey', width: 0.80, dash: 'solid' };
            traces[1].name        = '<b>'+ res.leg +'<b>';   // legend title
            traces[1].showlegend  = false;                   // Show legend

            
            if (DisplayData == "AI"   ) {

                T1 = ChannelList[ChNum].Results.SM_Parameters.T1_ai;   
                T2 = ChannelList[ChNum].Results.SM_Parameters.T2_ai;
                
                // update shapes 
                layout_update.shapes[0].x0 = T1;
                layout_update.shapes[0].x1 = T1;
                layout_update.shapes[0].y0 = 0;
                layout_update.shapes[0].y1 = 1;

                layout_update.shapes[1].x0 = T2;
                layout_update.shapes[1].x1 = T2;
                layout_update.shapes[1].y0 = 0;
                layout_update.shapes[1].y1 = 1;
            }

            let RD_Min = Min(traces[1].y).val;
            RD_Min -= (0.05 * RD_Min);

            layout_update.yaxis.title.text      = res.yTitle;   // This is the unit that user wants to see on the graph.
            layout_update.yaxis2.showticklabels = true;
            layout_update.yaxis2.title.text     = res.y2Title;

            layout_update.yaxis.autorange       = false;
            layout_update.yaxis.range           = [0, 1.05 * Max(traces[0].y).val];

            layout_update.yaxis2.autorange      = false;
            layout_update.yaxis2.matches        = null;
            layout_update.yaxis2.range          = [RD_Min,  1.05 * Max(traces[1].y).val];


            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = ChannelList[ChNum].Results.SM_Parameters.FiltPar.BaselineCorrection_String;

            // Show Filter_ID-Row in InfoBar
            FilterInfo  = ChannelList[ChNum].Results.SM_Parameters.FiltPar.FilterName_String;
            FilterInfo += "<br>" + ChannelList[ChNum].Results.SM_Parameters.FiltPar.FilterType_String;
            FilterInfo += " " + ChannelList[ChNum].Results.SM_Parameters.FiltPar.FilterBand;
            FilterInfo += "<br> Zero Phase: " + ChannelList[ChNum].Results.SM_Parameters.FiltPar.ZeroPhase;
            document.getElementById(FilterType_ID).innerHTML = FilterInfo;
            
            //
            Table = document.getElementById('Plotly_Stat_Table_' + ChannelList[ChNum].Unique_ID);
            Table.rows[3].cells[0].innerHTML = res.leg;
            Table.rows[3].cells[1].innerHTML = res.T_duration;

        }

    }
    else if (PageNo == 7) {
        // Drift 

        if (ChannelList[ChNum].Results.Drift.IsAnalysisCompleted) {

            // Get the data to plot 
            res  = Convert_Data_To_Graph_Unit_Drift(ChannelList[ChNum].Results.Drift.Drift, ChNum);   // Drift
            res2 = Convert_Data_To_Graph_Unit_Drift(ChannelList[ChNum].Results.Drift.Data,  ChNum);   // Displacment data channel-1 and channel-2

            // Get checkbox status
            IsFilter_CheckBox_Selected = document.getElementById(FilterResp_ID).checked;
            IsFFT_CheckBox_Selected    = document.getElementById(FilterFFT_ID).checked;

            if (!IsFilter_CheckBox_Selected && !IsFFT_CheckBox_Selected) {
                // Neither of the checkboxes are selected 

                traces[0].x           = ChannelList[ChNum].Results.Drift.time;
                traces[0].y           = res.Data;  // Drift data
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = '<b>Drift<b>';   // legend title
                traces[0].showlegend  = true;      // Show legend

                traces[1].x           = ChannelList[ChNum].Results.Drift.time;
                traces[1].y           = res2.Data[0];  // Data[0] - channel-1
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'red', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y1",
                traces[1].visible     = 'legendonly';
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'red', width: 0.80, dash: 'solid' };
                traces[1].name        = '<b>Channel-1<b>';   // legend title
                traces[1].showlegend  = true;      // Show legend

                traces[2].x           = ChannelList[ChNum].Results.Drift.time;
                traces[2].y           = res2.Data[1];  // Drift data
                traces[2].mode        = 'lines',
                traces[2].marker      = { color: 'green', size: 5, symbol: 'circle' },
                traces[2].yaxis       = "y1",
                traces[2].visible     = 'legendonly';
                traces[2].opacity     = 1.00;
                traces[2].line        = {color: 'green', width: 0.80, dash: 'solid' };
                traces[2].name        = '<b>Channel-2<b>';   // legend title
                traces[2].showlegend  = true;      // Show legend

                layout_update.yaxis.title.text      = res.yTitle;   // This is the unit that user wants to see on the graph.
                layout_update.yaxis2.showticklabels = false;
                layout_update.yaxis2.title.text     = '';

            }
            else if (IsFilter_CheckBox_Selected) {
                // Filter checkbox is selected 

                // Plot Filter Magnitude Response in trace[1]
                traces[0].x           = ChannelList[ChNum].Results.Drift.FiltPar.f;
                traces[0].y           = ChannelList[ChNum].Results.Drift.FiltPar.Mag;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = '<b>Magnitude<b>';      // legend title 
                traces[0].showlegend  = true;                   // Show legend

                // Plot Filter Phase Response in trace[2]
                traces[1].x           = ChannelList[ChNum].Results.Drift.FiltPar.f;
                traces[1].y           = ChannelList[ChNum].Results.Drift.FiltPar.Angle;
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'green', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y2",
                traces[1].visible     = true;
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'green', width: 1.50, dash: 'solid' };
                traces[1].name        = '<b>Phase<b>';          // legend title 
                traces[1].showlegend  = true;                   // Show legend

                layout_update.yaxis.title.text      = res.yTitle_FFT;
                layout_update.yaxis2.showticklabels = true;
                layout_update.yaxis2.title.text     = res.y2Title;

            }
            else if (IsFFT_CheckBox_Selected) {
                // FFT checkbox is selected 

                let Mag_FD   =[],  Angle_FD   =[],  f   =[];
                let Mag_FD_0 =[],  Angle_FD_0 =[],  f_0 =[];
                let Mag_FD_1 =[],  Angle_FD_1 =[],  f_1 =[];

                [Mag_FD, Angle_FD, f] = FourierSpec(res.Data,  ChannelList[ChNum].FSamp);  // FFT of Drift

                [Mag_FD_0, Angle_FD_0, f_0] = FourierSpec(res2.Data[0],  ChannelList[ChNum].FSamp);  // FFT of Channel-1 (disp)
                [Mag_FD_1, Angle_FD_1, f_1] = FourierSpec(res2.Data[1],  ChannelList[ChNum].FSamp);  // FFT of Channel-2 (disp)
                
                // trace[0]
                traces[0].x           = f;
                traces[0].y           = Mag_FD;
                traces[0].mode        = 'lines',
                traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
                traces[0].yaxis       = "y1",
                traces[0].visible     = true;
                traces[0].opacity     = 1.00;
                traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
                traces[0].name        = 'Drift';               // legend title
                traces[0].showlegend  = true;                  // Show legend

                // trace[1]
                traces[1].x           = f_0;
                traces[1].y           = Mag_FD_0;
                traces[1].mode        = 'lines',
                traces[1].marker      = { color: 'red', size: 5, symbol: 'circle' },
                traces[1].yaxis       = "y1",
                traces[1].visible     = 'legendonly';
                traces[1].opacity     = 1.00;
                traces[1].line        = {color: 'red', width: 0.80, dash: 'solid' };
                traces[1].name        = 'Channel-1';           // legend title
                traces[1].showlegend  = true;                  // Show legend

                // trace[1]
                traces[2].x           = f_1;
                traces[2].y           = Mag_FD_1;
                traces[2].mode        = 'lines',
                traces[2].marker      = { color: 'green', size: 5, symbol: 'circle' },
                traces[2].yaxis       = "y1",
                traces[2].visible     = 'legendonly';
                traces[2].opacity     = 1.00;
                traces[2].line        = {color: 'green', width: 0.80, dash: 'solid' };
                traces[2].name        = 'Channel-2';           // legend title
                traces[2].showlegend  = true;                  // Show legend

                layout_update.yaxis.title.text      = res.yTitle_FFT;
                layout_update.yaxis2.showticklabels = false;
                layout_update.yaxis2.title.text     = "";

            }

            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = ChannelList[ChNum].Results.Drift.FiltPar.BaselineCorrection_String;

            // Show Filter_ID-Row in InfoBar
            FilterInfo  = ChannelList[ChNum].Results.Drift.FiltPar.FilterName_String;
            FilterInfo += "<br>" + ChannelList[ChNum].Results.Drift.FiltPar.FilterType_String;
            FilterInfo += " " + ChannelList[ChNum].Results.Drift.FiltPar.FilterBand;
            FilterInfo += "<br> Zero Phase: " + ChannelList[ChNum].Results.Drift.FiltPar.ZeroPhase;
            document.getElementById(FilterType_ID).innerHTML = FilterInfo;

            // Update the Statistics of (RawData, Velocity, Displacement) in table - scaled to user-specified unit
            document.getElementById( Statictics_Peak_ID     ).innerHTML = res.Peak.toPrecision(4);
            document.getElementById( Statictics_Mean_ID     ).innerHTML = res.Mean.toPrecision(4);
            document.getElementById( Statictics_RMS_ID      ).innerHTML = res.RMS.toPrecision(4);
            document.getElementById( Statictics_Residual_ID ).innerHTML = res.Residual.toPrecision(4);
            
        }

    }
    else if (PageNo == 8) {
        // H/V Spectral Ratio

        if (ChannelList[ChNum].Results.HVSR.IsAnalysisCompleted) {

            traces[0].x           = ChannelList[ChNum].Results.HVSR.f;
            traces[0].y           = ChannelList[ChNum].Results.HVSR.HV;
            traces[0].mode        = 'lines',
            traces[0].marker      = { color: 'blue', size: 5, symbol: 'circle' },
            traces[0].yaxis       = "y1",
            traces[0].visible     = true;
            traces[0].opacity     = 1.00;
            traces[0].line        = {color: 'blue', width: 1.50, dash: 'solid' };
            traces[0].name        = '<b><b>';   // legend title
            traces[0].showlegend  = false;      // Show legend
            

            traces[1].x           = ChannelList[ChNum].Results.HVSR.f;
            traces[1].y           = ChannelList[ChNum].Results.HVSR.HV.map((v, i) => v * Math.exp(ChannelList[ChNum].Results.HVSR.Std[i]));
            traces[1].mode        = 'lines',
            traces[1].marker      = { color: 'black', size: 5, symbol: 'circle' },
            traces[1].yaxis       = "y1",
            traces[1].visible     = true;
            traces[1].opacity     = 1.00;
            traces[1].line        = {color: 'black', width: 0.80, dash: 'dash' };
            traces[1].name        = '<b><b>';   // legend title
            traces[1].showlegend  = false;      // Show legend

            traces[2].x           = ChannelList[ChNum].Results.HVSR.f;
            traces[2].y           = ChannelList[ChNum].Results.HVSR.HV.map((v, i) => v * Math.exp(-ChannelList[ChNum].Results.HVSR.Std[i]));
            traces[2].mode        = 'lines',
            traces[2].marker      = { color: 'black', size: 5, symbol: 'circle' },
            traces[2].yaxis       = "y1",
            traces[2].visible     = true;
            traces[2].opacity     = 1.00;
            traces[2].line        = {color: 'black', width: 0.80, dash: 'dash' };
            traces[2].name        = '<b><b>';   // legend title
            traces[2].showlegend  = false;      // Show legend

            layout_update.yaxis.type            = "linear";
            layout_update.yaxis.autorange       = false;
            layout_update.yaxis.range           = [0, 1.05 * Max(traces[1].y.slice(1)).val];

            layout_update.yaxis2.type            = "linear";
            layout_update.yaxis2.autorange       = false;
            layout_update.yaxis2.range           = [0, 1.05 * Max(traces[1].y.slice(1)).val];

            layout_update.xaxis.type            = "log";
            layout_update.xaxis.dtick           = 1;
            layout_update.xaxis.range           = [ChannelList[ChNum].Results.HVSR.f[1], ChannelList[ChNum].Results.HVSR.f.at(-1)];

            layout_update.yaxis.title.text      = '<b>H/V Amplitude</b>';
            layout_update.yaxis2.showticklabels = false;
            layout_update.yaxis2.title.text     = "";

            layout_update.xaxis.title.text      = '<b>Frequency (Hz)</b>';

            // Show Baseline-Row in InforBar
            document.getElementById(BaseLine_ID).innerHTML = ChannelList[ChNum].Results.HVSR.FiltPar.BaselineCorrection_String;

            // Show Filter_ID-Row in InfoBar
            FilterInfo  = ChannelList[ChNum].Results.HVSR.FiltPar.FilterName_String;
            FilterInfo += "<br>" + ChannelList[ChNum].Results.HVSR.FiltPar.FilterType_String;
            FilterInfo += " " + ChannelList[ChNum].Results.HVSR.FiltPar.FilterBand;
            FilterInfo += "<br> Zero Phase: " + ChannelList[ChNum].Results.HVSR.FiltPar.ZeroPhase;
            document.getElementById(FilterType_ID).innerHTML = FilterInfo;

            // Infor Table - first two lines
            


        }

    }

    // Update the graph
    Plotly.update(PlotArea_ID, traces, layout_update);

}
//-------------------------------------------------------------------------------------------------------------
async function Plotly_Create_Graph(Container_Id, Channel) {

    // This function create a new div-element that contains the Plotly-graph and the InforDiv 

    //Decleration of variables 
    let SubContainer, PlotlyDiv, InfoDiv, Span1, Span2, Table, Span2_Info

    // Create a div-element that will contain the following two  div-elements
    // The first div-element is to host Plotly-graph
    // The second div-element is to host Info-table
    SubContainer = document.createElement('div');
    SubContainer.setAttribute('id',    'Div_ID_' + Channel.Unique_ID);
    SubContainer.setAttribute('class', 'Plotly_Main_Container'      );

    // The first div-element Plotly graph
    PlotlyDiv           = document.createElement('div');
    PlotlyDiv.setAttribute('id',    'PlotArea_ID_' + Channel.Unique_ID);
    PlotlyDiv.setAttribute('class', 'Plotly_Div');

    // The second div-element Info-div
    InfoDiv = document.createElement('div');
    InfoDiv.setAttribute('class', 'Plotly_Info_Div');

    // InfoDiv contains 3 elements (2 spans element and 2 table element)
    Span1 = document.createElement('span');
    Span1.setAttribute('class', 'Plotly_Span_Title');
    Span1.innerHTML = Channel.FileName;

    Span2 = document.createElement('span');
    Span2.setAttribute('class', 'Plotly_Span_ChannelInfo');
    Span2_Info  = "(Ch-" + Channel.ChNum + ") (" +  Channel.Orientation  + ") (" + Number(Channel.FSamp).toFixed(3).toString() + " Hz)";
    Span2.innerHTML = Span2_Info;

    // Create a Table
    // wait for Info_Table() to be resolved and continue with the rest of the code
    Table = await Plotly_Info_Table(Channel);

    InfoDiv.appendChild(Span1);
    InfoDiv.appendChild(Span2);
    InfoDiv.appendChild(Table);

    // Append the two div-elements to SubContainer
    SubContainer.appendChild(PlotlyDiv);
    SubContainer.appendChild(InfoDiv);

    // Apend the SubContainer to the main-container
    document.getElementById(Container_Id).appendChild(SubContainer);
    
    // Create new Plotly graph using the PlotlyDiv element
    // wait for Create_NewGraph() to be resolved and continue with the rest of the code 
    await Plotly_NewGraph(PlotlyDiv, Channel);

    // Change the legend text for trace1 (raw-data)
    Plotly.restyle(PlotlyDiv, { name: Channel.FileName +'<br>'+ Span2_Info }, [0]);

}
//-------------------------------------------------------------------------------------------------------------
async function Plotly_NewGraph(Div_ID, Channel) {

    // This function creates a new Plotly grapgs with traces.
    // The first traces is filled with the raw-data of the channel.
    // Other two traces are left emty for use in other pages such as to plot filtered-data or itegrated-data (velocity, displacment)

    // Declaration of variables
    let traces=[], trace1, trace2, trace3, trace4, layout, config;
    let GraphTitle = '<b>' + Channel.FileName +' (Ch-' + Channel.ChNum + ') (' + Channel.Orientation  + ') (' + Number(Channel.FSamp).toFixed(3).toString() + ' Hz) <b>';
    let yLabel  = "<b>" + Channel.TypeString + "   [" + Channel.UnitString + "] <b>";
    let y2Label = "<b><b>";
    let xLabel  = "<b><b>";

    // Define Traces (Total of 3 Traces per graph)
    trace1 = {
        x             : Channel.time,
        y             : Multiply(Channel.ScaleFactor, Channel.data),
        mode          : 'lines',
        type          : 'scatter',
        yaxis         : "y1",
        name          : '<b>Raw Data<b>',                              // Legend name
        opacity       : 1.0,
        visible       : true,                                          // Show this trace
        line          : {color: 'blue', width: 1.00, dash: 'solid' },
        showlegend    : false,
    };
    trace2 = {
        x             : [],
        y             : [],
        mode          : 'lines',
        type          : 'scatter',
        yaxis         : "y1",
        name          : '<b>Filtered Data<b>',  // Legend name
        opacity       : 1.0,
        visible       : false,                                          // Hide this trace
        line          : {color: 'red', width: 1.00, dash: 'solid' },
        showlegend    : false,
    };
    trace3 = {
        x             : [],
        y             : [],
        mode          : 'lines',
        type          : 'scatter',
        yaxis         : "y2",
        name          : '<b>Phase<b>',                                  // Legend name
        opacity       : 1.0,
        visible       : false,                                          // Hide this trace
        line          : {color: 'green', width: 1.00, dash: 'solid' },
        showlegend    : false,
    };
    trace4 = {
        x             : [],
        y             : [],
        mode          : 'lines',
        type          : 'scatter',
        yaxis         : "y2",
        name          :  '',                                            // Legend name
        opacity       : 1.0,
        visible       : false,                                          // Hide this trace
        line          : {color: 'green', width: 1.00, dash: 'solid' },
        showlegend    : false,
    };
    traces.push(trace1);
    traces.push(trace2);
    traces.push(trace3);
    traces.push(trace4);

    // Define Layout
    layout = {
        title           : { text: GraphTitle, font: {size: 10 }, x: 0.5, xanchor: 'center', y: 0.98, yanchor: 'top'},
        xaxis           : { zeroline: false, automargin: true, tickfont: { size: 15 },                    linecolor: 'black', linewidth: 1, mirror: true, title: {text: xLabel,  standoff: 5, font: {family: "Arial", size: 17} }, autorange: true, type: 'linear' },
        yaxis           : { zeroline: true,  automargin: true, tickfont: { size: 15 }, tickformat: '.2e', linecolor: 'black', linewidth: 1, mirror: true, title: {text: yLabel,  standoff: 5, font: {family: "Arial", size: 17} }, autorange: true, rangemode: 'normal', type: 'linear' },
        yaxis2          : { zeroline: true,  automargin: true, tickfont: { size: 15 },                    linecolor: 'black', linewidth: 1, mirror: true, title: {text: y2Label, standoff: 5, font: {family: "Arial", size: 17} }, autorange: true, overlaying: 'y', side: 'right', showticklabels: false, matches: "y", rangemode: 'normal', type: 'linear'},
        plot_bgcolor    : '#ffffff', 
        paper_bgcolor   : '#ffffff',
        legend          : { x: 0.99, y:0.85, xanchor: 'right',  orientation: 'v', font: {size: 14, weight: 700}, bgcolor: '#f8f4f4', bordercolor: '#6c6a6a', borderwidth: 1.5, },
        autosize        : true,
        margin          : {t: 20, r:20, b:5, l:5},
        shapes          : [ { type: 'line',  x0: 0, x1: 0,  y0: 0, y1: 0,  line: { color: 'red', width: 1.5, dash: 'dash' }, yref: 'paper', layer: 'below'},
                            { type: 'line',  x0: 0, x1: 0,  y0: 0, y1: 0,  line: { color: 'red', width: 1.5, dash: 'dash' }, yref: 'paper', layer: 'below'},
                        ],
    };

    config = {
        edits                   : { legendPosition: true, }, // Allow moving/dragging the legend on plot
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
//-------------------------------------------------------------------------------------------------------------
async function Plotly_Info_Table(Channel) {

    let Tabel, Table_Body, row, cell, opt;
    let Statictics_Peak_ID, Statictics_Mean_ID, Statictics_RMS_ID, Statictics_Residual_ID;
    let FilterType_ID, FilterRow_ID, BaseLine_ID, BaseLineRow_ID, FilterResp_ID, GraphUnitRow_ID,FilterFFT_ID;
    let SDOF_Row_ID, SDOF_Plot_ID, SDOF_Cell_ID, SDOF_Method_ID, SDOF_Method_Row_ID;
    let Span, input, div1, div2, label, Unit_List, j, Unit_Cell_ID, Unit_Plot_ID, select;
    let ResSpec_Row_ID, ResSpec_Plot_ID, ResSpec_Cell_ID, Filter_Div_ID1,Filter_Div_ID2;
    let Spectrum_Row_ID, Spectrum_Plot_ID, Spectrum_Cell_ID;
    let SM_Par_Row_ID, SM_Par_Plot_ID, SM_Par_Cell_ID;

    Statictics_Peak_ID     = "Statictics_Peak_ID_" + Channel.Unique_ID;
    Statictics_Mean_ID     = "Statictics_Mean_ID_" + Channel.Unique_ID;
    Statictics_RMS_ID      = "Statictics_RMS_ID_" + Channel.Unique_ID;
    Statictics_Residual_ID = "Statictics_Residual_ID_" + Channel.Unique_ID;

    GraphUnitRow_ID        = "GraphUnitRow_ID_" + Channel.Unique_ID;
    Unit_Plot_ID           = "Unit_Plot_ID_" + Channel.Unique_ID;
    Unit_Cell_ID           = "Unit_Cell_ID_" + Channel.Unique_ID;
    SDOF_Row_ID            = "SDOF_Row_ID_" + Channel.Unique_ID;
    SDOF_Plot_ID           = "SDOF_Plot_ID_" + Channel.Unique_ID;
    SDOF_Cell_ID           = "SDOF_Cell_ID_" + Channel.Unique_ID;

    ResSpec_Row_ID         = "ResSpec_Row_ID_" + Channel.Unique_ID;
    ResSpec_Plot_ID        = "ResSpec_Plot_ID_" + Channel.Unique_ID;
    ResSpec_Cell_ID        = "ResSpec_Cell_ID_" + Channel.Unique_ID;

    Spectrum_Row_ID        = "Spectrum_Row_ID_" + Channel.Unique_ID;
    Spectrum_Plot_ID       = "Spectrum_Plot_ID_" + Channel.Unique_ID;
    Spectrum_Cell_ID       = "Spectrum_Cell_ID_" + Channel.Unique_ID;

    SM_Par_Row_ID          = "SM_Par_Row_ID_" + Channel.Unique_ID;
    SM_Par_Plot_ID         = "SM_Par_Plot_ID_" + Channel.Unique_ID;
    SM_Par_Cell_ID         = "SM_Par_Cell_ID_" + Channel.Unique_ID;

    SDOF_Method_ID         = "SDOF_Method_ID_" + Channel.Unique_ID;
    SDOF_Method_Row_ID     = "SDOF_Method_Row_ID_" + Channel.Unique_ID;
    BaseLine_ID            = "BaseLine_ID_" + Channel.Unique_ID;
    BaseLineRow_ID         = "BaseLineRow_ID_" + Channel.Unique_ID;
    FilterType_ID          = "FilterType_ID_" + Channel.Unique_ID;
    FilterRow_ID           = "FilterRow_ID_" + Channel.Unique_ID;
    FilterResp_ID          = "FilterResp_ID_" + Channel.Unique_ID;
    FilterFFT_ID           = "FilterFFT_ID_" + Channel.Unique_ID;
    Filter_Div_ID1         = "Filter_Div_ID1_" + Channel.Unique_ID;
    Filter_Div_ID2         = "Filter_Div_ID2_" + Channel.Unique_ID;

    
    // Create a table and tabelBody
    Tabel = document.createElement('table');
    Tabel.setAttribute('class', 'Plotly_Stat_Table');
    Tabel.setAttribute('id', 'Plotly_Stat_Table_'+Channel.Unique_ID);
    Tabel.style.tableLayout = 'fixed';
    Table_Body = document.createElement('tbody');
    Table_Body.setAttribute('class', '')
    Tabel.appendChild(Table_Body);

    // Create a new row for Peak-value
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_Peak_'+Channel.Unique_ID);
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    
    cell.style.width = '5.5rem';

    cell.innerHTML = "Peak";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.setAttribute('id', Statictics_Peak_ID);
    cell.innerHTML = Channel.Peak.toPrecision(4);

    // Create a new row for Mean-value
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_Mean_'+Channel.Unique_ID);
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Mean";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.setAttribute('id', Statictics_Mean_ID);
    cell.innerHTML = Channel.Mean.toPrecision(4);  

    // Create a new row for RMS-value
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_RMS_'+Channel.Unique_ID);
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "RMS";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.innerHTML = Channel.RMS.toPrecision(4);
    cell.setAttribute('id', Statictics_RMS_ID);

    // Create a new row for Residual-value
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_Residual_'+Channel.Unique_ID);
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Residual";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.innerHTML = Channel.Residual.toPrecision(4);
    cell.setAttribute('id', Statictics_Residual_ID);



    // Create a new row for SDOF Analysis Method
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_AnalysisMethod_'+Channel.Unique_ID);
    row.setAttribute('id', SDOF_Method_Row_ID);
    row.style.display = "none";
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Method";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.innerHTML = '';
    cell.id        = SDOF_Method_ID;



    // Create a new row for SDOF-List
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_Display_'+Channel.Unique_ID);
    row.setAttribute('id', SDOF_Row_ID);
    row.style.display = "none";

    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Display";

    select = document.createElement('select');
    select.setAttribute("id", SDOF_Plot_ID);
    select.setAttribute('class', 'form-custom');
    select.setAttribute('onchange', '');
    select.selectedIndex = 0;

    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.setAttribute('id', SDOF_Cell_ID);
    cell.appendChild(select);



    // Create a new row for ResSpec-List
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_ResSpec_Display_'+Channel.Unique_ID);
    row.setAttribute('id', ResSpec_Row_ID);
    row.style.display = "none";

    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Display";

    select = document.createElement('select');
    select.setAttribute("id", ResSpec_Plot_ID);
    select.setAttribute('class', 'form-custom');
    select.setAttribute('onchange', '');
    select.selectedIndex = 0;

    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.setAttribute('id', ResSpec_Cell_ID);
    cell.appendChild(select);


    // Create a new row for Spectrum-List
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_Spectrum_Display_'+Channel.Unique_ID);
    row.setAttribute('id', Spectrum_Row_ID);
    row.style.display = "none";

    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Display";

    select = document.createElement('select');
    select.setAttribute("id", Spectrum_Plot_ID);
    select.setAttribute('class', 'form-custom');
    select.setAttribute('onchange', '');
    select.selectedIndex = 0;

    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.setAttribute('id', Spectrum_Cell_ID);
    cell.appendChild(select);



    // Create a new row for SM_Parameters-List
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_SM_Par_Display_'+Channel.Unique_ID);
    row.setAttribute('id', SM_Par_Row_ID);
    row.style.display = "none";

    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Display";

    select = document.createElement('select');
    select.setAttribute("id", SM_Par_Plot_ID);
    select.setAttribute('class', 'form-custom');
    select.setAttribute('onchange', '');
    select.selectedIndex = 0;

    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.setAttribute('id', SM_Par_Cell_ID);
    cell.appendChild(select);


    // Create a new row for Graph-Unit
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_GraphUnit_'+Channel.Unique_ID);
    row.setAttribute('id', GraphUnitRow_ID);
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Unit";
    
    Unit_List  = List_Units(TypeAndUnit(Channel.TypeAndUnits).Unit);
    select     = document.createElement('select');
    select.setAttribute("id", Unit_Plot_ID);
    select.setAttribute('class', 'form-custom');
    for (j = 0; j < Unit_List.Units.length; j++) {
        opt = document.createElement("option");
        opt.text = Unit_List.Units[j];
        select.add(opt, null);
    }
    select.setAttribute('onchange', 'Plotly_Graph_Update('+ ChannelList_UniqueID(Channel.Unique_ID) + ')' );
    select.selectedIndex = Unit_List.UnitNum.indexOf(Channel.Unit) ;
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.setAttribute('id', Unit_Cell_ID);
    cell.appendChild(select);

    
    // Create a new row for Baseline
    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_Baseline_'+Channel.Unique_ID);
    row.setAttribute('id', BaseLineRow_ID);
    row.style.display = "none";
    cell = row.insertCell(0);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Left');
    cell.innerHTML = "Baseline";
    cell = row.insertCell(1);
    cell.setAttribute('class', 'Plotly_Stat_Body_Td_Right');
    cell.innerHTML = '';
    cell.id        = BaseLine_ID;
    
    // Create a new row for Filter and FFT
    Span = document.createElement('span');
    Span.setAttribute('class', 'Filter_Div_Span');
    Span.textContent = 'Filter';
    input = document.createElement('input');
    input.setAttribute('class', 'toggle');
    input.setAttribute('type', 'checkbox');
    input.setAttribute('id', FilterResp_ID);
    input.setAttribute('onclick', 'Plotly_Filter_Response(this)');
    label = document.createElement('label');
    label.setAttribute('for', FilterResp_ID);
    label.setAttribute('class', "toggle-label");
    div1 = document.createElement('div');
    div1.setAttribute('class', "Filter_Div");
    div1.setAttribute('id', Filter_Div_ID1);
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
    input.setAttribute('onclick', 'Plotly_Filter_Response(this)');
    label = document.createElement('label');
    label.setAttribute('for', FilterFFT_ID);
    label.setAttribute('class', "toggle-label");
    div2 = document.createElement('div');
    div2.setAttribute('class', "Filter_Div");
    div2.setAttribute('id', Filter_Div_ID2);
    div2.appendChild(Span);
    div2.appendChild(input);
    div2.appendChild(label);

    row = Tabel.insertRow(-1);
    row.setAttribute('id', 'Row_Filter_'+Channel.Unique_ID);
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

    // Adjust the width of the first column in the Table 
    for (let row of Tabel.rows) { row.cells[0].style.width = '8rem'; }

    return Tabel;
}
//-------------------------------------------------------------------------------------------------------------
function Plotly_Filter_Response(el) {

    // Declartion of variables 
    let ChNum;

    if (el.id.includes('FilterResp_ID_')) { 

        ChNum = ChannelList_UniqueID(el.id.replace('FilterResp_ID_', ""));
        document.getElementById('FilterFFT_ID_' + ChannelList[ChNum].Unique_ID).checked = false;    // uncheck FFT

    }
    else if (el.id.includes('FilterFFT_ID_')) { 

        ChNum  = ChannelList_UniqueID(el.id.replace('FilterFFT_ID_',  ""));
        document.getElementById('FilterResp_ID_' + ChannelList[ChNum].Unique_ID).checked = false;   // uncheck Filter Response

    }
    
    // Update Graph
    Plotly_Graph_Update(ChNum);
}
//-------------------------------------------------------------------------------------------------------------
function Plotly_Clear_Graph(ChNum) {

    // Return if no graph is created for this channel
    if (!ChannelList[ChNum].PlotGraph) { return; } 

    // Decleration of variables 
    let i, traces, layout_update, PlotArea_ID, SDOF_Plot_ID, Table;
    
    PlotArea_ID   = "PlotArea_ID_"   + ChannelList[ChNum].Unique_ID;

    // Get the existing traces in the PlotArea_ID
    traces = document.getElementById(PlotArea_ID).data;

    // Get the existing layout in the PlotArea_ID
    layout_update = document.getElementById(PlotArea_ID).layout;
    layout_update.yaxis.type       = "linear";
    layout_update.xaxis.type       = "linear";
    layout_update.yaxis.dtick      = null;
    layout_update.xaxis.dtick      = null;

    layout_update.xaxis.autorange  = true;
    layout_update.yaxis.autorange  = true;
    layout_update.yaxis2.autorange = true;

    layout_update.xaxis.rangemode  = 'normal';
    layout_update.yaxis.rangemode  = 'normal';
    layout_update.yaxis2.rangemode = 'normal';

    layout_update.yaxis2.matches   = 'y';

    layout_update.xaxis.title.text = '';


    // Delete shapes
    layout_update.shapes[0].x0 = 0;
    layout_update.shapes[0].x1 = 0;
    layout_update.shapes[0].y0 = 0;
    layout_update.shapes[0].y1 = 0;
    layout_update.shapes[1].x0 = 0;
    layout_update.shapes[1].x1 = 0;
    layout_update.shapes[1].y0 = 0;
    layout_update.shapes[1].y1 = 0;

    // Delete all traces
    for (i=0; i<traces.length; i++) {
        traces[i].type        = 'scatter';
        traces[i].x           = [];
        traces[i].y           = [];
        traces[i].mode        = 'lines',
        traces[i].marker      = { color: 'blue', size: 5, symbol: 'circle' },
        traces[i].yaxis       = "y1",
        traces[i].visible     = false;
        traces[i].opacity     = 0.35;
        traces[i].line        = {color: 'blue', width: 1.50, dash: 'solid',dash: 'solid' };
        traces[i].name        = '<b>Raw Data<b>';    // legend title
        traces[i].showlegend  = false;               // Show legend
    }

    // Get the Infor Table 
    Table = document.getElementById('Plotly_Stat_Table_' + ChannelList[ChNum].Unique_ID);

    Table.rows[3].cells[0].innerHTML = 'Residual';

    if (PageNo == 0) {
        // Load Data Page 
        Table.rows[0].style.display  = "table-row";  // Peak
        Table.rows[1].style.display  = "table-row";  // Mean
        Table.rows[2].style.display  = "table-row";  // RMS
        Table.rows[3].style.display  = "none";       // Residual
        Table.rows[4].style.display  = "none";       // Analysis Method
        Table.rows[5].style.display  = "none";       // SDOF Display
        Table.rows[6].style.display  = "none";       // ResSpec Display
        Table.rows[7].style.display  = "none";       // Spectrum Display
        Table.rows[8].style.display  = "none";       // SM_Par Display
        Table.rows[9].style.display  = "table-row";  // Graph Unit
        Table.rows[10].style.display = "none";       // Baseline 
        Table.rows[11].style.display = "none";       // Filter 

        // Two toggels on the Info Summary table for (FFT) and (Filter) checkboxes
        document.getElementById('Filter_Div_ID1_'+ChannelList[ChNum].Unique_ID).style.display = 'none';
        document.getElementById('Filter_Div_ID2_'+ChannelList[ChNum].Unique_ID).style.display = 'none';

    }
    else if (PageNo == 1) {
        // Filter Page 
        Table.rows[0].style.display  = "table-row";  // Peak
        Table.rows[1].style.display  = "table-row";  // Mean
        Table.rows[2].style.display  = "table-row";  // RMS
        Table.rows[3].style.display  = "none";       // Residual
        Table.rows[4].style.display  = "none";       // Analysis Method
        Table.rows[5].style.display  = "none";       // SDOF Display
        Table.rows[6].style.display  = "none";       // ResSpec Display
        Table.rows[7].style.display  = "none";       // Spectrum Display
        Table.rows[8].style.display  = "none";       // SM_Par Display
        Table.rows[9].style.display  = "table-row";  // Graph Unit
        Table.rows[10].style.display = "table-row";  // Baseline 
        Table.rows[11].style.display = "table-row";  // Filter
        
        // Two toggels on the Info Summary table for (FFT) and (Filter) checkboxes
        document.getElementById('Filter_Div_ID1_'+ChannelList[ChNum].Unique_ID).style.display = 'flex';
        document.getElementById('Filter_Div_ID2_'+ChannelList[ChNum].Unique_ID).style.display = 'flex';

    }
    else if (PageNo == 2) {
        // Integration page 
        Table.rows[0].style.display  = "table-row";  // Peak
        Table.rows[1].style.display  = "table-row";  // Mean
        Table.rows[2].style.display  = "table-row";  // RMS
        Table.rows[3].style.display  = "table-row";  // Residual
        Table.rows[4].style.display  = "none";       // Analysis Method
        Table.rows[5].style.display  = "none";       // SDOF Display
        Table.rows[6].style.display  = "none";       // ResSpec Display
        Table.rows[7].style.display  = "none";       // Spectrum Display
        Table.rows[8].style.display  = "none";       // SM_Par Display
        Table.rows[9].style.display  = "table-row";  // Graph Unit
        Table.rows[10].style.display = "table-row";  // Baseline 
        Table.rows[11].style.display = "table-row";  // Filter 

        // Two toggels on the Info Summary table for (FFT) and (Filter) checkboxes
        document.getElementById('Filter_Div_ID1_'+ChannelList[ChNum].Unique_ID).style.display = 'flex';
        document.getElementById('Filter_Div_ID2_'+ChannelList[ChNum].Unique_ID).style.display = 'flex';

    }
    else if (PageNo == 3) {
        // SDOF page 
        Table.rows[0].style.display  = "table-row";  // Peak
        Table.rows[1].style.display  = "table-row";  // Mean
        Table.rows[2].style.display  = "table-row";  // RMS
        Table.rows[3].style.display  = "none";       // Residual
        Table.rows[4].style.display  = "table-row";  // Analysis Method
        Table.rows[5].style.display  = "table-row";  // SDOF Display
        Table.rows[6].style.display  = "none";       // ResSpec Display
        Table.rows[7].style.display  = "none";       // Spectrum Display
        Table.rows[8].style.display  = "none";       // SM_Par Display
        Table.rows[9].style.display  = "table-row";  // Graph Unit
        Table.rows[10].style.display = "none";       // Baseline 
        Table.rows[11].style.display = "table-row";  // Filter 

        // Two toggels on the Info Summary table for (FFT) and (Filter) checkboxes
        document.getElementById('Filter_Div_ID1_'+ChannelList[ChNum].Unique_ID).style.display = 'none';
        document.getElementById('Filter_Div_ID2_'+ChannelList[ChNum].Unique_ID).style.display = 'flex';

    }
    else if (PageNo == 4) {
        // ResSpectrum page 
        Table.rows[0].style.display  = "none";       // Peak
        Table.rows[1].style.display  = "none";       // Mean
        Table.rows[2].style.display  = "none";       // RMS
        Table.rows[3].style.display  = "none";       // Residual
        Table.rows[4].style.display  = "table-row";  // Analysis Method
        Table.rows[5].style.display  = "none";       // SDOF Display
        Table.rows[6].style.display  = "table-row";  // ResSpec Display
        Table.rows[7].style.display  = "none";       // Spectrum Display
        Table.rows[8].style.display  = "none";       // SM_Par Display
        Table.rows[9].style.display  = "table-row";  // Graph Unit
        Table.rows[10].style.display = "table-row";  // Baseline 
        Table.rows[11].style.display = "table-row";  // Filter 
        
        // Two toggels on the Info Summary table for (FFT) and (Filter) checkboxes
        document.getElementById('Filter_Div_ID1_'+ChannelList[ChNum].Unique_ID).style.display = 'flex';
        document.getElementById('Filter_Div_ID2_'+ChannelList[ChNum].Unique_ID).style.display = 'none';
    }
    else if (PageNo == 5) {
        // Spectrum page 
        Table.rows[0].style.display  = "none";       // Peak
        Table.rows[1].style.display  = "none";       // Mean
        Table.rows[2].style.display  = "none";       // RMS
        Table.rows[3].style.display  = "none";       // Residual
        Table.rows[4].style.display  = "none";       // Analysis Method
        Table.rows[5].style.display  = "none";       // SDOF Display
        Table.rows[6].style.display  = "none";       // ResSpec Display
        Table.rows[7].style.display  = "table-row";  // Spectrum Display
        Table.rows[8].style.display  = "none";       // SM_Par Display
        Table.rows[9].style.display  = "table-row";  // Graph Unit
        Table.rows[10].style.display = "table-row";  // Baseline 
        Table.rows[11].style.display = "table-row";  // Filter 
        
        // Two toggels on the Info Summary table for (FFT) and (Filter) checkboxes
        document.getElementById('Filter_Div_ID1_'+ChannelList[ChNum].Unique_ID).style.display = 'none';
        document.getElementById('Filter_Div_ID2_'+ChannelList[ChNum].Unique_ID).style.display = 'none';

    }
    else if (PageNo == 6) {
        // SM Parameters page 
        Table.rows[0].style.display  = "none";       // Peak
        Table.rows[1].style.display  = "none";       // Mean
        Table.rows[2].style.display  = "none";       // RMS
        Table.rows[3].style.display  = "table-row";  // Significant Duration or Bracketed Duration (Residual)
        Table.rows[4].style.display  = "none";       // Analysis Method
        Table.rows[5].style.display  = "none";       // SDOF Display
        Table.rows[6].style.display  = "none";       // ResSpec Display
        Table.rows[7].style.display  = "none";       // Spectrum Display
        Table.rows[8].style.display  = "table-row";  // SM_Par Display
        Table.rows[9].style.display  = "none";       // Graph Unit
        Table.rows[10].style.display = "table-row";  // Baseline 
        Table.rows[11].style.display = "table-row";  // Filter 
        
        // Two toggels on the Info Summary table for (FFT) and (Filter) checkboxes
        document.getElementById('Filter_Div_ID1_'+ChannelList[ChNum].Unique_ID).style.display = 'none';
        document.getElementById('Filter_Div_ID2_'+ChannelList[ChNum].Unique_ID).style.display = 'none';
        
    }
    else if (PageNo == 7) {
        // Drift page 
        Table.rows[0].style.display  = "table-row";     // Peak
        Table.rows[1].style.display  = "table-row";     // Mean
        Table.rows[2].style.display  = "table-row";     // RMS
        Table.rows[3].style.display  = "table-row";     // Residual
        Table.rows[4].style.display  = "none";          // Analysis Method
        Table.rows[5].style.display  = "none";          // SDOF Display
        Table.rows[6].style.display  = "none";          // ResSpec Display
        Table.rows[7].style.display  = "none";          // Spectrum Display
        Table.rows[8].style.display  = "none";          // SM_Par Display
        Table.rows[9].style.display  = "table-row";     // Graph Unit
        Table.rows[10].style.display = "table-row";     // Baseline 
        Table.rows[11].style.display = "table-row";     // Filter 
        
        // Two toggels on the Info Summary table for (FFT) and (Filter) checkboxes
        document.getElementById('Filter_Div_ID1_'+ChannelList[ChNum].Unique_ID).style.display = 'flex';
        document.getElementById('Filter_Div_ID2_'+ChannelList[ChNum].Unique_ID).style.display = 'flex';
        
    }
    else if (PageNo == 8) {
        // H/V Spectral Ratio
        Table.rows[0].style.display  = "none";       // Peak
        Table.rows[1].style.display  = "none";       // Mean
        Table.rows[2].style.display  = "none";       // RMS
        Table.rows[3].style.display  = "none";       // Residual
        Table.rows[4].style.display  = "none";       // Analysis Method
        Table.rows[5].style.display  = "none";       // SDOF Display
        Table.rows[6].style.display  = "none";       // ResSpec Display
        Table.rows[7].style.display  = "none";       // Spectrum Display
        Table.rows[8].style.display  = "none";       // SM_Par Display
        Table.rows[9].style.display  = "none";       // Graph Unit
        Table.rows[10].style.display = "table-row";  // Baseline 
        Table.rows[11].style.display = "table-row";  // Filter 
        
        // Two toggels on the Info Summary table for FFT and Filter 
        document.getElementById('Filter_Div_ID1_'+ChannelList[ChNum].Unique_ID).style.display = 'none';
        document.getElementById('Filter_Div_ID2_'+ChannelList[ChNum].Unique_ID).style.display = 'none';
    }

    // Update the graph
    Plotly.update(PlotArea_ID, traces, layout_update);

}
