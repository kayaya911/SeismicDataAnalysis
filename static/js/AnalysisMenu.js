// stricter parsing and error handling
"use strict";

async function AnalysisMenu_Selection(a) {

    let i;

    if (a.id == "MainMenu_LoadData") {
      // Show-Hide Parameters Windows
      document.getElementById("Table_Channel_Div").style.display     = "flex";
      document.getElementById("Parameters_Filter").style.display     = "none";
      document.getElementById("Parameters_Integration").style.display= "none";
      document.getElementById("Parameters_SDOF").style.display       = "none";
      document.getElementById("Parameters_ResSpec").style.display    = "none";
      document.getElementById("Parameters_Newmark_Settings").style.display       = "none";
      document.getElementById("Parameters_Output_Unit_Settings").style.display   = "none";

      // Show-Hide Results-Containers
      document.getElementById("Graphs_LoadData_Container").style.display    = "flex";
      document.getElementById("Graphs_SDOF_Container").style.display        = "none";
      document.getElementById("Graphs_ResSpec_Container").style.display     = "none";

      // Show-Hide Calculation Button
      document.getElementById("Calculate_Div").style.display     = "none";

      // Update PageNo
      PageNo = 0;
      
      // Update all Graphs
      for (i=0; i<ChannelList.length; i++) { await Plot_Update(i); }

      // Adjust the height of the Parameter_Container
      document.getElementById("Parameter_Container").style.height    = "30%";
      document.getElementById("Parameter_Container").style.minHeight = "30%";
      document.getElementById("Parameter_Container").style.maxHeight = "fit-content";

    }
    else if (a.id == "MainMenu_Filter") {
      // Show-Hide Parameters Windows
      document.getElementById("Table_Channel_Div").style.display     = "none";
      document.getElementById("Parameters_Filter").style.display     = "flex";
      document.getElementById("Parameters_Integration").style.display= "none";
      document.getElementById("Parameters_SDOF").style.display       = "none";
      document.getElementById("Parameters_ResSpec").style.display    = "none";
      document.getElementById("Parameters_Newmark_Settings").style.display       = "none";
      document.getElementById("Parameters_Output_Unit_Settings").style.display   = "none";

      // Show-Hide Results-Containers
      document.getElementById("Graphs_LoadData_Container").style.display    = "flex";
      document.getElementById("Graphs_SDOF_Container").style.display        = "none";
      document.getElementById("Graphs_ResSpec_Container").style.display     = "none";

      // Show-Hide Calculation Button
      document.getElementById("Calculate_Div").style.display     = "flex";

      // innerHTML of Button
      document.getElementById("Calculate_Button").innerHTML = "Filter";

      // Update PageNo
      PageNo = 1;

      // Update all Graphs
      for (i=0; i<ChannelList.length; i++) { await Plot_Update(i); }

      // Adjust the height of the Parameter_Container
      document.getElementById("Parameter_Container").style.height    = "30%";
      document.getElementById("Parameter_Container").style.minHeight = "fit-content";
      document.getElementById("Parameter_Container").style.maxHeight = "fit-content";
      
    }
    else if (a.id == "MainMenu_Integral") {
      // Show-Hide Parameters Windows
      document.getElementById("Table_Channel_Div").style.display     = "none";
      document.getElementById("Parameters_Filter").style.display     = "flex";
      document.getElementById("Parameters_Integration").style.display= "flex";
      document.getElementById("Parameters_SDOF").style.display       = "none";
      document.getElementById("Parameters_ResSpec").style.display    = "none";
      document.getElementById("Parameters_Newmark_Settings").style.display       = "none";
      document.getElementById("Parameters_Output_Unit_Settings").style.display   = "none";

      // Show-Hide Results-Containers
      document.getElementById("Graphs_LoadData_Container").style.display    = "flex";
      document.getElementById("Graphs_SDOF_Container").style.display        = "none";
      document.getElementById("Graphs_ResSpec_Container").style.display     = "none";

      // Show-Hide Calculation Button
      document.getElementById("Calculate_Div").style.display     = "flex";

      // innerHTML of Button
      document.getElementById("Calculate_Button").innerHTML = "Integrate";

      // Update PageNo
      PageNo = 2;

      // Update all Graphs
      for (i=0; i<ChannelList.length; i++) { await Plot_Update(i); }

      // Adjust the height of the Parameter_Container
      document.getElementById("Parameter_Container").style.height    = "30%";
      document.getElementById("Parameter_Container").style.minHeight = "fit-content";
      document.getElementById("Parameter_Container").style.maxHeight = "fit-content";
    }
    else if (a.id == "MainMenu_SDOF") {
      // Show-Hide Parameters Windows
      document.getElementById("Table_Channel_Div").style.display     = "none";
      document.getElementById("Parameters_Filter").style.display     = "flex";
      document.getElementById("Parameters_Integration").style.display= "none";
      document.getElementById("Parameters_SDOF").style.display       = "flex";
      document.getElementById("Parameters_ResSpec").style.display    = "none";
      document.getElementById("Parameters_Newmark_Settings").style.display       = "none";
      document.getElementById("Parameters_Output_Unit_Settings").style.display   = "none";
      
      // Show-Hide Results-Containers
      document.getElementById("Graphs_LoadData_Container").style.display    = "none";
      document.getElementById("Graphs_SDOF_Container").style.display        = "flex";
      document.getElementById("Graphs_ResSpec_Container").style.display     = "none";

      // Show-Hide Calculation Button
      document.getElementById("Calculate_Div").style.display     = "flex";

      // innerHTML of Button
      document.getElementById("Calculate_Button").innerHTML = "SDOF Response";

      // Update PageNo
      PageNo = 3;

      // Adjust the height of the Parameter_Container
      document.getElementById("Parameter_Container").style.height    = "30%";
      document.getElementById("Parameter_Container").style.minHeight = "fit-content";
      document.getElementById("Parameter_Container").style.maxHeight = "fit-content";

    }
    else if (a.id == "MainMenu_ResSpec") {
      // Show-Hide Parameters Windows
      document.getElementById("Table_Channel_Div").style.display                 = "none";
      document.getElementById("Parameters_Filter").style.display                 = "flex";
      document.getElementById("Parameters_SDOF").style.display                   = "none";
      document.getElementById("Parameters_Integration").style.display            = "none";
      document.getElementById("Parameters_ResSpec").style.display                = "flex";
      document.getElementById("Parameters_Newmark_Settings").style.display       = "none";
      document.getElementById("Parameters_Output_Unit_Settings").style.display   = "none";

      // Show-Hide Results-Containers
      document.getElementById("Graphs_LoadData_Container").style.display    = "none";
      document.getElementById("Graphs_SDOF_Container").style.display        = "none";
      document.getElementById("Graphs_ResSpec_Container").style.display     = "flex";

      // Show-Hide Calculation Button
      document.getElementById("Calculate_Div").style.display     = "flex";

      // innerHTML of Button
      document.getElementById("Calculate_Button").innerHTML = "Response Spectrum";

      // Update PageNo
      PageNo = 4;

      // Adjust the height of the Parameter_Container
      document.getElementById("Parameter_Container").style.height    = "30%";
      document.getElementById("Parameter_Container").style.minHeight = "fit-content";
      document.getElementById("Parameter_Container").style.maxHeight = "fit-content";
    }
    else if ((a.id == "MainMenu_Settings") || (a.id == "Header_Settings")) {
      // Show-Hide Parameters Windows
      document.getElementById("Table_Channel_Div").style.display                 = "none";
      document.getElementById("Parameters_Filter").style.display                 = "none";
      document.getElementById("Parameters_Integration").style.display            = "none";
      document.getElementById("Parameters_SDOF").style.display                   = "none";
      document.getElementById("Parameters_ResSpec").style.display                = "none";
      document.getElementById("Parameters_Newmark_Settings").style.display       = "flex";
      document.getElementById("Parameters_Output_Unit_Settings").style.display   = "flex";

      // Show-Hide Results-Containers
      document.getElementById("Graphs_LoadData_Container").style.display    = "none";
      document.getElementById("Graphs_SDOF_Container").style.display        = "none";
      document.getElementById("Graphs_ResSpec_Container").style.display     = "none";

      // Show-Hide Calculation Button
      document.getElementById("Calculate_Div").style.display     = "none";

      // Update PageNo
      PageNo = 10;

      // Adjust the height of the Parameter_Container
      document.getElementById("Parameter_Container").style.height    = "30%";
      document.getElementById("Parameter_Container").style.minHeight = "fit-content";
      document.getElementById("Parameter_Container").style.maxHeight = "fit-content";

    }

    // Hide the Analysis_Menu on the screen
    document.getElementById("Analysis_Menu").style.display = "none";
    
}
function ShowHide_AnalysisMenu() {
    var x = document.getElementById("Analysis_Menu");
    if (x.style.display === "none") {
      x.style.display = "flex";
    } else {
      x.style.display = "none";
    }
}
