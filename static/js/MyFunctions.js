// stricter parsing and error handling
"use strict";

// Declaration of global variable
let ChannelList       = [];
let PageNo            = 0;
let Number_of_Graph   = 0;
let Max_Num_of_Graph  = 0;


//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
function OnLoad() {
    // Wecome message 
    console.log("Welcome to Website for Seismic Data Analysis...");
    
    //Hide the AnalysisMenu at start-up
    ShowHide_AnalysisMenu();

    // Add event listener to a button for uploading files
    document.getElementById('LoadInputFiles').addEventListener('change', function(e) {
        Load_Files(e); // Allows user to select files and uploads them
        AnalysisMenu_Selection(document.getElementById("MainMenu_LoadData")); // Take the screen to Load_Data page
    });

    // Close the Analysis_Menu on mouse-clicks anywhere outside of the Analysis_Menu
    document.body.addEventListener("click", () => { document.getElementById("Analysis_Menu").style.display = "none"; 
                                                    document.getElementById("Right_Click_Div").style.display = "none"; })
    document.getElementById("Analysis_Menu").addEventListener("click", (ev) => { ev.stopPropagation(); })
    document.getElementById("Right_Click_Div").addEventListener("click", (ev) => { ev.stopPropagation(); })
    document.getElementById("Main_Menu_Button").addEventListener("click", (ev) => { ev.stopPropagation(); })

    // Load_Data item from Analysis_Menu is selected 
    AnalysisMenu_Selection(document.getElementById("MainMenu_LoadData"));

    // Right Click Function
    Right_Click_ALL();

    // Filter Seetings 
    FilterName_Change();
    FilterType_Change();

}




