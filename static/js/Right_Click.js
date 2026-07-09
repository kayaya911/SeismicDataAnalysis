// stricter parsing and error handling
"use strict";

//------------------------------------------------------------------------------------------------------------------------------------------------
function Right_Click_ALL() {

    // Decleration of all variables 
    let TreeView, Panel1, Right_Click_Div, aa1, aa2, pp1, pp2;
    
    TreeView        = document.getElementById('sidebar');
    Panel1          = document.getElementById('panel1');
    Right_Click_Div = document.getElementById('Right_Click_Div');

    // Show context menu on right-click for sidebar
    TreeView.addEventListener('contextmenu', (e) => {
        e.preventDefault(); // Prevent default context menu
        Right_Click_Div.style.display = 'block';

        // Calculate the position of menu
        aa1 = e.pageY + Right_Click_Div.offsetHeight;
        aa2 = TreeView.offsetHeight+ TreeView.offsetTop;
        pp1 = e.pageX;
        if ( aa1 > aa2 ) { pp2 = aa2 - Right_Click_Div.offsetHeight - 10; } else { pp2 = e.pageY; }

        Right_Click_Div.style.left = `${pp1}px`;
        Right_Click_Div.style.top = `${pp2}px`;
    });

    // Show context menu on right-click for panel-1
    Panel1.addEventListener('contextmenu', (e) => {
        e.preventDefault(); // Prevent default context menu
        Right_Click_Div.style.display = 'block';

        // Calculate the position of menu
        aa1 = e.pageY + Right_Click_Div.offsetHeight;
        aa2 = TreeView.offsetHeight+ TreeView.offsetTop;
        pp1 = e.pageX;
        if ( aa1 > aa2 ) { pp2 = aa2 - Right_Click_Div.offsetHeight - 10; } else { pp2 = e.pageY; }

        Right_Click_Div.style.left = `${pp1}px`;
        Right_Click_Div.style.top = `${pp2}px`;
    });

    // Hide context menu when clicking elsewhere
    document.addEventListener('click', (e) => {
        if (!Right_Click_Div.contains(e.target)) {
            Right_Click_Div.style.display = 'none';
        }
    });

    // Hide on scroll
    document.addEventListener('scroll', () => {
        if (Right_Click_Div) {
            Right_Click_Div.style.display = 'none';
        }
    }, true);

    // Add <li> elements to Righ-Click Menu
    Add_li("Right_Click_ul", 'Right_Click_Select_All_Analysis',          'Select All Channels for Analysis', true);
    Add_li("Right_Click_ul", 'Right_Click_Select_All_Plotting',          'Select First ' + MaxPlotly_Graphs.toString() + ' Channels for Plotting', true);
    Add_li("Right_Click_ul", 'Right_Click_Select_All_Collapse_Sidebar',  'Expand Files in Sidebar', true);
    document.getElementById('Right_Click_ul').appendChild( document.createElement("br")  );
    
    Add_li("Right_Click_ul", 'Right_Click_Select_All_Collapse_Tables',   'Expand Files in Tables', true);
    Add_li("Right_Click_ul", 'Right_Click_Select_All_Single_Table',      'Switch Table View', true);
    
}
//------------------------------------------------------------------------------------------------------------------------------------------------
function Right_Click(el) {

    // Decleration of variables
    let checkboxes

    if (el.id === 'Right_Click_Select_All_Analysis') {
        Toggle_Sidebar_SelectAllChannels_For_Analysis();

    }
    else if (el.id === 'Right_Click_Select_All_Plotting') {
        Toggle_SidebBar_SelectAllChannels_For_Plotting();

    }
    else if (el.id === 'Right_Click_Select_All_Collapse_Sidebar') {
        Toggel_Sidebar_Collapse_Expand_Files();

    }
    else if (el.id === 'Right_Click_Select_All_Collapse_Tables') {
        Toggle_Tabels_Collapse_Files();

    }
    else if (el.id === 'Right_Click_Select_All_Single_Table') {
        Toggle_Tables_SingleTable();

    } else if (el.id.includes('FSamp_')) {

        // Decleration of variables 
        let i, status, result = [];

        // Get a string array of Sampling-Frequencies selected by user
        result = Get_UL_List_Checked('Right_Click_ul_Sampling', 'FSamp_');

        // Go through all channels 
        for (i=0; i<ChannelList.length; i++) { 
            if (result.includes(ChannelList[i].FSamp.toFixed(2))) { status = true;  } else { status = false; }
            ChannelList[i].Selected = status; 
            document.getElementById( 'FileTreeView_Checkbox_' + ChannelList[i].Unique_ID ).checked = status;
        }

        // Uncheck others
        Array.from(document.getElementById('Right_Click_ul_Type').querySelectorAll('li input[type="checkbox"]')).forEach( item => { item.checked = false; });
        Array.from(document.getElementById('Right_Click_ul_Azimuth').querySelectorAll('li input[type="checkbox"]')).forEach( item => { item.checked = false; });

        // Uncheck id='Right_Click_Select_All_Analysis'
        if (result.length != 0) { document.getElementById('Right_Click_Select_All_Analysis').checked = false; }
        

    } else if (el.id.includes('Type_')) {

        // Decleration of variables 
        let i, status, result = [];

        // Get a string array of Channel-Types selected by user
        result = Get_UL_List_Checked('Right_Click_ul_Type', 'Type_');

        // Go through all channels 
        for (i=0; i<ChannelList.length; i++) { 
            if (result.includes(ChannelList[i].TypeString)) { status = true;  } else { status = false; }
            ChannelList[i].Selected = status; 
            document.getElementById( 'FileTreeView_Checkbox_' + ChannelList[i].Unique_ID ).checked = status;
        }

        // Uncheck others
        Array.from(document.getElementById('Right_Click_ul_Sampling').querySelectorAll('li input[type="checkbox"]')).forEach( item => { item.checked = false; });
        Array.from(document.getElementById('Right_Click_ul_Azimuth').querySelectorAll('li input[type="checkbox"]')).forEach( item => { item.checked = false; });

        // Uncheck id='Right_Click_Select_All_Analysis'
        if (result.length != 0) { document.getElementById('Right_Click_Select_All_Analysis').checked = false; }

    } else if (el.id.includes('Azimuth_')) {

        // Decleration of variables 
        let i, status, result = [];

        // Get a string array of channel-orientations selected by user
        result = Get_UL_List_Checked('Right_Click_ul_Azimuth', 'Azimuth_');

        // Go through all channels 
        for (i=0; i<ChannelList.length; i++) { 
            if (result.includes(ChannelList[i].Orientation)) { status = true;  } else { status = false; }
            ChannelList[i].Selected = status; 
            document.getElementById( 'FileTreeView_Checkbox_' + ChannelList[i].Unique_ID ).checked = status;
        }

        // Uncheck others
        Array.from(document.getElementById('Right_Click_ul_Sampling').querySelectorAll('li input[type="checkbox"]')).forEach( item => { item.checked = false; });
        Array.from(document.getElementById('Right_Click_ul_Type').querySelectorAll('li input[type="checkbox"]')).forEach( item => { item.checked = false; });
        
        // Uncheck id='Right_Click_Select_All_Analysis'
        if (result.length != 0) { document.getElementById('Right_Click_Select_All_Analysis').checked = false; }
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------------
function Add_li(Parent_ID, Input_ID, TEXT_Content, status) {

    // <li> element with checkbox and label is created and appended it to Parrent_ID element

    // Decleration of variables 
    let li_el, input, label

    Parent_ID = document.getElementById(Parent_ID);

    li_el = document.createElement('li'    );
    input = document.createElement('input' );
    label = document.createElement('label' );

    input.setAttribute('class',      "form-check-input"  );
    input.setAttribute('type',       "checkbox"          );
    input.setAttribute('id',         Input_ID           );
    input.setAttribute('onclick',    "Right_Click(this)");
    input.checked = status;

    label.textContent = TEXT_Content;
    label.setAttribute('for',         Input_ID            );
    label.setAttribute('class',       'RightClick_Label'  );

    li_el.setAttribute('class',       "Right_Click_li"); 
    li_el.appendChild(input);
    li_el.appendChild(label);
    
    Parent_ID.appendChild(li_el);
}
//------------------------------------------------------------------------------------------------------------------------------------------------
function Get_UL_List(UL_ID, Opt) {

    // Returns a string array of unique sampling rate of all channels in ChannelList[]

    // Decleration of variables
    let Array_li, checkbox, result=[]
    
    // Get an array of <li> elements
    Array_li = Array.from( document.getElementById(UL_ID).querySelectorAll('li') );

    // For each <li> element in the array
    Array_li.forEach( li => {
        // Get the Checkbox element
        checkbox = li.querySelector('input[type="checkbox"]');

        if (Opt==1) {
            result.push(Number(checkbox.id.replace('FSamp_','')).toFixed(2));

        } else if (Opt == 2) {
            result.push(checkbox.id.replace('Type_',''));

        } else if (Opt == 3) {
            result.push(checkbox.id.replace('Azimuth_',''));

        }
        
    });
    // Retrun unique elements of result array
    return [...new Set(result)];
}
//------------------------------------------------------------------------------------------------------------------------------------------------
function Get_UL_List_Checked(UL_ID, Opt) {
    // Return a string array of one of the following attributes of Channel selected by user
    // Sampling Rate, Channle Type, or Azimuth
    let CheckedBoxes, result=[]
    CheckedBoxes = Array.from(document.getElementById(UL_ID).querySelectorAll('li input[type="checkbox"]')).filter(item => item.checked);
    CheckedBoxes.forEach( cb => {
        result.push(cb.id.replace(Opt,''));
    });
    return result;
}
//------------------------------------------------------------------------------------------------------------------------------------------------
function Order_Plotly_Graphs() {
    // This function re-orders the plotly graphs in the Panel2 as the order selected by the user in the sidebar

    // Decleration of variables 
    let CheckedBoxes, DivElements, DivID

    // All div-elements of Plotly graphs
    DivElements = Array.from( document.querySelectorAll('div[id*="Div_ID_"]') );

    // Delete the existing order of ploty graphs
    DivElements.forEach( item => { item.style.order = '';  });

    // All Plotly checkboxes that are user-selected
    CheckedBoxes = Array.from( document.querySelectorAll('input[type="checkbox"][class*="File_CheckBox_Ch_Green"]') ).filter(item => item.checked);

    // Re-order the plotly graphs
    CheckedBoxes.forEach( (item, i) => {
        DivID = item.id.replace("PlotChecbox_", "Div_ID_");
        document.getElementById(DivID).style.order = i; 
    });
}
//------------------------------------------------------------------------------------------------------------------------------------------------
function Get_Plotly_Graphs_Invisable() {
    // Returns an array of div elements of Plotly-Grapgh that are hidden (not visable)
    let AllDivs = Array.from(document.querySelectorAll('div[id*="Div_ID_"]'));
    return AllDivs.filter(div => window.getComputedStyle(div).display === 'none');
}
//------------------------------------------------------------------------------------------------------------------------------------------------