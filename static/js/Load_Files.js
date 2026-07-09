// await is a keyword in JavaScript that pauses execution of an async function until a Promise is settled

// stricter parsing and error handling
"use strict";

//-------------------------------------------------------------------------------------------------------------
class Channel {
    constructor() {

        let Unique_ID = Generate_Unique_ID();

        this.FileName             = undefined;         // FileName
        this.ChNum                = undefined;         // Channel number
        this.Unique_ID            = Unique_ID;         // Unique ID of CheckBox
        this.TableName            = undefined;         // TableName where the channel information will be listed. It is the name of the file where "-" and "." in the FileName are replaced by "_"
        this.FileListName         = undefined;         // FileListName where the channel information will be listed. It is the name of the file where "-" and "." in the FileName are replaced by "_" also _FileList indicator
        this.PlotGraph            = false;             // Whether this channel is selected by the user for plotting or not (true / false)
        this.Selected             = true;              // Whether this channel is selected by the user or not (true / false)
        this.NumSamples           = undefined;         // Number of digitized data points
        this.ScaleFactor          = undefined;         // User Scale Factor - Taken 1.0 if not specified
        this.Orientation          = undefined;         // Orientation
        this.DateTime             = undefined;         // Date&Time of the first sample in the record
        this.DateTime_End         = undefined;         // Date&Time of the last sample in the record
        this.Duration             = undefined;         // Total duration in seconds
        this.Lat                  = undefined;         // Latitude coordinate of the sensor
        this.Long                 = undefined;         // Longitude coordinate of the sensor
        this.FSamp                = undefined;         // Sampling Frequency in Hz
        this.delt                 = undefined;         // Time interval of the digitized record in seconds
        this.TypeAndUnits         = undefined;         // Type And Unit number - refer to the list
        this.Type                 = undefined;         // Type of the channel (0 - 9)
        this.TypeString           = undefined;         // String version of the type of the channel (e.g., acceleration, velocity, etc.)
        this.Unit                 = undefined;         // Unit of the channel (0 - 23)
        this.UnitString           = undefined;         // String version of the unit (e.g., g, m/s, etc.)
        this.IntervalTypeAndUnit  = undefined;         // Type And Unit number - refer to the list
        this.IntervalType         = undefined;         // (0-2)
        this.IntervalTypeString   = undefined;         // (Time, Period, etc.)
        this.IntervalUnit         = undefined;         // (0-4)
        this.IntervalUnitString   = undefined;         // (Second, UTC, etc.)
        this.RMS                  = undefined;         // RMS value of digitized raw data
        this.Peak                 = undefined;         // Peak value of digitized raw data
        this.Mean                 = undefined;         // Mean value of digitized raw data
        this.Residual             = undefined;         // Residual of digitized raw data  - last value of digitized raw data
        this.InstFreq             = undefined;         // Natural frequency of transducer (in Hz).
        this.InstPeriod           = undefined;         // Natural period of transducer (in seconds).
        this.InstDamp             = undefined;         // Damping ratio of transducer (fraction of critical).
        this.data                 = undefined;         // Digitized data
        this.time                 = undefined;         // Time array of the digitized data - starts from time zero.
        
        this.Results = { 
            Filter: [],
            Integral: [],
            SDOF: [],
            ResSpec: [],
            Spectrum: [],
            SM_Parameters: [],
            Drift: [],
            HVSR: [],
            P_Wave: [],
        };

        // Delete Results for this channes 
        this.DeleteResults = function() { 
            this.Results.Filter = [];
            this.Results.Integral = [];
            this.Results.SDOF = [];
            this.Results.ResSpec = [];
            this.Results.Spectrum = [];
            this.Results.SM_Parameters = [];
            this.Results.HVSR = [];
            this.Results.Drift = [];
            this.Results.P_Wave = [];
        }

        // Get data 
        this.GetData = function() {
            return Multiply(this.data, this.ScaleFactor);
        }

    }
}
//-------------------------------------------------------------------------------------------------------------
function IntervalTypeAndUnit(n, Opt) {
    // The supported type of Sampling-Interval and its units are given below
    // 
    //  n     ==>     Type                          Units
    //  ----------------------------------------------------------------
    //  1     ==>     0   ( Time )                  0 (Second)  [starts from zero and increases with constant delt intervals]
    //  2     ==>     0   ( Time )                  1 ('')      [contains Date And Time information]
    //  3     ==>     1   ( Period )                2 (Second)  [contains spectral periods in unit of seconds]
    
    if (Opt == null) { Opt = true; }

    if (Opt) {
        // TypeAndUnit number (n) is given
        if      (n == 1)   { return { 'IntervalTypeAndUnit': 1,   'Type': 0,   'Type_String': 'Time',   'Unit':  0,   'Unit_String': 's'   }; }   //  ( Time )   ( sec. )
        else if (n == 2)   { return { 'IntervalTypeAndUnit': 2,   'Type': 0,   'Type_String': 'Time',   'Unit':  1,   'Unit_String': ''    }; }   //  ( Time )   ( UTC  )
        else if (n == 3)   { return { 'IntervalTypeAndUnit': 3,   'Type': 1,   'Type_String': 'Period', 'Unit':  2,   'Unit_String': 's'   }; }   //  ( Period ) ( sec. )

    }
    else {
        // Unit number is given.
        if      (n == 0)   { return { 'IntervalTypeAndUnit': 1,   'Type': 0,   'Type_String': 'Time',   'Unit':  0,   'Unit_String': 's'   }; }   //  ( Time )   ( sec. )
        else if (n == 1)   { return { 'IntervalTypeAndUnit': 2,   'Type': 0,   'Type_String': 'Time',   'Unit':  1,   'Unit_String': ''    }; }   //  ( Time )   ( UTC  )
        else if (n == 2)   { return { 'IntervalTypeAndUnit': 3,   'Type': 1,   'Type_String': 'Period', 'Unit':  2,   'Unit_String': 's'   }; }   //  ( Period ) ( sec. )
    }

}
// Delete Records ---------------------------------------------------------------------------------------------
function DeleteRecordings() {
    
    if ( confirm('\n\nDELETE ALL RECORDINGS ?') ) {

        //Decleration of variables 
        let i, row, table, numOfRows;

        // Delete all channels 
        ChannelList = [];

        // Reset the Progress Bar
        ResetProgressBar("ProgressBar_LoadData");

        // Reset Graphs
        document.getElementById("Graphs_LoadData_Container").innerHTML = ""; 
        document.getElementById("FileListTree").innerHTML = "";

        // Clear RecordsTable
        table = document.getElementById("Table_Channels");
        numOfRows = table.rows.length;
        for (i=1; i<numOfRows-1; i++) { table.deleteRow(1); }        // delete all rows except the header and the last row
        row = table.insertRow(-1);                                   // add one empty raw to the end of the list
        row.setAttribute('class', 'Channels_Table_Body_Tr');
        table.deleteRow(1);                                          // delete the first row from the table

        // Clear InfoBar

        // Switch to the LoadData Page
        AnalysisMenu_Selection({id:"MainMenu_LoadData"});
        
    }
}
//-------------------------------------------------------------------------------------------------------------
async function Load_Files(ev) {
    // Appends the user-selected files as a Channel Object to the ChannelList[]
    // Files are read based on their file extension (e.g., VIF, V1, RAW, ASC, MSD, etc.)

    let array_files = ev.target.files;
    let i, len, delta, Indexs;

    len = array_files.length;

    if (len > 0) { 
        
        // Reset the ProgressBar
        document.getElementById("ProgressBar_Label").dataset.message = '0';
        await sleep(5);

        Indexs = LinStep(0, len-1, 1);
        delta  = 100 / len;
        
        for (i of Indexs) {

            // Skips this file if it is already uploaded.
            if (IsFileUploaded( array_files[i].name )) { 

                UpdateProgress(delta, 'ProgressBar_Label');
                await sleep(5);
                continue; 

            }

            // Reads the content of each file.
            await ReadFileContent( array_files[i] )

            // Waits for 5 miliseccond for secreen update
            await sleep(5);
        }
    }

    // Read the content of the file
    async function ReadFileContent(File_Blob) {

        // decleration of variables 
        let dataview, FileName, fileExt;

        // Convert Blob object to an ArrayBuffer of Promise (1 byte = 8 bits)
        File_Blob.arrayBuffer().then((result) => {
    
            // Create a dataView object to read the content of the ArrayBuffer
            dataview = new DataView(result);
            FileName = File_Blob.name;
    
            // Get the extension of the file
            fileExt = FileExtension( FileName ).toUpperCase();
            
            // Appends each channel in the file into the ChannelList[]
            // Also, appends each channel to the Table on the home-page
            if      (  fileExt == "VIF" )                                                   { Read_VIF_BCSIMS( FileName, delta, dataview ); }
            else if (( fileExt == "V1"  ) || (fileExt == "RAW"))                            { Read_V1_COSMOS(  FileName, delta, dataview ); }
            else if  ( fileExt == "V1C" )                                                   { Read_V1c_COSMOS( FileName, delta, dataview ); }
            else if (( fileExt == "MSD" ) || (fileExt == "MSEED"))                          { Read_MSD(        FileName, delta, dataview ); }
            else if (  fileExt == "TXT" )                                                   { Read_TXT(        FileName, delta, dataview ); }
            else if (  fileExt == "ASC" )                                                   { Read_ASC(        FileName, delta, dataview ); }
            else if (( fileExt == "AT2" ) || ( fileExt == "VT2" ) || ( fileExt == "DT2" ) ) { Read_PEER_Data(   FileName, delta, dataview); }
            else if (  fileExt == "DAT" ) {
                
                // Check if it matches the '=>' sync character from your existing Read_DAT
                let isType1 = dataview.getUint16(66, true) === 0x3E3D;

                if (isType1) { Read_DAT_GEOSIG(FileName, delta, dataview);  } 
                else         { Read_DAT_Free(  FileName, delta, dataview);  }
            }
        });
    }

    
    
}
//-------------------------------------------------------------------------------------------------------------
async function Add_To_Table(Channel) {

    // Appends single channel as a row entry to the Table
    // Each row consists of 11 columns as follow:
    // CheckBox - FileName - ChannelNumber - Duration (sec) - Sampling (Hz) - Type - Unit - Azimuth - Peak - Mean - RMS

    // Declare variables
    let table1, table2, numOfRows, numOfCols, chk, th, td, row, cell, select, Type_ID, opt, i, j;
    let Type_List, Unit_List, Unit_ID, Scale_ID, FSamp_ID, input;
    let UL_el, li_el, label, checkbox1, checkbox2, divCount

    // Get the list of TypeList and UnitList
    Type_List = List_TypeUnit().Types;
    Unit_List = List_Units(Channel.Unit).Units;

    // Number of row and columns
    table2    = document.getElementById('Table_FileNames');

    if (!IsFileUploaded(Channel.FileName)) {
        // This FileName does not exists (first time)
        // Create a row for FileName 
        row = table2.insertRow(-1);
        row.setAttribute('class', 'File-Name-Header');
        row.setAttribute('onClick', 'Toggle_Subtable(this)');
        td = document.createElement('td');
        td.setAttribute('class', 'File-Name');
        const ss = document.createElement('span');
        ss.innerHTML = Channel.FileName
        ss.style = 'margin-left: 0.5rem;';
        td.appendChild(FileIcon_SVG('OPEN'));
        td.appendChild(ss);

        row.appendChild(td);

        // Create the second row to host the subtable 
        row = table2.insertRow(-1);
        row.setAttribute('class', 'subtable-container expanded');
        td = document.createElement('td');
        td.setAttribute('colspan', '5');
        td.setAttribute('class', 'subtable-container-td');
        
        const table      = document.createElement('table');
        const thead      = document.createElement('thead');
        const headerRow  = document.createElement('tr');
        headerRow.setAttribute('class', 'Channels_Table_Thead');
        
        const Headers   = ['FileName', 'Ch',   'Length',   'Sampling',   'Type',   'Unit',   'Scale',   'Date & Time', 'Azimuth'];
        for (i=0; i<Headers.length; i++) {
            th = document.createElement('th');
            th.innerHTML = Headers[i];
            //th.setAttribute('class', '');
            headerRow.appendChild(th);
        }

        thead.appendChild(headerRow);
        table.appendChild(thead);
        const tbody = document.createElement('tbody');
        tbody.appendChild(document.createElement('tr'));
        table.appendChild(tbody);
        table.setAttribute('id', Channel.TableName);
        table.setAttribute('class', 'Channels_Table ');
        
        td.appendChild(table);
        row.appendChild(td);

        // Create an unsorted list for this file in the FileTreeView --------------------------------------
        UL_el = document.getElementById('FileListTree');
        li_el = document.createElement('li');
        li_el.setAttribute('class', 'TreeView_Files');
        li_el.setAttribute('id', Channel.FileListName);       //  FileList_GMT_2026_03_10_17_37_11_000_vif

        label                   = document.createElement('label');
        label.textContent       = '⮞';
        label.id                = 'FileTreeView_Arrow_' + Channel.TableName.toString();
        label.setAttribute('onclick', 'FileTreeView_CollapseFile(this)');
        label.setAttribute('class', 'LabelArrow');
        li_el.appendChild(label);

        li_el.appendChild(FileIcon_SVG('OPEN'));

        label             = document.createElement('label');
        label.textContent = Channel.FileName;
        label.id          = 'FileTreeView_FileLabel_' + Channel.TableName.toString();
        label.setAttribute('class', 'TreeView_File');
        label.setAttribute('onclick', 'FileTreeView_CollapseFile(this)');
        li_el.appendChild(label);
        UL_el.appendChild(li_el);
    } 

    // Add channel to the TreeView -------------------------------------------------------------------
    UL_el    = document.getElementById('FileListTree');
    li_el    = document.createElement('li');
    li_el.setAttribute('id', 'Ch-' + Channel.ChNum + '_' +  Channel.TableName);
    li_el.setAttribute('class', 'FileListTree_Channel_li'); 

    checkbox1 = document.createElement('input');
    checkbox2 = document.createElement('input');
    label    = document.createElement('label');

    checkbox1.className  = 'form-check-input File_CheckBox_Ch_Blue';
    checkbox1.type       = 'checkbox';
    checkbox1.id         = 'FileTreeView_Checkbox_' + Channel.Unique_ID;
    checkbox1.title      = 'Click to analyze channel';
    checkbox1.checked    = false;
    checkbox1.setAttribute('onclick',    'Channel_Click(this)' );

    checkbox2.className  = 'form-check-input File_CheckBox_Ch_Green';
    checkbox2.type       = 'checkbox';
    checkbox2.id         = 'PlotChecbox_' + Channel.Unique_ID;
    checkbox2.title      = 'Click to plot graph';
    checkbox2.checked    = false;
    checkbox2.setAttribute('onclick',   'Toggle_Sidebar_Checkbox_For_PlotGraph(this)' );
    
    label.textContent = "(Ch-" + Channel.ChNum + ") (" +  Channel.Orientation  + ") (" + Number(Channel.FSamp).toFixed(3).toString() + " Hz)";
    label.setAttribute('for', 'FileTreeView_Checkbox_' + Channel.Unique_ID);
    label.setAttribute('class', 'LabelClass');
    label.setAttribute('ondblclick', 'Channel_DoubleClick(this)' );

    li_el.appendChild(checkbox1);
    li_el.appendChild(checkbox2);
    li_el.appendChild(label);

    
    UL_el.insertBefore(li_el, Array.from(UL_el.querySelectorAll('li'))[findInsertionIndex(UL_el, Channel.FileListName)] || null);


    function FileIcon_SVG(Opt) {
        // Create SVG and PATH elements
        const main = document.createElement('div');
        main.innerHTML = '<svg><path></path></svg>';

        // Select the SVG element
        let Logo_svg = main.querySelector('svg');
        Logo_svg.setAttribute('class', 'TreeView_SVG'); 
        Logo_svg.setAttribute('onclick',   'FileTreeView_Select_Unselect_Channels(this)' );
        Logo_svg.setAttribute('user-data', 'OPEN' );
        Logo_svg.setAttribute('id', 'FileTreeView_SVG_' + Channel.TableName); 
        Logo_svg.setAttribute('viewBox', '0 0 50 50');

        // Select the PATH element 
        let path_1 = main.querySelector('path');
        if (Opt.toUpperCase() == 'OPEN') {
        
            // Open folder Icon
            path_1.setAttribute('d', 'M 3 4 C 1.355469 4 0 5.355469 0 7 L 0 43.90625 C -0.0625 44.136719 -0.0390625 44.378906 0.0625 44.59375 C 0.34375 45.957031 1.5625 47 3 47 L 42 47 C 43.492188 47 44.71875 45.875 44.9375 44.4375 C 44.945313 44.375 44.964844 44.3125 44.96875 44.25 C 44.96875 44.230469 44.96875 44.207031 44.96875 44.1875 L 45 44.03125 C 45 44.019531 45 44.011719 45 44 L 49.96875 17.1875 L 50 17.09375 L 50 17 C 50 15.355469 48.644531 14 47 14 L 47 11 C 47 9.355469 45.644531 8 44 8 L 18.03125 8 C 18.035156 8.003906 18.023438 8 18 8 C 17.96875 7.976563 17.878906 7.902344 17.71875 7.71875 C 17.472656 7.4375 17.1875 6.96875 16.875 6.46875 C 16.5625 5.96875 16.226563 5.4375 15.8125 4.96875 C 15.398438 4.5 14.820313 4 14 4 Z M 3 6 L 14 6 C 13.9375 6 14.066406 6 14.3125 6.28125 C 14.558594 6.5625 14.84375 7.03125 15.15625 7.53125 C 15.46875 8.03125 15.8125 8.5625 16.21875 9.03125 C 16.625 9.5 17.179688 10 18 10 L 44 10 C 44.5625 10 45 10.4375 45 11 L 45 14 L 8 14 C 6.425781 14 5.171875 15.265625 5.0625 16.8125 L 5.03125 16.8125 L 5 17 L 2 33.1875 L 2 7 C 2 6.4375 2.4375 6 3 6 Z M 8 16 L 47 16 C 47.5625 16 48 16.4375 48 17 L 43.09375 43.53125 L 43.0625 43.59375 C 43.050781 43.632813 43.039063 43.675781 43.03125 43.71875 C 43.019531 43.757813 43.007813 43.800781 43 43.84375 C 43 43.863281 43 43.886719 43 43.90625 C 43 43.917969 43 43.925781 43 43.9375 C 42.984375 43.988281 42.976563 44.039063 42.96875 44.09375 C 42.964844 44.125 42.972656 44.15625 42.96875 44.1875 C 42.964844 44.230469 42.964844 44.269531 42.96875 44.3125 C 42.84375 44.71875 42.457031 45 42 45 L 3 45 C 2.4375 45 2 44.5625 2 44 L 6.96875 17.1875 L 7 17.09375 L 7 17 C 7 16.4375 7.4375 16 8 16 Z');
        
        } else if (Opt.toUpperCase() == 'CLOSE') {
        
            // Closed folder icon
            path_1.setAttribute('d', 'M 5 4 C 3.346 4 2 5.346 2 7 L 2 13 L 3 13 L 47 13 L 48 13 L 48 11 C 48 9.346 46.654 8 45 8 L 18.044922 8.0058594 C 17.765922 7.9048594 17.188906 6.9861875 16.878906 6.4921875 C 16.111906 5.2681875 15.317 4 14 4 L 5 4 z M 3 15 C 2.448 15 2 15.448 2 16 L 2 43 C 2 44.657 3.343 46 5 46 L 45 46 C 46.657 46 48 44.657 48 43 L 48 16 C 48 15.448 47.552 15 47 15 L 3 15 z');

        } else {
            
            // Open folder icon
            path_1.setAttribute('d', 'M 3 4 C 1.355469 4 0 5.355469 0 7 L 0 43.90625 C -0.0625 44.136719 -0.0390625 44.378906 0.0625 44.59375 C 0.34375 45.957031 1.5625 47 3 47 L 42 47 C 43.492188 47 44.71875 45.875 44.9375 44.4375 C 44.945313 44.375 44.964844 44.3125 44.96875 44.25 C 44.96875 44.230469 44.96875 44.207031 44.96875 44.1875 L 45 44.03125 C 45 44.019531 45 44.011719 45 44 L 49.96875 17.1875 L 50 17.09375 L 50 17 C 50 15.355469 48.644531 14 47 14 L 47 11 C 47 9.355469 45.644531 8 44 8 L 18.03125 8 C 18.035156 8.003906 18.023438 8 18 8 C 17.96875 7.976563 17.878906 7.902344 17.71875 7.71875 C 17.472656 7.4375 17.1875 6.96875 16.875 6.46875 C 16.5625 5.96875 16.226563 5.4375 15.8125 4.96875 C 15.398438 4.5 14.820313 4 14 4 Z M 3 6 L 14 6 C 13.9375 6 14.066406 6 14.3125 6.28125 C 14.558594 6.5625 14.84375 7.03125 15.15625 7.53125 C 15.46875 8.03125 15.8125 8.5625 16.21875 9.03125 C 16.625 9.5 17.179688 10 18 10 L 44 10 C 44.5625 10 45 10.4375 45 11 L 45 14 L 8 14 C 6.425781 14 5.171875 15.265625 5.0625 16.8125 L 5.03125 16.8125 L 5 17 L 2 33.1875 L 2 7 C 2 6.4375 2.4375 6 3 6 Z M 8 16 L 47 16 C 47.5625 16 48 16.4375 48 17 L 43.09375 43.53125 L 43.0625 43.59375 C 43.050781 43.632813 43.039063 43.675781 43.03125 43.71875 C 43.019531 43.757813 43.007813 43.800781 43 43.84375 C 43 43.863281 43 43.886719 43 43.90625 C 43 43.917969 43 43.925781 43 43.9375 C 42.984375 43.988281 42.976563 44.039063 42.96875 44.09375 C 42.964844 44.125 42.972656 44.15625 42.96875 44.1875 C 42.964844 44.230469 42.964844 44.269531 42.96875 44.3125 C 42.84375 44.71875 42.457031 45 42 45 L 3 45 C 2.4375 45 2 44.5625 2 44 L 6.96875 17.1875 L 7 17.09375 L 7 17 C 7 16.4375 7.4375 16 8 16 Z');        
        
        }

        // Append PATH element tot SVG element
        Logo_svg.appendChild(path_1);
        
        // Return the SVG element
        return Logo_svg;
    }

    function findInsertionIndex(ulEl, fileListId) {
        const items = Array.from(ulEl.querySelectorAll('li'));
        
        // Find the index of the target FileList item
        const fileListIndex = items.findIndex(li => li.id === fileListId);
        if (fileListIndex === -1) {
            // FileList item not found
            throw new Error(`List item "${fileListId}" not found in Files & Channels`);
        }

        // From that point, find the last consecutive "Ch" item
        let lastChIndex = fileListIndex;
        for (let i = fileListIndex + 1; i < items.length; i++) {
            if (items[i].id.startsWith('Ch')) {
            lastChIndex = i;
            } else {
            break; // Stop at the first non-"Ch" item
            }
        }

        // Return the index AFTER the last "Ch" item
        return lastChIndex + 1;
    }
    
    // Add channel to Table ---------------------------------------------------------------------------
    // Number of row and columns
    table1    = document.getElementById(Channel.TableName);
    numOfRows = table1.rows.length;
    numOfCols = table1.rows[0].cells.length;

    if (table1.tBodies[0].children[0].children[0] == undefined) {
        // The first row of the table-body is empty
        // Use the empty row  
        row = table1.tBodies[0].children[0];
        numOfRows--;
    } else {
        // The first row of the table-body is not empty
        // Insert a row at the end of the table-body
        row = table1.insertRow(-1);
    }
    
    // Assign class to the row
    row.setAttribute('class', 'Channels_Table_Body_Tr');
    
    // Each column of the Table-row needs to be populated
    // 1st column (File Name) - Insert an empty row first and populate it -------------------------------------
    i = 0;
    cell = row.insertCell(i);
    cell.title     = Channel.FileName;
    cell.innerHTML = Channel.FileName;


    // 2nd column (Channel Number) - Insert an empty row first and populate it -------------------------------------
    i = 1;
    cell = row.insertCell(i);
    cell.title     = Channel.ChNum;
    cell.innerHTML = Channel.ChNum;

    // 3rd column (Duration of the records) - Insert an empty row first and populate it ----------------------------
    i = 2;
    cell = row.insertCell(i);
    cell.innerHTML = Number(Channel.Duration).toFixed(3);
    cell.title     = Number(Channel.Duration).toFixed(3);

    // 4th column (Sampling Frequency in Hz) - Insert an empty row first and populate it ---------------------------
    i = 3;
    // cell = row.insertCell(i);
    // cell.innerHTML = Number(Channel.FSamp).toFixed(3);
    // cell.title     = Number(Channel.FSamp).toFixed(3);


    cell   = row.insertCell(i);
    FSamp_ID        = "FSamp_ID_" + Channel.Unique_ID;
    input           = document.createElement('input');
    input.title     = "Sampling Frequency, Hz";
    input.type      = "number";
    input.style     = "width: 100%; text-align: right"
    input.value     = Channel.FSamp;
    input.id        = FSamp_ID;
    input.setAttribute('onchange', "ChannelFSamp_Update(this)" );
    input.setAttribute('onfocus', "this.oldValue = this.value;" );
    cell.appendChild(input);




    // 5th colum (Channel Type) - Insert an empty row first and populate it ----------------------------------------
    i = 4;
    cell   = row.insertCell(i);
    Type_ID = "Type_ID_" + Channel.Unique_ID;
    select  = document.createElement('select');
    select.setAttribute("id", Type_ID);
    select.setAttribute('class', 'form-select form-select-sm  Channels_Table_Select');
    select.setAttribute("onchange", 'ChannelType_Update(this, true)');
    for (j=0; j<Type_List.length; j++) {
        opt = document.createElement("option");
        opt.value = Type_List[j];
        opt.text = Type_List[j];
        select.add(opt, null);
    }
    select.selectedIndex = Channel.Type;
    cell.title = 'Select channel type';
    cell.appendChild(select);

    // 6th colum (Channel Unit) - Insert an empty row first and populate it ----------------------------------------
    i = 5;
    cell   = row.insertCell(i);
    Unit_ID = "Unit_ID_" + Channel.Unique_ID;
    select  = document.createElement('select');
    select.setAttribute("id", Unit_ID);
    select.setAttribute('class', 'form-select form-select-sm Channels_Table_Select');
    for (j = 0; j < Unit_List.length; j++) {
        opt = document.createElement("option");
        opt.text = Unit_List[j];
        select.add(opt, null);
    }
    select.setAttribute("onchange", "ChannelType_Update(this, false)");
    select.selectedIndex = Channel.Unit - List_Units(Channel.Unit).UnitNum[0];
    cell.title = 'Select channel unit';
    cell.appendChild(select);

    // 7th column (ScaleFactor) - Insert an empty row first and populate it ----------------------------------------
    i = 6;
    cell   = row.insertCell(i);
    Scale_ID        = "Scale_ID_" + Channel.Unique_ID;
    input           = document.createElement('input');
    input.title     = "User specified scale factor";
    input.type      = "number";
    input.style     = "width: 100%; text-align: right"
    input.value     = Channel.ScaleFactor;
    input.id        = Scale_ID;
    input.setAttribute('onchange', "ChannelScaleFactor_Update(this)" );
    cell.appendChild(input);

    // 8th column (Date& Time) - Insert an empty row first and populate it ----------------------------------------
    i = 7;
    cell = row.insertCell(i);
    cell.innerHTML = Channel.DateTime.replace("T", " "); 
    cell.title     = Channel.DateTime.replace("T", " "); 

    // 9th colum (Orientation/Azimuth) - Insert an empty row first and populate it ---------------------------------
    i = 8;
    cell   = row.insertCell(i);
    cell.innerHTML = Channel.Orientation;
    cell.title     = Channel.Orientation;


    // Add Info to Right-Click Menu -------------------------------------------------------------------------------
    if (!Get_UL_List('Right_Click_ul_Sampling',1).includes(Channel.FSamp.toFixed(2))) { 
        // Sampling Rate is not in the list - Add to the list
        Add_li('Right_Click_ul_Sampling',    'FSamp_' + Channel.FSamp.toFixed(2),    Channel.FSamp.toFixed(2) + ' Hz', false);
    }

    if (!Get_UL_List('Right_Click_ul_Type',2).includes(Channel.TypeString)) { 
        // TyprString is not in the list - Add to the list
        Add_li('Right_Click_ul_Type',    'Type_' + Channel.TypeString,    Channel.TypeString, false);
    }

    if (!Get_UL_List('Right_Click_ul_Azimuth',3).includes(Channel.Orientation)) { 
        // Orientation (Azimuth) is not in the list - Add to the list
        Add_li('Right_Click_ul_Azimuth',    'Azimuth_' + Channel.Orientation,    Channel.Orientation, false);
    }


    //--------------------------------------------------------------------------------------------------------------
    // Append to ChannelList
    ChannelList.push(Channel);


    // Now Add a Ployly Graph --------------------------------------------------------------------------------------
    if (Current_Plotly_Num < MaxPlotly_Graphs) {
        Plotly_Create_Graph('panel2', Channel);

        // Check for plotting Graph -  TRUE
        document.getElementById('PlotChecbox_' + Channel.Unique_ID).checked = true;
        Channel.PlotGraph = true;

        // Increase the number of Plotly-elememnts are in panel2
        Current_Plotly_Num++;
    }

    // Checkbox for analysis - TRUE --------------------------------------------------------------------------------
    document.getElementById('FileTreeView_Checkbox_' + Channel.Unique_ID).checked = true;

}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
async function Read_PEER_Data(FileName, delta, dataview) {

    // Reads a PEER NGA Strong Motion Database formatted file (.AT2  .VT2  .DT2)
    // File layout:
    //   Line 1 : "PEER NGA STRONG MOTION DATABASE RECORD"
    //   Line 2 : Event Name, Date, Station, Component  (comma separated)
    //   Line 3 : "ACCELERATION/VELOCITY/DISPLACEMENT TIME SERIES IN UNITS OF <unit>"
    //   Line 4 : "NPTS=  nnnn, DT=  .dddd SEC, ..."
    //   Remaining lines : digitized data values (free-format, multiple values per line)

    // Decleration of variables
    let Line1, Line2, Line3, Line4, temp, temp1;
    let EventInfo, Orientation, SeriesType, UnitOfMeasure, TypeAndUnitIndex, DateTime;
    let NPTS, DT, data=[], Time=[];
    let EndOfFile=false, Ind=0, TNumByte;

    TNumByte = dataview.byteLength;

    // --- Header lines ---------------------------------------------------------------------------
    Line1 = FGetL(dataview);    // "PEER NGA STRONG MOTION DATABASE RECORD"
    Line2 = FGetL(dataview);    // Event Name, Date, Station, Component
    Line3 = FGetL(dataview);    // "<TYPE> TIME SERIES IN UNITS OF <UNIT>"
    Line4 = FGetL(dataview);    // "NPTS=  nnnn, DT=  .dddd SEC, ..."

    // Event Name / Station / Orientation (last comma-separated field is the component/orientation)
    EventInfo   = Line2.split(',').map( s => s.trim() );
    Orientation = EventInfo.at(-1);

    // Series type (ACCELERATION / VELOCITY / DISPLACEMENT) and its unit
    SeriesType    = Line3.trim().split(/\s+/)[0].toUpperCase();
    UnitOfMeasure = Line3.slice( Line3.toUpperCase().indexOf('UNITS OF') + 8 ).trim().toUpperCase();

    // Map (Type, Unit) found in the header to the TypeAndUnit() index used by this webtool
    if      ( SeriesType == 'ACCELERATION' ) {
        if      ( UnitOfMeasure == 'G'      ) { TypeAndUnitIndex = 1;  }
        else if ( UnitOfMeasure == 'M/S2'   ) { TypeAndUnitIndex = 2;  }
        else if ( UnitOfMeasure == 'CM/S2'  ) { TypeAndUnitIndex = 3;  }
        else if ( UnitOfMeasure == 'MM/S2'  ) { TypeAndUnitIndex = 4;  }
        else                                  { TypeAndUnitIndex = 3;  } // default to cm/s² (PEER convention)
    }
    else if ( SeriesType == 'VELOCITY' ) {
        if      ( UnitOfMeasure == 'M/S'    ) { TypeAndUnitIndex = 12; }
        else if ( UnitOfMeasure == 'CM/S'   ) { TypeAndUnitIndex = 13; }
        else if ( UnitOfMeasure == 'MM/S'   ) { TypeAndUnitIndex = 14; }
        else                                  { TypeAndUnitIndex = 13; } // default to cm/s (PEER convention)
    }
    else if ( SeriesType == 'DISPLACEMENT' ) {
        if      ( UnitOfMeasure == 'M'      ) { TypeAndUnitIndex = 22; }
        else if ( UnitOfMeasure == 'CM'     ) { TypeAndUnitIndex = 23; }
        else if ( UnitOfMeasure == 'MM'     ) { TypeAndUnitIndex = 24; }
        else                                  { TypeAndUnitIndex = 23; } // default to cm (PEER convention)
    }
    else { TypeAndUnitIndex = 1; } // fallback to Acceleration (g)

    // NPTS and DT from Line4  (e.g. "NPTS=   1999, DT=   .0100 SEC, ...")
    NPTS = Number( ExtractFloats( Line4.slice( Line4.toUpperCase().indexOf('NPTS') ) )[0] );
    DT   = ExtractFloats( Line4.slice( Line4.toUpperCase().indexOf('DT=') ) )[0];

    // --- Data values ------------------------------------------------------------------------------
    // Values are free-format and may span several values per line, so read line-by-line until NPTS values are collected
    while ( ( data.length < NPTS ) && !EndOfFile ) {
        temp = ExtractFloats( FGetL(dataview) );
        data = data.concat(temp);
    }
    data     = data.slice(0, NPTS);
    Time     = data.map( (v,ii) => ii * DT );
    DateTime = convertDateString(EventInfo[1]);

    // --- Build the Channel object ------------------------------------------------------------------
    let res = new Channel();

    res.FileName         = FileName;
    res.TableName        = FileName.replace(/[-.]/g, '_');
    res.FileListName     = 'FileList_' + res.TableName;
    res.ChNum            = 0;
    res.NumSamples       = data.length;
    res.ScaleFactor      = 1;
    res.Orientation      = Orientation;
    res.DateTime         = DateTime;                          // PEER files do not include absolute date/time
    res.Duration         = Time.at(-1);
    res.DateTime_End     = ComputeEndDateTime(res.DateTime, res.Duration);
    res.Lat              = undefined;
    res.Long             = undefined;
    res.FSamp            = 1 / DT;
    res.delt             = DT;

    temp                  = TypeAndUnit( TypeAndUnitIndex );
    res.TypeAndUnits      = temp.TypeAndUnit;
    res.Type              = temp.Type;
    res.TypeString        = temp.Type_String;
    res.Unit              = temp.Unit;
    res.UnitString        = temp.Unit_String;

    temp1                    = IntervalTypeAndUnit(1);          // Time, starts from zero, constant delt
    res.IntervalTypeAndUnit  = temp1.IntervalTypeAndUnit;
    res.IntervalType         = temp1.Type;
    res.IntervalTypeString   = temp1.Type_String;
    res.IntervalUnit         = temp1.Unit;
    res.IntervalUnitString   = temp1.Unit_String;

    res.InstFreq          = undefined;
    res.InstPeriod        = undefined;
    res.InstDamp          = undefined;

    res.data              = data;
    res.time              = Time;

    // Statistics
    temp         = Statistics(res.data, res.ScaleFactor);
    res.Peak     = temp.Peak;
    res.Mean     = temp.Mean;
    res.RMS      = temp.RMS;
    res.Residual = res.data.at(-1) * res.ScaleFactor;

    // Add to the Main Table and Tree View
    await Add_To_Table( res );

    // Update the ProgressBar
    await UpdateProgress(delta, "ProgressBar_Label");
    await sleep(5);

    // --- Local helper functions --------------------------------------------------------------------

    function ExtractFloats(str) {
        // Extracts all floating point numbers (including exponential notation) from a string
        let regex, match, res=[];
        regex = /[+-]?\d*\.?\d+(?:[eE][+-]?\d+)?/g;
        while (match = regex.exec(str)) {
            res.push( Number(match[0]) );
        }
        return res;
    }

    function Read_Int8(dataview) {
        // Reads the next 1-byte and returns it

        // Return an empty string if EndOfFile
        if (EndOfFile) { return ''; }

        // Read the one-byte
        let content = dataview.getInt8( Ind );

        // Increase the Ind-byte index by 1
        Ind++;

        // Check EndOfFile
        if (Ind >= TNumByte) { EndOfFile = true; }

        // Return the content of the 1-byte read
        return content;
    }

    function FGetL(dataview) {
        // Reads the file byte-by-byte until the next Line Feed (\n)
        // Converts the bytes read to a string and returns it
        let dat=[], content;

        while (true) {
            // Break the while-loop if EndOfFile
            if (EndOfFile) { break; }

            // Read Signed Integer of one-byte
            content = Read_Int8( dataview );

            // If the content is Carriage Return (\r), then continue reading, but do not include it in the dat[] array
            if (content == 13) { continue; }

            // If the content is a Line Feed (\n), then end-of-line is reached, so break the while-loop and return the dat[] array
            if (content == 10) { break; }

            // Append the content to the dat[] array
            dat.push(content);
        };

        // Convert dat[] array into string and return
        return String.fromCharCode.apply(null, dat);
    }

    function convertDateString(a) {
        const [month, day, year] = a.split("/").map(Number);
        const mm = String(month).padStart(2, "0");
        const dd = String(day).padStart(2, "0");
        return `${year}-${mm}-${dd}T00:00:00.00`;
    }
}
//-------------------------------------------------------------------------------------------------------------
async function Read_DAT_GEOSIG(FileName, delta, dataview) {

    // Decleration of variables 
    let DateTime, DateTime_previous, DateTime_current, Data=[], NumSamp_previous, NumSamp_current, Flags, Karar, Opt={}
    let version, StationName, ChannelName, ChannelUnit, ChannelLSB, Flags1, FSamp, Reserved, Duration, TypeAndUnits, temp1, temp2
    let i=0

    // Read the File Header  (66 bytes total)
    version     = dataview.getUint16(i, 'littleEndian');   i += 2;  // Version            (2 bytes)
    StationName = GetCharacters(dataview, i, 6);           i += 6;  // Station name       (6 bytes)
    ChannelName = GetCharacters(dataview, i, 6);           i += 6;  // Channel Name       (6 bytes)
    ChannelUnit = GetCharacters(dataview, i, 6);           i += 6;  // channel units      (6 bytes)
    ChannelLSB  = dataview.getFloat64(i, 'littleEndian');  i += 8;  // Channel LSB        (8 bytes)
    FSamp       = dataview.getUint16(i, 'littleEndian');   i += 2;  // Sampling rate      (2 bytes)
    Flags1      = GetCharacters(dataview, i, 4);           i += 4;  // Flags              (4 bytes)
    Reserved    = GetCharacters(dataview, i, 32);          i += 32; // Reserved          (32 bytes)

    // Read the first block
    [DateTime_previous, NumSamp_previous, Flags, i] = Read_Dat_Block(dataview, i); 

    // DateTime of the first sample in the file
    DateTime =DateTime_previous;

    // Read the rest of the blocks in the file 
    if (DateTime_previous != -1) { Karar = true; } else { Karar = false; }
    
    while (Karar) {
        // Read the next data-block
        Opt = {Date: DateTime_previous, NumSamp: NumSamp_previous};
        [DateTime_current, NumSamp_current, Flags, i] = Read_Dat_Block(dataview, i, Opt); 

        // Check the SynchroChars in the first 2-byte of the block
        if (DateTime_current == -1) { Karar = false; break; }
        
        // Sway DataTime between current and previous blocks
        DateTime_previous = DateTime_current;
    }

    // Number of Samples in Data-array
    NumSamp_current = Data.length;

    Duration = NumSamp_current / FSamp;

    // Type and Unit of channel
    TypeAndUnits = 1; // Unknown, so it is assumed that this channel contains Acc readings of 'g' unit.
    temp1        = TypeAndUnit( TypeAndUnits );
    temp2        = IntervalTypeAndUnit( 1 );   // Type And Unit number - refer to the list (Time series)

    // Create new Channel Object and populate it
    let res = new Channel();

    res.FileName       = FileName;                              // FileName
    res.TableName      = FileName.replace(/[-.]/g, '_');
    res.FileListName   = 'FileList_' + res.TableName;
    res.ChNum          = ChannelName.replace(/\u0000/g, '');    // Channel number
    res.NumSamples     = NumSamp_current;          // Number of digitized data points
    res.ScaleFactor    = 1;                        // User defined Scale Factor - Taken 1.0 if not specified
    res.Orientation    = 'N/A';                    // Orientation of channel
    res.DateTime       = DateTime;                 // Date&Time of the first sample
    res.Duration       = Duration;                 // Total duration in seconds
    res.DateTime_End   = ComputeEndDateTime(res.DateTime, res.Duration); // Date&Time of the last sample in the record
    res.Lat            = undefined;                // Latitude coordinate of the sensor
    res.Long           = undefined;                // Longitude coordinate of the sensor
    res.FSamp          = FSamp;                    // Sampling Frequency in Hz
    res.delt           = 1/FSamp;                  // Time interval of the digitized record in seconds
    res.TypeAndUnits   = TypeAndUnits;             // Type And Unit number - refer to the list
    res.Type           = temp1.Type;               // Type of sensor reading - number
    res.TypeString     = temp1.Type_String;        // Type of sensor reading - string
    res.Unit           = temp1.Unit;               // Unit of sensor reading - number
    res.UnitString     = temp1.Unit_String;        // Unit of sensor reading - string

    res.IntervalTypeAndUnit = temp2.IntervalTypeAndUnit;  // Type And Unit number - refer to the list
    res.IntervalType        = temp2.Type;                 // (0-1)
    res.IntervalTypeString  = temp2.Type_String;          // (Time, Spectral Period, etc.)
    res.IntervalUnit        = temp2.Unit;                 // (0-2)
    res.IntervalUnitString  = temp2.Unit_String;          // (Second, DateTime, etc.)

    res.data           = Data;                              // Digitized data
    res.time           = Data.map((vv,ii) => ii / FSamp);   // Time array of the digitized data

    // Calculate Statictics
    temp1        = Statistics(res.data, res.ScaleFactor);
    res.Peak     = temp1.Peak;
    res.Mean     = temp1.Mean;
    res.RMS      = temp1.RMS;
    res.Residual = res.data.at(-1) * res.ScaleFactor;

    // Add to the Main Table and Tree View
    await Add_To_Table( res );

    // Update the ProgressBar per channel (Unlike other type of channels)
    await UpdateProgress(delta, "ProgressBar_Label");
    await sleep(5);



    function getSecondsDifference(datetime1, datetime2) {
        // Parse the datetime strings into Date objects
        const date1 = new Date(datetime1);
        const date2 = new Date(datetime2);
        
        // Calculate the difference in milliseconds
        const diffInMs = Math.abs(date2 - date1);
        
        // Convert to seconds
        const diffInSeconds = diffInMs / 1000;
        
        return diffInSeconds;
    }
    

    function Read_Dat_Block(dataview, i, Opt) {

        let SyncChar, NumSamp, Year_s, Month_s, DayOfWeek_s, Day_s, Hour_s, Minute_s
        let Second_s, MilliSecond_s, TimeQuality, byte, Flags=[], Data_Block=[], DateTime_Block
        let Time_Diff_sec, TolalNumSamples_InGap, Temp_Gap_Array=[]
        let j

        // Read the header of the block

        // The start of the block must be "=>" (0x3E3D) ----- If not, do not continue.
        if (i + 2 <= dataview.byteLength) { 
            SyncChar = dataview.getUint16(i, true); i += 2; // SynchroChars " => "  (2 bytes)
        } else {  
            return [-1, -1, -1, -1, i];       
        }
        NumSamp         = dataview.getUint32(i, 'littleEndian');   i += 4;  // # of data Samples in the block (4 bytes)
        Year_s          = dataview.getUint16(i, 'littleEndian');   i += 2;  // Year          (2 bytes)
        Month_s         = dataview.getUint16(i, 'littleEndian');   i += 2;  // Month         (2 bytes)
        DayOfWeek_s     = dataview.getUint16(i, 'littleEndian');   i += 2;  // DayOfWeek     (2 bytes)
        Day_s           = dataview.getUint16(i, 'littleEndian');   i += 2;  // Day           (2 bytes)
        Hour_s          = dataview.getUint16(i, 'littleEndian');   i += 2;  // Hour          (2 bytes)
        Minute_s        = dataview.getUint16(i, 'littleEndian');   i += 2;  // Minute        (2 bytes)
        Second_s        = dataview.getUint16(i, 'littleEndian');   i += 2;  // Second        (2 bytes)
        MilliSecond_s   = dataview.getUint16(i, 'littleEndian');   i += 2;  // MilliSecond   (2 bytes)
        TimeQuality     = dataview.getInt8(i, 'littleEndian');     i += 1;  // TimeQuality   (1 byte)
        
        byte      = dataview.getUint8(i);                    
        Flags[0] = (byte >> 0) & 1;  // bit 0   --> If this bit is true (1), then some data samples between this and the previous data block are lost
        Flags[1] = (byte >> 1) & 1;  // bit 1   --> If this bit is true (1), then the binary data in the block is 2-byte signed; otherwise, it is 4-byte signed.
        Flags[2] = (byte >> 2) & 1;  // bit 2   --> Not used.
        Flags[3] = (byte >> 2) & 1;  // bit 3   --> Not used.
        Flags[4] = (byte >> 2) & 1;  // bit 4   --> Not used.
        Flags[5] = (byte >> 2) & 1;  // bit 5   --> Not used.
        Flags[6] = (byte >> 2) & 1;  // bit 6   --> Not used.
        Flags[7] = (byte >> 2) & 1;  // bit 7   --> Not used.
        i += 1;  // Flags         (1 byte)

        // DateTime of the first sample in the block
        DateTime_Block =  DateTime_Construct(Year_s, Month_s, Day_s, Hour_s, Minute_s, Second_s, MilliSecond_s);

        // Check for data gab
        if ((Flags[0] == true) && (Opt != null)) {
            // There is a gab between this block and the previous block

            // Calculate the total time difference
            Time_Diff_sec = Math.abs(getSecondsDifference(Opt.Date, DateTime_Block));

            // Calculate the total number of samples between two blocks
            TolalNumSamples_InGap = Math.round(Time_Diff_sec) * FSamp - Opt.NumSamp;

            // Pre-allocate an array
            Data_Block = new Array(NumSamp + TolalNumSamples_InGap).fill(NaN);
        }
        else {
            // Pre-allocate an array of NumSamp-samples.
            Data_Block = new Array(NumSamp);

            TolalNumSamples_InGap = 0;
        }

        // Read the data in this block
        if (Flags[1]) {
            // The second-bit in the Flags-array is set to true; therefore, the data is 2-byte signed
            for (j = TolalNumSamples_InGap; j < NumSamp+TolalNumSamples_InGap; j++) {
                Data_Block[j] = dataview.getInt16(i, 'littleEndian');   i += 2; // 2-byte signed  (2 bytes)
            }
        }
        else {
            // The second-bit in the Flags-array is set to false; therefore, the data is 4-byte signed
            for (j = TolalNumSamples_InGap; j < NumSamp+TolalNumSamples_InGap; j++) {
                Data_Block[j] = dataview.getInt32(i, 'littleEndian');   i += 4; // 2-byte signed  (4 bytes)
            }
        }

        // push the Data_Block to Data-array (sample-by-sample to avoid stack overflow errors)
        for (j = 0; j < Data_Block.length; j++) { Data.push(Data_Block[j]); }

        // Return results
        return [DateTime_Block, NumSamp+TolalNumSamples_InGap, Flags, i];
    }

    function GetCharacters(dataview, offset, n) {
        //  Reads from dataview n-bytes and interprets them as characters.
        let j, stationName = '';
        for (j = 0; j < 6; j++) {
            stationName += String.fromCharCode(dataview.getUint8(offset + j));
        }
        return stationName
    }

    function toStrAndFill(val, n, Opt) {
        let Len, F, i;
        if (Opt==null) {Opt = true;}
        val = String(val);
        Len = val.length;
        if (Len < n) {
            F = n - Len;
            if (Opt) {for (i=0; i<F; i++) { val = '0' + val;  } }
            else     {for (i=0; i<F; i++) { val = val + '0';  } }
        }
        return val;
    }

    function DateTime_Construct(Year, Month, Day, Hour, Minute, Second, MilliSecond) {
        // DateTime of the first sample in the block
        return toStrAndFill(Year, 4) + "-" + toStrAndFill(Month, 2) + "-" + toStrAndFill(Day, 2) + "T" + toStrAndFill(Hour,2 ) + ":" + toStrAndFill(Minute, 2) + ":" + toStrAndFill(Second, 2) + '.' + toStrAndFill(MilliSecond, 3);
    }
}
//-------------------------------------------------------------------------------------------------------------
async function Read_DAT_Free(FileName, delta, dataview) {
    
    const decoder = new TextDecoder("utf-8");
    const uint8   = new Uint8Array(dataview.buffer);

    let start = 0;
    let Data = [];

    // 1. Scan the buffer for lines
    for (let i = 0; i < uint8.length; i++) {

        // Look for newline character (LF = 10)
        if (uint8[i] === 10 || i === uint8.length - 1) {
            // Extract just this line's bytes and convert to string
            let lineBytes = uint8.subarray(start, i);
            let line = decoder.decode(lineBytes).trim();
            
            if (line.length > 0) {
                // Split line into numbers
                let values = line.split(/\s+/).map(Number);
                Data.push(values);
            }
            start = i + 1;
        }
    }

    if (Data.length === 0) return;

    // 2. Map the row-based data into your Channel structures
    const numChannels  = Data[0].length;
    const numSamples   = Data.length;
    const FSamp        = 100;                         // Unknaow, so FSamp = 100Hz is assumed. 
    const TypeAndUnits = 1;                           // Unknown, so it is assumed that this channel contains Acc readings of 'g' unit.
    const temp1        = TypeAndUnit( TypeAndUnits );
    const temp2        = IntervalTypeAndUnit( 1 );    // Type And Unit number - refer to the list (Time series)
    const DateTime     = new Date().toISOString().replace("Z", "");

    for (let col = 0; col < numChannels; col++) {
        let res            = new Channel();
        res.FileName       = FileName;
        res.TableName      = FileName.replace(/[-.]/g, '_');
        res.FileListName   = 'FileList_' + res.TableName;
        res.ChNum          = col;
        res.NumSamples     = numSamples;                        // Number of digitized data points
        res.ScaleFactor    = 1;                                 // User defined Scale Factor - Taken 1.0 if not specified

        res.Orientation    = 'N/A';                             // Orientation of channel
        res.DateTime       = DateTime;                          // Date&Time of the first sample
        res.Duration       = numSamples / FSamp;                // Total duration in seconds
        res.DateTime_End   = ComputeEndDateTime(res.DateTime, res.Duration); // Date&Time of the last sample in the record
        res.Lat            = undefined;                         // Latitude coordinate of the sensor
        res.Long           = undefined;                         // Longitude coordinate of the sensor
        res.FSamp          = FSamp;                             // Sampling Frequency in Hz

        res.delt           = 1/FSamp;                  // Time interval of the digitized record in seconds
        res.TypeAndUnits   = TypeAndUnits;             // Type And Unit number - refer to the list
        res.Type           = temp1.Type;               // Type of sensor reading - number
        res.TypeString     = temp1.Type_String;        // Type of sensor reading - string
        res.Unit           = temp1.Unit;               // Unit of sensor reading - number
        res.UnitString     = temp1.Unit_String;        // Unit of sensor reading - string

        res.IntervalTypeAndUnit = temp2.IntervalTypeAndUnit;  // Type And Unit number - refer to the list
        res.IntervalType        = temp2.Type;                 // (0-1)
        res.IntervalTypeString  = temp2.Type_String;          // (Time, Spectral Period, etc.)
        res.IntervalUnit        = temp2.Unit;                 // (0-2)
        res.IntervalUnitString  = temp2.Unit_String;          // (Second, DateTime, etc.)

        res.data = [];
        for (let row = 0; row < numSamples; row++) { 
            res.data.push(Data[row][col]);                     // Digitized data
        } 
        res.time           =  res.data.map((vv,ii) => ii / FSamp);   // Time array of the digitized data
        
        // Calculate Statictics
        let temp3        = Statistics(res.data, res.ScaleFactor);
        res.Peak     = temp3.Peak;
        res.Mean     = temp3.Mean;
        res.RMS      = temp3.RMS;
        res.Residual = res.data.at(-1) * res.ScaleFactor;

        // Add to the Main Table and Tree View
        await Add_To_Table( res );

        // Update the ProgressBar per channel (Unlike other type of channels)
        await UpdateProgress(delta/numChannels, "ProgressBar_Label");
        await sleep(5);

    }
}
//-------------------------------------------------------------------------------------------------------------
async function Read_VIF_BCSIMS(FileName, delta, dataview) {
    // Read the header information
    // These variables only exist within this block {...}. Does not exist out of this function.
    let nVersion, Year, Month, Day, Hour, Minute, Second, MilliSec, NanoSec;
    let delt, NumOfChannel, NumSampPerCh, DateTime, DeviceName, OptComment, SF;
    let ScaleFactor=[], NumSamples=[], ChNum=[], FSamp=[], Orientation=[], ChName=[], ChOptChStr=[];
    let nDOFNode=[], fDOFDirX=[], fDOFDirY=[], fDOFDirZ=[], TypeAndUnits=[], Type=[], Unit=[], UnitString=[];
    let IntervalTypeAndUnits=[], IntervalType=[], IntervalTypeString=[], IntervalUnit=[], IntervalUnitString=[];
    let TypeString=[], OffsetFactor=[], data=[], time=[], Duration=[];
    let Ch, j, i, temp;
    let SSL = 128;
    let LSL = 1024;

    i = 0;
    nVersion      = dataview.getUint32(i,    'littleEndian');    i += 4;
    Year          = dataview.getUint32(i,    'littleEndian');    i += 4;
    Month         = dataview.getUint32(i,    'littleEndian');    i += 4;
    Day           = dataview.getUint32(i,    'littleEndian');    i += 4;
    Hour          = dataview.getUint32(i,    'littleEndian');    i += 4;
    Minute        = dataview.getUint32(i,    'littleEndian');    i += 4;
    Second        = dataview.getUint32(i,    'littleEndian');    i += 4;
    MilliSec      = dataview.getUint32(i,    'littleEndian');    i += 4;
    NanoSec       = dataview.getBigUint64(i, 'littleEndian');    i += 8; // ???
    delt          = dataview.getFloat64(i,   'littleEndian');    i += 8;
    NumOfChannel  = dataview.getUint32(i,    'littleEndian');    i += 4;
    NumSampPerCh  = dataview.getUint32(i,    'littleEndian');    i += 4;

    DateTime = DateAndTimeCal();

    // Pre-allocate this.data[] matrix
    data = new Array(NumOfChannel).fill().map(() => new Array(NumSampPerCh).fill(NaN));

    // Read the next 128-bytes
    DeviceName = ReadText(dataview, i, SSL);  i += 128;
    OptComment = ReadText(dataview, i, LSL);  i += 1024;

    for (Ch = 0; Ch < NumOfChannel; Ch++) {
        ScaleFactor[Ch]  = 1;
        NumSamples[Ch]   = NumSampPerCh;
        ChNum[Ch]        = Ch;
        FSamp[Ch]        = 1 /delt;
        Orientation[Ch]  = 'N/A';
        ChName[Ch]       = ReadText(dataview, i, SSL);  i += 128;
        ChOptChStr[Ch]   = ReadText(dataview, i, LSL);  i += 1024;
        nDOFNode[Ch]     = dataview.getUint32(i,     'littleEndian');  i += 4;
        fDOFDirX[Ch]     = dataview.getFloat32(i,    'littleEndian');  i += 4;
        fDOFDirY[Ch]     = dataview.getFloat32(i,    'littleEndian');  i += 4;
        fDOFDirZ[Ch]     = dataview.getFloat32(i,    'littleEndian');  i += 4;
        TypeAndUnits[Ch] = dataview.getUint32(i,     'littleEndian');  i += 4;
        temp            = TypeAndUnit( TypeAndUnits[Ch] );
        Type[Ch]        = temp.Type;
        TypeString[Ch]  = temp.Type_String;
        Unit[Ch]        = temp.Unit;
        UnitString[Ch]  = temp.Unit_String;

        SF               = dataview.getFloat64(i,'littleEndian');  i += 8;
        OffsetFactor[Ch] = dataview.getFloat64(i,'littleEndian');  i += 8;

        for (j = 0; j < NumSampPerCh; j++) {
            data[Ch][j] = (SF * dataview.getInt32(i, 'littleEndian')) + OffsetFactor[Ch];  i += 4;
        }

        time[Ch]      = data[0].map((vv,ii) => ii * delt);
        Duration[Ch]  = time[Ch].at(-1);

        temp                     = IntervalTypeAndUnit( 1 );   // Type And Unit number - refer to the list
        IntervalTypeAndUnits[Ch] = temp.IntervalTypeAndUnit;
        IntervalType[Ch]         = temp.Type;                  // (0-1)
        IntervalTypeString[Ch]   = temp.Type_String;           // (Time, Spectral Period, etc.)
        IntervalUnit[Ch]         = temp.Unit;                  // (0-2)
        IntervalUnitString[Ch]   = temp.Unit_String;           // (Second, DateTime, etc.)

        // Create new Channel Object and populate it
        let res = new Channel();

        res.FileName       = FileName;                 // FileName
        res.TableName      = FileName.replace(/[-.]/g, '_');
        res.FileListName   = 'FileList_' + res.TableName;
        res.ChNum          = Ch;                       // Channel number
        res.NumSamples     = NumSamples[Ch];           // Number of digitized data points
        res.ScaleFactor    = 1;                        // User defined Scale Factor - Taken 1.0 if not specified
        res.Orientation    = Orientation[Ch];          // Orientation of channel
        res.DateTime       = DateTime;                 // Date&Time of the first sample
        res.Duration       = Duration[Ch];             // Total duration in seconds
        res.DateTime_End   = ComputeEndDateTime(res.DateTime, res.Duration); // Date&Time of the last sample in the record
        res.Lat            = undefined;                // Latitude coordinate of the sensor
        res.Long           = undefined;                // Longitude coordinate of the sensor
        res.FSamp          = FSamp[Ch];                // Sampling Frequency in Hz
        res.delt           = delt;                     // Time interval of the digitized record in seconds
        res.TypeAndUnits   = TypeAndUnits[Ch];         // Type And Unit number - refer to the list
        res.Type           = Type[Ch];                 // Type of sensor reading - number
        res.TypeString     = TypeString[Ch];           // Type of sensor reading - string
        res.Unit           = Unit[Ch];                 // Unit of sensor reading - number
        res.UnitString     = UnitString[Ch];           // Unit of sensor reading - string

        res.IntervalTypeAndUnit = IntervalTypeAndUnits[Ch];    // Type And Unit number - refer to the list
        res.IntervalType        = IntervalType[Ch];           // (0-1)
        res.IntervalTypeString  = IntervalTypeString[Ch];     // (Time, Spectral Period, etc.)
        res.IntervalUnit        = IntervalUnit[Ch];           // (0-2)
        res.IntervalUnitString  = IntervalUnitString[Ch];     // (Second, DateTime, etc.)

        res.data           = data[Ch];              // Digitized data
        res.time           = time[Ch];              // Time array of the digitized data

        // Calculate Statictics
        temp         = Statistics(res.data, res.ScaleFactor);
        res.Peak     = temp.Peak;
        res.Mean     = temp.Mean;
        res.RMS      = temp.RMS;
        res.Residual = res.data.at(-1) * res.ScaleFactor;

        // Add to the Main Table and Tree View
        await Add_To_Table( res );

        // Update the ProgressBar per channel (Unlike other type of channels)
        await UpdateProgress(delta/NumOfChannel, "ProgressBar_Label");
        await sleep(5);
    }
    

    function DateAndTimeCal() {
        // These variables only exist within this block {...}. Does not exist out of this function.
        let YYYY, MM, DD, hh, mm, ss, ms;

        YYYY = toStrAndFill(Year, 4);
        MM   = toStrAndFill(Month, 2);
        DD   = toStrAndFill(Day, 2);
        hh   = toStrAndFill(Hour, 2);
        mm   = toStrAndFill(Minute, 2);
        ss   = toStrAndFill(Second, 2);
        ms   = toStrAndFill(MilliSec, 3);

        // Assign the DateTime variable of string type of ISO Format
        return YYYY + "-" + MM + "-" + DD + "T" + hh + ":" + mm + ":" + ss + "." + ms;
    }

    function toStrAndFill(val, n, Opt) {
        let Len, F, i;
        if (Opt==null) {Opt = true;}
        val = String(val);
        Len = val.length;
        if (Len < n) {
            F = n - Len;
            if (Opt) {for (i=0; i<F; i++) { val = '0' + val;  } }
            else     {for (i=0; i<F; i++) { val = val + '0';  } }
        }
        return val;
    }

    function ReadText(dataview, j, n) {
        // This function reads the next n-byte starting from the j(th) byte
        // Converts the bytes read to string and returns it.
        let dat=[], i, jj = j;

        for (i = 0; i < n; i++) { dat.push(dataview.getInt8(jj, 'littleEndian'));  jj++; };

        // Convert dat[] array into string and return it
        return String.fromCharCode.apply(null, dat);
    }
}
//-------------------------------------------------------------------------------------------------------------
async function Read_V1_COSMOS(FileName, delta, dataview) {
    let data=[],Time=[], ChannelNum=[], InstFreq=[], InstPeriod=[], InstDamp=[], delt=[], FSamp=[];
    let Lat=[], Long=[], StationNo=[], NumSamples=[], Azimuth, Orientation=[], DateTime=[], StationName=[], EqeTitleName=[];
    let ScaleFactor=[], IsVertical=[], TypeAndUnits=[], Type=[], Unit=[], UnitString=[], TypeString=[], Duration=[];
    let IntervalTypeAndUnits=[], IntervalType=[], IntervalTypeString=[], IntervalUnit=[], IntervalUnitString=[];
    let TNumByte, cN=0, EndOfFile=false, Ind=0;

    TNumByte = dataview.byteLength;

    // Read all channels - File usually includes 3 channels
    while (!EndOfFile) {
        ReadOneChannel( dataview, cN );
        cN++;
    }

    // Update the ProgressBar per channel (Unlike other type of channels)
    await UpdateProgress(delta, "ProgressBar_Label");
    await sleep(5);

    async function ReadOneChannel(dataview, ChNum) {
        // This function reads only one channel data from the file
        let header=[], IntNum=[], ReNum=[], Temp, Temp1=[], Temp2=[], Ind, Ind1, Ind2, Indx=0, i, j, n, m, d, t, temp;

        // Read the first 13 lines
        for (i = 0; i < 13; i++) {
            header = header.concat( FGetL(dataview) );
        }

        // Lines 14-20
        for (i = 0; i < 7; i++) {
            IntNum = IntNum.concat( Str2Int( FGetL(dataview) ) );
        }

        // Lines 21-27
        for (i = 0; i < 7; i++) {
            ReNum = ReNum.concat( Str2Float( FGetL(dataview), 10 ) )
        }

        ChannelNum[ChNum]   = IntNum[0];
        InstFreq[ChNum]     = ReNum[8];
        InstPeriod[ChNum]   = ReNum[0];
        InstDamp[ChNum]     = ReNum[1];
        delt[ChNum]         = ReNum[12] / 1000;    // ReNum[12] is in milliseconds
        FSamp[ChNum]        = 1 / delt[ChNum];
        Lat[ChNum]          = ReNum[28];
        Long[ChNum]         = ReNum[29];
        StationNo[ChNum]    = IntNum[13];
        NumSamples[ChNum]   = IntNum[27];

        Azimuth = IntNum[26];

        if      (Number(Azimuth) == 500) { Orientation[ChNum] = 'UP';  }
        else if (Number(Azimuth) == 600) { Orientation[ChNum] = 'DOWN';}
        else                             { Orientation[ChNum] = Number(Azimuth).toFixed(2); }

        DateTime[ChNum]     = ToString(IntNum[23],4) + "-" + ToString(IntNum[21], 2) + "-" + ToString(IntNum[22], 2) + "T" + ToString(IntNum[16], 2) + ":" + ToString(IntNum[17], 2) + ":" + ToString(IntNum[18], 2) + ".000";
        StationName[ChNum]  = header[5].slice(0, IntNum[29]);
        EqeTitleName[ChNum] = header[7].slice(0, IntNum[30]);
        ScaleFactor[ChNum]  = 1;

        // 500: Vertical Up       600: Vertical Down
        if ((IntNum[26] == 500) || (IntNum[26] = 600)) { IsVertical[ChNum] = true; }

        // Type and Unit of channel
        TypeAndUnits[ChNum] = 1; // Unknown, so it is assumed that this channel contains Acc readings of 'g' unit.
        temp                = TypeAndUnit( TypeAndUnits[ChNum] );
        Type[ChNum]         = temp.Type;
        TypeString[ChNum]   = temp.Type_String;
        Unit[ChNum]         = temp.Unit;
        UnitString[ChNum]   = temp.Unit_String;

        // Read the first line after the header
        Temp = FGetL( dataview );
        Ind  = Temp.indexOf(": (");

        if (Ind != -1) {
            // This is the information line before the recordings start
            // Format: (8f9.6)   n=8,  m=9,  d=6
            // 8f9.6 means 8 columns, each column includes 9 characters (numbers) with 6 decimal numbers
            n = parseInt( Temp[Ind + 3] );
            m = parseInt( Temp[Ind + 5] );
            d = parseInt( Temp[Ind + 7] );

            // Number of lines to read
            t = Math.ceil( NumSamples[ChNum] / n );

            data[ChNum] = new Array(NumSamples[ChNum]).fill();

            for (i = 0; i < t; i++) {
                Temp = Str2Float( FGetL(dataview), m);
                for (j = 0; j < Temp.length; j++) {
                    data[ChNum][Indx] = Temp[j];
                    Indx++;
                }
            }

            Time[ChNum]     = data[ChNum].map((v,i) => i * delt[ChNum] );
            Duration[ChNum] = Time[ChNum].at(-1);

            // End of record line
            FGetL(dataview);
        }
        else {
            // No information line in the file, so start reading the time-data pairs
            // 10 columns of Time-Data pairs
            data = [];
            Time = [];
            Ind1 = [0,2,4,6,8];
            Ind2 = [1,3,5,7,9];

            // Number of lines to read
            t = Math.floor( NumSamples[ChNum] / 5 );

            data[ChNum] = new Array(NumSamples[ChNum]).fill();
            Time[ChNum] = new Array(NumSamples[ChNum]).fill();

            // Temp is already populated above
            Temp = Str2Float( Temp, 7 );

            for (i = 0; i < t; i++) {
                Temp1 = GetByIndices(Temp, Ind2);
                Temp2 = GetByIndices(Temp, Ind1);
                for (j = 0; j < Temp.length/2; j++) {
                    data[ChNum][Indx] = Temp1[j];
                    Time[ChNum][Indx] = Temp2[j];
                    Indx++;
                }
                Temp  = Str2Float( FGetL(dataview), 7 );
            }

            // Read the last line of data - This line may contain less number of samples
            for (j = 0; j < Temp.length; j+=2) {
                data[ChNum][Indx] = Temp[j+1];
                Time[ChNum][Indx] = Temp[j];
                Indx++
            }

            Duration[ChNum] = Time[ChNum].at(-1);

            // End of record line
            FGetL(dataview);
        }

        temp                        = IntervalTypeAndUnit( 1 );   // Type And Unit number - refer to the list
        IntervalTypeAndUnits[ChNum] = temp.IntervalTypeAndUnit;
        IntervalType[ChNum]         = temp.Type;                  // (0-1)
        IntervalTypeString[ChNum]   = temp.Type_String;           // (Time, Spectral Period, etc.)
        IntervalUnit[ChNum]         = temp.Unit;                  // (0-2)
        IntervalUnitString[ChNum]   = temp.Unit_String;           // (Second, DateTime, etc.)


        // Create new Channel Object and populate it
        let res = new Channel();

        res.FileName       = FileName;             // FileName
        res.TableName      = FileName.replace(/[-.]/g, '_');
        res.FileListName   = 'FileList_' + res.TableName;
        res.ChNum          = ChNum;                // Channel number
        res.NumSamples     = NumSamples[ChNum];    // Number of digitized data points
        res.ScaleFactor    = 1;                    // User defined Scale Factor - Taken 1.0 if not specified
        res.Orientation    = Orientation[ChNum];   // Orientation of channel
        res.DateTime       = DateTime[ChNum];      // Date&Time of the first sample
        res.Duration       = Duration[ChNum];      // Total duration in seconds
        res.DateTime_End   = ComputeEndDateTime(res.DateTime, res.Duration); // Date&Time of the last sample in the record
        res.Lat            = Lat[ChNum];           // Latitude coordinate of the sensor
        res.Long           = Long[ChNum];          // Longitude coordinate of the sensor
        res.FSamp          = FSamp[ChNum];         // Sampling Frequency in Hz
        res.delt           = delt[ChNum];          // Time interval of the digitized record in seconds
        res.TypeAndUnits   = TypeAndUnits[ChNum];  // Type And Unit number - refer to the list
        res.Type           = Type[ChNum];          // Type of sensor reading - number
        res.TypeString     = TypeString[ChNum];    // Type of sensor reading - string
        res.Unit           = Unit[ChNum];          // Unit of sensor reading - number
        res.UnitString     = UnitString[ChNum];    // Unit of sensor reading - string
        res.IntervalTypeAndUnit = IntervalTypeAndUnits[ChNum];    // Type And Unit number - refer to the list
        res.IntervalType        = IntervalType[ChNum];           // (0-1)
        res.IntervalTypeString  = IntervalTypeString[ChNum];     // (Time, Spectral Period, etc.)
        res.IntervalUnit        = IntervalUnit[ChNum];           // (0-2)
        res.IntervalUnitString  = IntervalUnitString[ChNum];     // (Second, DateTime, etc.)
        res.InstFreq       = InstFreq[ChNum];      // Natural frequency of transducer (in Hz).
        res.InstPeriod     = InstPeriod[ChNum];    // Natural period of transducer (in seconds).
        res.InstDamp       = InstDamp[ChNum];      // Damping ratio of transducer (fraction of critical).
        res.data           = data[ChNum];          // Digitized data
        res.time           = Time[ChNum];          // Time array of the digitized data

        // Calculate Statictics
        temp         = Statistics(res.data, res.ScaleFactor);
        res.Peak     = temp.Peak;
        res.Mean     = temp.Mean;
        res.RMS      = temp.RMS;
        res.Residual = res.data.at(-1) * res.ScaleFactor;

        // Add to the Main Table and Tree View
        await Add_To_Table( res );
    }


    // Helper functions
    function Str2Float(Str, n) {
        // This function converts each consecutive n-character in Str[] array into an Float number
        let i, j, S=[], L=[];

        for (i = 0; i < Str.length; i+=n) {
            for (j=0; j<n; j++) { S+=Str[i+j];}
            L.push( parseFloat(S) );
            S=[];
        }
        return L;
    }

    function Str2Int(Str) {
        // This function converts each consecutive 5-character in Str[] array into an Integer number
        let i, j, S=[], L = [];

        for (i = 0; i < Str.length; i+=5) {
            for (j=0; j<5; j++) { S+=Str[i+j];}
            L.push(parseInt(S));
            S=[];
        }
        return L;
    }

    function Read_Int8(dataview) {
        // This function reads the next 1-byte and returns it

        // Return an empty string it is EndOfFile
        if (EndOfFile) { return ''; }

        // Read the one-byte
        let content = dataview.getInt8( Ind );

        // Increase the Ind-byte index by 1
        Ind++;

        // Check EndOfFile
        if (Ind >= TNumByte) { EndOfFile = true; }

        // Return the content of the 1-byte read
        return content;
    }

    function FGetL(dataview) {
        // This function reads the file byte-by-byte until the next Line Feed (\n)
        // Converts the bytes read to a string and returns it
        let dat=[], content;

        while (true) {
            // Break the while-loop if EndOfFile
            if (EndOfFile) { break; }

            // Read Signed Integer of one-byte
            content = Read_Int8( dataview );

            // If the content is Carriage Return (\r), then continue reading, but do not include it in the dat[] array
            if (content == 13) { continue; }

            // If the content is a Line Feed (\n), then EndOfFile is reached, so break the while-loop and return the dat[] array
            if (content == 10) { break; }

            // Append the content to the dat[] array
            dat.push(content);
        };

        // Convert dat[] array into string and return
        return String.fromCharCode.apply(null, dat);

    }

    function ToString(val, n, Opt) {
        let Len, F, i;

        if (Opt==null) {Opt = true;}

        val = String(val);
        Len = val.length;
        if (Len < n) {
             F = n - Len;
            if (Opt) {for (i=0; i<F; i++) { val = '0' + val;  } }
            else     {for (i=0; i<F; i++) { val = val + '0';  } }
        }
        return val;
    }
}
//-------------------------------------------------------------------------------------------------------------
async function Read_V1c_COSMOS( FileName, delta, dataview ) { 
    let TNumByte, EndOfFile=false;
    let i, Ind=0, header=[], IntHeaderNum=[], RealHeaderNum=[], data=[], Time=[];
    let H1, NumOfLines, NumCommentLines, NumOfDataValues;
    let count_Integer, widthChar_Integer, repeatCount, fieldWidth, decimalPlace;
    let Type, Scale, temp1, temp2, Duration;

    // Total number of bytes in the file 
    TNumByte = dataview.byteLength;

    // First Line in header
    H1         = FGetL(dataview);
    NumOfLines = Number(H1.substring( 46,  48)); 

    // Read the rest of header lines
    for (i = 0; i < NumOfLines-1; i++) {
        header = header.concat( FGetL(dataview) );
    }
    
    // Read Integer-header values
    H1                = FGetL(dataview);                   // First line in Integer-header string
    NumOfLines        = Number(H1.substring( 37,  40));    // Number of lines of Integer Header values that follow.
    [count_Integer, widthChar_Integer] = extractFormat(H1.substring( 56,  80));

    for (i = 0; i < NumOfLines; i++) {
        IntHeaderNum.push(parseFixedIntegers(FGetL(dataview), count_Integer, widthChar_Integer));
    }

    // Read Real-header values
    H1              = FGetL(dataview);                      // First line in Real-header string
    NumOfLines      = Number(H1.substring( 34,  37));       // Number of lines of Real Header values that follow.
   [repeatCount, fieldWidth, decimalPlace] = extractFormatNumbers(H1.substring( 53,  80));

    for (i = 0; i < NumOfLines; i++) {
        RealHeaderNum.push(parseFixedWidthFloats(FGetL(dataview), repeatCount, fieldWidth, decimalPlace));
    }

    // Comments 
    H1                = FGetL(dataview);
    NumCommentLines   = Number(H1.substring(  0,   4));     // Number of comment lines that follow
    for (i = 0; i < NumCommentLines; i++){ FGetL(dataview); }

    // Data Values Section 
    H1                = FGetL(dataview);
    NumOfDataValues   = Number(H1.substring(  0,   8));    // Number of data points following;
    [repeatCount, fieldWidth, decimalPlace] = extractFormatNumbers(H1.substring( 70,  80));
    for (i = 0; i < NumOfDataValues; i++){
        data.push(Number( FGetL(dataview) ));
    }

    // Date&Time
    let Year     = IntHeaderNum[3][9];
    let Month    = IntHeaderNum[4][1];
    let Day      = IntHeaderNum[4][2];
    let Hour     = IntHeaderNum[4][3];
    let Minute   = IntHeaderNum[4][4];
    let Second   = RealHeaderNum[5][4];
    let DateTime = ToString(Year,4) + "-" + ToString(Month, 2) + "-" + ToString(Day, 2) + "T" + ToString(Hour, 2) + ":" + ToString(Minute, 2) + ":" + Second.toFixed(3);

    // Date Type, Unit, Sacle factor
    let Physical_Parameter  = IntHeaderNum[0][1];    // Table-1
    let Units_Of_Data       = IntHeaderNum[0][2];    // Table 2
    let Type_Of_Record      = IntHeaderNum[0][4];    // Table 3

    [Type, Scale] = TypeUnit(Physical_Parameter, Units_Of_Data, Type_Of_Record);

    temp1 = TypeAndUnit(Type);
    temp2 = IntervalTypeAndUnit(1);

    let Lat     = RealHeaderNum[0][0];
    let Long    = RealHeaderNum[0][1];
    let Azimuth = AzimuthDet(RealHeaderNum[9][1]).toString();

    let delt  = RealHeaderNum[6][3];
    let FSamp = 1 / delt;
    Time      = data.map((v,i) => i * delt );
    Duration  = Time.at(-1);

    let SensorNaturalFrequency = RealHeaderNum[7][4];  // Sensor natural frequency (Hz)
    let SensorNaturalPeriod    = 1 / SensorNaturalFrequency; 
    let SensorDamping          = RealHeaderNum[8][0];  // Sensor damping (fraction of critical)
    let SensorSensitivity      = RealHeaderNum[8][1];  // Sensor sensitivity (for accelerometer, in volts/g (cm/g if film record);  for other sensors, in volts pe 0    motion unit given in Ihdr(3)).

    // Create new Channel Object and populate it
    let res = new Channel();
    
    res.FileName            = FileName;                                        // FileName
    res.TableName           = FileName.replace(/[-.]/g, '_');
    res.FileListName        = 'FileList_' + res.TableName;
    res.ChNum               = 0;                                               // Channel number
    res.NumSamples          = NumOfDataValues;                                 // Number of digitized data points
    res.ScaleFactor         = Scale;                                           // User defined Scale Factor - Taken 1.0 if not specified
    res.Orientation         = Azimuth;                                         // Orientation of channel
    res.DateTime            = DateTime;                                        // Date&Time of the first sample
    res.Duration            = Duration;                                        // Total duration in seconds
    res.DateTime_End        = ComputeEndDateTime(res.DateTime, res.Duration);  // Date&Time of the last sample in the record
    res.Lat                 = Lat;                                             // Latitude coordinate of the sensor
    res.Long                = Long;                                            // Longitude coordinate of the sensor
    res.FSamp               = FSamp;                                           // Sampling Frequency in Hz
    res.delt                = delt;                                            // Time interval of the digitized record in seconds
    res.TypeAndUnits        = Type;                                            // Type And Unit number - refer to the list
    res.Type                = temp1.Type;                                      // Type of sensor reading - number
    res.TypeString          = temp1.Type_String;                               // Type of sensor reading - string
    res.Unit                = temp1.Unit;                                      // Unit of sensor reading - number
    res.UnitString          = temp1.Unit_String;                               // Unit of sensor reading - string
    res.IntervalTypeAndUnit = temp2.IntervalTypeAndUnit;                       // Type And Unit number - refer to the list
    res.IntervalType        = temp2.Type;                                      // (0-1)
    res.IntervalTypeString  = temp2.Type_String;                               // (Time, Spectral Period, etc.)
    res.IntervalUnit        = temp2.Unit;                                      // (0-2)
    res.IntervalUnitString  = temp2.Unit_String;                               // (Second, DateTime, etc.)
    res.InstFreq            = SensorNaturalFrequency;                          // Natural frequency of transducer (in Hz).
    res.InstPeriod          = SensorNaturalPeriod;                             // Natural period of transducer (in seconds).
    res.InstDamp            = SensorDamping;                                   // Damping ratio of transducer (fraction of critical).
    res.data                = data;                                            // Digitized data
    res.time                = Time;                                            // Time array of the digitized data

    // Calculate Statictics
    temp1         = Statistics(res.data, res.ScaleFactor);
    res.Peak     = temp1.Peak;
    res.Mean     = temp1.Mean;
    res.RMS      = temp1.RMS;
    res.Residual = res.data.at(-1) * res.ScaleFactor;

    // Add to the Main Table and Tree View
    await Add_To_Table( res );

    // Update the ProgressBar per channel (Unlike other type of channels)
    await UpdateProgress(delta, "ProgressBar_Label");
    await sleep(5);

    // heper functions 
    function FGetL(dataview) {
        // This function reads the file byte-by-byte until the next Line Feed (\n)
        // Converts the bytes read to a string and returns it
        let dat=[], content;

        while (true) {
            // Break the while-loop if EndOfFile
            if (EndOfFile) { break; }

            // Read Signed Integer of one-byte
            content = Read_Int8( dataview );

            // If the content is Carriage Return (\r), then continue reading, but do not include it in the dat[] array
            if (content == 13) { continue; }

            // If the content is a Line Feed (\n), then EndOfFile is reached, so break the while-loop and return the dat[] array
            if (content == 10) { break; }

            // Append the content to the dat[] array
            dat.push(content);
        };

        // Convert dat[] array into string and return
        return String.fromCharCode.apply(null, dat);

    }
    function Read_Int8(dataview) {
        // This function reads the next 1-byte and returns it

        // Return an empty string it is EndOfFile
        if (EndOfFile) { return ''; }

        // Read the one-byte
        let content = dataview.getInt8( Ind );

        // Increase the Ind-byte index by 1
        Ind++;

        // Check EndOfFile
        if (Ind >= TNumByte) { EndOfFile = true; }

        // Return the content of the 1-byte read
        return content;
    }
    function extractFormat(str) {
        // str = (10I8)
        const match = str.match(/\((\d+)I(\d+)\)/);
        
        if (!match) { return [-1, -1]; }

        let count = parseInt(match[1]);   // 10
        let width = parseInt(match[2]);   // 8

        return [count, width]
    }
    function parseFixedIntegers(str, count, width) {
        const results = [];
        for (let i = 0; i < count; i++) {
            const chunk = str.slice(i * width, (i + 1) * width);
            results.push(parseInt(chunk, 10));
        }
        return results;
    }
    function extractFormatNumbers(str) {
        // Remove parentheses and letters, then split on '.' to handle decimal notation
        const match = str.match(/\((\d+)[A-Z](\d+)\.(\d+)\)/i);

        if (!match) { return [-1, -1, -1]; }
        
        let repeatCount  =  parseInt(match[1]); // repeat count → 5
        let fieldWidth   =  parseInt(match[2]); // field width  → 15
        let decimalPlace =  parseInt(match[3]); // decimal places → 6

        return [repeatCount, fieldWidth, decimalPlace];
    }
    function parseFixedWidthFloats(str, repeatCount, fieldWidth, decimalPlace) {
        const values = [];
        for (let i = 0; i < repeatCount; i++) {
            const field = str.substring(i * fieldWidth, (i + 1) * fieldWidth);
            values.push(parseFloat(field));
        }
        return values;
    }
    function ToString(val, n, Opt) {
        let Len, F, i;

        if (Opt==null) {Opt = true;}

        val = String(val);
        Len = val.length;
        if (Len < n) {
             F = n - Len;
            if (Opt) {for (i=0; i<F; i++) { val = '0' + val;  } }
            else     {for (i=0; i<F; i++) { val = val + '0';  } }
        }
        return val;
    }
    function TypeUnit(Physical_Parameter, Units_Of_Data, Type_Of_Record) {


        let Type  = 1;
        let Scale = 1;

        if      ([1, 10].includes(Physical_Parameter)) { 

            // Acceleration (1)
            // Angular acceleration (10)

            if (Units_Of_Data ==  2) { Type =  1; Scale = 1;     } 
            if (Units_Of_Data ==  4) { Type =  3; Scale = 1;     }
            if (Units_Of_Data ==  7) { Type =  3; Scale = 2.54;  }
            if (Units_Of_Data == 10) { Type =  3; Scale = 1;     }
            if (Units_Of_Data == 11) { Type =  1; Scale = 1e-3;  }
            if (Units_Of_Data == 12) { Type =  1; Scale = 1e-6;  }
            if (Units_Of_Data == 23) { Type =  1; Scale = 1;     }
            if (Units_Of_Data == 50) { Type =  1; Scale = 1;     }
            if (Units_Of_Data == 51) { Type =  1; Scale = 1;     }
            if (Units_Of_Data == 52) { Type =  1; Scale = 1;     }

        }
        else if ([2, 11].includes(Physical_Parameter)) { 

            // Velocity (2)
            // Angular velocity (11)
            
            if (Units_Of_Data ==  5) { Type =  13; Scale = 1;     } 
            if (Units_Of_Data ==  8) { Type =  13; Scale = 2.54;  } 
            if (Units_Of_Data == 24) { Type =  13; Scale = 1;     } 

        } 
        else if ([3, 4, 12].includes(Physical_Parameter)) { 

            // Displacement (absolute) 3
            // Displacement (relative) 4
            // Angular dispacement (12)

            if (Units_Of_Data ==  6) { Type =  23; Scale = 1;     } 
            if (Units_Of_Data ==  9) { Type =  23; Scale = 2.54;  } 
            if (Units_Of_Data == 25) { Type =  23; Scale = 1;     } 

        } 
        else if ([20, 21].includes(Physical_Parameter)) { 

            // Presure absolute (20)
            // Presure relative (21)

            if (Units_Of_Data == 60) { Type =  81; Scale = 6894.76; } 

        } 
        else if ([30, 31].includes(Physical_Parameter)) { 

            // Volumetric strain (30)
            // Linear Strain (31)

            if ([60].includes(Units_Of_Data)) { Type =  31; Scale = 1e-6; } 
        } 

        // return 
        return [Type, Scale];

    }
    function AzimuthDet(Az) {

        if ((Az >= 1) && (Az <=360)) { return Az; }

        if (Az == 400) { return "UP";       }
        if (Az == 401) { return "DOWN";     }
        if (Az == 402) { return "VERTICAL"; }

        if (Az == 500) { return "RADIAL";     }
        if (Az == 501) { return "TRANSVERSE"; }

        if (Az == 600) { return "LONG";     }
        if (Az == 601) { return "TANG";     }

        if (Az == 700) { return "H1";     }
        if (Az == 701) { return "H2";     }

        return "N/A"

    }
}
//-------------------------------------------------------------------------------------------------------------
async function Read_ASC(FileName, delta, dataview) {

    // Declaration of variables
    let EventName, Event_ID, Event_Date_YYYYMMDD, Event_Time_HHMMSS, Event_Latitude_Degree, Event_Longitude_Degree;
    let Event_Depth_KM, Hypocenter_Reference, Magnitude_W, Magnitude_W_Reference, Magnitude_L, Magnitude_L_Reference, Focal_Mechanism;
    let Network, Station_Code, Station_Name, Station_Latitude_Degree, Station_Longitude_Degree, Station_Elevation_m, Location;
    let Sensor_Depth_m, VS30_meterPerSecond, Site_Classification_EC8, Morphological_Classification, Epicentral_Distance, Earthquake_BackAzimuth_Degree;
    let DateTime_FirstSample_YYYYMMDD_HHMMSS, DateTime_FirstSample_Precision, Sampling_Interval_s, NumberOfSamples, Duration_s, Stream;
    let Units, Instrument, Instrument_Analog, Instrument_Frequency_Hz, Instrument_Damping, Full_Scale_g, N_Bit_Digital_Converter;
    let PGA_CmPerS2, Time_PGA_S, Baseline_Correction, FilterType, Filter_Order, Low_Cut_Frequency_Hz, High_Cut_Frequency_Hz, Late_Normal_Triggered;
    let Database_Version, Header_Format, Data_Type, Processing, Data_Timestamp_YYYYMMDD_HHMMSS, Data_License, Data_Citation, Data_Creator;
    let Original_Data_Mediator_Citation, Original_Data_Mediator, Original_Data_Creator_Citation, Original_Data_Creator;
    let i, data=[], Time=[], NumOfChannel, delt=[], FSamp=[], Duration=[], ScaleFactor=[], NumSamples=[], ChNum=[], temp=[], temp1=[];
    let Type=[], Unit=[], UnitString=[], TypeString=[], Orientation=[], DateTime=[];
    let IntervalTypeAndUnits=[], IntervalType=[], IntervalTypeString=[], IntervalUnit=[], IntervalUnitString=[];
    let EndOfFile=false, Ind=0, TNumByte, DT, InstFreq, InstPeriod, InstDamp;

    TNumByte = dataview.byteLength;

    EventName                             = RemoveAfter( FGetL(dataview) );
    Event_ID                              = RemoveAfter( FGetL(dataview) );
    Event_Date_YYYYMMDD                   = RemoveAfter( FGetL(dataview) );
    Event_Time_HHMMSS                     = RemoveAfter( FGetL(dataview) );
    Event_Latitude_Degree                 = RemoveAfter( FGetL(dataview) );
    Event_Longitude_Degree                = RemoveAfter( FGetL(dataview) );
    Event_Depth_KM                        = Number( RemoveAfter( FGetL(dataview) ));
    Hypocenter_Reference                  = RemoveAfter( FGetL(dataview) );
    Magnitude_W                           = RemoveAfter( FGetL(dataview) );
    Magnitude_W_Reference                 = RemoveAfter( FGetL(dataview) );
    Magnitude_L                           = RemoveAfter( FGetL(dataview) );
    Magnitude_L_Reference                 = RemoveAfter( FGetL(dataview) );
    Focal_Mechanism                       = RemoveAfter( FGetL(dataview) );
    Network                               = RemoveAfter( FGetL(dataview) );
    Station_Code                          = RemoveAfter( FGetL(dataview) );
    Station_Name                          = RemoveAfter( FGetL(dataview) );
    Station_Latitude_Degree               = RemoveAfter( FGetL(dataview) );
    Station_Longitude_Degree              = RemoveAfter( FGetL(dataview) );
    Station_Elevation_m                   = RemoveAfter( FGetL(dataview) );
    Location                              = RemoveAfter( FGetL(dataview) );
    Sensor_Depth_m                        = RemoveAfter( FGetL(dataview) );
    VS30_meterPerSecond                   = RemoveAfter( FGetL(dataview) );
    Site_Classification_EC8               = RemoveAfter( FGetL(dataview) );
    Morphological_Classification          = RemoveAfter( FGetL(dataview) );
    Epicentral_Distance                   = RemoveAfter( FGetL(dataview) );
    Earthquake_BackAzimuth_Degree         = RemoveAfter( FGetL(dataview) );
    DateTime_FirstSample_YYYYMMDD_HHMMSS  = RemoveAfter( FGetL(dataview) );
    DateTime_FirstSample_Precision        = RemoveAfter( FGetL(dataview) );
    Sampling_Interval_s                   = Number( RemoveAfter( FGetL(dataview) ));
    NumberOfSamples                       = Number( RemoveAfter( FGetL(dataview) ));
    Duration_s                            = RemoveAfter( FGetL(dataview) );
    Stream                                = RemoveAfter( FGetL(dataview) );
    Units                                 = RemoveAfter( FGetL(dataview) );
    Instrument                            = RemoveAfter( FGetL(dataview) );
    Instrument_Analog                     = RemoveAfter( FGetL(dataview) );
    Instrument_Frequency_Hz               = RemoveAfter( FGetL(dataview) );
    Instrument_Damping                    = RemoveAfter( FGetL(dataview) );
    Full_Scale_g                          = RemoveAfter( FGetL(dataview) );
    N_Bit_Digital_Converter               = RemoveAfter( FGetL(dataview) );
    PGA_CmPerS2                           = RemoveAfter( FGetL(dataview) );
    Time_PGA_S                            = RemoveAfter( FGetL(dataview) );
    Baseline_Correction                   = RemoveAfter( FGetL(dataview) );
    FilterType                            = RemoveAfter( FGetL(dataview) );
    Filter_Order                          = RemoveAfter( FGetL(dataview) );
    Low_Cut_Frequency_Hz                  = RemoveAfter( FGetL(dataview) );
    High_Cut_Frequency_Hz                 = RemoveAfter( FGetL(dataview) );
    Late_Normal_Triggered                 = RemoveAfter( FGetL(dataview) );
    Database_Version                      = RemoveAfter( FGetL(dataview) );
    Header_Format                         = RemoveAfter( FGetL(dataview) );
    Data_Type                             = RemoveAfter( FGetL(dataview) );
    Processing                            = RemoveAfter( FGetL(dataview) );
    Data_Timestamp_YYYYMMDD_HHMMSS        = RemoveAfter( FGetL(dataview) );
    Data_License                          = RemoveAfter( FGetL(dataview) );
    Data_Citation                         = RemoveAfter( FGetL(dataview) );
    Data_Creator                          = RemoveAfter( FGetL(dataview) );
    Original_Data_Mediator_Citation       = RemoveAfter( FGetL(dataview) );
    Original_Data_Mediator                = RemoveAfter( FGetL(dataview) );
    Original_Data_Creator_Citation        = RemoveAfter( FGetL(dataview) );
    Original_Data_Creator                 = RemoveAfter( FGetL(dataview) );

    // Skip the next 5-lines
    for (i = 0; i < 5; i++) { FGetL(dataview); }

    // Pre-allocation of data_array
    data[0] = [];

    // start reading data-points
    for (i = 0; i < NumberOfSamples; i++) { data[0].push( ExtractFloats( FGetL(dataview) )[0] ); }

    NumOfChannel = 1;

    for (i = 0; i < NumOfChannel; i++) {
        delt[i]          = Sampling_Interval_s;
        FSamp[i]         = 1 / delt[i];
        Time[i]          = data[i].map((v,ii) => ii * delt[i]);
        Duration[i]      = Time[i].at(-1);
        ScaleFactor[i]   = 1;
        NumSamples[i]    = NumberOfSamples;
        ChNum[i]         = i;
        temp             = TypeAndUnit(3);
        Type[i]          = temp.Type;
        TypeString[i]    = temp.Type_String;
        Unit[i]          = temp.Unit;
        UnitString[i]    = temp.Unit_String;
        Orientation[i]   = Stream;
        DT               = DateTime_FirstSample_YYYYMMDD_HHMMSS;
        DateTime[i]      = DT.slice(0,4) +'-'+  DT.slice(4,6) +'-'+ DT.slice(6,8) +'T'+ DT.slice(9,11) +':'+ DT.slice(11,13) +':'+ DT.slice(13,15) +'.000';

        InstFreq = Number( Instrument_Frequency_Hz );  if (InstFreq == 0) { InstFreq = undefined; InstPeriod = undefined } else { InstPeriod = 1 / InstFreq; }
        InstDamp = Number(Instrument_Damping);         if (InstDamp == 0) { InstDamp = undefined; }

        temp1                   = IntervalTypeAndUnit( 1 );   // Type And Unit number - refer to the list
        IntervalTypeAndUnits[i] = temp1.IntervalTypeAndUnit;
        IntervalType[i]         = temp1.Type;                  // (0-1)
        IntervalTypeString[i]   = temp1.Type_String;           // (Time, Spectral Period, etc.)
        IntervalUnit[i]         = temp1.Unit;                  // (0-2)
        IntervalUnitString[i]   = temp1.Unit_String;           // (Second, DateTime, etc.)

        // Create new Channel Object and populate it
        let res = new Channel();

        res.FileName       = FileName;                  // FileName
        res.TableName      = FileName.replace(/[-.]/g, '_');
        res.FileListName   = 'FileList_' + res.TableName;
        res.ChNum          = ChNum[i];                  // Channel number
        res.NumSamples     = NumSamples[i];             // Number of digitized data points
        res.ScaleFactor    = 1;                         // User defined Scale Factor - Taken 1.0 if not specified
        res.Orientation    = Orientation[i];            // Orientation of channel
        res.DateTime       = DateTime[i];               // Date&Time of the first sample
        res.Duration       = Duration[i];               // Total duration in seconds
        res.DateTime_End   = ComputeEndDateTime(res.DateTime, res.Duration); // Date&Time of the last sample in the record
        res.Lat            = Station_Latitude_Degree;   // Latitude coordinate of the sensor
        res.Long           = Station_Longitude_Degree;  // Longitude coordinate of the sensor
        res.FSamp          = FSamp[i];                  // Sampling Frequency in Hz
        res.delt           = delt[i];                   // Time interval of the digitized record in seconds
        res.TypeAndUnits   = temp.TypeAndUnit;          // Type And Unit number - refer to the list
        res.Type           = Type[i];                   // Type of sensor reading - number
        res.TypeString     = TypeString[i];             // Type of sensor reading - string
        res.Unit           = Unit[i];                   // Unit of sensor reading - number
        res.UnitString     = UnitString[i];             // Unit of sensor reading - string
        res.IntervalTypeAndUnit = IntervalTypeAndUnits[i];    // Type And Unit number - refer to the list
        res.IntervalType        = IntervalType[i];           // (0-1)
        res.IntervalTypeString  = IntervalTypeString[i];     // (Time, Spectral Period, etc.)
        res.IntervalUnit        = IntervalUnit[i];           // (0-2)
        res.IntervalUnitString  = IntervalUnitString[i];     // (Second, DateTime, etc.)
        res.InstFreq       = InstFreq;                  // Natural frequency of transducer (in Hz).
        res.InstPeriod     = InstPeriod;                // Natural period of transducer (in seconds).
        res.InstDamp       = InstDamp;                  // Damping ratio of transducer (fraction of critical).
        res.data           = data[i];                   // Digitized data
        res.time           = Time[i];                   // Time array of the digitized data

        // Calculate Statictics
        temp         = Statistics(res.data, res.ScaleFactor);
        res.Peak     = temp.Peak;
        res.Mean     = temp.Mean;
        res.RMS      = temp.RMS;
        res.Residual = res.data.at(-1) * res.ScaleFactor;

        // Add to the Main Table and Tree View
        await Add_To_Table( res );
    }

    // Update the ProgressBar per channel (Unlike other type of channels)
    await UpdateProgress(delta, "ProgressBar_Label");
    await sleep(5);

    function ExtractFloats(str) {
        let regex, match, res=[];
        regex = /[+-]?\d+(?:\.\d+)?/g;
        while (match = regex.exec(str)) {
          res.push( Number(match[0]) );
        }
        return res;
    }

    function RemoveAfter(str) {
        return str.slice( str.indexOf(":") + 1, ).trim();
    }

    function Read_Int8(dataview) {
        // This function reads the next 1-byte and returns it

        // Return an empty string it is EndOfFile
        if (EndOfFile) { return ''; }

        // Read the one-byte
        let content = dataview.getInt8( Ind );

        // Increase the Ind-byte index by 1
        Ind++;

        // Check EndOfFile
        if (Ind >= TNumByte) { EndOfFile = true; }

        // Return the content of the 1-byte read
        return content;
    }

    function FGetL(dataview) {
        // This function reads the file byte-by-byte until the next Line Feed (\n)
        // Converts the bytes read to a string and returns it
        let dat=[], content;

        while (true) {
            // Break the while-loop if EndOfFile
            if (EndOfFile) { break; }

            // Read Signed Integer of one-byte
            content = Read_Int8( dataview );

            // If the content is Carriage Return (\r), then continue reading, but do not include it in the dat[] array
            if (content == 13) { continue; }

            // If the content is a Line Feed (\n), then EndOfFile is reached, so break the while-loop and return the dat[] array
            if (content == 10) { break; }

            // Append the content to the dat[] array
            dat.push(content);
        };

        // Convert dat[] array into string and return
        return String.fromCharCode.apply(null, dat);

    }
}
//-------------------------------------------------------------------------------------------------------------
async function Read_TXT(FileName, delta, dataview) {

    // Declaration of variable
    let place, EarthquakeDate, EpicenterCoordinates, EarthquakeDepth, EarthquakeMagnitude;
    let StationID, StationCoordinates, Altitude, RecorderType, RecorderSerialNumber;
    let RecordDateTime, NumberOfSamples, delt_temp, NumOfChannel, DateTime, Lat, Long;
    let i, temp=[], temp1=[], Orientation=[], delt=[], FSamp=[], data=[], Time=[], Duration=[], ScaleFactor=[];
    let NumSamples=[],ChNum=[], Type=[], Unit=[], UnitString=[], TypeString=[];
    let IntervalTypeAndUnits=[], IntervalType=[], IntervalTypeString=[], IntervalUnit=[], IntervalUnitString=[];
    let EndOfFile=false, Ind=0, TNumByte;

    TNumByte = dataview.byteLength;

    // Skip the fist line;
    FGetL(dataview);

    place                  = RemoveAfter( FGetL(dataview) );
    EarthquakeDate         = RemoveAfter( FGetL(dataview) );
    EpicenterCoordinates   = RemoveAfter( FGetL(dataview) );
    EarthquakeDepth        = RemoveAfter( FGetL(dataview) );
    EarthquakeMagnitude    = RemoveAfter( FGetL(dataview) );
    StationID              = RemoveAfter( FGetL(dataview) );
    StationCoordinates     = RemoveAfter( FGetL(dataview) );
    Altitude               = Number( RemoveAfter( FGetL(dataview) ));
    RecorderType           = RemoveAfter( FGetL(dataview) );
    RecorderSerialNumber   = RemoveAfter( FGetL(dataview) );
    RecordDateTime         = RemoveAfter( FGetL(dataview) );
    DateTime               = DateTimeExtract(RecordDateTime)
    NumberOfSamples        = Number( RemoveAfter( FGetL(dataview) ));
    delt_temp              = Number( RemoveAfter( FGetL(dataview) ));
    NumOfChannel           = 3;

    // Skip the next 4 lines
    for (i=0; i<3; i++) { FGetL(dataview); }

    // Direction of the channels
    temp = FGetL(dataview).trim().split(' ').filter(Boolean);
    Orientation[0] = temp[0];
    Orientation[1] = temp[1];
    Orientation[2] = temp[2];

    // pre-allocation
    data[0] = [];
    data[1] = [];
    data[2] = [];

    for (i = 0; i < NumberOfSamples; i++) {
        temp = ExtractFloats( FGetL(dataview) );
        data[0].push( temp[0] );
        data[1].push( temp[1] );
        data[2].push( temp[2] );
    }

    // For each channel
    for (i=0; i<NumOfChannel; i++) {
        delt[i]        = delt_temp;
        FSamp[i]       = 1 / delt[i];
        Time[i]        = data[i].map((v,ii) => ii * delt[i]);
        Duration[i]    = Time[i].at(-1);
        ScaleFactor[i] = 1;
        NumSamples[i]  = NumberOfSamples;
        ChNum[i]       = i;
        temp           = TypeAndUnit(3);
        Type[i]        = temp.Type;
        TypeString[i]  = temp.Type_String;
        Unit[i]        = temp.Unit;
        UnitString[i]  = temp.Unit_String;

        Lat  = StationCoordinates.slice(0, StationCoordinates.indexOf('N')+1);
        Long = StationCoordinates.slice(StationCoordinates.indexOf('-')+1, StationCoordinates.indexOf('E')+1);

        temp1                  = IntervalTypeAndUnit( 1 );   // Type And Unit number - refer to the list
        IntervalTypeAndUnits[i] = temp1.IntervalTypeAndUnit;
        IntervalType[i]         = temp1.Type;                  // (0-1)
        IntervalTypeString[i]   = temp1.Type_String;           // (Time, Spectral Period, etc.)
        IntervalUnit[i]         = temp1.Unit;                  // (0-2)
        IntervalUnitString[i]   = temp1.Unit_String;           // (Second, DateTime, etc.)

        // Create new Channel Object and populate it
        let res = new Channel();

        res.FileName       = FileName;                  // FileName
        res.TableName      = FileName.replace(/[-.]/g, '_');
        res.FileListName   = 'FileList_' + res.TableName;
        res.ChNum          = ChNum[i];                  // Channel number
        res.NumSamples     = NumSamples[i];             // Number of digitized data points
        res.ScaleFactor    = 1;                         // User defined Scale Factor - Taken 1.0 if not specified
        res.Orientation    = Orientation[i];            // Orientation of channel
        res.DateTime       = DateTime;                  // Date&Time of the first sample
        res.Duration       = Duration[i];               // Total duration in seconds
        res.DateTime_End   = ComputeEndDateTime(res.DateTime, res.Duration); // Date&Time of the last sample in the record
        res.Lat            = Lat;                       // Latitude coordinate of the sensor
        res.Long           = Long;                      // Longitude coordinate of the sensor
        res.FSamp          = FSamp[i];                  // Sampling Frequency in Hz
        res.delt           = delt[i];                   // Time interval of the digitized record in seconds
        res.TypeAndUnits   = temp.TypeAndUnit;          // Type And Unit number - refer to the list
        res.Type           = Type[i];                   // Type of sensor reading - number
        res.TypeString     = TypeString[i];             // Type of sensor reading - string
        res.Unit           = Unit[i];                   // Unit of sensor reading - number
        res.UnitString     = UnitString[i];             // Unit of sensor reading - string
        res.IntervalTypeAndUnit = IntervalTypeAndUnits[i];   // Type And Unit number - refer to the list
        res.IntervalType        = IntervalType[i];           // (0-1)
        res.IntervalTypeString  = IntervalTypeString[i];     // (Time, Spectral Period, etc.)
        res.IntervalUnit        = IntervalUnit[i];           // (0-2)
        res.IntervalUnitString  = IntervalUnitString[i];     // (Second, DateTime, etc.)
        res.data           = data[i];                   // Digitized data
        res.time           = Time[i];                   // Time array of the digitized data

        // Calculate Statictics
        temp         = Statistics(res.data, res.ScaleFactor);
        res.Peak     = temp.Peak;
        res.Mean     = temp.Mean;
        res.RMS      = temp.RMS;
        res.Residual = res.data.at(-1) * res.ScaleFactor;

        // Add to the Main Table and Tree View
        await Add_To_Table( res );
    }

    // Update the ProgressBar per channel (Unlike other type of channels)
    await UpdateProgress(delta, "ProgressBar_Label");
    await sleep(5);

    function ExtractFloats(str) {
        let regex, match, res=[];
        regex = /[+-]?\d+(?:\.\d+)?/g;
        while (match = regex.exec(str)) {
          res.push( Number(match[0]) );
        }
        return res;
    }

    function DateTimeExtract(DT) {
        let Day, Month, Year, Hour, Minute, Second, DateTime;

        Day    = DT.slice(0, DT.indexOf("/")).trim();      temp  = DT.slice(DT.indexOf("/")+1, );
        Month  = temp.slice(0, temp.indexOf("/")).trim();  temp  = temp.slice(temp.indexOf("/")+1, );
        Year   = temp.slice(0, temp.indexOf(' '));         temp  = temp.slice(temp.indexOf(" ")+1, );
        Hour   = temp.slice(0, temp.indexOf(':'));         temp  = temp.slice(temp.indexOf(":")+1, );
        Minute = temp.slice(0, temp.indexOf(':'));         temp  = temp.slice(temp.indexOf(":")+1, );
        Second = temp.slice(0, temp.indexOf(' '));         temp  = temp.slice(temp.indexOf(":")+1, );
        Second = Number(Second).toFixed(3);
        DateTime = Year + "-" + ToString(Month, 2) + "-" + ToString(Day, 2) + "T" + ToString(Hour, 2) + ":" + ToString(Minute, 2) + ":" + Second;
        return DateTime;
    }

    function RemoveAfter(str) { return str.slice( str.indexOf(":") + 1, ).trim(); }

    function Read_Int8(dataview) {
        // This function reads the next 1-byte and returns it

        // Return an empty string it is EndOfFile
        if (EndOfFile) { return ''; }

        // Read the one-byte
        let content = dataview.getInt8( Ind );

        // Increase the Ind-byte index by 1
        Ind++;

        // Check EndOfFile
        if (Ind >= TNumByte) { EndOfFile = true; }

        // Return the content of the 1-byte read
        return content;
    }

    function FGetL(dataview) {
        // This function reads the file byte-by-byte until the next Line Feed (\n)
        // Converts the bytes read to a string and returns it
        let dat=[], content;

        while (true) {
            // Break the while-loop if EndOfFile
            if (EndOfFile) { break; }

            // Read Signed Integer of one-byte
            content = Read_Int8( dataview );

            // If the content is Carriage Return (\r), then continue reading, but do not include it in the dat[] array
            if (content == 13) { continue; }

            // If the content is a Line Feed (\n), then EndOfFile is reached, so break the while-loop and return the dat[] array
            if (content == 10) { break; }

            // Append the content to the dat[] array
            dat.push(content);
        };

        // Convert dat[] array into string and return
        return String.fromCharCode.apply(null, dat);

    }

    function ToString(val, n, Opt) {
        let Len, F, i;

        if (Opt==null) {Opt = true;}

        val = String(val);
        Len = val.length;
        if (Len < n) {
             F = n - Len;
            if (Opt) {for (i=0; i<F; i++) { val = '0' + val;  } }
            else     {for (i=0; i<F; i++) { val = val + '0';  } }
        }
        return val;
    }

}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
async function Read_MSD(FileName, delta, dataview) {
    // Declaration of variables
    let SequenceNumber=[], DataQualityIndicator=[], ReservedByte=[], StationIdentifierCode=[], LocationIdentifier=[], ChannelIdentifier=[], NetworkCode=[];
    let RecordStartTime_UTC_Epoch=[], RecordStartTimeISO=[], NumberOfSamples=[], SampRateFactor=[], SampRateMultiplier=[], SampleRate=[];
    let ActivityFlag=[], IO_Flag=[], DataQualityFlag=[], NumOfBlockettesFollow=[], TimeCorrection=[], OffsetBeginData=[], OffsetFirstBlockette=[];
    let Blocket=[];

    let be = false;    // Machine Format (true: Big Endian,  false: Little Endian)
    let File_Offset=0, Block_Offset=0, bN=0, i, j, Blockette_Code, BB, nibbles, x0, xn;
    let d=[], t=[];

    // Start reading each block in the file
    while (Block_Offset < dataview.byteLength) {

        // Read the Header (48-bytes) --------------------------------------------------------------------------
        ReadHeader(Block_Offset);

        // Initialize the data[] and time[] arrays
        d[bN] = new Array(NumberOfSamples[bN]).fill(NaN);
        t[bN] = new Array(NumberOfSamples[bN]).fill(NaN);

        // Create a Blocket for each block in the file
        Blocket[bN] = [];

        // Locate the pointer to the 1st Blockette in the file -------------------------------------------------
        Block_Offset += OffsetFirstBlockette[bN];

        // Read each Blockette one-by-one
        for (j = 0; j < NumOfBlockettesFollow[bN]; j++) {

            // Read Blockette_Code (2-byte)
            Blockette_Code = dataview.getUint16(Block_Offset, be);  Block_Offset += 2;

            if      (Blockette_Code === 1000) {
                // Read and store the content of the 1000-blockette (6-byte)
                Blocket[bN].push({
                    'BlocketteCode'        : Blockette_Code,
                    'OffsetNextBlockette'  : dataview.getUint16(Block_Offset,     be  ),  // 2-byte
                    'EncodingFormat'       : dataview.getUint8(Block_Offset  + 2, be  ),  // 1-byte
                    'WordOrder'            : dataview.getUint8(Block_Offset  + 3, be  ),  // 1-byte
                    'DataRecordLength'     : dataview.getUint8(Block_Offset  + 4, be  ),  // 1-byte
                    'Reserved'             : dataview.getUint8(Block_Offset  + 5, be  ),  // 1-byte
                });
                Block_Offset += 6;  // Update pointer
            }
            else if (Blockette_Code === 1001) {
                // Read and store the content of the 1001-blockette (6-byte)
                Blocket[bN].push({
                    'BlocketteCode'        : Blockette_Code,
                    'OffsetNextBlockette'  : dataview.getUint16(Block_Offset,    be  ),   // 2-byte
                    'TimingQuality'        : dataview.getUint8(Block_Offset + 2, be  ),   // 1-byte
                    'Micro_sec'            : dataview.getUint8(Block_Offset + 3, be  ),   // 1-byte
                    'Reserved'             : dataview.getUint8(Block_Offset + 4, be  ),   // 1-byte
                    'FrameCount'           : dataview.getUint8(Block_Offset + 5, be  ),   // 1-byte
                });
                Block_Offset += 6;  // Update pointer
            }
            else if (Blockette_Code === 100 ) {
                // Read and store the content of the 100-blockette (8-byte)
                Blocket[bN].push({
                    'BlocketteCode'        : Blockette_Code,
                    'OffsetNextBlockette'  : dataview.getUint16(Block_Offset,      be  ),  // 2-byte
                    'ActualSampleRate'     : dataview.getFloat32(Block_Offset + 2, be  ),  // 4-byte
                    'Flags'                : dataview.getUint8(Block_Offset   + 6, be  ),  // 1-byte
                    'Reserved'             : dataview.getUint8(Block_Offset   + 7, be  ),  // 1-byte
                });
                Block_Offset += 8;   // Update pointer
            }
            else {
                // Read and store the content of the unknown blockette (2-byte)
                Blocket[bN].push({
                    'BlocketteCode'        : Blockette_Code,
                    'OffsetNextBlockette'  : dataview.getUint16(Block_Offset, be ),  // 2-byte
                });
                Block_Offset += 2;  // Update pointer
            }
        }

        // Start reading the data streams ----------------------------------------------------------------------
        // Locate the pointer to the beginning of the data stream in the file
        Block_Offset = File_Offset + OffsetBeginData[bN];

        // Obtain the Blocket_1000
        BB = GetBlockette_1000(Blocket[bN]);

        // Read data streams based on the EncodingFormat
        if      (BB.EncodingFormat == 0) {} // Decoding format: ASCII text
        else if (BB.EncodingFormat == 1) {
            // Decoding format: 16-bit integers
            // Calculate total samples -  Each sample is 16-bit (2-byte) long
            let NumSample = (BB.DataRecordLength - OffsetBeginData[bN]) / 2;
            let Indx = 0;

            // Read each sample of 16-bits (2-byte) one-by-one
            for (i = 0; i < NumSample; i++) { d[bN][Indx] = dataview.getInt16(Block_Offset, be);  Block_Offset+=2;  Indx++; }

            // Trim data
            d[bN] = Truncate(d[bN], NumberOfSamples[bN]);
        }
        else if (BB.EncodingFormat == 2) {} // Decoding format: 24-bit integers
        else if (BB.EncodingFormat == 3) {
            // Decoding format: 32-bit integers
            // Calculate total samples -  Each sample is 32-bit (4-byte) long
            let NumSample = (BB.DataRecordLength - OffsetBeginData[bN]) / 4;
            let Indx = 0;

            // Read each sample of 32-bits (4-byte) one-by-one
            for (i = 0; i < NumSample; i++) { d[bN][Indx] = dataview.getInt32(Block_Offset, be);  Block_Offset+=4;  Indx++; }

            // Trim data
            d[bN] = Truncate(d[bN], NumberOfSamples[bN]);
        }
        else if (BB.EncodingFormat == 4) {
            // Decoding format: IEEE floating point - Float32
            // Calculate total samples -  Each sample is 32-bit (4-byte) long
            let NumSample = (BB.DataRecordLength - OffsetBeginData[bN]) / 4;
            let Indx = 0;

            // Read each sample of 32-bits (4-byte) one-by-one
            for (i = 0; i < NumSample; i++) { d[bN][Indx] = dataview.getFloat32(Block_Offset, be);  Block_Offset+=4;  Indx++; }

            // Trim data
            d[bN] = Truncate(d[bN], NumberOfSamples[bN]);
        }
        else if (BB.EncodingFormat == 5) {
            // Decoding format: IEEE double precision floating point - Float64
            // Calculate total samples -  Each sample is 64-bit (8-byte) long
            let NumSample = (BB.DataRecordLength - OffsetBeginData[bN]) / 8;
            let Indx = 0;

            // Read each sample of 64-bits (8-byte) one-by-one
            for (i = 0; i < NumSample; i++) { d[bN][Indx] = dataview.getFloat64(Block_Offset, be);  Block_Offset+=8;  Indx++; }

            // Trim data
            d[bN] = Truncate(d[bN], NumberOfSamples[bN]);
        }
        else if (BB.EncodingFormat == 10) {

            // Decoding format: STEIM-1
            // Calculate total number of frames in this block -  Each frame is 64-byte long
            let NumFrames = (BB.DataRecordLength - OffsetBeginData[bN]) / 64;

            let Indx = 0;

            // Read each frame one-by-one
            for (i = 0; i < NumFrames; i++) {

                // Get nibbles (sixteen 2-bit) : Nibbles contain the information about how the data-streams are stored
                nibbles = to32Bit(dataview.getUint32(Block_Offset, be));  Block_Offset += 4;

                // Loop over each nibbles (16 in total) - skip the first entry - Always a total of 16 nibbles
                for (j = 1; j < 16; j++) {
                    if      (nibbles[j] == 0) {
                        // This is not a time-series data - No Data
                        // If this is the first frame, then the following two 32-bit (4-byte) contains the first and last readings of the block
                        if ((i===0) && (j===1)) {
                            x0 = dataview.getInt32(Block_Offset,   be);
                            xn = dataview.getInt32(Block_Offset+4, be);
                        }
                        Block_Offset+=4;
                    }
                    else if (nibbles[j] == 1) {
                        // Four  (1-byte)
                        d[bN][Indx] = dataview.getInt8(Block_Offset, be); Block_Offset++; Indx++;
                        d[bN][Indx] = dataview.getInt8(Block_Offset, be); Block_Offset++; Indx++;
                        d[bN][Indx] = dataview.getInt8(Block_Offset, be); Block_Offset++; Indx++;
                        d[bN][Indx] = dataview.getInt8(Block_Offset, be); Block_Offset++; Indx++;
                    }
                    else if (nibbles[j] == 2) {
                        // Two  (2-byes)
                        d[bN][Indx] = dataview.getInt16(Block_Offset, be); Block_Offset+=2; Indx++;
                        d[bN][Indx] = dataview.getInt16(Block_Offset, be); Block_Offset+=2; Indx++;
                    }
                    else if (nibbles[j] == 3) {
                        // One  (4-byte)
                        d[bN][Indx] = dataview.getInt32(Block_Offset, be); Block_Offset+=4; Indx++;
                    }
                }

            }

            // The first value of each block needs to be adjusted if it is a cold-start
            if (x0 != d[bN][0]) { d[bN][0] = x0; }

            // Cumulative sum of the block over the differences
            d[bN] = CumSum( d[bN] );

        }
        else if (BB.EncodingFormat == 11) {

            // Decoding format: STEIM-2
            // Calculate total number of frames in this block -  Each frame is 64-byte long
            let NumFrames = (BB.DataRecordLength - OffsetBeginData[bN]) / 64;

            let nbit, FF, Indx = 0;

            // Read each frame one-by-one
            for (i = 0; i < NumFrames; i++) {

                // Get nibbles (sixteen 2-bit) : Nibbles contain the information about how the data streams are stored
                nibbles = to32Bit(dataview.getUint32(Block_Offset, be));  Block_Offset += 4;

                // Loop over each nibbles (16 in total) - skip the first entry - Always a total of 16 nibbles
                for (j = 1; j < 16; j++) {
                    if      (nibbles[j] == 0) {
                        // This is not a time-series data - No Data
                        // If this is the first frame, then the following two 32-bit (4-byte) contains the first and last readings of the block
                        if ((i==0) && (j==1)) {
                            x0 = dataview.getInt32(Block_Offset,   be);
                            xn = dataview.getInt32(Block_Offset+4, be);
                        }
                        Block_Offset+=4;
                    }
                    else if (nibbles[j] == 1) {
                        // Four  (1-byte)
                        d[bN][Indx] = dataview.getInt8(Block_Offset, be); Block_Offset++; Indx++;
                        d[bN][Indx] = dataview.getInt8(Block_Offset, be); Block_Offset++; Indx++;
                        d[bN][Indx] = dataview.getInt8(Block_Offset, be); Block_Offset++; Indx++;
                        d[bN][Indx] = dataview.getInt8(Block_Offset, be); Block_Offset++; Indx++;
                    }
                    else if (nibbles[j] == 2) {
                        // Type-A rules
                        nbit = to32Bit(dataview.getUint32(Block_Offset, be));
                        if (nbit[0] == 1) {
                            // 0b01 ==>  one (30-bit) difference in the frame
                            FF = toNBit(dataview.getUint32(Block_Offset, be), 30);  Block_Offset+=4;
                            d[bN][Indx] = FF[0];  Indx++;
                        }
                        else if (nbit[0] == 2) {
                            // 0b10 ==>  two (15-bit) difference in the frame
                            FF = toNBit(dataview.getUint32(Block_Offset, be), 15);  Block_Offset+=4;
                            d[bN][Indx] = FF[0];  Indx++;
                            d[bN][Indx] = FF[1];  Indx++;
                        }
                        else if (nbit[0] == 3) {
                            // 0b11 ==>  three (10-bit) difference in the frame
                            FF = toNBit(dataview.getUint32(Block_Offset, be), 10); Block_Offset+=4;
                            d[bN][Indx] = FF[0];  Indx++;
                            d[bN][Indx] = FF[1];  Indx++;
                            d[bN][Indx] = FF[2];  Indx++;
                        }
                    }
                    else if (nibbles[j] == 3) {
                        // Type-B rules
                        nbit = to32Bit(dataview.getUint32(Block_Offset, be));
                        if (nbit[0] == 0) {
                            // 0b01 ==>  five (6-bit) difference in the frame
                            FF = toNBit(dataview.getUint32(Block_Offset, be), 6); Block_Offset+=4;
                            d[bN][Indx] = FF[0];  Indx++;
                            d[bN][Indx] = FF[1];  Indx++;
                            d[bN][Indx] = FF[2];  Indx++;
                            d[bN][Indx] = FF[3];  Indx++;
                            d[bN][Indx] = FF[4];  Indx++;
                        }
                        else if (nbit[0] == 1) {
                            // 0b10 ==>  six (5-bit) difference in the frame
                            FF = toNBit(dataview.getUint32(Block_Offset, be), 5); Block_Offset+=4;
                            d[bN][Indx] = FF[0];  Indx++;
                            d[bN][Indx] = FF[1];  Indx++;
                            d[bN][Indx] = FF[2];  Indx++;
                            d[bN][Indx] = FF[3];  Indx++;
                            d[bN][Indx] = FF[4];  Indx++;
                            d[bN][Indx] = FF[5];  Indx++;
                        }
                        else if (nbit[0] == 2) {
                            // 0b11 ==>  seven (4-bit) difference in the frame
                            FF = toNBit(dataview.getUint32(Block_Offset, be), 4, 4); Block_Offset+=4;
                            d[bN][Indx] = FF[0];  Indx++;
                            d[bN][Indx] = FF[1];  Indx++;
                            d[bN][Indx] = FF[2];  Indx++;
                            d[bN][Indx] = FF[3];  Indx++;
                            d[bN][Indx] = FF[4];  Indx++;
                            d[bN][Indx] = FF[5];  Indx++;
                            d[bN][Indx] = FF[6];  Indx++;
                        }
                    }
                }

            }

            // The first value of each block needs to be adjusted if it is a cold-start
            if (x0 != d[bN][0]) { d[bN][0] = x0; }

            // Cumulative sum of the block over the differences
            d[bN] = CumSum( d[bN] );

        }

        // Create time vector for this block
        if (Int2String(ActivityFlag[bN], 8).at(6) == 0 ) {
            // Apply time correction
            RecordStartTime_UTC_Epoch[bN] += (TimeCorrection[bN] / (1000 / SampleRate[bN]));
            RecordStartTimeISO[bN]         = Epoch_To_ISO_Time( RecordStartTime_UTC_Epoch[bN] );
            t[bN] = t[bN].map((v, jj) =>  RecordStartTime_UTC_Epoch[bN] +  jj * 1000 / SampleRate[bN] );
        }

        File_Offset = Block_Offset;
        bN++;  // Block number
    }

    // Compile channel data
    MergeChannels();

    // Update the ProgressBar per channel (Unlike other type of channels)
    await UpdateProgress(delta, "ProgressBar_Label");
    await sleep(5);

    async function MergeChannels() {

        let i, j, jj, Ind, Ind1, temp, temp2, DateTime, DeltaSample, indices=[], data=[], time=[];
        let FSamp, delt;

        let IndexArray = findMatchingIndices(ChannelIdentifier, StationIdentifierCode, LocationIdentifier, NetworkCode);

        for (i=0; i < IndexArray.length; i++) {

            // Create new Channel Object and populate it
            let res = new Channel();

            indices = IndexArray[i];
            data    = [];
            time    = [];
            Ind      = indices[0];

            // this.ChannelIdentifier[i][0]; // General sampling rate and the response band of the instrument (H: High_Broad_Band, Sample Rate:  ≥ 80 to < 250, and Corner Period (sec): ≥ 10 sec.)
            // this.ChannelIdentifier[i][1]; // Instrument Code (N: Acceleration Sensor, H: High-gain seismometer, L: Low-gain seismometer)
            // this.ChannelIdentifier[i][2]; // Instrument orientation (Z, N, E: Traditional (Vertical, North-South, East-West)

            // Determine Channel Type and Unit
            if      (ChannelIdentifier[Ind][1] == 'H') { temp = TypeAndUnit(1);   }    // H: High-gain seismometer
            else if (ChannelIdentifier[Ind][1] == 'L') { temp = TypeAndUnit(1);   }    // L: Low-gain seismometer
            else if (ChannelIdentifier[Ind][1] == 'N') { temp = TypeAndUnit(1);   }    // N: Accelerometer
            else if (ChannelIdentifier[Ind][1] == 'A') { temp = TypeAndUnit(122); }    // A: Tilt
            else if (ChannelIdentifier[Ind][1] == 'D') { temp = TypeAndUnit(81);  }    // D: Pressure
            else if (ChannelIdentifier[Ind][1] == 'I') { temp = TypeAndUnit(71);  }    // I: Humidity
            else if (ChannelIdentifier[Ind][1] == 'K') { temp = TypeAndUnit(61);  }    // K: Temperature
            else if (ChannelIdentifier[Ind][1] == 'S') { temp = TypeAndUnit(31);  }    // S: Linear Strain
            else if (ChannelIdentifier[Ind][1] == 'V') { temp = TypeAndUnit(31);  }    // V: Volumetric Strain
            else if ((ChannelIdentifier[Ind][1] == 'W') && (ChannelIdentifier[Ind][2] == 'S')) { temp = TypeAndUnit(51); }    // W: Wind - Speed
            else if ((ChannelIdentifier[Ind][1] == 'W') && (ChannelIdentifier[Ind][2] == 'D')) { temp = TypeAndUnit(41); }    // D: Wind - Direction
            else { temp = TypeAndUnit(1); }     // N: Accelerometer

            // Start merging block-by-block
            for (jj=0; jj<d[Ind].length; jj++) { data.push(d[Ind][jj]); }

            // Assume Sampling rate is same across all blocks
            FSamp = SampleRate[Ind];
            delt  = 1 / FSamp;

            // DateTime
            DateTime = RecordStartTimeISO[Ind]

            for (j=1; j < indices.length; j++) {

                Ind  = indices[j];
                Ind1 = indices[j-1];

                // Calculate the number of sample points (in gap) between this block and previous block
                DeltaSample = Math.floor( Math.abs( (t[Ind].at(0) - t[Ind1].at(-1)) / (1000 / SampleRate[Ind]) - 1 ) );

                if (DeltaSample == 0) {
                    // Number of samples points (in gap) between this block and previous block is zero
                    // This means that there is no gap between these two data blocks

                    // Merge this block to this.data
                    for (jj=0; jj<d[Ind].length; jj++) { data.push(d[Ind][jj]); }

                }
                else {
                    // Number of samples points (in gap) between this block and previous block is more than 1 sample.
                    // Therefore, the gap in this.data[i] will be filled with NaN
                    // The gap does not need to be filled with NaN; instead, it can be filled with data samples interpolated between two blocks.
                    // Potential improvement for future versions.
                    for (jj=0; jj<DeltaSample; jj++) { data.push(NaN); }

                    for (jj=0; jj<d[Ind].length; jj++) { data.push(d[Ind][jj]); }
                }
            }
            // time vector
            time = data.map((v,ii) => ii * delt);

            temp2  = IntervalTypeAndUnit( 1 );   // Type And Unit number - refer to the list

            // Populate channel
            res.FileName            = FileName;
            res.TableName           = FileName.replace(/[-.]/g, '_');
            res.FileListName        = 'FileList_' + res.TableName;
            res.DateTime            = DateTime;
            res.Orientation         = ChannelIdentifier[Ind];
            res.ScaleFactor         = 1;
            res.ChNum               = i;
            res.TypeAndUnits        = temp.TypeAndUnit;
            res.Type                = temp.Type;
            res.TypeString          = temp.Type_String;
            res.Unit                = temp.Unit;
            res.UnitString          = temp.Unit_String;
            res.IntervalTypeAndUnit = temp2.IntervalTypeAndUnit;    // Type And Unit number - refer to the list
            res.IntervalType        = temp2.Type;           // (0-1)
            res.IntervalTypeString  = temp2.Type_String;     // (Time, Spectral Period, etc.)
            res.IntervalUnit        = temp2.Unit;           // (0-2)
            res.IntervalUnitString  = temp2.Unit_String;     // (Second, DateTime, etc.)
            res.Duration            = time.at(-1);
            res.DateTime_End        = ComputeEndDateTime(res.DateTime, res.Duration); // Date&Time of the last sample in the record
            res.NumSamples          = data.length;
            res.FSamp               = FSamp;
            res.delt                = delt;
            res.data                = data;
            res.time                = time;

            // Calculate Statictics
            temp         = Statistics(res.data, res.ScaleFactor);
            res.Peak     = temp.Peak;
            res.Mean     = temp.Mean;
            res.RMS      = temp.RMS;
            res.Residual = res.data.at(-1) * res.ScaleFactor;

            // Add to the Main Table and Tree View
            await Add_To_Table( res );
            await sleep(5);
        }
    }

    function findMatchingIndices(ChannelIdentifier, StationIdentifierCode, LocationIdentifier, NetworkCode) {
        const groups = {};

        for (let i = 0; i < ChannelIdentifier.length; i++) {
            const key = `${ChannelIdentifier[i]}|${StationIdentifierCode[i]}|${LocationIdentifier[i]}|${NetworkCode[i]}`;

            if (!groups[key]) {
                groups[key] = [];
            }
            groups[key].push(i);
        }

        // Return only groups that appear more than once (actual duplicates)
        return Object.values(groups).filter(indices => indices.length > 1);
    }

    function toNBit(num, nbit, skipBit) {
        // This function converts unsigned-32bit integer number to nbit-bit numbers - Skips the leading skipBit bits

        if (skipBit == null) {skipBit=2;}

        let str, i, res;

        res = [];
        str = Int2String(num, 32); // converts unsigned-32bit integer number to a binary strings of N-bit long

        // Skip the first skipBit bits
        for (i = skipBit; i+nbit < 33; i+=nbit) {
            if (str.slice(i, i+1) == '1') {
                // This is a negative number, so first Flip all bits (0 ==> 1, and 1 ==> 0), and then parse it to integer and add one.
                // This known as Two's compliment method to store negative numbers as binary
                res.push(-1 * parseInt(flipBits(str.slice(i, i+nbit)) ,2) - 1);
            }
            else {
                // This is a positive number, so no need to flipBits
                res.push(parseInt(str.slice(i+1, i+nbit), 2));
            }
        }
        return res;
    }

    function flipBits(str) { return str.split('').map((v, i) => {return (1-v).toString()}).join(''); }

    function to32Bit(num) {
        // Returns nibbles which contains the information about how the data-streams are stored
        let str, i, res;

        // Nibbles are stored in res=[] - always has a length of 16
        res = new Array(16);

        // Converts the unsigned-32bit integer number to binary string of 32-bit long
        str = Int2String(num, 32);

        for (i = 0; i < 16; i++) {
            let a1 = str[i*2];
            let a2 = str[i*2+1];
            if ((a1==0) && (a2==0)) { res[i] = 0; }
            if ((a1==0) && (a2==1)) { res[i] = 1; }
            if ((a1==1) && (a2==0)) { res[i] = 2; }
            if ((a1==1) && (a2==1)) { res[i] = 3; }
        }
        return res;
    }

    function Int2String(num, N) {
        // This function converts unsigned-32bit integer number to a binary strings of N-bit long AND returns it
        let i, n, str;
        if (num >= 0) { str = num.toString(2); } else { str = (~num).toString(2); }
        if (str.length < N) { n = N - str.length; for (i = 0; i < n; i++) { str = "0" + str; }; }
        return str
    }

    function GetBlockette_1000(A) {
        // Declaration variables
        let i, Len

        // Returns the Blockette-1000
        Len = A.length;
        for (i = 0; i < Len; i++) {
            if (A[i].BlocketteCode === 1000) {
                return {
                    'DataRecordLength'  :   Math.pow(2, A[i].DataRecordLength),  // usually 4096, but could be different.
                    'EncodingFormat'    :   A[i].EncodingFormat,
                    'WordOrder'         :   A[i].WordOrder
                }
            }
        }
        // Create a Blockette-1000 and return it.
        return {
            'DataRecordLength'  :   Math.pow(2, 12),  // 4096 - It is assumed that each block of data will be 4096 byte.
            'EncodingFormat'    :   10,
            'WordOrder'         :   1
        }
    }

    function ReadHeader(i) {
        // This function reads the header (48 bytes) starting from the i(th)-byte

        let Year, Day, Hour, Minute, Second, MilliSecond;

        SequenceNumber[bN]        = ReadNextBytes(i,  6); i += 6;  // Sequence number
        DataQualityIndicator[bN]  = ReadNextBytes(i,  1); i += 1;  // Data Header quality indicator
        ReservedByte[bN]          = ReadNextBytes(i,  1); i += 1;  // Reserved byte
        StationIdentifierCode[bN] = ReadNextBytes(i,  5); i += 5;  // Station code
        LocationIdentifier[bN]    = ReadNextBytes(i,  2); i += 2;  // Location identifier
        ChannelIdentifier[bN]     = ReadNextBytes(i,  3); i += 3;  // Channel identifier
        NetworkCode[bN]           = ReadNextBytes(i,  2); i += 2;  // Network code

        // B_TIME (10-byte)
        [Year, Day, Hour, Minute, Second, MilliSecond] = Read_B_Time(i);   i += 10; // B_Time

        RecordStartTime_UTC_Epoch[bN] = Date.UTC(Year, 0, Day, Hour, Minute, Second, (MilliSecond / 10) );
        RecordStartTimeISO[bN]        = Epoch_To_ISO_Time( RecordStartTime_UTC_Epoch[bN] );
        NumberOfSamples[bN]           = dataview.getUint16(i, be);   i += 2;  // Number of Samples
        SampRateFactor[bN]            = dataview.getInt16(i,  be);   i += 2;  // Sample Rate Factor
        SampRateMultiplier[bN]        = dataview.getInt16(i,  be);   i += 2;  // Sample rate multiplier

        // Calculate Sample Rate
        if (SampRateFactor[bN] > 0) {
            if (SampRateMultiplier[bN] >= 0) { SampleRate[bN] =      SampRateFactor[bN] * SampRateMultiplier[bN]; }
            else                             { SampleRate[bN] = -1 * SampRateFactor[bN] / SampRateMultiplier[bN]; }
        }
        else {
            if (SampRateMultiplier[bN] >= 0) { SampleRate[bN] = -1 * SampRateMultiplier[bN] / SampRateFactor[bN];   }
            else                             { SampleRate[bN] =  1 / (SampRateFactor[bN] * SampRateMultiplier[bN]); }
        }

        ActivityFlag[bN]          = dataview.getUint8(i,  be);   i += 1;  // Activity flags
        IO_Flag[bN]               = dataview.getUint8(i,  be);   i += 1;  // I/O Flags
        DataQualityFlag[bN]       = dataview.getUint8(i,  be);   i += 1;  // Data Quality flags
        NumOfBlockettesFollow[bN] = dataview.getUint8(i,  be);   i += 1;  // Number of blockettes that follow
        TimeCorrection[bN]        = dataview.getInt32(i,  be);   i += 4;  // Time correction in milliseconds
        OffsetBeginData[bN]       = dataview.getUint16(i, be);   i += 2;  // Offset to beginning of data
        OffsetFirstBlockette[bN]  = dataview.getUint16(i, be);   i += 2;  // Offset to beginning of first blockette
    }

    function ReadNextBytes(j, n) {
        // This function reads the next n-byte starting from the jj(th) byte in the file
        // Converts the bytes read to a string

        let dat=[], i, jj = j;

        for (i = 0; i < n; i++) {
            // Read the next n-byte from the file and append it to the dat[] array
            // Machine Format is not important for 8-bit
            dat.push(dataview.getInt8(jj, be));
            jj++;
        };

        // Convert dat[] array into string and return it
        return String.fromCharCode.apply(null, dat);
    }

    function Read_B_Time(Indx) {

        let Year, Day, Hour, Minute, Second, MilliSecond, Data_alignment, SwapFlag;

        // B_TIME
        Year           = dataview.getUint16(Indx    );  // 2-byte Year
        Day            = dataview.getUint16(Indx + 2);  // 2-byte Day of a year  (January 1st is 1)
        Hour           = dataview.getUint8(Indx  + 4);  // 1-byte Hour of day (0-23)
        Minute         = dataview.getUint8(Indx  + 5);  // 1-byte Minutes of day (0-59)
        Second         = dataview.getUint8(Indx  + 6);  // 1-byte Seconds of day (0—59, 60 for leap seconds)
        Data_alignment = dataview.getUint8(Indx  + 7);  // 1-byte Unused for data (required for alignment)
        MilliSecond    = dataview.getUint16(Indx + 8);  // 2-byte .0001 seconds (0—9999)

        // Automatic detection of little/big-endian encoding
        // % Automatic detection of little/big-endian encoding
        // -- by F. Beauducel, March 2014 --
        //
        // If the 2-byte day is >= 512, the file is not opened in the correct
        // endianness. If the day is 1 or 256, there is a possible byte-swap and we
        // need to check also the year; but we need to consider what is a valid year:
        // - years from 1801 to 2047 are OK (swapbytes >= 2312)
        // - years from 2048 to 2055 are OK (swapbytes <= 1800)
        // - year 2056 is ambiguous (swapbytes = 2056)
        // - years from 2057 to 2311 are OK (swapbytes >= 2312)
        // - year 1799 is ambiguous (swapbytes = 1799)
        // - year 1800 is suspicious (swapbytes = 2055)
        //
        // Thus, the only cases for which we are 'sure' there is a byte-swap, are:
        // - day >= 512
        // - (day == 1 or day == 256) and (year < 1799 or year > 2311)
        //
        // Note: in IRIS libmseed, the test is only year>2050 or year<1920.

        if ((Day >= 512) || ([1, 256].includes(Day) && ((Year > 2311) || (Year < 1799)))) {
            be           = true;
            Year         = Byte_swap16(Year);
            Day          = Byte_swap16(Day);
            MilliSecond  = Byte_swap16(MilliSecond);
        }
        else {
            be = false;
        }
        return [Year, Day, Hour, Minute, Second, MilliSecond]
    }

    function Byte_swap16(val) {
        return ((val & 0xFF) << 8)
               | ((val >> 8) & 0xFF);
    }

    function Epoch_To_ISO_Time(Ep) {
        let TempDate    = new Date( Ep );
        let Year        = TempDate.getUTCFullYear();
        let Month       = TempDate.getUTCMonth() + 1;
        let Day         = TempDate.getUTCDate();
        let Hour        = TempDate.getUTCHours();
        let Minute      = TempDate.getUTCMinutes();
        let Second      = TempDate.getUTCSeconds();
        let MilliSecond = TempDate.getUTCMilliseconds();
        return ToString(Year, 4) + '-' + ToString(Month, 2) + '-' + ToString(Day, 2) + 'T' + ToString(Hour, 2) + ':' + ToString(Minute, 2) + ':' + ToString(Second, 2) + '.' + ToString(MilliSecond, 4, false);
    }

    function ToString(val, n, Opt) {
        let Len, F, i;

        if (Opt==null) {Opt = true;}

        val = String(val);
        Len = val.length;
        if (Len < n) {
             F = n - Len;
            if (Opt) {for (i=0; i<F; i++) { val = '0' + val;  } }
            else     {for (i=0; i<F; i++) { val = val + '0';  } }
        }
        return val;
    }


}
//-------------------------------------------------------------------------------------------------------------
function ComputeEndDateTime(startDateTime, lengthSeconds) {
  const start = new Date(startDateTime);
  const end = new Date(start.getTime() + lengthSeconds * 1000);

  const fmt = (d) =>
    `${d.getFullYear()}-${String(d.getMonth() + 1).padStart(2, "0")}-${String(d.getDate()).padStart(2, "0")} ` +
    `${String(d.getHours()).padStart(2, "0")}:${String(d.getMinutes()).padStart(2, "0")}:` +
    `${String(d.getSeconds()).padStart(2, "0")}.${String(d.getMilliseconds()).padStart(3, "0")}`;

  return fmt(end).replace(" ", "T");
}
//-------------------------------------------------------------------------------------------------------------
function FileExtension(FileName) {
    // Returns the extension of a file as string
    let Ind = FileName.indexOf('.');
    if (Ind == -1) { return ""; } else { return FileName.split('.').pop(); }
}
//-------------------------------------------------------------------------------------------------------------
function IsFileUploaded(FileName) {
    // Return TRUE if the FileName is included in the ChannelList already; otherwise, returns FALSE.
    for (let j=0; j<ChannelList.length; j++) {
        if (ChannelList[j].FileName == FileName) { return true; }
    }
    return false;
}
//-------------------------------------------------------------------------------------------------------------
// File Loading Progress Bar
async function UpdateProgress(delta, Div_ID) {
    // decleration of variables
    let ProgressBar_Label, currentLevel;

    if (Div_ID == null) { Div_ID = "ProgressBar_Label"; }

    ProgressBar_Label = document.getElementById(Div_ID);

    // Get the current progress_level from data-message
    currentLevel = Number(ProgressBar_Label.dataset.message);

    // Increase the progress_lebel by delta
    currentLevel += delta; 

    if (currentLevel.toFixed(0) < 100) { 
        // Update the progressBar
        ProgressBar_Label.innerHTML = 'Loading Files (Please Wait): ' + currentLevel.toFixed(0).toString() + "%";
        ProgressBar_Label.dataset.message = currentLevel.toString();
        ProgressBar_Label.style.color = 'red';
    }
    else {
        // Update the progressBar
        ProgressBar_Label.innerHTML = "Files Loaded: 100%";
        ProgressBar_Label.dataset.message = '100';
        ProgressBar_Label.style.color = 'black';
    }
}
//-------------------------------------------------------------------------------------------------------------
// Download Files to client browser 
function StartNewWorkBook() {
    let WorkBook = XLSX.utils.book_new();
    WorkBook.Props = {
        Title:   "SGM Seismic Data Analysis",
        Subject: "SGM Seismic Data Analysis",
        Author:  "SGM Seismic Data Analysis",
    };
    return WorkBook;
}
//-------------------------------------------------------------------------------------------------------------
function StartNewWorkSheet() {
    // Start an empty worksheet
    return XLSX.utils.aoa_to_sheet([[]]);
}
//-------------------------------------------------------------------------------------------------------------
function ColumnStyle(WorkSheet, range, columnConfig) {
    if (!range || !WorkSheet || !Array.isArray(columnConfig)) return;

    // 1. Set column widths
    WorkSheet['!cols'] = columnConfig.map(col => ({ wch: col.width }));

    // 2. Apply styles column by column
    columnConfig.forEach((col, colIndex) => {
        const align = col.align || {};
        const alignment = {
            horizontal: align.horizontal || 'left',
            vertical:   align.vertical   || 'center',
            wrapText:   align.wrapText   || false
        };

        for (let row = 0; row <= range.e.r; row++) {
            const cellAddress = XLSX.utils.encode_cell({ c: colIndex, r: row });

            // Ensure the cell exists (critical for empty columns or missing headers)
            if (!WorkSheet[cellAddress]) {
                WorkSheet[cellAddress] = { t: 's', v: '' };  // empty string cell
            }

            // Initialize style object if missing
            if (!WorkSheet[cellAddress].s) {
                WorkSheet[cellAddress].s = {};
            }

            // Apply alignment to all cells
            WorkSheet[cellAddress].s.alignment = { ...alignment };

            // Apply font: bold + larger for header, normal + smaller for body
            if (row === 0) {
                // Header row
                WorkSheet[cellAddress].s.font = { bold: true, sz: 14 };
            } else {
                // Body rows
                WorkSheet[cellAddress].s.font = { bold: false, sz: 12 };
            }
        }
    });
}
//-------------------------------------------------------------------------------------------------------------
function AddDataToWorkSheet(WorkSheet, Header, Pos1, Data, Pos2, Data2=null, Pos3=null, Data3=null, Pos4=null) {
    // Add data to workSheet
    XLSX.utils.sheet_add_aoa(WorkSheet, Header, {origin: Pos1});
    XLSX.utils.sheet_add_aoa(WorkSheet, Data,   {origin: Pos2});
    if (Data2 != null) {
        XLSX.utils.sheet_add_aoa(WorkSheet, Data2,   {origin: Pos3});
    }
    if (Data3 != null) {
        XLSX.utils.sheet_add_aoa(WorkSheet, Data3,   {origin: Pos4});
    }
}
//-------------------------------------------------------------------------------------------------------------
async function DonwloadExcel_LoadDataPage() {
    
    if (IsOnHelpPage) { DisableButtons(true);   window.print();   DisableButtons(false); return; }

    // Disable all buttons
    DisableButtons(true);

    // Decleration of Variables 
    let i, j, FileName, WorkBook, WorkSheet, header, data, data2, data3, rN=1, range, columnConfig, temp, temp2, Indx, plotCount, NumOfResults;
    let CEA1, CEA2;
  
    FileName = "SDA_Results.xlsx";
  
    // Retrun if no file is uploaded
    if (ChannelList.length == 0) { 
        ProgressBar_Update('Please upload data files for processing !', 'red');
        DisableButtons(false);  
        return; 
    }

    // Disable Download Button during processing (if applicable)
    // Prevents user from triggering multiple simultaneous download operations
    document.getElementById("Header_Download").disabled = true;

    // Update the ProgressBar
    ProgressBar_Update('Please wait -- Downloading SDA_Results.xlsx File', 'red');
    await sleep(50);
    
    // Start new workbook
    WorkBook = StartNewWorkBook();

    // Create new workSheet for the list of channel with details
    WorkSheet = StartNewWorkSheet();

    // List_Of_Records
    header = ["No", "File Name", "Channel", "Duration (s)", "Sampling (Hz)", "Type", "Unit", "Scale Factor", "Date&Time", "Azimuth", "Peak", "Mean", "RMS"];
    XLSX.utils.sheet_add_aoa(WorkSheet, [header], {origin: "A1"});
    for (i=0; i<ChannelList.length; i++) {
        XLSX.utils.sheet_add_aoa(WorkSheet, [[i+1]],                                        {origin: {r:rN, c:0}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].FileName]],                    {origin: {r:rN, c:1}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].ChNum]],                       {origin: {r:rN, c:2}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].Duration]],                    {origin: {r:rN, c:3}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].FSamp]],                       {origin: {r:rN, c:4}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].TypeString]],                  {origin: {r:rN, c:5}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].UnitString]],                  {origin: {r:rN, c:6}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].ScaleFactor]],                 {origin: {r:rN, c:7}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].DateTime.replace('T','  ')]],  {origin: {r:rN, c:8}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].Orientation]],                 {origin: {r:rN, c:9}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].Peak.toPrecision(4)]],         {origin: {r:rN, c:10}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].Mean.toPrecision(4)]],         {origin: {r:rN, c:11}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].RMS.toPrecision(4)]],          {origin: {r:rN, c:12}});
        rN++;
    }

    // All column formatting
    columnConfig = [
        { width: 10, align: { horizontal: 'center', vertical: 'center' } },  // Col 0: No
        { width: 25, align: { horizontal: 'left',   vertical: 'center' } },  // Col 1: FileName
        { width: 10, align: { horizontal: 'center', vertical: 'center' } },  // Col 2: Channel Number
        { width: 15, align: { horizontal: 'right',  vertical: 'center' } },  // Col 3: Duration (s)
        { width: 17, align: { horizontal: 'right',  vertical: 'center' } },  // Col 4: Sampling (Hz)
        { width: 15, align: { horizontal: 'center', vertical: 'center' } },  // Col 5: Type
        { width: 10, align: { horizontal: 'center', vertical: 'center' } },  // Col 6: Unit
        { width: 15, align: { horizontal: 'right',  vertical: 'center' } },  // Col 7: Scale Factor
        { width: 28, align: { horizontal: 'center', vertical: 'center' } },  // Col 8: Date&Time
        { width: 10, align: { horizontal: 'center', vertical: 'center' } },  // Col 9: Azimuth
        { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 10: Peak
        { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 11: Mean
        { width: 12, align: { horizontal: 'right',  vertical: 'center' } }   // Col 12: RMS
    ];

    range = XLSX.utils.decode_range(WorkSheet['!ref']);
    ColumnStyle(WorkSheet, range, columnConfig) ;

    // Append Worksheet to Workbook
    XLSX.utils.book_append_sheet(WorkBook, WorkSheet, "List_Of_Records");

    // Loop over each channel in ChannelList[]
    for (i = 0; i < ChannelList.length; i++) {
  
        // Create new workSheet for the channel
        WorkSheet = StartNewWorkSheet();
      
        if (PageNo == 0) {
            // File name 
            FileName = "SDA_Results_Rawdata.xlsx";

            // Populate Time-RawData-FilteredData to WorkSheet
            header = ["Time (s)", "Raw Data (" + ChannelList[i].UnitString + ")"];
            data = Transpose([ChannelList[i].time,  Multiply(ChannelList[i].data, ChannelList[i].ScaleFactor)]);
            AddDataToWorkSheet(WorkSheet, [header], "A1", data, "A2");
            
            // Populate FileName, ChannelNumber, Duration, Sampling Frequency, Data Type,  Data Unit, Scale Factor, Date&Time, Orientation, Peak, Mean, and RMS
            // Statistics are calculated on raw data with scale factor applied. 
            header = [["File Name"],             ["Channel Number"],     ["Duration (s)"],          ["Scale Factor"],             ["Orientation"],              ["Sampling Frequency (Hz)"], ["Peak (" + ChannelList[i].UnitString +")"],  ["Mean (" + ChannelList[i].UnitString +")"],  ["RMS (" + ChannelList[i].UnitString +")"]];   
            data   = [[ChannelList[i].FileName], [ChannelList[i].ChNum], [ChannelList[i].Duration], [ChannelList[i].ScaleFactor], [ChannelList[i].Orientation], [ChannelList[i].FSamp],      [ChannelList[i].Peak],                        [ChannelList[i].Mean],                        [ChannelList[i].RMS]];
            AddDataToWorkSheet(WorkSheet, header, "D1", data, "E1");

            // All column formatting
            columnConfig = [
                { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 0: Time (s)
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 1: Raw Data 
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 2: EMPTY COLUMN
                { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 3: Atribute
                { width: 35, align: { horizontal: 'left',   vertical: 'center' } },  // Col 4: Value
            ];
            range = XLSX.utils.decode_range(WorkSheet['!ref']);
            ColumnStyle(WorkSheet, range, columnConfig);

        }
        else if (PageNo == 1) {
            // File name 
            FileName = "SDA_Results_Filter.xlsx";

            // Skip this channel if Filter Results are nor ready
            if (!ChannelList[i].Results.Filter.IsAnalysisCompleted) { continue; }

            // Populate Time-RawData-FilteredData to WorkSheet
            header = [["Time (s)", "Raw Data (" + ChannelList[i].UnitString + ")", "Filtered Data (" + ChannelList[i].UnitString + ")"]];
            data   = [ChannelList[i].time, Multiply(ChannelList[i].data, ChannelList[i].ScaleFactor)];
            
            // Add Filtered Data if it is exists
            // Filtered data is already scaled by Scale factor
            if (ChannelList[i].Results.Filter.IsAnalysisCompleted) {
                data.push(ChannelList[i].Results.Filter.FilteredData);
            }
            AddDataToWorkSheet(WorkSheet, header, "A1", Transpose(data), "A2");

            // Populate Filter Response to WorkSheet
            if (ChannelList[i].Results.Filter.IsAnalysisCompleted) {
                header = [["Freq (Hz)", "Magnitude", "Phase (rad)"]];
                data = Transpose([ChannelList[i].Results.Filter.f, 
                                    ChannelList[i].Results.Filter.Mag, 
                                    ChannelList[i].Results.Filter.Angle,
                                ]);
                AddDataToWorkSheet(WorkSheet, header, "D1", data, "D2");
            }
            
            // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
            header = [["File Name"],             ["Channel Number"],     ["Duration (s)"],          ["Scale Factor"],             ["Orientation"],              ["Sampling Frequency (Hz)"]];   
            data   = [[ChannelList[i].FileName], [ChannelList[i].ChNum], [ChannelList[i].Duration], [ChannelList[i].ScaleFactor], [ChannelList[i].Orientation], [ChannelList[i].FSamp]];
            AddDataToWorkSheet(WorkSheet, header, "H1", data, "I1");

            // Populate Peak, Mean, RMS
            AddDataToWorkSheet(WorkSheet, [["Raw Data (" + ChannelList[i].UnitString +")"]], "I8", [["Filtered Data (" + ChannelList[i].UnitString +")"]], "J8");
            
            header = [ ["Peak"],  ["Mean"],  ["RMS"]];
            data  = [[ChannelList[i].Peak], [ChannelList[i].Mean], [ChannelList[i].RMS]];  // Statistics of raw data and filtered data already include ScaleFactor
            data2 = [[ChannelList[i].Results.Filter.Peak], [ChannelList[i].Results.Filter.Mean], [ChannelList[i].Results.Filter.RMS]];
            AddDataToWorkSheet(WorkSheet, header, "H9", data, "I9", data2, "J9");

            // Populate Filter Settings to WorkSheet
            header = [["Filter Settings"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, "H13", data, "I13");
            

            header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
            data = [ [ChannelList[i].Results.Filter.BaselineCorrection_String],
                     [ChannelList[i].Results.Filter.FilterName_String],
                     [ChannelList[i].Results.Filter.FilterType_String],
                     [ChannelList[i].Results.Filter.FilterOrder],
                     [ChannelList[i].Results.Filter.FilterBand],
                     [ChannelList[i].Results.Filter.ZeroPhase],
                     [ChannelList[i].Results.Filter.IsStable],
                   ];
            AddDataToWorkSheet(WorkSheet, header, "H14", data, "I14");

            // Populate Filter Coefficients to WorkSheet
            header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
            data = [ChannelList[i].Results.Filter.a, ChannelList[i].Results.Filter.b];
            AddDataToWorkSheet(WorkSheet, header, "H22", data, "I22");

            // All column formatting
            columnConfig = [
                { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 0: Time (s)
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 1: Raw Data 
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 2: Filtered Data
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 3: Frequency (Hz)
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 4: Magnitude
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 5: Phase (rad)
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 5: Empty
                { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 5: Atribute
                { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 5: Value
                { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 6: Value
            ];
            range = XLSX.utils.decode_range(WorkSheet['!ref']);
            ColumnStyle(WorkSheet, range, columnConfig);

            // Make some cells bold and larger fonts size 
            ['I8', 'J8', 'H9', 'H10', 'H11', 'H13', 'H22', 'H23'].forEach(cell => {
                WorkSheet[cell].s = { 
                    ...WorkSheet[cell].s, 
                    font: { bold: true, sz: 14 } 
                };
            });


        }
        else if (PageNo == 2) {
            // File name 
            FileName = "SDA_Results_Integral.xlsx";

            // Skip this channel if Integral Results are nor ready
            if (!ChannelList[i].Results.Integral.IsAnalysisCompleted) { continue; }

            // Integral Units
            temp = TypeAndUnit(ChannelList[i].TypeAndUnits);
            temp.IntegrationUnits.FirstIntegral.Unit_String;

            // Populate Time, Rawdata, Velocity and Diplacement 
            header = ["Time (s)", "Raw Data (" + ChannelList[i].UnitString + ")", "Velocity (" + temp.IntegrationUnits.FirstIntegral.Unit_String + ")", "Displacement (" + temp.IntegrationUnits.SecondIntegral.Unit_String + ")"];
            data   = [ChannelList[i].time, Multiply(ChannelList[i].data, ChannelList[i].ScaleFactor)];

            if (ChannelList[i].Results.Integral.IsAnalysisCompleted) {
                data.push( ChannelList[i].Results.Integral.Vel  );
                data.push( ChannelList[i].Results.Integral.Disp );
            }
            AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");

            // Populate Filter Response to WorkSheet
            if (ChannelList[i].Results.Integral.IsAnalysisCompleted) {
                header = [["Freq (Hz)", "Magnitude", "Phase (rad)"]];
                data = Transpose([ ChannelList[i].Results.Integral.f, 
                                   ChannelList[i].Results.Integral.Mag, 
                                   ChannelList[i].Results.Integral.Angle,
                                ]);
                AddDataToWorkSheet(WorkSheet, header, "F1", data, "F2");
            }
            
            // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
            header = [["File Name"],             ["Channel Number"],     ["Duration (s)"],          ["Scale Factor"],             ["Orientation"],              ["Sampling Frequency (Hz)"]];   
            data   = [[ChannelList[i].FileName], [ChannelList[i].ChNum], [ChannelList[i].Duration], [ChannelList[i].ScaleFactor], [ChannelList[i].Orientation], [ChannelList[i].FSamp]];
            AddDataToWorkSheet(WorkSheet, header, "J1", data, "K1");
            
            // Populate Peak, Mean, RMS
            AddDataToWorkSheet(WorkSheet, [["Raw Data (" + ChannelList[i].UnitString + ")"]], "K8", [["Velocity (" + temp.IntegrationUnits.FirstIntegral.Unit_String + ")"]], "L8", [["Displacement (" + temp.IntegrationUnits.SecondIntegral.Unit_String + ")"]], "M8");
            
            header = [ ["Peak"],  ["Mean"],  ["RMS"], ["Residual"]];
            data  = [[ChannelList[i].Peak], [ChannelList[i].Mean], [ChannelList[i].RMS], [ChannelList[i].Residual]];  // Statistics of raw data and filtered data already include ScaleFactor
            data2 = [[ChannelList[i].Results.Integral.Peak_Vel], [ChannelList[i].Results.Integral.Mean_Vel], [ChannelList[i].Results.Integral.RMS_Vel], [ChannelList[i].Results.Integral.Residual_Vel]];
            data3 = [[ChannelList[i].Results.Integral.Peak_Disp], [ChannelList[i].Results.Integral.Mean_Disp], [ChannelList[i].Results.Integral.RMS_Disp], [ChannelList[i].Results.Integral.Residual_Disp]];
            AddDataToWorkSheet(WorkSheet, header, "J9", data, "K9", data2, "L9", data3, "M9");

            // Populate Filter Settings to WorkSheet
            header = [["Filter Settings"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, "J14", data, "K14");
            

            header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
            data = [ [ChannelList[i].Results.Integral.BaselineCorrection_String],
                     [ChannelList[i].Results.Integral.FilterName_String],
                     [ChannelList[i].Results.Integral.FilterType_String],
                     [ChannelList[i].Results.Integral.FilterOrder],
                     [ChannelList[i].Results.Integral.FilterBand],
                     [ChannelList[i].Results.Integral.ZeroPhase],
                     [ChannelList[i].Results.Integral.IsStable],
                   ];
            AddDataToWorkSheet(WorkSheet, header, "J15", data, "K15");

            // Populate Filter Coefficients to WorkSheet
            header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
            data = [ChannelList[i].Results.Integral.a, ChannelList[i].Results.Integral.b];
            AddDataToWorkSheet(WorkSheet, header, "J23", data, "K23");

            // All column formatting
            columnConfig = [
                { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 0:  Time (s)
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 1:  Raw Data 
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 2:  Velocity
                { width: 23, align: { horizontal: 'right',  vertical: 'center' } },  // Col 3:  Displacement 
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 4:  Empty
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 5:  Frequency
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 6:  Magnitude
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 7:  Phase (rad)
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 8:  Empty
                { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 9:  Atribute
                { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 10: Value
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 11: Value
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 12: Value
            ];
            range = XLSX.utils.decode_range(WorkSheet['!ref']);
            ColumnStyle(WorkSheet, range, columnConfig);
            
            // Make some cells bold and larger fonts size 
            ['K8', 'L8', 'M8', 'J9', 'J10', 'J11', 'J12', 'J14', 'J23', 'J24'].forEach(cell => {
                WorkSheet[cell].s = { 
                    ...WorkSheet[cell].s, 
                    font: { bold: true, sz: 14 } 
                };
            });

        }
        else if (PageNo == 3) {
            // File name 
            FileName = "SDA_Results_SDOF.xlsx";

            // Skip this channel if SDOF Results are nor ready
            if (!ChannelList[i].Results.SDOF.IsAnalysisCompleted) { continue; }

            // Skip if this channels has no SDOF results
            if (ChannelList[i].Results.SDOF                     == []    ) { continue; }
            if (ChannelList[i].Results.SDOF.IsAnalysisCompleted == false ) { continue; }

            // Integral Units
            temp = TypeAndUnit(ChannelList[i].TypeAndUnits);

            Indx = ChannelList[i].Results.SDOF.AnalysisMethod;

            if      (Indx == 0) { 

                // Populate Time, Rawdata, Relative Acc, totoal Acc, Velocity, Displacement, Forces and Energy
                header = [  "Time (s)", 
                            "Velocity ("               + temp.IntegrationUnits.FirstIntegral.Unit_String + ")",
                            "Displacement ("           + temp.IntegrationUnits.SecondIntegral.Unit_String + ")", 
                            "Kinetic Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Damping Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Strain Energy ("          + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Input Energy ("           + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                        ];

                data   = [  ChannelList[i].Results.SDOF.Time, 
                            ChannelList[i].Results.SDOF.Vel,
                            ChannelList[i].Results.SDOF.Disp,
                            ChannelList[i].Results.SDOF.Ek,
                            ChannelList[i].Results.SDOF.Ed,
                            ChannelList[i].Results.SDOF.Es,
                            ChannelList[i].Results.SDOF.Ei
                        ];
                        
                AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");


                // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
                header = [["File Name"],             ["Analysis Method"],                                  ["Duration (s)"],          ["Sampling Interval (sec)"]];   
                data   = [[ChannelList[i].FileName], [ChannelList[i].Results.SDOF.AnalysisMethod_string],  [ChannelList[i].Duration], [ChannelList[i].delt]      ];
                AddDataToWorkSheet(WorkSheet, header, "I1", data, "J1");
                
                // SDOF Information
                header = [["SDOF"]];
                data = "";
                AddDataToWorkSheet(WorkSheet, header, "I6", data, "J6");

                header = [["Frequency"],                   ["Damping Ratio"],                  ["Initial Displacement (cm)"],        ["Initial Velocity (cm/s)"],         ];   
                data   = [[ChannelList[i].Results.SDOF.f], [ChannelList[i].Results.SDOF.ksi],  [ChannelList[i].Results.SDOF.InDisp], [ChannelList[i].Results.SDOF.InVel], ];
                AddDataToWorkSheet(WorkSheet, header, "I7", data, "J7");


                // Populate Peak, Mean, RMS
                AddDataToWorkSheet(WorkSheet, 
                    [["Peak"  ]],   "J15",
                    [["Mean"  ]],   "K15", 
                    [["RMS"   ]],   "L15"
                );

                header = [  ["Velocity ("              + temp.IntegrationUnits.FirstIntegral.Unit_String + ")"],
                            ["Displacement ("          + temp.IntegrationUnits.SecondIntegral.Unit_String + ")"], 
                            ["Kinetic Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Damping Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Strain Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Input Energy ("          + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                        ];

                data  = [   [ChannelList[i].Results.SDOF.FiltPar.Peak_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ei],
                        ];
                data2 = [   [ChannelList[i].Results.SDOF.FiltPar.Mean_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ei],
                        ];
                data3 = [   [ChannelList[i].Results.SDOF.FiltPar.RMS_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ei],
                        ];
                AddDataToWorkSheet(WorkSheet, header, "I16", data, "J16", data2, "K16", data3, "L16");

                

                // All column formatting
                columnConfig = [
                    { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 0:  Time (s)
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 1:  Velocity
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 2:  Displacement
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 3:  Kinetic Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 4:  Damping Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 5:  Strain Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 6:  Input Energy
                    { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 7:  Empty
                    { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 8:  Value
                    { width: 30, align: { horizontal: 'right',  vertical: 'center' } },  // Col 9:  Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 10: Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 11: Value
                ];
                range = XLSX.utils.decode_range(WorkSheet['!ref']);
                ColumnStyle(WorkSheet, range, columnConfig);

                // Make some cells bold and larger fonts size 
                ['J15', 'K15', 'L15', 'I6'].forEach(cell => {
                    WorkSheet[cell].s = { 
                        ...WorkSheet[cell].s, 
                        font: { bold: true, sz: 14 } 
                    };
                });

            }
            else if (Indx == 1) { 

                // Populate Time, Rawdata, Relative Acc, totoal Acc, Velocity, Displacement, Forces and Energy
                header = [  "Time (s)", 
                            "Velocity ("               + temp.IntegrationUnits.FirstIntegral.Unit_String + ")",
                            "Displacement ("           + temp.IntegrationUnits.SecondIntegral.Unit_String + ")", 
                            "Steady-State ("           + temp.IntegrationUnits.SecondIntegral.Unit_String + ")", 
                            "Transient Displacement (" + temp.IntegrationUnits.SecondIntegral.Unit_String + ")", 
                            "Harmonic Force ("         + 'Mass • ' + ChannelList[i].UnitString + ")",
                            "Kinetic Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Damping Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Strain Energy ("          + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Input Energy ("           + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                        ];

                data   = [  ChannelList[i].time, 
                            ChannelList[i].Results.SDOF.Vel,
                            ChannelList[i].Results.SDOF.Disp,
                            ChannelList[i].Results.SDOF.Uc,
                            ChannelList[i].Results.SDOF.Up,
                            ChannelList[i].Results.SDOF.HarForce,
                            ChannelList[i].Results.SDOF.Ek,
                            ChannelList[i].Results.SDOF.Ed,
                            ChannelList[i].Results.SDOF.Es,
                            ChannelList[i].Results.SDOF.Ei
                        ];
                        
                AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");


                // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
                header = [["File Name"],             ["Analysis Method"],                                  ["Duration (s)"],          ["Sampling Interval (sec)"]];   
                data   = [[ChannelList[i].FileName], [ChannelList[i].Results.SDOF.AnalysisMethod_string],  [ChannelList[i].Duration], [ChannelList[i].delt]      ];
                AddDataToWorkSheet(WorkSheet, header, "L1", data, "M1");
                
                // SDOF Information
                header = [["SDOF"]];
                data = "";
                AddDataToWorkSheet(WorkSheet, header, "L6", data, "M6");

                header = [["Frequency"],                   ["Damping Ratio"],                  ["Initial Displacement (cm)"],        ["Initial Velocity (cm/s)"],         ["Harmonic Excitation Amplitude"],    ["Harmonic Excitation Frequency [Hz]"]   ];   
                data   = [[ChannelList[i].Results.SDOF.f], [ChannelList[i].Results.SDOF.ksi],  [ChannelList[i].Results.SDOF.InDisp], [ChannelList[i].Results.SDOF.InVel], [ChannelList[i].Results.SDOF.HarAmp], [ChannelList[i].Results.SDOF.HarF]       ];
                AddDataToWorkSheet(WorkSheet, header, "L7", data, "M7");


                // Populate Peak, Mean, RMS
                AddDataToWorkSheet(WorkSheet, 
                    [["Peak"  ]],   "M15",
                    [["Mean"  ]],   "N15", 
                    [["RMS"   ]],   "O15"
                );

                header = [  ["Velocity ("              + temp.IntegrationUnits.FirstIntegral.Unit_String + ")"],
                            ["Displacement ("          + temp.IntegrationUnits.SecondIntegral.Unit_String + ")"], 
                            ["Steady-State ("          + temp.IntegrationUnits.SecondIntegral.Unit_String + ")"], 
                            ["Transient-State ("       + temp.IntegrationUnits.SecondIntegral.Unit_String + ")"], 
                            ["Harmonic Force ("        + 'Mass • ' + ChannelList[i].UnitString + ")"],
                            ["Kinetic Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Damping Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Strain Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Input Energy ("          + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                        ];

                data  = [   [ChannelList[i].Results.SDOF.FiltPar.Peak_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Uc],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Up],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_HarFor],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ei],
                        ];
                data2 = [   [ChannelList[i].Results.SDOF.FiltPar.Mean_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Uc],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Up],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_HarFor],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ei],
                        ];
                data3 = [   [ChannelList[i].Results.SDOF.FiltPar.RMS_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Uc],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Up],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_HarFor],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ei],
                        ];
                AddDataToWorkSheet(WorkSheet, header, "L16", data, "M16", data2, "N16", data3, "O16");

                

                // All column formatting
                columnConfig = [
                    { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 0:  Time (s)
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 1:  Velocity
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 2:  Displacement
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 3:  Steady-State Displacement
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 4:  Transient Displacment 
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 5:  Harmonic Force
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 6:  Kinetic Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 7:  Damping Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 8:  Strain Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 9:  Input Energy
                    { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 10: Empty
                    { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 11: Value
                    { width: 30, align: { horizontal: 'right',  vertical: 'center' } },  // Col 12: Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 13: Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 14: Value
                ];
                range = XLSX.utils.decode_range(WorkSheet['!ref']);
                ColumnStyle(WorkSheet, range, columnConfig);

                // Make some cells bold and larger fonts size 
                ['M15', 'N15', 'O15', 'L6'].forEach(cell => {
                    WorkSheet[cell].s = { 
                        ...WorkSheet[cell].s, 
                        font: { bold: true, sz: 14 } 
                    };
                });

            }
            else if ((Indx == 2) || (Indx == 3)) { 

                // Populate Time, Rawdata, Relative Acc, totoal Acc, Velocity, Displacement, Forces and Energy
                header = [  "Time (s)", 
                            "Raw Data ("              + ChannelList[i].UnitString + ")", 
                            "Relative Acceleration (" + ChannelList[i].UnitString + ")",
                            "Total Acceleration ("    + ChannelList[i].UnitString + ")",
                            "Velocity ("              + temp.IntegrationUnits.FirstIntegral.Unit_String + ")",
                            "Displacement ("          + temp.IntegrationUnits.SecondIntegral.Unit_String + ")", 
                            "Kinetic Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Damping Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Strain Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Input Energy ("          + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                        ];

                data   = [  ChannelList[i].time, 
                            Multiply(ChannelList[i].data, ChannelList[i].ScaleFactor),
                            ChannelList[i].Results.SDOF.acc,
                            ChannelList[i].Results.SDOF.Acc,
                            ChannelList[i].Results.SDOF.Vel,
                            ChannelList[i].Results.SDOF.Disp,
                            ChannelList[i].Results.SDOF.Ek,
                            ChannelList[i].Results.SDOF.Ed,
                            ChannelList[i].Results.SDOF.Es,
                            ChannelList[i].Results.SDOF.Ei
                        ];
                AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");

                // Filter Response 
                header = [["Freq (Hz)", "Magnitude", "Phase (rad)"]];
                data = Transpose([ ChannelList[i].Results.SDOF.FiltPar.f, 
                                    ChannelList[i].Results.SDOF.FiltPar.Mag, 
                                    ChannelList[i].Results.SDOF.FiltPar.Angle,
                                    ]);
                AddDataToWorkSheet(WorkSheet, header, "L1", data, "L2");

                // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
                header = [["File Name"],             ["Analysis Method"],                                  ["Channel Number"],     ["Duration (s)"],          ["Scale Factor"],             ["Orientation"],              ["Sampling Frequency (Hz)"]];   
                data   = [[ChannelList[i].FileName], [ChannelList[i].Results.SDOF.AnalysisMethod_string],  [ChannelList[i].ChNum], [ChannelList[i].Duration], [ChannelList[i].ScaleFactor], [ChannelList[i].Orientation], [ChannelList[i].FSamp]];
                AddDataToWorkSheet(WorkSheet, header, "P1", data, "Q1");
                
                // SDOF Information
                // Populate Filter Settings to WorkSheet
                header = [["SDOF"]];
                data = "";
                AddDataToWorkSheet(WorkSheet, header, "P9", data, "Q9");

                header = [["Frequency"],                   ["Damping Ratio"],                  ["Initial Displacement (cm)"],        ["Initial Velocity (cm/s)"]           ];   
                data   = [[ChannelList[i].Results.SDOF.f], [ChannelList[i].Results.SDOF.ksi],  [ChannelList[i].Results.SDOF.InDisp], [ChannelList[i].Results.SDOF.InVel]   ];
                AddDataToWorkSheet(WorkSheet, header, "P10", data, "Q10");


                // Populate Peak, Mean, RMS
                AddDataToWorkSheet(WorkSheet, 
                    [["Peak"  ]],   "Q17",
                    [["Mean"  ]],   "R17", 
                    [["RMS"   ]],   "S17"
                );

                header = [  ["Raw Data ("              + ChannelList[i].UnitString + ")"], 
                            ["Relative Acceleration (" + ChannelList[i].UnitString + ")"],
                            ["Total Acceleration ("    + ChannelList[i].UnitString + ")"],
                            ["Velocity ("              + temp.IntegrationUnits.FirstIntegral.Unit_String + ")"],
                            ["Displacement ("          + temp.IntegrationUnits.SecondIntegral.Unit_String + ")"], 
                            ["Kinetic Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Damping Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Strain Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Input Energy ("          + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                        ];

                data  = [[ChannelList[i].Peak], 
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_acc], 
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Acc],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ei],
                        ];
                data2 = [[ChannelList[i].Mean], 
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_acc], 
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Acc],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ei],
                        ];
                data3 = [[ChannelList[i].RMS], 
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_acc], 
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Acc],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ei],
                        ];
                AddDataToWorkSheet(WorkSheet, header, "P18", data, "Q18", data2, "R18", data3, "S18");

                // Populate Filter Settings to WorkSheet
                header = [["Filter Settings"]];
                data = "";
                AddDataToWorkSheet(WorkSheet, header, "P31", data, "Q31");

                header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
                data = [[ChannelList[i].Results.SDOF.FiltPar.BaselineCorrection_String],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterName_String],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterType_String],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterOrder],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterBand],
                        [ChannelList[i].Results.SDOF.FiltPar.ZeroPhase],
                        [ChannelList[i].Results.SDOF.FiltPar.IsStable],
                    ];
                AddDataToWorkSheet(WorkSheet, header, "P32", data, "Q32");

                // Populate Filter Coefficients to WorkSheet
                header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
                data = [ChannelList[i].Results.SDOF.FiltPar.a, ChannelList[i].Results.SDOF.FiltPar.b];
                AddDataToWorkSheet(WorkSheet, header, "P40", data, "Q40");

                // All column formatting
                columnConfig = [
                    { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 0:  Time (s)
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 1:  Raw Data 
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 2:  Relative Acceleration
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 3:  Total Acceleration
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 4:  Velocity
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 5:  Displacement
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 6:  Kinetic Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 7: Damping Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 8: Strain Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 9: Input Energy
                    { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 10: Empty
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 11: Frequency
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 12: Magnitude
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 13: Phase
                    { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 14: Empty
                    { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 15: Value
                    { width: 30, align: { horizontal: 'right',  vertical: 'center' } },  // Col 16: Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 17: Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 18: Value
                ];
                range = XLSX.utils.decode_range(WorkSheet['!ref']);
                ColumnStyle(WorkSheet, range, columnConfig);

                // Make some cells bold and larger fonts size 
                ['Q17', 'R17', 'S17', 'P9', 'P31', 'P40', 'P41'].forEach(cell => {
                    WorkSheet[cell].s = { 
                        ...WorkSheet[cell].s, 
                        font: { bold: true, sz: 14 } 
                    };
                });

                
            }
            else if (Indx == 4) { 

                // Populate Time, Rawdata, Relative Acc, totoal Acc, Velocity, Displacement, Forces and Energy
                header = [  "Time (s)", 
                            "Raw Data ("              + ChannelList[i].UnitString + ")", 
                            "Relative Acceleration (" + ChannelList[i].UnitString + ")",
                            "Total Acceleration ("    + ChannelList[i].UnitString + ")",
                            "Velocity ("              + temp.IntegrationUnits.FirstIntegral.Unit_String + ")",
                            "Displacement ("          + temp.IntegrationUnits.SecondIntegral.Unit_String + ")", 
                            "Spring Force ("          + 'Mass • ' + ChannelList[i].UnitString + ")", 
                            "Damping Force ("         + 'Mass • ' + ChannelList[i].UnitString + ")", 
                            "Inertia Force ("         + 'Mass • ' + ChannelList[i].UnitString + ")", 
                            "Kinetic Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Damping Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Strain Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                            "Input Energy ("          + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )", 
                        ];

                data   = [  ChannelList[i].time, 
                            Multiply(ChannelList[i].data, ChannelList[i].ScaleFactor),
                            ChannelList[i].Results.SDOF.acc,
                            ChannelList[i].Results.SDOF.Acc,
                            ChannelList[i].Results.SDOF.Vel,
                            ChannelList[i].Results.SDOF.Disp,
                            ChannelList[i].Results.SDOF.Fs,
                            ChannelList[i].Results.SDOF.Fc,
                            ChannelList[i].Results.SDOF.Fi,
                            ChannelList[i].Results.SDOF.Ek,
                            ChannelList[i].Results.SDOF.Ed,
                            ChannelList[i].Results.SDOF.Es,
                            ChannelList[i].Results.SDOF.Ei
                        ];
                AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");

                // Filter Response 
                header = [["Freq (Hz)", "Magnitude", "Phase (rad)"]];
                data = Transpose([ ChannelList[i].Results.SDOF.FiltPar.f, 
                                    ChannelList[i].Results.SDOF.FiltPar.Mag, 
                                    ChannelList[i].Results.SDOF.FiltPar.Angle,
                                    ]);
                AddDataToWorkSheet(WorkSheet, header, "O1", data, "O2");

                // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
                header = [["File Name"],             ["Analysis Method"],                                  ["Channel Number"],     ["Duration (s)"],          ["Scale Factor"],             ["Orientation"],              ["Sampling Frequency (Hz)"]];   
                data   = [[ChannelList[i].FileName], [ChannelList[i].Results.SDOF.AnalysisMethod_string],  [ChannelList[i].ChNum], [ChannelList[i].Duration], [ChannelList[i].ScaleFactor], [ChannelList[i].Orientation], [ChannelList[i].FSamp]];
                AddDataToWorkSheet(WorkSheet, header, "S1", data, "T1");
                
                // SDOF Information
                // Populate Filter Settings to WorkSheet
                header = [["SDOF"]];
                data = "";
                AddDataToWorkSheet(WorkSheet, header, "S9", data, "T9");

                header = [["Frequency"],                   ["Damping Ratio"],                  ["Initial Displacement (cm)"],        ["Initial Velocity (cm/s)"]           ];   
                data   = [[ChannelList[i].Results.SDOF.f], [ChannelList[i].Results.SDOF.ksi],  [ChannelList[i].Results.SDOF.InDisp], [ChannelList[i].Results.SDOF.InVel]   ];
                AddDataToWorkSheet(WorkSheet, header, "S10", data, "T10");


                // Populate Peak, Mean, RMS
                AddDataToWorkSheet(WorkSheet, 
                    [["Peak"  ]],   "T17",
                    [["Mean"  ]],   "U17", 
                    [["RMS"   ]],   "V17"
                );

                header = [  ["Raw Data ("              + ChannelList[i].UnitString + ")"], 
                            ["Relative Acceleration (" + ChannelList[i].UnitString + ")"],
                            ["Total Acceleration ("    + ChannelList[i].UnitString + ")"],
                            ["Velocity ("              + temp.IntegrationUnits.FirstIntegral.Unit_String + ")"],
                            ["Displacement ("          + temp.IntegrationUnits.SecondIntegral.Unit_String + ")"], 
                            ["Spring Force ("          + 'Mass • ' + ChannelList[i].UnitString + ")"], 
                            ["Damping Force ("         + 'Mass • ' + ChannelList[i].UnitString + ")"], 
                            ["Inertia Force ("         + 'Mass • ' + ChannelList[i].UnitString + ")"], 
                            ["Kinetic Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Damping Energy ("        + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Strain Energy ("         + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                            ["Input Energy ("          + 'Mass • (' + temp.IntegrationUnits.FirstIntegral.Unit_String + ")² )"], 
                        ];

                data  = [[ChannelList[i].Peak], 
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_acc], 
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Acc],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Fs],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Fc],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Fi],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Ei],
                        ];
                data2 = [[ChannelList[i].Mean], 
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_acc], 
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Acc],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Fs],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Fc],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Fi],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Ei],
                        ];
                data3 = [[ChannelList[i].RMS], 
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_acc], 
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Acc],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Fs],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Fc],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Fi],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ek],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ed],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Es],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Ei],
                        ];
                AddDataToWorkSheet(WorkSheet, header, "S18", data, "T18", data2, "U18", data3, "V18");

                // Populate Filter Settings to WorkSheet
                header = [["Filter Settings"]];
                data = "";
                AddDataToWorkSheet(WorkSheet, header, "S31", data, "T31");

                header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
                data = [[ChannelList[i].Results.SDOF.FiltPar.BaselineCorrection_String],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterName_String],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterType_String],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterOrder],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterBand],
                        [ChannelList[i].Results.SDOF.FiltPar.ZeroPhase],
                        [ChannelList[i].Results.SDOF.FiltPar.IsStable],
                    ];
                AddDataToWorkSheet(WorkSheet, header, "S32", data, "T32");

                // Populate Filter Coefficients to WorkSheet
                header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
                data = [ChannelList[i].Results.SDOF.FiltPar.a, ChannelList[i].Results.SDOF.FiltPar.b];
                AddDataToWorkSheet(WorkSheet, header, "S40", data, "T40");

                // All column formatting
                columnConfig = [
                    { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 0:  Time (s)
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 1:  Raw Data 
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 2:  Relative Acceleration
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 3:  Total Acceleration
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 4:  Velocity
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 5:  Displacement
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 6:  Spring Force
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 7:  Damping Force
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 8:  Inertia Force
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 9:  Kinetic Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 10: Damping Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 11: Strain Energy
                    { width: 40, align: { horizontal: 'right',  vertical: 'center' } },  // Col 12: Input Energy
                    { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 13: Empty
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 14: Frequency
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 15: Magnitude
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 16: Phase
                    { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 17: Empty
                    { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 18: Value
                    { width: 30, align: { horizontal: 'right',  vertical: 'center' } },  // Col 19: Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 20: Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 21: Value
                ];
                range = XLSX.utils.decode_range(WorkSheet['!ref']);
                ColumnStyle(WorkSheet, range, columnConfig);

                // Make some cells bold and larger fonts size 
                ['T17', 'U17', 'V17', 'S9', 'S31', 'S40', 'S41'].forEach(cell => {
                    WorkSheet[cell].s = { 
                        ...WorkSheet[cell].s, 
                        font: { bold: true, sz: 14 } 
                    };
                });
                
            }
            else if (Indx == 5) { 

                // Populate Time, Rawdata, Relative Acc, totoal Acc, Velocity, Displacement, Forces and Energy
                header = [  "Time (s)", 
                            "Raw Data ("              + ChannelList[i].UnitString + ")", 
                            "Relative Acceleration (" + ChannelList[i].UnitString + ")",
                            "Total Acceleration ("    + ChannelList[i].UnitString + ")",
                            "Velocity ("              + temp.IntegrationUnits.FirstIntegral.Unit_String + ")",
                            "Displacement ("          + temp.IntegrationUnits.SecondIntegral.Unit_String + ")", 
                            "Spring Force ("          + 'Mass • ' + ChannelList[i].UnitString + ")", 
                            "Damping Force ("         + 'Mass • ' + ChannelList[i].UnitString + ")", 
                            "Inertia Force ("         + 'Mass • ' + ChannelList[i].UnitString + ")", 
                        ];

                data   = [  ChannelList[i].time, 
                            Multiply(ChannelList[i].data, ChannelList[i].ScaleFactor),
                            ChannelList[i].Results.SDOF.acc,
                            ChannelList[i].Results.SDOF.Acc,
                            ChannelList[i].Results.SDOF.Vel,
                            ChannelList[i].Results.SDOF.Disp,
                            ChannelList[i].Results.SDOF.Fs,
                            ChannelList[i].Results.SDOF.Fc,
                            ChannelList[i].Results.SDOF.Fi,
                        ];
                AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");

                // Filter Response 
                header = [["Freq (Hz)", "Magnitude", "Phase (rad)"]];
                data = Transpose([ ChannelList[i].Results.SDOF.FiltPar.f, 
                                    ChannelList[i].Results.SDOF.FiltPar.Mag, 
                                    ChannelList[i].Results.SDOF.FiltPar.Angle,
                                    ]);
                AddDataToWorkSheet(WorkSheet, header, "K1", data, "K2");

                // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
                header = [["File Name"],             ["Analysis Method"],                                  ["Channel Number"],     ["Duration (s)"],          ["Scale Factor"],             ["Orientation"],              ["Sampling Frequency (Hz)"]];   
                data   = [[ChannelList[i].FileName], [ChannelList[i].Results.SDOF.AnalysisMethod_string],  [ChannelList[i].ChNum], [ChannelList[i].Duration], [ChannelList[i].ScaleFactor], [ChannelList[i].Orientation], [ChannelList[i].FSamp]];
                AddDataToWorkSheet(WorkSheet, header, "O1", data, "P1");
                
                // SDOF Information
                header = [["SDOF"]];
                data = "";
                AddDataToWorkSheet(WorkSheet, header, "O9", data, "P9");

                header = [["Frequency"],                   ["Damping Ratio"],                  ["Initial Displacement (cm)"],        ["Initial Velocity (cm/s)"],         ["Yield Force"],                          ["Post-Yield Hardening Ratio"]     ];   
                data   = [[ChannelList[i].Results.SDOF.f], [ChannelList[i].Results.SDOF.ksi],  [ChannelList[i].Results.SDOF.InDisp], [ChannelList[i].Results.SDOF.InVel], [ChannelList[i].Results.SDOF.YieldForce], [ChannelList[i].Results.SDOF.del]  ];
                AddDataToWorkSheet(WorkSheet, header, "O10", data, "P10");


                // Populate Peak, Mean, RMS
                AddDataToWorkSheet(WorkSheet, 
                    [["Peak"  ]],   "P17",
                    [["Mean"  ]],   "Q17", 
                    [["RMS"   ]],   "R17"
                );

                header = [  ["Raw Data ("              + ChannelList[i].UnitString + ")"], 
                            ["Relative Acceleration (" + ChannelList[i].UnitString + ")"],
                            ["Total Acceleration ("    + ChannelList[i].UnitString + ")"],
                            ["Velocity ("              + temp.IntegrationUnits.FirstIntegral.Unit_String + ")"],
                            ["Displacement ("          + temp.IntegrationUnits.SecondIntegral.Unit_String + ")"], 
                            ["Spring Force ("          + 'Mass • ' + ChannelList[i].UnitString + ")"], 
                            ["Damping Force ("         + 'Mass • ' + ChannelList[i].UnitString + ")"], 
                            ["Inertia Force ("         + 'Mass • ' + ChannelList[i].UnitString + ")"], 
                        ];

                data  = [[ChannelList[i].Peak], 
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_acc], 
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Acc],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Fs],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Fc],
                            [ChannelList[i].Results.SDOF.FiltPar.Peak_Fi],
                        ];
                data2 = [[ChannelList[i].Mean], 
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_acc], 
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Acc],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Fs],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Fc],
                            [ChannelList[i].Results.SDOF.FiltPar.Mean_Fi],
                        ];
                data3 = [[ChannelList[i].RMS], 
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_acc], 
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Acc],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Vel],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Disp],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Fs],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Fc],
                            [ChannelList[i].Results.SDOF.FiltPar.RMS_Fi],
                        ];
                AddDataToWorkSheet(WorkSheet, header, "O18", data, "P18", data2, "Q18", data3, "R18");

                // Populate Filter Settings to WorkSheet
                header = [["Filter Settings"]];
                data = "";
                AddDataToWorkSheet(WorkSheet, header, "O31", data, "P31");

                header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
                data = [[ChannelList[i].Results.SDOF.FiltPar.BaselineCorrection_String],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterName_String],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterType_String],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterOrder],
                        [ChannelList[i].Results.SDOF.FiltPar.FilterBand],
                        [ChannelList[i].Results.SDOF.FiltPar.ZeroPhase],
                        [ChannelList[i].Results.SDOF.FiltPar.IsStable],
                    ];
                AddDataToWorkSheet(WorkSheet, header, "O32", data, "P32");

                // Populate Filter Coefficients to WorkSheet
                header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
                data = [ChannelList[i].Results.SDOF.FiltPar.a, ChannelList[i].Results.SDOF.FiltPar.b];
                AddDataToWorkSheet(WorkSheet, header, "O40", data, "P40");

                // All column formatting
                columnConfig = [
                    { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 0:  Time (s)
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 1:  Raw Data 
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 2:  Relative Acceleration
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 3:  Total Acceleration
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 4:  Velocity
                    { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 5:  Displacement
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 6:  Spring Force
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 7:  Damping Force
                    { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Col 8:  Inertia Force
                    { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 9: Empty
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 10: Frequency
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 11: Magnitude
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 12: Phase
                    { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 13: Empty
                    { width: 30, align: { horizontal: 'left',   vertical: 'center' } },  // Col 14: Value
                    { width: 30, align: { horizontal: 'right',  vertical: 'center' } },  // Col 15: Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 16: Value
                    { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Col 17: Value
                ];
                range = XLSX.utils.decode_range(WorkSheet['!ref']);
                ColumnStyle(WorkSheet, range, columnConfig);

                // Make some cells bold and larger fonts size 
                ['P17', 'Q17', 'R17', 'O9', 'O31', 'O40', 'O41'].forEach(cell => {
                    WorkSheet[cell].s = { 
                        ...WorkSheet[cell].s, 
                        font: { bold: true, sz: 14 } 
                    };
                });
            }

        }
        else if (PageNo == 4) {
            // File name 
            FileName = "SDA_Results_ResponseSpectrum.xlsx";

            // Skip this channel if Response Spectrum Results are nor ready
            if (!ChannelList[i].Results.ResSpec.IsAnalysisCompleted) { continue; }

            Indx = ChannelList[i].Results.ResSpec.AnalysisMethod;

            // Populate Time and Rawdata
            header = [  "Time (s)", 
                        "Raw Data ("              + ChannelList[i].UnitString + ")", 
                    ];

            data   = [  ChannelList[i].time, 
                        Multiply(ChannelList[i].data, ChannelList[i].ScaleFactor),
                    ];
            AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");
            
            // Populate period values 
            AddDataToWorkSheet(WorkSheet, [["Period (s)"]], "D1", Transpose([ChannelList[i].Results.ResSpec.T]), "D2");


            if      (Indx == 0) { plotCount = ChannelList[i].Results.ResSpec.DampingRatioCount; temp = 'ksi_';  NumOfResults = 6;  }
            else if (Indx == 1) { plotCount = ChannelList[i].Results.ResSpec.DuctilityCount;    temp = 'mu_';   NumOfResults = 4;  }
            else if (Indx == 2) { plotCount = ChannelList[i].Results.ResSpec.DuctilityCount;    temp = 'mu_';   NumOfResults = 4;  }
            
            CEA1 = "E1";
            CEA2 = "E2";

            // loop over each DampingRatio or Constant-Ductility
            for (j=0; j<plotCount; j++) {

                temp2 = temp + (j+1).toString();

                header = [  "SA ("  + temp2 + " = " + ChannelList[i].Results.ResSpec[temp2] + ") (" + ChannelList[i].UnitString + ")",  
                            "Sa ("  + temp2 + " = " + ChannelList[i].Results.ResSpec[temp2] + ") (" + ChannelList[i].UnitString + ")",
                            "SV ("  + temp2 + " = " + ChannelList[i].Results.ResSpec[temp2] + ") (" + ChannelList[i].UnitString + ")",
                            "SD ("  + temp2 + " = " + ChannelList[i].Results.ResSpec[temp2] + ") (" + ChannelList[i].UnitString + ")",
                        ];
                data   = [  ChannelList[i].Results.ResSpec.SA[j], 
                            ChannelList[i].Results.ResSpec.Sa[j],
                            ChannelList[i].Results.ResSpec.SV[j],
                            ChannelList[i].Results.ResSpec.SD[j],
                        ];

                AddDataToWorkSheet(WorkSheet, [header], CEA1, Transpose(data), CEA2); 
                CEA1 = ShiftExcellAddress(CEA1, 4, 0);
                CEA2 = ShiftExcellAddress(CEA1, 0, 1); 
                
                if (Indx==0) {
                    header = [  "SPa (" + temp2 + " = " + ChannelList[i].Results.ResSpec[temp2] + ") (" + ChannelList[i].UnitString + ")",
                                "SPv (" + temp2 + " = " + ChannelList[i].Results.ResSpec[temp2] + ") (" + ChannelList[i].UnitString + ")",
                            ];
                    data   = [  ChannelList[i].Results.ResSpec.SPa[j],
                                ChannelList[i].Results.ResSpec.SPv[j],
                            ];
                    AddDataToWorkSheet(WorkSheet, [header], CEA1, Transpose(data), CEA2); 
                    CEA1 = ShiftExcellAddress(CEA1, 2, 0);
                    CEA2 = ShiftExcellAddress(CEA1, 0, 1); 
                }

            }
            CEA1 = ShiftExcellAddress(CEA1, 1, 0);
            CEA2 = ShiftExcellAddress(CEA1, 1, 0);
            
            // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
            header = [["File Name"],             ["Analysis Method"],                                     ["Channel Number"],     ["Duration (s)"],          ["Scale Factor"],             ["Orientation"],              ["Sampling Frequency (Hz)"]];   
            data   = [[ChannelList[i].FileName], [ChannelList[i].Results.ResSpec.AnalysisMethod_string],  [ChannelList[i].ChNum], [ChannelList[i].Duration], [ChannelList[i].ScaleFactor], [ChannelList[i].Orientation], [ChannelList[i].FSamp]];
            AddDataToWorkSheet(WorkSheet, header, CEA1, data, CEA2);

            CEA1 = ShiftExcellAddress(CEA1, 0, 8);
            CEA2 = ShiftExcellAddress(CEA2, 0, 8);

            // Spectrum Information
            header = [["Spectrum"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, CEA1, data, CEA2);

            CEA1 = ShiftExcellAddress(CEA1, 0, 1);
            CEA2 = ShiftExcellAddress(CEA2, 0, 1);

            header = [["Minimum Period"],                     ["Period Step"],                          ["Maximum Period"]                      ];   
            data   = [[ChannelList[i].Results.ResSpec.T_Min], [ChannelList[i].Results.ResSpec.T_Step],  [ChannelList[i].Results.ResSpec.T_Max]  ];
            AddDataToWorkSheet(WorkSheet, header, CEA1, data, CEA2);
            
           
            if (Indx == 1) { 
                CEA1 = ShiftExcellAddress(CEA1, 0, 3);
                CEA2 = ShiftExcellAddress(CEA2, 0, 3);

                AddDataToWorkSheet(WorkSheet, [["Post Yield Hardening Ratio"]], CEA1, [[ChannelList[i].Results.ResSpec.PostYieldHard]], CEA2); 
                
                CEA1 = ShiftExcellAddress(CEA1, 0, 2);
                CEA2 = ShiftExcellAddress(CEA2, 0, 2);
            }
            else if (Indx == 2) { 
                CEA1 = ShiftExcellAddress(CEA1, 0, 3);
                CEA2 = ShiftExcellAddress(CEA2, 0, 3);

                AddDataToWorkSheet(WorkSheet, [["Post Yield Hardening Ratio"], ["Unloading Stiffness Degradation"]], CEA1, [[ChannelList[i].Results.ResSpec.PostYieldHard], [ChannelList[i].Results.ResSpec.Stiff_Deg]], CEA2); 
                
                CEA1 = ShiftExcellAddress(CEA1, 0, 3);
                CEA2 = ShiftExcellAddress(CEA2, 0, 3);
            }
            else {
                CEA1 = ShiftExcellAddress(CEA1, 0, 4);
                CEA2 = ShiftExcellAddress(CEA2, 0, 4);
            }
            

            // Populate Filter Settings to WorkSheet
            header = [["Filter Settings"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, CEA1, data, CEA2);

            CEA1 = ShiftExcellAddress(CEA1, 0, 1);
            CEA2 = ShiftExcellAddress(CEA2, 0, 1);

            header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
            data = [[ChannelList[i].Results.ResSpec.FiltPar.BaselineCorrection_String],
                    [ChannelList[i].Results.ResSpec.FiltPar.FilterName_String],
                    [ChannelList[i].Results.ResSpec.FiltPar.FilterType_String],
                    [ChannelList[i].Results.ResSpec.FiltPar.FilterOrder],
                    [ChannelList[i].Results.ResSpec.FiltPar.FilterBand],
                    [ChannelList[i].Results.ResSpec.FiltPar.ZeroPhase],
                    [ChannelList[i].Results.ResSpec.FiltPar.IsStable],
                ];
            AddDataToWorkSheet(WorkSheet, header, CEA1, data, CEA2);

            CEA1 = ShiftExcellAddress(CEA1, 0, 8);
            CEA2 = ShiftExcellAddress(CEA2, 0, 8);

            // Populate Filter Coefficients to WorkSheet
            header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
            data = [ChannelList[i].Results.ResSpec.FiltPar.a, ChannelList[i].Results.ResSpec.FiltPar.b];
            AddDataToWorkSheet(WorkSheet, header, CEA1, data, CEA2);
            

            // All column formatting
            columnConfig = [
                { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 0:  Time (s)
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Col 1:  Raw Data 
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Col 13: Empty
                { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Col 2:  Period
            ];
            for (j=0; j<plotCount*NumOfResults; j++) {
                columnConfig.push( { width: 30, align: { horizontal: 'right',  vertical: 'center' } }  );  // spectral values 
            }
            columnConfig.push( { width:  5, align: { horizontal: 'left',  vertical: 'center' } }  );  // Empty
            columnConfig.push( { width: 30, align: { horizontal: 'left',  vertical: 'center' } }  );  // Empty
            columnConfig.push( { width: 30, align: { horizontal: 'left',  vertical: 'center' } }  );  // Empty

            range = XLSX.utils.decode_range(WorkSheet['!ref']);
            ColumnStyle(WorkSheet, range, columnConfig);

            // Make some cells bold and larger fonts size
                                                         WorkSheet[CEA1].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            CEA1 = ShiftExcellAddress(CEA1, 0,  1);      WorkSheet[CEA1].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            CEA1 = ShiftExcellAddress(CEA1, 0,-10);      WorkSheet[CEA1].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            CEA1 = ShiftExcellAddress(CEA1, 0, -5);      WorkSheet[CEA1].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};

        }
        else if (PageNo == 5) { 
            // File name 
            FileName = "SDA_Results_Spectrum.xlsx";

            // Skip this channel if Spectrum Results are nor ready
            if (!ChannelList[i].Results.Spectrum.IsAnalysisCompleted) { continue; }

            // Populate Time-RawData-AriasIntensity-CAV to WorkSheet
            header = [  "Time (s)",
                        "Raw Data (" + ChannelList[i].UnitString + ")",
                    ];
            data =  [   ChannelList[i].time,  
                        Multiply(ChannelList[i].data, ChannelList[i].ScaleFactor),
                    ];
            AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");


            // Populate Time-RawData-AriasIntensity-CAV to WorkSheet
            header = [  "Frequency (Hz)",
                        "Fourier Amplitude (" + ChannelList[i].UnitString + ")",
                        "Phase (rad)",
                        "Power Spectrum (" + ChannelList[i].UnitString + ")",
                        "PSD  (" + ChannelList[i].UnitString + ")² / Hz",
                    ];
            data =  [   ChannelList[i].Results.Spectrum.Freq_vector,  
                        ChannelList[i].Results.Spectrum.FFT,
                        ChannelList[i].Results.Spectrum.FFTAngle,
                        ChannelList[i].Results.Spectrum.PowerSpectrum,
                        ChannelList[i].Results.Spectrum.PSD,
                    ];
            AddDataToWorkSheet(WorkSheet, [header], "D1", Transpose(data), "D2");

            // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
            header = [["File Name"],             ["Channel Number"],     ["Duration (s)"],          ["Scale Factor"],             ["Orientation"],              ["Sampling Frequency (Hz)"]];   
            data   = [[ChannelList[i].FileName], [ChannelList[i].ChNum], [ChannelList[i].Duration], [ChannelList[i].ScaleFactor], [ChannelList[i].Orientation], [ChannelList[i].FSamp]];
            AddDataToWorkSheet(WorkSheet, header, "J1", data, "K1");
            
            // Populate Filter Settings to WorkSheet
            header = [["Spectrum Settings"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, "J8", data, "K8");
            

            header = [["Window Length (s)"],                          ["Overlap Ratio"],                               ["Smoothing Window Type"]                                 ];   
            data   = [[ChannelList[i].Results.Spectrum.WindowLength], [ChannelList[i].Results.Spectrum.OverlapRatio],  [ChannelList[i].Results.Spectrum.SmoothingWindow_string]  ];
            AddDataToWorkSheet(WorkSheet, header, "J9", data, "K9");
            
            // Populate Filter Settings to WorkSheet
            header = [["Filter Settings"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, "J13", data, "K13");

            header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
            data = [[ChannelList[i].Results.Spectrum.FiltPar.BaselineCorrection_String],
                    [ChannelList[i].Results.Spectrum.FiltPar.FilterName_String],
                    [ChannelList[i].Results.Spectrum.FiltPar.FilterType_String],
                    [ChannelList[i].Results.Spectrum.FiltPar.FilterOrder],
                    [ChannelList[i].Results.Spectrum.FiltPar.FilterBand],
                    [ChannelList[i].Results.Spectrum.FiltPar.ZeroPhase],
                    [ChannelList[i].Results.Spectrum.FiltPar.IsStable],
                ];
            AddDataToWorkSheet(WorkSheet, header, "J14", data, "K14");

            // Populate Filter Coefficients to WorkSheet
            header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
            data = [ChannelList[i].Results.Spectrum.FiltPar.a, ChannelList[i].Results.Spectrum.FiltPar.b];
            AddDataToWorkSheet(WorkSheet, header, "J22", data, "K22");

            // All column formatting
            columnConfig = [
                { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Time (s)
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Raw Data 
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Empty
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Frequency
                { width: 30, align: { horizontal: 'right',  vertical: 'center' } },  // Fourier Amplitude 
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Phase
                { width: 30, align: { horizontal: 'right',  vertical: 'center' } },  // Power Spectrum
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Power Spectral Density
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Empty
                { width: 35, align: { horizontal: 'left',   vertical: 'center' } },  // 
                { width: 25, align: { horizontal: 'left',   vertical: 'center' } },  // 
            ];

            range = XLSX.utils.decode_range(WorkSheet['!ref']);
            ColumnStyle(WorkSheet, range, columnConfig);

            WorkSheet["J8" ].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["J13"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["J22"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["J23"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};

            

        }
        else if (PageNo == 6) {
            // File name 
            FileName = "SDA_Results_SM_Parameters.xlsx";
            
            // Skip this channel if Strong Motion Results are nor ready
            if (!ChannelList[i].Results.SM_Parameters.IsAnalysisCompleted) { continue; }

            temp = TypeAndUnit(ChannelList[i].TypeAndUnits)

            // Populate Time-RawData-AriasIntensity-CAV to WorkSheet
            header = [  "Time (s)",
                        "Raw Data (" + ChannelList[i].UnitString + ")",
                        "Arias Intensity (" + temp.IntegrationUnits.FirstIntegral.Unit_String + ")",
                        "CAV (" + temp.IntegrationUnits.FirstIntegral.Unit_String + ")",
                    ];
            data =  [   ChannelList[i].time,  
                        Multiply(ChannelList[i].data, ChannelList[i].ScaleFactor),
                        ChannelList[i].Results.SM_Parameters.AI,
                        ChannelList[i].Results.SM_Parameters.CAV,
                    ];
            AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");

            // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
            header = [["File Name"],             ["Channel Number"],     ["Duration (s)"],          ["Scale Factor"],             ["Orientation"],              ["Sampling Frequency (Hz)"]];   
            data   = [[ChannelList[i].FileName], [ChannelList[i].ChNum], [ChannelList[i].Duration], [ChannelList[i].ScaleFactor], [ChannelList[i].Orientation], [ChannelList[i].FSamp]];
            AddDataToWorkSheet(WorkSheet, header, "F1", data, "G1");

            // Populate SM Parameters to WorkSheet
            header = [["Strong Motion Parameters"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, "F8", data, "G8");

            header = [  ["Peak Ground Acceleration (" + ChannelList[i].UnitString + ")" ],
                        ["Root Mean Square ("         + ChannelList[i].UnitString + ")" ],
                        ["Maximum Arias Intensity (" + temp.IntegrationUnits.FirstIntegral.Unit_String + ")" ],
                        ["Significant Duration (s)"],
                        ["Bracketed Duration (s)"],
                        ["Effective pseudo-spectral acceleration (" + ChannelList[i].UnitString + ")" ],
                        ["Effective pseudo-spectral velocity ("     + temp.IntegrationUnits.FirstIntegral.Unit_String  + ")" ],
                        ["Housner Spectral Intensity ("             + temp.IntegrationUnits.FirstIntegral.Unit_String + ")" ],
                        ["Katayama Spectral Intensity ("            + temp.IntegrationUnits.FirstIntegral.Unit_String + ")" ],
                    ];   

            data   = [ [ChannelList[i].Peak.toPrecision(2)],
                       [ChannelList[i].RMS.toPrecision(2)],
                       [ChannelList[i].Results.SM_Parameters.AI_MaxVal.toPrecision(2)], 
                       [ChannelList[i].Results.SM_Parameters.Ts.toPrecision(2)],  
                       [ChannelList[i].Results.SM_Parameters.Td.toPrecision(2)],
                       [ChannelList[i].Results.SM_Parameters.EPa.toPrecision(2)],
                       [ChannelList[i].Results.SM_Parameters.EPv.toPrecision(2)],
                       [ChannelList[i].Results.SM_Parameters.HSI.toPrecision(2)],
                       [ChannelList[i].Results.SM_Parameters.kSI.toPrecision(2)],
                    ];
            AddDataToWorkSheet(WorkSheet, header, "F9", data, "G9");

            // Populate Filter Settings to WorkSheet
            header = [["Filter Settings"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, "F19", data, "G19");

            header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
            data = [[ChannelList[i].Results.SM_Parameters.FiltPar.BaselineCorrection_String],
                    [ChannelList[i].Results.SM_Parameters.FiltPar.FilterName_String],
                    [ChannelList[i].Results.SM_Parameters.FiltPar.FilterType_String],
                    [ChannelList[i].Results.SM_Parameters.FiltPar.FilterOrder],
                    [ChannelList[i].Results.SM_Parameters.FiltPar.FilterBand],
                    [ChannelList[i].Results.SM_Parameters.FiltPar.ZeroPhase],
                    [ChannelList[i].Results.SM_Parameters.FiltPar.IsStable],
                ];
            AddDataToWorkSheet(WorkSheet, header, "F20", data, "G20");

            // Populate Filter Coefficients to WorkSheet
            header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
            data = [ChannelList[i].Results.SM_Parameters.FiltPar.a, ChannelList[i].Results.SM_Parameters.FiltPar.b];
            AddDataToWorkSheet(WorkSheet, header, "F28", data, "G28");

            // All column formatting
            columnConfig = [
                { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Time (s)
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Raw Data 
                { width: 30, align: { horizontal: 'right',  vertical: 'center' } },  // Arias Intensity
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // CAV
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Empty
                { width: 45, align: { horizontal: 'left',   vertical: 'center' } },  // SM Parameter Key
                { width: 20, align: { horizontal: 'left',   vertical: 'center' } },  // SM Parameter Value
            ];

            range = XLSX.utils.decode_range(WorkSheet['!ref']);
            ColumnStyle(WorkSheet, range, columnConfig);

            WorkSheet["F8" ].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["F19"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["F28"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["F29"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
        }
        else if (PageNo == 7) {

            // Skip all channels expect the first one 
            // Drift is always stored in the frst channel
            if (i!=0) { continue; }

            // Skip this channel if Drift Results are not ready
            if (!ChannelList[i].Results.Drift.IsAnalysisCompleted) { continue; }

            FileName = "SDA_Results_Drift.xlsx";

            // Populate Time, Channel-1, Channel-2, Drift to WorkSheet
            header = [  "Time (s)",
                        "Channel-1 (" + ChannelList[i].Results.Drift.Unit_String + ")",
                        "Channel-2 (" + ChannelList[i].Results.Drift.Unit_String + ")",
                        "Drift ("     + ChannelList[i].Results.Drift.Unit_String + ")",
                    ];
            data =  [   ChannelList[i].Results.Drift.Drift.map((vv,ii) => ii / ChannelList[i].FSamp),
                        ChannelList[i].Results.Drift.Data[0],
                        ChannelList[i].Results.Drift.Data[1],
                        ChannelList[i].Results.Drift.Drift,
                    ];
            AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");
            
            // Populate Filter Response to WorkSheet
            header = [  "Freq (Hz)", 
                        "Magnitude", 
                        "Phase (rad)"
                    ];
            data =  [   ChannelList[i].Results.Drift.FiltPar.f, 
                        ChannelList[i].Results.Drift.FiltPar.Mag, 
                        ChannelList[i].Results.Drift.FiltPar.Angle,
                    ];
            AddDataToWorkSheet(WorkSheet, [header], "F1", Transpose(data), "F2");

            // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
            header = [["Channel-1"],
                      ["Channel-2"],
                      ["Sampling Frequency (Hz)"]
                    ];   
            data =  [["File Name"],
                      [ChannelList[i].Results.Drift.FileName[0]],
                      [ChannelList[i].Results.Drift.FileName[1]],
                      [ChannelList[i].FSamp]
                    ];
            AddDataToWorkSheet(WorkSheet, header, "J2", data, "K1");

            header = [["ChNum"]];   
            data   = [[ChannelList[i].Results.Drift.ChNum[0]], [ChannelList[i].Results.Drift.ChNum[1]]];
            AddDataToWorkSheet(WorkSheet, header, "L1", data, "L2");

            header = [["Orientation"]];   
            data   = [[ChannelList[i].Results.Drift.Azimuth[0]], [ChannelList[i].Results.Drift.Azimuth[1]]];
            AddDataToWorkSheet(WorkSheet, header, "M1", data, "M2");

            header = [["Synced Date&Time (Start)"], ["Synced Date&Time (End)"], ["Overlapped Segment Length (s)"]];
            data   = [[ChannelList[i].Results.Drift.Sync_Date_Time_Start], [ChannelList[i].Results.Drift.Sync_Date_Time_End], [ChannelList[i].Results.Drift.OverlappedSegment_Length]];
            AddDataToWorkSheet(WorkSheet, header, "J5", data, "K5");

            // Populate Filter Settings to WorkSheet
            header = [["Filter Settings"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, "J9", data, "K9");

            header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
            data = [[ChannelList[i].Results.Drift.FiltPar.BaselineCorrection_String],
                    [ChannelList[i].Results.Drift.FiltPar.FilterName_String],
                    [ChannelList[i].Results.Drift.FiltPar.FilterType_String],
                    [ChannelList[i].Results.Drift.FiltPar.FilterOrder],
                    [ChannelList[i].Results.Drift.FiltPar.FilterBand],
                    [ChannelList[i].Results.Drift.FiltPar.ZeroPhase],
                    [ChannelList[i].Results.Drift.FiltPar.IsStable],
                ];
            AddDataToWorkSheet(WorkSheet, header, "J10", data, "K10");

            // Populate Filter Coefficients to WorkSheet
            header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
            data = [ChannelList[i].Results.Drift.FiltPar.a, ChannelList[i].Results.Drift.FiltPar.b];
            AddDataToWorkSheet(WorkSheet, header, "J18", data, "K18");

            // All column formatting
            columnConfig = [
                { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Time (s)
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Channel-1
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Channel-1
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Drift
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Empty
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Frequency
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Magnitude
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Phase
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Empty
                { width: 35, align: { horizontal: 'left',   vertical: 'center' } },  // Parameters-1
                { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Parameters-2
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Parameters-3
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Parameters-4
            ];

            range = XLSX.utils.decode_range(WorkSheet['!ref']);
            ColumnStyle(WorkSheet, range, columnConfig);

            WorkSheet["J9"].s  = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["J18"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["J19"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            

        }
        else if (PageNo == 8) {

            // Skip all channels expect the first one 
            // HVSR is always stored in the frst channel
            if (i!=0) { continue; }

            // Skip this channel if HVSR Results are not ready
            if (!ChannelList[i].Results.HVSR.IsAnalysisCompleted) { continue; }

            FileName = "SDA_Results_HVSR.xlsx";

            let HVSR_Status = HVSR_Table_Check();   

            // Populate Time, Horizontal-1, Horizontal-2, Vertical to WorkSheet
            header = [  "Time (s)",
                        "Horizontal-1 (" + HVSR_Status.Unit_String[0] + ")",
                        "Horizontal-2 (" + HVSR_Status.Unit_String[0] + ")",
                        "Vertical ("     + HVSR_Status.Unit_String[0] + ")",
                    ];
            data =  [   ChannelList[i].Results.HVSR.Data[0].map((vv,ii) => ii / ChannelList[i].FSamp),
                        ChannelList[i].Results.HVSR.Data[0],
                        ChannelList[i].Results.HVSR.Data[1],
                        ChannelList[i].Results.HVSR.Data[2],
                    ];
            AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");

            // Populate frequency and H/V amplitude
            header = [  "Frequency (Hz)",
                        "HVSR_Mean Amplitude",
                        "HVSR_Std_logNormal",
                    ];
            data =  [   ChannelList[i].Results.HVSR.f,  
                        ChannelList[i].Results.HVSR.HV,
                        ChannelList[i].Results.HVSR.Std,
                    ];
            AddDataToWorkSheet(WorkSheet, [header], "F1", Transpose(data), "F2");

            
            // Populate FileName, ChannelNumber, Duration (s), Data Type, Data Unit, Scale Factor, Orientation and Sampling Frequency
            header = [["Horizontal-1"],
                      ["Horizontal-2"],
                      ["Vertical"],
                      ["Window length (samples)"],
                      ["Overlap length (samples)"],
                      ["Sampling Frequency (Hz)"]
                    ];   
            data   = [["File Name"],
                      [ChannelList[i].Results.HVSR.FileName[0]],
                      [ChannelList[i].Results.HVSR.FileName[1]],
                      [ChannelList[i].Results.HVSR.FileName[2]],
                      [ChannelList[i].Results.HVSR.RWL],
                      [ChannelList[i].Results.HVSR.OVS],
                      [ChannelList[i].FSamp]
                    ];
            AddDataToWorkSheet(WorkSheet, header, "J2", data, "K1");

            header = [["ChNum"]];   
            data   = [[ChannelList[i].Results.HVSR.ChNum[0]], [ChannelList[i].Results.HVSR.ChNum[1]], [ChannelList[i].Results.HVSR.ChNum[2]]];
            AddDataToWorkSheet(WorkSheet, header, "L1", data, "L2");

            header = [["Orientation"]];   
            data   = [[ChannelList[i].Results.HVSR.Azimuth[0]], [ChannelList[i].Results.HVSR.Azimuth[1]], [ChannelList[i].Results.HVSR.Azimuth[2]]];
            AddDataToWorkSheet(WorkSheet, header, "M1", data, "M2");
            
            header = [["Synced Date&Time (Start)"], ["Synced Date&Time (End)"], ["Overlapped Segment Length (s)"]];
            data   = [[ChannelList[i].Results.HVSR.Sync_Date_Time_Start], [ChannelList[i].Results.HVSR.Sync_Date_Time_End], [ChannelList[i].Results.HVSR.OverlappedSegment_Length]];
            AddDataToWorkSheet(WorkSheet, header, "J8", data, "K8");


            // Populate Filter Settings to WorkSheet
            header = [["Filter Settings"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, "J12", data, "K12");

            header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
            data = [[ChannelList[i].Results.HVSR.FiltPar.BaselineCorrection_String],
                    [ChannelList[i].Results.HVSR.FiltPar.FilterName_String],
                    [ChannelList[i].Results.HVSR.FiltPar.FilterType_String],
                    [ChannelList[i].Results.HVSR.FiltPar.FilterOrder],
                    [ChannelList[i].Results.HVSR.FiltPar.FilterBand],
                    [ChannelList[i].Results.HVSR.FiltPar.ZeroPhase],
                    [ChannelList[i].Results.HVSR.FiltPar.IsStable],
                ];
            AddDataToWorkSheet(WorkSheet, header, "J13", data, "K13");

            // Populate Filter Coefficients to WorkSheet
            header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
            data = [ChannelList[i].Results.HVSR.FiltPar.a, ChannelList[i].Results.HVSR.FiltPar.b];
            AddDataToWorkSheet(WorkSheet, header, "J21", data, "K21");

            header = [["Upper bound = HVSR * EXP(Std)"],["Lower bound = HVSR * EXP(-Std)"]];
            data = "";
            AddDataToWorkSheet(WorkSheet, header, "J25", data, "K25");

            // All column formatting
            columnConfig = [
                { width: 12, align: { horizontal: 'right',  vertical: 'center' } },  // Time (s)
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Horizontal-1
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Horizontal-1
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Vertical
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Empty
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Frequency
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // H/V Amplitude
                { width: 25, align: { horizontal: 'right',  vertical: 'center' } },  // Standard Deviation
                { width:  5, align: { horizontal: 'right',  vertical: 'center' } },  // Empty
                { width: 35, align: { horizontal: 'left',   vertical: 'center' } },  // Parameters-1
                { width: 35, align: { horizontal: 'right',  vertical: 'center' } },  // Parameters-2
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Parameters-3
                { width: 20, align: { horizontal: 'right',  vertical: 'center' } },  // Parameters-4
            ];

            range = XLSX.utils.decode_range(WorkSheet['!ref']);
            ColumnStyle(WorkSheet, range, columnConfig);

            WorkSheet["J12"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["J21"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["J22"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["J25"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};
            WorkSheet["J26"].s = {font: { bold: true, sz: 14 }, align: { horizontal: 'left',  vertical: 'center' }};

        }

        // Add the workSheet to WorkBook
        XLSX.utils.book_append_sheet(WorkBook, WorkSheet, (i+1).toString());
        
        // Update the progressBar
        // Calculate the Increment for the Progress_Bar 
        ProgressBar_Update('Please wait -- ' + (i/ChannelList.length*100).toFixed(0).toString() + '% Done', 'red');
        await sleep(50);
    }


    // 
    Add_Disclaimer();

  
    // Inform user on progress (Wait for screen update)
    ProgressBar_Update("Please wait -- Writing XLSX File", 'red');
    await sleep(50);
  
    // Write workBook to XLSX file
    XLSX.writeFile(WorkBook, FileName, { compression: true });
  
    // Inform user on progress
    ProgressBar_Update('Downloading Done! ---> ' + FileName.toString().trim(), 'black');

    // Enable all buttons
    DisableButtons(false);

    // Helper functions
    function ShiftExcellAddress(row_excel, col, row) {
        // Strip optional $ anchors (e.g. "$B$3" → "B3")
        const clean = row_excel.replace(/\$/g, "");

        // Split into column letters and row number
        const match = clean.match(/^([A-Z]+)(\d+)$/i);
        if (!match) throw new Error(`Invalid cell address: "${row_excel}"`);

        const [, colStr, rowStr] = match;

        // --- Column: letters → 1-based number ---
        const colToNum = (letters) =>
            letters.toUpperCase().split("").reduce((acc, ch) => acc * 26 + (ch.charCodeAt(0) - 64), 0);

        // --- Column: 1-based number → letters ---
        const numToCol = (n) => {
            let s = "";
            while (n > 0) {
            const rem = (n - 1) % 26;
            s = String.fromCharCode(65 + rem) + s;
            n = Math.floor((n - 1) / 26);
            }
            return s;
        };

        const newCol = colToNum(colStr) + col;
        const newRow = parseInt(rowStr, 10) + row;

        if (newCol < 1) throw new Error(`Column shift goes out of bounds (result: ${newCol})`);
        if (newRow < 1) throw new Error(`Row shift goes out of bounds (result: ${newRow})`);

        return numToCol(newCol) + newRow;
    }
    function Add_Disclaimer() {
        
        // Create new workSheet for Disclaimer
        WorkSheet = StartNewWorkSheet();

        AddDataToWorkSheet(WorkSheet, [["Disclaimer"]], "A1", [[]], "A2");

        const paragraphsList = Array.from(document.querySelectorAll('.Disclaimer_Body p')).map(p => p.textContent.trim());
        let currentRow = 2;

        for (const tt of paragraphsList) {
            let headerCell = "A" + currentRow;
            let dataCell   = "A" + (currentRow + 1);
            AddDataToWorkSheet(WorkSheet, [[tt]], headerCell, [[]], dataCell);
            currentRow += 1;
        }

        // Initialize the !merges array if it doesn't exist yet
        if (!WorkSheet['!merges']) WorkSheet['!merges'] = [];

        // Initialize the !rows array if it doesn't exist yet
        if (!WorkSheet['!rows']) WorkSheet['!rows'] = [];

        // Helper function to set heights safely (Row 1 is index 0, Row 3 is index 2, etc.)
        const setRowHeight = (rowIndex, height) => {
            WorkSheet['!rows'][rowIndex] = { hpt: height };
        };

        // Set row heights in points
        setRowHeight(0, 40);  // Row 1
        setRowHeight(1, 70);  // Row 3
        setRowHeight(2, 95);  // Row 5
        setRowHeight(3, 70);  // Row 7
        setRowHeight(4, 55);  // Row 9
        setRowHeight(5, 30); // Row 11
        
        // Add the merge range for A1 to E1 (0-indexed: Row 0, Col 0 to Row 0, Col 4)
        WorkSheet['!merges'].push(XLSX.utils.decode_range('A1:E1'));
        WorkSheet['!merges'].push(XLSX.utils.decode_range('A2:E2'));
        WorkSheet['!merges'].push(XLSX.utils.decode_range('A3:E3'));
        WorkSheet['!merges'].push(XLSX.utils.decode_range('A4:E4'));
        WorkSheet['!merges'].push(XLSX.utils.decode_range('A5:E5'));
        WorkSheet['!merges'].push(XLSX.utils.decode_range('A6:E6'));

        if (!WorkSheet['!cols']) WorkSheet['!cols'] = [];
        WorkSheet['!cols'][0] = { wch: 25 }; // Width of Column A
        WorkSheet['!cols'][1] = { wch: 20 }; // Width of Column B
        WorkSheet['!cols'][2] = { wch: 20 }; // Width of Column C
        WorkSheet['!cols'][3] = { wch: 20 }; // Width of Column D
        WorkSheet['!cols'][4] = { wch: 25 }; // Width of Column E

        // Apply styles to A1 (the top-left cell of the merge holds the value and style)
        WorkSheet["A1"].s  = {font: { bold: true, sz: 22 }, alignment: { horizontal: 'left', vertical: 'center', wrapText: true } };
        WorkSheet["A2"].s  = {font: { bold: true, sz: 14 }, alignment: { horizontal: 'left', vertical: 'center', wrapText: true } };
        WorkSheet["A3"].s  = {font: { bold: true, sz: 14 }, alignment: { horizontal: 'left', vertical: 'center', wrapText: true } };
        WorkSheet["A4"].s  = {font: { bold: true, sz: 14 }, alignment: { horizontal: 'left', vertical: 'center', wrapText: true } };
        WorkSheet["A5"].s  = {font: { bold: true, sz: 14 }, alignment: { horizontal: 'left', vertical: 'center', wrapText: true } };
        WorkSheet["A6"].s  = {font: { bold: true, sz: 14 }, alignment: { horizontal: 'left', vertical: 'center', wrapText: true } };

        // Add the Disclaimer workSheet to WorkBook
        XLSX.utils.book_append_sheet(WorkBook, WorkSheet, 'Disclaimer');

    }


}
//-----------------------------------------------------------------------------------------------
async function LoadServerSampleFilesDirectly() {
    var sampleFiles = [
        './static/data/CE12951.V1',
        './static/data/CE13873.V1',
        './static/data/20111023104043_0401.txt',
        './static/data/20111023104118_1302.txt',
        './static/data/TWF025-n.7000m9g4.HLN.10.V1c',
        './static/data/TWF025-n.7000m9g4.HLZ.10.V1c',
        './static/data/TWF025-n.7000m9g4.HLE.10.V1c',
    ];

    try {
        if (typeof ProgressBar_Update === "function") {
            ProgressBar_Update('Fetching sample files from server...', 'blue');
        }
        
        var dataTransfer = new DataTransfer();

        for (var i = 0; i < sampleFiles.length; i++) {
            var url = sampleFiles[i];
            var fileName = url.substring(url.lastIndexOf('/') + 1);
            
            var response = await fetch(url);
            if (!response.ok) {
                throw new Error(`HTTP network error! status: ${response.status} for file: ${fileName}`);
            }
            
            var fileBlob = await response.blob();
            var fileObject = new File([fileBlob], fileName, { type: fileBlob.type });
            dataTransfer.items.add(fileObject);
        }

        // 🟢 FIXED: Target the actual hidden input element that holds the listener
        var fileInputElement = document.getElementById("LoadInputFiles");
        
        if (fileInputElement) {
            fileInputElement.files = dataTransfer.files;

            // This now successfully wakes up Load_Files(e) 
            var changeEvent = new Event('change', { bubbles: true });
            fileInputElement.dispatchEvent(changeEvent);
            
            if (typeof ProgressBar_Update === "function") {
                ProgressBar_Update('Sample files successfully uploaded to browser!', 'green');
            }
        } else {
            console.error("Target element id='LoadInputFiles' was not found in the DOM.");
        }

    } catch (error) {
        console.error("Failed to upload server files:", error);
        if (typeof ProgressBar_Update === "function") {
            ProgressBar_Update('Error loading sample files from server.', 'red');
        }
    }
}
//-----------------------------------------------------------------------------------------------
function Init_DragAndDrop_Upload() {

    let Overlay, DragCounter;

    Overlay     = document.getElementById('DragDrop_Overlay');
    DragCounter = 0;

    function Is_File_Drag(e) {
        return e.dataTransfer && Array.from(e.dataTransfer.types || []).includes('Files');
    }

    window.addEventListener('dragover', (e) => {
        if (!Is_File_Drag(e)) { return; }
        e.preventDefault();
        e.dataTransfer.dropEffect = 'copy';
    });

    window.addEventListener('dragenter', (e) => {
        if (!Is_File_Drag(e)) { return; }
        e.preventDefault();
        DragCounter++;
        Overlay.style.display = 'flex';
    });

    window.addEventListener('dragleave', (e) => {
        if (!Is_File_Drag(e)) { return; }
        e.preventDefault();
        DragCounter--;
        if (DragCounter <= 0) {
            DragCounter = 0;
            Overlay.style.display = 'none';
        }
    });

    window.addEventListener('drop', (e) => {
        if (!Is_File_Drag(e)) { return; }
        e.preventDefault();

        DragCounter = 0;
        Overlay.style.display = 'none';

        let DroppedFiles = e.dataTransfer.files;
        if (DroppedFiles && DroppedFiles.length > 0) {
            Load_Files({ target: { files: DroppedFiles } });
            AnalysisMenu_Selection(document.getElementById("MainMenu_LoadData"));
        }
    });
}
//-----------------------------------------------------------------------------------------------

