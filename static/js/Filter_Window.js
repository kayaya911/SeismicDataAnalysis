// stricter parsing and error handling
"use strict";


// -------------------------------------------------------------------------
// Calculation Button ------------------------------------------------------
async function Calculation() {
    if      (PageNo == 1) { await Filter__Data();    }
    else if (PageNo == 2) { await Integrate__Data(); }
}
// -------------------------------------------------------------------------
// Filter Calculation ------------------------------------------------------
function BaselineCorrection_Change() {
    // Declaration of variables
    let Indx = document.getElementById('BaselineCorrection').selectedIndex;
}
function FilterName_Change() {
    // Declaration of variables
    let FilterTable, Indx;

    // Get the index number of the FilterName
    Indx = document.getElementById('FilterName').selectedIndex;

    // Disable the FilterTable rows if no filtering is selected 
    if (Indx === 0){
        document.getElementById('FilterType').disabled     = true;
        document.getElementById('FilterOrder').disabled     = true;
        document.getElementById('F1').disabled              = true;
        document.getElementById('F2').disabled              = true;
        document.getElementById('RippleSize').disabled      = true;
        document.getElementById('ZeroPhase').disabled       = true;
    }
    else {
        document.getElementById('FilterType').disabled      = false;
        document.getElementById('FilterOrder').disabled     = false;
        document.getElementById('F1').disabled              = false;
        document.getElementById('F2').disabled              = false;
        document.getElementById('RippleSize').disabled      = false;
        document.getElementById('ZeroPhase').disabled       = false;
    }

    // Display Maximum Ripplesize if Chebyshev-filter type is selected
    FilterTable = document.getElementById('FilterTable');
    if (Indx === 2) {
        FilterTable.rows[3].style.display = "table-row";
    }
    else {
        FilterTable.rows[3].style.display = "none";
    }
}
function FilterType_Change() {
    // Declaration of variables
    let FilterTable, Indx, F1, F2;

    // Get all filter paraeters
    Indx = document.getElementById('FilterType').selectedIndex;

    // Show F2 if Band-Pass is selected; otherwise, hide F2.
    FilterTable = document.getElementById('FilterTable');
    if      (Indx === 0) { FilterTable.rows[2].style.display = "none";      } 
    else if (Indx === 1) { FilterTable.rows[2].style.display = "none";      }
    else if (Indx === 2) { 
        F1 = Number(document.getElementById('F1').value);
        F2 = Number(document.getElementById('F2').value);
        if (F1>F2) { 
            document.getElementById('F1').value = F2;
            document.getElementById('F2').value = F1;
        }
        else if (F1 == F2) {
            document.getElementById('F2').value = F2 + 0.5;
        }
        FilterTable.rows[2].style.display = "table-row";
    }
}
function FilterOrder_Change() {
    // Filter order must be between zero and 8   (0 < FilterOrder <=8)
    // No changes are made in case of invalid entry.

    // Declaration of variables
    let x = document.getElementById('FilterOrder');

    if (Number(x.value) <= 0 || Number(x.value) > 12) {
        // Invalid entry.
        document.getElementById('FilterOrder').value = x.oldValue;
    }
    else {
        // Valid entry.
        document.getElementById('FilterOrder').value = String(Math.ceil(Number(x.value)));
    }
}
function F1_Change() {
    // Declaration of variables
    let x          = document.getElementById('F1');
    let F2         = Number(document.getElementById('F2').value);
    let FilterType = document.getElementById('FilterType').selectedIndex;

    // Make sure that the F1 cut-off frequency is between 0 < F1 < F2 < FNyquist
    if (Number(x.value) <= 0) {
        document.getElementById('F1').value = x.oldValue;
    }
    else if ((FilterType == 2) && (Number(x.value) >= F2)) {
        // Band-Pass filter 
        document.getElementById('F1').value = x.oldValue;
    }
    else {
        document.getElementById('F1').value        = String(Number(x.value));
        document.getElementById('F1').defaultValue = String(Number(x.value));
    }
}
function F2_Change() {

    // Declaration of variables
    let F1 = Number(document.getElementById('F1').value);
    let x  = document.getElementById('F2');

    // Make sure that the F2 cut-off corner frequency is between:  0 < F1 < F2 < FNyquist
    if (Number(x.value) <= 0) {
        document.getElementById('F2').value = x.oldValue;
    }
    else if ( Number(x.value) <= F1 ) {
        document.getElementById('F2').value = x.oldValue;
    }
    else {
        document.getElementById('F2').value        = String(Number(x.value));
        document.getElementById('F2').defaultValue = String(Number(x.value));
    }
}
function RippleSize_Change() {
    // Ripple Size (dB) is used for Chebyshev tyype of filter only
    // Ripple size must be greater than zero.
    let x = document.getElementById('RippleSize');

    if (Number(x.value) <= 0) {
        // Invalid entry
        document.getElementById('RippleSize').value = x.defaultValue;
    }
    else {
        // Valid entry
        document.getElementById('RippleSize').value        = String(Number(x.value));
        document.getElementById('RippleSize').defaultValue = String(Number(x.value));
    }
}
function ZeroPhase_Change() {

    if (document.getElementById("ZeroPhase").checked){
        // Do something
    }
    else {
        // Do something else
    }
}
async function GetFilterParameters() {

    // Declaration of variables
    let BaselineCorrection, FilterName, FilterType, FilterOrder, F1, F2, FilterBand, RippleSize, ZeroPhase, ShowHideFilterResp;
    let BaselineCorrection_String, FilterName_String, FilterType_String;

    // Filter parameters 
    BaselineCorrection = document.getElementById('BaselineCorrection').selectedIndex;
    FilterName         = document.getElementById('FilterName').selectedIndex;
    FilterType         = document.getElementById('FilterType').selectedIndex;
    FilterOrder        = Number(document.getElementById('FilterOrder').value);
    F1                 = Number(document.getElementById('F1').value);
    F2                 = Number(document.getElementById('F2').value);
    RippleSize         = Number(document.getElementById('RippleSize').value);
    ZeroPhase          = document.getElementById('ZeroPhase').checked;  

    if      (BaselineCorrection == 0 ) { BaselineCorrection_String = "None";        }
    else if (BaselineCorrection == 1 ) { BaselineCorrection_String = "Remove Mean"; }
    else if (BaselineCorrection == 2 ) { BaselineCorrection_String = "Linear";      }
    else if (BaselineCorrection == 3 ) { BaselineCorrection_String = "Quadratic";   }
    else if (BaselineCorrection == 4 ) { BaselineCorrection_String = "Qubic";       }

    if      (FilterName == 0) { FilterName_String = "None";        }
    else if (FilterName == 1) { FilterName_String = "Butterworth"; }
    else if (FilterName == 2) { FilterName_String = "Chebyshev";   }
    else if (FilterName == 3) { FilterName_String = "Bessel";      }

    if      (FilterType == 0) { FilterType_String = "Lowpass";   FilterBand = "[" + F1.toString() + " Hz]";  }
    else if (FilterType == 1) { FilterType_String = "Highpass";  FilterBand = "[" + F1.toString() + " Hz]";  }
    else if (FilterType == 2) { FilterType_String = "Bandpass";  FilterBand = "[" + F1.toString() + " - " + F2.toString() + " Hz]";    }

    if (PageNo == 1) {

        // Return for Filtering 
        return {
            FilteredData                : undefined,
            BaselineCorrection          : BaselineCorrection,
            BaselineCorrection_String   : BaselineCorrection_String,
            FilterName                  : FilterName,
            FilterName_String           : FilterName_String,
            FilterType                  : FilterType,
            FilterType_String           : FilterType_String,
            FilterOrder                 : FilterOrder,
            F1                          : F1,
            F2                          : F2,
            FilterBand                  : FilterBand,
            RippleSize                  : RippleSize,
            ZeroPhase                   : ZeroPhase,
            ZeroPhase                   : ZeroPhase,
            Peak                        : undefined,
            Mean                        : undefined,
            RMS                         : undefined,
            a                           : undefined,
            b                           : undefined,
            zf                          : undefined,
            H                           : undefined,
            f                           : undefined,
            Mag                         : undefined,
            Angle                       : undefined,
            IsStable                    : undefined,                // True if filter is estimated to be stable; false otherwise. 
            ErrorMessage                : undefined,                // Error Message
            CutOff_Freq_Check           : undefined,                // Check if cuf-off frequencies are correct with respect to Nyquist frequency
            Unit_Ind                    : 0,                        // Index of the Unit on Graph 
        }
    }
    else if (PageNo == 2) {
        
        // Return for Integration
        return {
            Vel                         : undefined,
            Disp                        : undefined,
            BaselineCorrection          : BaselineCorrection,
            BaselineCorrection_String   : BaselineCorrection_String,
            FilterName                  : FilterName,
            FilterName_String           : FilterName_String,
            FilterType                  : FilterType,
            FilterType_String           : FilterType_String,
            FilterOrder                 : FilterOrder,
            F1                          : F1,
            F2                          : F2,
            FilterBand                  : FilterBand,
            RippleSize                  : RippleSize,
            ZeroPhase                   : ZeroPhase,
            ZeroPhase                   : ZeroPhase,
            Peak_vel                    : undefined,
            Mean_Vel                    : undefined,
            RMS_Vel                     : undefined,
            Peak_Disp                   : undefined,
            Mean_Disp                   : undefined,
            RMS_Disp                    : undefined,
            a                           : undefined,
            b                           : undefined,
            zf                          : undefined,
            H                           : undefined,
            f                           : undefined,
            Mag                         : undefined,
            Angle                       : undefined,
            IsStable                    : undefined,                // True if filter is estimated to be stable; false otherwise. 
            ErrorMessage                : undefined,                // Error Message
            CutOff_Freq_Check           : undefined,                // Check if cuf-off frequencies are correct with respect to Nyquist frequency 
            Unit_Ind                    : 0,                        // Index of the Unit on Graph -- Acceleration
            Unit_Ind_v                  : 0,                        // Index of the Unit on Graph -- Velocity 
            Unit_Ind_d                  : 0,                        // Index of the Unit on Graph -- Displacement 
        }
    }
    
}
async function Filter_Cut_Off_Check(Channel, FiltPar) {

    // Decleration of variables 
    let FilterType  = Number(FiltPar.FilterType);
    let F1          = Number(FiltPar.F1);
    let F2          = Number(FiltPar.F2);
    let Nyquist     = Channel.FSamp / 2;

    // Check Filter Cut-Off Frequency
    if (FilterType == 0) { 
        // Low-Pass Filter
        if ( F1 >= Nyquist) { 
            FiltPar.CutOff_Freq_Check = false;
            FiltPar.ErrorMessage      = "Check Filter Seetings! (F1)";     
            return FiltPar; 
        } 
    }
    else if (FilterType == 1) { 
        // High-Pass Filter 
        if ( F1 >= Nyquist) { 
            FiltPar.CutOff_Freq_Check = false;
            FiltPar.ErrorMessage      = "Check Filter Seetings! (F1)";     
            return FiltPar; 
        } 
    }
    else if (FilterType == 2) { 
        // Band-Pass Filter 
        if ((F1 >= Nyquist) || (F2 >= Nyquist) || (F1 >= F2)) { 
            FiltPar.CutOff_Freq_Check = false;
            FiltPar.ErrorMessage      = "Check Filter Seetings! (F1- F2)"; 
            return FiltPar; 
        } 
    }
    FiltPar.CutOff_Freq_Check = true;
    return FiltPar;
}
async function Filter_Stability(Channel, FiltPar) {

    // Decleration of variables 
    let FR, FRZ;
    let FilterType  = Number(FiltPar.FilterType);
    let FilterOrder = Number(FiltPar.FilterOrder);
    let F1          = Number(FiltPar.F1);
    let F2          = Number(FiltPar.F2);

    // Calculate Filter (a) and (b) coefficients 
    if (FilterType == 0) {
        // Lowpass filter 
        FR =  Butterworth_LowPass(FilterOrder, F1, Channel.FSamp);
    }
    else if (FilterType == 1) {
        // Highpass filter 
        FR = Butterworth_HighPass(FilterOrder, F1, Channel.FSamp);
    }
    else if (FilterType == 2) {
        // Bandpass filter 
        FR = Butterworth_BandPass(FilterOrder, F1, F2, Channel.FSamp);
    }

    FiltPar.a     = FR.a;
    FiltPar.b     = FR.b;
    FiltPar.zf    = FR.zf;

    // Check Filter stability using (a) and (b) coefficients
    if ((FiltPar.a != undefined) || (FiltPar.b != undefined)) { 
        
        // Filter is Stable 
        FiltPar.IsStable     = true; 

        // Calculate Filter Response 
        FRZ           = Freqz(FR.b, FR.a, 512, Channel.FSamp, "ONESIDED");
        FiltPar.H     = FRZ.H;
        FiltPar.f     = FRZ.f;
        FiltPar.Mag   = Abs(FRZ.H);
        FiltPar.Angle = Angle(FRZ.H);
        return FiltPar;
    } 
    else { 

        // Filter is Unstable
        FiltPar.IsStable = false; 
        FiltPar.ErrorMessage = "Unstable Filter!"; 
        return FiltPar;
    }
}
async function Filter__Data() {

    let FiltPar = [];

    // Disable CALCULATE Button
    document.getElementById("Calculate_Button").disabled = true;

    // Inform user on progress
    document.getElementById("ProgressBar_Label").innerHTML = "<b>Filtering Selected Channels</b>";

    // Reset the ProgressBar
    await ResetProgressBar("ProgressBar_LoadData");
    await sleep(5);

    let delta = 100 / ChannelList.length;

    // Loop over each channel
    for (let i=0; i<ChannelList.length; i++) {

        // Skip this Channel if it is not selcted for analysis
        if (!ChannelList[i].Selected) { await IncreaseProgressBar( delta, "ProgressBar_LoadData"); await sleep(5); continue; }

        // Get new set of filter parameters   
        FiltPar = await GetFilterParameters();

        // Store the default filter-settings
        ChannelList[i].Results.Filter = FiltPar;
        
        // Apply Baseline Correction
        if (!FiltPar.BaselineCorrection == 0) {  FiltPar = await BaseLineCorrection(ChannelList[i], FiltPar);  }
        
        // Filter the Channel 
        FiltPar = await FilterChannelData(ChannelList[i], FiltPar);
        
        // Calculate Statistics on filtered data
        if ((!FiltPar.BaselineCorrection == 0) || FiltPar.IsStable) { FiltPar = await Statistics_Channel(FiltPar); };
        
        // Store final Results 
        ChannelList[i].Results.Filter = FiltPar;

        // Update Graph
        await Plot_Update(i);
        
        // Update ProgressBar
        await IncreaseProgressBar( delta, "ProgressBar_LoadData")
        await sleep(5);
    }
    // Enable CALCULATE Button
    document.getElementById("Calculate_Button").disabled = false;


    async function BaseLineCorrection(Channel, FiltPar) {

        // Get the RawData for the channel and apply Scaling Factor 
        let FilteredData = Mult(Channel.data, Channel.ScaleFactor);

        // Apply Detrend 
        if      (FiltPar.BaselineCorrection == 1) { FilteredData = Detrend(FilteredData, 0);  }
        else if (FiltPar.BaselineCorrection == 2) { FilteredData = Detrend(FilteredData, 1);  }
        else if (FiltPar.BaselineCorrection == 3) { FilteredData = Detrend(FilteredData, 2);  }
        else if (FiltPar.BaselineCorrection == 4) { FilteredData = Detrend(FilteredData, 3);  }

        // Store the Baseline Corrected Data 
        FiltPar.FilteredData = FilteredData;

        return FiltPar;
    }
    async function FilterChannelData(Channel, FiltPar) {

        if (FiltPar.FilterName != 0) { 

            // Check for Cut-off frequencies with respect to Nyquist Frequency 
            FiltPar = await Filter_Cut_Off_Check(Channel, FiltPar);

            // Check if FiltPar.CutOff_Freq_Check is TRUE
            if (FiltPar.CutOff_Freq_Check) {

                // Check if Filter is stable
                FiltPar = await Filter_Stability(Channel, FiltPar); 

                // Filter the Channel 
                if (FiltPar.IsStable) {  FiltPar = await Filter_RawData(Channel, FiltPar); };
            }
        }
        return FiltPar;
    }
    async function Filter_RawData(Channel, FiltPar) {

        // Check if FilteredData exist, if not get the RawData to filter 
        if (FiltPar.FilteredData == undefined) { FiltPar.FilteredData = Mult(Channel.data, Channel.ScaleFactor); }

        let FiltResp;

        // Apply filter
        if (FiltPar.ZeroPhase) { FiltResp = FiltFilt(FiltPar.b,   FiltPar.a,   FiltPar.FilteredData);   }
        else                   { FiltResp = Filter(  FiltPar.b,   FiltPar.a,   FiltPar.FilteredData);   }

        FiltPar.FilteredData = FiltResp.y;
        FiltPar.zf           = FiltResp.zf;

        return FiltPar;
    }
    async function Statistics_Channel(FiltPar) {

        // Update Statistics of Filtered Data
        FiltPar.Peak      =  Max( Abs(FiltPar.FilteredData) ).val;
        FiltPar.Mean      =  Mean( FiltPar.FilteredData );
        FiltPar.RMS       =  Rms( FiltPar.FilteredData );

        return FiltPar;
    }    
}
// -------------------------------------------------------------------------
// Integration -------------------------------------------------------------
async function AccVellDisp_Select(a) {
    
    // Decleration of variables 
    let i, flag = false;

    if (a.id == "Int_Acceleration") {

        if (document.getElementById("Int_Acceleration").checked) { flag = true; }

        document.getElementById("Int_Acceleration").checked  = true;
        document.getElementById("Int_Velocity").checked      = false;
        document.getElementById("Int_Displacement").checked  = false;
        
        // Update Graphs
        if (flag) { for (i=0; i<ChannelList.length; i++) { await Plot_Update(i);} }
        
    }
    else if (a.id == "Int_Velocity") {

        if (document.getElementById("Int_Velocity").checked) { flag = true; }

        document.getElementById("Int_Acceleration").checked  = false;
        document.getElementById("Int_Velocity").checked      = true;
        document.getElementById("Int_Displacement").checked  = false;

        // Update Graphs
        if (flag) { for (i=0; i<ChannelList.length; i++) { await Plot_Update(i);} }
    }
    else if (a.id == "Int_Displacement") {

        if (document.getElementById("Int_Displacement").checked) { flag = true; }

        document.getElementById("Int_Acceleration").checked  = false;
        document.getElementById("Int_Velocity").checked      = false;
        document.getElementById("Int_Displacement").checked  = true;

        // Update Graphs
        if (flag) { for (i=0; i<ChannelList.length; i++) { await Plot_Update(i);} }
    }
}
async function Integrate__Data() {
    
    let i, IntegralPar=[], RawData;

    // Disable CALCULATE Button
    document.getElementById("Calculate_Button").disabled = true;

    // Inform user on progress
    document.getElementById("ProgressBar_Label").innerHTML = "<b>Integrating Selected Channels</b>";

    // Reset the ProgressBar
    await ResetProgressBar("ProgressBar_LoadData");
    await sleep(5);

    let delta = 100 / ChannelList.length;

    // Loop over each channel
    for (i=0; i<ChannelList.length; i++) {

        // Skip this Channel if it is not selected for analysis
        if ((!ChannelList[i].Selected) || (ChannelList[i].Type !=0)) {await IncreaseProgressBar( delta, "ProgressBar_LoadData"); await Plot_Update(i); await sleep(5); continue; }

        // Get new set of filter parameters   
        IntegralPar = await GetFilterParameters();

        // Check Filter Stability
        IntegralPar = await Filter_Check(ChannelList[i], IntegralPar);

        // Continue Integration if Filter is Stable 
        if (IntegralPar.IsStable) {

            // Get the RawData for the channel
            RawData = Mult(ChannelList[i].data, ChannelList[i].ScaleFactor);

            // Integrate RawData to estimate Velocity
            IntegralPar.Vel = Filter([0.5/ChannelList[i].FSamp,   0.5/ChannelList[i].FSamp], [1, -1],   await BLC_Filter(RawData, IntegralPar)).y;
            IntegralPar.Vel = await BLC_Filter(IntegralPar.Vel, IntegralPar);

            // Integrate Velocity to estimate Displacement
            IntegralPar.Disp = Filter([0.5/ChannelList[i].FSamp,  0.5/ChannelList[i].FSamp], [1, -1],   IntegralPar.Vel).y;

            // Calculate Statistics on Vel and Disp
            IntegralPar = await Statistics_VelDisp(IntegralPar);
        }
        
        // Store final Results 
        ChannelList[i].Results.Integral = IntegralPar;

        // Update Graph
        await Plot_Update(i);
    
        // Update ProgressBar
        await IncreaseProgressBar( delta, "ProgressBar_LoadData")
        await sleep(5);
    }
    // Enable CALCULATE Button
    document.getElementById("Calculate_Button").disabled = false;


    async function Filter_Check(Channel, IntegralPar) {
        if (IntegralPar.FilterName != 0) { 

            // Check for Cut-off frequencies with respect to Nyquist Frequency 
            IntegralPar = await Filter_Cut_Off_Check(Channel, IntegralPar);

            // Check if IntegralPar.CutOff_Freq_Check is TRUE
            if (IntegralPar.CutOff_Freq_Check) {

                // Check if Filter is stable
                IntegralPar = await Filter_Stability(Channel, IntegralPar); 

            }
        }
        return IntegralPar;
    }
    async function BLC_Filter(Data, IntegralPar) {
        
        let FilteredData = Copy(Data);

        // Apply Baseline Correction
        if (!IntegralPar.BaselineCorrection == 0) {  
            // Apply Detrend 
            if      (IntegralPar.BaselineCorrection == 1) { FilteredData = Detrend(FilteredData, 0);  }
            else if (IntegralPar.BaselineCorrection == 2) { FilteredData = Detrend(FilteredData, 1);  }
            else if (IntegralPar.BaselineCorrection == 3) { FilteredData = Detrend(FilteredData, 2);  }
            else if (IntegralPar.BaselineCorrection == 4) { FilteredData = Detrend(FilteredData, 3);  }
        }

        // Apply filtering
        if (IntegralPar.ZeroPhase) { FilteredData = FiltFilt( IntegralPar.b,   IntegralPar.a,   FilteredData).y;   }
        else                       { FilteredData =   Filter( IntegralPar.b,   IntegralPar.a,   FilteredData).y;   }

        return FilteredData;

    }
    async function Statistics_VelDisp(IntegralPar) {

        // Update Statistics for Vel and Disp
        IntegralPar.Peak_Vel   =  Max( Abs(IntegralPar.Vel) ).val;
        IntegralPar.Mean_Vel   =  Mean( IntegralPar.Vel );
        IntegralPar.RMS_Vel    =  Rms( IntegralPar.Vel );

        IntegralPar.Peak_Disp  =  Max( Abs(IntegralPar.Disp) ).val;
        IntegralPar.Mean_Disp  =  Mean( IntegralPar.Disp );
        IntegralPar.RMS_Disp   =  Rms( IntegralPar.Disp );

        return IntegralPar;
    }    
}
// -------------------------------------------------------------------------
// PWave -------------------------------------------------------------------
function GetPWaveParameters() {

    // Declaration of variables
    let PWave_Method;

    // PWave parameters 
    PWave_Method = document.getElementById('PWave_Method').selectedIndex;
    BinNumber = Number(document.getElementById('PWave_HistBinNum').value);
    Period = Number(document.getElementById('PWave_Period').value);
    Damping = Number(document.getElementById('PWave_ksi').value);

    return {
        PWave_Method    : PWave_Method,
        numBins         : BinNumber,
        Period          : Period,
        Damping         : Damping,
        PWaveArrival    : undefined,
        Ed              : undefined,
        RawData         : undefined,
        FilteredData    : undefined,
        level           : undefined,
    }

}
function Enable_Disable_PWave(val) {
    document.getElementById('PWave_Method').disabled = val; 
    document.getElementById('PWave_HistBinNum').disabled = val; 
    document.getElementById('PWave_Period').disabled = val; 
    document.getElementById('PWave_ksi').disabled = val; 
}
