// stricter parsing and error handling
"use strict";


//----------------------------------------------------------------------------------------------------------------------
// Collapse all functions in Visual Studio Code 
// Ctrl + K  and  Ctrl + 0

//----------------------------------------------------------------------------------------------------------------------
// PWave detection
function PhasePicker_Hist(Data, delt, Type, ksi, Tn, numBins) {
    
    // Decleration of variables 
    let Result, Res, Td, Ed, Mode, locs, PWaveArrival=undefined;

    // Check agruments
    if (ksi == null)     { ksi = 0.60; }
    if (Tn == null)      { if (delt<=0.01) {Tn = 100;} else {Tn = 10;} } 
    if (numBins == null) { numBins = Math.round(2/delt); }

    // Normalize and Filter Data
    Result = Normalize_Filter(Data, Type);
    
    // Calculate the reponse of SDOF to FilteredData.
    Res = PieceWiseLin(Tn, ksi, delt, 0, 0, Result.Data);

    // Calculate the rate of change (power) of damping energy dissipated for a unit mass
    Td = Tn / Math.sqrt(1-ksi*ksi);
    Ed  = Mult(Mult_El(Res.Vel, Res.Vel), 2*ksi*(2*Math.PI/Td));

    // Calculate the low-state level, which is the mode of the highest bin in lower-histogram
    // The low-state corresponds to the PWave-phase
    Mode = StateLevel(Ed, numBins);

    // Find the index numbers where Ed > Mode
    locs = Find(CheckRelation(Ed, ">", Mode));

    if (locs[0] != undefined) {
        // locs[] is not an empty array
        if (locs[0] > 0) {
            // locs[0] must be bigger than zero
            // Find the zero corssings on the FilteredData before the P-wave phase arrival
            locs = Find(CheckRelation(Mult_El(ExtractSubset(Result.Data, 0, locs[0]-1), ExtractSubset(Result.Data, 1, locs[0])), "<", 0));  

            if (locs[0] != undefined) {
                // Find the last zero crossing and calculate the time of the PWave onset
                PWaveArrival = locs.slice(-1)[0] * delt;

                if (PWaveArrival == 0.0) { PWaveArrival = undefined; }
            }
        }
    }

    return {
        PWaveArrival: PWaveArrival,
        FilterOrder : Result.FilterOrder,
        F1          : Result.F1,
        F2          : Result.F2,
        Period      : Tn,
        Damping     : ksi,
        numBins     : numBins,
        Mode        : Mode,
    };
   
    function Normalize_Filter(Data, Type) {
         // Decleration of variables 
        let a, b, zf, FilterFlag, Ind, FilterOrder, F1, F2, Result=[];

        // Remove mean
        Data = Detrend(Data, 0);

        // Normalize Data to prevent numerical instability from very low amplitudes
        Data = Divide(Data, Max(Abs(Data)).val);

        // Select FilterOrder, F1, and F2 based on Type of Record (e.g, Acc and Vel)
        if (Type == 0) {
            // Acceleration records - strongMotion
            FilterOrder = 4;
            F1          = 0.1;
            F2          = 20;
            FilterFlag  = true;
        }
        else if (Type == 1) {
            // Velocity Records - weakMotion
            FilterOrder = 4;
            F1          = 0.1;
            F2          = 10;
            FilterFlag  = true;
        }
        else {
            FilterFlag = false;
        }

        if (FilterFlag) {
            // Desing 4th order Butterworth BandPass filer 
            [b, a, zf] = Butterworth_BandPass(FilterOrder, F1, F2, 1/delt);

            // Filter the normalized Data
            [Data, zf] = FiltFilt(b, a, Data);
        }
        
        // Apply Baseline Correction
        Data = Detrend(Data, 1);
        
        // Find the index of absoulute peak value 
        Ind = Find(CheckRelation(Abs(Data), "==",  Max(Abs(Data)).val))[0]; 

        // Take segment of waveform from beginning to absolute peak value (recommended for fast processing)
        Data = ExtractSubset(Data, 0, Ind);

        return {
            Data        : Data,
            FilterOrder : FilterOrder,
            F1          : F1,
            F2          : F2,
        };
    }

    function StateLevel(Ed, N) {
            
        let i, yMax, yMin, idx, H, Ry, dy, bins;
        let iLow, iHigh, lLow, lHigh, Mode;

        yMax = Max(Ed).val;
        yMin = Min(Ed).val - Number.EPSILON;

        // Compute Histogram
        idx = Mult(Subtract(Ed, yMin), (N/(yMax-yMin))).map((v,i) => Math.floor(v));
        idx.map((v,i) => {if ((v<0)||(v>N-1)) {idx.splice(i,1);}});
        H = Zeros(N);
        for (i=0; i<idx.length; i++) { H[idx[i]] = H[idx[i]] + 1;}

        // Compute Center of Each Bin
        yMin = Min(Ed).val;
        Ry = yMax - yMin;
        dy = Ry / N;
        bins = Add(Mult(Subtract(LinSpace(1, N, N), 0.50), dy), yMin)

        // Compute State Levels
        for (i=0; i<H.length; i++) { if (H[i]>0) {iLow = i; break;} }
        for (i=H.length-1; i>=0; i--) { if (H[i]>0) {iHigh = i; break;} }

        lLow  = Copy(iLow);
        lHigh = iLow + Math.floor((iHigh - iLow)/2);

        // Calculate the mode of the hiest bin in lower-histogram
        Mode = HistMode(H,  lLow, lHigh, dy);

        return Mode;
    }

    function HistMode(H, lLow, lHigh, dy) {

        let Ind, X1=[], X2=[], Y1=[], Y2=[], a1, a2, b1, b2, R;

        // Index of the hiest bin of lower histograms
        Ind = lLow + Max(ExtractSubset(H, lLow, lHigh)).Indx;

        if (Ind == 0) {
            X1 = [0,    1    ];
            Y1 = [0,    H[0] ];
            X2 = [0,    1    ];
            Y2 = [H[0], H[1] ];
        }
        else if (Ind == lHigh) {
            X1 = [0,        1      ];
            Y1 = [H[Ind],   0      ];
            X2 = [0,        1      ];
            Y2 = [H[Ind-1], H[Ind] ];
        }
        else {
            X1 = [0,        1        ];
            Y1 = [H[Ind],   H[Ind+1] ];
            X2 = [0,        1        ];
            Y2 = [H[Ind-1], H[Ind]   ];
        }
        
        R = PolyFit(X1, Y1, 1);
        a1 = R.P[0] / R.Std;     
        b1 = R.P[1] - (R.P[0] / R.Std * R.Mean);

        R = PolyFit(X2, Y2, 1);
        a2 = R.P[0] / R.Std;     
        b2 = R.P[1] - (R.P[0] / R.Std * R.Mean);

        return (Ind + (b2-b1)/(a1-a2)) * dy;
    }
}

//----------------------------------------------------------------------------------------------------------------------
// PWave detection
function Eq_Interpolation(X, M) {
    // Interpolating a time series is an important task in ground-motion processing and in computing response spectra.
    // There are many time-domain resampling methods (for example, curve fitting) in the literature.
    // An acceleration time series sampled at N samples-per-second (sps) should have no energy beyond the Nyquist frequency (fNyquist).
    // Ideally, any interpolation method should not introduce energy at frequencies beyond fNyquist of the original time series.
    // Linear interpolating to a higher sampling rate in time domain with an anti-aliasing filter violates this condition, 
    // and such resampled records are expected to lead to errors in the Fourier spectra and response spectra.
    // Other approaches for interpolation include zero-padding in the frequency-domain, windowed sinc interpolation in the time-domain
    // and resampling using poly-phase filter implementation.
    //
    // This function implements the frequency-domain zero padding (FDZP) approach for resampling from N samples to (N×M) samples. 
    // The procedure consists of steps defined in Lyons, 2014
    // Lyons, R.G., 2014, Understanding digital signal processing, third edition: Indiana, Prentice Hall, 954 p.

    let N, a, K;
    let Re, Im, Ind, Re1=[], Im1=[];

    // Length of data array
    N = X.length;

    // Round the Scaling factor, M
    M = Math.round(Abs(M));

    // Calculate FFT of data array
    [Re, Im] = FFT(X);

    // New Array of Interpolated Data 
    K = N * M;
    Re1 = new Array(K).fill(0);
    Im1 = new Array(K).fill(0);

    // The last index of positive frequencies 
    if (IsEven(N)) { a=(N/2); } else { a=(N-1)/2+1; }

    // Zero-pad FFT frequency in the middle to preserve the conjugate symmetry in the transform domain
    for (i=0; i<a; i++) {
        Re1[i] = Re[i];
        Im1[i] = Im[i];
    }
    Ind = K - a;
    for (i=a; i<N; i++) {
        Re1[Ind] = Re[i];
        Im1[Ind] = Im[i];
        Ind++;
    }

    // Compute IFFT
    [Re1, Im1] = IFFT(Re1, Im1);
    
    // Scale the Real part of IFFT by M to compansate for the 1/M amplitude loss 
    return Mult(Re1, M);
}

function ZeroCrossing(Data) {
    // List of indexes of zero crossings in Data[] array

    // Decleration of variables
    let i, Zp, Zc, result=[];

    // Starting Zone
    Zp = 0;
    for (i=0; i<Data.length; i++) { if (Sign(Data[i]) != 0) { Zp = Sign(Data[i]); break; } }
    
    // Determine all zero crossings
    if (Zp != 0) { 
        for (i=0; i<Data.length; i++) { 
            Zc = Sign(Data[i]);
            if ((Zc == Zp) || (Zc == 0)) { continue; } else { result.push(i-1); Zp=Zc; }
        }
    }
    return {
        All: result,
        First: result[0],
        Last: result.slice(-1)[0],
    }
}

function Cosine_Taper(Data, Nt, Ne) {
    
    // Nt : The width of the full cosine taper, Ntaper (in samples), is taken as the number of samples from the
    // beginning of the record to the last zero crossing before the event onset.
    //
    // Ne : the window length of the cosine taper at the end of the Data signal.  (e.g., 3 seconds)
    // The taper length typically is one half of the duration of the pre-event for the begining of the Data
    
    // Nt must be an even number 
    if (!IsEven(Nt)) { Nt = Nt + 1; }

    // length of Data
    N = Data.length;
    
    // Create cosine taper 
    w = Ones(N);
    for (i=0;      i<Nt/2; i++) { w[i] = 0.5 * (1-Math.cos(2*Math.PI*i/Nt));               }
    for (i=N-Nt/2; i<N;    i++) { w[i] = 0.5 * (1+Math.cos(2*Math.PI/Nt*(i-(N-Nt/2)+1)));  }
    return Mult_El(Data, w);
}
//----------------------------------------------------------------------------------------------------------------------
// Numerical Integration and Differentiation

function Cumtrapz(Y, delt) {
    // Computes the approximate cumulative integral of Y[] via the trapezoidal method with unit spacing of delt
    // This function is cross checked against Matlab Q = cumtrapz(delt, Y) built-in function

    // Check default values of input arguments
    if (delt == null) { delt = 1; };

    // Check if the input arguments are correct; otherwise, throw an error.
    //if (delt <= 0)         {throw new Error("delt (time interval) cannot be equal or less than zero.");}
    //if (!Array.isArray(Y)) {throw new Error("Y[] must be an array.");};

    let n = Y.length,  sum = 0;
    delt /= 2;
    return new Array(n).fill(0).map((v,i) => {if (i>0){sum+=(Y[i]+Y[i-1])*delt; return sum;} else {return 0;} });
}

function Trapz(Y, delt) {

    if (Array.isArray(delt) && typeof(delt[0].length) == "undefined") {
        // delt[] is an 1D array
        // Computes the approximate integral of Y[] array via the trapezoidal method with respect to the coordinates
        // specified by delt[] array.
        // Y[] and delt[] arrays must have the same length

        // Check if the input arguments are correct; otherwise, throw an error.
        if (!Array.isArray(Y))       {throw new Error("Y[] must be an array.");};
        if (!Array.isArray(delt))    {throw new Error("delt[] must be an array.");};
        if (Y.length != delt.length) {throw new Error("Dimensions of Y[] and delt[] arrays are inconsistent");};

        var result = 0;

        for (i = 1; i < Y.length; i++) {
            result += ((Y[i - 1] + Y[i]) / 2) * (delt[i] - delt[i - 1]);
        }

        return result;

    } else {
        // delt is a scalar
        // Computes the approximate integral of Y via the trapezoidal method with unit spacing of delt

        // Check default values of input arguments
        if (delt == null) { delt = 1; };

        // Check if the input arguments are correct; otherwise, throw an error.
        if (delt <= 0)            {throw new Error("delt (time interval) cannot be equal or less than zero.");}
        if (!Array.isArray(Y))    {throw new Error("Y[] must be an array.");};

        var n = Y.length;
        return (Statistics_Sum(Y) - (Y[0] + Y[n-1]) / 2) * delt;
    }

}

function Compute_VelDisp(Acc, FSamp, PWaveArrival) {
    // Check the sampling rate of the acceleration time series.
    // If the acceleration has low time resolution (less than about 200 samples-per-second [sps]), then the time series is
    // interpolated to 200 sps by resampling in the frequency domain prior to processing in order to reduce numerical noise (Kalkan and Stephens, 2017).

    // Decleration of variables
    let M, RawData, delt, Type, ksi, Tn, numBins, PWaveArrival_Ind
    let Vel, Vel_Trend, th
    
    // STEP-0 
    // Interpolate the Data if Sampling rate is less than 200 Hz
    if (FSamp < 200) {
        // Scaling Factor
        M = Math.ceil(200 / FSamp);
        RawData = Eq_Interpolation(Acc, M);
        FSamp = M * FSamp;
    }
    else {
        RawData = Copy(Acc);
    }

    // Sampling Interval
    delt = 1 / FSamp;

    // STEP-1  (Determine Event Onset) --------------------------------------------------------------------------------
    if (PWaveArrival == null) {
        // Determine Event Onset 
        Type    = 0;
        ksi     = 0.60;
        Tn      = 100;
        numBins = 200;
        PWaveArrival = PhasePicker_Hist(Acc, delt, Type, ksi, Tn, numBins).PWaveArrival;

        if (PWaveArrival == undefined) { 
            message = "PWave Onset is not detected! Calculation is Terminated!";
            return;
        }
    }

    // Index of PWave onset
    PWaveArrival_Ind = Math.floor(PWaveArrival / delt);

    // Apply Buffer Interval of 1 second to pre-event onset if the length of pre-event greater than 3 seconds.
    // Buffer Interval reduces the uncertainty in the onset-time. 
    if (PWaveArrival_Ind >= 3*FSamp) { PWaveArrival_Ind -= 1 * FSamp;}

    // Pre-Event interval
    Pre_Event       = ExtractSubset(RawData, 0, PWaveArrival_Ind);
    Pre_Event_Mean  = Mean(Pre_Event);
    Pre_Event_Trend = BestFitTrend(Pre_Event, true); 

    // STEP-2  (Remove pre-event mean from RawData) -------------------------------------------------------------------
    RawData = Subtract(RawData, Pre_Event_Mean);

    // STEP-3  (Integrate to Velocity) --------------------------------------------------------------------------------
    Vel = Cumtrapz(RawData, delt);

    // STEP-4  (Compute the best Fit-Trend in Velocity) ---------------------------------------------------------------
    Vel_Trend = BestFitTrend(Vel, true);

    // STEP-5  (Remove derivative of bestFit-Trend from Accelertation)
    if (Vel_Trend.length == 2) {
        // Best fit is a 1st order line; therefore, remove the mean of the trend from RawData[] array
        RawData = Subtract(RawData, Vel_Trend[1]);
    }
    else {
        // Best fit is a 2nd order polynomial; therefore, remove the derivative of the trend from RawData[] array
        RawData = RawData.map((v,i) => v - Vel_Trend[1] - 2*Vel_Trend[2]*(i+1));
    }
    
    // Integrate Acceleration to Velocity
    Vel = Cumtrapz(RawData, delt);

    // STEP-6  (Quality Check on Velocity)
    th = 0.01;
    if (QualityControl(Vel, th)) {
        
        // The width of cosine taper (Nt) is taken as the number of samples from beginning of the record
        // to the last zero crossing (a2) before the event onset. 
        a2 = ZeroCrossing(Pre_Event).Last;
        
        // Compute the width of cosine taper at each end of Velocity record as the minimum of following:
        //   - Half length of pre-event interval from begining to the last zero crossing 
        //   - 2 seconds
        Nt = Min(Math.floor(a2/2), 2*FSamp).val;

        // Cosine Taper of Velocity 
        Ne = 3 * FSamp;
        RawData = Cosine_Taper(RawData, Nt, Ne);

        // Zero-padding of Velocity at either end 
        // Vector of Zero padding (Number of samples) on each of the Velocity
        Tp  = Zeros(Math.ceil(1.5 * 4 / 0.1 * FSamp / 2));
        RawData = Concat(Flip(Concat(Flip(RawData),Tp)), Tp);

        // Acasual BandPass Filter
        N = 4;
        fc1 = 0.1;
        fc2 = (FSamp / 2) * 0.80;
        [b, a, zf] = Butterworth_BandPass(N, fc1, fc2, FSamp);
        [RawData, zf] = FiltFilt(b, a, RawData);

        // Integrate to Velocity and Displacement 
        Vel = Cumtrapz(RawData, delt);
        Disp = Cumtrapz(Vel, delt);

        console.log("QC OK ***************************")

        return [RawData, Vel, Disp];

    }
    else {
        // Run Adaptive Baseline Correction 
        console.log("QC NOT OK")
    }



    // Quality Control
    function QualityControl(Vel, th) {

        N = Vel.length;
        WinLen = 100;
        w1 = WinLen - 1;
        w2 = N - WinLen; 

        // Find first and last zero crossings for leading interval
        a1 = ZeroCrossing(ExtractSubset(Vel, 0, w1)).First;
        if (a1 == undefined) { return false; }

        // Find first and last zero crossings for trailing interval
        b2 = ZeroCrossing(ExtractSubset(Vel, w2, N-1)).Last;
        if (b2 == undefined) { return false; }

        // Mean of leading and trailing intervals 
        M1 = Mean(ExtractSubset(Vel, 0, a1));
        M2 = Mean(ExtractSubset(Vel, b2, N-1));

        if ((M1 <= th) && (M2 <=th)) { return true; } else { return false; }

    }

    // Best-FitTrend
    function BestFitTrend(X, Opt) {

        // Choose best regression based on Root Mean Square Deviation (RMSD)
        // Decleration of variables
        let R1, R2, R3;

        // Opt == true  ===>> Use 1st and 2nd order polynomials 
        // Opt == false ===>> Use 1st, 2nd and 3rd order polynomilas 
        if (Opt == null) { Opt = false; }

        if (Opt) {
            // Calculate regressions of 1st and 2nd order polynomials  
            R1 = Detrend(X, 1, true);
            R2 = Detrend(X, 2, true);

            if (R1.RMDS < R2.RMDS) { return R1.Coeff; } else { return R1.Coeff; }
        }
        else {
            // Calculate regressions of 1st, 2nd and 3rd order polynomials 
            R1 = Detrend(X, 1, true);
            R2 = Detrend(X, 2, true);
            R2 = Detrend(X, 3, true);

            if      ((R1.RMDS < R2.RMDS) && (R1.RMDS < R3.RMDS)) { return R1.Coeff; }
            else if ((R2.RMDS < R1.RMDS) && (R2.RMDS < R3.RMDS)) { return R2.Coeff; }
            else if ((R3.RMDS < R1.RMDS) && (R3.RMDS < R2.RMDS)) { return R3.Coeff; }
        }
    }

}

//----------------------------------------------------------------------------------------------------------------------
// Ground Motion Analysis including Ambient Vibration

function Arias(data, delt, g) {
    // Returns the Arias Intensity defined by the Arias in 1970
    // Unit of the Arias intensity is m/s  provided that the data[] array has the units of m/s**2
    // Earth's gravity is taken as g = 9.80665 m/s**2  if not given

    //  data[]:  Data vector
    //    delt:  Sampling interval
    //       g:  Earth's gravity
    //
    // Outputs:
    //         AI:  Cumulative Arias Intensity
    //  AI_MaxVal:  Maximum value of Arias Intensity
    //         T1:  Time  5% of AI_MaxVal is exceeded
    //         T2:  Time 95% of AI_MaxVal is exceeded
    //         Td:  Significant Duration

    // Check default values of input arguments
    if (delt == null) { delt = 1;       };
    if (   g == null) {    g = 9.80665; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if (delt <= 0)            {throw new Error("delt (time interval) cannot be equal or less than zero.");}
    if (!Array.isArray(data)) {throw new Error("data[] must be an array.");};

    // Calculate Arias Intensity
    AI        = Vector_Multiply(Cumtrapz(Vector_Pow(data, 2), delt), (Math.PI/2/g));
    AI_MaxVal = Statistics_Max(AI).maxVal;

    A1 = 0.05 * AI_MaxVal;
    A2 = 0.95 * AI_MaxVal;

    flag1 = true;
    flag2 = false;
    T1 = 0;
    T2 = 0;

    AI.map( (v, i) => {
        if ((v >= A1) && (flag1)) { T1 = i * delt; T2 = T1; flag1 = false; flag2 = true;  }
        if ((v >= A2) && (flag2)) { T2 = i * delt;                         flag2 = false; }
    });

    Td = T2 - T1;

    return [AI, AI_MaxVal, T1, T2, Td];
}

function Cav(Y, delt) {
    // Returns the Cumulative absolute velocity
    // Y[] array contains the velocity readings with unit spacing of delt.

    // Check default values of input arguments
    if (delt == null) { delt = 1; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if (delt <= 0)         {throw new Error("delt (time interval) cannot be equal or less than zero.");}
    if (!Array.isArray(Y)) {throw new Error("Y[] must be an array.");};

    return Cumtrapz(Vector_Abs(Y), delt)

}

function HVSR(EW, NS, UD, RWL, OVS, FSamp, NFFT, Option) {
    // This function returns the average H/V spectral ratio using 1D arrays of EW[], NS[], UD[].
    // The each of these arrays are divided into multiple segments of RWL-samples and overlapped by OVS sample.
    // Each segment is multiplied by Hamming window to prevent spectral leakage before taking the NFFT-point of FFT.
    // No smoothing is performed on the calculated FAS.
    // The FFT of the two horizontal components (EW[] and NS[] arrays) are combined before H/V calculation using one
    // of the 4 Options as follow:
    //
    //  Inputs:
    //      EW      : East-West recording - 1D array
    //      NS      : North-South recording - 1D array
    //      UD      : Vertical recording - 1D array
    //      RWL     : Number of samples in each data-segment
    //      OVS     : Number of samples that two data segments are overlapped
    //      NFFT    : Number of discrete Fourier Transform points in FFT estimate for each data segment
    //      FSamp   : Sampling frequency
    //      Option  : Determines how to combine two horizontal components of EW[] and NS[] arrays
    //                1: Geometric Mean
    //                2: Vector summation
    //                3: Quadratic Mean
    //                4: Arithmetic mean
    //
    //  Outputs:
    //      HV      : H/V Ratio Spectrum
    //      Std     : Standard deviation for each discrete frequency of the H/V ratio spectrum
    //      f       : Frequency vector over which the H/V ratio spectrum is calculated
    //
    //  by Dr. Yavuz Kaya, P.Eng.;
    //  Created on      17.Mar.2017
    //  Last modified   01.Jly.2023
    //  Yavuz.Kaya.ca@gmail.com

    // Check if the input arguments are correct; otherwise, throw an error.
    if (!Array.isArray(EW)) {throw new Error("EW[] must be an array.");};
    if (!Array.isArray(NS)) {throw new Error("NS[] must be an array.");};
    if (!Array.isArray(UD)) {throw new Error("UD[] must be an array.");};

    // Check length of EW[], NS[], and UD[] arrays
    if ((EW.length != NS.length) || (EW.length != UD.length)) {throw new Error("Dimensions of EW[], NS[] and UD[] arrays are inconsistent.");};

    // Check default values of input arguments
    if (RWL    == null) { RWL    = Math.max(2, Math.floor(EW.length / 8.00)); };
    if (OVS    == null) { OVS    = Math.floor(RWL * 0.25); };
    if (FSamp  == null) { FSamp  = 1; };
    if (NFFT   == null) { NFFT   = NextPow2(RWL); };
    if (Option == null) { Option = 3; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if ((RWL <= 0) || (RWL > EW.length) || (RWL > NFFT) || (OVS >= RWL) || (OVS < 0) || (NFFT < 2) || (FSamp <= 0)) {
        {throw new Error("Input arguments are inconsistent (RWL <= 0) || (RWL > EW.Length) || (RWL > NFFT) || (OVS >= RWL) || (OVS < 0) || (NFFT < 2) || (FSamp <= 0)");};
    }

    // Declare variables
    var HV1, HV2, Std, DW, Win, K, a1, a2, EW1, NS1, UD1, FS1, FS2, FS3, Mag1, Mag2, Mag3, df, f

    // Pre-allocation
    HV1 = new Array(NFFT).fill(0);
    HV2 = new Array(NFFT).fill(0);
    Std = new Array(NFFT).fill(0);

    // Calculate initial parameters
    DW = RWL - OVS;

    // Create a window function of RWL-point and normalize it so that the power is 1 Watt.
    Win = Hamming(RWL);
    Win = Vector_Multiply(Win, Math.sqrt(RWL / Statistics_Sum(Vector_Pow(Win, 2)) ) );

    // Update Indexes
    K = 0;
    a1 = 0;
    a2 = RWL - 1;

    while (a2 < EW.length) {

        EW1 = Vector_GetRange(EW, a1, a2);
        NS1 = Vector_GetRange(NS, a1, a2);
        UD1 = Vector_GetRange(UD, a1, a2);

        // Skip the this segment if NaN exists in any of the EW1[], NS1[], or UD1[] arrays
        if (IsNaN(EW1) || IsNaN(NS1) || IsNaN(UD1)) { a1 += DW; a2 += DW; continue; }

        // Extract the segment of the EW[], NS[] and UD[] arrays starting from Index a1 to Index a2
        // Then calculate the Fourier Amplitude Spectrum (FAS) of each windowed-segment
        [Re1, Im1] = FFT(Vector_Multiply(EW1, Win), null, NFFT);
        [Re2, Im2] = FFT(Vector_Multiply(NS1, Win), null, NFFT);
        [Re3, Im3] = FFT(Vector_Multiply(UD1, Win), null, NFFT);

        // Calculate Magnitude Spectrum
        [Mag1, ] = FFT_MagPhase(Re1, Im1);
        [Mag2, ] = FFT_MagPhase(Re2, Im2);
        [Mag3, ] = FFT_MagPhase(Re3, Im3);

        // Fourier Amplitude Spectrum need to be smoothed.
        // However, Smoothing is not applied because the Konno & Ohmachi smoothing
        // is very time consuming to calculate; therefore, it is not performed.
        // Create Konno & Ohmachi window for smoothing.

        // Combine two horizontal recordings
        if      (Option == 1) {
            // Geometric Mean is recommended by the SESAME project (2004) and also
            // adopted by Picozzi et al. (2005), Haghshenas et al. (2008), Pileggi et al. (2011)
            Mag1 = new Array(NFFT).fill().map((v,i) => Math.sqrt(Mag1[i] * Mag2[i]));
        }
        else if (Option == 2) {
            // Vector summation used by Sauriau et al. (2007) and Puglia et al. (2011),
            Mag1 = new Array(NFFT).fill().map((v,i) => Math.sqrt(Mag1[i]*Mag1[i] + Mag2[i]*Mag2[i]));
        }
        else if (Option == 3) {
            // Quadratic mean is considered by Bonnefoy-Claudet et al. (2006, 2008) and Fäh et al. (2001)
            Mag1 = new Array(NFFT).fill().map((v,i) => Math.sqrt((Mag1[i]*Mag1[i] + Mag2[i]*Mag2[i])/2));
        }
        else if (Option == 4) {
            // Arithmetic mean is considered by Chavez-Garcia et al. (2007)
            Mag1 = new Array(NFFT).fill().map((v,i) => (Mag1[i] + Mag2[i])/2);
        }

        // Calculate H/V ratios
        Mag3 = new Array(NFFT).fill().map((v,i) => Mag1[i] / Mag3[i] );         //  Calculate H/V

        // Sum all H/V ratios
        HV1 = new Array(NFFT).fill().map((v,i) => HV1[i] + Mag3[i] );          //  Sum of H/V ratios
        HV2 = new Array(NFFT).fill().map((v,i) => HV2[i] + Mag3[i]*Mag3[i] );  //  Sum of (H/V)^2

        // Update the Indexes
        a1 += DW;
        a2 += DW;
        K++;
    }

    if (K > 1) {

        // Average the individual HVSR
        HV1 = Vector_Divide(HV1, K);

        // Normalize HV2 by K
        HV2 = Vector_Divide(HV2, (K*K));

        // Calculate standard deviation of the averaged-HVSR for each discrete frequency
        Std = new Array(NFFT).fill().map((v,i) => Math.sqrt( (HV2[i]-HV1[i]*HV1[i]/K) / (K-1) )  );

    }

    // Frequency vector
    df = FSamp / NFFT;
    f  = Vector_LinSpace(0, (NFFT-1)*df, NFFT);

    // Return OnceSided Spectrum Only
    if (IsEven(NFFT)) { NumSamp   = (NFFT / 2) + 1; } else { NumSamp = (NFFT + 1) / 2; }

    HV1 = Vector_GetRange(HV1, 0, NumSamp-1);
    Std = Vector_GetRange(Std, 0, NumSamp-1);
    f   = Vector_GetRange(f,   0, NumSamp-1);

    return [HV1, Std, f];
}

function RotateRecord(x, y, theta) {
    // This function rotates records of x[] and y[] array in the x-y plane counterclockwise through an angle of theta
    // with respect to the x axis about the origin of a two-dimensional Cartesian coordinate system.
    //
    //    R is the Rotation matrix
    //    R = |  cos(theta)  -sin(theta) |   *   | x |    =    | x*cos(theta)-y*sin(theta) |
    //        |  sin(theta)   cos(theta) |       | y |         | x*sin(theta)+y*cos(theta) |
    //

    // Check if the input arguments are correct; otherwise, throw an error.
    if (!Array.isArray(x)) {throw new Error("x[] must be an array.");};
    if (!Array.isArray(y)) {throw new Error("y[] must be an array.");};

    var c, s, Ar, xh, yh;

    Ar = Deg2Rad(theta);
    c  = Math.cos(Ar);
    s  = Math.sin(Ar);

    xh = new Array(n).fill().map((v,i) => x[i]*c - y[i]*s );
    yh = new Array(n).fill().map((v,i) => x[i]*s + y[i]*c );

    return [xh, yh];
}

function SpecIntensity(data, delt, ksi, T) {
    // Returns the Spectral Intensity defined by the Housner in 1952
    // Total area underneath the Spectral Velocity between 0.1s and 2.5s is calculated.
    // The total area is then dived by 2.4

    let SD, SV, SA, Sa;

    // Check default values of input arguments
    if (delt == null) { delt = 1;    };
    if (ksi  == null) { ksi  = 0.20; };
    if (T    == null) { T    = [0.10, 0.15, 0.20, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5]; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if (delt <= 0)               {throw new Error("delt (time interval) cannot be equal or less than zero.");};
    if ((ksi < 0) || (ksi >= 1)) {throw new Error("ksi (damping ratio) cannot be less than zero or greater than or equal to 1.");};
    if (!Array.isArray(data))    {throw new Error("data[] must be an array.");};
    if (!Array.isArray(T))       {throw new Error("T[] must be an array.");};

    // Calculate Response Spectrum
    [SD, SV, SA, Sa] = SDOF_ResSpectrum(data, delt, ksi, T);

    // Return Spectral Intensity as the area underneath the velocity spectrum
    return Trapz(SV, T);

}

function BracketedDuration(data, delt, th) {

    // Returns the bracketed duration of the data[] vector
    //  data[]:  Data vector
    //    delt:  Sampling interval
    //      th:  Threshold
    //
    // Outputs:
    //      T1:  Time of the first exceedance of threshold in the data[] vector
    //      T2:  Time of the last exceedance of threshold in the data[] vector
    //      Td:  Bracketed duration (T2 - T1) in seconds

    flag = true;
    T1 = 0;
    T2 = 0;

    data.map( (v, i) => {
        if      ((Math.abs(v) >= th) && (flag)) { T1 = i * delt; T2 = T1; flag = false; }
        else if  (Math.abs(v) >= th)            { T2 = i * delt;                        }
    });

    Td = T2 - T1;
    return [T1, T2, Td];
}

//----------------------------------------------------------------------------------------------------------------------
// System of Linear Equations

function LinearEquations_Hessenberg(data) {
    // Hessenberg form of matrix
    // H = Hessenberg(data) finds H, the Hessenberg form of matrix data[,].
    // [P,H] = Hessenberg(data) produces a Hessenberg matrix H and a unitary matrix P so that A = P*H*P'
    // data[,] matrix must be square.
    // This function is cross-checked against Matlab [P,H]=hess(A) built-in function 2023-06-12

    // A transformation that reduces a general matrix to Hessenberg form will reduce a Hermitian matrix to tridiagonal form.
    // Tridiagonal matrix is a band matrix that has nonzero elements only on the main diagonal, the subdiagonal/lower diagonal (the first diagonal below this),
    // and the supradiagonal/upper diagonal (the first diagonal above the main diagonal)

    var H, n, z, v, K, P, i, j, k;

    // Throws error if the matrix is not square...
    if (data.length != data[0].length) {throw new Error("Matrix must be square.");}

    // Make a copy of the data[,] matrix
    H = Matrix_Copy(data);
    n = H[0].length;
    P = Matrix_Eye(n);

    for (k = 0; k < n - 2; k++) {

        // Find the Householder reflections
        z = Vector_GetRange(Matrix_GetColumn(H, k), k + 1, n - 1);
        v = Vector_Multiply(Vector_Copy(z), -1);
        v[0] = -Math.sign(z[0]) * Vector_Norm(z) - z[0];
        v = Vector_Divide(v, Math.sqrt(Vector_Dot(v, v)));

        K = Matrix_Eye(n);

        for (i = 0; i < v.length; i++) {
            for (j = i; j < v.length; j++) {
                K[i + k + 1][j + k + 1] -= 2 * v[i] * v[j];
                K[j + k + 1][i + k + 1] = K[i + k + 1][j + k + 1];
            }
        }

        // Unitary matrix P
        P = Matrix_Multiply(P, K);

        // T[n,n] = V[n] * V[n]
        // K[n,n] = I[n,n] - 2 * T[n,n]
        // H = K * H * K
        H = Matrix_Multiply(Matrix_Multiply(K, H), K);
    }
    return [P, H];
}

function GaussElimination(A, y, Opt) {

    // Gauss Elimination to solve systems of linear equation of (A * x) = y
    // If Opt is false, then no row exchange is applied unless the leading term is equal to zero.
    // This function is double checked for correctness 2023-06-12

    // Check default values of input arguments
    if (Opt == null) { Opt = true; };

    // Check if the dimensions of input arguments are correct; otherwise, throw an error.
    //if ((A.length != A[0].length) || (A.length != y.length)) { throw new Error("Dimensions of A[,] and y[] are inconsistent."); }
    
    var IA, n, U, i, j, k, pivot, PP;
    
    n  = A.length;
    IA = infoArray(y);

    if (IA.is1D)  { 

        // Copy from A[,] array to U[,] array
        U = Copy(A);

        // Forward Substitution with optional row-exchange
        for (i = 0; i < n; i++) {
            // Row-exchange
            if ((U[i][i] == 0) || Opt) {
                // Determine the maximum value and the Index number of each column of the U[,] array
                PP = Max(Abs(ExtractSubset(GetColumn(U, i), i, n - 1))).Indx;
                U = SwapRow(U, [PP + i], [i]);
                y = SwapRow(y, [PP + i], [i]);
            }
            
            for (j = i + 1; j < n; j++) {
                pivot = U[j][i] / U[i][i];
                y[j] -= pivot * y[i];

                for (k = i; k < n; k++) {
                    U[j][k] -= pivot * U[i][k];
                }
            }
        }

        // Backward Substitution - No row exchange is required
        for (i = n - 1; i > 0; i--) {
            // if U[i, i] is equal to ZERO, Divide by Zero error occurs.
            // Then there is no solution.  This code needs to be improved to avoid division by zero
            y[i] /= U[i][i];
            U[i][i] = 1;

            for (j = 0; j < i; j++) {
                y[j] -= U[j][i] * y[i];
            }
        }

        y[0] /= U[0][0];

        return y;
    }
    else if (IA.is2D) { 
        for (i=0; i<y[0].length; i++) { ReplaceCol(y, i, GaussElimination(A, GetColumn(y, i) , Opt)); }
        return y;
    }
}

//let A = [[2, 4], [ 8,  5]]; 
//let y = [[7,14], [-5,-10]];
//let temp1 = GaussElimination(A, y);
//console.log(temp1)

//----------------------------------------------------------------------------------------------------------------------
// Polynomials

function Conv(a, b) {
    // This function returns the convolution of a[] and b[] arrays of real-values.
    // Resulting array has the length of (a.length + b.length - 1)
    // If a[] and b[] are arrays of polynomial coefficients, convolving them is equivalent to multiplying the two polynomials.
    // This function cross checked against Matlab w=conv(u,v) built-in function 2023-06-19

    // Check if the input arguments are correct; otherwise, throw an error.
    if (!Array.isArray(a)) {throw new Error("a[] must be an array.");};
    if (!Array.isArray(b)) {throw new Error("b[] must be an array.");};

    //
    let N1, N2, N, N2FFT, Re, Im;

    // Calculate convolution of a[] and b[] arrays
    N1 = a.length;
    N2 = b.length;
    N = N1 + N2 - 1;

    N2FFT = NextPow2(N);

    // pre-allocate variables
    let Re1 = new Array(N2FFT).fill(0);
    let Im1 = new Array(N2FFT).fill(0);
    let Re2 = new Array(N2FFT).fill(0);
    let Im2 = new Array(N2FFT).fill(0);

    [Re1, Im1] = FFT(a, null, N2FFT);
    [Re2, Im2] = FFT(b, null, N2FFT);

    for (i = 0; i < N2FFT; i++) {

        Re = Re1[i] * Re2[i] - Im1[i] * Im2[i];
        Im = Re1[i] * Im2[i] + Im1[i] * Re2[i];

        Re1[i] = Re;
        Im1[i] = Im;
    };

    [Re1, Im1] = IFFT(Re1, Im1, N2FFT);

    return Vector_Truncate(Re1, N);
}

function Detrend_With_PolyFit(X, N) {

    // Removes a polynomial trend of order N from the vector X[].
    // If N = 0, removes the mean for the data in X
    //
    // Y = Detrend(X) removes the best straight-line fit linear trend from the data in X[] array and
    // returns the residual in r[] array.

    // Check default values of input arguments
    if (N == null) { N = 1; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if (N < 0) {throw new Error("N cannot be less than zero.");};

    if (Array.isArray(X) && typeof(X[0].length) == "number") {
        // X[][] is a matrix

        let nr, nc, time, data, i, r, result, PF;

        nr   = X.length;
        nc   = X[0].length;
        time = Vector_LinSpace(0, (nr-1), nr);
        result = new Array(nr).fill().map((v,i) => new Array(nc).fill(0));
        for (i =0; i < nc; i++) {
            data = Matrix_GetColumn(X, i);
            PF = PolyFit(time, data, N);
            Matrix_ReplaceColumn(result, i, PF.Residual);
        }
        return result;
    }
    else if (Array.isArray(X) && typeof(X[0].length) == "undefined") {
        // X[] is a 1D array

        let n, mu, time, r;

        if (N == 0){
            n = X.length;
            mu = Statistics_Mean(X);
            return new Array(n).fill().map((v,i) => X[i] - mu );
        }
        else {
            n    = X.length;
            time = Vector_LinSpace(0, (n-1), n);
            PF = PolyFit(time, X, N);
            return PF.Residual;
        }
    }
    else { return 0; }
}

function Detrend(X, N, Opt) {

    // Removes a polynomial trend of order N from the vector X[].
    // If N = 0, removes the mean for the data in X
    // Y = Detrend(X) removes the best fit curve of order N, in least-square sense, from the data in X[]

    // Decleration of variables 
    let n, mu, b1, b2, b3, b4, i, a1, a2, a3, a4, a5, a6, A=[], B=[], C=[], Result=[], rmsd;

    // Check default values of N argument
    if (N == null) { N = 1; };
    if (Opt == null) { Opt = false; }

    // Length of data.
    n = X.length;

    if (N == 0) {
        // Removes mean from X[] array
        mu = Mean(X);
        Result = new Array(n).fill().map((v,i) => X[i] - mu );
    }
    else if (N == 1) {
        // Remove 1st order straight line from X[] array
        b1 = Sum(X);
        b2 = 0;
        for (i=0; i<n; i++) { b2 += X[i] * (i+1); }

        a1 = n * (n + 1) / 2;
        a2 = n * (n + 1) * (2*n + 1) / 6;

        A = [[n, a1], [a1, a2]];
        B = [[b1], [b2]];
        C = Mult(Inv(A), B);  //  f(x)=C[1]*x + C[0]   ==>   C[1]: mean,   //C[0]:y-intercept

        // Detrend applied array 
        Result = new Array(n).fill().map((v, ii) => X[ii] - C[0] - C[1]*(ii+1));
    }
    else if (N == 2) {
        // Remove 2nd order parabola curve from X[] array
        b1 = Sum(X);
        b2 = 0;
        b3 = 0;
        for (i=0; i<n; i++) { b2 += X[i] * (i+1); b3 += X[i] * (i+1)**2; }

        a1 = n * (n + 1) / 2;
        a2 = n * (n + 1) * (2*n + 1) / 6;
        a3 = a1 * a1;
        a4 = (n * (n + 1) * (2*n + 1) * (3*n*n + 3*n - 1)) / 30;

        A = [[n, a1, a2], [a1, a2, a3], [a2, a3, a4]];
        B = [[b1], [b2], [b3]];
        C = Mult(Inv(A), B);

        // Detrend applied array 
        Result = new Array(n).fill().map((v,ii) => X[ii] - C[0] - C[1]*(ii+1) - C[2]*(ii+1)**2);
    }
    else if (N == 3) {
        // Remove 3rd order Cubic curve from X[] array
        b1 = Sum(X);
        b2 = 0;
        b3 = 0;
        b4 = 0;
        for (i=0; i<n; i++) { b2 += X[i] * (i+1); b3 += X[i] * (i+1)**2; b4 += X[i] * (i+1)**3; }

        a1 = n * (n + 1) / 2;
        a2 = n * (n + 1) * (2*n + 1) / 6;
        a3 = a1 * a1;
        a4 = (n * (n + 1) * (2*n + 1) * (3*n*n + 3*n - 1)) / 30;
        a5 = (n**2 * (n+1)**2 * (2*n*n + 2*n - 1)) / 12;
        a6 = (n * (n+1) * (2*n + 1) * (3*n**4 + 6*n**3 - 3*n +1)) / 42;

        A = [[n, a1, a2, a3], [a1, a2, a3, a4], [a2, a3, a4, a5], [a3, a4, a5, a6]];
        B = [[b1], [b2], [b3], [b4]];
        C = Mult(Inv(A), B);

        // Detrend applied array 
        Result = new Array(n).fill().map((v,ii) => X[ii] - C[0] - C[1]*(ii+1) - C[2]*(ii+1)**2 - C[3]*(ii+1)**3);
    }

    // Return results
    if (Opt) {
        // Return additional information about the regression 
        rmsd = 0;
        for (i=0; i<n; i++) { rmsd += Result[i]**2; }
        rmsd = Math.sqrt(rmsd / n);
        return {
            "Coeff"     :   C,
            "RMSD"      :   rmsd,
            "Result"    :   Result,
        };
    }
    else {
        // Return the result only
        return Result;
    }
}

function Poly(x) {
    // Returns the coefficients of the polynomial whose roots are the elements of x[] array
    // Converts complex-poles or real-poles in the x[] array to polynomial coefficients of a[] array
    
    // If x is a number, covert it to 1D array
    if (isNumber(x)) { x = [x]; }

    let i, j, n, c, c1, c2, c3;
    n = x.length; 
    c = new Array(n+1).fill().map((x) => new ComplexNum(0, 0));
    c[0].Re = 1;
    c[0].Im = 0;
    for (i = 0; i < n; i++) {
        c1 = ExtractSubset(c, 1, i + 1);
        c2 = ExtractSubset(c, 0, i);
        c2 = c2.map((v,j) => Mult(c2[j], x[i]) );
        c3 = new Array(i+1).fill().map((v,j) => c1[j].Subtract(c2[j]));
        for (j = 0; j < c3.length; j++) {
            c[j + 1] = c3[j];
        }
    }
    return c;
}

function PolyFit(X, Y, N) {

    // This function return coefficients for a polynomial p(X) of degree N that is a best fit (in a least-squares sense)
    // for the data in Y[] array. The length of p(X) is N+1 as follow:
    //
    //  p(X) = p(1) * x^N  +  p(2) * x^(N-1)  +   ...   +  P(N) * X^1  + P(N+1)
    //
    //  X    : Array of x-coordinates of data in Y[] array
    //  Y    : Array of data
    //  N    : Degree of polynomial
    //
    //  p    : Polynomial coefficients in descending order for each (X-mu)/std
    //  r    : Residual vector
    //  mu   : Mean value of data in X[] array
    //  stdv : Standard deviation of data in X[] array
    //  R    :
    //  df   : Degree of freedom
    //  Norm : Norm of residuals
    //
    //  by Dr. Yavuz Kaya, P.Eng.;
    //  Created on      03.Mar.2019
    //  Last modified   24.Jun.2023
    //  Yavuz.Kaya.ca@gmail.com

    var mu, stdv, V=[], p=[], r=[], df, Res;

    // Calculate the mean and the standard deviation of X[] array
    mu   = Mean(X);
    stdv = Std(X);

    // Remove mean and then divide by standard deviation
    X = X.map((v,i) => (X[i]-mu)/stdv);

    // Construct the Vandermonde matrix as follow:
    //
    //  |  X[1]^N   X[2]^(N-1)   ...   X[1]^2   X[1]   1   |
    //  |  X[2]^N   X[2]^(N-1)   ...   X[2]^2   X[2]   1   |
    //  |  X[3]^N   X[3]^(N-1)   ...   X[3]^2   X[3]   1   |
    //  |    .          .        ...      .      .     .   |
    //  |    .          .        ...      .      .     .   |
    //  |    .          .        ...      .      .     .   |
    //  |  X[.]^N   X[.]^(N-1)   ...   X[.]^2   X[.]   1   |
    //

    V = new Array(X.length).fill().map((vv,j) => new Array(N+1).fill().map((v,i) => Math.pow(X[j], N-i)));

    // QR Decomposition - return economy-size of Q[,] and R[,] matrices
    Res = QR(V, 1); 

    // Gauss-elimination to solve for the polynomial coefficients
    p = GaussElimination(Res.R, Mult(Transpose(Res.Q), Y));

    //  Residual
    r = Subtract(Y, Mult(V, p));

    // df:  Degree of Freedom
    df = Max([0, Y.length - (N+1)]).val;

    return {
        P           : p,
        Residual    : r,
        Mean        : mu,
        Std         : stdv,
        R           : Res.R,
        df          : df,
        Norm_res    : Norm(r)
    }
}

function Smoothing(data, Win) {
    // Smooths the data[] array using the window function in Win[] array


    // Check default values of input arguments
    if (Win == null) { Win = Rectwin(3); };

    // Check if the input arguments are correct; otherwise, throw an error.
    if (Win.length >= data.length)  {throw new Error("Length of smoothing window cannot be greater than number of samples in data[] array");};
    if (!Array.isArray(data))       {throw new Error("data[] must be an array.");};
    if (!Array.isArray(Win))        {throw new Error("Win[] must be an array.");};

    // Make sure that the length of smoothing window is an odd number
    // Throw error otherwise.
    if (Win.length %2 == 0) {throw new Error("Length of smoothing window must be an odd number.");}

    m = data.length;
    n = Win.length;

    // Normalize the Window function
    Factor = Math.sqrt(n / Win.reduce((acc, val, i) => acc + Win[i] * Win[i], 0));
    Win    = Vector_Divide(Win, Factor);

    // Remove the mean from the data[] array
    mu   = Statistics_Mean(data);
    data = Vector_Subtract(data, mu);

    // Calculate the scale factor
    SF     = (1/n) * (n/Statistics_Sum(Win));
    lc     = (n - 1) / 2;
    return Vector_GetRange(Conv(data, Win), lc, m+lc-1).map((v,i) => v * SF + mu);
}

function Xcorr(a, b, ScaleOpt) {
    // This function returns the cross-correlation of two discrete-time sequences of a[] and b[] arrays.
    // Cross-correlation measures the similarity between an array of a[] and shifted (lagged) copies of an array b[]
    // as a function of the lag.
    //
    // ScaleOpt variable specifies a normalization option for the cross-correlation or auto-correlation.
    // Any option other than 'none' (the default) requires a[] and b[] to have the same length.
    //
    // ScaleOpt = 'none'        : Unscaled cross-correlation
    // ScaleOpt = 'biased'      : Biased estimate of the cross-correlation
    // ScaleOpt = 'unbiased'    : Unbiased estimate of the cross-correlation
    // ScaleOpt = 'normalized'  : Normalizes the sequence so that the auto-correlations at zero lag equal 1
    //
    //  by Dr. Yavuz Kaya, P.Eng.;
    //  Created on      03.Mar.2019
    //  Last modified   25.Jun.2023
    //  Yavuz.Kaya.ca@gmail.com

    // Check default values of input arguments
    if (ScaleOpt == null) { ScaleOpt = 'none'; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if (!Array.isArray(a)) {throw new Error("a[] must be an array.");};
    if (!Array.isArray(b)) {throw new Error("b[] must be an array.");};

    N1 = a.length;
    N2 = b.length;
    N  = 2 * Math.max(N1, N2) - 1;

    // ScaleOpt must be one of the
    if ((ScaleOpt != 'none') && (N1 != N2)) {throw new Error("a[] and b[] arrays must be the same length.");};

    // Insert zeros at the beginning of a[] array
    if (N1 < N2)  { a = Vector_Concatenate(new Array(N2-1).fill(0), a); }
    else          { a = Vector_Concatenate(new Array(N1-1).fill(0), a); }

    // Calculate the FFT length
    N2FFT = NextPow2(N);

    // Lag vector
    lag = Vector_LinSpace(-(N-1)/2, (N-1)/2, N);

    // Take FFT of both arrays
    [Re1, Im1] = FFT(a, null, N2FFT);
    [Re2, Im2] = FFT(b, null, N2FFT);

    // Calculate product of FFT1 and conjugate(FFT2)
    for (i = 0; i < N2FFT; i++) {
        Re =  Re1[i] * Re2[i] + Im1[i] * Im2[i];
        Im = -Re1[i] * Im2[i] + Im1[i] * Re2[i];

        Re1[i] = Re;
        Im1[i] = Im;
    }

    // Cross-correction: first N element of IFFT(FFT1).Re
    [Re, ] = IFFT(Re1, Im1, N2FFT);
    r = Vector_Truncate(Re, N);

    if (ScaleOpt == 'biased')     { r = Vector_Divide(r, N1); }
    if (ScaleOpt == 'unbiased')   { r = r.map((v,i) => v / (N1-Math.abs(lag[i])) ); }
    if (ScaleOpt == 'normalized') { r1=Vector_Dot(a,a); r2=Vector_Dot(b,b); r = Vector_Divide(r, Math.sqrt(r1*r2)); }

    return [r, lag];
}
//----------------------------------------------------------------------------------------------------------------------
// Statistics

function Hist(data, NumBins) {
    
    // Declaration of variables
    let i, maxVal, minVal, BinEdges, BinWidth, hist;

    // Check default values of input arguments
    if (NumBins == null) { 
        // The Freedman-Diaconis’ Rule
        BinWidth = 2 * IQR(data) / Math.pow(data.length, 1/3);
        NumBins = Math.ceil((Max(data).val - Min(data).val) / BinWidth); 
    };

    maxVal   = Max(data).val;
    minVal   = Min(data).val;
    BinEdges = LinSpace(minVal, maxVal, NumBins+1);
    BinWidth = BinEdges[1] - BinEdges[0];

    hist = new Array(NumBins).fill(0);

    for (i = 0; i < data.length; i++) { hist[BinNumber(data[i], BinEdges)]++; }

    return {
        'Values':    hist,
        'BinEdges':  BinEdges,
        'BinWidth':  BinWidth,
        'NumBins':   NumBins,
        'maxVal':    maxVal,
        'minVal':    minVal
    };

    // Find the bin number for x
    function BinNumber(x, BinEdges) {
        for (let i = 0; i < BinEdges.length - 1; i++) {
            if (i != BinEdges.length - 2) { if ((x >= BinEdges[i]) && (x < BinEdges[i+1])) { return i; }  }
            else { if ((x >= BinEdges[i]) && (x <= BinEdges[i+1])) { return i; } }
        }
    }
}

function IQR(data) {
    // Interquartile range of data[] calculated as the difference between the 75th and the 25th percentiles of the data contained in data[]
    let pr = Percentile(data, [25, 75]);
    return pr[1] - pr[0];
}

//----------------------------------------------------------------------------------------------------------------------
// Windows

function Blackman(M) {
    // Returns a Blackman Window of M-pont

    // Check if the input arguments are correct; otherwise, throw an error.
    if (M <= 0) {throw new Error("Number of points in window cannot be equal or less than zero.");}

    var i, a, h;

    if (M==1) {return new Array(M).fill(1);}

    h = new Array(M).fill(0);
    for (i = 0; i < M; i++) {
        a = 2 * Math.PI * i / (M - 1);

        h[i] = 0.42 - 0.5 * Math.cos(a) + 0.08 * Math.cos(2 * a);
    }
    h[0] = 0;
    h[M - 1] = 0;
    return h;
}

function GaussWin(M, alpha) {
    // Returns a Gaussian Window of M-point
    // alpha is width factor, and it is specified as a positive real scalar.
    // alpha is inversely proportional to the width of the window

    // Check if the input arguments are correct; otherwise, throw an error.
    if (M <= 0)    {throw new Error("Number of points in window cannot be equal or less than zero.");}
    if (alpha < 0) {throw new Error("alpha cannot be less than zero.");}

    var sigma, v;

    sigma =  (M - 1) / 2 / alpha;
    v     = -(M - 1) / 2 - 1;

    return new Array(M).fill().map(() => {v++; return Math.exp(-Math.pow(v, 2) / 2 / sigma / sigma); } );
}

function Hamming(M) {
    // Returns a Hamming Window of M-point

    // Check if the input arguments are correct; otherwise, throw an error.
    if (M <= 0) {throw new Error("Number of points in window cannot be equal or less than zero.");}

    return new Array(M).fill().map((v,i) => 0.54 - 0.46 * Math.cos(2 * Math.PI * i / (M - 1)) );
}

function Hann(M) {
    // Returns a hanning Window of M-point

    // Check if the input arguments are correct; otherwise, throw an error.
    if (M <= 0) {throw new Error("Number of points in window cannot be equal or less than zero.");}

    return new Array(M).fill().map((v,i) => 0.5 * (1 - Math.cos(2 * Math.PI * i / (M - 1))) );
}

function Triang(M) {
    // Returns a triangular Window of M-point

    // Check if the input arguments are correct; otherwise, throw an error.
    if (M <= 0) {throw new Error("Number of points in window cannot be equal or less than zero.");}

    return new Array(M).fill().map((v,i) =>  1 - Math.abs(2 * (i - 0.5 * (M - 1)) / (M + 1)) );
}

function Rectwin(M) {
    // Reruns a Rectangular window of M-point
    return new Array(M).fill(1);
}

//----------------------------------------------------------------------------------------------------------------------
// Signal processing tools

function Cpsd(Input, Output, RWL, OVS, FSamp, NFFT, OneSided) {
    // This function calculates the Cross-Power Spectral Density (Sxy) between a 1D Input array and 2D Output matrix
    // using Welch's method. Each column in the Output matrix is a separate channel.
    // Sampling frequency for the Input and Output is assumed be same.
    // Frequency Response Function (Hxy) between Input and Output matrices is calculated.
    // Impulse Response Function h(t) as a function of time is calculated as the Inverse Fourier Transform of the Hxy.
    //
    // Sxy and Hxy is calculated for each column of the Output matrix.
    // Both the Input array and each column of Output matrix is divided into equal length of segments of RWL samples
    // simultaneously, and each segment is overlapped by OVS samples. Sxy is calculated for each successive segment
    // and averaged across all segments to form the Transfer function. Hxy is calculated as the ratio of Sxy and Sx.
    //
    // Cross power spectral density is the Fourier Transform of the cross-correlation function (Xcorr).
    // Cross-correlation function is a function that defines the relationship between two random signals.
    //
    // Cross power spectral density is a spectral analysis that compares two signals.
    // It gives the total noise power spectral density of two signals. The only condition is that there should be some
    // phase difference or time delay between these two signals.  CPSD analysis is most suitable for studying the effect
    // of stationary, but stochastic signals.
    //
    //  Input:
    //      Input    : 1D Input array
    //      Output   : 2D Output matrix (each column is a separate channel)
    //      RWL      : Number of samples in each data-segment
    //      OVS      : Number of samples that two data segments are overlapped
    //      FSamp    : Sampling frequency (assumed to be same both for Input array and Output matrix
    //      NFFT     : Number of discrete Fourier Transform points in FFT estimate for each data segment
    //      OneSided : 1-sided or 2-sided spectrum
    //                 0 for one-sided  :::  1 for two-sided
    //                 Impulse Response h(t) is a function of time, so h(t) is not affected by OneSided argument
    //
    //  Output:
    //      Sxy     : Cross-power spectral density (same size as Output matrix)
    //      Sx      : Power spectral density of Input array
    //      Sy      : Power spectral density of Output matrix for each column
    //      Cxy     : Magnitude Squared Coherence
    //      Hxy     : Frequency Response Function (same size as Output matrix)
    //      h       : Impulse Response Function of the System (same size as Output matrix)
    //      f       : Frequency vector over which the spectrum is calculated
    //
    //  by Dr. Yavuz Kaya, P.Eng.;
    //  Created on      03.Mar.2019
    //  Last modified   21.jun.2023
    //  Yavuz.Kaya.ca@gmail.com

    var row, col, len, DW1, Win, SF, kk, kk2, df, f, SxyR, SxyI, Sy, Sx, Cxy, HxyR, HxyI, h, tt, i, j, a1, a2, K
    var TempArray, Re_O, Re_I, Ka_I, Ka_O, Re1, Re2, Im1, Im2, R1, I1, hR

    // Determine the size of the Output[] array
    row = Output.length;
    col = Output[0].length;
    len = Input.length;

    // Check default values of input arguments
    if (RWL      == null) { RWL      = Math.max(2, Math.floor(row / 8.00)); };
    if (OVS      == null) { OVS      = Math.floor(RWL * 0.25); };
    if (FSamp    == null) { FSamp    = 1; };
    if (NFFT     == null) { NFFT     = RWL; };
    if (OneSided == null) { OneSided = 1; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if ((OneSided != true) && (OneSided != false)) {throw new Error("OneSided argument must be either true or false."); }
    if ((RWL <= 0) || (RWL > len) || (RWL > NFFT) || (OVS >= RWL) || (OVS < 0) || (NFFT < 2) || (FSamp <= 0)) {
        {throw new Error("Input arguments are inconsistent (RWL <= 0) || (RWL > Input.Length) || (RWL > NFFT) || (OVS >= RWL) || (OVS < 0) || (NFFT < 2) || (FSamp <= 0)");};
    }
    if (!Array.isArray(Input))  {throw new Error("Input[] must be 1D array.");};
    if (!Array.isArray(Output)) {throw new Error("Output[] must be 2D array.");};
    if (row < 2)                {throw new Error("Number of elements in each column of Output[] matrix cannot be less than 2.");};
    if (len < 2)                {throw new Error("Number of elements in Input[] array cannot be less than 2.");};
    if (row != len)             {throw new Error("Dimensions of Input[] array and Output[] matrix is inconsistent.");};

    // Calculate initial parameters
    DW1 = RWL - OVS;

    // Create hamming window of RWL-length and normalize it
    Win = Hamming(RWL);
    Win = Vector_Divide(Win, Statistics_Sum(Win));

    // Scale Factor for the PSD normalization
    // Units of PSD is "Hertz" because it is scaled by the sampling frequency to obtain the PSD
    SF = Statistics_Sum(Vector_Pow(Win, 2)) * FSamp;

    // Length of spectrum vector based on OneSided variable (one-sider or two-sided)
    if (OneSided) {
        // One-sided spectrum
        if (NFFT % 2 == 0) { kk = NFFT / 2 + 1;  kk2 = 2; } else { kk = (NFFT + 1) / 2; kk2 = 1;}

        // Adjust the scale factor for one-sided spectrum
        SF /= 2;

        // Frequency vector
        df = FSamp / NFFT;
        f  = Vector_LinSpace(0, (kk-1)*df, kk);

    } else {
        // Two-sided spectrum
        kk = NFFT;

        // Frequency vector
        df = FSamp / NFFT;
        f  = Vector_LinSpace(0, (NFFT-1)*df, NFFT);
    }

    // Pre-allocation for Power Spectral Density (PSD)
    SxyR = new Array(kk).fill(0).map((v,i)   => new Array(col).fill(0));
    SxyI = new Array(kk).fill(0).map((v,i)   => new Array(col).fill(0));
    Sy   = new Array(kk).fill(0).map((v,i)   => new Array(col).fill(0));
    Sx   = new Array(kk).fill(0);
    Cxy  = new Array(kk).fill(0).map((v,i)   => new Array(col).fill(0));
    HxyR = new Array(kk).fill(0).map((v,i)   => new Array(col).fill(0));
    HxyI = new Array(kk).fill(0).map((v,i)   => new Array(col).fill(0));
    h    = new Array(NFFT).fill(0).map((v,i) => new Array(col).fill(0));

    for (j = 0; j < col; j++) {

        a1 = 0;
        a2 = RWL - 1;
        K = 0;

        // Reset the Sx array
        for (i = 0; i < kk; i++) { Sx[i] = 0; }

        // Obtain the data in each column
        TempArray = Matrix_GetColumn(Output, j);

        while (a2 < TempArray.length) {
            // Extract the segment of the Output[] array starting from Index a1 to Index a2
            Re_O = Vector_GetRange(TempArray, a1, a2);
            Ka_O = IsNaN(Re_O);

            Re_I = Vector_GetRange(Input, a1, a2);
            Ka_I = IsNaN(Re_I);

            // Skip one iteration in the for-loop if Output[] array OR Input[] array contains NaN
            if (Ka_O || Ka_I) { a1 += DW1; a2 += DW1; continue; };

            // Fourier Amplitude Spectrum (FAS) of the windowed-segment
            [Re1, Im1] = FFT(Vector_Multiply(Re_O, Win), null, NFFT);
            [Re2, Im2] = FFT(Vector_Multiply(Re_I, Win), null, NFFT);

            for (tt = 0; tt < kk; tt++) {
               SxyR[tt][j] +=  Re2[tt] * Re1[tt] + Im2[tt] * Im1[tt];
               SxyI[tt][j] += -Re2[tt] * Im1[tt] + Im2[tt] * Re1[tt];
               Sx[tt]      +=  Re2[tt] * Re2[tt] + Im2[tt] * Im2[tt];
               Sy[tt][j]   +=  Re1[tt] * Re1[tt] + Im1[tt] * Im1[tt];
            }

            // Update the Indexes
            a1 += DW1;
            a2 += DW1;
            K += 1;
        }

        // Normalization of the averaged Power Spectrum to calculate PSD
        // Units of PSD is Hertz because it is scaled by the sampling frequency to obtain the psd
        Sx = Vector_Divide(Sx, (SF * K));

        for (tt = 0; tt < kk; tt++) {
            SxyR[tt][j] /= (SF * K);
            SxyI[tt][j] /= (SF * K);
            Sy[tt][j]   /= (SF * K);

            Cxy[tt][j] = (SxyR[tt][j] * SxyR[tt][j] + SxyI[tt][j] * SxyI[tt][j]) / (Sx[tt] * Sy[tt][j]);

            HxyR[tt][j] = SxyR[tt][j] / Sx[tt];
            HxyI[tt][j] = SxyI[tt][j] / Sx[tt];
        }

        if (OneSided) {
            SxyR[0][j] /= 2;
            SxyI[0][j] /= 2;
            Sy[0][j]   /= 2;
            Sx[0]      /= 2;

            if (NFFT % 2 == 0) {
                // Odd
                SxyR[kk-1][j] /= 2;
                SxyI[kk-1][j] /= 2;
                Sy[kk-1][j]   /= 2;
                Sx[kk-1]      /= 2;
            }
            R1 = Vector_Concatenate(Matrix_GetColumn(HxyR, j), Vector_Flip(Vector_GetRange(Matrix_GetColumn(HxyR, j), 1, kk-kk2)));
            I1 = Vector_Concatenate(Matrix_GetColumn(HxyI, j), Vector_Multiply(Vector_Flip(Vector_GetRange(Matrix_GetColumn(HxyI, j), 1, kk-kk2)), -1));
        } else {
            R1 = Matrix_GetColumn(HxyR, j);
            I1 = Matrix_GetColumn(HxyI, j);
        }

        [hR, ] = IFFT(R1, I1);
        hR = Vector_Concatenate([hR[0]], Vector_Flip(Vector_GetRange(hR, 1, NFFT - 1)));

        for (tt = 0; tt < NFFT; tt++) {
            h[tt][j] = hR[tt];
        }
    }

    return [SxyR, SxyI, Sy, Sx, Cxy, HxyR, HxyI, h, f];

}

function Filter(b, a, x, zi) {
    // This function filters the x[] array to create the filtered data of y[] array using the filter coefficients
    // described in b[] and a[] arrays.
    //
    // The filter is a "Direct Form II Transposed" implementation of the standard difference equation as follow:
    //
    //   a[1]*y[n] = b[1]*x[n] + b[2]*x[n-1] + ... + b[nb+1]*x[n-nb]
    //                         - a[2]*y[n-1] - ... - a[na+1]*y[n-na]
    //
    //  b[], a[]    : Filter coefficients of numerator and denominator, respectively
    //       x[]    : Raw data to be filtered
    //        zi    : Filter initial conditions
    //
    //       y[]    : Filtered data
    //        zf    : Filter final conditions
    //
    // This function is cross checked with Matlab [y, zf] = FILTER(b, a, x, zi) built-in function. 202306-20

    //
    var i, j, y, tmp, zf, N, N1, N2, N3;

    N = x.length;
    N1 = a.length;
    N2 = b.length;
    N3 = Math.max(N1, N2) - 1;

    // Check default values of input arguments
    if (zi == null) { zi = new Array(N3).fill(0); };

    // Check if the input arguments are correct; otherwise, throw an error.
    if (!Array.isArray(a))  {throw new Error("Filter coefficients a[] must be an array.");};
    if (!Array.isArray(b))  {throw new Error("Filter coefficients b[] must be an array.");};
    if (!Array.isArray(x))  {throw new Error("Raw data x[] must be an array.");};
    if (!Array.isArray(zi)) {throw new Error("Filter initial conditions zi[] must be an array.");};
    if (zi.length != N3)    {throw new Error("Dimension of initial conditions zi[] array is inconsistent.");}

    // Calculate the filtered data
    y = new Array(N).fill(0);
    for (i = 0; i < N; ++i) {
        tmp = 0.0;
        for (j = 0; j < N2; ++j) {
            if (i - j < 0) {continue;}
            tmp += b[j] * x[i - j];
        }
        for (j = 1; j < N1; ++j) {
            if (i - j < 0) {continue;}
            tmp -= a[j] * y[i - j];
        }

        // apply filter initial conditions
        if (i < N3) { tmp += zi[i]; }

        tmp /= a[0];
        y[i] = tmp;
    }

    // Calculate filter final conditions
    zf = FilterFinalCondition(b, a, x, y);

    return {"y": y,   "zf": zf};
}

function FiltFilt(b, a, data) {

    // This function performs zero-phase digital filtering by processing the input array of data[] in both the forward
    // and reverse directions. After filtering the data in the forward direction, the function reverses the filtered
    // sequence and runs it back through the filter. The result has these characteristics:
    //    Zero phase distortion.
    //    A filter transfer function equal to the squared magnitude of the original filter transfer function.
    //    A filter order that is double the order of the filter specified by b and a.
    //
    // FiltFilt minimizes start-up and ending transients by matching initial conditions.
    // Do not use FiltFilt with differentiator and Hilbert FIR filters, because the operation of these filters
    // depends heavily on their phase response.
    //

    // Check if the input arguments are correct; otherwise, throw an error.
    if (!Array.isArray(a)) {throw new Error("Filter coefficients a[] must be an array.");};
    if (!Array.isArray(b)) {throw new Error("Filter coefficients b[] must be an array.");};

    var i, z, zi, N, M, x, K1, K2, y;

    // Calculate filter initial conditions for zero-phase response
    zi = Filter_Zi_Initials(b, a);

    //Pad the data with 3 times the filter length on both sides by rotating the data by 180°
    N = 3 * (a.length - 1);
    M = N + data.length;
    x = new Array(data.length + 2 * N).fill(0);

    for (i = 0; i < data.length; i++) { x[N+i] = data[i]; }

    K1 = 2 * x[N];
    K2 = 2 * x[M - 1];

    for (i = 0; i < N; i++) {
        x[N - i - 1] = K1 - x[N + i + 1];
        x[M + i] = K2 - x[M - i - 2];
    }

    //Calculate the initial values for the given sequence
    z = Mult_El(zi, x[0]);

    // Forward filtering
    y = Filter(b, a, x, z).y;

    // Revers the filtered data
    y = Flip(y);

    // Calculate the initial values for the filtered-sequence
    z = Mult_El(zi, y[0]);

    // Backward filtering
    y= Filter(b, a, y, z).y;

    // Revers the filtered data
    y = Flip(y);
    y = ExtractSubset(y, N, M - 1);

    return {"y": y,   "zf": zi};

}

function FilterFinalCondition(b, a, x, y) {
    // Calculates the filter digital-final conditions
    var i, n, ny, N, t, zf, say, sum;

    ny = y.length;
    N = b.length;
    t = 1;

    zf = new Array(N-1).fill(0);

    for (i = 0; i < N - 1; i++) {
        say = 0;
        sum = 0;
        for (n = ny - 1; n > (ny - N + t - 1); n--) {
            if (n < 0) { continue; }
            sum += x[n] * b[say + t] - y[n] * a[say + t];
            say++;
        }
        t++;
        zf[i] = sum;
    }
    return zf;
}

function Filter_Zi_Initials(b, a) {
    // Calculates the filter initial conditions for Filt-Filt filter for zero-phase response

    var i, j, N, A1, A2, A3, A4, zi, sum;

    N = a.length - 1;
    A1 = Eye(N);
    A2 = new Array(N).fill(0).map(() => new Array(N).fill(0));
    A3 = new Array(N).fill(0);
    zi = new Array(N).fill(0);
    sum = 0;

    for (i = 0; i < N; i++) {
       A2[i][0] = -a[i + 1];
    }
    for (i = 0; i < N - 1; i++) {
       A2[i][i + 1] = 1;
    }
    for (i = 0; i < N; i++) {
        A3[i] = b[i + 1] - b[0] * a[i + 1];
    }
    A4 =  Inv(Subtract(A1, A2));
    for (i = 0; i < N; i++) {
        sum = 0;
        for (j = 0; j < N; j++) {
            sum += A4[i][j] * A3[j];
        }
        zi[i] = sum;
    }
    return zi;
}

function Freq(b, a, f0, FSamp) {
    // Frequency response of digital filter at f=f0 Hz
    var dig_w, s, i, na, sum1, sum2;
    dig_w = (2 * Math.PI * f0 / FSamp);
    s = Exp(new ComplexNum(0, dig_w));
    na = a.length;

    sum1 = new ComplexNum(0, 0);
    sum2 = new ComplexNum(0, 0);

    for (i = 0; i < na; i++) {
        sum1 = Add(sum1, Mult(Pow(s, na-1-i), b[i]));
        sum2 = Add(sum2, Mult(Pow(s, na-1-i), a[i]));
    }
    return Divide(sum1, sum2);
}

function Freqz(b, a, N, FSamp, option) {
    //   FREQZ Frequency response of digital filter
    //   [H, W] = FREQZ(B, A, N) returns the N-point complex frequency response
    //   vector H and the N-point frequency vector W in radians/sample of the filter:
    //
    //                jw               -jw              -jmw
    //        jw  B(e)     b(1) + b(2)e + .... + b(m+1)e
    //     H(e) = ----- = ------------------------------------
    //               jw               -jw              -jnw
    //            A(e)     a(1) + a(2)e + .... + a(n+1)e
    //
    //   given numerator and denominator coefficients in vectors B and A.
    //
    //   [H, F] = FREQZ(B, A, N, FSamp) returns the N-point complex frequency response
    //   vector H and the N-point frequency vector F (in Hz) of the filter.
    //   FSamp is the sampling frequency in Hz.
    //
    //   Option : "WHOLE" or "TWOSIDED" uses N-points around the whole unit circle
    //   (Results are validated using MATLAB)
    let i, j, H, na, f, w, s, sum1, sum2;

    // Check arguments for default values 
    if (N == null) { N = 512; }
    if (FSamp == null) { FSamp = 1; }
    if (option == null) { option = "ONESIDED"; }

    H = new Array(N).fill().map((x) => new ComplexNum(0,0));
    na = a.length;

    // Remove whitespace and convert to UpperCase
    option = option.replace(/\s/g, "").toUpperCase()

    if ((option == "WHOLE") || (option == "TWOSIDED")) {
        f = LinSpace(0, (FSamp - FSamp / N), N);
        w = new Array(f.length).fill().map((v,i) => f[i] * (2 * Math.PI / FSamp));
    } else {
        FSamp /= 2;
        f = LinSpace(0, (FSamp - FSamp / N), N);
        w = new Array(f.length).fill().map((v,i) => f[i] * (Math.PI / FSamp));
    }
    for (j = 0; j < N; j++) {
        s = Exp(new ComplexNum(0,w[j]));

        sum1 = new ComplexNum(0, 0);
        sum2 = new ComplexNum(0, 0);

        for (i = 0; i < na; i++) {
            sum1 = Add(sum1, Mult_El(Pow(s, na - 1 - i), b[i]));
            sum2 = Add(sum2, Mult_El(Pow(s, na - 1 - i), a[i]));
        }
        H[j] = Divide(sum1, sum2);
    }
    return {"H": H,  "f": f};
}

function tf2zp() {

    let a = [1,2,3,4,5];
    let b = [3,4,5];

    //[b, a] = tfchk(b, a);

    m = b.length;
    n = a.length;

    if (a.length == 0) { coeff = a[0]; } else { coef = 1; }
    if (coef == 0) { console.log("ERROR-1");  return; }



    return [b, a];



    function tfchk(b, a) {
        let m, n;
        m = b.length;
        n = a.length;
        if (n-m > 0) { b = new Array(n-m).fill(0).concat(b); }
        return [b, a];
    }

}

//----------------------------------------------------------------------------------------------------------------------
// Butterworth Filter Design
function Butterworth_BandPass(N, fc1, fc2, FSamp) {
    let k, theta, b, a, a1, a2, zf, K, Fc1, Fc2, BW_hz, F0, f0, alpha, beta, h, t1, t2;
    let pa = new Array(2 * N).fill();
    let pl = new Array(2 * N).fill();
    let q  = new Array(2 * N).fill();
    let p  = new Array(2 * N).fill().map((x) => new ComplexNum(0,0));

    // continuous pre-warped frequency
    Fc1 = FSamp / Math.PI * Math.tan(Math.PI * fc1 / FSamp);
    Fc2 = FSamp / Math.PI * Math.tan(Math.PI * fc2 / FSamp);

    // Hz prewarped - 3 dB bandwidth
    BW_hz = Fc2 - Fc1;

    // Hz geometric mean frequency
    F0 = Math.sqrt(Fc1 * Fc2);

    t1 = new ComplexNum(2 * Math.PI * F0, 0);
    t2 = new ComplexNum(0, 1);

    for (k = 0; k < N; k++) {
        theta = (2 * k + 1) * Math.PI / (2 * N);

        // pa are poles of analog filter with cutoff frequency of 1 rad/s
        pl = new ComplexNum(-Math.sin(theta), Math.cos(theta));

        alpha = Mult(BW_hz / F0 * 0.5, pl);
        beta = Sqrt(Add(Mult(Pow(alpha, 2), -1), 1));

        pa[k] = Mult(Add(Mult(beta, t2), alpha), t1);
        pa[k + N] = Mult(Add(Mult(Mult(beta, t2), -1), alpha), t1);

        // poles and zeros in the z plane
        a1 = Divide(pa[k], 2 * FSamp);
        a2 = Divide(pa[k + N], 2 * FSamp);

        p[k] = Divide(Add(a1, 1), Add(Mult(a1, -1), 1));
        p[k + N] = Divide(Add(a2, 1), Add(Mult(a2, -1), 1));

        q[k] = -1;
        q[k + N] = 1;
    }
    
    // convert poles to polynomial coefficients, a[]
    a = Real(Poly(p));

    // convert zeros to polynomial coefficients, b[]
    b = Real(Poly(q));

    // Scale coefficients so that amplitude is 1.0 at f0
    f0 = Math.sqrt(fc1 * fc2);   /// Fc1 or fc1  not sure ???????????????????????????????????????

    h = Freq(b, a, f0, FSamp); 

    // Calculate amplitude scale factor
    b = Divide(b, Abs(h)); 

    // Empty array - final conditions of the filter
    zf = new Array(b.length - 1).fill(0);

    return {"b": b, "a": a, "zf": zf};
}

function Butterworth_HighPass(N, fc, FSamp) {
    let i, k, pa, theta, b, a, a1, sum1, sum2, K, Fc, zf;
    let q = new Array(N).fill();
    let p = new Array(N).fill().map((x) => new ComplexNum(0, 0));

    // continuous pre-warped frequency
    Fc = FSamp / Math.PI * Math.tan(Math.PI * fc / FSamp);

    for (k = 0; k < N; k++) {
        theta = (2 * k + 1) * Math.PI / (2 * N);

        // pa are poles of analog filter with cutoff frequency of 1 rad/s
        pa = new ComplexNum(-Math.sin(theta), Math.cos(theta));

        // scale poles (in frequency) by (2 * pi * Fc) to transform from analog LP to analog HP
        pa = Divide(new ComplexNum(2 * Math.PI * Fc, 0), pa)

        // poles and zeros in the z plane
        a1 = Divide(pa, 2 * FSamp)

        p[k] = Divide(Add(1, a1), Subtract(1, a1));
        q[k] = 1;
    }

    // convert poles to polynomial coefficients, a[]
    a = Real(Poly(p));

    // convert zeros to polynomial coefficients, b[]
    b = Real(Poly(q));

    // Amplitude scale factor for gain = 1 at f = fs/2  (z = -1)
    sum1 = 0;
    sum2 = 0;
    for (i = 0; i < N + 1; i++) {
        sum1 += Math.pow(-1, i) * a[i];
        sum2 += Math.pow(-1, i) * b[i];
    }
    K = sum1 / sum2;

    b = Mult(b, K); 

    // Empty array - final conditions of the filter
    zf = new Array(b.length - 1).fill(0);

    return {"b": b, "a": a, "zf": zf};

}

function Butterworth_LowPass(N, fc, FSamp) {

    let k, pa, b, a, zf, a1, K, Fc, theta;
    let q = new Array(N).fill();
    let p = new Array(N).fill().map((x) => new ComplexNum(0, 0));

    // continuous pre-warped frequency
    Fc = FSamp / Math.PI * Math.tan(Math.PI * fc / FSamp);

    for (k = 0; k < N; k++) {
        theta = (2 * k + 1) * Math.PI / (2 * N);

        // pa are poles of analog filter with cutoff frequency of 1 rad/s
        pa = new ComplexNum(-Math.sin(theta), Math.cos(theta));
        
        // scale poles (in frequency) by 2 * pi * Fc
        pa = Mult(pa, 2 * Math.PI * Fc);

        // poles and zeros in the z plane
        a1 = Divide(pa, 2 * FSamp);

        p[k] = Divide(Add(1, a1), Subtract(1, a1));
        q[k] = -1;
    }

    // convert poles p[] to polynomial coefficients, a[]
    a = Real(Poly(p));

    // convert zeros q[] to polynomial coefficients, b[]
    b = Real(Poly(q));
   
    // Calculate amplitude scale factor
    K = Divide(Sum(a), Sum(b)); 

    b = Mult(b, K); 
  
    // Empty array - final conditions of the filter
    zf = new Array(b.length - 1).fill(0);

    return {"b": b, "a": a, "zf": zf};
}

function Cheby1(Order, Rp, Wn) {

    let z, p, k, a,b,c,d;

    //  Get N-th order Chebyshev type-I lowpass analog prototype
    [z, p, k] = cheb1ap(Order, Rp);

    // Transform to state-space
    let tt  = zpk2ss(z,p,k);
    
    return [a,b,c,d];

    // N-th order Chebyshev type-I lowpass analog prototype
    function cheb1ap(n, Rp) {
        
        // Decleration of variables 
        let epsilon, mu, p, realp, imagp, k, z=[];

        epsilon = Math.sqrt( 10 ** (0.1 * Rp) - 1);
        mu      = Math.asinh(1 / epsilon) / n;
        p       = Exp(Num2Complex(0, Add( Mult( LinStep(1, 2*n-1, 2), Math.PI / 2 / n), Math.PI/2 )));
        realp   = Real(p);    realp = Divide(Add(realp, Flip(realp)),2);
        imagp   = Imag(p);    imagp = Divide(Subtract(imagp, Flip(imagp)),2);  
        p       = Num2Complex( Mult(realp, Math.sinh(mu)),  Mult(imagp, Math.cosh(mu)));
        k       = Real(Prod(Mult(p, -1)));
        if (IsEven(n)) {k = Divide(k, Math.sqrt(1+ epsilon*epsilon)); }
        return [z, p, k];
    }

}

function zpk2ss(z, p, k) {

    // Decleration of variables 
    let a, b, c, d, num, np, nz, den, wn, t, i, index, a1, b1, c1, d1, ma1, na2;

    p = GetRange(p, Find(IsFinite(p)));
    z = GetRange(z, Find(IsFinite(z)));
    
    // Group into complex pairs
    np = p.length;
    nz = z.length;
    
    try {
        z = cplxpair(z);
        p = cplxpair(p);
    } finally {
        z = cplxpair(z, 1e6*nz*Norm(z)*Number.EPSILON+Number.EPSILON);
        p = cplxpair(p, 1e6*nz*Norm(z)*Number.EPSILON+Number.EPSILON);
    }

    // Initialize state-space matrices for running series
    a=[]; b=[]; c=[]; d=1;

    //  If odd number of poles AND zeros, convert the pole and zero at the end into state-space.
    //  H(s) = (s-z1)/(s-p1) = (s + num(2)) / (s + den(2))
    if (!IsEven(np) && !IsEven(nz)) {
        console.log("BIR")
        a = p[np-1];    b = 1;    c = p[np-1] - z[nz-1];    d = 1;
        np = np - 1; 
        nz = nz - 1;
    }
    
    // If odd number of poles only, convert the pole at the end into state-space.
    //  H(s) = 1/(s-p1) = 1/(s + den(2)) 
    if (!IsEven(np)) {
        console.log("IKI")
        a = p[np-1];    b = 1;    c = 1;    d = 0;
        np = np - 1;
    }
    
    //  If odd number of zeros only, convert the zero at the end, along with a pole-pair into state-space.
    //   H(s) = (s+num(2))/(s^2+den(2)s+den(3)) 
     if (!IsEven(nz)) {
        console.log("UC")
        num = Real(Poly(z[nz-1])); 
        den = Real(Poly(ExtractSubset(p, np-2, np-1)));
        wn  = Math.sqrt(Prod(Abs(ExtractSubset(p, np-2, np-1))));
        if (wn == 0) { wn = 1; } 
        t   = Diag( [1, 1/wn] ); // Balancing transformation
        a = Mult(GaussElimination(t, [[-den[1], -den[2]], [1, 0]]), t); 
        b = Transpose( GaussElimination(t, [1, 0]) );    
        c = Mult( [1, num[1]], t);
        d = 0;
        nz = nz - 1;
        np = np - 2;
    }
    
    // Now we have an even number of poles and zeros, although not necessarily the same number - there may be more poles.
    //   H(s) = (s^2+num(2)s+num(3))/(s^2+den(2)s+den(3))
    // Loop through rest of pairs, connecting in series to build the model.
    i = 0;
    while (i < nz) {
        
        index = LinStep(i, i+1, 1);
        num   = Real(Poly( GetRange(z, index) )); 
        den   = Real(Poly( GetRange(p, index) )); 
        wn    = Math.sqrt(Prod(Abs( GetRange(p, index) )));
        if (wn == 0) { wn = 1; } 
        t     = Diag([1, 1/wn]);  // Balancing transformation
        a1    = Mult(GaussElimination(t, [[-den[1], -den[2]], [1, 0]]), t);
        b1    = Transpose( GaussElimination(t, [1, 0]) ); 
        c1    = Mult( [num[1]-den[1], num[2]-den[2]], t);
        d1    = 1;
        // Next lines perform series connection 
        ma1   = a.length;
        na2   = a1[0].length;
        if   (isNumber(a)) { a = Concat(Concat(a, Zeros(na2),     1),   Concat(Mult(b1,c),   a1, 1), 0); } 
        else               { a = Concat(Concat(a, Zeros(ma1,na2), 1),   Concat(Mult(b1,[c]), a1, 1), 0); }
        b     = Concat(b, Mult(b1, d), 0);
        c     = Concat(Mult(c, d1), c1, 1);
        d     = Mult(d1, d);
        i     = i + 2;
        
    }

    // Take care of any left over unmatched pole pairs.
    //   H(s) = 1/(s^2+den(2)s+den(3))
    while (i < np) {
        index = LinStep(i, i+1, 1);
        den   = Real(Poly( GetRange(p, index) ));
        wn    = Math.sqrt(Prod(Abs(  GetRange(p, index) )));
        if (wn == 0) { wn = 1; } 
        t     = Diag([1, 1/wn]);  // Balancing transformation
        a1    = Mult(GaussElimination(t, [[-den[1], -den[2]], [1, 0]]), t);
        b1    = Transpose( GaussElimination(t, [1, 0]) );
        c1    = Mult( [0, 1], t);
        d1    = 0;
        // Next lines perform series connection 
        ma1   = a.length;
        na2   = a1[0].length;
        if   (isNumber(a)) { a = Concat(Concat(a, Zeros(na2),     1),   Concat(Mult(b1,c),   a1, 1), 0); } 
        else               { a = Concat(Concat(a, Zeros(ma1,na2), 1),   Concat(Mult(b1,[c]), a1, 1), 0); }
        b     = Concat(b, Mult(b1, d), 0);
        c     = Concat(Mult(c, d1), c1, 1);
        d     = Mult(d1, d);
        i     = i + 2;
    }

    // Apply gain k:
    c = Mult(c, k);
    d = Mult(d, k);

    return {a,b,c,d}
}

let z = [1,2,30,54,98];
let p = [ new ComplexNum(10.23,  90.45),
          new ComplexNum(10.23, -90.45),
          new ComplexNum(20.88,  80.55),
          new ComplexNum(20.88, -80.55),
          new ComplexNum(30.99,  70.99),
          new ComplexNum(30.99, -70.99),
          new ComplexNum(40.15,  60.15),
          new ComplexNum(40.15, -60.15),   12,25,3,18
        ];
let k = 0.3548;

let ss= zpk2ss(z, p, k)
console.log(ss.a)
console.log(ss.b)
console.log(ss.c)
console.log(ss.d)

function cplxpair(p, tol) {

    // Decleration of variables
    let i, j, complexNumbers=[],  realNumbers=[];
    let sortedComplex= [], used, current, conjugateFound, realDiff, imDiff, isConjugate, Message;

    // Check input variable
    if (tol == null) { tol = 100 * Number.EPSILON; }

    // Separate complex and real numbers
    p.forEach(num => {
        if (typeof(num) === 'object' && Math.abs(num.Im) > tol) { complexNumbers.push(num); } else { realNumbers.push(num); }
    });

    // Check for even number of complex numbers
    if (complexNumbers.length % 2 !== 0) { 
        Message = "Number of complex numbers must be even to form conjugate pairs";
        throw new Error(Message);
    }

    // Sort complex numbers into conjugate pairs
    used = new Set();

    for (i = 0; i < complexNumbers.length; i++) {
        if (used.has(i)) continue;

        current        = complexNumbers[i];
        conjugateFound = false;

        // Find conjugate pair
        for (j = i + 1; j < complexNumbers.length; j++) {
            if (used.has(j)) continue;

            // Check if numbers are conjugates within tolerance
            realDiff    = Math.abs(current.Re - complexNumbers[j].Re);
            imDiff      = Math.abs(Math.abs(current.Im) - Math.abs(complexNumbers[j].Im));
            isConjugate = realDiff < tol && imDiff < tol && Math.sign(current.Im) !== Math.sign(complexNumbers[j].Im);

            if (isConjugate) {
                // Add pair with negative imaginary part first
                if (current.Im < 0) { sortedComplex.push(current, complexNumbers[j]); } 
                else                { sortedComplex.push(complexNumbers[j], current); }
                used.add(i);
                used.add(j);
                conjugateFound = true;
                break;
            }
        }
        if (!conjugateFound) {
            Message = "No conjugate pair found for complex number: " + JSON.stringify(current); 
            throw new Error(Message);
        }
    }
    
    // Check if sortedComplex is empty
    if (sortedComplex.length > 0) { 
        let Ind       = Sort(Real(sortedComplex)).Ind;  
        sortedComplex = GetRange(sortedComplex, Ind);
    } else {
        sortedComplex = [];
    }

    // Combine complex pairs and real numbers
    return Concat(sortedComplex, Sort(realNumbers, "ASC").data);
}

function Chebyshev_Filter() {

    let n, Rp, fc, Wn, ftype;
    let fs, u;

    n     = 6;
    Rp    = 10;
    fs    = 1000;
    Wn    = [0.6, 0.8];
    btype = 2;

    u     = [];

    Bw = [];
    // Convert to low-pass prototype estimate
    if (btype == 1)	{
        // lowpass
        fs = 2;
        u  = 2 * fs * Math.tan(Math.PI * Wn / fs);
        Wn = u;
    }
    else if (btype == 2) {
        // bandpass
        fs   = 2;
        u[0] = 2 * fs * Math.tan(Math.PI * Wn[0] / fs);
        u[1] = 2 * fs * Math.tan(Math.PI * Wn[1] / fs);

        Bw = u[2] - u[1];
        Wn = Math.sqrt( u[0] * u[1] );	// center frequency
    }
    else if (btype == 3) {
        // highpass
        fs = 2;
        u  = 2 * fs * Math.tan(Math.PI * Wn / fs);
        Wn = u;
    }


    function cheb1ap(n, rp) {
        // Chebyshev Type I analog lowpass filter prototype
        // Returns the zeros, poles, and gain of an N-th order normalized analog prototype Chebyshev Type I lowpass filter
        // with Rp decibels of ripple in the pass-band.  Chebyshev Type I filters are maximally flat in the stop-band.

        // Cast to enforce precision rules
        // n = double(n);
        // rp = double(rp);

        epsilon = Math.sqrt( Math.pow(10, (0.1 * rp)) - 1 );
        mu = Math.asinh(1 / epsilon) / n;

        let yy = Vector_LinIncrement(1, 2*n-1, 2);
        yy = Vector_Multiply(yy, Math.PI/(2*n));
        yy = Vector_Add(yy, Math.PI / 2);

        //new Complex(0, yy)

        //p = Math.exp( 1i *(Math.PI * (1:2:2 * n-1) / (2*n) + Math.PI / 2) ).';
        //realp = real(p); realp = (realp + flipud(realp))./2;
        //imagp = imag(p); imagp = (imagp - flipud(imagp))./2;
        //p = complex(sinh(mu).*realp , cosh(mu).*imagp);
        //z = [];
        //k = real(prod(-p));
        //if ~rem(n,2)	% n is even so patch k
        //    k = k/sqrt((1 + epsilon^2));
        //end


    }

}

function Chebyshev_Filter_OLD() {

    let A, B, TA, TB;

    FC = 0.6;  // Enter cutoff frequency (0 to .5)
    LH = 0;    // Enter 0 for Low-Pass  OR enter 1 for High-Pass filter
    PR = 1;    // Enter percent ripple (0 to 29): ", PR
    N  = 6;    // Enter number of poles (2,4,...20): ", NP

    A = new Array(22).fill(0);    A[2] = 1;
    B = new Array(22).fill(0);    B[2] = 1;

    TA = new Array(22).fill(0);
    TB = new Array(22).fill(0);

    //LOOP FOR EACH POLE-PAIR
    for (P = 1; P < N/2; P++) {
        [A0, A1, A2, B1, B2] = ABC(FC, LH, PR, P);

        // Add coefficients to the cascade
        for (I = 0; I < 22; I++) {
            TA[I] = A[I];
            TB[I] = B[I];
        }

        for (I = 2; I < 22; I++) {
            A[I] = A0 * TA[I] + A1 * TA[I - 1] + A2 * TA[I - 2];
            B[I] = TB[I] - B1 * TB[I - 1] - B2 * TB[I - 2];
        }
    }

    B[2] = 0; // Finish combining coefficients
    for (I = 0; I < 20; I++){
        A[I] =  A[I + 2];
        B[I] = -B[I + 2];
    }

    SA = 0;  //NORMALIZE THE GAIN
    SB = 0;
    for (I = 0; I < 20; I++){
        if (LH = 0) { SA += A[I];                      SB +=B[I]; }
        if (LH = 1) { SA += A[I] * Math.pow(-1, I);    SB += B[I] * Math.pow(-1, I);}
    }

    GAIN = SA / (1 - SB);
    for (I = 0; I < 20; I++) { A[I] /= GAIN; }

    return [B, A];


    function ABC(FC, LH, PR, P) {
        RP = -Math.cos(Math.PI / (N * 2) + (P - 1) * Math.PI / N);
        IP =  Math.sin(Math.PI / (N * 2) + (P - 1) * Math.PI / N);

        //Warp from a circle to an ellipse
        if (PR != 0) {
            ES = Math.sqrt( Math.pow(100 / (100 - PR), 2) - 1 );
            VX = (1 / N) * Math.log( (1 / ES) + Math.sqrt( (1/ES*ES) + 1) );
            KX = (1 / N) * Math.log( (1 / ES) + Math.sqrt( (1/ES*ES) - 1) );
            KX = (Math.exp(KX) + Math.exp(-KX)) / 2;
            RP = RP * ( (Math.exp(VX) - Math.exp(-VX) ) /2 ) / KX;
            IP = IP * ( (Math.exp(VX) + Math.exp(-VX) ) /2 ) / KX;
        }

        T = 2 * Math.tan(1/2);
        W = 2 * Math.PI * FC;
        M = RP * RP + IP * IP;
        D = 4 - (4 * RP * T) + (M * T * T);
        X0 = (T * T) / D;
        X1 = 2 * T * T / D;
        X2 = T * T / D ;
        Y1 = (8 - 2 * M * T *T) / D;
        Y2 = (-4 - (4 * RP * T) - (M * T * T)) / D;

        if (LH = 1) { K = -Math.cos(W / 2 + 1/2) / Math.cos(W / 2 - 1/2); }
        if (LH = 0) { K =  Math.sin(1/2 - W / 2) / Math.sin(1/2 + W / 2); }

        D = 1 + (Y1 * K) - (Y2 * K * K);
        A0 = (X0 - X1 * K + X2 * K * K) / D;
        A1 = ((-2 * X0 * K) + X1 + (X1 * K * K) - (2 * X2 * K)) / D;
        A2 = ((X0 * K * K) - (X1 * K) + X2) / D;
        B1 = (2*K + Y1 + (Y1 * K * K) - (2 * Y2 * K)) / D;
        B2 = (-(K * K) - (Y1 * K) + Y2) / D;
        if (LH = 1) { A1 = -A1; B1 = -B1}

        return [A0, A1, A2, B1, B2];
    }



}
//----------------------------------------------------------------------------------------------------------------------
// Fourier Transform

function FFT(REX, IMX, NFFT) {

    let N = REX.length
    let Re=[], Im=[], wR=[], wI=[];
    let i, Ind, N2FFT;

    // NFFT is the length of frequency vecor and equals to the length of REX[] if not given
    if (NFFT == null) { NFFT = N; };

    // Create Re[] and Im[] arrays of length NFFT
    Re = new Array(NFFT).fill(0);
    Im = new Array(NFFT).fill(0);

    // Last Index to be coppied 
    if (N < NFFT) {Ind=N;} else {Ind=NFFT;}

    // Copy from REX[] and IMX[] into Re[] and Im[]
    for (i=0; i<Ind; i++) { Re[i]=REX[i]; }
    if  (IMX != null)     { for (i=0; i<Ind; i++) { Im[i]=IMX[i]; } }

    // Determine the next power of 2
    N2FFT = NextPow2(NFFT);

    // If NFFT is power of 2, then comupte in-place FFT and return.
    if (NFFT == N2FFT) { FFT_Algorithm(Re, Im, NFFT); return [Re, Im]; }

    // Twiddle Factors 
    [wR, wI] = TwiddleFactors(N);

    // Use decimation in time with mixed-split-RADIX algorithms to calculate DFT
    return DFT(Re, Im);

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    function FFT_Algorithm(REX, IMX, NFFT) {

        let NM1 = NFFT - 1;
        let ND2 = NFFT / 2;
        let M = parseInt( Math.log(NFFT) / Math.log(2) );
        let J = ND2;
        let I, L, TR, TI, K, UR, UI, LE, LE2, SR, SI, JM1, IP;
    
        // Start Bit reversal
        for (I = 1; I <= NFFT -2; I++) {
            if (I < J) {
                TR = REX[J];
                TI = IMX[J];
                REX[J] = REX[I];
                IMX[J] = IMX[I]
                REX[I] = TR;
                IMX[I] = TI;
            }
            K = ND2;
            while (K <= J){
                J = J - K;
                K = K / 2;
            }
            J = J + K;
        }
        for (L = 1; L <= M; L ++) {
            LE = parseInt(Math.pow(2, L));
            LE2 = LE / 2;
            UR = 1;
            UI = 0;
            SR = Math.cos(Math.PI/LE2);
            SI = -Math.sin(Math.PI/LE2);
    
            for (J = 1; J <= LE2; J++){
                JM1 = J - 1;
                for (I = JM1; I <= NM1; I+=LE) {
                    IP = I + LE2;
                    TR = REX[IP] * UR - IMX[IP] * UI;
                    TI = REX[IP] * UI + IMX[IP] * UR;
                    REX[IP] = REX[I] - TR;
                    IMX[IP] = IMX[I] - TI;
                    REX[I] = REX[I] + TR;
                    IMX[I] = IMX[I] + TI;
                }
                TR = UR;
                UR = TR * SR - UI * SI;
                UI = TR * SI + UI * SR;
            }
        }
    }

    function DFT(Re, Im) {
        
        // Return if length of Re[] array is equal to one (1)
        if (Re.length == 1) { return [Re, Im]; }

        if ((Re.length % 2) == 0) {
            // GO TO RADIX-2
            return RADIX2(Re, Im);
        }
        if ((Re.length % 3) == 0) {
            // GOTO RADIX-3
            return RADIX3(Re, Im)
        }
        if ((Re.length % 5) == 0) {
            // GOTO RADIX-5
            return RADIX5(Re, Im)
        }
        else {
            // GOTO COOLEY-TUKEY
            return DFT_Algorithm(Re, Im);
        }
    }

    function DFT_Algorithm(Re, Im) {

        let NRe = Re.length;
        let i, b=new Array(NRe);
        let L, M, facL=[], facM=[];

        // Index numbers
        for (i=0; i<NRe; i++) { b[i]=i; }

        // Prime-Factors
        [L, M, facL, facM] = PrimeFactors(Factors(NRe), NRe);
        
        // Compute DFT
        if (L==1 || M==1) {
            return Bluestein_Algorithm(Re, Im, b);
        }
        else {
            return CooleyTukey(Re, Im, b, L, M, facL, facM);
        }
    }

    function CooleyTukey(Re, Im, b, L, M, facL, facM) {

        // Variables 
        let i, j, Ind, Idx=0, Flag, L1, M1, st, N1=b.length;
		let bb=new Array(M), facLL=[], facMM=[], d=new Array(N1);
		let Xr=[], Xi=[];

        // Flag
        if (facM.length==1) {Flag=true;} else {Flag=false;}

        for (i=0; i<L; i++) {
			// Temp vectors
            bb=[]; Xr=[]; Xi=[];
			Ind=0; for (j=i; j<N1; j+=L) { bb[Ind]=b[j];   d[Idx]=b[j];   Idx++; Ind++; }
			
            // Compute the M-point DFT of each row
            if (Flag) { 
                // If M is bigger than 103, use Bluestein_Algorithm
                if (M >= 101) { [Xr, Xi] = Bluestein_Algorithm(Re, Im, bb); }
                else          { [Xr, Xi] = DTFT(Re, Im, bb);                }
            }
            else      { 
				[L1, M1, facLL, facMM] = PrimeFactors(facM, M);
				[Xr, Xi] = CooleyTukey(Re, Im, bb, L1, M1, facLL, facMM); 
			}
		 	
			// Multiply the resulting array by the phase factors W^(l*q)
			st = N/N1;
			
            for (j=0; j<M; j++) {
                Ind = (j * i * st);
                Re[bb[j]] = (wR[Ind] * Xr[j]) - (wI[Ind] * Xi[j]);
                Im[bb[j]] = (wR[Ind] * Xi[j]) + (wI[Ind] * Xr[j]);
            }
        }
        
        // clear Xr, Xi, bb ---------------------------------------------------------
        bb=new Array(L);
		
        // Flag
        if (facL.length==1) {Flag=true;} else {Flag=false;}
		
        for (i=0; i<M; i++) {
			// Temp vectors
			Ind=0; for (j=i*L; j<(i+1)*L; j++) { bb[Ind]=b[j]; Ind++; }
			
            // Compute the L-point DFT of each column
            if (Flag) { 
                // If L is bigger than 103, use Bluestein_Algorithm
                if (L >= 101) { [Xr, Xi] = Bluestein_Algorithm(Re, Im, bb); }
                else          { [Xr, Xi] = DTFT(Re, Im, bb);                }
            }
            else      { 
				[L1, M1, facLL, facMM] = PrimeFactors(facL, L);
				[Xr, Xi] = CooleyTukey(Re, Im, bb, L1, M1, facLL, facMM); 
			}
            
            for (j=0; j<L; j++) {
                Re[bb[j]] = Xr[j];
                Im[bb[j]] = Xi[j];
            }
        }
        // Read the resulting array row-wise and return
        return [d.map((v,i)=>Re[v]), d.map((v,i)=>Im[v])];
        //return [GetRange(Re, d), GetRange(Im, d)];
    }

    function Bluestein_Algorithm(Re, Im, bb) {

        // Discrete Time Fourier Transform for prime-length signals.
        // Matlab code of Bluestein Algorithm
        // x = [1,12,3,14,5,6,17,8,9,10,11,12,13];
        // N = length(x);
        // L = 2^ceil(log((2*N-1))/log(2));
        // n = 0:N-1;
        // W = exp(-j*pi*n.*n/N);
        // FW = fft([conj(W), zeros(1,L-2*N+1), conj(W(N:-1:2))], L);
        // y = ifft(FW .* fft(x.*W, L));
        // y = y(1:N) .* W;
        // "https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Fast_Fourier_Transforms_(Burrus)/04%3A_The_DFT_as_Convolution_or_Filtering/4.03%3A_The_Chirp_Z-Transform_or_Bluestein's_Algorithm"
    
        let NN = bb.length; 
        let wRR=[], wII=[], Mr=[], Mi=[], Fr=[], Fi=[];
        let a, b, L;

        // length of FFT - power of 2
        L  = Math.pow(2, Math.ceil(Math.log((2*NN-1))/Math.log(2)));

        // Twiddle Factors for Bluestein
        if (NN % 2 == 0) {
            // NN is even number 
            [wRR, wII] = TwiddleFactors_Bluestein(NN);
        }
        else {
            // NN is odd number 
            wRR = new Array(NN);
            wII = new Array(NN);
            wRR[0] = wR[0];
            wII[0] = wI[0];

            a = N / NN;

            for (i=1; i<NN; i+=2) {
                b = (((i*i)+NN)*a/2) % (N);
                wRR[i]    =  wR[b];
                wII[i]    =  wI[b];
                wRR[NN-i] = -wRR[i];
                wII[NN-i] = -wII[i];
            }
            for (i=2; i<NN; i+=2) {
                b = ((i*i)*a/2) % (N);
                wRR[i]    =  wR[b];
                wII[i]    =  wI[b];
                wRR[NN-i] = -wRR[i];
                wII[NN-i] = -wII[i];
            }
        }
        
        //
        Mr = new Array(L).fill(0);
        Mi = new Array(L).fill(0);
        for (i=0; i<NN; i++) {
            a = Re[bb[i]]*wRR[i] - Im[bb[i]]*wII[i];
            b = Re[bb[i]]*wII[i] + Im[bb[i]]*wRR[i];
            Mr[i] = a;
            Mi[i] = b;
        }
        
        // Zero padding
        Fr = new Array(L).fill(0);
        Fi = new Array(L).fill(0);
        Fr[0] = wRR[0];
        Fi[0] = wII[0];
        for (i=1; i<NN; i++) {
            Fr[i] = wRR[i];
            Fi[i] = -wII[i];
            Fr[L-i] = Fr[i];
            Fi[L-i] = Fi[i];
        }
        
        // FFT
        FFT_Algorithm(Fr, Fi, L);
        FFT_Algorithm(Mr, Mi, L);
    
        for (i=0; i<L; i++) {
            a = Fr[i]*Mr[i] - Fi[i]*Mi[i];
            b = Fr[i]*Mi[i] + Fi[i]*Mr[i];
            Fr[i] = a;
            Fi[i] = b;
        }
        
        // Inverse IFT
        [Fr, Fi] = IFFT(Fr, Fi, L);
        
        for (i=0; i<NN; i++) {
            a = Fr[i]*wRR[i] - Fi[i]*wII[i];
            b = Fr[i]*wII[i] + Fi[i]*wRR[i];
            wRR[i] = a;
            wII[i] = b;
        }
        return [wRR,  wII]
    }

    // RADIX-2 FFT 
    function RADIX2(Re, Im) {

        let NN = Re.length;
        let k, idx, Ind1, Ind2;
        let tt = NN/2;

        let Re1 = new Array(tt),   Im1 = new Array(tt);
        let Re2 = new Array(tt),   Im2 = new Array(tt);
        
        // Decimating the time domain sequence
        for (k=0; k<tt; k++) {
            Ind1 = 2*k;
            Ind2 = Ind1 + 1;

            Re1[k] = Re[Ind1];    Im1[k] = Im[Ind1];
            Re2[k] = Re[Ind2];    Im2[k] = Im[Ind2];
        }
        
        // Performing DFT on each sub-sequence recursively
        [Re1, Im1] = DFT(Re1, Im1);
        [Re2, Im2] = DFT(Re2, Im2);

        Ind1 = N / NN;

        // Combining the separately calculated DFTs
        for (k=0; k<NN; k++) {
            idx = k % (tt);
            Ind2 = k * Ind1;

            Re[k] = Re1[idx] + (Re2[idx]*wR[Ind2]) - (Im2[idx]*wI[Ind2]);
            Im[k] = Im1[idx] + (Re2[idx]*wI[Ind2]) + (Im2[idx]*wR[Ind2]);
        }
        return [Re, Im]
    }

    function RADIX3(Re, Im) {

        let NN = Re.length;
        let k, idx, Ind1, Ind2, Ind3;
        let tt = NN/3;

        let Re1 = new Array(tt),   Im1 = new Array(tt);
        let Re2 = new Array(tt),   Im2 = new Array(tt);
        let Re3 = new Array(tt),   Im3 = new Array(tt);
        
        // Decimating the time domain sequence
        for (k=0; k<tt; k++) {
            Ind1 = 3*k;
            Ind2 = Ind1 + 1;
            Ind3 = Ind2 + 1;

            Re1[k] = Re[Ind1];    Im1[k] = Im[Ind1];
            Re2[k] = Re[Ind2];    Im2[k] = Im[Ind2];
            Re3[k] = Re[Ind3];    Im3[k] = Im[Ind3];
        }
        
        // Performing DFT on each sub-sequence recursively
        [Re1, Im1] = DFT(Re1, Im1);
        [Re2, Im2] = DFT(Re2, Im2);
        [Re3, Im3] = DFT(Re3, Im3);
        
        Ind1 = N / NN;

        // Combining the separately calculated DFTs
        for (k=0; k<NN; k++) {
            idx  = k % (tt);
            Ind2 = k * Ind1;
            Ind3 = (2 * Ind2) % (N);

            Re[k] = Re1[idx] + (Re2[idx]*wR[Ind2]) - (Im2[idx]*wI[Ind2]) + (Re3[idx]*wR[Ind3]) - (Im3[idx]*wI[Ind3]);
            Im[k] = Im1[idx] + (Re2[idx]*wI[Ind2]) + (Im2[idx]*wR[Ind2]) + (Re3[idx]*wI[Ind3]) + (Im3[idx]*wR[Ind3]);
        }
        return [Re, Im]
    }

    // RADIX-5 FFT 
    function RADIX5(Re, Im) {

        let NN = Re.length;
        let k, idx, Ind1, Ind2, Ind3, Ind4, Ind5;
        let tt = NN/5;
        
        let Re1 = new Array(tt),   Im1 = new Array(tt);
        let Re2 = new Array(tt),   Im2 = new Array(tt);
        let Re3 = new Array(tt),   Im3 = new Array(tt);
        let Re4 = new Array(tt),   Im4 = new Array(tt);
        let Re5 = new Array(tt),   Im5 = new Array(tt);
        
        // Decimating the time domain sequence
        for (k=0; k<tt; k++) {
            Ind1 = 5*k;
            Ind2 = Ind1 + 1;
            Ind3 = Ind2 + 1;
            Ind4 = Ind3 + 1;
            Ind5 = Ind4 + 1;

            Re1[k] = Re[Ind1];    Im1[k] = Im[Ind1];
            Re2[k] = Re[Ind2];    Im2[k] = Im[Ind2];
            Re3[k] = Re[Ind3];    Im3[k] = Im[Ind3];
            Re4[k] = Re[Ind4];    Im4[k] = Im[Ind4];
            Re5[k] = Re[Ind5];    Im5[k] = Im[Ind5];
        }
        
        // Performing DFT on each sub-sequence recursively
        [Re1, Im1] = DFT(Re1, Im1);
        [Re2, Im2] = DFT(Re2, Im2);
        [Re3, Im3] = DFT(Re3, Im3);
        [Re4, Im4] = DFT(Re4, Im4);
        [Re5, Im5] = DFT(Re5, Im5);
        
        Ind1 = N / NN;

        // Combining the separately calculated DFTs
        for (k=0; k<NN; k++) {
            idx = k % (tt);
            Ind2 = k * Ind1;
            Ind3 = (2 * Ind2) % (N);
            Ind4 = (3 * Ind2) % (N);
            Ind5 = (4 * Ind2) % (N);

            Re[k] = Re1[idx] + (Re2[idx]*wR[Ind2]) - (Im2[idx]*wI[Ind2]) + (Re3[idx]*wR[Ind3]) - (Im3[idx]*wI[Ind3]) + (Re4[idx]*wR[Ind4]) - (Im4[idx]*wI[Ind4]) + (Re5[idx]*wR[Ind5]) - (Im5[idx]*wI[Ind5]);
            Im[k] = Im1[idx] + (Re2[idx]*wI[Ind2]) + (Im2[idx]*wR[Ind2]) + (Re3[idx]*wI[Ind3]) + (Im3[idx]*wR[Ind3]) + (Re4[idx]*wI[Ind4]) + (Im4[idx]*wR[Ind4]) + (Re5[idx]*wI[Ind5]) + (Im5[idx]*wR[Ind5]);

        }
        return [Re, Im]
    }

    // Discrete Time Fourier Transform
    function DTFT(Re, Im, b) {
        let N4=b.length;
        let Indx, k, n, rr, ii, Ind, Ind2=N4-1;
        let Xrr=new Array(N4), Xii=new Array(N4);
        let st = N/N4;

        for (k=0; k<N4; k++) {
            // for n=0
            rr   = 0;
            ii   = 0;
            Indx = 0;
    
            for (n=0; n<N4; n++) {
                Ind = (Indx * st);
                rr += (Re[b[n]] * wR[Ind]) - (Im[b[n]] * wI[Ind]);
                ii += (Re[b[n]] * wI[Ind]) + (Im[b[n]] * wR[Ind]);
                Indx += k;
                Indx = Indx % N4;
                //if (Indx > Ind2) { Indx -= N4; }
            }
            Xrr[k] = rr;
            Xii[k] = ii;
        }
        return [Xrr, Xii]
    }

    function PrimeFactors(fac, N) {
        
        // Decleration of variables
        let i, a=Math.sqrt(N), res=[];
        let r=1, r1, r2, r3, r4, r5, r6;
        let f1=[], f2=Copy(fac);

        // Prime Factors
        for (i=0; i<fac.length; i++) {
            r *= fac[i];
            f1.push(fac[i]); // add the end of f1[] array
            f2.shift();      // remove the first element from f2[] array
            if (r>=a) { 
                r1 = r/fac[i];   r2 = N/r1;   r3 = r1+r2;
                r4 = r;          r5 = N/r4;   r6 = r4+r5; 
                if (r3<r6){
                    f1.pop();           // remove last element from f1[] array
                    f2.unshift(fac[i])  // add to the first element of f2[] array
                    if (r2>r1) {return [r2, r1, f2, f1];} else {return [r1, r2, f1, f1];}
                }
                else {
                    if (r5>r4) {return [r5, r4, f2, f1];} else {return [r4, r5, f1, f2];}
                }
            }
        }
        return res;
    }
        
    function TwiddleFactors(N) {

        // Decleration of variables
        let a, b, c, i, Ind1, Type;
    
        // Pre-allocation of 2D array of w[]
        let wR = new Array(N);
        let wI = new Array(N);
    
        a = 2 * Math.PI / N;
    
        if (N % 2 == 0) {
            // N is Even
            Ind1 = (N / 2) - 1;
            Type = 1;
        } else {
            // N is Odd
            Ind1 = (N - 1) / 2;
            Type = 0;
        }
        wR[0] = 1;
        wI[0] = 0;
        for (i = 1; i <= Ind1; i++) {
            b = Math.cos(a*i);
            c = Math.sin(a*i);
            wR[i] = b;
            wI[i] = -c;
            wR[N-i] = b;
            wI[N-i] = c;
        }
        if (Type == 1) {
            wR[Ind1+1] = Math.cos(a*(Ind1+1));
            wI[Ind1+1] = -Math.sin(a*(Ind1+1));
        };
        return [wR, wI]; 
    } 

    function TwiddleFactors_Bluestein(N) {

        // Decleration of variables
        let a, b, c, d, i, Ind1;
    
        // Pre-allocation of 2D array of w[]
        let wR = new Array(N);
        let wI = new Array(N);
        wR[0] = 1;
        wI[0] = 0;
    
        a = Math.PI / N;

        if (N % 2 == 0) {
            // N is Even
            Ind1 = (N / 2) - 1;
            for (i=1; i<=Ind1; i++) {
                d = a*i*i;
                b = Math.cos(d);
                c = Math.sin(d);
                wR[i] = b;
                wI[i] = -c;
                wR[N-i] = b;
                wI[N-i] = -c;
            }
            wR[Ind1+1] = Math.cos(a*(Ind1+1)*(Ind1+1));
            wI[Ind1+1] = -Math.sin(a*(Ind1+1)*(Ind1+1));
        } else {
            // N is Odd
            Ind1 = (N - 1) / 2;
            for (i=1; i<=Ind1; i++) {
                d = a*i*i;
                b = Math.cos(d);
                c = Math.sin(d);
                wR[i] = b;
                wI[i] = -c;
                wR[N-i] = -b;
                wI[N-i] = c;
            }
        }
        return [wR, wI]; 
    }
}

function GCD(a, b) {
    // Convert 'a' and 'b' to numbers.
    a = +a;
    b = +b;
    // Check if 'a' or 'b' is NaN (Not a Number).
    if (a !== a || b !== b) {
        return [NaN, NaN, NaN];
    }

    // Check if 'a' or 'b' is Infinity or -Infinity.
    if (a === Infinity || a === -Infinity || b === Infinity || b === -Infinity) {
        return [Infinity, Infinity, Infinity];
    }

    // Check if 'a' or 'b' are decimals.
    if ((a % 1 !== 0) || (b % 1 !== 0)) {
        return false;
    }

    // Initialize variables for signs and quotient.
    var signX = (a < 0) ? -1 : 1
    let signY = (b < 0) ? -1 : 1
    let x=0, y=1, u=1, v=0, q, r, m, n;

    // Get the absolute values of 'a' and 'b'.
    a = Math.abs(a);
    b = Math.abs(b);

    // Implement Euclid's algorithm to find the GCD.
    while (a !== 0) {
        q = Math.floor(b / a);
        r = b % a;
        m = x - u * q;
        n = y - v * q;
        b = a;
        a = r;
        x = u;
        y = v;
        u = m;
        v = n;
    }

    // Return an array containing the GCD, and the coefficients for 'a' and 'b'.
    return [b, signX * x, signY * y];
}

function Factorization(N) {

    // N must be an integer number 
    if (!Number.isInteger(N)) { return []; }

    // Return if N==1
    if (N==1) {return [1, 1];}

    // Factors of N
    let fac = Factors(N);

    let result_1 = PrimeFactors(fac, N);
    let result_2 = RelativelyPrimeFactors(fac, N);

    return ChoseMethod(result_1, result_2, fac, N);



    function ChoseMethod(result_1, result_2, fac, N) {

        // Decleration of variables
        let sum1=Sum(result_1), sum2=Sum(result_2);

        // Chose Good-Thomas is sum1==sum2
        if (sum1 == sum2) { 
            return { 
                "CRT"       : CRT(result_2, N), 
                "LM"        : result_2,  
                "Method"    : "GTM",
            }; 
        }

        // Determine the lowest sum
        Ind = Min([sum1, sum2]).Indx;

        if (Ind == 0) {
            // Cooley-Tukey is the most efficient method
            return { 
                "LM"        : result_1,  
                "Method"    : "CTM",
            };
        }
        else {
            // Good-Thomas is the most efficient method
            return { 
                "CRT"       : CRT(result_2, N), 
                "LM"        : result_2,  
                "Method"    : "GTM",
            };
        }
    }

    function PrimeFactors(fac, N) {
        
        // Decleration of variables
        let i, a=Math.sqrt(N), res=[];
        let r=1, r1, r2, r3, r4, r5, r6;

        // Prime Factors
        for (i=0; i<fac.length; i++) {
            r *= fac[i];
            if (r>=a) { 
                r1 = r/fac[i];   r2 = N/r1;   r3 = r1+r2;
                r4 = r;          r5 = N/r4;   r6 = r4+r5; 
                if (r3<r6){
                    if (r2>r1) {res=[r2, r1];} else {res=[r1, r2];}
                }
                else {
                    if (r5>r4) {res=[r5, r4];} else {res=[r4, r5];}
                }
                return res;
            }
        }
        return res;
    }

    function RelativelyPrimeFactors(fac, N) {

        // Decleration of variables
        let UFac = [...new Set(fac)];
        let M = UFac.length;
        let res = [];
        let a=Math.sqrt(N);
        
        // Relatively PrimeFactors
        for (i=0; i<M; i++) {
            res.push(Pow(UFac[i], CountUnique(fac, UFac[i])));
        }
        if (res.length==1) {res.push(1);};
        if (res.length > 2) { res = PrimeFactors(res, N); }
        return res;
    }

    function CountUnique(a, b) {
        c=0; 
        a.map((v,i) => {if (v==b) {c+=1;}});
        return c;
    }

    function NumOfOperation(N1, N2, Type) {

        // Decleration of variables
        let t1, t2;

        if (Type == 0) {
            // Cooley-Tukey Method

            // N1 times N2-point DFT
            t1 = Mult(DFT_NumOperation(N1), N2);

            // (N2 times N1-point DFT) + (4*N2 real-number multiplication) + (2*N2 real-number additions)
            t2 = Mult(Add(DFT_NumOperation(N2),  [4*N2, 2*N2]), N1);

            return [t1, t2];
        }
        else if (Type == 1) {
            // Good-Thomas Method
            // N1 times N2-point DFT
            t1 = Mult(DFT_NumOperation(N1), N2);
            
            // (N2 times N1-point DFT)
            t2 = Mult(DFT_NumOperation(N2), N1);

            return [t1, t2];
        }
        else {
            // FFT
            return (N1*N2) * Math.log2(N1*N2);
        }
    }

    function DFT_NumOperation(N) {
        // return total number of real-number multiplications and rel-number additions
        return [4*N*N, 2*N*(2*N-1)];
    }

    function CRT(m, N) {
        // Decleration of varables 
        let i, Mi, g, ni, Ni, c=[];
    
        for (i=0; i<m.length; i++) {
            Mi          = N / m[i];
            [g, ni, Ni] = GCD(m[i], Mi);
            c[i]        = Mi * Ni;
        }
        return c;
    }
}

function Factors(n) {

    if (!Number.isInteger(n)) { return []; }

    let factors=[], divisor = 2;
    
    while (n >= 2) {
        if (n % divisor == 0) {
        factors.push(divisor);
        n = n / divisor;
        } else {
        divisor++;
        }
    }
    return factors;
}

function IFFT(REX, IMX, NFFT) {

    // Decleration of variables 
    let i;

    // Check default values of input arguments
    if (IMX  == null)  { IMX = new Array(REX.length).fill(0); };
    if (NFFT == null)  { NFFT = REX.length; };

    // Decleration of variables
    let Re=[], Im = [];

    // Reverse the Real and Imaginary part and take FFT
    [Re, Im] = FFT(IMX, REX, NFFT);

    // Divide the results by NFFT and revers Real and Imaginary arrays
    for (i=0; i<NFFT; i++) { Re[i] /= NFFT;  Im[i] /= NFFT;   }

    return [Im, Re];
}

function FFT_MagPhase(REX, IMX, NFFT) {
    // Check default values of input arguments
    if (NFFT == null) { NFFT = REX.length; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if (!Array.isArray(REX))      {throw new Error("REX[] must be an array.");};
    if (!Array.isArray(IMX))      {throw new Error("IMX[] must be an array.");};
    if (REX.length != IMX.length) {throw new Error("REX[] and IMX[] must be the same length.");};

    var Mag =  new Array(NFFT).fill(0);
    var Phase =  new Array(NFFT).fill(0);

    for (var i = 0; i < NFFT; i++) {
       Mag[i]   = Math.sqrt(REX[i] * REX[i] + IMX[i] * IMX[i]);
       Phase[i] = Math.atan(IMX[i] / REX[i]);
    }

    return [Mag, Phase];
}

function FFT_Angle(REX, IMX, NFFT) {
    // Check default values of input arguments
    if (NFFT == null) { NFFT = REX.length; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if (!Array.isArray(REX))      {throw new Error("REX[] must be an array.");};
    if (!Array.isArray(IMX))      {throw new Error("IMX[] must be an array.");};
    if (REX.length != IMX.length) {throw new Error("REX[] and IMX[] must be the same length.");};

    var Angle =  new Array(NFFT).fill(0);

    // Return four quadrant arc-tangent of a complex number REX[] and IMX[]
    for (var i = 0; i < NFFT; i++) {
        if (REX[i] > 0) {
            Angle[i] = Math.atan(IMX[i] / REX[i]);
        }
        else if ((REX[i] < 0) && (IMX[i] >= 0)) {
            Angle[i] = Math.atan(IMX[i] / REX[i]) + Math.PI;
        }
        else if ((REX[i] < 0) && (IMX[i] < 0)) {
            Angle[i] = Math.atan(IMX[i] / REX[i]) - Math.PI;
        }
        else if ((REX[i] == 0) && (IMX[i] > 0)) {
            Angle[i] = Math.PI / 2;
        }
        else if ((REX[i] == 0) && (IMX[i] < 0)) {
            Angle[i] = -Math.PI / 2;
        }
        else if ((REX[i] == 0) && (IMX[i] == 0)) {
            // this is actually undefined
            Angle[i] = 0;
        }
    }
    return Angle;
}

//----------------------------------------------------------------------------------------------------------------------
// Spectrum Estimation

function Hilbert(f) {
    // Hilbert transform is a mathematical process performed on a rel signal of f(t), yielding a new real signal H(t)
    // that is a 90-degree phase-shifted version of f(t).
    // This function calculates the Hilbert transform of a real-valued f[] array.
    // Hilbert transform H[t] adds a phase shift of ±90° (π⁄2 radians) to every frequency component of f[] array.
    // The sign of the phase shift depends on the sign of the frequency (positive or negative).
    // The Hilbert transform is a component of the analytic representation of a real-valued signal, u(t)
    //
    //     u(t) = f(t) + (1i) * H(t) =  A(t) * exp(i * P(t))
    //
    //  Instantaneous amplitude     : A(t) = Sqrt( f(t) * f(t)  +  H(t) * H(t) )   (Also, envelope function)
    //  Instantaneous phase         : P(t) = atan(H(t)) / f(t))
    //  Instantaneous frequency     : w(t) = d( P(t) ) / dt       (time derivative of phase)
    //  Instantaneous energy        : e(t) = abs(A(t)) ^2
    //
    // Impulse response of Hilbert transform is h(n)
    //  h(n) = FSamp * (1- cos(pi*n)) / (pi*n)  using L'Hopital's Rule, for n=0, h(0)=0
    //  n     : discrete time-domain integer index (...,-3,-2,-1,0,1,2,3...)
    //  FSamp : sample rate measured in samples/second
    //
    // Some properties of Hilbert Transform:
    //
    //  1- Hilbert transform of a constant signal is equal to zero:
    //     f(t)  <==>  H(t)
    //     f(t) = [c, c, c, ... d], then  H( f(t) ) = 0
    //
    //  2- Time shifting and time-dilation:
    //     f(t)        <==>  H(t)
    //     f(t - to)   <==>  H(t- to)
    //     f(at)       <==>  sign(a) H(at)
    //
    //  3- The Hilbert transform of the derivative of a signal is the derivative of the Hilbert transform
    //     f(t)  <==>  H(t)
    //     H( df(t)/dt ) = dH(t)/ dt
    //
    //  4- A signal f(t) and its Hilbert transform H(t) are orthogonal to each other.
    //     f(t)  <==>  H(t)
    //     Integral( f(t) x H(t) ) = 0
    //
    //  5- Hilbert transform of a hilbert transform equals to negative of the original signal.
    //     f(t)         <==>  H(t)
    //     H( H(t) )  =  (-1) * f(t)
    //
    //  6- Signal f(t) and its hilbert transform H(t) has the same auto-correlation function
    //     f(t)         <==>  H(t)
    //     Xcorr( f(t) )  =  Xcorr( H(f) )
    //
    //  7- Magnitude spectrum of f(t) and H(t) are same
    //     f(t)       <==>  H(t)
    //     abs( F(w) ) =  abs( H(w) )
    //
    //  Output of this function is a complex values signal u(t) = f(t) + (1i) * H(t)
    //  Real values equals the f(t) array itself.
    //  Imaginary values equals the Hilbert transform of the f(t) array.
    //

    // Check if the input arguments are correct; otherwise, throw an error.
    if (!Array.isArray(data)) {throw new Error("data[] must be an array.");};

    // Length of data[] array
    n = f.length;

    // determine if the n is even or odd number
    if (IsEven(n)) { m = (n / 2) + 1; } else { m = (n + 1) / 2; }

    // Fourier transform of data[] array
    [Re, Im] = FFT(f);

    // Pre-allocation
    h = new Array(n).fill(0);
    for (i = 1; i < m - 1; i++) {
        Re[i] = 2 * Re[i];
        Im[i] = 2 * Im[i];
    }
    for (i = m; i < n; i++) {
        Re[i] = 0;
        Im[i] = 0;
    }

    // Inverse FFT
    [Re, Im] = IFFT(Re, Im);

    return [Re, Im];
}

function FourierSpec(data, FSamp, NFFT, OneSided) {

    // This function calculates the Fourier Amplitude Spectrum of the data[] array.
    // Fourier Magnitude Spectrum is normalized by the length of the data[] array
    // If OnSided is true, then half of the spectrum is returned.
    //
    //  by Dr. Yavuz Kaya, P.Eng.;
    //  Created on      22.Jul.2023
    //  Last modified   22.Jul.2023
    //  Yavuz.Kaya.ca@gmail.com
    //

    // Decleration of Variables
    let n, df, NumSamp, f=[], REX=[], IMX=[], Mag=[], Angle=[];

    // Check default values of input arguments
    if (FSamp    == null)  { FSamp    = 1;    };
    if (NFFT     == null)  { NFFT = data.length; };
    if (OneSided == null)  { OneSided = true; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if ((OneSided != true) && (OneSided != false)) {throw new Error("OneSided argument must be either true or false."); }
    if (FSamp <= 0)           {throw new Error("FSamp (Sampling frequency) cannot be equal or less than zero."); }
    if (!Array.isArray(data)) {throw new Error("data[] must be an array.");};

    // Calculate Fourier Amplitude Spectrum
    n = NFFT;
    [REX, IMX] = FFT(data, null, n);
    [Mag, ]    = FFT_MagPhase(REX, IMX);
    Angle      = FFT_Angle(REX, IMX);

    // Normalize the Fourier Magnitude spectrum
    Mag = Mag.map((v, i) => v / n);

    // Frequency vector
    df = FSamp / n;
    f  = LinSpace(0, (n-1)*df, n);

    if (OneSided) {
        // Up to the Nyquist Frequency

        // Double the amplitude of the Fourier Magnitude spectrum
        Mag = Mag.map((v, i) => 2 * v);

        if (IsEven(n)) {
            NumSamp   = (n / 2) + 1;
            Mag[0]   /= 2;
            Mag[n/2] /= 2;
        }
        else {
            NumSamp = (n + 1) / 2;
            Mag[0] /= 2;
        }

        Mag   = ExtractSubset(Mag,   0, NumSamp-1);
        Angle = ExtractSubset(Angle, 0, NumSamp-1);
        f     = ExtractSubset(f,     0, NumSamp-1);
    }
    return [Mag, Angle, f]
}

function FFT_Average(data, RWL, OVS, FSamp, NFFT, WinOption, OneSided) {

    //  This function returns average Fourier Amplitude Spectrum (FAS) and Fourier Phase Spectrum (FPS) of data[] array.
    //  The data[] array is divided into multiple equal-segments of RWL-points and overlapped by OVS sample.
    //  If Option is True, then each data-segment is multiplied by Hamming window of RWL-points before taking
    //  the NFFT-length of FFT. FFTs of all segments are averaged to form the FAS and FPS for the entire data[] array.
    //
    //  Input Arguments:
    //     data   : 1D array for which FAS and FPS is calculated
    //     RWL    : Number of samples in each data-segment
    //     OVS    : Number of samples that two data segments are overlapped
    //     NFFT   : Number of discrete Fourier Transform points in FFT estimate for each data segment
    //     FSamp  : Sampling frequency
    //  WinOption : TRUE/FALSE (True: Hamming window,  False: No Windowing )
    //              Each data-segment will be multiplied by a window function before FFT.
    //  OneSided  : 1-sided or 2-sided spectrum
    //              0 for one-sided  :::  1 for two-sided
    //
    //  Output:
    //      FAS  : Averaged Fourier Amplitude Spectrum
    //      FPS  : Averaged Fourier Phase Spectrum
    //        f  : Frequency vector over which the FAS and FPS are calculated
    //
    //  by Dr. Yavuz Kaya, P.Eng.;
    //  Created on      17.Mar.2017
    //  Last modified   25.Jun.2023
    //  Yavuz.Kaya.ca@gmail.com

    // Declare variables
    var Real, Win, DW, K, df, f, a1, a2, kk, SF, Mag, FAS, FS1, FS2, FF;

    // Check default values of input arguments
    if (RWL       == null) { RWL       = Math.max(2, Math.floor(row / 8.00)); };
    if (OVS       == null) { OVS       = Math.floor(RWL * 0.25); };
    if (FSamp     == null) { FSamp     = 1; };
    if (NFFT      == null) { NFFT      = RWL; };
    if (WinOption == null) { WinOption = true; };
    if (OneSided  == null) { OneSided  = false; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if ((OneSided != true) && (OneSided != false)) {throw new Error("OneSided argument must be either true or false."); }
    if ((RWL <= 0) || (RWL > data.length) || (RWL > NFFT) || (OVS >= RWL) || (OVS < 0) || (NFFT < 2) || (FSamp <= 0)) {
        {throw new Error("Input arguments are inconsistent (RWL <= 0) || (RWL > data.Length) || (RWL > NFFT) || (OVS >= RWL) || (OVS < 0) || (NFFT < 2) || (FSamp <= 0)");};
    }
    if (!Array.isArray(data))  {throw new Error("data[] must be 1D array.");};

    // Pre-allocation of variables
    FAS = new Array(NFFT).fill(0);
    FS1 = new Array(NFFT).fill(0);
    FS2 = new Array(NFFT).fill(0);

    // Calculate initial parameters
    DW = RWL - OVS;

    // Create a window function of RWL-point and normalize it so that the power is 1 Watt.
    if (WinOption) { Win = Hamming(RWL); } else { Win = Rectwin(RWL); }
    Win = Vector_Multiply(Win, Math.sqrt(RWL / Statistics_Sum(Vector_Pow(Win, 2))));

    // Icg  = sum(Win.^ 2) / RWL;    % Incoherent Power Gain of the window Function
    // Cg   = (sum(Win) / RWL) ^ 2;  % Coherent Power Gain of the window Function
    // ENBW = Icg / Cg;              % Effective Noise Bandwidth of the Window

    // Assuming NFFT >= RWL, which is already satisfied above
    // Scaling is done as follow:
    // Each segment is divided by the number of bin(NFFT)
    // Scale each segment by a factor of(NFFT/ RWL) in case NFFT >= RWL
    // Window correction factor(RWL/ sum(Win)) is applied to each segment
    SF = (1 / NFFT) * (NFFT / RWL) * (RWL / Statistics_Sum(Win));

    K = 0;
    a1 = 0;
    a2 = RWL - 1;

    while (a2 < data.length) {
        // Extract the segment of the data[] array starting from Index a1 to Index a2
        Real = Vector_GetRange(data, a1, a2);

        // Skip the this segment if data[] array contains NaN
        if (IsNaN(Real)) { a1 += DW; a2 += DW; continue; }

        // Fourier Amplitude Spectrum (FAS) of the windowed-segment
        [Re, Im] = FFT(Vector_Multiply(Real, Win), null, NFFT);

        // Sum the Fourier Magnitude response of each segment
        [Mag, ] = FFT_MagPhase(Re, Im);
        FAS     = FAS.map((v,i) => v + Mag[i])

        // Fourier Amplitude Spectrum (FAS) of the windowed-segment
        [Re, Im] = FFT(Real, null, NFFT);

        // Continue adding the Real and Imaginary
        FS1 = FS1.map((v,i) => v + Re[i]);
        FS2 = FS2.map((v,i) => v + Im[i]);

        // Update the Indexes
        a1 += DW;
        a2 += DW;
        K++;
    }
    // Normalization of the averaged FFT Spectrum
    FAS = Vector_Divide(FAS, K / SF);
    FS1 = FFT_Angle(Vector_Divide(FS1, K), Vector_Divide(FS2, K));

    // Frequency vector
    df = FSamp / NFFT;
    f  = Vector_LinSpace(0, (NFFT-1)*df, NFFT);

    return [FAS, FS1, f];
}

function PowerSpectralDensity(data, RWL, OVS, FSamp, NFFT, WinOption, OneSided) {
    // This function calculates Power Spectral Density (PSD) for 1D data[] array.
    //
    // The periodogram is not a consistent estimator of the true power spectral density of a wide-sense stationary
    // process. Welch’s technique to reduce the variance of the periodogram breaks the time series into segments,
    // usually overlapping.
    //
    // Welch’s method computes a modified periodogram (power-spectrum) for each segment and then averages these
    // estimates to produce the estimate of the power spectral density. Because the process is wide-sense stationary
    // and Welch’s method uses PSD estimates of different segments of the time series, the modified periodograms
    // represent approximately uncorrelated estimates of the true PSD and averaging reduces the variability.
    //
    // The segments are typically multiplied by a window function, such as a Hamming window, so that Welch’s method
    // amounts to averaging modified periodograms. Because the segments usually overlap, data values at the beginning
    // and end of the segment tapered by the window in one segment, occur away from the ends of adjacent segments.
    // This guards against the loss of information caused by windowing.
    //
    // This method (Welch's method) relies on the averaging the periodograms of each segment to reduce the variance of
    // the Power Spectrum estimator.  K-independent average reduces the variance by a factor of K as long as the
    // averaged segments are independent. Resolution of the PSD increases as the length of the window increases.
    //
    // The power spectral density (PSD) indicates the strength of the energy as a function of frequency.
    // In other words, it shows at which frequencies energy variations are strong and at which frequencies they are weak.
    // The unit of PSD is energy per bandwidth and energy can be obtained within a specific frequency range by
    // integrating PSD within that frequency range.
    //
    // If data[] array is in units of "g", then the PSD is in units of (g•s)²/Hz
    //
    // The integral of the PSD over a frequency interval gives the mean-square value of the content of the signal
    // within the frequency interval. (This is known as Parseval’s Theorem.)
    //
    //  Input:
    //      Data      : 1D array of data.
    //      RWL       : Number of samples in each data-segment
    //      OVS       : Number of samples that two data segments are overlapped
    //      FSamp     : Sampling frequency in Hz.
    //      NFFT      : Number of discrete Fourier Transform points in FFT estimate for each data segment
    //      WinOption : true/false (true: Hamming window,  false: No Windowing)
    //                  Each data-segment will be multiplied by a window function before FFT.
    //                  WinOption == true will run Pwelch algorithm with Hamming window,
    //                  WinOption == false will run Pbartlett algorithm with Rectangular window
    //      OneSided  : 1-sided or 2-sided spectrum
    //                  0 returns one-sided spectrum :::  1 returns two-sided spectrum
    //
    //  Output:
    //      PSD     : Averaged Power Spectral Density Spectrum
    //        fq    : Frequency vector over which the power spectral density is calculated
    //
    //  by Dr. Yavuz Kaya, P.Eng.;
    //  Created on      18.Mar.2017
    //  Last modified   03.Feb.2024
    //  Yavuz.Kaya.ca@gmail.com
    //

    // Check default values of input arguments
    if (RWL       == null) { RWL       = Math.max(2, Math.floor(data.length / 8.00)); };
    if (OVS       == null) { OVS       = Math.floor(RWL * 0.25); };
    if (FSamp     == null) { FSamp     = 1; };
    if (NFFT      == null) { NFFT      = NextPow2(RWL); };
    if (WinOption == null) { WinOption = true; };
    if (OneSided  == null) { OneSided  = true; };

    // Check if the input arguments are correct; otherwise, throw an error.
    if ((OneSided != true) && (OneSided != false)) {throw new Error("OneSided argument must be either true or false."); }
    if ((RWL <= 0) || (RWL > data.length) || (RWL > NFFT) || (OVS >= RWL) || (OVS < 0) || (NFFT < 2) || (FSamp <= 0)) {
        {throw new Error("Input arguments are inconsistent (RWL <= 0) || (RWL > data.Length) || (RWL > NFFT) || (OVS >= RWL) || (OVS < 0) || (NFFT < 2) || (FSamp <= 0)");};
    }
    if (!Array.isArray(data))  {throw new Error("data[] must be 1D array.");};

    //
    let PSD, DW, Win, SF, K, a1, a2, Re, Im, df, fq;

    // Pre-allocate
    PSD = new Array(NFFT).fill(0);

    // Calculate initial parameters
    DW = RWL - OVS;

    // Create a window function of RWL-point and normalize it
    if (WinOption) { Win = Hamming(RWL); } else { Win = Rectwin(RWL); }
    Win = Vector_Divide(Win, Statistics_Sum(Win));

    // Scale Factor for the PSD normalization
    // Units of PSD is "Hertz" because it is scaled by the sampling frequency to obtain the PSD
    SF = Statistics_Sum(Vector_Pow(Win, 2)) * FSamp

    K = 0;
    a1 = 0;
    a2 = RWL - 1;

    while (a2 < data.length) {
        // Extract the segment of the data[] array starting from Index a1 to Index a2
        Re = Vector_GetRange(data, a1, a2);

        // Skip one iteration in the for-loop if data[] array contains NaN
        if (IsNaN(Re)) { a1 += DW; a2 += DW; continue; }

        // Fourier Amplitude Spectrum (FAS) of the windowed-segment
        [Re, Im] = FFT(Vector_Multiply(Re, Win), null, NFFT);

        // Magnitude and Phase of the FAS
        [Re, _] = FFT_MagPhase(Re, Im);

        // Continue adding the Power Spectrum
        PSD = PSD.map((v, i) => v + Re[i] * Re[i] );

        // Update the Indexes
        a1 += DW;
        a2 += DW;
        K++;
    }

    // Normalization of the averaged Power Spectrum to calculate PSD
    // Units of PSD is Hertz because it is scaled by the sampling frequency to obtain the psd
    PSD = Vector_Divide(PSD, (SF * K));

    // Frequency vector
    df = FSamp / NFFT;
    fq = Vector_LinSpace(0, (NFFT-1)*df, NFFT);

    if (OneSided) {
        // Return one-sided spectrum
        PSD = PSD.map((v, i) => v * 2);

        if (NFFT % 2 == 0) {
            // NFFT is even
            PSD[0]      /= 2;
            PSD[NFFT/2] /= 2;
            return [Vector_Truncate(PSD, NFFT/2+1), Vector_Truncate(fq, NFFT/2+1)];
        }
        else {
            // NFFT is odd
            PSD[0]          /= 2;
            return [Vector_Truncate(PSD, (NFFT+1)/2), Vector_Truncate(fq, (NFFT+1)/2)];
        }
    }
    return [PSD, fq];
}

function TransferFunction(Input, Output, RWL, FSamp, OLR) {

    // This function calculates the Transfer Function between 1D Input[] array and and 2D Output[] matrix.
    // Sampling frequency for Input[] and Output[] must be the same.
    //
    // Transfer function is calculated for each column of the 2D Output[] matrix.
    // The raw data in Input and Output matrices is divided into equal length of segments of RWL-point.
    // Each segment is overlapped by OLR ratio. Transfer function is calculated for each successive segment and averaged
    // across segments to form the Transfer function.
    //
    // Input:
    //   Input   : Input matrix of single column.
    //   Output  : Output matrix of multiple columns
    //   RWL     : Running window length (seconds)
    //   FSamp   : Sampling frequency of the Input and Output matrices
    //   OLR     : Overlap ratio between two successive segments
    //
    // Output:
    //   TF   : Transfer function - Same size as Output matrix
    //     f  : Frequency vector over which the TF is calculated
    //
    // by Dr. Yavuz Kaya, P.Eng.;
    // Created on      02.March.2019
    // Last modified   03.March.2019

}

function Data_V() {

    return [
        [-1.22696384381571	,  1.07804587481937	    , -1.52196836203740	,  0.549349147991792	,  1.80448013956929],
        [-0.327239140494225	,  0.308178759981126	, 0.621011307663875 , 	0.467634317562715	,  -0.632060379431905],
        [0.891646400708002	,  0.299638542702110	, -1.50752258512374	,  0.191457003009039	,  1.31645820719660],
        [0.288185578130522	,  -0.197218455541783	, -1.67940213835796	,  -0.229804955766629	,  1.55159205849066],
        [2.26519543357233	,  -0.146422432942845	, 0.788976200730389	,  -0.579216839923128	,  -1.46890810350677],
        [-0.0478903304637657	,  -0.103070869778914	, -0.654282360534659, 	0.480481177166826	,  0.176938356086965],
        [-1.55186406853224	,  -2.79898194583242	, 1.24490301031691	,   -0.386832641866730	,  3.46625991439170],
        [0.444083312031626	,  0.393275176045361	, -1.29226441786498	,   0.421615541369374	,  -0.214626724445032],
        [-0.911818876809533	,  0.990247152651875	, -0.614352473622596, 	1.08770585526875	,  0.486282563419722],
        [0.0494346475387195	,  -1.29763291067202	, 0.241695211355152	,   -2.24933189844400	,  0.330885119922978]
    ];

}

function Data500() {
    return [
        0.0844603174995124         ,
        1.58489881951931           ,
        1.14191473571310           ,
        -1.35238543142068          ,
        -0.161949703835975         ,
        -1.22773381992825          ,
        -0.808465693280332         ,
        0.357118124475001          ,
        0.602699894904838          ,
        0.250279042639578          ,
        0.220123212754366          ,
        2.71336478210263           ,
        -0.282417558950304         ,
        0.886748365919865          ,
        0.331160874951649          ,
        1.97683976656366           ,
        -0.592803621717091         ,
        0.691035310906601          ,
        -2.01267551063306          ,
        -1.39440286294597          ,
        -0.900861802981435         ,
        -0.332096863213568         ,
        0.647135629496336          ,
        0.00571869556480775        ,
        1.32620784539887           ,
        -1.28066820801125          ,
        -1.48890876789754          ,
        0.0406288538947070         ,
        0.182261943151885          ,
        0.207474438868625          ,
        0.256851706288613          ,
        -0.613609946295178         ,
        0.0111314109387508         ,
        -1.08718231688222          ,
        -0.712838296912265         ,
        0.735694152231514          ,
        -0.362444063967712         ,
        0.139300253852541          ,
        -0.287166050813671         ,
        -1.18285914119541          ,
        -0.645753330465427         ,
        0.461895536354062          ,
        -0.0499236131354404        ,
        0.340460553296415          ,
        -0.561811966875045         ,
        0.786128343195693          ,
        -0.0850711050018647        ,
        -0.0854545522344731        ,
        -0.452883039599746         ,
        1.10407864791520           ,
        -0.823045556266867         ,
        1.63854699583821           ,
        0.565146161655363          ,
        -0.385904976706853         ,
        -0.248886363657854         ,
        -0.147972865120699         ,
        0.456650149696909          ,
        0.771789829912884          ,
        -0.265244762646484         ,
        1.02473644259262           ,
        2.03254087370090           ,
        -0.110671616550137         ,
        0.587622041677598          ,
        -0.603574029970699         ,
        0.743376151981725          ,
        -1.26287771896791          ,
        -1.14789098815766          ,
        -0.502960245348216         ,
        1.38815094833161           ,
        0.0340621985023298         ,
        0.331627332353053          ,
        0.749882026270739          ,
        0.184411291772378          ,
        -0.133147451274978         ,
        -0.126917070737305         ,
        -0.508294862416309         ,
        0.822546698728909          ,
        1.80818414975520           ,
        -0.727543590917313         ,
        -0.919579185444817         ,
        0.538286148450088          ,
        -1.82926721080124          ,
        2.12403498762405           ,
        -1.04139230556933          ,
        -0.843095558085898         ,
        -0.576964952073649         ,
        -1.24870710595472          ,
        0.0886233309130156         ,
        -1.11980366747178          ,
        1.55603658727517           ,
        -0.269023452811481         ,
        0.226436955216374          ,
        -1.07584079808804          ,
        -0.381094602536626         ,
        -0.414955102491428         ,
        -1.21565399415450          ,
        -0.133232316097266         ,
        0.231609250903858          ,
        -0.234627669743434         ,
        0.463927413378933          ,
        -0.484680544879962         ,
        0.737342904985498          ,
        1.31764050593737           ,
        0.610845522451412          ,
        1.06027141932350           ,
        -0.271025532665868         ,
        0.480599767654982          ,
        -1.15123014576584          ,
        -0.691460165945098         ,
        -0.792938439984827         ,
        -0.472164667772703         ,
        0.614933503786697          ,
        -1.28348442687361          ,
        1.47263598840625           ,
        0.879090034648008          ,
        1.92340391603222           ,
        0.121297697280397          ,
        0.452122487725470          ,
        1.51374311908824           ,
        0.113565735264003          ,
        -0.222393140931744         ,
        -0.706992648161236         ,
        -0.469177256142784         ,
        0.778856892459288          ,
        0.0464883755516863         ,
        0.139522521211399          ,
        -1.38081942593971          ,
        0.182805477412276          ,
        1.87857040722629           ,
        0.754804170266752          ,
        0.949789160464642          ,
        -0.0267954233894763        ,
        -1.71504042199688          ,
        -2.33471004097728          ,
        0.685877753845871          ,
        -0.983726287768541         ,
        0.440093078102060          ,
        0.755712928297109          ,
        0.595477311414294          ,
        1.03359863732638           ,
        1.20521540196228           ,
        -0.287072663457633         ,
        0.870399103594100          ,
        1.87725126470526           ,
        -2.97957712564628          ,
        -1.24916324064870          ,
        1.07297224792392           ,
        -0.432679372256279         ,
        -1.06793713159694          ,
        0.423024635030762          ,
        0.188891220128056          ,
        2.29957855593735           ,
        -0.432722290100711         ,
        0.513137536121674          ,
        -1.68707768409984          ,
        -0.687904306634916         ,
        -1.04524349890623          ,
        1.34305048741140           ,
        0.114161070368347          ,
        -0.488921411756530         ,
        -1.02992971074612          ,
        0.679733092291711          ,
        -1.78168512252772          ,
        1.10841551947173           ,
        0.240379781553238          ,
        -0.0555494762342227        ,
        -1.35705482833232          ,
        -2.33795081145851          ,
        0.266751883053658          ,
        -2.13307662186752          ,
        -0.0302561498147974        ,
        -0.0617570528346412        ,
        0.0251442212354898         ,
        0.0445340315003364         ,
        -0.664670327755326         ,
        -0.159025782987608         ,
        -0.243195483880197         ,
        0.784496756564593          ,
        -0.512573071272982         ,
        -1.15223820177729          ,
        -1.28179030121212          ,
        -1.37126374089868          ,
        1.01204549147543           ,
        -0.704766516611475         ,
        -0.465840483749377         ,
        -0.677349453625763         ,
        -0.500957906878932         ,
        -0.804844807073594         ,
        0.934208475807774          ,
        -0.431460357364588         ,
        0.901340812559861          ,
        -0.510979192621600         ,
        -0.265642105725161         ,
        0.658107576861553          ,
        1.64017979054324           ,
        0.976986656472511          ,
        1.39237780192702           ,
        -0.0772003313956825        ,
        -0.688470530383388         ,
        0.261892005769165          ,
        -1.12861587995265          ,
        0.965646248007031          ,
        -1.36178293882454          ,
        0.318969512380783          ,
        -0.178769247162191         ,
        0.00706035868074639        ,
        1.62607607507216           ,
        0.651769346907816          ,
        0.0622587274866813         ,
        1.04554911563484           ,
        0.269444204147106          ,
        0.233288482076063          ,
        0.507118585280466          ,
        -0.449760897224305         ,
        0.143105599733476          ,
        -0.770045852787616         ,
        -1.97764586739574          ,
        -0.840246919983134         ,
        -0.677175549021427         ,
        0.899075014575140          ,
        -0.389596201270305         ,
        0.999054087627100          ,
        0.797063307579231          ,
        -0.468338444480064         ,
        -0.921062699495156         ,
        -0.401401246872318         ,
        -1.21124180636743          ,
        1.22033083796495           ,
        0.489694114684168          ,
        -0.0376231226272725        ,
        1.99843406269389           ,
        0.235672922676430          ,
        1.13066649378719           ,
        -0.509653857682602         ,
        -0.425166990170222         ,
        -0.659651847519786         ,
        -1.65684536140666          ,
        -0.549747918864291         ,
        -0.723683269080522         ,
        -0.208409141354759         ,
        -0.439816301778245         ,
        1.22739201128968           ,
        -0.631893942589536         ,
        -0.688405825813722         ,
        1.68944418672528           ,
        -0.286531922158002         ,
        0.126133199250934          ,
        1.43337122601953           ,
        -0.290861308952875         ,
        1.54254641763595           ,
        -0.188386181361206         ,
        0.856138544618607          ,
        -1.43554894862110          ,
        -0.0161771389276517        ,
        0.798219759324313          ,
        -0.686601015888339         ,
        -1.14138425589188          ,
        0.173006977924011          ,
        0.141391367008740          ,
        -1.57773285983536          ,
        -1.36296186508612          ,
        1.10141738727393           ,
        0.579610843239650          ,
        0.0703982735740168         ,
        0.476287333455335          ,
        0.760722356138319          ,
        0.404250882746384          ,
        1.18093630067466           ,
        0.480354594985207          ,
        -2.07485916370505          ,
        -0.387995204047038         ,
        0.576167090196590          ,
        -0.616026122447960         ,
        -1.09288371332933          ,
        1.38684767121464           ,
        0.835835605703340          ,
        -0.512107642852344         ,
        -1.07457796022912          ,
        1.79046977009799           ,
        0.659876580222513          ,
        -0.203521157273675         ,
        -1.29979361247396          ,
        0.916346328316059          ,
        1.66995977981522           ,
        -0.518806870589207         ,
        0.235339555450697          ,
        0.687921281535145          ,
        -2.36741268795267          ,
        0.477633958412511          ,
        0.171681846875096          ,
        -0.358260527169860         ,
        1.92936600849405           ,
        0.495298019587797          ,
        1.69343001278806           ,
        1.81105191522048           ,
        0.0120369562479433         ,
        -1.35587031834132          ,
        1.05700712315876           ,
        0.110558275361712          ,
        0.0251172880598136         ,
        2.02935343620480           ,
        0.620757658235027          ,
        -0.549396880625381         ,
        2.02964413221311           ,
        0.393263629274566          ,
        0.187767227047191          ,
        -0.438321152080045         ,
        0.715674275634101          ,
        -2.67416611312864          ,
        -0.189848687267985         ,
        0.427719529869678          ,
        -0.572371832328432         ,
        -1.49288841594528          ,
        0.648755423579811          ,
        -1.16259938929584          ,
        -0.351513693454996         ,
        0.149800447307552          ,
        -0.620307695936027         ,
        0.934192330045269          ,
        -0.718416727309709         ,
        1.26612054474156           ,
        -0.559792568035452         ,
        -0.285938350780082         ,
        -0.854151746051820         ,
        -0.728698377443024         ,
        0.858980829971565          ,
        -1.32299935394893          ,
        1.70530147215672           ,
        -1.21248264849629          ,
        -0.751638722621229         ,
        0.859947807684077          ,
        1.07734867444062           ,
        1.66649064143326           ,
        0.409368936409304          ,
        -0.298259082670436         ,
        -1.34237251926560          ,
        0.0567690240582586         ,
        0.674469946627434          ,
        0.378563596204813          ,
        1.89525260030615           ,
        2.36871200116630           ,
        0.424509988937863          ,
        1.17040574119634           ,
        0.846226557460673          ,
        -0.215994300293751         ,
        -2.25276435706569          ,
        -0.375113950621171         ,
        -2.12525478281230          ,
        -0.540678177385886         ,
        0.531709709552514          ,
        -0.208395481957775         ,
        1.36195716259413           ,
        -3.06758436563190          ,
        0.891450691996614          ,
        -1.46644479253244          ,
        2.17192597349202           ,
        0.0253583858314756         ,
        0.0132258875971179         ,
        1.35642186664452           ,
        -1.23377685216833          ,
        -1.51444563040441          ,
        0.791411108317160          ,
        2.47053643689727           ,
        0.942575001816039          ,
        0.558140823835601          ,
        0.508065326145460          ,
        -0.303896568598906         ,
        -0.177427512715403         ,
        0.181846806994695          ,
        2.39334570658494           ,
        -1.71536929982915          ,
        0.408959972384774          ,
        0.324688859105424          ,
        -1.07132465381623          ,
        0.329504989997233          ,
        0.208794833525136          ,
        -0.774576798768013         ,
        -0.0568302900924565        ,
        -1.28105010751899          ,
        -0.835290002591060         ,
        0.250303449659204          ,
        0.535190297144180          ,
        -1.06817411013383          ,
        -1.95416443843637          ,
        -0.615917505887273         ,
        -0.981513732349280         ,
        0.519718627051927          ,
        0.773009597890718          ,
        -0.270269010627654         ,
        -1.40104011654006          ,
        -0.290199658234517         ,
        -0.0209496149083120        ,
        -0.0258419034181870        ,
        -0.206935361840776         ,
        -0.852258176593703         ,
        2.99509332206469           ,
        0.100373435329947          ,
        0.468569377235733          ,
        0.0426769723396555         ,
        1.36809942050754           ,
        0.108211389246283          ,
        0.264507406199269          ,
        1.59805237151102           ,
        -1.62672904913716          ,
        -0.0870039046856149        ,
        -0.695667403654184         ,
        -0.115854930499358         ,
        -0.0751020592755356        ,
        -0.359421978654797         ,
        0.748644595979229          ,
        -0.699422862559458         ,
        -0.936703315408536         ,
        1.81200845550490           ,
        1.41536574955107           ,
        1.15943166548398           ,
        2.58329080048668           ,
        1.92809236241687           ,
        -0.491577816448374         ,
        0.863429166675993          ,
        0.595162653876330          ,
        -0.00698142378955546       ,
        1.31560149496539           ,
        1.05087796709405           ,
        1.11421936787212           ,
        -1.13915719773180          ,
        0.482076856317395          ,
        -1.29505342487988          ,
        -0.370672796578085         ,
        -0.438477397008251         ,
        0.0912340589656165         ,
        -0.642797330030394         ,
        1.86073692461133           ,
        -1.15337146901165          ,
        0.428963773101857          ,
        0.731183627561607          ,
        -0.494760187913981         ,
        0.366297354432319          ,
        0.121565323754563          ,
        -1.67071543565691          ,
        -0.128262130477708         ,
        2.05793903367714           ,
        -1.00711638693528          ,
        -1.42793256466421          ,
        0.114410897091231          ,
        -0.525012553949627         ,
        1.72284086868408           ,
        -0.474782059938199         ,
        -1.57173068052326          ,
        0.417779833203497          ,
        -0.504241350678695         ,
        -2.04243009624114          ,
        0.602132488521351          ,
        -0.0819766631435109        ,
        0.000108929764751487       ,
        -1.01869796118760          ,
        0.753489924822795          ,
        -0.0842281626030491        ,
        0.737873033982041          ,
        -1.15230745495358          ,
        1.99692073682035           ,
        -0.449275586391165         ,
        0.0611206439733801         ,
        -0.778902920880378         ,
        -2.59508655511599          ,
        2.67740663798458           ,
        0.0732079666256241         ,
        -1.03470421681874          ,
        -0.0223440923786126        ,
        -0.404853430083130         ,
        0.331029021236492          ,
        -0.787448559736036         ,
        -0.650747723262264         ,
        -0.647302275288246         ,
        -1.09283702174396          ,
        -0.945528145362311         ,
        -0.113898342155809         ,
        1.07793519441805           ,
        -1.81112062789967          ,
        0.685248667196579          ,
        1.66814939766873           ,
        -0.609182627190437         ,
        0.518914513588425          ,
        -0.205884036495407         ,
        -2.62345313722010          ,
        0.331926649515461          ,
        -0.974676325621187         ,
        -0.710677701935702         ,
        -1.28712698789762          ,
        0.0616655910573507         ,
        0.600250883538885          ,
        -1.08181009701912          ,
        -0.228794374223121         ,
        -0.517896904880999         ,
        0.713311888286190          ,
        0.498610988319403          ,
        -0.282849494396830         ,
        0.265920819509207          ,
        -0.622814504274104         ,
        0.367843086304039          ,
        -0.931664824296078         ,
    ];
}