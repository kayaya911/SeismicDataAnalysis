// stricter parsing and error handling
"use strict";


//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------

class Channel {
    constructor() {

        let Unique_ID = Generate_Unique_ID();

        this.FileName             = undefined;         // FileName
        this.ChNum                = undefined;         // Channel number
        this.Unique_ID            = Unique_ID;         // Unique ID of CheckBox
        this.Selected             = true;              // Weather this channle is slected by user or NOT (true / false)
        this.NumSamples           = undefined;         // Number of digitized data points
        this.ScaleFactor          = undefined;         // User Scale Factor - Taken 1.0 if not specified
        this.Orientation          = undefined;         // Orientation
        this.DateTime             = undefined;         // Date&Time of the first sample in the records
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
            P_Wave: [],
        };

    }

    // Statisctis update
    async Statistics_update() {
        let Temp_Data  =  Mult(this.data, this.ScaleFactor);
        this.Peak      =  Max( Abs(Temp_Data) ).val;
        this.Mean      =  Mean( Temp_Data );
        this.RMS       =  Rms( Temp_Data );
    }

    // PWave this.data
    async PWave(PWavePar) {

        // Decleration of variables 
        let Data, Result, i, Res, numBins, message;

        // Copy number of histogram-bins in a separate variable
        numBins = Copy(PWavePar.numBins);

        // Get the RawData using the scaling factor 
        Data = await this.#Get_ChannelData(this);

        if (PWavePar.PWave_Method == 0) {
            // PhasePicker using Histogram Method 
            Result = PhasePicker_Hist(Data, this.delt, this.Type, PWavePar.ksi, PWavePar.Period, numBins);

            // If PWave Not Detected, try for a maximun of 2 more times with different numBins
            if (Result.PWaveArrival == undefined) {

                message = this.FileName + "\xA0 (Ch: " + this.ChNum + ")" + "\xA0 (FSamp: " + this.FSamp + ") \xA0 (Phase-Picker) \xA0 (numBins=" + numBins + ") P-Wave Not Detected!";
                InfoBar("WARNING", AnalysisType, message);
                
                for (i=0; i < 2; i++) {
                    numBins = Math.ceil(numBins / 2);  
                    if (numBins >= 10 ) { Res = PhasePicker_Hist(Data, this.delt, this.Type, PWavePar.ksi, PWavePar.Period, numBins); }
                    if (Res.PWaveArrival != undefined) { 
                        message = this.FileName + "\xA0 (Ch: " + this.ChNum + ")" + "\xA0 (FSamp: " + this.FSamp + ") \xA0 (Phase-Picker) \xA0 (numBins=" + numBins + ") - Successful.";
                        InfoBar("", AnalysisType, message);
                        Result=Res; 
                        break; 
                    }
                    else {
                        message = this.FileName + "\xA0 (Ch: " + this.ChNum + ")" + "\xA0 (FSamp: " + this.FSamp + ") \xA0 (Phase-Picker) \xA0 (numBins=" + numBins + ") P-Wave Not Detected!";
                        InfoBar("WARNING", AnalysisType, message);
                    }
                }
            }
            else {
                message = this.FileName + "\xA0 (Ch: " + this.ChNum + ")" + "\xA0 (FSamp: " + this.FSamp + ") \xA0 (Phase-Picker) \xA0 - Successful.";
                InfoBar("", AnalysisType, message);
            }
        }
        else if (PWavePar.PWave_Method == 1) {
            // LTA-STA method 
        }
        else if (PWavePar.PWave_Method == 2) {
            // AIC method 
        }
        
        // Store results
        this.Results.P_Wave = {
            PWave_Method    : 0,
            RawData         : Detrend(Data, 0),
            PWaveArrival    : Result.PWaveArrival,
            FilterOrder     : Result.FilterOrder,
            F1              : Result.F1,
            F2              : Result.F2,
            Period          : Result.Period,
            Damping         : Result.Damping,
            numBins         : Result.numBins,
        };

    }

    // Private methods
    async #Get_ChannelData() {
        // Returns the RawData after applying the user-defined scaling factor for that channel
        return Mult(this.data, this.ScaleFactor);
    }

}

//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
// Type and Unit functions

function IntervalTypeAndUnit(n, Opt) {
    // The supported Interval types and units are given below
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
function TypeAndUnit(n, Opt) {
    // The supported measurement types and units are given below
    //
    //  n     ==>     Type                          Units
    //  ----------------------------------------------------------------
    //  1     ==>     0   ( Acceleration   )        0   (g)
    //  2     ==>     0   ( Acceleration   )        1   (m/s²)
    //  3     ==>     0   ( Acceleration   )        2   (cm/s²)
    //  4     ==>     0   ( Acceleration   )        3   (mm/s²)
    //  11    ==>     1   ( Velocity       )        4   (gs)
    //  12    ==>     1   ( Velocity       )        5   (m/s)
    //  13    ==>     1   ( Velocity       )        6   (cm/s)
    //  14    ==>     1   ( Velocity       )        7   (mm/s)
    //  21    ==>     2   ( Displacement   )        8   (gs²)
    //  22    ==>     2   ( Displacement   )        9   (m)
    //  23    ==>     2   ( Displacement   )        10  (cm)
    //  24    ==>     2   ( Displacement   )        11  (mm)
    //  31    ==>     3   ( Strain         )        12  (Strain)
    //  41    ==>     4   ( Wind Direction )        13  (Degrees)
    //  51    ==>     5   ( Wind Speed     )        14  (Mile/h)
    //  52    ==>     5   ( Wind Speed     )        15  (km/h)
    //  53    ==>     5   ( Wind Speed     )        16  (m/s)
    //  54    ==>     5   ( Wind Speed     )        17  (Beaufort)
    //  61    ==>     6   ( Temperature    )        18  (°C)
    //  62    ==>     6   ( Temperature    )        19  (°F)
    //  71    ==>     7   ( Humidity       )        20  (%)
    //  81    ==>     8   ( Pore Pressure  )        21  (N/m²)
    //  121   ==>     9   ( Tilt           )        22  (Degrees)
    //  122   ==>     9   ( Tilt           )        23  (Volts)

    if (Opt == null) { Opt = true; }

    if (Opt) {
        // TypeAndUnit number (n) is given
        if      (n == 1)   { return { 'TypeAndUnit': 1  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 0 ,   'Unit_String': 'g'        }; }   //  ( Acceleration   )  (g)
        else if (n == 2)   { return { 'TypeAndUnit': 2  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 1 ,   'Unit_String': 'm/s²'     }; }   //  ( Acceleration   )  (m/s²)
        else if (n == 3)   { return { 'TypeAndUnit': 3  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 2 ,   'Unit_String': 'cm/s²'    }; }   //  ( Acceleration   )  (cm/s²)
        else if (n == 4)   { return { 'TypeAndUnit': 4  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 3 ,   'Unit_String': 'mm/s²'    }; }   //  ( Acceleration   )  (mm/s²)
        else if (n == 11)  { return { 'TypeAndUnit': 11 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 4 ,   'Unit_String': 'gs'       }; }   //  ( Velocity       )  (gs)
        else if (n == 12)  { return { 'TypeAndUnit': 12 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 5 ,   'Unit_String': 'm/s'      }; }   //  ( Velocity       )  (m/s)
        else if (n == 13)  { return { 'TypeAndUnit': 13 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 6 ,   'Unit_String': 'cm/s'     }; }   //  ( Velocity       )  (cm/s)
        else if (n == 14)  { return { 'TypeAndUnit': 14 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 7 ,   'Unit_String': 'mm/s'     }; }   //  ( Velocity       )  (mm/s)
        else if (n == 21)  { return { 'TypeAndUnit': 21 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 8 ,   'Unit_String': 'gs²'      }; }   //  ( Displacement   )  (gs²)
        else if (n == 22)  { return { 'TypeAndUnit': 22 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 9 ,   'Unit_String': 'm'        }; }   //  ( Displacement   )  (m)
        else if (n == 23)  { return { 'TypeAndUnit': 23 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 10,   'Unit_String': 'cm'       }; }   //  ( Displacement   )  (cm)
        else if (n == 24)  { return { 'TypeAndUnit': 24 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 11,   'Unit_String': 'mm'       }; }   //  ( Displacement   )  (mm)
        else if (n == 31)  { return { 'TypeAndUnit': 31 ,   'Type': 3,   'Type_String': 'Strain',         'Unit': 12,   'Unit_String': 'Strain'   }; }   //  ( Strain         )  (Strain)
        else if (n == 41)  { return { 'TypeAndUnit': 41 ,   'Type': 4,   'Type_String': 'Wind Direction', 'Unit': 13,   'Unit_String': 'Degrees'  }; }   //  ( Wind Direction )  (Degrees)
        else if (n == 51)  { return { 'TypeAndUnit': 51 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 14,   'Unit_String': 'Mile/h'   }; }   //  ( Wind Speed     )  (Mile/h)
        else if (n == 52)  { return { 'TypeAndUnit': 52 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 15,   'Unit_String': 'km/h'     }; }   //  ( Wind Speed     )  (km/h)
        else if (n == 53)  { return { 'TypeAndUnit': 53 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 16,   'Unit_String': 'm/s'      }; }   //  ( Wind Speed     )  (m/s)
        else if (n == 54)  { return { 'TypeAndUnit': 54 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 17,   'Unit_String': 'Beaufort' }; }   //  ( Wind Speed     )  (Beaufort)
        else if (n == 61)  { return { 'TypeAndUnit': 61 ,   'Type': 6,   'Type_String': 'Temperature',    'Unit': 18,   'Unit_String': '°C'       }; }   //  ( Temperature    )  (°C)
        else if (n == 62)  { return { 'TypeAndUnit': 62 ,   'Type': 6,   'Type_String': 'Temperature',    'Unit': 19,   'Unit_String': '°F'       }; }   //  ( Temperature    )  (°F)
        else if (n == 71)  { return { 'TypeAndUnit': 71 ,   'Type': 7,   'Type_String': 'Humidity',       'Unit': 20,   'Unit_String': '%'        }; }   //  ( Humidity       )  (%)
        else if (n == 81)  { return { 'TypeAndUnit': 81 ,   'Type': 8,   'Type_String': 'Pore Pressure',  'Unit': 21,   'Unit_String': 'N/m²'     }; }   //  ( Pore Pressure  )  (N/m²)
        else if (n == 121) { return { 'TypeAndUnit': 121,   'Type': 9,   'Type_String': 'Tilt',           'Unit': 22,   'Unit_String': 'Degrees'  }; }   //  ( Tilt           )  (Degrees)
        else if (n == 122) { return { 'TypeAndUnit': 122,   'Type': 9,   'Type_String': 'Tilt',           'Unit': 23,   'Unit_String': 'Volts'    }; }   //  ( Tilt           )  (Volts)
    }
    else {
        // Unit number is given.
        if       ( n == 0  ) { return { 'TypeAndUnit': 1  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 0 ,   'Unit_String': 'g'        }; }   //  ( Acceleration   )  (g)
        else if  ( n == 1  ) { return { 'TypeAndUnit': 2  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 1 ,   'Unit_String': 'm/s²'     }; }   //  ( Acceleration   )  (m/s²)
        else if  ( n == 2  ) { return { 'TypeAndUnit': 3  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 2 ,   'Unit_String': 'cm/s²'    }; }   //  ( Acceleration   )  (cm/s²)
        else if  ( n == 3  ) { return { 'TypeAndUnit': 4  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 3 ,   'Unit_String': 'mm/s²'    }; }   //  ( Acceleration   )  (mm/s²)
        else if  ( n == 4  ) { return { 'TypeAndUnit': 11 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 4 ,   'Unit_String': 'gs'       }; }   //  ( Velocity       )  (gs)
        else if  ( n == 5  ) { return { 'TypeAndUnit': 12 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 5 ,   'Unit_String': 'm/s'      }; }   //  ( Velocity       )  (m/s)
        else if  ( n == 6  ) { return { 'TypeAndUnit': 13 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 6 ,   'Unit_String': 'cm/s'     }; }   //  ( Velocity       )  (cm/s)
        else if  ( n == 7  ) { return { 'TypeAndUnit': 14 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 7 ,   'Unit_String': 'mm/s'     }; }   //  ( Velocity       )  (mm/s)
        else if  ( n == 8  ) { return { 'TypeAndUnit': 21 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 8 ,   'Unit_String': 'gs²'      }; }   //  ( Displacement   )  (gs²)
        else if  ( n == 9  ) { return { 'TypeAndUnit': 22 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 9 ,   'Unit_String': 'm'        }; }   //  ( Displacement   )  (m)
        else if  ( n == 10 ) { return { 'TypeAndUnit': 23 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 10,   'Unit_String': 'cm'       }; }   //  ( Displacement   )  (cm)
        else if  ( n == 11 ) { return { 'TypeAndUnit': 24 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 11,   'Unit_String': 'mm'       }; }   //  ( Displacement   )  (mm)
        else if  ( n == 12 ) { return { 'TypeAndUnit': 31 ,   'Type': 3,   'Type_String': 'Strain',         'Unit': 12,   'Unit_String': 'Strain'   }; }   //  ( Strain         )  (Strain)
        else if  ( n == 13 ) { return { 'TypeAndUnit': 41 ,   'Type': 4,   'Type_String': 'Wind Direction', 'Unit': 13,   'Unit_String': 'Degrees'  }; }   //  ( Wind Direction )  (Degrees)
        else if  ( n == 14 ) { return { 'TypeAndUnit': 51 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 14,   'Unit_String': 'Mile/h'   }; }   //  ( Wind Speed     )  (Mile/h)
        else if  ( n == 15 ) { return { 'TypeAndUnit': 52 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 15,   'Unit_String': 'km/h'     }; }   //  ( Wind Speed     )  (km/h)
        else if  ( n == 16 ) { return { 'TypeAndUnit': 53 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 16,   'Unit_String': 'm/s'      }; }   //  ( Wind Speed     )  (m/s)
        else if  ( n == 17 ) { return { 'TypeAndUnit': 54 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 17,   'Unit_String': 'Beaufort' }; }   //  ( Wind Speed     )  (Beaufort)
        else if  ( n == 18 ) { return { 'TypeAndUnit': 61 ,   'Type': 6,   'Type_String': 'Temperature',    'Unit': 18,   'Unit_String': '°C'       }; }   //  ( Temperature    )  (°C)
        else if  ( n == 19 ) { return { 'TypeAndUnit': 62 ,   'Type': 6,   'Type_String': 'Temperature',    'Unit': 19,   'Unit_String': '°F'       }; }   //  ( Temperature    )  (°F)
        else if  ( n == 20 ) { return { 'TypeAndUnit': 71 ,   'Type': 7,   'Type_String': 'Humidity',       'Unit': 20,   'Unit_String': '%'        }; }   //  ( Humidity       )  (%)
        else if  ( n == 21 ) { return { 'TypeAndUnit': 81 ,   'Type': 8,   'Type_String': 'Pore Pressure',  'Unit': 21,   'Unit_String': 'N/m²'     }; }   //  ( Pore Pressure  )  (N/m²)
        else if  ( n == 22 ) { return { 'TypeAndUnit': 121,   'Type': 9,   'Type_String': 'Tilt',           'Unit': 22,   'Unit_String': 'Degrees'  }; }   //  ( Tilt           )  (Degrees)
        else if  ( n == 23 ) { return { 'TypeAndUnit': 122,   'Type': 9,   'Type_String': 'Tilt',           'Unit': 23,   'Unit_String': 'Volts'    }; }   //  ( Tilt           )  (Volts)
    }
}
function UnitConversion(In, Out, Opt) {
    // Opt: true  ==> TypeAndUnit (n) is given
    // Opt: false ==> Unit number is given

    let iU, oU;

    if (Opt == null) { Opt = true; }

    iU = TypeAndUnit(In,  Opt);
    oU = TypeAndUnit(Out, Opt);

    if      (iU.Unit == 0) {
        // Acceleration, g
        if      (oU.Unit == 0) { return 1;      }  // g
        else if (oU.Unit == 1) { return 9.81;   }  // m/s²
        else if (oU.Unit == 2) { return 981;    }  // cm/s²
        else if (oU.Unit == 3) { return 9810;   }  // mm/s²
    }
    else if (iU.Unit == 1) {
        // Acceleration, m/s²
        if      (oU.Unit == 0) { return 1/9.81; }  // g
        else if (oU.Unit == 1) { return 1;      }  // m/s²
        else if (oU.Unit == 2) { return 100;    }  // cm/s²
        else if (oU.Unit == 3) { return 1000;   }  // mm/s²
    }
    else if (iU.Unit == 2) {
        // Acceleration, cm/s²
        if      (oU.Unit == 0) { return 1/981;  }  // g
        else if (oU.Unit == 1) { return 0.01;   }  // m/s²
        else if (oU.Unit == 2) { return 1;      }  // cm/s²
        else if (oU.Unit == 3) { return 10;     }  // mm/s²
    }
    else if (iU.Unit == 3) {
        // Acceleration, mm/s²
        if      (oU.Unit == 0) { return 1/9810; }  // g
        else if (oU.Unit == 1) { return 0.001;  }  // m/s²
        else if (oU.Unit == 2) { return 0.1;    }  // cm/s²
        else if (oU.Unit == 3) { return 1;      }  // mm/s²
    }
    else if (iU.Unit == 4) {
        // Velocity, gs
        if      (oU.Unit == 4) { return 1;      }  // gs
        else if (oU.Unit == 5) { return 9.81;   }  // m/s
        else if (oU.Unit == 6) { return 981;    }  // cm/s
        else if (oU.Unit == 7) { return 9810;   }  // mm/s
    }
    else if (iU.Unit == 5) {
        // Velocity, m/s
        if      (oU.Unit == 4) { return 1/9.81; }  // g
        else if (oU.Unit == 5) { return 1;      }  // m/s
        else if (oU.Unit == 6) { return 100;    }  // cm/s
        else if (oU.Unit == 7) { return 1000;   }  // mm/s
    }
    else if (iU.Unit == 6) {
        // Velocity, cm/s
        if      (oU.Unit == 4) { return 1/981;  }  // g
        else if (oU.Unit == 5) { return 0.01;   }  // m/s
        else if (oU.Unit == 6) { return 1;      }  // cm/s
        else if (oU.Unit == 7) { return 10;     }  // mm/s
    }
    else if (iU.Unit == 7) {
        // Velocity, mm/s
        if      (oU.Unit == 4) { return 1/9810; }  // g
        else if (oU.Unit == 5) { return 0.001;  }  // m/s
        else if (oU.Unit == 6) { return 0.1;    }  // cm/s
        else if (oU.Unit == 7) { return 1;      }  // mm/s
    }
    else if (iU.Unit == 8) {
        // Displacement, gs²
        if      (oU.Unit == 8)  { return 1;      }  // s²
        else if (oU.Unit == 9)  { return 9.81;   }  // m
        else if (oU.Unit == 10) { return 981;    }  // cm
        else if (oU.Unit == 11) { return 9810;   }  // mm
    }
    else if (iU.Unit == 9) {
        // Displacement, m
        if      (oU.Unit == 8)  { return 1/9.81; }  // gs²
        else if (oU.Unit == 9)  { return 1;      }  // m
        else if (oU.Unit == 10) { return 100;    }  // cm
        else if (oU.Unit == 11) { return 1000;   }  // mm
    }
    else if (iU.Unit == 10) {
        // Displacement, cm
        if      (oU.Unit == 8)  { return 1/981;  }  // gs²
        else if (oU.Unit == 9)  { return 0.01;   }  // m
        else if (oU.Unit == 10) { return 1;      }  // cm
        else if (oU.Unit == 11) { return 10;     }  // mm
    }
    else if (iU.Unit == 11) {
        // Displacement, mm
        if      (oU.Unit == 8)  { return 1/9810; }  // gs²
        else if (oU.Unit == 9)  { return 0.001;  }  // m
        else if (oU.Unit == 10) { return 0.1;    }  // cm
        else if (oU.Unit == 11) { return 1;      }  // mm
    }
    else if (iU.Unit == 12) {
        // Strain, strain
        if      (oU.Unit == 12)  { return 1; }  // strain
    }
    else if (iU.Unit == 13) {
        // Wind direction, degrees
        if      (oU.Unit == 13)  { return 1; }  // degrees
    }
    else if (iU.Unit == 14) {
        // Wind Speed, Mile/h
        if      (oU.Unit == 14) { return 1;         }  // Mile/h
        else if (oU.Unit == 15) { return 1.609344;  }  // km/h
        else if (oU.Unit == 16) { return 44.704;    }  // cm/s
        else if (oU.Unit == 17) { return 0;         }  // Beaufort
    }
    else if (iU.Unit == 15) {
        // Wind Speed, km/h
        if      (oU.Unit == 14) { return 0.621371;  }  // Mile/h
        else if (oU.Unit == 15) { return 1.0;       }  // km/h
        else if (oU.Unit == 16) { return 27.7778;   }  // cm/s
        else if (oU.Unit == 17) { return 0;         }  // Beaufort
    }
    else if (iU.Unit == 16) {
        // Wind Speed, cm/s
        if      (oU.Unit == 14) { return 0.0223694; }  // Mile/h
        else if (oU.Unit == 15) { return 0.036;     }  // km/h
        else if (oU.Unit == 16) { return 1;         }  // cm/s
        else if (oU.Unit == 17) { return 0;         }  // Beaufort
    }
    else if (iU.Unit == 17) {
        // Wind Speed, cm/s
        if      (oU.Unit == 14) { return 0.0;       }  // Mile/h
        else if (oU.Unit == 15) { return 0.0;       }  // km/h
        else if (oU.Unit == 16) { return 0;         }  // cm/s
        else if (oU.Unit == 17) { return 1;         }  // Beaufort
    }
    else if (iU.Unit == 18) {
        // temperature, °C
        if      (oU.Unit == 18) { return 1.0;       }  // °C
        else if (oU.Unit == 19) { return 0.0;       }  // °F
    }
    else if (iU.Unit == 19) {
        // temperature, °F
        if      (oU.Unit == 18) { return 0.0;       }  // °C
        else if (oU.Unit == 19) { return 1.0;       }  // °F
    }
    else if (iU.Unit == 20) {
        // Humidity, %
        if      (oU.Unit == 20) { return 1.0;       }  // %
    }
    else if (iU.Unit == 21) {
        // Pore Pressure, N/m²
        if      (oU.Unit == 21) { return 1.0;       }  // N/m²
    }
    else if (iU.Unit == 22) {
        // Tilt, Degrees
        if      (oU.Unit == 22) { return 1.0;       }  // Degrees
        else if (oU.Unit == 23) { return 0.0;       }  // Volts
    }
    else if (iU.Unit == 23) {
        // Tilt, Volts
        if      (oU.Unit == 22) { return 0.0;       }  // Degrees
        else if (oU.Unit == 23) { return 1.0;       }  // Volts
    }
}
function TypeUnitList() {
    // Returns a list of Types and Units
    // N : Number of measurement types
    let N=24, i, temp, Type = [], Unit =[];
    for (i=0; i<N; i++) { temp = TypeAndUnit(i, false); Type.push(temp.Type_String); Unit.push(temp.Unit_String); }
    // Remove duplicates in the Type and Unit arrays
    Type = [...new Set(Type)];
    Unit = [...new Set(Unit)];
    return { 'Types': Type, 'Units': Unit }
}
function UnitList(n, Opt) {
    // N : Number of measurement types

    if (Opt == null) { Opt = true; }

    let N=24, i, temp, temp2, Unit=[], UnitNum=[];

    if (Opt) {
        // Return a list of all Units of a specific unit
        // n : Unit number
        temp = TypeAndUnit(n, false).Type_String;
        for (i=0; i<N; i++) { temp2 = TypeAndUnit(i, false); if (temp == temp2.Type_String) { Unit.push(temp2.Unit_String), UnitNum.push(temp2.Unit); } }
        return { 'Units': Unit, 'UnitNum': UnitNum }

    } else {
        // Return a list of all Units of a specific Type
        // n : Type number
        for (i=0; i<N; i++) { temp2 = TypeAndUnit(i, false); if (n == temp2.Type) { Unit.push(temp2.Unit_String), UnitNum.push(temp2.Unit); } }
        return { 'Units': Unit, 'UnitNum': UnitNum }
    }

}
function GraphUnits(ChNum) {
    
    // Decleration of variables 
    let Or_Data, Gr_Data, ScaleFactor, Ind, result, Unit_Plot_ID, Unit_ID;
    let yLabel_vel="", yLabel_disp="";

    // Original Data Units 
    Or_Data = TypeAndUnit(ChannelList[ChNum].TypeAndUnits);
    
    // User sepecified Unit 
    Unit_Plot_ID = "Unit_Plot_ID_" + ChannelList[ChNum].Unique_ID;
    Unit_ID      = "Unit_ID_" + ChannelList[ChNum].Unique_ID;
    if (PageNo == 0) { Ind = document.getElementById(Unit_ID).selectedIndex;      }
    else             { Ind = document.getElementById(Unit_Plot_ID).selectedIndex; }

    // Graph Units 
    if      (Or_Data.Type == 0) { 
        // Acceleration -- Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 1; } 
        else if (Ind == 1) { Gr_Data = 2; }  
        else if (Ind == 2) { Gr_Data = 3; }
        else if (Ind == 3) { Gr_Data = 4; }

        yLabel_vel  = "<b>Velocity" + "  [" + UnitList(4).Units[Ind] + "] <b>";
        yLabel_disp = "<b>Displacement" + "  [" + UnitList(8).Units[Ind] + "] <b>";
        
    }
    else if (Or_Data.Type == 1) { 
        // Velocity -- Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 11; }
        else if (Ind == 1) { Gr_Data = 12; }
        else if (Ind == 2) { Gr_Data = 13; }
        else if (Ind == 3) { Gr_Data = 14; }
    }
    else if (Or_Data.Type == 2) { 
        // Displacement -- Get the selected Unit Index
        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 21; }
        else if (Ind == 1) { Gr_Data = 22; }
        else if (Ind == 2) { Gr_Data = 23; }
        else if (Ind == 3) { Gr_Data = 24; }
    }
    else if (Or_Data.Type == 3) { 
        // Strain -- Get the selected Unit Index
        // Find the TypeAndUnit number
        Gr_Data = 31;

    }
    else if (Or_Data.Type == 4) { 
        // Wind Direction -- Get the selected Unit Index
        // Find the TypeAndUnit number
        Gr_Data = 41;

    }
    else if (Or_Data.Type == 5) { 
        // Wind Speed -- Get the selected Unit Index
        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 51; }
        else if (Ind == 1) { Gr_Data = 52; }
        else if (Ind == 2) { Gr_Data = 53; }
        else if (Ind == 3) { Gr_Data = 54; }

    }
    else if (Or_Data.Type == 6) { 
        // Temperature -- Get the selected Unit Index
        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 61; }
        else if (Ind == 1) { Gr_Data = 62; }

    }
    else if (Or_Data.Type == 7) { 
        // Humidity -- Get the selected Unit Index
        // Find the TypeAndUnit number
        Gr_Data = 71;

    }
    else if (Or_Data.Type == 8) {
        // PorePressure -- Get the selected Unit Index
        // Find the TypeAndUnit number
        Gr_Data = 81;

     }
    else if (Or_Data.Type == 9) { 
        // Tilt -- Get the selected Unit Index
        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 71; }
        else if (Ind == 1) { Gr_Data = 72; }

    }
    
    ScaleFactor = UnitConversion(ChannelList[ChNum].TypeAndUnits, Gr_Data);
    Gr_Data     = TypeAndUnit(Gr_Data);

    result = {
        Or_Data     : Or_Data, 
        Gr_Data     : Gr_Data, 
        ScaleFactor : ScaleFactor,
        yLabel      : "<b>" + Gr_Data.Type_String + "  [" + Gr_Data.Unit_String + "] <b>",
        yLabel_vel  : yLabel_vel,
        yLabel_disp : yLabel_disp,
        yLabel_FFT  : "<b>Normalized Magnitude<b>",
    };
   
    return result;
}
function GraphUnits_old(ChNum) {
    
    // Decleration of variables 
    let Or_Data, Gr_Data, ScaleFactor, Ind, result;

    // Original Data Units 
    Or_Data = TypeAndUnit(ChannelList[ChNum].TypeAndUnits);

    // Graph Units 
    if      (Or_Data.Type == 0) { 
        
        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_Acceleration").selectedIndex ;
        
        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 1; }
        else if (Ind == 1) { Gr_Data = 2; }
        else if (Ind == 2) { Gr_Data = 3; }
        else if (Ind == 3) { Gr_Data = 4; }

    }
    else if (Or_Data.Type == 1) { 

        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_Velocity").selectedIndex ;

        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 11; }
        else if (Ind == 1) { Gr_Data = 12; }
        else if (Ind == 2) { Gr_Data = 13; }
        else if (Ind == 3) { Gr_Data = 14; }
        
    }
    else if (Or_Data.Type == 2) { 

        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_Displacement").selectedIndex ;

        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 21; }
        else if (Ind == 1) { Gr_Data = 22; }
        else if (Ind == 2) { Gr_Data = 23; }
        else if (Ind == 3) { Gr_Data = 24; }

    }
    else if (Or_Data.Type == 3) { 

        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_Strain").selectedIndex ;

        // Find the TypeAndUnit number
        Gr_Data = 31;

    }
    else if (Or_Data.Type == 4) { 

        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_WindDir").selectedIndex ;

        // Find the TypeAndUnit number
        Gr_Data = 41;

    }
    else if (Or_Data.Type == 5) { 

        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_WindSpeed").selectedIndex ;

        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 51; }
        else if (Ind == 1) { Gr_Data = 52; }
        else if (Ind == 2) { Gr_Data = 53; }
        else if (Ind == 3) { Gr_Data = 54; }

    }
    else if (Or_Data.Type == 6) { 
        
        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_Temperature").selectedIndex ;

        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 61; }
        else if (Ind == 1) { Gr_Data = 62; }

    }
    else if (Or_Data.Type == 7) { 

        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_Humidity").selectedIndex ;

        // Find the TypeAndUnit number
        Gr_Data = 71;

    }
    else if (Or_Data.Type == 8) {

        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_PorePressure").selectedIndex ;

        // Find the TypeAndUnit number
        Gr_Data = 81;

     }
    else if (Or_Data.Type == 9) { 

        // Get the selected Unit Index
        Ind = document.getElementById("Graph_Output_Tilt").selectedIndex ;

        // Find the TypeAndUnit number
        if      (Ind == 0) { Gr_Data = 71; }
        else if (Ind == 1) { Gr_Data = 72; }

    }
    
    ScaleFactor = UnitConversion(ChannelList[ChNum].TypeAndUnits, Gr_Data);

    result = {
        Or_Data     : Or_Data, 
        Gr_Data     : TypeAndUnit(Gr_Data), 
        ScaleFactor : ScaleFactor
    };
   
    return result;
}
function IntegralUnits(Indx) {

    let Ind_Acc, Ind_Vel, Ind_Disp, Or_Data;

    // Get the user-specified units for Acc, Vel, Disp
    Ind_Acc  = document.getElementById("Graph_Output_Acceleration").selectedIndex ;
    Ind_Vel  = document.getElementById("Graph_Output_Velocity"    ).selectedIndex ;
    Ind_Disp = document.getElementById("Graph_Output_Displacement").selectedIndex ;

    // Original Data Units 
    Or_Data = TypeAndUnit( ChannelList[Indx].TypeAndUnits );

    if      (Or_Data.Type == 0) {
        return {
            Acc_Type   : "Acceleration",
            Acc_Unit   : document.getElementById("Graph_Output_Acceleration").value,

            Vel_Type   : "Velocity",
            Vel_Unit   : document.getElementById("Graph_Output_Velocity").value,

            Disp_Type  : "Displacement",
            Disp_Unit  : document.getElementById("Graph_Output_Displacement").value,
        }
    }
    else if (Or_Data.Type == 1) {
        return {
            Acc_Type   : "Velocity",
            Acc_Unit   : document.getElementById("Graph_Output_Velocity").value,

            Vel_Type   : "Displacement",
            Vel_Unit   : document.getElementById("Graph_Output_Displacement").value,

            Disp_Type  : "",
            Disp_Unit  : "",
        }
    }

}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
// Delete Records
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
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
async function Load_Files(ev) {
    // Appends the user-selected files as a Channel Object to the ChannelList[]
    // Files are read based on their file extension (e.g., VIF, V1, RAW, ASC, MSD, etc.)

    let array_files = ev.target.files;
    let i, len, delta, Indexs, message;

    len = array_files.length;

    if (len > 0) { 

        // Reset the ProgressBar
        await ResetProgressBar("ProgressBar_LoadData");
        await sleep(5);

        Indexs = LinSpace(0, len-1, len);
        delta  = 100 / len;
        
        for (i of Indexs) {

            // Skips this file if it is already uploaded.
            if (IsFileUploaded( array_files[i].name )) { 
                IncreaseProgressBar(delta, "ProgressBar_LoadData");  // Update the ProgressBar
                await sleep(5);
                continue; 
            }

            // Reads the content of each file.
            await ReadFileContent( array_files[i] )

            // Waits for 5 miliseccond for secreen update
            await sleep(5);
        }
    }

    function ReadFileContent(File_Blob) {
        return new Promise((resolve) => {
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
                if      ( fileExt == "VIF" )                          { Read_VIF( FileName, delta, dataview ); }
                else if ((fileExt == "V1"  ) || (fileExt == "RAW"))   { Read_V1(  FileName, delta, dataview ); }
                else if ((fileExt == "MSD" ) || (fileExt == "MSEED")) { Read_MSD( FileName, delta, dataview ); }
                else if ( fileExt == "TXT" )                          { Read_TXT( FileName, delta, dataview ); }
                else if ( fileExt == "ASC" )                          { Read_ASC( FileName, delta, dataview ); }
            });
            resolve();
        });    
    }
    
}

function Add_To_FileTree(Channel) {

    // Add New File if not exists
    if (document.getElementById(Channel.FileName) == null) { File_Entry(Channel.FileName, Channel.FileName); }
    
    // Unique Index of Channle
    let Indx = ChannelList_UniqueID(Channel.Unique_ID);
    Indx = "("+ (Indx+1).toString().padStart(3, '0') + ")";

    // Add Channel
    let li_Name = Indx + "   (Ch-" + Channel.ChNum + ") (" +  Channel.Orientation  + ") (" + Number(Channel.FSamp).toFixed(3).toString() + " Hz)";
    Channel_Entry(li_Name);

    // Create Plotly Graph for this channel
    Create_Plotly_Graph("Graphs_LoadData_Container", Channel);
    
    // Update Right_Click Menu
    Update_RightClick_Menu();

    function File_Entry(label_id, TextCont) {
        let UL_el = document.getElementById('FileListTree');

        let li_el = document.createElement('li');
        let label = document.createElement('label');
        
        label.textContent = TextCont;
        label.id          = label_id;
        li_el.className   = "Files";
        li_el.appendChild(Create_SVG());
        li_el.appendChild(label);
        UL_el.appendChild(li_el);
    }

    function Channel_Entry(TextCont) {
        let UL_el = document.getElementById('FileListTree');

        let li_el    = document.createElement('li');
        let checkbox = document.createElement('input');
        let label    = document.createElement('label');

        checkbox.className  = "Parameter-check-input form-check-input File_CheckBox_Ch";
        checkbox.type       = "checkbox";
        checkbox.id         = Channel.Unique_ID;
        checkbox.checked    = true;
        checkbox.setAttribute('onclick',   "Channel_Click(this)" );
 
        label.textContent = TextCont;
        label.setAttribute('for', Channel.Unique_ID);
        label.setAttribute('class', 'LabelClass');
        
        li_el.className  = "Files_Ch"; 
        li_el.appendChild(checkbox);
        li_el.appendChild(label);
        
        UL_el.appendChild(li_el);
    }

    function Create_SVG() {

        //<svg width="1320.000000pt" height="920.000000pt" viewBox="0 0 1320.000000 920.000000" preserveAspectRatio="xMidYMid meet">
        //   <g transform="translate(0.000000,920.000000) scale(0.100000,-0.100000)" fill="#000000" stroke="none">
        //    <path d="M120 4605 l0 -4445 6450 0 6450 0 0 4445 0 4445 -6450 0 -6450 0 0 -4445z m12710 0 l0 -4255 -6260 0 -6260 0 0 4255 0 4255 6260 0 6260 0 0 -4255z"/>
        //    <path d="M5562 8349 c-123 -21 -238 -111 -274 -215 -10 -27 -29 -137 -43 -244 -13 -107 -45 -357 -70 -555 -84 -662 -224 -1764 -310 -2437 -19 -150 -64 -507 -100 -793 -36 -286 -75 -587 -85 -670 -11 -82 -20 -156 -20 -164 0 -31 -12 0 -25 66 -13 59 -164 849 -411 2138 -84 438 -92 464 -179 544 -73 66 -146 94 -246 94 -116 0 -226 -56 -288 -146 -16 -23 -136 -267 -268 -542 -132 -275 -315 -659 -409 -852 l-169 -353 -632 0 -633 0 0 -365 0 -366 782 3 783 3 57 23 c31 13 79 43 107 68 58 52 77 88 336 629 169 354 195 404 195 377 0 -12 62 -338 349 -1832 88 -459 205 -1067 259 -1350 79 -408 105 -526 126 -568 114 -225 414 -264 585 -77 86 95 67 -11 241 1355 22 168 66 521 100 785 34 264 78 617 100 785 21 168 59 467 84 665 26 198 73 571 106 828 33 258 60 476 60 484 0 50 15 -26 60 -303 27 -170 83 -518 125 -774 111 -686 197 -1218 265 -1640 34 -206 65 -391 70 -410 16 -58 67 -136 113 -174 75 -61 130 -81 227 -81 76 0 92 3 152 33 96 47 145 106 203 244 92 216 184 436 345 818 89 212 165 388 169 393 7 7 82 -210 445 -1278 81 -236 158 -458 172 -492 34 -82 88 -141 166 -184 60 -32 68 -34 163 -34 88 0 107 3 158 27 120 56 94 11 850 1445 l122 233 778 0 777 0 0 365 0 365 -889 0 c-559 0 -908 -4 -942 -10 -77 -15 -137 -47 -194 -105 -41 -42 -90 -127 -284 -495 -303 -575 -299 -567 -305 -550 -3 8 -65 188 -137 400 -72 212 -157 462 -189 555 -32 94 -118 345 -191 559 -143 419 -159 453 -245 518 -65 50 -133 72 -219 71 -103 0 -175 -31 -250 -107 l-60 -61 -224 -535 c-124 -294 -227 -539 -231 -543 -4 -4 -26 113 -49 260 -125 803 -645 3904 -663 3950 -27 74 -104 155 -181 192 -71 34 -138 43 -215 30z"/>
        //    </g>
        //</svg>
        
        // Create SVG element
        let Logo_svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        Logo_svg.setAttribute('class', 'File_SVG'); 
        Logo_svg.setAttribute('id', ''); 
        Logo_svg.setAttribute('width', '1.5rem');
        Logo_svg.setAttribute('height', '1.0rem');
        Logo_svg.setAttribute('viewBox', '0 0 1320.000000 920.000000');
        Logo_svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');

        let g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        g.setAttribute('transform', 'translate(0.000000,920.000000) scale(0.10000,-0.10000)');
        g.setAttribute('fill', '#000000');
        g.setAttribute('stroke', 'none');

        let path_1 = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        path_1.setAttribute('d', 'M120 4605 l0 -4445 6450 0 6450 0 0 4445 0 4445 -6450 0 -6450 0 0 -4445z m12710 0 l0 -4255 -6260 0 -6260 0 0 4255 0 4255 6260 0 6260 0 0 -4255z');
        
        let path_2 = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        path_2.setAttribute('d', 'M5562 8349 c-123 -21 -238 -111 -274 -215 -10 -27 -29 -137 -43 -244 -13 -107 -45 -357 -70 -555 -84 -662 -224 -1764 -310 -2437 -19 -150 -64 -507 -100 -793 -36 -286 -75 -587 -85 -670 -11 -82 -20 -156 -20 -164 0 -31 -12 0 -25 66 -13 59 -164 849 -411 2138 -84 438 -92 464 -179 544 -73 66 -146 94 -246 94 -116 0 -226 -56 -288 -146 -16 -23 -136 -267 -268 -542 -132 -275 -315 -659 -409 -852 l-169 -353 -632 0 -633 0 0 -365 0 -366 782 3 783 3 57 23 c31 13 79 43 107 68 58 52 77 88 336 629 169 354 195 404 195 377 0 -12 62 -338 349 -1832 88 -459 205 -1067 259 -1350 79 -408 105 -526 126 -568 114 -225 414 -264 585 -77 86 95 67 -11 241 1355 22 168 66 521 100 785 34 264 78 617 100 785 21 168 59 467 84 665 26 198 73 571 106 828 33 258 60 476 60 484 0 50 15 -26 60 -303 27 -170 83 -518 125 -774 111 -686 197 -1218 265 -1640 34 -206 65 -391 70 -410 16 -58 67 -136 113 -174 75 -61 130 -81 227 -81 76 0 92 3 152 33 96 47 145 106 203 244 92 216 184 436 345 818 89 212 165 388 169 393 7 7 82 -210 445 -1278 81 -236 158 -458 172 -492 34 -82 88 -141 166 -184 60 -32 68 -34 163 -34 88 0 107 3 158 27 120 56 94 11 850 1445 l122 233 778 0 777 0 0 365 0 365 -889 0 c-559 0 -908 -4 -942 -10 -77 -15 -137 -47 -194 -105 -41 -42 -90 -127 -284 -495 -303 -575 -299 -567 -305 -550 -3 8 -65 188 -137 400 -72 212 -157 462 -189 555 -32 94 -118 345 -191 559 -143 419 -159 453 -245 518 -65 50 -133 72 -219 71 -103 0 -175 -31 -250 -107 l-60 -61 -224 -535 c-124 -294 -227 -539 -231 -543 -4 -4 -26 113 -49 260 -125 803 -645 3904 -663 3950 -27 74 -104 155 -181 192 -71 34 -138 43 -215 30z');
        
        // Append elements
        g.appendChild(path_1);
        g.appendChild(path_2);
        Logo_svg.appendChild(g);
        
        return Logo_svg;
    }

}
 
function Add_To_Table(Channel) {

    // Appends single channel as a row entry to the Table
    // Each row consists of 11 columns as follow:
    // CheckBox - FileName - ChannelNumber - Duration (sec) - Sampling (Hz) - Type - Unit - Azimuth - Peak - Mean - RMS

    // Declare variables
    let table1, numOfRows, numOfCols, row, cell, select, Type_ID, opt, i, j;
    let Type_List, Unit_List, Unit_ID, Scale_ID, input;

    // Get the list of TypeList and UnitList
    Type_List = TypeUnitList().Types;
    Unit_List = UnitList(Channel.Unit).Units;

    // Number of row and columns
    table1    = document.getElementById('Table_Channels');
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
    // 1st column - (Checkbox) Insert an empty row first and populate it -------------------------------------------
    i = 0;
    cell = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Right');
    cell.setAttribute('onClick', "RePosition_Graph(" + (numOfRows-1).toString() + ")");
    cell.innerHTML = numOfRows;

    // 2nd column (FileName) - Insert an empty row first and populate it -------------------------------------------
    i = 1;
    cell = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Left');
    cell.setAttribute('onClick', "RePosition_Graph(" + (numOfRows-1).toString() + ")");
    cell.innerHTML = Channel.FileName;  

    // 3rd column (Channel Number) - Insert an empty row first and populate it -------------------------------------
    i = 2;
    cell = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Right');
    cell.setAttribute('onClick', "RePosition_Graph(" +  (numOfRows-1).toString() + ")");
    cell.innerHTML = Channel.ChNum;

    // 4th column (Duration of the records) - Insert an empty row first and populate it ----------------------------
    i = 3;
    cell = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Right');
    cell.innerHTML = Number(Channel.Duration).toFixed(3);
    cell.setAttribute('onClick', "RePosition_Graph(" +  (numOfRows-1).toString() + ")");

    // 5th column (Sampling Frequency in Hz) - Insert an empty row first and populate it ---------------------------
    i = 4;
    cell = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Right');
    cell.innerHTML = Number(Channel.FSamp).toFixed(3);
    cell.setAttribute('onClick', "RePosition_Graph(" +  (numOfRows-1).toString() + ")");

    // 6th colum (Channel Type) - Insert an empty row first and populate it ----------------------------------------
    i = 5;
    cell   = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Right');
    Type_ID = "Type_ID_" + Channel.Unique_ID;
    select  = document.createElement('select');
    select.setAttribute("id", Type_ID);
    select.setAttribute('class', 'form-select form-select-sm');
    select.setAttribute("onchange", "ChannelType_Update(this, true)");
    for (j=0; j<Type_List.length; j++) {
        opt = document.createElement("option");
        opt.value = Type_List[j];
        opt.text = Type_List[j];
        select.add(opt, null);
    }
    select.selectedIndex = Channel.Type;
    cell.appendChild(select);

    // 7th colum (Channel Unit) - Insert an empty row first and populate it ----------------------------------------
    i = 6;
    cell   = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Right');
    Unit_ID = "Unit_ID_" + Channel.Unique_ID;
    select  = document.createElement('select');
    select.setAttribute("id", Unit_ID);
    select.setAttribute('class', 'form-select form-select-sm');
    for (j = 0; j < Unit_List.length; j++) {
        opt = document.createElement("option");
        opt.value = Unit_List[j];
        opt.text = Unit_List[j];
        select.add(opt, null);
    }
    select.setAttribute("onchange", "ChannelType_Update(this, false)");
    select.selectedIndex = Channel.Unit - UnitList(Channel.Unit).UnitNum[0];
    cell.appendChild(select);

    // 8th column (Multiply-ScaleFactor) - Insert an empty row first and populate it ----------------------------------------
    i = 7;
    cell   = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Right');
    Scale_ID        = "Scale_ID_" + Channel.Unique_ID;
    input           = document.createElement('input');
    input.title     = "User specified scale factor";
    input.type      = "number";
    input.style     = "width: 100%;"
    input.value     = Channel.ScaleFactor;
    input.id        = Scale_ID;
    //input.onfocus   = "this.oldValue = this.value;";
    input.setAttribute('onchange', "ChannelScaleFactor_Update(this)" );
    cell.appendChild(input);

    // 8th column (Date& Time) - Insert an empty row first and populate it ----------------------------------------
    i = 8;
    cell = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Left');
    cell.setAttribute('onClick', "RePosition_Graph(" +  (numOfRows-1).toString() + ")");
    cell.innerHTML = Channel.DateTime.replace("T", " "); 

    // 9th colum (Orientation/Azimuth) - Insert an empty row first and populate it ---------------------------------
    i = 9;
    cell   = row.insertCell(i);
    cell.setAttribute('class', 'Channels_Table_Body_Td_Left');
    cell.innerHTML = Channel.Orientation;
    cell.setAttribute('onClick', "RePosition_Graph(" +  (numOfRows-1).toString() + ")");

    // Append to ChannelList
    Append_ChList( Channel );

    // Add to the File Tree
    Add_To_FileTree(Channel);
}

function Append_ChList(a) {
    let i, IA = infoArray(a);
    if (!IA.is1D && !IA.is2D) {  ChannelList.push(a);              }
    else if (IA.is1D)         {  a.map((v,i) => Append_ChList(v)); }
    else if (IA.is2D)         {  a.map((v,i) => Append_ChList(v)); }
}

async function ChannelScaleFactor_Update(el) {

    // Updates the ScaleFactor of a channel in the ChannelList[]
    // Updates the graph as required.

    // Declaration of variables
    let ChNum, Scale_Fac;

    // Obtain the channel Number
    ChNum = ChannelList_UniqueID( el.id.replace("Scale_ID_", "") );

    // Ge the user-specified scale factor
    Scale_Fac = Number(document.getElementById(el.id).value);

    // Updates the ScaleFactor
    ChannelList[ChNum].ScaleFactor = Scale_Fac;

    // Update Peak, Mean and RMS for this channel
    await ChannelList[ChNum].Statistics_update();

    // Update Graph
    await Plot_Update(ChNum);
}

async function ChannelType_Update(el, Opt) {

    // Declaration of variables
    let N=24, j, Indx, Indx2, Indx3, opt, temp, TempList, ChNum, Unit_ID;

    if (Opt) {
        // Channel type has been altered

        // Selected Channel-Type index
        Indx = document.getElementById(el.id).selectedIndex;

        // Get Channel Number 
        ChNum = ChannelList_UniqueID(el.id.replace("Type_ID_", ""));

        // Create the corresponding list of units
        TempList = UnitList(Indx, false).Units;

        // Update the Channel-Units based on the selected Channel-Type
        Unit_ID = el.id.replace("Type_ID", "Unit_ID");
        document.getElementById(Unit_ID).innerText = null;
        for (j = 0; j < TempList.length; j++) {
            opt = document.createElement("option");
            opt.value = TempList[j];
            opt.text = TempList[j];
            document.getElementById(Unit_ID).options[j] = opt;
        }

        // Select the unit number of the fist item in the unit-list
        Indx2 = UnitList(Indx, false).UnitNum[0];

    } else {
        // Channel unit has been altered
        
        // Selected Channel-Type index
        Indx = document.getElementById(el.id.replace("Unit_ID", "Type_ID")).selectedIndex;

        // Update the ChannelList[ChNum] accordingly.
        ChNum = ChannelList_UniqueID(el.id.replace("Unit_ID_", ""));

        // Index of the selected channel unit
        Indx3 = document.getElementById(el.id).selectedIndex; 

        // Obtain the unit number of the fist item in the unit-list
        Indx2 = UnitList(Indx, false).UnitNum[Indx3];
    }

    // Update the ChannelList[ChNum] accordingly.
    for (j=0; j<N; j++) {

        temp = TypeAndUnit(j, false);
        if ((temp.Type == Indx) && (temp.Unit == Indx2)) {

            temp = TypeAndUnit(temp.TypeAndUnit);

            ChannelList[ChNum].Type         = temp.Type;
            ChannelList[ChNum].TypeAndUnits = temp.TypeAndUnit;
            ChannelList[ChNum].TypeString   = temp.Type_String;
            ChannelList[ChNum].Unit         = temp.Unit;
            ChannelList[ChNum].UnitString   = temp.Unit_String;
            break;
        }
    }

    // Update Graph
    await Plot_Update(ChNum);

    // Update Right_Click Menu
    Update_RightClick_Menu();
}

async function Read_VIF(FileName, delta, dataview) {
    // Read the header information
    // These variables only exist within this block {...}. Does not exist out of this function.
    let nVersion, Year, Month, Day, Hour, Minute, Second, MilliSec, NanoSec;
    let delt, NumOfChannel, NumSampPerCh, DateTime, DeviceName, OptComment, AA, SF;
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


        //  return { 'IntervalTypeAndUnit': 1,   'Type': 0,   'Type_String': 'Time',   'Unit':  0,   'Unit_String': 's'   };

        temp                     = IntervalTypeAndUnit( 1 );   // Type And Unit number - refer to the list
        IntervalTypeAndUnits[Ch] = temp.IntervalTypeAndUnit;
        IntervalType[Ch]         = temp.Type;                  // (0-1)
        IntervalTypeString[Ch]   = temp.Type_String;           // (Time, Spectral Period, etc.)
        IntervalUnit[Ch]         = temp.Unit;                  // (0-2)
        IntervalUnitString[Ch]   = temp.Unit_String;           // (Second, DateTime, etc.)

        // Create new Channel Object and populate it
        let res = new Channel();

        res.FileName       = FileName;                 // FileName
        res.ChNum          = Ch;                       // Channel number
        res.NumSamples     = NumSamples[Ch];           // Number of digitized data points
        res.ScaleFactor    = 1;                        // User defined Scale Factor - Taken 1.0 if not specified
        res.Orientation    = Orientation[Ch];          // Orientation of channel
        res.DateTime       = DateTime;                 // Date&Time of the first sample
        res.Duration       = Duration[Ch];             // Total duration in seconds
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
        await res.Statistics_update();

        // Add to the Main Table and Tree View
        Add_To_Table( res );

        // Update the ProgressBar per channel (Unlike other type of channels)
        document.getElementById("ProgressBar_Label").innerHTML = "<b>Loading Your Files</b>";
        IncreaseProgressBar(delta/NumOfChannel, "ProgressBar_LoadData");
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

async function Read_V1(FileName, delta, dataview) {
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

    // Update the ProgressBar
    document.getElementById("ProgressBar_Label").innerHTML = "<b>Loading Your Files</b>";
    IncreaseProgressBar(delta, "ProgressBar_LoadData");

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

        DateTime[ChNum]     = IntNum[23] + "-" + ToString(IntNum[22], 2) + "-" + ToString(IntNum[21], 2) + "T" + ToString(IntNum[16], 2) + ":" + ToString(IntNum[17], 2) + ":" + ToString(IntNum[18], 2) + ".000";
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
            // 8f9.6 means 8 columns, each column includes 9 numbers with 6 decimal numbers
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
                Temp1 = GetRange(Temp, Ind2);
                Temp2 = GetRange(Temp, Ind1);
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

        temp                     = IntervalTypeAndUnit( 1 );   // Type And Unit number - refer to the list
        IntervalTypeAndUnits[ChNum] = temp.IntervalTypeAndUnit;
        IntervalType[ChNum]         = temp.Type;                  // (0-1)
        IntervalTypeString[ChNum]   = temp.Type_String;           // (Time, Spectral Period, etc.)
        IntervalUnit[ChNum]         = temp.Unit;                  // (0-2)
        IntervalUnitString[ChNum]   = temp.Unit_String;           // (Second, DateTime, etc.)


        // Create new Channel Object and populate it
        let res = new Channel();

        res.FileName       = FileName;             // FileName
        res.ChNum          = ChNum;                // Channel number
        res.NumSamples     = NumSamples[ChNum];    // Number of digitized data points
        res.ScaleFactor    = 1;                    // User defined Scale Factor - Taken 1.0 if not specified
        res.Orientation    = Orientation[ChNum];   // Orientation of channel
        res.DateTime       = DateTime[ChNum];      // Date&Time of the first sample
        res.Duration       = Duration[ChNum];      // Total duration in seconds
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
        await res.Statistics_update();

        // Add to the Main Table and Tree View
        Add_To_Table( res );
        
    }

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
        res.ChNum          = ChNum[i];                  // Channel number
        res.NumSamples     = NumSamples[i];             // Number of digitized data points
        res.ScaleFactor    = 1;                         // User defined Scale Factor - Taken 1.0 if not specified
        res.Orientation    = Orientation[i];            // Orientation of channel
        res.DateTime       = DateTime[i];               // Date&Time of the first sample
        res.Duration       = Duration[i];               // Total duration in seconds
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
        await res.Statistics_update();

        // Add to the Main Table and Tree View
        Add_To_Table( res );

    }

    // Update the ProgressBar
    document.getElementById("ProgressBar_Label").innerHTML = "<b>Loading Your Files</b>";
    IncreaseProgressBar(delta, "ProgressBar_LoadData");

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
        res.ChNum          = ChNum[i];                  // Channel number
        res.NumSamples     = NumSamples[i];             // Number of digitized data points
        res.ScaleFactor    = 1;                         // User defined Scale Factor - Taken 1.0 if not specified
        res.Orientation    = Orientation[i];            // Orientation of channel
        res.DateTime       = DateTime;                  // Date&Time of the first sample
        res.Duration       = Duration[i];               // Total duration in seconds
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
        await res.Statistics_update();

        // Add to the Main Table and Tree View
        Add_To_Table( res );
    }

    // Update the ProgressBar
    document.getElementById("ProgressBar_Label").innerHTML = "<b>Loading Your Files</b>";
    IncreaseProgressBar(delta, "ProgressBar_LoadData");

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
            d[bN] = Vector_Truncate(d[bN], NumberOfSamples[bN]);
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
    MergeChannels()

    // Update the ProgressBar
    document.getElementById("ProgressBar_Label").innerHTML = "<b>Loading Your Files</b>";
    IncreaseProgressBar(delta, "ProgressBar_LoadData");

    async function MergeChannels() {

        let i, j, jj, Ind, Ind1, temp, temp2, DateTime, DeltaSample, Azimuth=[], indices=[], data=[], time=[], NumOfChannel;
        let FSamp, delt;

        Azimuth      = [...new Set(ChannelIdentifier)];
        NumOfChannel = Azimuth.length;

        for (i=0; i < NumOfChannel; i++) {

            // Create new Channel Object and populate it
            let res = new Channel();

            indices = [];
            data    = [];
            time    = [];

            // this.Azimuth[i][0]; // General sampling rate and the response band of the instrument (H: High_Broad_Band, Sample Rate:  ≥ 80 to < 250, and Corner Period (sec): ≥ 10 sec.)
            // this.Azimuth[i][1]; // Instrument Code (N: Acceleration Sensor, H: High-gain seismometer, L: Low-gain seismometer)
            // this.Azimuth[i][2]; // Instrument orientation (Z, N, E: Traditional (Vertical, North-South, East-West)

            // Determine Channel Type and Unit
            if      (Azimuth[i][1] == 'H') { temp = TypeAndUnit(1);   }    // H: High-gain seismometer
            else if (Azimuth[i][1] == 'L') { temp = TypeAndUnit(1);   }    // L: Low-gain seismometer
            else if (Azimuth[i][1] == 'N') { temp = TypeAndUnit(1);   }    // N: Accelerometer
            else if (Azimuth[i][1] == 'A') { temp = TypeAndUnit(122); }    // A: Tilt
            else if (Azimuth[i][1] == 'D') { temp = TypeAndUnit(81);  }    // D: Pressure
            else if (Azimuth[i][1] == 'I') { temp = TypeAndUnit(71);  }    // I: Humidity
            else if (Azimuth[i][1] == 'K') { temp = TypeAndUnit(61);  }    // K: Temperature
            else if (Azimuth[i][1] == 'S') { temp = TypeAndUnit(31);  }    // S: Linear Strain
            else if (Azimuth[i][1] == 'V') { temp = TypeAndUnit(31);  }    // V: Volumetric Strain
            else if ((Azimuth[i][1] == 'W') && (Azimuth[i][2] == 'S')) { temp = TypeAndUnit(51); }    // W: Wind - Speed
            else if ((Azimuth[i][1] == 'W') && (Azimuth[i][2] == 'D')) { temp = TypeAndUnit(41); }    // D: Wind - Direction
            else { temp = TypeAndUnit(1); }     // N: Accelerometer

            // find indices of each channel
            ChannelIdentifier.filter( function (v, ii) { if (v == Azimuth[i] ) { indices.push(ii); } });

            // Start merging block-by-block
            Ind = indices[0];
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
            res.FileName       = FileName;
            res.DateTime       = DateTime;
            res.Orientation    = Azimuth[i];
            res.ScaleFactor    = 1;
            res.ChNum          = i;
            res.TypeAndUnits   = temp.TypeAndUnit;
            res.Type           = temp.Type;
            res.TypeString     = temp.Type_String;
            res.Unit           = temp.Unit;
            res.UnitString     = temp.Unit_String;
            res.IntervalTypeAndUnit = temp2.IntervalTypeAndUnit;    // Type And Unit number - refer to the list
            res.IntervalType        = temp2.Type;           // (0-1)
            res.IntervalTypeString  = temp2.Type_String;     // (Time, Spectral Period, etc.)
            res.IntervalUnit        = temp2.Unit;           // (0-2)
            res.IntervalUnitString  = temp2.Unit_String;     // (Second, DateTime, etc.)
            res.Duration       = time.at(-1);
            res.NumSamples     = data.length;
            res.FSamp          = FSamp;
            res.delt           = delt;
            res.data           = data;
            res.time           = time;

            // Calculate Statictics
            await res.Statistics_update();

            // Add to the Main Table and Tree View
            Add_To_Table( res );
        }
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

function FileExtension(FileName) {
    // Returns the extension of a file as string
    let Ind = FileName.indexOf('.');
    if (Ind == -1) { return ""; } else { return FileName.split('.').pop(); }
}

function IsFileUploaded(FileName) {
    // Return TRUE if the FileName is included in the ChannelList already; otherwise, returns FALSE.
    for (let j=0; j<ChannelList.length; j++) {
        if (ChannelList[j].FileName == FileName) { return true; }
    }
    return false;
}

//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
// File Loading Progress Bar
async function IncreaseProgressBar(delta, Div_ID) {

  // decleration of variables
  let ProgressBar, currentProgress;

  if (Div_ID == null) { Div_ID = "ProgressBar"; }

  ProgressBar = document.getElementById(Div_ID);

  // Get the current progress
  currentProgress = Number(ProgressBar.style.width.replace("%", ""));

  // Increase the progressBar by delta
  currentProgress += delta;   

  if (currentProgress.toFixed(0) < 100) { 
    // Update the progressBar
    ProgressBar.style.width = currentProgress.toString() + "%";
    ProgressBar.innerHTML = currentProgress.toFixed(0) + "%";
  }
  else {
    // Update the progressBar
    ProgressBar.style.width = '100%';
    ProgressBar.innerHTML = 'Done!';
  }
  return;

}
async function ResetProgressBar(Div_ID) {
  // decleration of variables
  let ProgressBar;

  if (Div_ID == null) { Div_ID = "ProgressBar"; }

  ProgressBar = document.getElementById(Div_ID);
  
  // Update the progressBar
  ProgressBar.style.width = '0%';
  ProgressBar.innerHTML = '';
  return;
}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
// Download Files to client browser 
function StartNewWorkBook() {
    let WorkBook = XLSX.utils.book_new();
    WorkBook.Props = {
        Title:   "SGM Seismic Data Analysis",
        Subject: "",
        Author:  "",
    };
    return WorkBook;
}
function StartNewWorkSheet() {
    // Start an empty worksheet
    return XLSX.utils.aoa_to_sheet([[]]);
}
function AddDataToWorkSheet(WorkSheet, Header, Pos1, Data, Pos2, Data2=null, Pos3=null) {
    // Add data to workSheet
    XLSX.utils.sheet_add_aoa(WorkSheet, Header, {origin: Pos1});
    XLSX.utils.sheet_add_aoa(WorkSheet, Data,   {origin: Pos2});
    if (Data2 != null) {
        XLSX.utils.sheet_add_aoa(WorkSheet, Data2,   {origin: Pos3});
    }
}
async function DonwloadExcel_LoadDataPage() {
  
    // Decleration of Variables 
    let i, FileName, WorkBook, WorkSheet, header, data, data2, res, delta, rN=1;
  
    FileName = "RawData.xlsx";
  
    // Retrun if no file is uploaded
    if (ChannelList.length == 0) { resolve(); return; }

    // Disable CALCULATE Button
    document.getElementById("Calculate_Button").disabled = true;
  
    // Reset the ProgressBar
    document.getElementById("ProgressBar_Label").innerHTML = "<b>Downloading XLSX File</b>";
    await ResetProgressBar("ProgressBar_LoadData");
    await sleep(50);
  
    // Increment for Progress Bar 
    delta = 100 / ChannelList.length;
  
    // Start new workbook
    WorkBook = StartNewWorkBook();

    // Create new workSheet for the channel
    WorkSheet = StartNewWorkSheet();

    // List_Of_Records
    header = ["No", "File Name", "Channel Number", "Duration (s)", "Sampling (Hz)", "Type", "Unit", "Scale Factor", "Date&Time", "Azimuth", "Peak", "Mean", "RMS"];
    XLSX.utils.sheet_add_aoa(WorkSheet, [header], {origin: "A1"});
    for (i=0; i<ChannelList.length; i++) {
        XLSX.utils.sheet_add_aoa(WorkSheet, [[i]],                                     {origin: {r:rN, c:0}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].FileName]],               {origin: {r:rN, c:1}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].ChNum]],                  {origin: {r:rN, c:2}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].Duration]],               {origin: {r:rN, c:3}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].FSamp]],                  {origin: {r:rN, c:4}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].TypeString]],             {origin: {r:rN, c:5}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].UnitString]],             {origin: {r:rN, c:6}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].ScaleFactor]],            {origin: {r:rN, c:7}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].DateTime]],               {origin: {r:rN, c:8}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].Orientation]],            {origin: {r:rN, c:9}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].Peak.toPrecision(4)]],    {origin: {r:rN, c:10}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].Mean.toPrecision(4)]],    {origin: {r:rN, c:11}});
        XLSX.utils.sheet_add_aoa(WorkSheet, [[ChannelList[i].RMS.toPrecision(4)]],     {origin: {r:rN, c:12}});
        rN++;
    }
    XLSX.utils.book_append_sheet(WorkBook, WorkSheet, "List_Of_Records");

    // Loop over each channel in ChannelList[]
    for (i = 0; i < ChannelList.length; i++) {
  
      // Get the user specified unit
      res = GraphUnits(i);
  
      // Create new workSheet for the channel
      WorkSheet = StartNewWorkSheet();
      
      if (PageNo == 0) {
        // Populate Time-RawData-FilteredData to WorkSheet
        header = ["Time (s)", "Raw Data (" + res.Gr_Data.Unit_String + ")"];
        data = Transpose([ChannelList[i].time, 
                          Mult(Mult(ChannelList[i].data, ChannelList[i].ScaleFactor), res.ScaleFactor)
                        ]);
        AddDataToWorkSheet(WorkSheet, [header], "A1", data, "A2");

        // Populate FileName, ChannelNumber, Recording Type, Orientation and Sampling Frequency
        header = [["File Name"],  ["Channel Number"], ["Recording Type"], ["Orientation"], ["Sampling Frequency (Hz)"], ["Peak (" + res.Gr_Data.Unit_String +")"], ["Mean (" + res.Gr_Data.Unit_String +")"], ["RMS (" + res.Gr_Data.Unit_String +")"]];   
        data = [[ChannelList[i].FileName], [ChannelList[i].ChNum], [ChannelList[i].TypeString], [ChannelList[i].Orientation], [ChannelList[i].FSamp], [ChannelList[i].Peak * res.ScaleFactor], [ChannelList[i].Mean * res.ScaleFactor], [ChannelList[i].RMS * res.ScaleFactor]];
        AddDataToWorkSheet(WorkSheet, header, "D1", data, "E1");

      }
      else if (PageNo == 1) {

        // Populate Time-RawData-FilteredData to WorkSheet
        header = ["Time (s)", "Raw Data (" + res.Gr_Data.Unit_String + ")", "Filtered Data (" + res.Gr_Data.Unit_String + ")"];
        data   = [ChannelList[i].time, 
                  Mult(Mult(ChannelList[i].data, ChannelList[i].ScaleFactor), res.ScaleFactor)];
        
        // Add Filtered Data if it is calculated
        if (ChannelList[i].Results.Filter.FilteredData != undefined) {
            data.push(Mult(ChannelList[i].Results.Filter.FilteredData, res.ScaleFactor));
        }
        AddDataToWorkSheet(WorkSheet, [header], "A1", Transpose(data), "A2");

        // Populate Filter Response to WorkSheet
        if (ChannelList[i].Results.Filter.f != undefined) {
            header = ["Freq (Hz)", "Magnitude", "Phase (rad)"];
            data = Transpose([ChannelList[i].Results.Filter.f, 
                              ChannelList[i].Results.Filter.Mag, 
                              ChannelList[i].Results.Filter.Angle,
                            ]);
            AddDataToWorkSheet(WorkSheet, [header], "E1", data, "E2");
        }
        

        // Populate FileName, ChannelNumber, Recording Type, Orientation and Sampling Frequency
        header = [["File Name"],  ["Channel Number"], ["Recording Type"], ["Orientation"], ["Sampling Frequency (Hz)"]];   
        data = [[ChannelList[i].FileName], [ChannelList[i].ChNum], [ChannelList[i].TypeString], [ChannelList[i].Orientation], [ChannelList[i].FSamp]];
        AddDataToWorkSheet(WorkSheet, header, "I1", data, "J1");

        // Populate Peak, Mean, RMS
        AddDataToWorkSheet(WorkSheet, [["Raw Data"]], "J7", [["Filtered Data"]], "K7");
        
        header = [["Peak (" + res.Gr_Data.Unit_String +")"],  ["Mean (" + res.Gr_Data.Unit_String +")"], ["RMS (" + res.Gr_Data.Unit_String +")"]];   
        data  = [[ChannelList[i].Peak * res.ScaleFactor], [ChannelList[i].Mean * res.ScaleFactor], [ChannelList[i].RMS * res.ScaleFactor]];
        data2 = [[ChannelList[i].Results.Filter.Peak * res.ScaleFactor], [ChannelList[i].Mean * res.ScaleFactor], [ChannelList[i].Results.Filter.RMS * res.ScaleFactor]];
        AddDataToWorkSheet(WorkSheet, header, "I8", data, "J8", data2, "K8");

        // Populate Filter Settings to WorkSheet
        header = ["Filter Settings"];
        data = "";
        AddDataToWorkSheet(WorkSheet, [header], "I12", data, "J12");

        header = [["Baseline Correction"],  ["Filter Name"], ["Filter Type"], ["Filter Order"], ["Cut-Off Frequency (Hz)"], ["Zero Phase"], ["Filter Stable"]];
        data = [[ChannelList[i].Results.Filter.BaselineCorrection_String],
                [ChannelList[i].Results.Filter.FilterName_String],
                [ChannelList[i].Results.Filter.FilterType_String],
                [ChannelList[i].Results.Filter.FilterOrder],
                [ChannelList[i].Results.Filter.FilterBand],
                [ChannelList[i].Results.Filter.ZeroPhase],
                [ChannelList[i].Results.Filter.IsStable],
            ];
        AddDataToWorkSheet(WorkSheet, header, "I13", data, "J13");

        // Populate Filter Coefficients to WorkSheet
        header = [["Filter a_Coefficients"], ["Filter b_Coefficients"]];
        data = [ChannelList[i].Results.Filter.a, ChannelList[i].Results.Filter.b];
        AddDataToWorkSheet(WorkSheet, header, "I20", data, "J20");

      }
      else if (PageNo == 2) {
        // Populate Time-RawData-FilteredData to WorkSheet
        header = ["Time (s)", "Raw Data (" + res.Gr_Data.Unit_String + ")", "Velocity (" + res.Gr_Data.Unit_String + " • s)", "Displacement (" + res.Gr_Data.Unit_String + " • s²)"];
        data = Transpose([ChannelList[i].time, 
                          Mult(Mult(ChannelList[i].data, ChannelList[i].ScaleFactor), res.ScaleFactor),
                          Mult(ChannelList[i].Results.Integral.Vel,  res.ScaleFactor),
                          Mult(ChannelList[i].Results.Integral.Disp, res.ScaleFactor),
                         ]);
        AddDataToWorkSheet(WorkSheet, [header], "A1", data, "A2");

        
      }
      
      // Add the workSheet to WorkBook
      XLSX.utils.book_append_sheet(WorkBook, WorkSheet, i.toString());
  
      // Update the progressBar
      IncreaseProgressBar(delta, "ProgressBar_LoadData");
  
      // Wait for screen update
      await sleep(50);
    }
  
    // Inform user on progress
    document.getElementById("ProgressBar_LoadData").innerHTML = "Writing XLSX File";
  
    // Wait for screen update
    await sleep(100);
  
    // Write workBook to XLSX file
    XLSX.writeFile(WorkBook, FileName, { compression: true });
  
    // Inform user on progress
    document.getElementById("ProgressBar_LoadData").innerHTML = "Done!";

    // Enable CALCULATE Button
    document.getElementById("Calculate_Button").disabled = false;
  
    // return the promise
    return;
    
}
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// Right Click Menu
function Right_Click_ALL() {

    // Decleration of all variables 
    let TreeView, Right_Click_Div, aa1, aa2, pp1, pp2;
    
    TreeView        = document.getElementById('TreeView');
    Right_Click_Div = document.getElementById('Right_Click_Div');

    // Show context menu on right-click
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

    // Hide context menu when clicking elsewhere
    document.addEventListener('click', (e) => {
        if (!Right_Click_Div.contains(e.target) && e.target !== TreeView) {
            Right_Click_Div.style.display = 'none';
        }
    });

}
function Update_RightClick_Menu() {

    // Decleration of variables 
    let List_Or, List_FSamp, List_Type;
    
    [List_Or, List_FSamp, List_Type] = Get_Orientation();
  
    // Clear the content of the UL elements
    document.getElementById('Right_Click_ul_First').innerHTML = '';
    document.getElementById('Right_Click_ul_Second').innerHTML = '';
    document.getElementById('Right_Click_ul_Third').innerHTML = '';
    document.getElementById('Right_Click_ul_Forth').innerHTML = '';
    
    // Hide the Right_Click_Menu
    document.getElementById("Right_Click_Div").style.display = "none";

    // Polulate Select and Unselect options
    Fist_Level("SelectAll_CheckBox",   "Select All",   "Right_Click_ul_First", "1");

    // Populate the selections based on the orientation of channel 
    if (List_Or.length > 0) {
      document.getElementById("Right_Click_Br_label_1").style.display = "block";
      document.getElementById("Right_Click_Br_label_2").style.display = "block";
      document.getElementById("Right_Click_Br_label_3").style.display = "block";
    }

    for (let i=0; i<List_Or.length; i++) {
      Fist_Level(List_Or[i],   List_Or[i],   "Right_Click_ul_Second", "2");
    } 

    for (let i=0; i<List_FSamp.length; i++) {
      Fist_Level("FSamp_"+Number(List_FSamp[i]).toFixed(2), Number(List_FSamp[i]).toFixed(2)+" Hz",   "Right_Click_ul_Third", "3");
    }

    for (let i=0; i<List_Type.length; i++) {
      Fist_Level("Type_" + List_Type[i], List_Type[i],   "Right_Click_ul_Forth", "4");
    }


    function Fist_Level(Input_ID, TEXT_Content, UL_ID, UserData) {

      // Decleratin of variables
      let UL_el, li_el, input, label
    
      UL_el = document.getElementById(UL_ID);

      li_el = document.createElement('li');
      input = document.createElement('input');
      label = document.createElement('label');

      input.setAttribute('className', "Parameter-check-input form-check-input");
      input.setAttribute('type',      "checkbox"         );
      input.setAttribute('id',         Input_ID          );
      input.setAttribute('userspec',   UserData);
      input.setAttribute('onclick',    "Right_Click(this)");

      label.textContent = TEXT_Content;
      label.setAttribute('for',         Input_ID            );
      label.setAttribute('class',       'RightClick_Label'  );

      li_el.setAttribute('className', "Right_Click_li"); 
      li_el.appendChild(input);
      li_el.appendChild(label);
      
      UL_el.appendChild(li_el);
    }

    function Get_Orientation() {
      // Decleratin of variables
      let List_1, List_2, List_3;

      List_1 = [];
      List_2 = [];
      List_3 = [];
      
      for (let i=0; i<ChannelList.length; i++) {
        List_1[i] = ChannelList[i].Orientation;
        List_2[i] = ChannelList[i].FSamp;
        List_3[i] = ChannelList[i].TypeString;
      }
      return [ [...new Set(List_1)], [...new Set(List_2)], [...new Set(List_3)] ];  // removes duplications
    }
}
function Right_Click(aa) {

    let i, K1, K2, Div_ID, value = aa.getAttribute('userspec');
    
    if (value === "1") {
      // SELECT ALL CHANNELS 
      
      UnSelect_By_UL("Right_Click_ul_Second", false);
      UnSelect_By_UL("Right_Click_ul_Third", false);
      UnSelect_By_UL("Right_Click_ul_Forth", false);

      if (aa.checked) {  
          UnSelect_By_UL("FileListTree", true);
          K1 = true;
          K2 = 'flex';
      }
      else {
          UnSelect_By_UL("FileListTree", false);
          K1 = false;
          K2 = 'none';
      }
      // Filter channels
      for (i=0; i<ChannelList.length; i++) { 
          ChannelList[i].Selected = K1; 
          Div_ID = "Div_ID_" + ChannelList[i].Unique_ID;
          document.getElementById(Div_ID).style.display = K2;
      }
    }
    else if (value === "2"){
      // SELECT BY AZIMUTH
      document.getElementById("SelectAll_CheckBox").checked   = false;
      UnSelect_By_UL("Right_Click_ul_Third", false);
      UnSelect_By_UL("Right_Click_ul_Forth", false);
      Filter_by_Azimuth();
    }
    else if (value === "3") {
      // SELECT BY SAMPLING FREQUENCY
      document.getElementById("SelectAll_CheckBox").checked   = false;
      UnSelect_By_UL("Right_Click_ul_Second", false);
      UnSelect_By_UL("Right_Click_ul_Forth", false);
      Filter_by_FSamp();
    }
    else if (value === "4") {
      // SELECT BY TYPE
      document.getElementById("SelectAll_CheckBox").checked   = false;
      UnSelect_By_UL("Right_Click_ul_Second", false);
      UnSelect_By_UL("Right_Click_ul_Third", false);
      Filter_by_Type();
    }
}
function UnSelect_By_UL(UL_ID, Karar) {
    let i, CB, List_li = document.getElementById(UL_ID).children;
    for (i=0; i<List_li.length; i++) {
        CB = List_li[i].querySelector('input[type="checkbox"]');
        if (CB != null) { CB.checked = Karar; }
    }
}
function Filter_by_Azimuth() {
    // Get unique list of Azimuth
    let i, CB, Div_ID, Az = [];
    let List_li = document.getElementById("Right_Click_ul_Second").children;
    for (i=0; i<List_li.length; i++) {
        CB = List_li[i].querySelector('input[type="checkbox"]');
        if (CB != null) { if (CB.checked) {Az.push(CB.id)} }
    }

    // Filter channels 
    for (i=0; i<ChannelList.length; i++) { 
        Div_ID = "Div_ID_" + ChannelList[i].Unique_ID;

        if (Az.includes(ChannelList[i].Orientation)) {
          ChannelList[i].Selected = true; 
          document.getElementById(Div_ID).style.display             = 'flex';
          document.getElementById(ChannelList[i].Unique_ID).checked = true;
        }
        else {
          ChannelList[i].Selected = false; 
          document.getElementById(Div_ID).style.display             = 'none';
          document.getElementById(ChannelList[i].Unique_ID).checked = false;
        }
    }
}
function Filter_by_FSamp() {
    // Get unique list of Sampling Frequency
    let i, CB, Div_ID, FS = [];
    let List_li = document.getElementById("Right_Click_ul_Third").children;
    for (i=0; i<List_li.length; i++) {
        CB = List_li[i].querySelector('input[type="checkbox"]');
        if (CB != null) { if (CB.checked) {FS.push(Number(CB.id.replace("FSamp_", "", 1)))} }
    }
    
    // Filter channels 
    for (i=0; i<ChannelList.length; i++) { 
        Div_ID = "Div_ID_" + ChannelList[i].Unique_ID;

        if (FS.includes(Number(ChannelList[i].FSamp.toFixed(2)))) {
          ChannelList[i].Selected = true; 
          document.getElementById(Div_ID).style.display             = 'flex';
          document.getElementById(ChannelList[i].Unique_ID).checked = true;
        }
        else {
          ChannelList[i].Selected = false; 
          document.getElementById(Div_ID).style.display             = 'none';
          document.getElementById(ChannelList[i].Unique_ID).checked = false;
        }
    }
}
function Filter_by_Type() {
    // Get unique list of Azimuth
    let i, CB, Div_ID, Ty = [];
    let List_li = document.getElementById("Right_Click_ul_Forth").children;
    for (i=0; i<List_li.length; i++) {
        CB = List_li[i].querySelector('input[type="checkbox"]');
        if (CB != null) { if (CB.checked) {Ty.push(CB.id.replace("Type_", "", 1))} }
    }
    
    // Filter channels 
    for (i=0; i<ChannelList.length; i++) { 
        Div_ID = "Div_ID_" + ChannelList[i].Unique_ID;

        if (Ty.includes(ChannelList[i].TypeString)) {
          ChannelList[i].Selected = true; 
          document.getElementById(Div_ID).style.display             = 'flex';
          document.getElementById(ChannelList[i].Unique_ID).checked = true;
        }
        else {
          ChannelList[i].Selected = false; 
          document.getElementById(Div_ID).style.display             = 'none';
          document.getElementById(ChannelList[i].Unique_ID).checked = false;
        }
    }
}