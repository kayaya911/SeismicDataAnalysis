// stricter parsing and error handling
"use strict";

//-------------------------------------------------------------------------------------------------------------
function TypeAndUnit(n, Opt) {

    // Returns measurement type and unit information
    //
    //  Parameters:
    //    n   : Measurment index (if Opt=true) or Unit number (if Opt=false)
    //    Opt : Boolen - If true, n is Measurment-index; if false, n is Unit-number
    //
    //
    //  Retruns:
    //    Object with TypeAndUnit, Type, Type_String, Unit, Unit_String properties
    //
    //
    // Supported measurement types and units are given below
    //  Measurment
    //  Index          Type                          Unit
    //  ----------     ----------------------        -----------------
    //  1              0   ( Acceleration   )        0   (g)
    //  2              0   ( Acceleration   )        1   (m/s²)
    //  3              0   ( Acceleration   )        2   (cm/s²)
    //  4              0   ( Acceleration   )        3   (mm/s²)
    //  11             1   ( Velocity       )        4   (gs)
    //  12             1   ( Velocity       )        5   (m/s)
    //  13             1   ( Velocity       )        6   (cm/s)
    //  14             1   ( Velocity       )        7   (mm/s)
    //  21             2   ( Displacement   )        8   (gs²)
    //  22             2   ( Displacement   )        9   (m)
    //  23             2   ( Displacement   )        10  (cm)
    //  24             2   ( Displacement   )        11  (mm)
    //  31             3   ( Strain         )        12  (Strain)
    //  41             4   ( Wind Direction )        13  (Degrees)
    //  51             5   ( Wind Speed     )        14  (Mile/h)
    //  52             5   ( Wind Speed     )        15  (km/h)
    //  53             5   ( Wind Speed     )        16  (m/s)
    //  54             5   ( Wind Speed     )        17  (Beaufort)
    //  61             6   ( Temperature    )        18  (°C)
    //  62             6   ( Temperature    )        19  (°F)
    //  71             7   ( Humidity       )        20  (%)
    //  81             8   ( Pore Pressure  )        21  (N/m²)
    //  121            9   ( Tilt           )        22  (Degrees)
    //  122            9   ( Tilt           )        23  (Volts)
    //
    // Author   : Dr. Yavuz Kaya, P.Eng.
    // Modified : 27.Jan.2026
    
    // Check default variables 
    if (Opt == null) { Opt = true; }


    if (Opt) {
        // TypeAndUnit Index (n) is given
        if      (n == 1)   { return { 'TypeAndUnit': 1  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 0 ,   'Unit_String': 'g',          'IntegrationUnits': IntegrationUnits(1  ) }; }   //  ( Acceleration   )  (g)
        else if (n == 2)   { return { 'TypeAndUnit': 2  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 1 ,   'Unit_String': 'm/s²',       'IntegrationUnits': IntegrationUnits(2  ) }; }   //  ( Acceleration   )  (m/s²)
        else if (n == 3)   { return { 'TypeAndUnit': 3  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 2 ,   'Unit_String': 'cm/s²',      'IntegrationUnits': IntegrationUnits(3  ) }; }   //  ( Acceleration   )  (cm/s²)
        else if (n == 4)   { return { 'TypeAndUnit': 4  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 3 ,   'Unit_String': 'mm/s²',      'IntegrationUnits': IntegrationUnits(4  ) }; }   //  ( Acceleration   )  (mm/s²)
        else if (n == 11)  { return { 'TypeAndUnit': 11 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 4 ,   'Unit_String': 'gs',         'IntegrationUnits': IntegrationUnits(11 ) }; }   //  ( Velocity       )  (gs)
        else if (n == 12)  { return { 'TypeAndUnit': 12 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 5 ,   'Unit_String': 'm/s',        'IntegrationUnits': IntegrationUnits(12 ) }; }   //  ( Velocity       )  (m/s)
        else if (n == 13)  { return { 'TypeAndUnit': 13 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 6 ,   'Unit_String': 'cm/s',       'IntegrationUnits': IntegrationUnits(13 ) }; }   //  ( Velocity       )  (cm/s)
        else if (n == 14)  { return { 'TypeAndUnit': 14 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 7 ,   'Unit_String': 'mm/s',       'IntegrationUnits': IntegrationUnits(14 ) }; }   //  ( Velocity       )  (mm/s)
        else if (n == 21)  { return { 'TypeAndUnit': 21 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 8 ,   'Unit_String': 'gs²',        'IntegrationUnits': IntegrationUnits(21 ) }; }   //  ( Displacement   )  (gs²)
        else if (n == 22)  { return { 'TypeAndUnit': 22 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 9 ,   'Unit_String': 'm',          'IntegrationUnits': IntegrationUnits(22 ) }; }   //  ( Displacement   )  (m)
        else if (n == 23)  { return { 'TypeAndUnit': 23 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 10,   'Unit_String': 'cm',         'IntegrationUnits': IntegrationUnits(23 ) }; }   //  ( Displacement   )  (cm)
        else if (n == 24)  { return { 'TypeAndUnit': 24 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 11,   'Unit_String': 'mm',         'IntegrationUnits': IntegrationUnits(24 ) }; }   //  ( Displacement   )  (mm)
        else if (n == 31)  { return { 'TypeAndUnit': 31 ,   'Type': 3,   'Type_String': 'Strain',         'Unit': 12,   'Unit_String': 'Strain',     'IntegrationUnits': IntegrationUnits(31 ) }; }   //  ( Strain         )  (Strain)
        else if (n == 41)  { return { 'TypeAndUnit': 41 ,   'Type': 4,   'Type_String': 'Wind Direction', 'Unit': 13,   'Unit_String': 'Degrees',    'IntegrationUnits': IntegrationUnits(41 ) }; }   //  ( Wind Direction )  (Degrees)
        else if (n == 51)  { return { 'TypeAndUnit': 51 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 14,   'Unit_String': 'Mile/h',     'IntegrationUnits': IntegrationUnits(51 ) }; }   //  ( Wind Speed     )  (Mile/h)
        else if (n == 52)  { return { 'TypeAndUnit': 52 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 15,   'Unit_String': 'km/h',       'IntegrationUnits': IntegrationUnits(52 ) }; }   //  ( Wind Speed     )  (km/h)
        else if (n == 53)  { return { 'TypeAndUnit': 53 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 16,   'Unit_String': 'm/s',        'IntegrationUnits': IntegrationUnits(53 ) }; }   //  ( Wind Speed     )  (m/s)
        else if (n == 54)  { return { 'TypeAndUnit': 54 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 17,   'Unit_String': 'Beaufort',   'IntegrationUnits': IntegrationUnits(54 ) }; }   //  ( Wind Speed     )  (Beaufort)
        else if (n == 61)  { return { 'TypeAndUnit': 61 ,   'Type': 6,   'Type_String': 'Temperature',    'Unit': 18,   'Unit_String': '°C',         'IntegrationUnits': IntegrationUnits(61 ) }; }   //  ( Temperature    )  (°C)
        else if (n == 62)  { return { 'TypeAndUnit': 62 ,   'Type': 6,   'Type_String': 'Temperature',    'Unit': 19,   'Unit_String': '°F',         'IntegrationUnits': IntegrationUnits(62 ) }; }   //  ( Temperature    )  (°F)
        else if (n == 71)  { return { 'TypeAndUnit': 71 ,   'Type': 7,   'Type_String': 'Humidity',       'Unit': 20,   'Unit_String': '%',          'IntegrationUnits': IntegrationUnits(71 ) }; }   //  ( Humidity       )  (%)
        else if (n == 81)  { return { 'TypeAndUnit': 81 ,   'Type': 8,   'Type_String': 'Pore Pressure',  'Unit': 21,   'Unit_String': 'N/m²',       'IntegrationUnits': IntegrationUnits(81 ) }; }   //  ( Pore Pressure  )  (N/m²)
        else if (n == 121) { return { 'TypeAndUnit': 121,   'Type': 9,   'Type_String': 'Tilt',           'Unit': 22,   'Unit_String': 'Degrees',    'IntegrationUnits': IntegrationUnits(121) }; }   //  ( Tilt           )  (Degrees)
        else if (n == 122) { return { 'TypeAndUnit': 122,   'Type': 9,   'Type_String': 'Tilt',           'Unit': 23,   'Unit_String': 'Volts',      'IntegrationUnits': IntegrationUnits(122) }; }   //  ( Tilt           )  (Volts)
        else { return TypeAndUnit(1, Opt); } // fallback 
    }
    else {
        // Unit number is given.
        if       ( n == 0  ) { return { 'TypeAndUnit': 1  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 0 ,   'Unit_String': 'g',         'IntegrationUnits': IntegrationUnits(1  )}; }   //  ( Acceleration   )  (g)
        else if  ( n == 1  ) { return { 'TypeAndUnit': 2  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 1 ,   'Unit_String': 'm/s²',      'IntegrationUnits': IntegrationUnits(2  )}; }   //  ( Acceleration   )  (m/s²)
        else if  ( n == 2  ) { return { 'TypeAndUnit': 3  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 2 ,   'Unit_String': 'cm/s²',     'IntegrationUnits': IntegrationUnits(3  )}; }   //  ( Acceleration   )  (cm/s²)
        else if  ( n == 3  ) { return { 'TypeAndUnit': 4  ,   'Type': 0,   'Type_String': 'Acceleration',   'Unit': 3 ,   'Unit_String': 'mm/s²',     'IntegrationUnits': IntegrationUnits(4  )}; }   //  ( Acceleration   )  (mm/s²)
        else if  ( n == 4  ) { return { 'TypeAndUnit': 11 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 4 ,   'Unit_String': 'gs',        'IntegrationUnits': IntegrationUnits(11 )}; }   //  ( Velocity       )  (gs)
        else if  ( n == 5  ) { return { 'TypeAndUnit': 12 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 5 ,   'Unit_String': 'm/s',       'IntegrationUnits': IntegrationUnits(12 )}; }   //  ( Velocity       )  (m/s)
        else if  ( n == 6  ) { return { 'TypeAndUnit': 13 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 6 ,   'Unit_String': 'cm/s',      'IntegrationUnits': IntegrationUnits(13 )}; }   //  ( Velocity       )  (cm/s)
        else if  ( n == 7  ) { return { 'TypeAndUnit': 14 ,   'Type': 1,   'Type_String': 'Velocity',       'Unit': 7 ,   'Unit_String': 'mm/s',      'IntegrationUnits': IntegrationUnits(14 )}; }   //  ( Velocity       )  (mm/s)
        else if  ( n == 8  ) { return { 'TypeAndUnit': 21 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 8 ,   'Unit_String': 'gs²',       'IntegrationUnits': IntegrationUnits(21 )}; }   //  ( Displacement   )  (gs²)
        else if  ( n == 9  ) { return { 'TypeAndUnit': 22 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 9 ,   'Unit_String': 'm',         'IntegrationUnits': IntegrationUnits(22 )}; }   //  ( Displacement   )  (m)
        else if  ( n == 10 ) { return { 'TypeAndUnit': 23 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 10,   'Unit_String': 'cm',        'IntegrationUnits': IntegrationUnits(23 )}; }   //  ( Displacement   )  (cm)
        else if  ( n == 11 ) { return { 'TypeAndUnit': 24 ,   'Type': 2,   'Type_String': 'Displacement',   'Unit': 11,   'Unit_String': 'mm',        'IntegrationUnits': IntegrationUnits(24 )}; }   //  ( Displacement   )  (mm)
        else if  ( n == 12 ) { return { 'TypeAndUnit': 31 ,   'Type': 3,   'Type_String': 'Strain',         'Unit': 12,   'Unit_String': 'Strain',    'IntegrationUnits': IntegrationUnits(31 )}; }   //  ( Strain         )  (Strain)
        else if  ( n == 13 ) { return { 'TypeAndUnit': 41 ,   'Type': 4,   'Type_String': 'Wind Direction', 'Unit': 13,   'Unit_String': 'Degrees',   'IntegrationUnits': IntegrationUnits(41 )}; }   //  ( Wind Direction )  (Degrees)
        else if  ( n == 14 ) { return { 'TypeAndUnit': 51 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 14,   'Unit_String': 'Mile/h',    'IntegrationUnits': IntegrationUnits(51 )}; }   //  ( Wind Speed     )  (Mile/h)
        else if  ( n == 15 ) { return { 'TypeAndUnit': 52 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 15,   'Unit_String': 'km/h',      'IntegrationUnits': IntegrationUnits(52 )}; }   //  ( Wind Speed     )  (km/h)
        else if  ( n == 16 ) { return { 'TypeAndUnit': 53 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 16,   'Unit_String': 'm/s',       'IntegrationUnits': IntegrationUnits(53 )}; }   //  ( Wind Speed     )  (m/s)
        else if  ( n == 17 ) { return { 'TypeAndUnit': 54 ,   'Type': 5,   'Type_String': 'Wind Speed',     'Unit': 17,   'Unit_String': 'Beaufort',  'IntegrationUnits': IntegrationUnits(54 )}; }   //  ( Wind Speed     )  (Beaufort)
        else if  ( n == 18 ) { return { 'TypeAndUnit': 61 ,   'Type': 6,   'Type_String': 'Temperature',    'Unit': 18,   'Unit_String': '°C',        'IntegrationUnits': IntegrationUnits(61 )}; }   //  ( Temperature    )  (°C)
        else if  ( n == 19 ) { return { 'TypeAndUnit': 62 ,   'Type': 6,   'Type_String': 'Temperature',    'Unit': 19,   'Unit_String': '°F',        'IntegrationUnits': IntegrationUnits(62 )}; }   //  ( Temperature    )  (°F)
        else if  ( n == 20 ) { return { 'TypeAndUnit': 71 ,   'Type': 7,   'Type_String': 'Humidity',       'Unit': 20,   'Unit_String': '%',         'IntegrationUnits': IntegrationUnits(71 )}; }   //  ( Humidity       )  (%)
        else if  ( n == 21 ) { return { 'TypeAndUnit': 81 ,   'Type': 8,   'Type_String': 'Pore Pressure',  'Unit': 21,   'Unit_String': 'N/m²',      'IntegrationUnits': IntegrationUnits(81 )}; }   //  ( Pore Pressure  )  (N/m²)
        else if  ( n == 22 ) { return { 'TypeAndUnit': 121,   'Type': 9,   'Type_String': 'Tilt',           'Unit': 22,   'Unit_String': 'Degrees',   'IntegrationUnits': IntegrationUnits(121)}; }   //  ( Tilt           )  (Degrees)
        else if  ( n == 23 ) { return { 'TypeAndUnit': 122,   'Type': 9,   'Type_String': 'Tilt',           'Unit': 23,   'Unit_String': 'Volts',     'IntegrationUnits': IntegrationUnits(122)}; }   //  ( Tilt           )  (Volts)
        else { return TypeAndUnit(0, Opt); } // fallback 
    }
}
//-------------------------------------------------------------------------------------------------------------
function IntegrationUnits(n) {

    // Returns the first and second integral units for a given measurement type
    //
    //  Parameters:
    //    n   : Measurment index
    //
    //  Returns:
    //    Object with FirstIntegral and SecondIntegral properties, each containing:
    //      - Type_String : String description of the integral type
    //      - Unit_String : Unit of measurement after integration
    //
    //
    // Integration chain for supported measurement types:
    //  Index (n)      Original Type/Unit                   First Integral                           Second Integral
    //  ----------     --------------------------------     ------------------------------------     ------------------------------------
    //  1              ( Acceleration   )  (g)              ( Velocity        )  (gs)                ( Displacement     )  (gs²)
    //  2              ( Acceleration   )  (m/s²)           ( Velocity        )  (m/s)               ( Displacement     )  (m)
    //  3              ( Acceleration   )  (cm/s²)          ( Velocity        )  (cm/s)              ( Displacement     )  (cm)
    //  4              ( Acceleration   )  (mm/s²)          ( Velocity        )  (mm/s)              ( Displacement     )  (mm)
    //  11             ( Velocity       )  (gs)             ( Displacement    )  (gs²)               (∫∫ Velocity       )  (gs³)
    //  12             ( Velocity       )  (m/s)            ( Displacement    )  (m)                 (∫∫ Velocity       )  (m·s)
    //  13             ( Velocity       )  (cm/s)           ( Displacement    )  (cm)                (∫∫ Velocity       )  (cm·s)
    //  14             ( Velocity       )  (mm/s)           ( Displacement    )  (mm)                (∫∫ Velocity       )  (mm·s)
    //  21             ( Displacement   )  (gs²)            (∫ Displacement   )  (gs³)               (∫∫ Displacement   )  (gs⁴)
    //  22             ( Displacement   )  (m)              (∫ Displacement   )  (m·s)               (∫∫ Displacement   )  (m·s²)
    //  23             ( Displacement   )  (cm)             (∫ Displacement   )  (cm·s)              (∫∫ Displacement   )  (cm·s²)
    //  24             ( Displacement   )  (mm)             (∫ Displacement   )  (mm·s)              (∫∫ Displacement   )  (mm·s²)
    //  31             ( Strain         )  (Strain)         (∫ Strain         )  (Strain·s)          (∫∫ Strain         )  (Strain·s²)
    //  41             ( Wind Direction )  (Degrees)        (∫ Wind Direction )  (Degrees·s)         (∫∫ Wind Direction )  (Degrees·s²)
    //  51             ( Wind Speed     )  (Mile/h)         (∫ Wind Speed     )  (Mile/h·s)          (∫∫ Wind Speed     )  (Mile/h·s²)
    //  52             ( Wind Speed     )  (km/h)           (∫ Wind Speed     )  (km/h·s)            (∫∫ Wind Speed     )  (km/h·s²)
    //  53             ( Wind Speed     )  (m/s)            (∫ Wind Speed     )  (m)                 (∫∫ Wind Speed     )  (m·s)
    //  54             ( Wind Speed     )  (Beaufort)       (∫ Wind Speed     )  (Beaufort·s)        (∫∫ Wind Speed     )  (Beaufort·s²)
    //  61             ( Temperature    )  (°C)             (∫ Temperature    )  (°C·s)              (∫∫ Temperature    )  (°C·s²)
    //  62             ( Temperature    )  (°F)             (∫ Temperature    )  (°F·s)              (∫∫ Temperature    )  (°F·s²)
    //  71             ( Humidity       )  (%)              (∫ Humidity       )  (%·s)               (∫∫ Humidity       )  (%·s²)
    //  81             ( Pore Pressure  )  (N/m²)           (∫ Pore Pressure  )  (N/m²·s)            (∫∫ Pore Pressure  )  (N/m²·s²)
    //  121            ( Tilt           )  (Degrees)        (∫ Tilt           )  (Degrees·s)         (∫∫ Tilt           )  (Degrees·s²)
    //  122            ( Tilt           )  (Volts)          (∫ Tilt           )  (Volts·s)           (∫∫ Tilt           )  (Volts·s²)
    //
    // Author   : Dr. Yavuz Kaya, P.Eng.
    // Modified : 27.Jan.2026

    if      (n == 1)   { return { FirstIntegral: { 'Type_String': 'Velocity',          'Unit_String': 'gs'         }, SecondIntegral: { 'Type_String': 'Displacement',       'Unit_String': 'gs²'         } } }
    else if (n == 2)   { return { FirstIntegral: { 'Type_String': 'Velocity',          'Unit_String': 'm/s'        }, SecondIntegral: { 'Type_String': 'Displacement',       'Unit_String': 'm'           } } }
    else if (n == 3)   { return { FirstIntegral: { 'Type_String': 'Velocity',          'Unit_String': 'cm/s'       }, SecondIntegral: { 'Type_String': 'Displacement',       'Unit_String': 'cm'          } } }
    else if (n == 4)   { return { FirstIntegral: { 'Type_String': 'Velocity',          'Unit_String': 'mm/s'       }, SecondIntegral: { 'Type_String': 'Displacement',       'Unit_String': 'mm'          } } }
    else if (n == 11)  { return { FirstIntegral: { 'Type_String': 'Displacement',      'Unit_String': 'gs²'        }, SecondIntegral: { 'Type_String': '∫∫ Velocity',        'Unit_String': 'gs³'         } } }
    else if (n == 12)  { return { FirstIntegral: { 'Type_String': 'Displacement',      'Unit_String': 'm'          }, SecondIntegral: { 'Type_String': '∫∫ Velocity',        'Unit_String': 'm·s'         } } }
    else if (n == 13)  { return { FirstIntegral: { 'Type_String': 'Displacement',      'Unit_String': 'cm'         }, SecondIntegral: { 'Type_String': '∫∫ Velocity',        'Unit_String': 'cm·s'        } } }
    else if (n == 14)  { return { FirstIntegral: { 'Type_String': 'Displacement',      'Unit_String': 'mm'         }, SecondIntegral: { 'Type_String': '∫∫ Velocity',        'Unit_String': 'mm·s'        } } }
    else if (n == 21)  { return { FirstIntegral: { 'Type_String': '∫ Displacement',    'Unit_String': 'gs³'        }, SecondIntegral: { 'Type_String': '∫∫ Displacement',    'Unit_String': 'gs⁴'         } } }
    else if (n == 22)  { return { FirstIntegral: { 'Type_String': '∫ Displacement',    'Unit_String': 'm·s'        }, SecondIntegral: { 'Type_String': '∫∫ Displacement',    'Unit_String': 'm·s²'        } } }
    else if (n == 23)  { return { FirstIntegral: { 'Type_String': '∫ Displacement',    'Unit_String': 'cm·s'       }, SecondIntegral: { 'Type_String': '∫∫ Displacement',    'Unit_String': 'cm·s²'       } } }
    else if (n == 24)  { return { FirstIntegral: { 'Type_String': '∫ Displacement',    'Unit_String': 'mm·s'       }, SecondIntegral: { 'Type_String': '∫∫ Displacement',    'Unit_String': 'mm·s²'       } } }
    else if (n == 31)  { return { FirstIntegral: { 'Type_String': '∫ Strain',          'Unit_String': 'Strain·s'   }, SecondIntegral: { 'Type_String': '∫∫ Strain',          'Unit_String': 'Strain·s²'   } } }
    else if (n == 41)  { return { FirstIntegral: { 'Type_String': '∫ Wind Direction',  'Unit_String': 'Degrees·s'  }, SecondIntegral: { 'Type_String': '∫∫ Wind Direction',  'Unit_String': 'Degrees·s²'  } } }
    else if (n == 51)  { return { FirstIntegral: { 'Type_String': '∫ Wind Speed',      'Unit_String': 'Mile/h·s'   }, SecondIntegral: { 'Type_String': '∫∫ Wind Speed',      'Unit_String': 'Mile/h·s²'   } } }
    else if (n == 52)  { return { FirstIntegral: { 'Type_String': '∫ Wind Speed',      'Unit_String': 'km/h·s'     }, SecondIntegral: { 'Type_String': '∫∫ Wind Speed',      'Unit_String': 'km/h·s²'     } } }
    else if (n == 53)  { return { FirstIntegral: { 'Type_String': '∫ Wind Speed',      'Unit_String': 'm'          }, SecondIntegral: { 'Type_String': '∫∫ Wind Speed',      'Unit_String': 'm·s'         } } }
    else if (n == 54)  { return { FirstIntegral: { 'Type_String': '∫ Wind Speed',      'Unit_String': 'Beaufort·s' }, SecondIntegral: { 'Type_String': '∫∫ Wind Speed',      'Unit_String': 'Beaufort·s²' } } }
    else if (n == 61)  { return { FirstIntegral: { 'Type_String': '∫ Temperature',     'Unit_String': '°C·s'       }, SecondIntegral: { 'Type_String': '∫∫ Temperature',     'Unit_String': '°C·s²'       } } }
    else if (n == 62)  { return { FirstIntegral: { 'Type_String': '∫ Temperature',     'Unit_String': '°F·s'       }, SecondIntegral: { 'Type_String': '∫∫ Temperature',     'Unit_String': '°F·s²'       } } }
    else if (n == 71)  { return { FirstIntegral: { 'Type_String': '∫ Humidity',        'Unit_String': '%·s'        }, SecondIntegral: { 'Type_String': '∫∫ Humidity',        'Unit_String': '%·s²'        } } }
    else if (n == 81)  { return { FirstIntegral: { 'Type_String': '∫ Pore Pressure',   'Unit_String': 'N/m²·s'     }, SecondIntegral: { 'Type_String': '∫∫ Pore Pressure',   'Unit_String': 'N/m²·s²'     } } }
    else if (n == 121) { return { FirstIntegral: { 'Type_String': '∫ Tilt',            'Unit_String': 'Degrees·s'  }, SecondIntegral: { 'Type_String': '∫∫ Tilt',            'Unit_String': 'Degrees·s²'  } } }
    else if (n == 122) { return { FirstIntegral: { 'Type_String': '∫ Tilt',            'Unit_String': 'Volts·s'    }, SecondIntegral: { 'Type_String': '∫∫ Tilt',            'Unit_String': 'Volts·s²'    } } }

}
//-------------------------------------------------------------------------------------------------------------
function List_TypeUnit() {

    // Returns lists of all unique measurement types and units supported by the system
    //
    //  Returns:
    //    Object with Types and Units properties:
    //      - Types : Array of unique measurement type strings   (e.g., ['Acceleration', 'Velocity', ...])
    //      - Units : Array of unique unit strings               (e.g., ['g', 'm/s²', 'cm/s²', ...])
    //
    // Author   : Dr. Yavuz Kaya, P.Eng.
    // Modified : 27.Jan.2026

    // N : Number of measurement types
    let N=24, i, temp, Type = [], Unit =[];
    for (i=0; i<N; i++) { 
        temp = TypeAndUnit(i, false); 
        Type.push(temp.Type_String); 
        Unit.push(temp.Unit_String); 
    }
    
    // Remove duplicates in the Type and Unit arrays
    Type = [...new Set(Type)];
    Unit = [...new Set(Unit)];
    
    return { Types: Type, Units: Unit }
}
//-------------------------------------------------------------------------------------------------------------
function List_Units(n, Opt) {

    // Returns a list of all available units for a specific measurement type
    //
    //  Parameters:
    //    n   : Unit-Number (if Opt=true) or Type-Number (if Opt=false)
    //    Opt : Boolean (True)   n is unit-number
    //          Boolean (False)  n is Type-index
    //
    //  Returns:
    //    Object with Units and UnitNum properties:
    //      - Units    : Array of unit strings for the specified type  (e.g., ['g', 'm/s²', 'cm/s²', 'mm/s²'])
    //      - UnitNum  : Array of corresponding unit-numbers           (e.g., [0, 1, 2, 3])
    //
    // Author   : Dr. Yavuz Kaya, P.Eng.
    // Modified : 27.Jan.2026

    // Check default variables
    if (Opt == null) { Opt = true; }

    // N : Number of measurement types
    let N=24, i, temp, temp2, Unit=[], UnitNum=[];

    if (Opt) {
        // Return a list of all Units of a specific unit-number
        // n : Unit number
        temp = TypeAndUnit(n, false).Type_String;

        // Loop over unit numbers from 0 to 23
        for (i=0; i<N; i++) { 
            temp2 = TypeAndUnit(i, false); 
            
            if (temp == temp2.Type_String) { 
                Unit.push(temp2.Unit_String); 
                UnitNum.push(temp2.Unit); 
            } 

        }
        return { Units: Unit, UnitNum: UnitNum }

    } else {
        // Return a list of all Units of a specific Type-Index
        // n : Type Index

        // Loop over unit numbers from 0 to 23
        for (i=0; i<N; i++) { 
            temp2 = TypeAndUnit(i, false); 
            
            if (n == temp2.Type) { 
                Unit.push(temp2.Unit_String); 
                UnitNum.push(temp2.Unit); 
            } 
        }
        return { Units: Unit, UnitNum: UnitNum }
    }

}
//-------------------------------------------------------------------------------------------------------------
function Convert_Units_Data(Data, In, Out, Opt) {

    // Converts data array from one unit to another within the same measurement type
    //
    //  Parameters:
    //    Data : Array of numerical values to be converted
    //    In   : Input Measurment index  (if Opt=true) or 
    //           Unit number             (if Opt=false)
    //    Out  : Output Measurment index (if Opt=true) or 
    //           Unit number             (if Opt=false)
    //    Opt  : Boolean - if true,  In and Out are Measurment indices; 
    //                     if false, In and Out are Unit numbers
    //
    //  Returns:
    //    Object with Data and Unit properties:
    //      - Data  : Array of converted values in the output unit
    //      - Unit : String representation of the output unit       (e.g., 'g', 'm/s²', '°C', etc.)
    //
    //  Conversion Method:
    //    Most conversions: OutputValue = (InputValue × Scale) + Offset
    //    Special cases:
    //      - Temperature (°C ↔ °F): Uses scale and offset for proper conversion
    //      - Beaufort scale: Uses lookup tables for non-linear conversion
    //      - Same unit: Returns data unchanged
    //
    // Author   : Dr. Yavuz Kaya, P.Eng.
    // Modified : 27.Jan.2026

    let iU, oU, Scale=0, Offset=0, Unit='';

    // Check default variables
    if (Opt == null) { Opt = true; }

    iU = TypeAndUnit(In,  Opt);
    oU = TypeAndUnit(Out, Opt);

    if (iU.Unit == 0) {
        // Acceleration, g
        if      (oU.Unit == 0) { Scale = 1;          Offset = 0;  Unit = 'g';                                 }  // g
        else if (oU.Unit == 1) { Scale = 9.81;       Offset = 0;  Unit = 'm/s²';                              }  // m/s²
        else if (oU.Unit == 2) { Scale = 981;        Offset = 0;  Unit = 'cm/s²';                             }  // cm/s²
        else if (oU.Unit == 3) { Scale = 9810;       Offset = 0;  Unit = 'mm/s²';                             }  // mm/s²
    }                                                                                                          
    else if (iU.Unit == 1) {                                                                                   
        // Acceleration, m/s²                                                                                  
        if      (oU.Unit == 0) { Scale = 1/9.81;     Offset = 0;  Unit = 'g';                                 }  // g
        else if (oU.Unit == 1) { Scale = 1;          Offset = 0;  Unit = 'm/s²';                              }  // m/s²
        else if (oU.Unit == 2) { Scale = 100;        Offset = 0;  Unit = 'cm/s²';                             }  // cm/s²
        else if (oU.Unit == 3) { Scale = 1000;       Offset = 0;  Unit = 'mm/s²';                             }  // mm/s²
    }                                                                                                          
    else if (iU.Unit == 2) {                                                                                   
        // Acceleration, cm/s²                                                                                 
        if      (oU.Unit == 0) { Scale = 1/981;      Offset = 0;  Unit = 'g';                                 }  // g
        else if (oU.Unit == 1) { Scale = 0.01;       Offset = 0;  Unit = 'm/s²';                              }  // m/s²
        else if (oU.Unit == 2) { Scale = 1;          Offset = 0;  Unit = 'cm/s²';                             }  // cm/s²
        else if (oU.Unit == 3) { Scale = 10;         Offset = 0;  Unit = 'mm/s²';                             }  // mm/s²
    }                                                                                                          
    else if (iU.Unit == 3) {                                                                                   
        // Acceleration, mm/s²                                                                                 
        if      (oU.Unit == 0) { Scale = 1/9810;     Offset = 0;  Unit = 'g';                                 }  // g
        else if (oU.Unit == 1) { Scale = 0.001;      Offset = 0;  Unit = 'm/s²';                              }  // m/s²
        else if (oU.Unit == 2) { Scale = 0.1;        Offset = 0;  Unit = 'cm/s²';                             }  // cm/s²
        else if (oU.Unit == 3) { Scale = 1;          Offset = 0;  Unit = 'mm/s²';                             }  // mm/s²
    }                                                                                                          
    else if (iU.Unit == 4) {                                                                                   
        // Velocity, gs                                                                                        
        if      (oU.Unit == 4) { Scale = 1;          Offset = 0;  Unit = 'gs';                                }  // gs
        else if (oU.Unit == 5) { Scale = 9.81;       Offset = 0;  Unit = 'm/s';                               }  // m/s
        else if (oU.Unit == 6) { Scale = 981;        Offset = 0;  Unit = 'cm/s';                              }  // cm/s
        else if (oU.Unit == 7) { Scale = 9810;       Offset = 0;  Unit = 'mm/s';                              }  // mm/s
    }                                                                                                          
    else if (iU.Unit == 5) {                                                                                   
        // Velocity, m/s                                                                                       
        if      (oU.Unit == 4) { Scale = 1/9.81;     Offset = 0;  Unit = 'gs';                                }  // gs
        else if (oU.Unit == 5) { Scale = 1;          Offset = 0;  Unit = 'm/s';                               }  // m/s
        else if (oU.Unit == 6) { Scale = 100;        Offset = 0;  Unit = 'cm/s';                              }  // cm/s
        else if (oU.Unit == 7) { Scale = 1000;       Offset = 0;  Unit = 'mm/s';                              }  // mm/s
    }                                                                                                          
    else if (iU.Unit == 6) {                                                                                   
        // Velocity, cm/s                                                                                      
        if      (oU.Unit == 4) { Scale = 1/981;      Offset = 0;  Unit = 'gs';                                }  // gs
        else if (oU.Unit == 5) { Scale = 0.01;       Offset = 0;  Unit = 'm/s';                               }  // m/s
        else if (oU.Unit == 6) { Scale = 1;          Offset = 0;  Unit = 'cm/s';                              }  // cm/s
        else if (oU.Unit == 7) { Scale = 10;         Offset = 0;  Unit = 'mm/s';                              }  // mm/s
    }                                                                                                          
    else if (iU.Unit == 7) {                                                                                   
        // Velocity, mm/s                                                                                      
        if      (oU.Unit == 4) { Scale = 1/9810;     Offset = 0;  Unit = 'gs';                                }  // gs
        else if (oU.Unit == 5) { Scale = 0.001;      Offset = 0;  Unit = 'm/s';                               }  // m/s
        else if (oU.Unit == 6) { Scale = 0.1;        Offset = 0;  Unit = 'cm/s';                              }  // cm/s
        else if (oU.Unit == 7) { Scale = 1;          Offset = 0;  Unit = 'mm/s';                              }  // mm/s
    }                                                                                                          
    else if (iU.Unit == 8) {                                                                                   
        // Displacement, gs²                                                                                   
        if      (oU.Unit == 8)  { Scale = 1;         Offset = 0;  Unit = 'gs²';                               }  // gs²
        else if (oU.Unit == 9)  { Scale = 9.81;      Offset = 0;  Unit = 'm';                                 }  // m
        else if (oU.Unit == 10) { Scale = 981;       Offset = 0;  Unit = 'cm';                                }  // cm
        else if (oU.Unit == 11) { Scale = 9810;      Offset = 0;  Unit = 'mm';                                }  // mm
    }                                                                                                          
    else if (iU.Unit == 9) {                                                                                   
        // Displacement, m                                                                                     
        if      (oU.Unit == 8)  { Scale = 1/9.81;    Offset = 0;  Unit = 'gs²';                               }  // gs²
        else if (oU.Unit == 9)  { Scale = 1;         Offset = 0;  Unit = 'm';                                 }  // m
        else if (oU.Unit == 10) { Scale = 100;       Offset = 0;  Unit = 'cm';                                }  // cm
        else if (oU.Unit == 11) { Scale = 1000;      Offset = 0;  Unit = 'mm';                                }  // mm
    }                                                                                                          
    else if (iU.Unit == 10) {                                                                                  
        // Displacement, cm                                                                                    
        if      (oU.Unit == 8)  { Scale = 1/981;     Offset = 0;  Unit = 'gs²';                               }  // gs²
        else if (oU.Unit == 9)  { Scale = 0.01;      Offset = 0;  Unit = 'm';                                 }  // m
        else if (oU.Unit == 10) { Scale = 1;         Offset = 0;  Unit = 'cm';                                }  // cm
        else if (oU.Unit == 11) { Scale = 10;        Offset = 0;  Unit = 'mm';                                }  // mm
    }                                                                                                          
    else if (iU.Unit == 11) {                                                                                  
        // Displacement, mm                                                                                    
        if      (oU.Unit == 8)  { Scale = 1/9810;    Offset = 0;  Unit = 'gs²';                               }  // gs²
        else if (oU.Unit == 9)  { Scale = 0.001;     Offset = 0;  Unit = 'm';                                 }  // m
        else if (oU.Unit == 10) { Scale = 0.1;       Offset = 0;  Unit = 'cm';                                }  // cm
        else if (oU.Unit == 11) { Scale = 1;         Offset = 0;  Unit = 'mm';                                }  // mm
    }                                                                                                          
    else if (iU.Unit == 12) {                                                                                  
        // Strain, strain                                                                                      
        if      (oU.Unit == 12) { return { Data: Data, Unit: 'Strain' };                                      }  // strain
    }                                                                                                          
    else if (iU.Unit == 13) {                                                                                  
        // Wind direction, degree                                                                              
        if      (oU.Unit == 13) { return { Data: Data, Unit: 'Degrees' };                                     }  // degrees
    }                                                
    else if (iU.Unit == 14) {                        
        // Wind Speed, Mile/h                        
        if      (oU.Unit == 14) { Scale = 1;         Offset = 0;                          Unit = 'Mile/h';    }  // Mile/h
        else if (oU.Unit == 15) { Scale = 1.609344;  Offset = 0;                          Unit = 'km/h';      }  // km/h
        else if (oU.Unit == 16) { Scale = 44.704;    Offset = 0;                          Unit = 'm/s';       }  // m/s
        else if (oU.Unit == 17) { return { Data: kmhToBeaufort(Multiply(Data, 1.609344)), Unit: 'Beaufort'};  }  // Beaufort
    }                                                
    else if (iU.Unit == 15) {                        
        // Wind Speed, km/h                          
        if      (oU.Unit == 14) { Scale = 0.621371;  Offset = 0;      Unit = 'Mile/h';                        }  // Mile/h
        else if (oU.Unit == 15) { Scale = 1.0;       Offset = 0;      Unit = 'km/h';                          }  // km/h
        else if (oU.Unit == 16) { Scale = 27.7778;   Offset = 0;      Unit = 'm/s';                           }  // m/s
        else if (oU.Unit == 17) { return { Data: kmhToBeaufort(Data), Unit: 'Beaufort'};                      }  // Beaufort
    }                                                                              
    else if (iU.Unit == 16) {                                                      
        // Wind Speed, m/s                                                         
        if      (oU.Unit == 14) { Scale = 0.0223694; Offset = 0;                     Unit = 'Mile/h';         }  // Mile/h
        else if (oU.Unit == 15) { Scale = 0.036;     Offset = 0;                     Unit = 'km/h';           }  // km/h
        else if (oU.Unit == 16) { Scale = 1.0;       Offset = 0;                     Unit = 'm/s';            }  // m/s
        else if (oU.Unit == 17) { return { Data: kmhToBeaufort(Multiply(Data, 3.6)), Unit: 'Beaufort' };      }  // Beaufort
    }                                                
    else if (iU.Unit == 17) {                        
        // Wind Speed, Beaufort                          
        if      (oU.Unit == 14) { return { Data: Multiply(beaufortToKmh(Data), 0.621371), Unit: 'Mile/h'   }; }  // Mile/h
        else if (oU.Unit == 15) { return { Data: beaufortToKmh(Data),                     Unit: 'km/h'     }; }  // km/h
        else if (oU.Unit == 16) { return { Data: Multiply(beaufortToKmh(Data), (1/3.6)),  Unit: 'm/s'      }; }  // m/s
        else if (oU.Unit == 17) { return { Data: Data,                                    Unit: 'Beaufort' }; }  // Beaufort
    }                                                
    else if (iU.Unit == 18) {                        
        // temperature, °C                           
        if      (oU.Unit == 18) { return { Data: Data,                         Unit : '°C' };                 }  // °C
        else if (oU.Unit == 19) { return { Data: Add(Multiply(Data, 9/5), 32), Unit : '°F' };                 }  // °F
    }                                                
    else if (iU.Unit == 19) {                        
        // temperature, °F                           
        if      (oU.Unit == 18) { return { Data: Add(Multiply(Data, 5/9), -160/9), Unit : '°C' };             }  // °C
        else if (oU.Unit == 19) { return { Data: Data,                             Unit : '°F' };             }  // °F
    }                                                
    else if (iU.Unit == 20) {                        
        // Humidity, %                               
        if      (oU.Unit == 20) { return { Data: Data,                             Unit : '%' };              }  // %
    }                                                                                   
    else if (iU.Unit == 21) {                                                           
        // Pore Pressure, N/m²                                                          
        if      (oU.Unit == 21) { return { Data: Data,                             Unit : 'N/m²' };           }  // N/m²
    }                                                
    else if (iU.Unit == 22) {                        
        // Tilt, Degrees                             
        if      (oU.Unit == 22) { Scale = 1.0;       Offset = 0;                                               }  // Degrees
        else if (oU.Unit == 23) { Scale = 1.0;       Offset = 0;                                               }  // Volts
    }                                                
    else if (iU.Unit == 23) {                        
        // Tilt, Volts                               
        if      (oU.Unit == 22) { Scale = 1.0;       Offset = 0;                                               }  // Degrees
        else if (oU.Unit == 23) { Scale = 1.0;       Offset = 0;                                               }  // Volts
    }
    
    // return results
    return { 
        Data    : Add(Multiply(Data, Scale), Offset), 
        Unit    : Unit 
    };

    // Helper functions
    function kmhToBeaufort(kmhArray) {
        // Converts wind speed from km/h to Beaufort scale (0-12)
        // Uses standard Beaufort scale thresholds
        return kmhArray.map(kmh => {
            if (kmh < 1)   return 0;
            if (kmh <= 5)  return 1;
            if (kmh <= 11) return 2;
            if (kmh <= 19) return 3;
            if (kmh <= 28) return 4;
            if (kmh <= 38) return 5;
            if (kmh <= 49) return 6;
            if (kmh <= 61) return 7;
            if (kmh <= 74) return 8;
            if (kmh <= 88) return 9;
            if (kmh <= 102) return 10;
            if (kmh <= 117) return 11;
            return 12;
        });
    }
    function beaufortToKmh(beaufortArray) {
        // Converts Beaufort scale (0-12) to wind speed in km/h
        // Uses midpoint values for each Beaufort number
        const beaufortToKmhMap = {
            0: 0,      // Calm: < 1 km/h
            1: 3,      // Light air: 1-5 km/h
            2: 8,      // Light breeze: 6-11 km/h
            3: 15,     // Gentle breeze: 12-19 km/h
            4: 24,     // Moderate breeze: 20-28 km/h
            5: 33,     // Fresh breeze: 29-38 km/h
            6: 44,     // Strong breeze: 39-49 km/h
            7: 55,     // Near gale: 50-61 km/h
            8: 68,     // Gale: 62-74 km/h
            9: 81,     // Strong gale: 75-88 km/h
            10: 95,    // Storm: 89-102 km/h
            11: 110,   // Violent storm: 103-117 km/h
            12: 120    // Hurricane: ≥ 118 km/h
        };
        return beaufortArray.map(beaufort => beaufortToKmhMap[beaufort]);
    }
}
//-------------------------------------------------------------------------------------------------------------
function Convert_Data_To_Graph_Unit(Data, ChNum) {

    let Or_Data, LU, Ind, temp, temp_Stats, Stats;
    let Indx_Acc, Indx_Vel, Indx_Disp;

    // Original data using Measurment-unit
    Or_Data = TypeAndUnit(ChannelList[ChNum].TypeAndUnits); 

    // List of Units of the original data using Type-Number
    LU = List_Units(ChannelList[ChNum].Type, false);

    // Indexnumber of User-sepecified Unit on the Plotly Graph
    Ind = document.getElementById("Unit_Plot_ID_" + ChannelList[ChNum].Unique_ID).selectedIndex;

    // Convert Data to user-specified unit on Plotly Graph  ---  Using Unit-Number
    temp = Convert_Units_Data(Data,   Or_Data.Unit,   LU.UnitNum[Ind],   false);

    // Original statistical values are already scaled by Scale Factor of the channel 
    // Therefore, we just need to conver the units.
    if        (PageNo == 0) {
        
        temp_Stats = [ChannelList[ChNum].Peak, ChannelList[ChNum].Mean, ChannelList[ChNum].RMS];

        Stats = Convert_Units_Data(temp_Stats, Or_Data.Unit, LU.UnitNum[Ind], false);
    
    } else if (PageNo == 1) {
        if (ChannelList[ChNum].Results.Filter.IsAnalysisCompleted) {
            temp_Stats = [ChannelList[ChNum].Results.Filter.Peak, ChannelList[ChNum].Results.Filter.Mean, ChannelList[ChNum].Results.Filter.RMS];
            Stats = Convert_Units_Data(temp_Stats, Or_Data.Unit, LU.UnitNum[Ind], false);
        } else {
            Stats = { Data : ['','',''] };
        }
    } else if (PageNo == 2) {
        if (ChannelList[ChNum].Results.Integral.IsAnalysisCompleted) {
            // Find out which results to plot (Acc, Vel or Disp)
            Indx_Acc  = document.getElementById("Int_Acceleration").checked;
            Indx_Vel  = document.getElementById("Int_Velocity").checked;
            Indx_Disp = document.getElementById("Int_Displacement").checked;

            if      (Indx_Acc)  { temp_Stats = [ChannelList[ChNum].Peak,                       ChannelList[ChNum].Mean,                       ChannelList[ChNum].RMS,                       ChannelList[ChNum].Residual                       ];  }
            else if (Indx_Vel)  { temp_Stats = [ChannelList[ChNum].Results.Integral.Peak_Vel,  ChannelList[ChNum].Results.Integral.Mean_Vel,  ChannelList[ChNum].Results.Integral.RMS_Vel,  ChannelList[ChNum].Results.Integral.Residual_Vel  ];  }
            else if (Indx_Disp) { temp_Stats = [ChannelList[ChNum].Results.Integral.Peak_Disp, ChannelList[ChNum].Results.Integral.Mean_Disp, ChannelList[ChNum].Results.Integral.RMS_Disp, ChannelList[ChNum].Results.Integral.Residual_Disp ];  }

            Stats = Convert_Units_Data(temp_Stats, Or_Data.Unit, LU.UnitNum[Ind], false);

            temp.Residual    = Stats.Data[3];

        } else {
            Stats = { Data : ['','','',''] };
        }
    } else if (PageNo == 3) {
        if (ChannelList[ChNum].Results.SDOF.IsAnalysisCompleted) {

            Stats = { Data : [1,1,1] };
        } else {
            Stats = { Data : ['','',''] };
        }
    }

    //  temp.Data
    //  temp.Unit
    temp.Vel_Units   = LU.Units.map( item => item + ' • s');
    temp.Disp_Units  = LU.Units.map( item => item + ' • s²');
    
    temp.yTitle      = '<b>' + Or_Data.Type_String + ' [' + temp.Unit  + ']<b>';
    temp.yTitle_Vel  = '<b>' + Or_Data.IntegrationUnits.FirstIntegral.Type_String  + ' [' + temp.Vel_Units[Ind]   + ']<b>';
    temp.yTitle_Disp = '<b>' + Or_Data.IntegrationUnits.SecondIntegral.Type_String + ' [' + temp.Disp_Units[Ind]  + ']<b>';
    temp.yTitle_FFT  = '<b>Magnitude<b>';
    temp.y2Title     = '<b>Phase<b>';
    
    temp.Or_Units    = Select_Element(LU.Units);
    temp.Vel_Select  = Select_Element(temp.Vel_Units);
    temp.Disp_Select = Select_Element(temp.Disp_Units);

    temp.Peak        = Stats.Data[0];
    temp.Mean        = Stats.Data[1];
    temp.RMS         = Stats.Data[2];
    

    return temp;
    
    // Creates new Select-Element for Plotly Info-table (Graph-Unit List)
    function Select_Element(Unit_List) {
        
        // Decleration of variables
        let j, opt, select, ID;

        ID = 'Unit_Plot_ID_' + ChannelList[ChNum].Unique_ID;

        // Create select element and populate it 
        select = document.createElement('select');
        select.setAttribute('id', ID);
        select.setAttribute('class', 'form-select form-select-sm');
        select.setAttribute('onchange', 'Plotly_Graph_Update(' + ChannelList_UniqueID(ChannelList[ChNum].Unique_ID) + ')');
        
        // Options for the select element 
        for (j = 0; j < Unit_List.length; j++) {
            opt = document.createElement("option");
            opt.value = Unit_List[j];
            opt.text = Unit_List[j];
            select.add(opt, null);
        }
        select.selectedIndex = document.getElementById(ID).selectedIndex;
        return select;
    } 

}
//-------------------------------------------------------------------------------------------------------------
function Convert_Data_To_Graph_Unit_SDOF(Data, ChNum) {
    
    let Or_Data, LU, Ind, Indx, temp, ss;

    // SDOF-data using Measurment-index
    Or_Data = TypeAndUnit(ChannelList[ChNum].Results.SDOF.TypeAndUnits); 

    // List of Units of the SDOF-data using Type-Number
    LU = List_Units(Or_Data.Type, false);

    // Indexnumber of User-sepecified Unit on the Plotly Graph
    Ind = document.getElementById("Unit_Plot_ID_" + ChannelList[ChNum].Unique_ID).selectedIndex;

    // Convert Data to user-specified unit on Plotly Graph  ---  Using Unit-Number
    temp = Convert_Units_Data(Data,   Or_Data.Unit,   LU.UnitNum[Ind],   false);
    
    // Double convert the data if applicable 
    if (ChannelList[ChNum].Results.SDOF.DisplayData == "E") {
        temp = Convert_Units_Data(temp.Data,   Or_Data.Unit,   LU.UnitNum[Ind],   false);

        Indx      = document.getElementById('Unit_Plot_ID_'+ChannelList[ChNum].Unique_ID).selectedIndex;
        temp.Unit = document.getElementById('Unit_Plot_ID_'+ChannelList[ChNum].Unique_ID)[Indx].innerHTML;

        Or_Data.Type_String = 'Energy';  
        temp.Unit           = 'Mass • ' + temp.Unit;

    }

    if      (ChannelList[ChNum].Results.SDOF.DisplayData == 'Fs'  )   { Or_Data.Type_String = 'Spring Force';          temp.Unit = 'Mass • ' + temp.Unit; }
    else if (ChannelList[ChNum].Results.SDOF.DisplayData == 'Fc'  )   { Or_Data.Type_String = 'Damping Force';         temp.Unit = 'Mass • ' + temp.Unit; }
    else if (ChannelList[ChNum].Results.SDOF.DisplayData == 'Fi'  )   { Or_Data.Type_String = 'Inertia Force';         temp.Unit = 'Mass • ' + temp.Unit; }
    else if (ChannelList[ChNum].Results.SDOF.DisplayData == 'acc' )   { Or_Data.Type_String = 'Relative Acceleration';                                    }
    else if (ChannelList[ChNum].Results.SDOF.DisplayData == 'Acc' )   { Or_Data.Type_String = 'Total Acceleration';                                       }
    else if (ChannelList[ChNum].Results.SDOF.DisplayData == 'Hyst')   { Or_Data.Type_String = 'Spring Force';          temp.Unit = 'Mass • ' + temp.Unit; }
    else if (ChannelList[ChNum].Results.SDOF.DisplayData == 'HarFor') { Or_Data.Type_String = 'Harmonic Force';        temp.Unit = 'Mass • ' + temp.Unit; }

    // Original statistical values are already scaled by Scale Factor of the channel 
    // Therefore, we just need to conver the units.
    if (ChannelList[ChNum].Results.SDOF.IsAnalysisCompleted) {
        ss    = Statistics(temp.Data);
    } else {
        ss = { Peak:'', Mean:'', RMS:'' };
    }

    //  temp
    temp.yTitle      = '<b>' + Or_Data.Type_String + '  [' + temp.Unit  + ']<b>';
    temp.yTitle_FFT  = '<b>Magnitude<b>';
    temp.y2Title     = '<b>Phase<b>';

    temp.Peak        = ss.Peak;
    temp.Mean        = ss.Mean;
    temp.RMS         = ss.RMS;

    return temp;

}
//-------------------------------------------------------------------------------------------------------------
function Convert_Data_To_Graph_Unit_ResSpec(Data, ChNum) {
    
    let Or_Data, LU, Ind, Indx, temp, ss, arr;

    // SDOF-data using Measurment-index
    Or_Data = TypeAndUnit(ChannelList[ChNum].Results.ResSpec.TypeAndUnits); 

    // List of Units of the SDOF-data using Type-Number
    LU = List_Units(Or_Data.Type, false);

    // Indexnumber of User-sepecified Unit on the Plotly Graph
    Ind = document.getElementById("Unit_Plot_ID_" + ChannelList[ChNum].Unique_ID).selectedIndex;

    // Convert Data to user-specified unit on Plotly Graph  ---  Using Unit-Number
    temp = Convert_Units_Data(Data,   Or_Data.Unit,   LU.UnitNum[Ind],   false);

    if      (ChannelList[ChNum].Results.ResSpec.DisplayData == 'acc' )   { Or_Data.Type_String = 'Relative Acceleration';    }
    else if (ChannelList[ChNum].Results.ResSpec.DisplayData == 'Acc' )   { Or_Data.Type_String = 'Total Acceleration';       }
    else if (ChannelList[ChNum].Results.ResSpec.DisplayData == 'Disp')   { Or_Data.Type_String = 'Displacement Spectrum';    }
    else if (ChannelList[ChNum].Results.ResSpec.DisplayData == 'Vel' )   { Or_Data.Type_String = 'Velocity Spectrum';        }
    else if (ChannelList[ChNum].Results.ResSpec.DisplayData == 'PAcc')   { Or_Data.Type_String = 'Pseudo-Acceleration';      }
    else if (ChannelList[ChNum].Results.ResSpec.DisplayData == 'PVel')   { Or_Data.Type_String = 'Pseudo-Velocity';          }

    // Original statistical values are already scaled by Scale Factor of the channel 
    // Therefore, we just need to conver the units.
    if (ChannelList[ChNum].Results.ResSpec.IsAnalysisCompleted) {
        ss    = Statistics(temp.Data);
    } else {
        ss = { Peak:'', Mean:'', RMS:'' };
    }

    //  temp
    temp.yTitle      = '<b>' + Or_Data.Type_String + '  [' + temp.Unit  + ']<b>';
    temp.yTitle_FFT  = '<b>Magnitude<b>';
    temp.y2Title     = '<b>Phase<b>';

    temp.Peak        = ss.Peak;
    temp.Mean        = ss.Mean;
    temp.RMS         = ss.RMS;
    
    return temp;

}
//-------------------------------------------------------------------------------------------------------------
function Convert_Data_To_Graph_Unit_SM_Par(Data, ChNum) {
    
    let Or_Data, LU, Ind, Indx, temp, ss, arr, leg, T_duration, Type_String;

    // Data using Measurment-index
    Or_Data = TypeAndUnit(ChannelList[ChNum].Results.SM_Parameters.TypeAndUnits);

    // List of Units of the SDOF-data using Type-Number
    LU = List_Units(Or_Data.Type, false);

    // Indexnumber of User-sepecified Unit on the Plotly Graph
    Ind = document.getElementById("Unit_Plot_ID_" + ChannelList[ChNum].Unique_ID).selectedIndex;

    // Convert Data to user-specified unit on Plotly Graph  ---  Using Unit-Number
    temp      = Convert_Units_Data([1],   Or_Data.Unit,   LU.UnitNum[Ind],   false);
    temp.Data = Data;

    if      (ChannelList[ChNum].Results.SM_Parameters.DisplayData == 'AI'  )   { Type_String = 'Arias Intensity';  leg = 'Significant Duration, s';  T_duration = ChannelList[ChNum].Results.SM_Parameters.Ts.toFixed(2); }
    else if (ChannelList[ChNum].Results.SM_Parameters.DisplayData == 'CAV' )   { Type_String = 'CAV';              leg = 'Bracketed Duration, s';    T_duration = ChannelList[ChNum].Results.SM_Parameters.Td.toFixed(2); }

    Or_Data = TypeAndUnit(ChannelList[ChNum].TypeAndUnits);

    //  temp
    temp.yTitle      = '<b>' + Type_String + '  [' + Or_Data.IntegrationUnits.FirstIntegral.Unit_String  + ']<b>';
    temp.leg         = leg;
    temp.T_duration  = T_duration;
    temp.y2Title     = '<b>' + ChannelList[ChNum].TypeString + '  [' + ChannelList[ChNum].UnitString  + ']<b>';

    return temp;

}
//-------------------------------------------------------------------------------------------------------------
function Convert_Data_To_Graph_Unit_Drift(Data, ChNum) {
    
    let Or_Data, LU, Ind, Indx, temp, ss, arr, U;

    Or_Data = TypeAndUnit(ChannelList[ChNum].Results.Drift.TypeAndUnits); 

    // List of Units of the SDOF-data using Type-Number
    LU = List_Units(Or_Data.Type, false);

    // Indexnumber of User-sepecified Unit on the Plotly Graph
    Ind = document.getElementById("Unit_Plot_ID_" + ChannelList[ChNum].Unique_ID).selectedIndex;

    // Convert Data to user-specified unit on Plotly Graph  ---  Using Unit-Number
    temp = Convert_Units_Data(Data,   Or_Data.Unit,   LU.UnitNum[Ind],   false);
    
    // Original statistical values are already scaled by Scale Factor of the channel 
    // Therefore, we just need to conver the units.
    if (ChannelList[ChNum].Results.Drift.IsAnalysisCompleted) {
        ss    = Statistics(temp.Data);
    } else {
        ss = { Peak:'', Mean:'', RMS:'', Residual:'' };
    }

    //  temp
    temp.yTitle      = '<b> Displacement [' + temp.Unit + ']<b>';
    temp.yTitle_FFT  = '<b>Magnitude<b>';
    temp.y2Title     = '<b>Phase<b>';

    temp.Peak        = ss.Peak;
    temp.Mean        = ss.Mean;
    temp.RMS         = ss.RMS;
    temp.Residual    = temp.Data.at(-1);
    
    return temp;

}
//-------------------------------------------------------------------------------------------------------------
function Convert_Data_To_Graph_Unit_Spectrum(Data, ChNum) {
    
    let Or_Data, LU, Ind, Indx, temp, ss, arr, leg, T_duration, aaa;

    // SDOF-data using Measurment-index
    Or_Data = TypeAndUnit(ChannelList[ChNum].Results.SM_Parameters.TypeAndUnits); 

    // List of Units of the SDOF-data using Type-Number
    LU = List_Units(Or_Data.Type, false);

    // Indexnumber of User-sepecified Unit on the Plotly Graph
    Ind = document.getElementById("Unit_Plot_ID_" + ChannelList[ChNum].Unique_ID).selectedIndex;

    // Convert Data to user-specified unit on Plotly Graph  ---  Using Unit-Number
    temp = Convert_Units_Data(Data,   Or_Data.Unit,   LU.UnitNum[Ind],   false);

    if      (ChannelList[ChNum].Results.Spectrum.DisplayData == 'FFT' )   { Or_Data.Type_String = 'Fourier Amplitude'; aaa = '';        }
    else if (ChannelList[ChNum].Results.Spectrum.DisplayData == 'POW' )   { Or_Data.Type_String = 'Power Spectrum';    aaa = ' ²'       }
    else if (ChannelList[ChNum].Results.Spectrum.DisplayData == 'PSD' )   { Or_Data.Type_String = 'Spectral Density';  aaa = ' ² / Hz'  }
    else if (ChannelList[ChNum].Results.Spectrum.DisplayData == 'HET' )   { Or_Data.Type_String = 'Spectrogram';       aaa = ' ²'       }

    // Original statistical values are already scaled by Scale Factor of the channel 
    // Therefore, we just need to conver the units.
    if (ChannelList[ChNum].Results.Spectrum.IsAnalysisCompleted) {
        ss    = Statistics(temp.Data);
    } else {
        ss = { Peak:'', Mean:'', RMS:'' };
    }

    //  temp
    temp.yTitle      = '<b>' + Or_Data.Type_String + '  [ ' + temp.Unit  + ' ]' +aaa+'<b>';

    if (ChannelList[ChNum].Results.Spectrum.DisplayData == 'HET' )   { temp.yTitle      = '<b>Frequency (Hz)<b>'; }

    temp.Peak        = ss.Peak;
    temp.Mean        = ss.Mean;
    temp.RMS         = ss.RMS;
    
    return temp;

}
//-------------------------------------------------------------------------------------------------------------
function Update_Units_infoTable(i) {

    // retrun if no graph to plot
    if (!ChannelList[i].PlotGraph) { return; }

    // Declaration of varibalers 
    let AM, II, Units_SelectElement, Ind, Type, DisplayData;
    let SDOF_Plot_ID, Unit_Cell_ID;

    // Defaults
    Units_SelectElement = Select_Element(List_Units(8).Units);                    
    Type                = 2;  
    DisplayData         = 'Disp';

    // SDOF_ID
    SDOF_Plot_ID = "SDOF_Plot_ID_"  + ChannelList[i].Unique_ID;
    Unit_Cell_ID = "Unit_Cell_ID_"  + ChannelList[i].Unique_ID;

    AM   = document.getElementById('SDOF_Analysis').selectedIndex;   // Analysis Method (Newmark, piesewise, etc) from SDOF Parameetrs Window
    II   = document.getElementById(SDOF_Plot_ID).selectedIndex;      // Index of the SDOF_ResultsDisplay in InfoTable
    
    if (II == -1) { return; }

    if      (AM == 0) {  
        // Free Vibration
        if (II == 0) {  Units_SelectElement = Select_Element(List_Units(8).Units);                    Type=2;  DisplayData='Disp';   } // Displacement
        if (II == 1) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='Vel';    } // Velocity
        if (II == 2) {  Units_SelectElement = Select_Element(['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']);  Type=1;  DisplayData='E';      } // E  ['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']
    } 
    else if (AM == 1) {  
        // Forced Vibration
        if (II == 0) {  Units_SelectElement = Select_Element(List_Units(8).Units);                    Type=2;  DisplayData='Disp';   } // Displacement
        if (II == 1) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='Vel';    } // Velocity
        if (II == 2) {  Units_SelectElement = Select_Element(['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']);  Type=1;  DisplayData='E';      } // E  ['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']
        if (II == 3) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='HarFor'; } // Harmonic Force
    } 
    else if (AM == 2) {  
        // Piece-Wise Exact
        if (II == 0) {  Units_SelectElement = Select_Element(List_Units(8).Units);                    Type=2;  DisplayData='Disp';   } // Displacement
        if (II == 1) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='Vel';    } // Velocity
        if (II == 2) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='acc';    } // Relative acceleration
        if (II == 3) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Acc';    } // Total acceleration
        if (II == 4) {  Units_SelectElement = Select_Element(['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']);  Type=1;  DisplayData='E';      } // E  ['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']
    }
    else if (AM == 3) {  
        // Central Difference Method 
        if (II == 0) {  Units_SelectElement = Select_Element(List_Units(8).Units);                    Type=2;  DisplayData='Disp';   } // Displacement
        if (II == 1) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='Vel';    } // Velocity
        if (II == 2) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='acc';    } // Relative acceleration
        if (II == 3) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Acc';    } // Total acceleration
        if (II == 4) {  Units_SelectElement = Select_Element(['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']);  Type=1;  DisplayData='E';      } // E  ['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']
    } 
    else if (AM == 4) {  
        // Newmark Linear 
        if (II == 0) {  Units_SelectElement = Select_Element(List_Units(8).Units);                    Type=2;  DisplayData='Disp';   } // Displacement
        if (II == 1) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='Vel';    } // Velocity
        if (II == 2) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='acc';    } // Relative acceleration
        if (II == 3) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Acc';    } // Total acceleration
        if (II == 4) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Fs';     } // Spring (restoring) force per unit mass
        if (II == 5) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Fc';     } // Damping force per unit mass
        if (II == 6) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Fi';     } // Inertia force per unit mass
        if (II == 7) {  Units_SelectElement = Select_Element(['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']);  Type=1;  DisplayData='E';      } // E  ['g²s²', 'm²/s²', 'cm²/s²', 'mm²/s²']
    } 
    else if (AM == 5) {  
        // Newmark non-Linear
        if (II == 0) {  Units_SelectElement = Select_Element(List_Units(8).Units);                    Type=2;  DisplayData='Disp';   } // Displacement
        if (II == 1) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='Vel';    } // Velocity
        if (II == 2) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='acc';    } // Relative acceleration
        if (II == 3) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Acc';    } // Total acceleration
        if (II == 4) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Fs';     } // Spring (restoring) force per unit mass
        if (II == 5) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Fc';     } // Damping force per unit mass
        if (II == 6) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Fi';     } // Inertia force per unit mass
        if (II == 7) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Hyst';   } // Fs-Disp  Inertia force per unit mass versus Displacement
    } 

    // Assign select-element to cell-element in table
    document.getElementById(Unit_Cell_ID).innerHTML = "";
    document.getElementById(Unit_Cell_ID).appendChild(Units_SelectElement);

    // Index number of the Unit on the Data-Page (PageNo=0)
    Ind = document.getElementById("Unit_ID_" + ChannelList[i].Unique_ID).selectedIndex;

    // Measurment-Index of the user-selected Display-and-Unit on InfoTable
    ChannelList[i].Results.SDOF.TypeAndUnits = TypeAndUnit(List_Units(Type, false).UnitNum[Ind], false).TypeAndUnit;
    ChannelList[i].Results.SDOF.DisplayData  = DisplayData;

    // Update Graph
    Plotly_Graph_Update(i);

    // Creates new Select-Element for Plotly Info-table (Graph-Unit List)
    function Select_Element(Unit_List) {
        
        // Decleration of variables
        let j, opt, select, ID;

        ID = 'Unit_Plot_ID_' + ChannelList[i].Unique_ID;

        // Create select element and populate it 
        select = document.createElement('select');
        select.setAttribute('id', ID);
        select.setAttribute('class', 'form-select form-select-sm');
        select.setAttribute('onchange', 'Plotly_Graph_Update(' + ChannelList_UniqueID(ChannelList[i].Unique_ID) + ')');
        
        // Options for the select element 
        for (j = 0; j < Unit_List.length; j++) {
            opt = document.createElement("option");
            opt.value = Unit_List[j];
            opt.text = Unit_List[j];
            select.add(opt, null);
        }
        select.selectedIndex = document.getElementById(ID).selectedIndex;
        return select;
    }

}
//-------------------------------------------------------------------------------------------------------------
function Update_Units_infoTable_ResSpec(i) {

    // retrun if no graph to plot
    if (!ChannelList[i].PlotGraph) { return; }

    // Declaration of varibalers 
    let AM, II, Units_SelectElement, Ind, Type, DisplayData;
    let ResSpec_Plot_ID, Unit_Cell_ID;

    // Defaults
    Units_SelectElement = Select_Element(List_Units(8).Units);                    
    Type                = 2;  
    DisplayData         = 'Disp';

    // ResSpec_ID
    ResSpec_Plot_ID = "ResSpec_Plot_ID_"  + ChannelList[i].Unique_ID;
    Unit_Cell_ID    = "Unit_Cell_ID_"     + ChannelList[i].Unique_ID;

    AM   = document.getElementById('ResSpec_AnalysisMethod').selectedIndex;   // Analysis Method (Elastic Spectra, Inelestic spectra, etc) from ResSpec-Parameetrs Window
    II   = document.getElementById(ResSpec_Plot_ID).selectedIndex;            // Index of the ResSpec_ResultsDisplay in InfoTable
    
    if (II == -1) { return; }

    if      (AM == 0) {  
        // Elastic Spectra
        if (II == 0) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Acc';    } // Total Acceleration
        if (II == 1) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='acc';    } // Relative Acceleration
        if (II == 2) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='Vel';    } // Velocity
        if (II == 3) {  Units_SelectElement = Select_Element(List_Units(8).Units);                    Type=2;  DisplayData='Disp';   } // Displacement
        if (II == 4) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='PAcc';   } // Pseudo-Acceleration
        if (II == 5) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='PVel';   } // Pseudo-Velocity
    } 
    else if (AM == 1) {  
        // Constant Ductility Inelastic Spectra (Bilinear Hysteretic model)
        if (II == 0) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Acc';    } // Total Acceleration
        if (II == 1) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='acc';    } // Relative Acceleration
        if (II == 2) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='Vel';    } // Velocity
        if (II == 3) {  Units_SelectElement = Select_Element(List_Units(8).Units);                    Type=2;  DisplayData='Disp';   } // Displacement
    } 
    else if (AM == 2) {  
        // Constant Ductility Inelastic Spectra (Clough Bilinear model)
        if (II == 0) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='Acc';    } // Total Acceleration
        if (II == 1) {  Units_SelectElement = Select_Element(List_Units(0).Units);                    Type=0;  DisplayData='acc';    } // Relative Acceleration
        if (II == 2) {  Units_SelectElement = Select_Element(List_Units(4).Units);                    Type=1;  DisplayData='Vel';    } // Velocity
        if (II == 3) {  Units_SelectElement = Select_Element(List_Units(8).Units);                    Type=2;  DisplayData='Disp';   } // Displacement
    }

    // Assign select-element to cell-element in table
    document.getElementById(Unit_Cell_ID).innerHTML = "";
    document.getElementById(Unit_Cell_ID).appendChild(Units_SelectElement);

    // Index number of the Unit on the Data-Page (PageNo=0)
    Ind = document.getElementById("Unit_ID_" + ChannelList[i].Unique_ID).selectedIndex;

    // Measurment-Index of the user-selected Display-and-Unit on InfoTable
    ChannelList[i].Results.ResSpec.TypeAndUnits = TypeAndUnit(List_Units(Type, false).UnitNum[Ind], false).TypeAndUnit;
    ChannelList[i].Results.ResSpec.DisplayData  = DisplayData;

    // Update Graph
    Plotly_Graph_Update(i);

    // Creates new Select-Element for Plotly Info-table (Graph-Unit List)
    function Select_Element(Unit_List) {
        
        // Decleration of variables
        let j, opt, select, ID;

        ID = 'Unit_Plot_ID_' + ChannelList[i].Unique_ID;

        // Create select element and populate it 
        select = document.createElement('select');
        select.setAttribute('id', ID);
        select.setAttribute('class', 'form-select form-select-sm');
        select.setAttribute('onchange', 'Plotly_Graph_Update(' + ChannelList_UniqueID(ChannelList[i].Unique_ID) + ')');
        
        // Options for the select element 
        for (j = 0; j < Unit_List.length; j++) {
            opt = document.createElement("option");
            opt.value = Unit_List[j];
            opt.text = Unit_List[j];
            select.add(opt, null);
        }
        select.selectedIndex = document.getElementById(ID).selectedIndex;
        return select;
    }

}
//-------------------------------------------------------------------------------------------------------------
function Update_Units_infoTable_SM_Par(i) {

    // retrun if no graph to plot
    if (!ChannelList[i].PlotGraph) { return; }

    // Declaration of varibalers 
    let II, Units_SelectElement, Ind, Type, DisplayData;
    let SM_Par_Plot_ID, Unit_Cell_ID;

    // Defaults
    Units_SelectElement = Select_Element(List_Units(8).Units);                    
    Type                = 2;  
    DisplayData         = 'Disp';

    // ResSpec_ID
    SM_Par_Plot_ID = "SM_Par_Plot_ID_"  + ChannelList[i].Unique_ID;
    Unit_Cell_ID   = "Unit_Cell_ID_"    + ChannelList[i].Unique_ID;

    II   = document.getElementById(SM_Par_Plot_ID).selectedIndex;            // Index of the SM_Par_ResultsDisplay in InfoTable
    
    if (II == -1) { return; }

    if      (II == 0) { Units_SelectElement = Select_Element(List_Units(4).Units);   Type=1;    DisplayData='AI';   } // Arias Intensity ( velocity unit)
    else if (II == 1) { Units_SelectElement = Select_Element(List_Units(4).Units);   Type=1;    DisplayData='CAV';  } // CAV  ( velocity unit)
    

    // Assign select-element to cell-element in table
    document.getElementById(Unit_Cell_ID).innerHTML = "";
    document.getElementById(Unit_Cell_ID).appendChild(Units_SelectElement);

    // Index number of the Unit on the Data-Page (PageNo=0)
    Ind = document.getElementById("Unit_ID_" + ChannelList[i].Unique_ID).selectedIndex;

    // Measurment-Index of the user-selected Display-and-Unit on InfoTable
    ChannelList[i].Results.SM_Parameters.TypeAndUnits = TypeAndUnit(List_Units(Type, false).UnitNum[Ind], false).TypeAndUnit;
    ChannelList[i].Results.SM_Parameters.DisplayData  = DisplayData;

    // Update Graph
    Plotly_Graph_Update(i);

    // Creates new Select-Element for Plotly Info-table (Graph-Unit List)
    function Select_Element(Unit_List) {
        
        // Decleration of variables
        let j, opt, select, ID;

        ID = 'Unit_Plot_ID_' + ChannelList[i].Unique_ID;

        // Create select element and populate it 
        select = document.createElement('select');
        select.setAttribute('id', ID);
        select.setAttribute('class', 'form-select form-select-sm');
        select.setAttribute('onchange', 'Plotly_Graph_Update(' + ChannelList_UniqueID(ChannelList[i].Unique_ID) + ')');
        
        // Options for the select element 
        for (j = 0; j < Unit_List.length; j++) {
            opt = document.createElement("option");
            opt.value = Unit_List[j];
            opt.text = Unit_List[j];
            select.add(opt, null);
        }
        select.selectedIndex = document.getElementById(ID).selectedIndex;
        return select;
    }

}
//-------------------------------------------------------------------------------------------------------------
function Update_Units_infoTable_Drift(i) {

    // retrun if no graph to plot
    if (!ChannelList[i].PlotGraph) { return; }

    // Declaration of varibalers 
    let AM, II, Units_SelectElement, Ind, Type, DisplayData;
    let ResSpec_Plot_ID, Unit_Cell_ID;

    // Defaults
    Units_SelectElement = Select_Element(List_Units(8).Units);                    
    Type                = 2;  

    // ResSpec_ID
    ResSpec_Plot_ID = "ResSpec_Plot_ID_"  + ChannelList[i].Unique_ID;
    Unit_Cell_ID    = "Unit_Cell_ID_"     + ChannelList[i].Unique_ID;

    // Assign select-element to cell-element in table
    document.getElementById(Unit_Cell_ID).innerHTML = "";
    document.getElementById(Unit_Cell_ID).appendChild(Units_SelectElement);

    // Index number of the Unit on the Data-Page (PageNo=0)
    Ind = document.getElementById("Unit_ID_" + ChannelList[i].Unique_ID).selectedIndex;

    // Measurment-Index of the user-selected Display-and-Unit on InfoTable
    ChannelList[i].Results.Drift.TypeAndUnits = TypeAndUnit(List_Units(Type, false).UnitNum[Ind], false).TypeAndUnit;

    // Update Graph
    Plotly_Graph_Update(i);

    // Creates new Select-Element for Plotly Info-table (Graph-Unit List)
    function Select_Element(Unit_List) {
        
        // Decleration of variables
        let j, opt, select, ID;

        ID = 'Unit_Plot_ID_' + ChannelList[i].Unique_ID;

        // Create select element and populate it 
        select = document.createElement('select');
        select.setAttribute('id', ID);
        select.setAttribute('class', 'form-select form-select-sm');
        select.setAttribute('onchange', 'Plotly_Graph_Update(' + ChannelList_UniqueID(ChannelList[i].Unique_ID) + ')');
        
        // Options for the select element 
        for (j = 0; j < Unit_List.length; j++) {
            opt = document.createElement("option");
            opt.value = Unit_List[j];
            opt.text = Unit_List[j];
            select.add(opt, null);
        }
        select.selectedIndex = document.getElementById(ID).selectedIndex;
        return select;
    }

}
//-------------------------------------------------------------------------------------------------------------
function Update_Units_infoTable_Spectrum(i) {

    // retrun if no graph to plot
    if (!ChannelList[i].PlotGraph) { return; }

    // Declaration of varibalers 
    let II, Units_SelectElement, Ind, Type, DisplayData;
    let Spectrum_Plot_ID, Unit_Cell_ID;

    // Defaults
    Units_SelectElement = Select_Element(List_Units(8).Units);                    
    Type                = 2;  
    DisplayData         = 'Disp';

    // ResSpec_ID
    Spectrum_Plot_ID = "Spectrum_Plot_ID_"  + ChannelList[i].Unique_ID;
    Unit_Cell_ID      = "Unit_Cell_ID_"     + ChannelList[i].Unique_ID;

    II   = document.getElementById(Spectrum_Plot_ID).selectedIndex;            // Index of the Spectrum_ResultsDisplay in InfoTable
    
    if (II == -1) { return; }

    if      (II == 0) { Units_SelectElement = Select_Element(List_Units(ChannelList[i].Unit).Units);   Type = ChannelList[i].Type;    DisplayData='FFT';   } // Fourier Amplitude Spectrum
    else if (II == 1) { Units_SelectElement = Select_Element(List_Units(ChannelList[i].Unit).Units);   Type = ChannelList[i].Type;    DisplayData='POW';   } // Power Spectrum
    else if (II == 2) { Units_SelectElement = Select_Element(List_Units(ChannelList[i].Unit).Units);   Type = ChannelList[i].Type;    DisplayData='PSD';   } // Spectral Density Spectrum
    else if (II == 3) { Units_SelectElement = Select_Element(List_Units(ChannelList[i].Unit).Units);   Type = ChannelList[i].Type;    DisplayData='HET';   } // Spectrogram
    

    // Assign select-element to cell-element in table
    document.getElementById(Unit_Cell_ID).innerHTML = "";
    document.getElementById(Unit_Cell_ID).appendChild(Units_SelectElement);

    // Index number of the Unit on the Data-Page (PageNo=0)
    Ind = document.getElementById("Unit_ID_" + ChannelList[i].Unique_ID).selectedIndex;

    // Measurment-Index of the user-selected Display-and-Unit on InfoTable
    ChannelList[i].Results.Spectrum.TypeAndUnits = TypeAndUnit(List_Units(Type, false).UnitNum[Ind], false).TypeAndUnit;
    ChannelList[i].Results.Spectrum.DisplayData  = DisplayData;

    // Update Graph
    Plotly_Graph_Update(i);

    // Creates new Select-Element for Plotly Info-table (Graph-Unit List)
    function Select_Element(Unit_List) {
        
        // Decleration of variables
        let j, opt, select, ID;

        ID = 'Unit_Plot_ID_' + ChannelList[i].Unique_ID;

        // Create select element and populate it 
        select = document.createElement('select');
        select.setAttribute('id', ID);
        select.setAttribute('class', 'form-select form-select-sm');
        select.setAttribute('onchange', 'Plotly_Graph_Update(' + ChannelList_UniqueID(ChannelList[i].Unique_ID) + ')');
        
        // Options for the select element 
        for (j = 0; j < Unit_List.length; j++) {
            opt = document.createElement("option");
            opt.value = Unit_List[j];
            opt.text = Unit_List[j];
            select.add(opt, null);
        }
        select.selectedIndex = document.getElementById(ID).selectedIndex;
        return select;
    }

}
//-------------------------------------------------------------------------------------------------------------


