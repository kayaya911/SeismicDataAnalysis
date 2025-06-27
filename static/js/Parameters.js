// stricter parsing and error handling
"use strict";

// This function hides/shows parameters DIV
function Parameters_Button() {
    var x = document.getElementById("Parameter_Container");
    if (x.style.display === "none") {
      x.style.display = "flex";
    } else {
      x.style.display = "none";
    }
}

// 