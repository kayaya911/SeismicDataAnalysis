// stricter parsing and error handling
"use strict";


function sleep(t) {
    return new Promise((resolve) => setTimeout(resolve, t));
}


function Generate_Unique_ID() {
  let id_1, temp;
  do {
    temp = crypto.randomUUID();
    id_1 = `${temp}`;
  } while (document.getElementById(id_1));
  
  // replace all "-" with "_" before returning 
  return id_1.replace(/-/gi,"_");
}
