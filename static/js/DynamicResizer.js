// stricter parsing and error handling
"use strict";

//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
// Vertical resizer
const verticalResizer  = document.getElementById('Vertical_Divider');
const treeMenu         = document.getElementById('File_Tree_Containe');
let isResizingVertical = false;

verticalResizer.addEventListener('mousedown', (e) => {
    isResizingVertical = true;
    document.addEventListener('mousemove', resizeVertical);
    document.addEventListener('mouseup', stopResizeVertical);
});

function resizeVertical(e) {
    if (!isResizingVertical) return;
    // The clientX property of the MouseEvent object provides the horizontal coordinate of the mouse cursor (in pixels) 
    // relative to the viewport (the visible area of the browser window), not the entire page.
    const newWidth = e.clientX;
    if (newWidth > 100 && newWidth < window.innerWidth - 100) {
        treeMenu.style.flex = `0 0 ${newWidth}px`;
    }
}

function stopResizeVertical() {
    isResizingVertical = false;
    document.removeEventListener('mousemove', resizeVertical);
    document.removeEventListener('mouseup', stopResizeVertical);
}


//-------------------------------------------------------------------------------------------------------------
// Horizontal resizer
const horizontalResizer   = document.getElementById('Horizontal_Divider');
const resultsParameters   = document.getElementById('Parameter_Container');   
let isResizingHorizontal = false;

horizontalResizer.addEventListener('mousedown', (e) => {
    isResizingHorizontal = true;
    document.addEventListener('mousemove', resizeHorizontal);
    document.addEventListener('mouseup', stopResizeHorizontal);
});

function resizeHorizontal(e) {
    if (!isResizingHorizontal) return;
    const resultsMenu = document.querySelector('.Results_Container');
    const newHeight = e.clientY - resultsMenu.offsetTop;
    if (newHeight > 100 && newHeight < resultsMenu.offsetHeight - 100) {
        resultsParameters.style.flex = `0 0 ${newHeight}px`;
    }
}

function stopResizeHorizontal() {
    isResizingHorizontal = false;
    document.removeEventListener('mousemove', resizeHorizontal);
    document.removeEventListener('mouseup', stopResizeHorizontal);
}
