// stricter parsing and error handling
"use strict";

//------------------------------------------------------------------------------------------------------------------------------------------------
// DISCLAIMER MODAL
//
// Shows a one-time acknowledgment modal on first visit. Stores acceptance in
// localStorage so it will not show again unless DISCLAIMER_VERSION is bumped
// (e.g. after editing the disclaimer text) or the user clears site data.
//------------------------------------------------------------------------------------------------------------------------------------------------

const DISCLAIMER_STORAGE_KEY = 'SeismicWebtool_DisclaimerAccepted';
const DISCLAIMER_VERSION     = '1.0'; // bump this string to force re-acceptance from all users

//------------------------------------------------------------------------------------------------------------------------------------------------
function Disclaimer_Init() {

    // Decleration of variables
    let Overlay, Checkbox, AcceptBtn, stored, record;

    Overlay   = document.getElementById('Disclaimer_Overlay');
    Checkbox  = document.getElementById('Disclaimer_Checkbox');
    AcceptBtn = document.getElementById('Disclaimer_Accept_Btn');

    if (!Overlay) { return; }

    // Check if user has already accepted this version
    stored = localStorage.getItem(DISCLAIMER_STORAGE_KEY);
    if (stored) {
        try {
            record = JSON.parse(stored);
            if (record.version === DISCLAIMER_VERSION) {
                return; // already accepted current version, do not show
            }
        } catch (e) {
            // malformed record in storage, fall through and show modal
        }
    }

    // Show the modal
    Overlay.classList.add('Disclaimer_Visible');

    // Wire up checkbox -> enables/disables Accept button
    Checkbox.addEventListener('change', () => {
        AcceptBtn.classList.toggle('Disclaimer_Enabled', Checkbox.checked);
    });

    // Wire up Accept button
    AcceptBtn.addEventListener('click', () => {
        if (!Checkbox.checked) { return; }
        Disclaimer_Accept();
    });
}
//------------------------------------------------------------------------------------------------------------------------------------------------
function Disclaimer_Accept() {

    // Decleration of variables
    let record, Overlay;

    record = {
        version:    DISCLAIMER_VERSION,
        acceptedAt: new Date().toISOString()
    };

    localStorage.setItem(DISCLAIMER_STORAGE_KEY, JSON.stringify(record));

    Overlay = document.getElementById('Disclaimer_Overlay');
    Overlay.classList.remove('Disclaimer_Visible');
}
//------------------------------------------------------------------------------------------------------------------------------------------------
document.addEventListener('DOMContentLoaded', Disclaimer_Init);
//------------------------------------------------------------------------------------------------------------------------------------------------