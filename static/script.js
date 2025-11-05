/* ===================================================================
   C.ORIGAMI VIEWER — script.js
   (Plain JS, no CSS or HTML)
   =================================================================== */
   "use strict";

   /*************************************************************
    * Constants and Global Variables
    *************************************************************/
   const WINDOW_WIDTH = 2097152; // 2 Mb window
   window.chromosomeLengths = {};
   window.initComplete = false;
   
   /*************************************************************
    * Local Storage Helpers
    *************************************************************/
   function storeFormFields() {
     // Save first region and other fields (exclude region_end1 here!)
     const fields = [
       "region_chr1", "region_start1",
       "del_start", "del_width",
       "perturb_width",
       "model_select", "atac_bw_path",
       "ctcf_bw_path", "peaks_file_path"
     ];
     fields.forEach(field => {
       const elem = document.getElementById(field);
       if (elem) localStorage.setItem(field, elem.value);
     });
   
     // Special handling for region_end1: only store it if user actually typed something
     const endField1 = document.getElementById("region_end1");
     if (endField1 && endField1.value.trim()) {
       localStorage.setItem("region_end1", endField1.value);
     } else {
       localStorage.removeItem("region_end1");
     }
   
     // Same idea for region_end2:
     const endField2 = document.getElementById("region_end2");
     if (endField2 && endField2.value.trim()) {
       localStorage.setItem("region_end2", endField2.value);
     } else {
       localStorage.removeItem("region_end2");
     }
   
     // Save second region fields (chr + start):
     const secondFields = ["region_chr2", "region_start2"];
     secondFields.forEach(field => {
       const elem = document.getElementById(field);
       if (elem) localStorage.setItem(field, elem.value);
     });
   
     // Save ds_option
     const dsRadios = document.getElementsByName("ds_option");
     for (let radio of dsRadios) {
       if (radio.checked) {
         localStorage.setItem("ds_option", radio.value);
         break;
       }
     }
   }
   
   function restoreFormFields() {
     const fields = [
       "region_chr1", "region_start1",
       "del_start", "del_width",
       "perturb_width",
       "model_select", "atac_bw_path",
       "ctcf_bw_path", "peaks_file_path",
       "genome_select"
     ];
     fields.forEach(field => {
       const storedVal = localStorage.getItem(field);
       const elem = document.getElementById(field);
       if (storedVal && elem) elem.value = storedVal;
     });
   
     // region_end1
     const storedEnd1 = localStorage.getItem("region_end1");
     const endElem1 = document.getElementById("region_end1");
     if (endElem1) {
       if (storedEnd1) endElem1.value = storedEnd1;
       else endElem1.value = "";
     }
   
     // region_end2
     const storedEnd2 = localStorage.getItem("region_end2");
     const endElem2 = document.getElementById("region_end2");
     if (endElem2) {
       if (storedEnd2) endElem2.value = storedEnd2;
       else endElem2.value = "";
     }
   
     // region 2 basics
     ["region_chr2", "region_start2"].forEach(field => {
       const storedVal = localStorage.getItem(field);
       const elem = document.getElementById(field);
       if (storedVal && elem) elem.value = storedVal;
     });
   
     const dsOption = localStorage.getItem("ds_option");
     if (dsOption) {
       const radioElem = document.getElementById("ds_" + dsOption);
       if (radioElem) radioElem.checked = true;
     }
   }
   
   /*************************************************************
    * Form Behavior Helpers
    *************************************************************/
   function populateChromosomeDropdown() {
     const genomeSelect = document.getElementById("genome_select");
     let genome = genomeSelect ? genomeSelect.value : "hg38";
     let url = (genome === "mm10") ? "static/mm10_chr_lengths.json" : "static/hg38_chr_lengths.json";
   
     return fetch(url)
       .then(response => response.json())
       .then(data => {
         window.chromosomeLengths = data;
   
         const chrSelect1 = document.getElementById("region_chr1");
         const chrSelect2 = document.getElementById("region_chr2");
         [chrSelect1, chrSelect2].forEach(select => {
           if (!select) return;
           select.innerHTML = "";
           for (const chr in data) {
             const cleanChr = chr.replace(/^chr/i, "");
             const option = document.createElement("option");
             option.value = chr;     // e.g. "chr1"
             option.text = cleanChr; // e.g. "1", "2", "X"
             select.appendChild(option);
           }
         });
         console.log("Chromosome dropdown populated for", genome);
       })
       .catch(error => {
         console.error("Error fetching chromosome lengths:", error);
       });
   }
   
   function updateEndPosition() {
     // Region 1
     const startField1 = document.getElementById("region_start1");
     const endField1   = document.getElementById("region_end1");
     const startVal1   = parseInt(startField1.value, 10);
     if (!isNaN(startVal1)) {
       const computedEnd1 = startVal1 + WINDOW_WIDTH;
       if (endField1.disabled) {
         // deletion mode: bake computed end into value
         endField1.value = computedEnd1;
         endField1.placeholder = "";
       } else if (!endField1.value || isNaN(parseInt(endField1.value, 10))) {
         endField1.placeholder = String(computedEnd1);
       } else {
         endField1.placeholder = "";
       }
     } else {
       endField1.placeholder = "";
     }
   
     // Region 2
     const secondRegionDiv = document.getElementById("second_chr_fields");
     if (secondRegionDiv && secondRegionDiv.style.display !== "none") {
       const startField2 = document.getElementById("region_start2");
       const endField2   = document.getElementById("region_end2");
       const startVal2   = parseInt(startField2.value, 10);
       if (!isNaN(startVal2)) {
         const computedEnd2 = startVal2 + WINDOW_WIDTH;
         if (!endField2.value || isNaN(parseInt(endField2.value, 10))) {
           endField2.placeholder = String(computedEnd2);
         } else {
           endField2.placeholder = "";
         }
       } else {
         endField2.placeholder = "";
       }
     }
   }
   
   /*************************************************************
    * Validation Helpers
    *************************************************************/
   function validateRegion1Bounds() {
     const chrSelect = document.getElementById("region_chr1");
     const startElem = document.getElementById("region_start1");
     const endElem   = document.getElementById("region_end1");
   
     // If something’s not loaded, bail out
     if (!chrSelect || !window.chromosomeLengths) return null;
   
     const selectedChr = chrSelect.value;
     const maxLength   = window.chromosomeLengths[selectedChr];
   
     if (typeof maxLength === "undefined") return null;
   
     const start = parseInt(startElem.value, 10);
     let end     = parseInt(endElem.value, 10);
     if (isNaN(end)) {
       end = start + WINDOW_WIDTH; // default to 2 Mb from start
     }
   
     // Out-of-bounds checks
     if (start < 0 || start > maxLength) {
       return `Region 1: Start position (${start}) is out of bounds for ${selectedChr}, which has length ${maxLength}.`;
     }
     if (end < 0 || end > maxLength) {
       return `Region 1: End position (${end}) is out of bounds for ${selectedChr}, which has length ${maxLength}.`;
     }
   
     // Logical checks
     if (end < start) {
       return "Region 1: Chromosome start position cannot be greater than end";
     }
     if (end - start > WINDOW_WIDTH * 10) {
       return "Region 1: window must be under 20 Mb in length";
     }
   
     return null;
   }
   
   function validateRegion2Bounds() {
     // Only validate if the second region is visible
     const secondChrDiv = document.getElementById("second_chr_fields");
     if (!secondChrDiv || secondChrDiv.style.display === "none") {
       return null;
     }
   
     const chrSelect2 = document.getElementById("region_chr2");
     const startElem2 = document.getElementById("region_start2");
     const endElem2   = document.getElementById("region_end2");
   
     if (!chrSelect2 || !window.chromosomeLengths) return null;
   
     const selectedChr2 = chrSelect2.value;
     const maxLength2   = window.chromosomeLengths[selectedChr2];
     if (typeof maxLength2 === "undefined") return null;
   
     const start2 = parseInt(startElem2.value, 10);
     let end2     = parseInt(endElem2.value, 10);
     if (isNaN(end2)) {
       end2 = start2 + WINDOW_WIDTH;
     }
   
     // Out-of-bounds checks
     if (start2 < 0 || start2 > maxLength2) {
       return `Region 2: Start position (${start2}) is out of bounds for ${selectedChr2}, which has length ${maxLength2}.`;
     }
     if (end2 < 0 || end2 > maxLength2) {
       return `Region 2: End position (${end2}) is out of bounds for ${selectedChr2}, which has length ${maxLength2}.`;
     }
   
     // Logical checks
     if (end2 < start2) {
       return "Region 2: Chromosome start position cannot be greater than end";
     }
     if (end2 - start2 > WINDOW_WIDTH * 10) {
       return "Region 2: window must be under 20 Mb in length";
     }
   
     return null;
   }
   
   function validateDeletionAreaErrors() {
     const dsOption = document.querySelector('input[name="ds_option"]:checked');
     if (!dsOption || dsOption.value !== "deletion") {
       return null;
     }
   
     const regionStartElem = document.getElementById("region_start1");
     const regionEndElem   = document.getElementById("region_end1");
     const delStartElem    = document.getElementById("del_start");
     const delWidthElem    = document.getElementById("del_width");
   
     const regionStart   = parseInt(regionStartElem.value, 10);
     const regionEnd     = regionEndElem.value ? parseInt(regionEndElem.value, 10) : regionStart + WINDOW_WIDTH;
     const deletionStart = parseInt(delStartElem.value, 10);
     const deletionWidth = parseInt(delWidthElem.value, 10);
     const deletionEnd   = deletionStart + deletionWidth;
   
     if (deletionStart < regionStart || deletionEnd > regionEnd) {
       return "Deletion area is out of bounds for the selected region.";
     }
     return null;
   }
   
   /**
    * Additional check for combined Region 1 + Region 2 total <= 20 Mb
    */
   function validateTotalRegionSize() {
     const start1 = parseInt(document.getElementById("region_start1").value, 10);
     let end1     = parseInt(document.getElementById("region_end1").value, 10);
     if (isNaN(end1)) {
       end1 = start1 + WINDOW_WIDTH;
     }
     const size1 = isNaN(start1) ? 0 : (end1 - start1);
   
     let size2 = 0;
     const secondChrDiv = document.getElementById("second_chr_fields");
     if (secondChrDiv && secondChrDiv.style.display !== "none") {
       const start2 = parseInt(document.getElementById("region_start2").value, 10);
       let end2     = parseInt(document.getElementById("region_end2").value, 10);
       if (isNaN(end2)) {
         end2 = start2 + WINDOW_WIDTH;
       }
       size2 = isNaN(start2) ? 0 : (end2 - start2);
     }
   
     if ((size1 + size2) > 20000000) {
       return "Combined region sizes must not exceed 20 Mb in total.";
     }
     return null;
   }
   
   /*************************************************************
    * Aggregator for All Errors
    *************************************************************/
   function collectAndDisplayErrors() {
     if (!window.initComplete) {
       console.log("Skipping validation because init is not complete.");
       return;
     }
     const errorContainer = document.getElementById("global-error-module");
     if (!errorContainer) return;
   
     const errors = [];
     const warnings = [];
   
     // Existing checks
     const eReg1 = validateRegion1Bounds();
     if (eReg1) errors.push(eReg1);
   
     const eReg2 = validateRegion2Bounds();
     if (eReg2) errors.push(eReg2);
   
     const eDel = validateDeletionAreaErrors();
     if (eDel) errors.push(eDel);
   
     const eTotal = validateTotalRegionSize();
     if (eTotal) errors.push(eTotal);
   
     // New: screening edges warning
     const eScreenEdges = checkScreeningRegionEdges();
     if (eScreenEdges) {
       if (eScreenEdges.isWarning) {
         warnings.push(eScreenEdges.message);
       } else if (eScreenEdges.isError) {
         errors.push(eScreenEdges.message);
       }
     }
   
     // New: minimum region size warning
     const eMinSize = checkMinimumRegionSize();
     if (eMinSize && eMinSize.isWarning) {
       warnings.push(eMinSize.message);
     }
   
     // Clear + render
     errorContainer.innerHTML = "";
   
     if (errors.length === 0 && warnings.length === 0) {
       errorContainer.style.display = "none";
     } else {
       errorContainer.style.display = "block";
   
       if (errors.length > 0) {
         const ulErrors = document.createElement("ul");
         ulErrors.classList.add("global-error-list");
         errors.forEach(msg => {
           const li = document.createElement("li");
           li.textContent = msg;
           ulErrors.appendChild(li);
         });
         errorContainer.appendChild(ulErrors);
       }
   
       if (warnings.length > 0) {
         const ulWarn = document.createElement("ul");
         ulWarn.classList.add("global-warning-list");
         warnings.forEach(msg => {
           const li = document.createElement("li");
           li.textContent = msg;
           ulWarn.appendChild(li);
         });
         errorContainer.appendChild(ulWarn);
       }
     }
   
     const submitBtn = document.querySelector('input[type="submit"].submit-button');
     if (submitBtn) {
       if (errors.length > 0 || !checkRequiredFieldsFilled()) {
         submitBtn.disabled = true;
         submitBtn.classList.add("disabled");
       } else {
         submitBtn.disabled = false;
         submitBtn.classList.remove("disabled");
       }
     }
   }
   
   /*************************************************************
    * Toggles & Additional Behaviors
    *************************************************************/
   function showSecondChr() {
     const secondChrDiv  = document.getElementById('second_chr_fields');
     const togglePlus    = document.getElementById('toggle-second-chr');
     const chimericInput = document.getElementById('chimeric_active');
     const peaksInput    = document.getElementById("peaks_file_upload");
     const peaksDropdown = document.getElementById("peaks_file_path");
   
     if (secondChrDiv) secondChrDiv.style.display = "block";
     if (togglePlus) togglePlus.classList.add("active");
     if (chimericInput) chimericInput.value = "true";
     localStorage.setItem("chimeric_active", "true");
   
     if (peaksInput) {
       peaksInput.disabled = true;
     }
   
     if (peaksDropdown) {
       if (!peaksDropdown.dataset.originalOptions) {
         peaksDropdown.dataset.originalOptions = peaksDropdown.innerHTML;
       }
       peaksDropdown.innerHTML = '<option value="auto_generated">Automatically generate</option>';
       peaksDropdown.value = "auto_generated";
       peaksDropdown.disabled = true;
     }
   
     collectAndDisplayErrors();
     updateNormLocks();
   }
   
   function hideSecondChr() {
     const secondChrDiv  = document.getElementById('second_chr_fields');
     const togglePlus    = document.getElementById('toggle-second-chr');
     const chimericInput = document.getElementById('chimeric_active');
     const peaksInput    = document.getElementById("peaks_file_upload");
     const peaksDropdown = document.getElementById("peaks_file_path");
   
     if (secondChrDiv) secondChrDiv.style.display = "none";
     if (togglePlus) togglePlus.classList.remove("active");
     if (chimericInput) chimericInput.value = "false";
     localStorage.setItem("chimeric_active", "false");
   
     if (peaksInput) {
       peaksInput.disabled = false;
     }
   
     if (peaksDropdown && peaksDropdown.dataset.originalOptions) {
       peaksDropdown.innerHTML = peaksDropdown.dataset.originalOptions;
       peaksDropdown.disabled = false;
     }
   
     collectAndDisplayErrors();
     updateNormLocks();
   }
   
   function toggleOptionalFields() {
     const dsOption = document.querySelector('input[name="ds_option"]:checked');
     const delFields = document.getElementById("deletion-fields");
     const scrFields = document.getElementById("screening-fields-2");
     const peaksContainer = document.getElementById("peaks_file_container");
     const togglePlusElemContainer = document.getElementById("toggle-second-chr-container");
     const endField1 = document.getElementById("region_end1");
   
     if (!dsOption) return;
   
     if (dsOption.value === "deletion") {
       if (delFields) delFields.style.display = "block";
       if (scrFields) scrFields.style.display = "none";
       if (peaksContainer) peaksContainer.style.display = "none";
       hideSecondChr(); // always hide second region in deletion mode
       if (togglePlusElemContainer) togglePlusElemContainer.style.display = "none";
   
       if (endField1) {
         const startVal = parseInt(document.getElementById("region_start1").value, 10) || 0;
         const computedEnd = startVal + WINDOW_WIDTH;
         endField1.value = computedEnd;
         endField1.disabled = true;
         endField1.classList.add("disabled-field");
       }
     } else if (dsOption.value === "screening") {
       if (scrFields) scrFields.style.display = "block";
       if (delFields) delFields.style.display = "none";
       if (peaksContainer) peaksContainer.style.display = "block";
       if (togglePlusElemContainer) togglePlusElemContainer.style.display = "block";
   
       if (endField1) {
         endField1.disabled = false;
         endField1.classList.remove("disabled-field");
       }
     } else {
       if (delFields) delFields.style.display = "none";
       if (scrFields) scrFields.style.display = "none";
       if (peaksContainer) peaksContainer.style.display = "none";
       if (togglePlusElemContainer) togglePlusElemContainer.style.display = "block";
   
       if (endField1) {
         endField1.disabled = false;
         endField1.classList.remove("disabled-field");
       }
     }
   }
   
   function toggleSecondChr() {
     const dsOption = document.querySelector('input[name="ds_option"]:checked');
     if (dsOption && dsOption.value === "deletion") {
       hideSecondChr();
     } else {
       const secondChrDiv = document.getElementById('second_chr_fields');
       if (secondChrDiv && (secondChrDiv.style.display === "none" || secondChrDiv.style.display === "")) {
         showSecondChr();
         updateEndPosition();
       } else {
         hideSecondChr();
       }
     }
   }
   
   function checkFormRequirements() {
     const dsOption = document.querySelector('input[name="ds_option"]:checked');
     if (!dsOption) return true;
   
     if (dsOption.value === "deletion") {
       if (!document.getElementById("del_start").value.trim() ||
           !document.getElementById("del_width").value.trim()) {
         alert("Please provide both Deletion Start and Deletion Width for Deletion mode.");
         return false;
       }
     } else if (dsOption.value === "screening") {
       if (!document.getElementById("perturb_width").value.trim()) {
         alert("Please provide both Perturb Width and Step Size for Screening mode.");
         return false;
       }
     }
     return true;
   }
   
   /*************************************************************
    * File Dropdown Handling
    *************************************************************/
   function populateDropdownFromServer(selectId, fileType) {
     const selectElem = document.getElementById(selectId);
     return fetch("/list_uploads?file_type=" + encodeURIComponent(fileType))
       .then(response => response.json())
       .then(serverFiles => {
         serverFiles.forEach(file => {
           const alreadyExists = Array.from(selectElem.options)
             .some(opt => opt.value === file.value);
           if (!alreadyExists) {
             const newOption = document.createElement("option");
             newOption.value = file.value;
             newOption.text = file.name;
             selectElem.appendChild(newOption);
           }
         });
         console.log("Dropdown", selectId, "populated with:", serverFiles);
       })
       .catch(err => console.error(`Error populating ${fileType} dropdown:`, err));
   }
   
   /*************************************************************
    * AJAX File Upload Helper
    *************************************************************/
   function ajaxUploadFile(fileInputId, fileType, dropdownId) {
     const fileInput = document.getElementById(fileInputId);
     if (!fileInput.files || fileInput.files.length === 0) {
       console.log("No file selected for", fileType);
       return;
     }
   
     const file = fileInput.files[0];
   
     // 1) Enforce maximum file size of 5 GB
     const maxSizeBytes = 5 * 1024 * 1024 * 1024; // 5 GB
     if (file.size > maxSizeBytes) {
       alert("File is too large! Must be under 5 GB.");
       fileInput.value = "";
       return;
     }
   
     // 2) Enforce correct file extension based on fileType
     const ext = file.name.toLowerCase().split('.').pop();
     if (fileType === "atac" || fileType === "ctcf") {
       if (ext !== "bw" && ext !== "bigwig") {
         alert("File extension not recognized as .bw or .bigWig. Please upload a bigWig file.");
         fileInput.value = "";
         return;
       }
     } else if (fileType === "peaks") {
       const validPeaksExt = ["bed", "narrowpeak", "broadpeak"];
       if (!validPeaksExt.includes(ext)) {
         alert("Please upload a peaks file (.bed, .narrowPeak, or .broadPeak).");
         fileInput.value = "";
         return;
       }
     }
   
     console.log("Selected file for " + fileType + ":", file.name);
   
     const uploadButton = document.querySelector("label[for='" + fileInputId + "']");
     if (uploadButton) {
       uploadButton.textContent = 'Uploading...';
       uploadButton.classList.add('loading');
     }
   
     // Construct FormData for POST
     const formData = new FormData();
     if (fileType === "peaks") {
       formData.append("peaks_file", file);
     } else {
       formData.append(fileType + "_bw_file", file);
     }
   
     fetch("/upload_file?file_type=" + encodeURIComponent(fileType), {
       method: "POST",
       body: formData
     })
       .then(response => response.json())
       .then(data => {
         if (data.error) {
           console.error("Upload error for " + fileType + ":", data.error);
           alert("Error uploading file: " + data.error);
         } else {
           console.log("Upload successful for " + fileType + ":", data);
           const dropdown = document.getElementById(dropdownId);
           const newOption = document.createElement("option");
           newOption.value = data.saved_path;
           newOption.text  = data.display_name;
           dropdown.appendChild(newOption);
           dropdown.value = data.saved_path;
           localStorage.setItem(dropdownId, data.saved_path);
         }
       })
       .catch(error => {
         console.error("Error uploading file for " + fileType + ":", error);
         alert("Error uploading file. Check console for details.");
       })
       .finally(() => {
         if (uploadButton) {
           uploadButton.textContent = 'Upload';
           uploadButton.classList.remove('loading');
         }
         collectAndDisplayErrors(); // re-validate after upload
         updateNormLocks();
       });
   }
   
   /*************************************************************
    * Helper to Execute Inline Scripts After AJAX Update
    *************************************************************/
   function executeScripts(container) {
     console.log("Executing scripts from container:", container);
     const scripts = container.querySelectorAll("script");
     scripts.forEach(script => {
       const newScript = document.createElement("script");
       if (script.src) {
         newScript.src = script.src;
         newScript.onload = () => console.log("Loaded external script:", script.src);
       } else {
         newScript.text = script.textContent;
         console.log("Executing inline script:", script.textContent.substring(0, 100));
       }
       document.head.appendChild(newScript);
     });
   }
   
   /*************************************************************
    * Screening Helpers
    *************************************************************/
   function runScreening() {
     console.log("runScreening() is called");
   
     const regionChr = document.getElementById("region_chr1").value;
     const regionStartVal = parseInt(document.getElementById("region_start1").value, 10);
     let regionEndVal   = parseInt(document.getElementById("region_end1").value, 10);
   
     if (isNaN(regionEndVal)) {
       regionEndVal = regionStartVal + WINDOW_WIDTH;
       document.getElementById("region_end1").value = regionEndVal;
     }
   
     const screeningContainerElem = document.getElementById("screening_chart");
     if (screeningContainerElem) {
       screeningContainerElem.innerHTML = "<div class='loader' style='display:block;margin:0 auto;'></div>";
     }
   
     const paramsObj = {
       region_chr:   regionChr,
       model_select: document.getElementById("model_select").value,
       region_start: regionStartVal,
       region_end:   regionEndVal,
       perturb_width: document.getElementById("perturb_width").value,
       ctcf_bw_path:  document.getElementById("ctcf_bw_path").value,
       atac_bw_path:  document.getElementById("atac_bw_path").value,
       output_dir:    document.getElementById("output_dir").value,
       peaks_file:    document.getElementById("peaks_file_path").value
     };
     const queryParams = new URLSearchParams(paramsObj).toString();
     console.log("Sending screening request with parameters:", paramsObj);
   
     const xhr = new XMLHttpRequest();
     xhr.onreadystatechange = function() {
       if (xhr.readyState === 4) {
         if (xhr.status === 200) {
           console.log("Screening request successful. Response:", xhr.responseText);
           const response = JSON.parse(xhr.responseText);
           if (response.screening_config) {
             const screeningConfig = response.screening_config;
             if (screeningContainerElem) {
               screeningContainerElem.innerHTML = "";
               // drawColumnChart defined in d3_plots.js
               drawColumnChart('#screening_chart', screeningConfig);
             }
           } else {
             if (screeningContainerElem) {
               screeningContainerElem.innerHTML = "<p>No screening config returned.</p>";
             }
           }
         } else {
           console.error("Error generating screening plot. Status:", xhr.status);
           if (screeningContainerElem) {
             screeningContainerElem.innerHTML = "<p>Error generating screening plot.</p>";
           }
         }
       }
     };
     xhr.open("GET", "/run_screening?" + queryParams, true);
     xhr.send();
   }
   
   /*************************************************************
    * Checking Required Fields
    *************************************************************/
   function checkRequiredFieldsFilled() {
     const dsOption = document.querySelector('input[name="ds_option"]:checked');
     if (!dsOption) {
       console.log("No ds_option selected; disabling button");
       return false;
     }
   
     // Basic required fields
     const reqIds = [
       "region_chr1",
       "region_start1",
       "model_select",
       "genome_select",
       "atac_bw_path"
     ];
   
     if (dsOption.value === "deletion") {
       reqIds.push("del_start", "del_width");
     } else if (dsOption.value === "screening") {
       reqIds.push("perturb_width");
     }
   
     for (let id of reqIds) {
       const elem = document.getElementById(id);
       if (!elem) {
         console.log(`Required field #${id} not found in DOM`);
         return false;
       }
       if (!elem.value || !elem.value.trim()) {
         console.log(`Required field #${id} is empty -> '${elem.value}'`);
         return false;
       }
     }
   
     console.log("All required fields are filled!");
     return true;
   }
   
   /*************************************************************
    * Document Ready
    *************************************************************/
   document.addEventListener('DOMContentLoaded', function() {
     console.log("DOM fully loaded");
   
     // 1) Kick off all fetches in parallel
     const pChr  = populateChromosomeDropdown(); // Chromosome
     const pAtac = populateDropdownFromServer("atac_bw_path", "atac");
     const pCtcf = populateDropdownFromServer("ctcf_bw_path", "ctcf");
     const pPeaks= populateDropdownFromServer("peaks_file_path", "peaks");
   
     // 2) Wait for *all* promises to resolve
     Promise.all([pChr, pAtac, pCtcf, pPeaks])
       .then(() => {
         // Now region_chr1 / region_chr2 are filled,
         // and atac_bw_path, ctcf_bw_path, peaks_file_path are also filled
   
         // Only now do we do final init:
         restoreFormFields();
         updateNormalizationLabels();
         toggleOptionalFields();
         updateEndPosition();
         updateNormLocks();
   
         // Mark init complete & validate
         window.initComplete = true;
         collectAndDisplayErrors();
         updateNormLocks();
   
         console.log("All dropdowns done. Initialization complete!");
       })
       .catch(err => {
         console.error("Populate calls failed:", err);
       });
   
     // Restore second region on reload
     const dsOption = document.querySelector('input[name="ds_option"]:checked');
     const chimericVal = localStorage.getItem("chimeric_active");
     if (dsOption && dsOption.value === "deletion") {
       hideSecondChr();
     } else {
       if (chimericVal === "true") {
         showSecondChr();
       } else {
         hideSecondChr();
       }
     }
   
     // Region 1 changes
     const rs1 = document.getElementById("region_start1");
     if (rs1) {
       rs1.addEventListener("input", () => {
         storeFormFields();
         updateEndPosition();
         collectAndDisplayErrors();
       });
     }
     const re1 = document.getElementById("region_end1");
     if (re1) {
       re1.addEventListener("input", () => {
         updateEndPosition();
         collectAndDisplayErrors();
       });
     }
     const rc1 = document.getElementById("region_chr1");
     if (rc1) {
       rc1.addEventListener("change", () => {
         storeFormFields();
         collectAndDisplayErrors();
       });
     }
   
     // Region 2 changes
     const startField2 = document.getElementById("region_start2");
     if (startField2) {
       startField2.addEventListener("input", () => {
         storeFormFields();
         updateEndPosition();
         collectAndDisplayErrors();
       });
     }
     const re2 = document.getElementById("region_end2");
     if (re2) {
       re2.addEventListener("input", () => {
         updateEndPosition();
         collectAndDisplayErrors();
       });
     }
   
     // Deletion fields
     const delStart = document.getElementById("del_start");
     if (delStart) {
       delStart.addEventListener("input", () => {
         storeFormFields();
         collectAndDisplayErrors();
       });
     }
     const delWidth = document.getElementById("del_width");
     if (delWidth) {
       delWidth.addEventListener("input", () => {
         storeFormFields();
         collectAndDisplayErrors();
       });
     }
   
     // Screening fields
     const pw = document.getElementById("perturb_width");
     if (pw) pw.addEventListener("input", storeFormFields);
   
     // Model changes
     const modelSel = document.getElementById("model_select");
     if (modelSel) {
       modelSel.addEventListener("change", () => {
         updateNormalizationLabels();
         updateNormLocks();
         const ctcfDD = document.getElementById("ctcf_bw_path");
         if (ctcfDD) ctcfDD.dispatchEvent(new Event("change", { bubbles: true }));
         collectAndDisplayErrors();
         storeFormFields();
       });
     }
   
     // DS option toggles
     const dsRadios = document.getElementsByName("ds_option");
     dsRadios.forEach(radio => {
       radio.addEventListener("change", () => {
         storeFormFields();
         toggleOptionalFields();
         collectAndDisplayErrors();
         updateNormLocks();
       });
     });
   
     // Second chromosome toggle
     const togglePlusElemContainer = document.getElementById("toggle-second-chr-container");
     if (togglePlusElemContainer) {
       togglePlusElemContainer.addEventListener("click", toggleSecondChr);
     } else {
       console.error("Element with id 'toggle-second-chr-container' not found");
     }
   
     // File uploads
     const atacFile = document.getElementById("atac_bw_file");
     if (atacFile) atacFile.addEventListener("change", () => {
       ajaxUploadFile("atac_bw_file", "atac", "atac_bw_path");
     });
     const ctcfFile = document.getElementById("ctcf_bw_file");
     if (ctcfFile) ctcfFile.addEventListener("change", () => {
       ajaxUploadFile("ctcf_bw_file", "ctcf", "ctcf_bw_path");
     });
     const peaksFile = document.getElementById("peaks_file_upload");
     if (peaksFile) peaksFile.addEventListener("change", () => {
       ajaxUploadFile("peaks_file_upload", "peaks", "peaks_file_path");
     });
   
     // Normalization labels / locks (initial sync)
     updateNormalizationLabels();
     updateNormLocks();
   
     /********** AWS/NYU guard **********/
     // Server-side flag rendered into <body data-aws-enabled="...">
     const serverFlag = document.body && document.body.dataset
       ? document.body.dataset.awsEnabled === "true"
       : false;
   
     // Hard guard: treat NYU hosts as non-AWS no matter what the server flag says
     const onNYU =
       /^10\.185\./.test(location.hostname) ||
       location.hostname.endsWith(".nyumc.org");
   
     // Optional emergency kill switch (set via DevTools: localStorage.setItem('DISABLE_WORKER','1'))
     const killSwitch = localStorage.getItem("DISABLE_WORKER") === "1";
   
     // Final truth: only AWS when server says so AND we're not on NYU AND no kill switch
     const awsEnabled = serverFlag && !onNYU && !killSwitch;
   
     /********** Submit handler wiring **********/
     const formElem = document.getElementById("corigami-form");
     let runInProgress = false;
   
     if (formElem) {
       formElem.addEventListener("submit", async (e) => {
         e.preventDefault();
         if (runInProgress) return;
         if (!checkFormRequirements()) return;
   
         runInProgress = true;
         const container = document.getElementById("output-container");
         if (container) {
           container.innerHTML = '<div class="loader" style="display:block;margin:0 auto;"></div>';
         }
   
         // wipe 2nd-region inputs if the toggle is off
         const toggle = document.getElementById("toggle-second-chr");
         if (!toggle || !toggle.classList.contains("active")) {
           ["region_chr2", "region_start2", "region_end2"].forEach(id => {
             const el = document.getElementById(id);
             if (el) el.value = "";
           });
         }
   
         let job_id;
         try {
           const payload = Object.fromEntries(new FormData(formElem));
           const resp = await fetch("/api/predict", {
             method: "POST",
             headers: { "Content-Type": "application/json" },
             body: JSON.stringify(payload),
           });
           if (resp.status !== 202) throw new Error(`HTTP ${resp.status}`);
           const j = await resp.json();
           job_id = j.job_id;
         } catch (err) {
           console.error("predict launch failed:", err);
           if (container) container.innerHTML = "<p style='color:red;'>Could not start job.</p>";
           runInProgress = false;
           return;
         }
   
         async function waitUntilDone(id) {
           // Poll status until finished/failed
           while (true) {
             const res = await fetch(`/api/job/${id}`);
             if (!res.ok) throw new Error(`HTTP ${res.status}`);
             const { status } = await res.json();
             if (status === "failed") throw new Error("worker failed");
             if (status === "finished") return;
             // wait 1.5s
             await new Promise(r => setTimeout(r, 1500));
           }
         }
   
         try {
           await waitUntilDone(job_id);
         } catch (err) {
           console.error("job failed:", err);
           if (container) container.innerHTML = "<p style='color:red;'>Prediction job failed. Check logs.</p>";
           runInProgress = false;
           return;
         }
   
         try {
           const html = await fetch(`/api/job/${job_id}/html`).then(r => r.text());
           if (container) {
             container.innerHTML = html;
             executeScripts(container);
           }
         } catch (err) {
           console.error("partial fetch failed:", err);
           if (container) container.innerHTML = "<p style='color:red;'>Could not fetch results.</p>";
         } finally {
           runInProgress = false;
         }
       });
     } else {
       console.error("Form #corigami-form not found; submit handler not attached.");
     }
   
     // ───────── AWS “Connect to server” button logic (honors awsEnabled) ─────────
     (function workerButtonBlock() {
       const btn = document.getElementById("worker-btn");
       // Helper: add/remove CSS class that makes the button visible
       const showBtn = (show) => {
         if (!document.body) return;
         document.body.classList.toggle("show-worker", !!show);
       };
   
       if (!btn) return; // If there is no button in DOM, nothing to do.
   
       if (!awsEnabled) {
         // Not on AWS → remove button and bail (prevents polling + console spam)
         btn.remove();
         showBtn(false);
         return;
       }
   
       // On AWS: wire up polling logic
       const BOOT_TIMEOUT_MS = 10 * 60 * 1000; // 10 minutes
       let startInitiatedAt = null;
   
       // Let CSS control visibility; stop using hidden attribute for control
       btn.hidden = false;
       btn.setAttribute("aria-hidden", "false");
   
       function setStatus(color /* 'red'|'orange'|'green' */, text, disabled) {
         btn.classList.remove("status-red", "status-orange", "status-green");
         btn.classList.add(`status-${color}`);
         btn.textContent = text;
         btn.disabled = !!disabled;
       }
   
       async function checkOnce() {
         try {
           const res = await fetch("/api/worker-status");
           if (!res.ok) {
             // Non-200: keep button visible but don't crash the loop
             showBtn(true);
             setStatus("red", "Status error", false);
             return 10000;
           }
           const { running, starting } = await res.json();
           const now = Date.now();
   
           if (running) {
             // server up → hide button, slow poll
             showBtn(false);
             return 60000;
           }
   
           // not running → show button
           showBtn(true);
   
           if (starting || (startInitiatedAt && now - startInitiatedAt < BOOT_TIMEOUT_MS)) {
             setStatus("orange", "Server starting…", true);
             return 5000;
           }
   
           if (startInitiatedAt && now - startInitiatedAt >= BOOT_TIMEOUT_MS) {
             startInitiatedAt = null;
             setStatus("red", "Connection error. Retry?", false);
             return 10000;
           }
   
           setStatus("red", "Connect to server", false);
           return 10000;
         } catch (e) {
           // Network/JSON error → show button but stay quiet-ish
           console.error("status poll failed", e);
           showBtn(true);
           setStatus("red", "Status error", false);
           return 10000;
         }
       }
   
       async function poll() {
         const next = await checkOnce();
         setTimeout(poll, next);
       }
       poll();
   
       btn.addEventListener("click", async () => {
         if (btn.disabled) return;
         startInitiatedAt = Date.now();
         showBtn(true);
         setStatus("orange", "Server starting…", true);
         try {
           const resp = await fetch("/api/start-worker", { method: "POST" });
           if (![200, 202].includes(resp.status)) throw new Error(`status ${resp.status}`);
           // polling loop will update UI
         } catch (e) {
           console.error("start-worker failed", e);
           startInitiatedAt = null;
           setStatus("red", "Start failed", false);
         }
       });
     })();
   }); // DOMContentLoaded end
   
   /*************************************************************
    * Screening warnings
    *************************************************************/
   function checkScreeningRegionEdges() {
     const dsOption = document.querySelector('input[name="ds_option"]:checked');
     if (!dsOption || dsOption.value !== "screening") {
       return null;
     }
   
     const secondChrDiv = document.getElementById("second_chr_fields");
     if (secondChrDiv && secondChrDiv.style.display !== "none") {
       return {
         isWarning: true,
         message: "When providing two regions, screening will skip the first and last 1 Mb of the combined region."
       };
     }
   
     const chrSelect = document.getElementById("region_chr1");
     const startElem = document.getElementById("region_start1");
     const endElem   = document.getElementById("region_end1");
     if (!chrSelect || !window.chromosomeLengths) return null;
   
     const selectedChr = chrSelect.value;
     const chrLength   = window.chromosomeLengths[selectedChr];
     if (typeof chrLength === "undefined") return null;
   
     const start = parseInt(startElem.value, 10);
     const defaultEnd = start + WINDOW_WIDTH;
     const end   = endElem.value ? parseInt(endElem.value, 10) : defaultEnd;
   
     const buffer = WINDOW_WIDTH / 2; // 1 Mb
     if (start < buffer || end > chrLength - buffer) {
       return {
         isWarning: true,
         message: "Screening will not occur within 1 Mb of the start or end of the chromosome."
       };
     }
     return null;
   }
   
   function checkMinimumRegionSize() {
     const TWO_MB = 2097152;
   
     const start1 = parseInt(document.getElementById("region_start1").value, 10);
     let end1   = parseInt(document.getElementById("region_end1").value, 10);
     if (isNaN(start1)) return null;
   
     if (isNaN(end1)) {
       end1 = start1 + TWO_MB;
     }
     const size1 = end1 - start1;
   
     const secondChrDiv = document.getElementById("second_chr_fields");
     if (secondChrDiv && secondChrDiv.style.display !== "none") {
       const start2 = parseInt(document.getElementById("region_start2").value, 10);
       let end2   = parseInt(document.getElementById("region_end2").value, 10);
       if (isNaN(start2)) return null;
   
       if (isNaN(end2)) {
         end2 = start2 + TWO_MB;
       }
       const size2 = end2 - start2;
       if (size1 < TWO_MB || size2 < TWO_MB) {
         return {
           isWarning: true,
           message: "For best results, please ensure each provided region is greater than 2 Mb in length."
         };
       }
     } else {
       if (size1 < TWO_MB) {
         return {
           isWarning: true,
           message: "For best results, please ensure each provided region is greater than 2 Mb in length."
         };
       }
     }
     return null;
   }
   
   /*************************************************************
    * Normalization UI helpers
    *************************************************************/
   function updateNormalizationLabels() {
     const modelVal = document.getElementById("model_select").value;
     const atacLabel = document.getElementById("apply_atac_norm_label");
     const ctcfLabel = document.getElementById("apply_ctcf_norm_label");
     const ctcfRow = document.getElementById("ctcf_norm_row");
     const reqDiv = document.getElementById("model-requirements-list");
     const predictCtcfRow = document.getElementById("predict_ctcf_row");
     const ctcfDropdown = document.getElementById("ctcf_bw_path");
     const noneOption = ctcfDropdown ? ctcfDropdown.querySelector('option[value="none"]') : null;
   
     if (!atacLabel || !ctcfLabel || !ctcfRow || !reqDiv || !predictCtcfRow || !ctcfDropdown || !noneOption) return;
   
     // Basic label/row toggles
     if (modelVal === "IMR90") {
       ctcfRow.style.display = "none";
       predictCtcfRow.style.display = "none";
       atacLabel.innerHTML = "Apply <strong>log norm</strong> to my raw ATAC file";
       reqDiv.innerHTML = "This model requires a <strong>log norm</strong> ATAC file and a <strong>log2FC norm</strong> CTCF file.";
       noneOption.style.display = "none";
       if (ctcfDropdown.value === "none") {
         ctcfDropdown.value = "./corigami_data/data/hg38/imr90/genomic_features/ctcf_log2fc.bw";
       }
     } else if (modelVal === "BALL") {
       ctcfRow.style.display = "";
       predictCtcfRow.style.display = "";
       atacLabel.innerHTML = "Apply <strong>minmax norm</strong> to my raw ATAC file";
       ctcfLabel.innerHTML = "Apply <strong>minmax norm</strong> to my raw CTCF file";
       reqDiv.innerHTML = "This model is trained on expected/observed data and requires <strong>minmax normalized</strong> ATAC and CTCF files.";
       noneOption.style.display = "";
     } else {
       ctcfRow.style.display = "";
       predictCtcfRow.style.display = "";
       atacLabel.innerHTML = "Apply normalization to my raw ATAC file";
       ctcfLabel.innerHTML = "Apply normalization to my raw CTCF file";
       reqDiv.innerHTML = "";
       noneOption.style.display = "";
     }
   
     // Helper: pick the first option that has data-model = modelVal
     function fallbackToMatchingPreset(dropdown, modelValue) {
       const matching = Array.from(dropdown.options).find(
         o => o.getAttribute("data-model") &&
              o.getAttribute("data-model").toLowerCase() === modelValue.toLowerCase()
       );
       if (matching) {
         matching.style.display = "";
         dropdown.value = matching.value;
       } else {
         dropdown.value = "none";
       }
     }
   
     // Hide unrelated preset options in ATAC, CTCF, and peaks
     const atacDropdown  = document.getElementById("atac_bw_path");
     const peaksDropdown = document.getElementById("peaks_file_path");
     [atacDropdown, ctcfDropdown, peaksDropdown].forEach(dd => {
       if (!dd) return;
       const opts = Array.from(dd.options);
       opts.forEach(opt => {
         const dm = opt.getAttribute("data-model");
         if (!dm) {
           opt.style.display = ""; // universal
         } else if (dm.toLowerCase() === modelVal.toLowerCase()) {
           opt.style.display = "";
         } else {
           opt.style.display = "none";
           if (dd.value === opt.value) {
             // fallback to the first matching preset or "none"
             fallbackToMatchingPreset(dd, modelVal);
           }
         }
       });
     });
   }
   
   /*************************************************************
    * Normalization lock helpers
    *************************************************************/
   function _selectedIsPreset(selectEl) {
     if (!selectEl) return false;
     const opt = selectEl.options[selectEl.selectedIndex];
     if (!opt) return false;
   
     // 1) official signal
     if (opt.hasAttribute("data-model")) return true;
   
     // 2) value-based fallback for presets shipped with the app
     const v = String(opt.value || "").toLowerCase();
     return v.includes("/corigami_data/") || v.startsWith("./corigami_data/");
   }
   
   function _lock(wrapperEl, checkboxEl, shouldCheck) {
     if (!wrapperEl || !checkboxEl) return;
     checkboxEl.checked = !!shouldCheck;
     wrapperEl.classList.toggle("checked", !!shouldCheck);
     wrapperEl.classList.add("disabled-toggle", "locked"); // visual + block clicks
   }
   function _unlock(wrapperEl) {
     if (!wrapperEl) return;
     wrapperEl.classList.remove("disabled-toggle", "locked");
   }
   
   function updateNormLocks() {
     const model  = document.getElementById("model_select")?.value || "IMR90";
   
     const atacDD = document.getElementById("atac_bw_path");
     const ctcfDD = document.getElementById("ctcf_bw_path");
   
     const atacCB = document.getElementById("apply_atac_norm");
     const atacW  = document.getElementById("apply_atac_norm_wrapper");
   
     const ctcfCB = document.getElementById("apply_ctcf_norm");
     const ctcfW  = document.getElementById("apply_ctcf_norm_wrapper");
   
     const autoCB = document.getElementById("predict_ctcf");
     const autoW  = document.getElementById("predict_ctcf_wrapper");
   
     // reset → then apply locks
     _unlock(atacW); _unlock(ctcfW); _unlock(autoW);
   
     const atacIsPreset = _selectedIsPreset(atacDD);
     const ctcfIsPreset = _selectedIsPreset(ctcfDD);
     const ctcfIsAuto   = !!(autoCB && autoCB.checked) || (ctcfDD && ctcfDD.value === "none");
   
     // ATAC
     if (atacIsPreset && model === "IMR90") {
       // IMR-90 + preset ATAC → log-norm = CHECKED + locked
       _lock(atacW, atacCB, true);
     } else if (atacIsPreset && model === "BALL") {
       // B-ALL + preset ATAC → minmax = UNCHECKED + locked
       _lock(atacW, atacCB, false);
     }
   
     // CTCF (B-ALL rules)
     if (model === "BALL" && ctcfIsPreset) {
       _lock(ctcfW, ctcfCB, false); // preset CTCF → minmax unchecked + locked
     }
     if (model === "BALL" && ctcfIsAuto) {
       _lock(autoW, autoCB, true);  // auto-generate locked ON
       _lock(ctcfW, ctcfCB, true);  // norm locked ON
     }
   }
   