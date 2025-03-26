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
  // Save first region and other fields
  const fields = [
    "region_chr1", "region_start1", 
    /* remove region_end1 from this list! */
    "del_start", "del_width", 
    "perturb_width", "step_size", 
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

  const normAtacRadio = document.querySelector('input[name="norm_atac"]:checked');
  if (normAtacRadio) localStorage.setItem("norm_atac", normAtacRadio.value);

  const normCtcfRadio = document.querySelector('input[name="norm_ctcf"]:checked');
  if (normCtcfRadio) localStorage.setItem("norm_ctcf", normCtcfRadio.value);

  const dsRadios = document.getElementsByName("ds_option");
  for (let radio of dsRadios) {
    if (radio.checked) {
      localStorage.setItem("ds_option", radio.value);
      break;
    }
  }
}

function restoreFormFields() {
  // Fields for first region and other parameters (excluding region_end1!)
  const fields = [
    "region_chr1", "region_start1",
    "del_start", "del_width", 
    "perturb_width", "step_size", 
    "model_select", "atac_bw_path", 
    "ctcf_bw_path", "peaks_file_path", 
    "genome_select"
  ];
  fields.forEach(field => {
    const storedVal = localStorage.getItem(field);
    const elem = document.getElementById(field);
    if (storedVal && elem) {
      elem.value = storedVal;
    }
  });

  // Restore region_end1 ONLY if user explicitly typed a value
  const storedEnd1 = localStorage.getItem("region_end1");
  const endElem1 = document.getElementById("region_end1");
  if (storedEnd1 && endElem1) {
    endElem1.value = storedEnd1;
  } else {
    // user never typed it, so remain blank and show placeholder
    endElem1.value = "";
  }

  // Same idea for region_end2
  const storedEnd2 = localStorage.getItem("region_end2");
  const endElem2 = document.getElementById("region_end2");
  if (storedEnd2 && endElem2) {
    endElem2.value = storedEnd2;
  } else {
    endElem2.value = "";
  }

  // Restore region_chr2 and region_start2
  const secondFields = ["region_chr2", "region_start2"];
  secondFields.forEach(field => {
    const storedVal = localStorage.getItem(field);
    const elem = document.getElementById(field);
    if (storedVal && elem) {
      elem.value = storedVal;
    }
  });

  // Norm fields
  const normAtac = localStorage.getItem("norm_atac");
  if (normAtac) {
    const radio = document.querySelector(`input[name="norm_atac"][value="${normAtac}"]`);
    if (radio) radio.checked = true;
  }
  const normCtcf = localStorage.getItem("norm_ctcf");
  if (normCtcf) {
    const radio = document.querySelector(`input[name="norm_ctcf"][value="${normCtcf}"]`);
    if (radio) radio.checked = true;
  }

  // DS option
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
  // This now *only* populates chromosomes and returns the promise,
  // no calls to restoreFormFields() or collectAndDisplayErrors().
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
          let cleanChr = chr.replace(/^chr/i, "");
          let option = document.createElement("option");
          option.value = chr; // e.g. "chr1"
          option.text = `${cleanChr} (${data[chr]})`;
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
  const startVal1   = parseInt(startField1.value);
  if (!isNaN(startVal1)) {
    const computedEnd1 = startVal1 + WINDOW_WIDTH;
    if (!endField1.value || isNaN(parseInt(endField1.value))) {
      endField1.placeholder = computedEnd1;
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
    const startVal2   = parseInt(startField2.value);
    if (!isNaN(startVal2)) {
      const computedEnd2 = startVal2 + WINDOW_WIDTH;
      if (!endField2.value || isNaN(parseInt(endField2.value))) {
        endField2.placeholder = computedEnd2;
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

  if (!chrSelect || !window.chromosomeLengths) return null;

  const selectedChr = chrSelect.value;
  const maxLength = window.chromosomeLengths[selectedChr];
  if (typeof maxLength === "undefined") return null;

  const start = parseInt(startElem.value);
  const end   = endElem.value ? parseInt(endElem.value) : start + WINDOW_WIDTH;

  if (end < start) {
    return "Region 1: Chromosome start position cannot be greater than end";
  }
  if (start < 0 || end < 0 || start > maxLength || end > maxLength) {
    return `Region 1: ${selectedChr} must be between 0 and ${maxLength}`;
  }
  if (end - start > WINDOW_WIDTH * 10) {
    return "Region 1: window must be under 20 Mb in length";
  }
  return null;
}

function validateRegion2Bounds() {
  const secondChrDiv = document.getElementById("second_chr_fields");
  if (!secondChrDiv || secondChrDiv.style.display === "none") {
    return null;
  }
  const chrSelect2  = document.getElementById("region_chr2");
  const startElem2  = document.getElementById("region_start2");
  const endElem2    = document.getElementById("region_end2");
  if (!chrSelect2 || !window.chromosomeLengths) return null;

  const selectedChr2 = chrSelect2.value;
  const maxLength2   = window.chromosomeLengths[selectedChr2];
  if (typeof maxLength2 === "undefined") return null;

  const start2 = parseInt(startElem2.value);
  const end2   = endElem2.value ? parseInt(endElem2.value) : start2 + WINDOW_WIDTH;

  if (end2 < start2) {
    return "Region 2: Chromosome start position cannot be greater than end";
  }
  if (start2 < 0 || end2 < 0 || start2 > maxLength2 || end2 > maxLength2) {
    return `Region 2: ${selectedChr2} must be between 0 and ${maxLength2}`;
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

  const regionStart   = parseInt(regionStartElem.value);
  const regionEnd     = regionEndElem.value ? parseInt(regionEndElem.value) : regionStart + WINDOW_WIDTH;
  const deletionStart = parseInt(delStartElem.value);
  const deletionWidth = parseInt(delWidthElem.value);
  const deletionEnd   = deletionStart + deletionWidth;

  if (deletionStart < regionStart || deletionEnd > regionEnd) {
    return "Deletion area is out of bounds for the selected region.";
  }
  return null;
}

function checkCtcfAtacRequirements() {
  const ctcfFile = document.getElementById("ctcf_bw_path").value;
  const atacState = document.querySelector('input[name="norm_atac"]:checked').value;
  if (ctcfFile === "none" && (atacState === "log" || atacState === "minmax")) {
    return "Raw ATAC BigWig is needed when CTCF is not provided.";
  }
  return null;
}

function checkNormalizationMismatch() {
  const atacState = document.querySelector('input[name="norm_atac"]:checked').value;
  const ctcfState = document.querySelector('input[name="norm_ctcf"]:checked').value;
  const ctcfFile  = document.getElementById("ctcf_bw_path").value;

  if (ctcfFile !== "none" &&
      atacState !== "none" &&
      ctcfState !== "none" &&
      atacState !== ctcfState) {
    return "ATAC & CTCF files cannot already have two different normalization methods.";
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

  // -- Already existing checks --
  const eReg1 = validateRegion1Bounds();
  if (eReg1) errors.push(eReg1);

  const eReg2 = validateRegion2Bounds();
  if (eReg2) errors.push(eReg2);

  const eDel = validateDeletionAreaErrors();
  if (eDel) errors.push(eDel);

  const eCtcf = checkCtcfAtacRequirements();
  if (eCtcf) errors.push(eCtcf);

  const eNormMismatch = checkNormalizationMismatch();
  if (eNormMismatch) errors.push(eNormMismatch);

  const eTotal = validateTotalRegionSize();
  if (eTotal) errors.push(eTotal);

  // -- NEW: Check screening region edges --
  const eScreenEdges = checkScreeningRegionEdges(); 
  if (eScreenEdges) {
    if (eScreenEdges.isWarning) {
      warnings.push(eScreenEdges.message);
    } else if (eScreenEdges.isError) {
      errors.push(eScreenEdges.message);
    }
  }

  // -- NEW function for region minimum size
  const eMinSize = checkMinimumRegionSize();
  if (eMinSize && eMinSize.isWarning) {
    warnings.push(eMinSize.message);
  }

  // Clear old content
  errorContainer.innerHTML = "";

  if (errors.length === 0 && warnings.length === 0) {
    errorContainer.style.display = "none";
  } else {
    errorContainer.style.display = "block";

    // 1) Show blocking errors
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

    // 2) Show non-blocking warnings
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
  if (errors.length > 0 || !checkRequiredFieldsFilled()) {
    submitBtn.disabled = true;
    submitBtn.classList.add("disabled");
  } else {
    submitBtn.disabled = false;
    submitBtn.classList.remove("disabled");
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

  secondChrDiv.style.display = "block";
  togglePlus.classList.add("active");
  chimericInput.value = "true";
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
}

function hideSecondChr() {
  const secondChrDiv  = document.getElementById('second_chr_fields');
  const togglePlus    = document.getElementById('toggle-second-chr');
  const chimericInput = document.getElementById('chimeric_active');
  const peaksInput    = document.getElementById("peaks_file_upload");
  const peaksDropdown = document.getElementById("peaks_file_path");

  secondChrDiv.style.display = "none";
  togglePlus.classList.remove("active");
  chimericInput.value = "false";
  localStorage.setItem("chimeric_active", "false");

  if (peaksInput) {
    peaksInput.disabled = false;
  }

  if (peaksDropdown && peaksDropdown.dataset.originalOptions) {
    peaksDropdown.innerHTML = peaksDropdown.dataset.originalOptions;
    peaksDropdown.disabled = false;
  }

  collectAndDisplayErrors();
}

function toggleOptionalFields() {
  const dsOption        = document.querySelector('input[name="ds_option"]:checked');
  const delFields       = document.getElementById("deletion-fields");
  const scrFields       = document.getElementById("screening-fields-2");
  const peaksContainer  = document.getElementById("peaks_file_container");
  const togglePlusElemContainer = document.getElementById("toggle-second-chr-container");

  if (!dsOption) return;

  if (dsOption.value === "deletion") {
    if (delFields) delFields.style.display = "block";
    if (scrFields) scrFields.style.display = "none";
    if (peaksContainer) peaksContainer.style.display = "none";
    hideSecondChr(); // always hide in deletion
    if (togglePlusElemContainer) togglePlusElemContainer.style.display = "none";
  } else if (dsOption.value === "screening") {
    if (scrFields) scrFields.style.display = "block";
    if (delFields) delFields.style.display = "none";
    if (peaksContainer) peaksContainer.style.display = "block";
    if (togglePlusElemContainer) togglePlusElemContainer.style.display = "block";
  } else {
    if (delFields) delFields.style.display = "none";
    if (scrFields) scrFields.style.display = "none";
    if (peaksContainer) peaksContainer.style.display = "none";
    if (togglePlusElemContainer) togglePlusElemContainer.style.display = "block";
  }
}

function toggleSecondChr() {
  const dsOption = document.querySelector('input[name="ds_option"]:checked');
  if (dsOption && dsOption.value === "deletion") {
    hideSecondChr();
  } else {
    const secondChrDiv = document.getElementById('second_chr_fields');
    if (secondChrDiv.style.display === "none" || secondChrDiv.style.display === "") {
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
    if (!document.getElementById("perturb_width").value.trim() || 
        !document.getElementById("step_size").value.trim()) {
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
  console.log("Selected file for " + fileType + ":", file.name);

  const uploadButton = document.querySelector("label[for='" + fileInputId + "']");
  if (uploadButton) {
    uploadButton.textContent = 'Uploading...';
    uploadButton.classList.add('loading');
  }

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
    } else {
      console.log("Upload successful for " + fileType + ":", data);
      const dropdown = document.getElementById(dropdownId);
      const newOption = document.createElement("option");
      newOption.value = data.saved_path;
      newOption.text = data.display_name;
      dropdown.appendChild(newOption);
      dropdown.value = data.saved_path;
      localStorage.setItem(dropdownId, data.saved_path);
    }
  })
  .catch(error => console.error("Error uploading file for " + fileType + ":", error))
  .finally(() => {
    if (uploadButton) {
      uploadButton.textContent = 'Upload';
      uploadButton.classList.remove('loading');
    }
    collectAndDisplayErrors();
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
 * Screening and Normalization Helpers
 *************************************************************/
function runScreening() {
  console.log("runScreening() is called");

  let regionChr = document.getElementById("region_chr1").value;
  let regionStartVal = parseInt(document.getElementById("region_start1").value, 10);
  let regionEndVal   = parseInt(document.getElementById("region_end1").value, 10);

  if (isNaN(regionEndVal)) {
    regionEndVal = regionStartVal + 2097152;
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
    step_size:     document.getElementById("step_size").value,
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
          const screeningConfig = JSON.parse(response.screening_config);
          if (screeningContainerElem) {
            screeningContainerElem.innerHTML = "";
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

function updateCtcfNormalization() {
  const ctcfSelect = document.getElementById('ctcf_bw_path');
  const ctcfSwitch = document.getElementById('ctcf-switch');
  const ctcfNoNorm = document.getElementById('ctcf-no-norm');
  const ctcfLog    = document.getElementById('ctcf-log');
  const ctcfMinmax = document.getElementById('ctcf-minmax');

  if (ctcfSelect.value === 'none') {
    ctcfNoNorm.disabled = true;
    ctcfLog.disabled    = true;
    ctcfMinmax.disabled = true;
    ctcfSwitch.classList.add('disabled-toggle');
  } else {
    ctcfNoNorm.disabled = false;
    ctcfLog.disabled    = false;
    ctcfMinmax.disabled = false;
    ctcfSwitch.classList.remove('disabled-toggle');
  }
}

function updateTrainingNormField() {
  const ctcfFile   = document.getElementById("ctcf_bw_path").value;
  const atacState  = document.querySelector('input[name="norm_atac"]:checked').value;
  const ctcfState  = document.querySelector('input[name="norm_ctcf"]:checked').value;
  const trainingContainer = document.getElementById("training-norm-container");

  if (ctcfFile === "none") {
    document.getElementById("training-minmax").checked = true;
    document.getElementById("training-log").disabled   = true;
    trainingContainer.classList.add("disabled-toggle");
  } else if (atacState !== "none") {
    if (atacState === "log") {
      document.getElementById("training-log").checked = true;
    } else {
      document.getElementById("training-minmax").checked = true;
    }
    document.getElementById("training-log").disabled   = true;
    document.getElementById("training-minmax").disabled= true;
    trainingContainer.classList.add("disabled-toggle");
  } else if (ctcfState !== "none") {
    if (ctcfState === "log") {
      document.getElementById("training-log").checked = true;
    } else {
      document.getElementById("training-minmax").checked = true;
    }
    document.getElementById("training-log").disabled   = true;
    document.getElementById("training-minmax").disabled= true;
    trainingContainer.classList.add("disabled-toggle");
  } else {
    document.getElementById("training-log").disabled    = false;
    document.getElementById("training-minmax").disabled = false;
    trainingContainer.classList.remove('disabled-toggle');
  }

  collectAndDisplayErrors();
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
    reqIds.push("perturb_width", "step_size");
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
      toggleOptionalFields();
      updateCtcfNormalization();
      updateTrainingNormField();
      updateEndPosition();

      // Mark init complete & validate
      window.initComplete = true;
      collectAndDisplayErrors();

      console.log("All dropdowns done. Initialization complete!");
    })
    .catch(err => {
      console.error("Populate calls failed:", err);
    });

  // ========== Restore second region on reload (shortened) ==========
  // (We do this after ds_option is set in localStorage,
  //  but if you rely on region_chr1 <option>, remember
  //  they are only guaranteed after the fetch. 
  //  That's why we do final checks in .then(...).)
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
  // ================================================================

  // Region 1 changes
  document.getElementById("region_start1").addEventListener("input", () => {
    storeFormFields();
    updateEndPosition();
    collectAndDisplayErrors();
  });
  document.getElementById("region_end1").addEventListener("input", () => {
    updateEndPosition();
    collectAndDisplayErrors();
  });
  document.getElementById("region_chr1").addEventListener("change", () => {
    storeFormFields();
    collectAndDisplayErrors();
  });

  // Region 2 changes
  const startField2 = document.getElementById("region_start2");
  if (startField2) {
    startField2.addEventListener("input", () => {
      storeFormFields();
      updateEndPosition();
      collectAndDisplayErrors();
    });
  }
  document.getElementById("region_end2").addEventListener("input", () => {
    updateEndPosition();
    collectAndDisplayErrors();
  });

  // Deletion fields
  document.getElementById("del_start").addEventListener("input", () => {
    storeFormFields();
    collectAndDisplayErrors();
  });
  document.getElementById("del_width").addEventListener("input", () => {
    storeFormFields();
    collectAndDisplayErrors();
  });

  // Screening fields
  document.getElementById("perturb_width").addEventListener("input", storeFormFields);
  document.getElementById("step_size").addEventListener("input", storeFormFields);

  // Model changes
  document.getElementById("model_select").addEventListener("change", storeFormFields);

  // DS option toggles
  const dsRadios = document.getElementsByName("ds_option");
  dsRadios.forEach(radio => {
    radio.addEventListener("change", () => {
      storeFormFields();
      toggleOptionalFields();
      collectAndDisplayErrors();
    });
  });

  // Second chromosome toggle
  const togglePlusElemContainer = document.getElementById("toggle-second-chr-container");
  if (togglePlusElemContainer) {
    togglePlusElemContainer.addEventListener("click", toggleSecondChr);
  } else {
    console.error("Element with id 'toggle-second-chr-container' not found");
  }

  // CTCF selection & normalization
  const ctcfSelectElem = document.getElementById("ctcf_bw_path");
  ctcfSelectElem.addEventListener("change", () => {
    storeFormFields();
    updateCtcfNormalization();
    updateTrainingNormField();
  });

  // ATAC & CTCF norm radios
  document.getElementsByName("norm_atac").forEach(radio => {
    radio.addEventListener("change", () => {
      storeFormFields();
      updateTrainingNormField();
    });
  });
  document.getElementsByName("norm_ctcf").forEach(radio => {
    radio.addEventListener("change", () => {
      storeFormFields();
      updateTrainingNormField();
    });
  });

  // File uploads
  document.getElementById("atac_bw_file").addEventListener("change", () => {
    ajaxUploadFile("atac_bw_file", "atac", "atac_bw_path");
  });
  document.getElementById("ctcf_bw_file").addEventListener("change", () => {
    ajaxUploadFile("ctcf_bw_file", "ctcf", "ctcf_bw_path");
  });
  document.getElementById("peaks_file_upload").addEventListener("change", () => {
    ajaxUploadFile("peaks_file_upload", "peaks", "peaks_file_path");
  });

  // Instructions modal
  const modal = document.getElementById("instructions-modal");
  const btn   = document.getElementById("openModal");
  const span  = document.getElementsByClassName("close")[0];
  if (btn) {
    btn.onclick = (e) => {
      e.preventDefault();
      modal.style.display = "block";
    };
  }
  if (span) {
    span.onclick = () => {
      modal.style.display = "none";
    };
  }
  window.onclick = (event) => {
    if (event.target === modal) {
      modal.style.display = "none";
    }
  };

  /*************************************************************
   * Handle "Cancel" run logic
   *************************************************************/
  let runInProgress = false;
  const submitBtn = document.querySelector('input[type="submit"].submit-button');
  const formElem  = document.getElementById("corigami-form");

  // function cancelCurrentRun() {
  //   fetch("/cancel_run", {
  //     method: "POST",
  //     credentials: "same-origin"
  //   })
  //     .then(res => res.json())
  //     .then(data => {
  //       console.log("Cancel run response:", data);
  //       runInProgress = false;
  //       if (submitBtn) submitBtn.value = "Generate plots";
  //     })
  //     .catch(err => console.error("Error canceling run:", err));
  // }

  formElem.addEventListener("submit", function(e) {
    if (runInProgress) {
      e.preventDefault();
      return false;
    }

    e.preventDefault();
    if (!checkFormRequirements()) return false;

    const container = document.getElementById("output-container");
    container.innerHTML = '<div class="loader" style="display: block; margin: 0 auto;"></div>';

    // If second-chr toggle is off, clear the second region fields
    if (!document.getElementById('toggle-second-chr').classList.contains("active")) {
      document.getElementById("region_chr2").value   = "";
      document.getElementById("region_start2").value = "";
      document.getElementById("region_end2").value   = "";
    }

    const formData = new FormData(formElem);
    fetch("/", {
      method: "POST",
      body: formData,
      headers: { "X-Requested-With": "XMLHttpRequest" }
    })
    .then(response => {
      runInProgress = false;
      if (submitBtn) submitBtn.value = "Generate plots";
      if (!response.ok) throw new Error("Network response was not ok");
      return response.text();
    })
    .then(html => {
      container.innerHTML = html;
      executeScripts(container);
    })
    .catch(error => {
      console.error("Error during form submission:", error);
      container.innerHTML = "<p style='color:red;'>Error generating plots.</p>";
    });
  });

  // If server set screening_mode = true, auto-run screening
  if (typeof screening_mode !== "undefined" && screening_mode === true) {
    runScreening();
  }
});

/*************************************************************
 * checkScreeningRegionEdges + checkMinimumRegionSize (duplicate)
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
  console.log("checkMinimumRegionSize() is called");
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
