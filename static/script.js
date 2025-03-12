/*************************************************************
 * Constants and Global Variables
 *************************************************************/
const WINDOW_WIDTH = 2097152; // 2 Mb window
window.chromosomeLengths = {};
/*************************************************************
 * Local Storage Helpers
 *************************************************************/
function storeFormFields() {
  // Save first region and other fields
  const fields = [
    "region_chr", "region_start", "region_end",
    "del_start", "del_width", "perturb_width",
    "step_size", "model_select", "atac_bw_path",
    "ctcf_bw_path", "peaks_file_path"
  ];
  fields.forEach(field => {
    const elem = document.getElementById(field);
    if (elem) localStorage.setItem(field, elem.value);
  });

  // Save second region fields as well
  const secondFields = ["region_chr2", "region_start2", "region_end2"];
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
  // Fields for first region and other parameters
  const fields = [
    "region_chr", "region_start", "region_end",
    "del_start", "del_width", "perturb_width",
    "step_size", "model_select", "atac_bw_path",
    "ctcf_bw_path", "peaks_file_path", "genome_select"
  ];
  fields.forEach(field => {
    const storedVal = localStorage.getItem(field);
    const elem = document.getElementById(field);
    if (storedVal && elem) elem.value = storedVal;
  });

  // Fields for second region
  const secondFields = ["region_chr2", "region_start2", "region_end2"];
  secondFields.forEach(field => {
    const storedVal = localStorage.getItem(field);
    const elem = document.getElementById(field);
    if (storedVal && elem) elem.value = storedVal;
  });

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

  // Get both dropdown elements
  const chrSelect1 = document.getElementById("region_chr");
  const chrSelect2 = document.getElementById("region_chr2");
  const chrSelects = [chrSelect1, chrSelect2];

  fetch(url)
    .then(response => response.json())
    .then(data => {
      window.chromosomeLengths = data; // store for later validation
      chrSelects.forEach(select => {
        if (!select) return;
        select.innerHTML = ""; // Clear previous options

        // Iterate over each chromosome key in data
        for (const chr in data) {
          // Remove any leading "chr" (case-insensitive) for display only.
          let cleanChr = chr.replace(/^chr/i, "");
          let option = document.createElement("option");
          option.value = chr; // Preserve the original value (e.g., "chr1")
          option.text = `${cleanChr} (${data[chr]})`;
          select.appendChild(option);
        }
      });
      // Validate region bounds in case start/end values are already set
      validateRegionBounds();
    })
    .catch(error => {
      console.error("Error fetching chromosome lengths:", error);
    });
}


function validateRegionBounds() {
  // Validate region 1
  const errorElem1 = document.getElementById("region-bound-error");
  let errorMessage1 = "";
  const chrSelect = document.getElementById("region_chr");
  const regionStartElem = document.getElementById("region_start");
  const regionEndElem = document.getElementById("region_end");
  
  if (chrSelect && window.chromosomeLengths) {
    const selectedChr = chrSelect.value;
    const maxLength = window.chromosomeLengths[selectedChr];
    if (typeof maxLength === "undefined") {
      errorMessage1 = "Region must be within 0 and unknown bounds.";
    } else {
      const start = parseInt(regionStartElem.value);
      let end = (!regionEndElem.value || isNaN(parseInt(regionEndElem.value)))
                  ? start + WINDOW_WIDTH
                  : parseInt(regionEndElem.value);
      
      // Bound check first
      if (start < 0 || start > maxLength || end < 0 || end > maxLength) {
        errorMessage1 = `Region must be within 0 and ${maxLength} bounds.`;
      }
      // Then length check (only if no bound error)
      else if (end - start > WINDOW_WIDTH * 10) {
        errorMessage1 = "Region must be less than 20mb in length.";
      }
    }
  }
  
  if (errorMessage1) {
    errorElem1.textContent = errorMessage1;
    errorElem1.style.display = "block";
  } else {
    errorElem1.style.display = "none";
  }
  
  // Validate region 2 (if toggled on)
  const secondChrDiv = document.getElementById("second_chr_fields");
  const errorElem2 = document.getElementById("region-bound-error-2");
  let errorMessage2 = "";
  
  if (secondChrDiv && secondChrDiv.style.display !== "none") {
    const chrSelect2 = document.getElementById("region_chr2");
    const regionStartElem2 = document.getElementById("region_start2");
    const regionEndElem2 = document.getElementById("region_end2");
    
    if (chrSelect2 && window.chromosomeLengths) {
      const selectedChr2 = chrSelect2.value;
      const maxLength2 = window.chromosomeLengths[selectedChr2];
      if (typeof maxLength2 === "undefined") {
        errorMessage2 = "Region must be within 0 and unknown bounds.";
      } else {
        const start2 = parseInt(regionStartElem2.value);
        let end2 = (!regionEndElem2.value || isNaN(parseInt(regionEndElem2.value)))
                   ? start2 + WINDOW_WIDTH
                   : parseInt(regionEndElem2.value);
        
        // Bound check first
        if (start2 < 0 || start2 > maxLength2 || end2 < 0 || end2 > maxLength2) {
          errorMessage2 = `Region must be within 0 and ${maxLength2} bounds.`;
        }
        // Then length check (only if no bound error)
        else if (end2 - start2 > WINDOW_WIDTH * 10) {
          errorMessage2 = "Region must be less than 20mb in length.";
        }
      }
    }
  }
  
  if (errorElem2) {
    if (errorMessage2) {
      errorElem2.textContent = errorMessage2;
      errorElem2.style.display = "block";
    } else {
      errorElem2.style.display = "none";
    }
  }
}





function updateEndPosition() {
  // First region
  const startField1 = document.getElementById("region_start");
  const endField1 = document.getElementById("region_end");
  const startVal1 = parseInt(startField1.value);
  
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
  
  // Second region (only if toggled on)
  const secondRegionDiv = document.getElementById("second_chr_fields");
  if (secondRegionDiv && secondRegionDiv.style.display !== "none") {
    const startField2 = document.getElementById("region_start2");
    const endField2 = document.getElementById("region_end2");
    if (startField2 && endField2) {
      let startVal2 = parseInt(startField2.value);
      // If no value is entered, assign a default value (adjust as needed)
      if (isNaN(startVal2)) {
        startVal2 = 1500000;  // default for second region
        startField2.value = startVal2; // optionally set the field's value
      }
      const computedEnd2 = startVal2 + WINDOW_WIDTH;
      if (!endField2.value || isNaN(parseInt(endField2.value))) {
        endField2.placeholder = computedEnd2;
      } else {
        endField2.placeholder = "";
      }
    }
  }
}


function toggleOptionalFields() {
  const dsOption = document.querySelector('input[name="ds_option"]:checked');
  const delFields = document.getElementById("deletion-fields");
  const scrFields = document.getElementById("screening-fields-2");
  const peaksContainer = document.getElementById("peaks_file_container");
  // Get the plus toggle element for adding a second chromosome
  const togglePlusElem = document.getElementById("toggle-second-chr");

  if (dsOption) {
    if (dsOption.value === "deletion") {
      if (delFields) delFields.style.display = "block";
      if (scrFields) scrFields.style.display = "none";
      if (peaksContainer) peaksContainer.style.display = "none";
      // Hide the plus toggle when in deletion mode
      if (togglePlusElem) togglePlusElem.style.display = "none";
    } else if (dsOption.value === "screening") {
      if (scrFields) scrFields.style.display = "block";
      if (delFields) delFields.style.display = "none";
      if (peaksContainer) peaksContainer.style.display = "block";
    } else {
      if (delFields) delFields.style.display = "none";
      if (scrFields) scrFields.style.display = "none";
      if (peaksContainer) peaksContainer.style.display = "none";
      if (togglePlusElem) togglePlusElem.style.display = "block";
    }
  }
}


function checkFormRequirements() {
  const dsOption = document.querySelector('input[name="ds_option"]:checked');
  if (!dsOption) return true;
  if (dsOption.value === "deletion") {
    if (!document.getElementById("del_start").value.trim() || !document.getElementById("del_width").value.trim()) {
      alert("Please provide both Deletion Start and Deletion Width for Deletion mode.");
      return false;
    }
  } else if (dsOption.value === "screening") {
    if (!document.getElementById("perturb_width").value.trim() || !document.getElementById("step_size").value.trim()) {
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
  const xhr = new XMLHttpRequest();
  xhr.open("GET", "/list_uploads?file_type=" + encodeURIComponent(fileType));
  xhr.onload = function() {
    if (xhr.status === 200) {
      const serverFiles = JSON.parse(xhr.responseText);
      for (let i = selectElem.options.length - 1; i >= 0; i--) {
        const option = selectElem.options[i];
        if (option.value !== "none" && !option.value.startsWith("./corigami_data")) {
          selectElem.remove(i);
        }
      }
      serverFiles.forEach(file => {
        const exists = Array.from(selectElem.options).some(opt => opt.value === file.value);
        if (!exists) {
          const newOption = document.createElement("option");
          newOption.value = file.value;
          newOption.text = file.name;
          selectElem.appendChild(newOption);
        }
      });
      const storedVal = localStorage.getItem(selectId);
      if (storedVal) {
        for (let option of selectElem.options) {
          if (option.value === storedVal) {
            selectElem.value = storedVal;
            break;
          }
        }
      }
      console.log("Dropdown", selectId, "populated with:", serverFiles);
    }
  };
  xhr.send();
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

  // 1) Read regionChr from the front-end (or from window.screening_params if you prefer).
  let regionChr = document.getElementById("region_chr").value;

  // 2) Parse region_start
  let regionStartVal = parseInt(document.getElementById("region_start").value, 10);
  if (isNaN(regionStartVal)) {
    // Fallback: if user left start blank or typed something invalid, 
    // let's set 1,048,576 or 0 or some default:
    regionStartVal = 1048576; 
  }

  // 3) Parse region_end
  let regionEndVal = parseInt(document.getElementById("region_end").value, 10);
  if (isNaN(regionEndVal)) {
    // If the user left end blank or typed non-numeric => fallback
    regionEndVal = regionStartVal + 2097152; // default to 2 Mb beyond start
    // Optionally update the form field so the user sees the new value
    document.getElementById("region_end").value = regionEndVal;
  }

  // 4) Example check: If standard mode (chr != chrCHIM) & regionStart < 1048576 => warn
  if (regionChr !== "chrCHIM" && regionStartVal < 1048576) {
    alert("Cannot run screening if start < 1,048,576 in standard mode!");
    return;
  }

  // (If you want to handle chimeric lengths automatically, e.g. if regionChr === "chrCHIM" 
  //  you might override regionEndVal with the actual chim_len from window.screening_params
  //  but that's optional.)

  // 5) Show a loader in the screening chart container
  const screeningContainerElem = document.getElementById("screening_chart");
  if (screeningContainerElem) {
    screeningContainerElem.innerHTML = "<div class='loader' style='display:block;margin:0 auto;'></div>";
  }

  // 6) Build params
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

  // 7) AJAX to /run_screening
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
  console.log("Sending screening request to /run_screening?" + queryParams);
  xhr.send();
}


function updateCtcfNormalization() {
  const ctcfSelect = document.getElementById('ctcf_bw_path');
  const ctcfSwitch = document.getElementById('ctcf-switch');
  const ctcfNoNorm = document.getElementById('ctcf-no-norm');
  const ctcfLog = document.getElementById('ctcf-log');
  const ctcfMinmax = document.getElementById('ctcf-minmax');
  if (ctcfSelect.value === 'none') {
    ctcfNoNorm.disabled = true;
    ctcfLog.disabled = true;
    ctcfMinmax.disabled = true;
    ctcfSwitch.classList.add('disabled-toggle');
  } else {
    ctcfNoNorm.disabled = false;
    ctcfLog.disabled = false;
    ctcfMinmax.disabled = false;
    ctcfSwitch.classList.remove('disabled-toggle');
  }
}

function updateTrainingNormField() {
  const ctcfFile = document.getElementById("ctcf_bw_path").value;
  const atacState = document.querySelector('input[name="norm_atac"]:checked').value;
  const ctcfState = document.querySelector('input[name="norm_ctcf"]:checked').value;
  const trainingContainer = document.getElementById("training-norm-container");
  if (ctcfFile === "none") {
    document.getElementById("training-minmax").checked = true;
    document.getElementById("training-log").disabled = true;
    trainingContainer.classList.add("disabled-toggle");
  } else if (atacState !== "none") {
    if (atacState === "log") {
      document.getElementById("training-log").checked = true;
    } else if (atacState === "minmax") {
      document.getElementById("training-minmax").checked = true;
    }
    document.getElementById("training-log").disabled = true;
    document.getElementById("training-minmax").disabled = true;
    trainingContainer.classList.add("disabled-toggle");
  } else if (ctcfState !== "none") {
    if (ctcfState === "log") {
      document.getElementById("training-log").checked = true;
    } else if (ctcfState === "minmax") {
      document.getElementById("training-minmax").checked = true;
    }
    document.getElementById("training-log").disabled = true;
    document.getElementById("training-minmax").disabled = true;
    trainingContainer.classList.add("disabled-toggle");
  } else {
    document.getElementById("training-log").disabled = false;
    document.getElementById("training-minmax").disabled = false;
    trainingContainer.classList.remove('disabled-toggle');
  }
}

/*************************************************************
 * Validate Deletion Area (Single Definition)
 *************************************************************/
function validateDeletionArea() {
  // Check if we're in deletion mode
  const dsOption = document.querySelector('input[name="ds_option"]:checked');
  if (!dsOption || dsOption.value !== "deletion") {
    // Optionally hide the error message if not in deletion mode
    const errorElem = document.getElementById("deletion-error");
    if (errorElem) errorElem.style.display = "none";
    return;
  }

  const regionStartElem = document.getElementById("region_start");
  const regionEndElem = document.getElementById("region_end");
  const delStartElem = document.getElementById("del_start");
  const delWidthElem = document.getElementById("del_width");
  const errorElem = document.getElementById("deletion-error");
  
  if (!regionStartElem || !delStartElem || !delWidthElem) return;
  
  const regionStart = parseInt(regionStartElem.value);
  const regionEnd = (regionEndElem && regionEndElem.value && !isNaN(parseInt(regionEndElem.value)))
                      ? parseInt(regionEndElem.value)
                      : regionStart + WINDOW_WIDTH;
  const deletionStart = parseInt(delStartElem.value);
  const deletionWidth = parseInt(delWidthElem.value);
  const deletionEnd = deletionStart + deletionWidth;
  
  if (deletionStart < regionStart || deletionEnd > regionEnd) {
    errorElem.style.display = "block";
    errorElem.textContent = "Deletion area is out of bounds";
  } else {
    errorElem.style.display = "none";
  }
}


/*************************************************************
 * Window onload: Attach Event Listeners and Initialize
 *************************************************************/
function toggleSecondChr() {
  const dsOption = document.querySelector('input[name="ds_option"]:checked');
  console.log("toggleSecondChr called with dsOption:", dsOption);
  // Hide second chromosome fields if in deletion or screening mode
  if (dsOption && (dsOption.value === "deletion")) {
    const secondChrDiv = document.getElementById('second_chr_fields');
    const togglePlus = document.getElementById('toggle-second-chr');
    const chimericInput = document.getElementById('chimeric_active');
    secondChrDiv.style.display = "none";
    togglePlus.classList.remove("active");
    chimericInput.value = "false";
    return;  // exit early
  }
  
  // Otherwise, proceed with the normal toggle behavior
  const secondChrDiv = document.getElementById('second_chr_fields');
  const togglePlus = document.getElementById('toggle-second-chr');
  const chimericInput = document.getElementById('chimeric_active');
  
  if (secondChrDiv.style.display === "none" || secondChrDiv.style.display === "") {
    // Show second chromosome fields and set chimeric flag to true
    secondChrDiv.style.display = "block";
    togglePlus.classList.add("active");
    chimericInput.value = "true";
    
    // Optionally set a default value for second region start if empty
    const regionStart2 = document.getElementById("region_start2");
    if (regionStart2 && !regionStart2.value) {
      regionStart2.value = "1500000";
    }
    updateEndPosition();
  } else {
    // Hide second chromosome fields and set chimeric flag to false
    secondChrDiv.style.display = "none";
    togglePlus.classList.remove("active");
    chimericInput.value = "false";
  }
}




document.addEventListener('DOMContentLoaded', function() {
  console.log("DOM fully loaded");
  
  // Your initialization code:
  populateChromosomeDropdown();
  restoreFormFields();
  updateEndPosition();
  toggleOptionalFields();
  updateCtcfNormalization();
  updateTrainingNormField();
  validateDeletionArea();
  validateRegionBounds();

  document.getElementById("region_start").addEventListener("input", () => {
    storeFormFields();
    updateEndPosition();
    validateDeletionArea();
    validateRegionBounds();
  });
  document.getElementById("region_end").addEventListener("input", () => {
    validateDeletionArea();
    validateRegionBounds();
    updateEndPosition();
  });
  document.getElementById("region_end2").addEventListener("input", () => {
    validateRegionBounds();
    updateEndPosition();
  });
  document.getElementById("region_chr").addEventListener("change", () => {
    storeFormFields();
    validateRegionBounds();
  });
  document.getElementById("del_start").addEventListener("input", () => {
    storeFormFields();
    validateDeletionArea();
  });
  document.getElementById("del_width").addEventListener("input", () => {
    storeFormFields();
    validateDeletionArea();
  });
  document.getElementById("perturb_width").addEventListener("input", storeFormFields);
  document.getElementById("step_size").addEventListener("input", storeFormFields);
  document.getElementById("region_chr").addEventListener("change", storeFormFields);
  document.getElementById("model_select").addEventListener("change", storeFormFields);

  const dsRadios = document.getElementsByName("ds_option");
  dsRadios.forEach(radio => {
    radio.addEventListener("change", () => {
      storeFormFields();
      toggleOptionalFields();
      validateDeletionArea();
    });
  });

  // Attach event listener for the toggle-plus element
  const togglePlusElem = document.getElementById("toggle-second-chr");
  if (togglePlusElem) {
    togglePlusElem.addEventListener("click", toggleSecondChr);
    console.log("toggle-plus event listener attached");
  } else {
    console.error("Element with id 'toggle-second-chr' not found");
  }

  const startField2 = document.getElementById("region_start2");
  if (startField2) {
    startField2.addEventListener("input", () => {
      storeFormFields();
      updateEndPosition();
      validateRegionBounds();
    });
  }
  const ctcfSelectElem = document.getElementById("ctcf_bw_path");
  ctcfSelectElem.addEventListener("change", () => {
    storeFormFields();
    updateCtcfNormalization();
    updateTrainingNormField();
  });

  document.getElementsByName("norm_atac").forEach(radio => {
    radio.addEventListener("change", updateTrainingNormField);
  });
  document.getElementsByName("norm_ctcf").forEach(radio => {
    radio.addEventListener("change", updateTrainingNormField);
  });

  // Populate dropdowns
  populateDropdownFromServer("atac_bw_path", "atac");
  populateDropdownFromServer("ctcf_bw_path", "ctcf");
  populateDropdownFromServer("peaks_file_path", "peaks");

  // Attach immediate AJAX upload for file inputs
  document.getElementById("atac_bw_file").addEventListener("change", () => {
    ajaxUploadFile("atac_bw_file", "atac", "atac_bw_path");
  });
  document.getElementById("ctcf_bw_file").addEventListener("change", () => {
    ajaxUploadFile("ctcf_bw_file", "ctcf", "ctcf_bw_path");
  });
  document.getElementById("peaks_file_upload").addEventListener("change", () => {
    ajaxUploadFile("peaks_file_upload", "peaks", "peaks_file_path");
  });

  // Modal handling for instructions
  const modal = document.getElementById("instructions-modal");
  const btn = document.getElementById("openModal");
  const span = document.getElementsByClassName("close")[0];
  if (btn) {
    btn.onclick = (e) => {
      e.preventDefault();
      modal.style.display = "block";
    };
  }
  if (span) {
    span.onclick = () => { modal.style.display = "none"; };
  }
  window.onclick = (event) => {
    if (event.target === modal) modal.style.display = "none";
  };

  // Form submission via AJAX
  const formElem = document.getElementById("corigami-form");
  formElem.addEventListener("submit", function(e) {
    const togglePlusElem = document.getElementById("toggle-second-chr");
    if (!togglePlusElem.classList.contains("active")) {
      // Clear second region fields if the plus toggle is off.
      document.getElementById("region_chr2").value = "";
      document.getElementById("region_start2").value = "";
      document.getElementById("region_end2").value = "";
    }
    e.preventDefault();
    if (!checkFormRequirements()) return false;
    const container = document.getElementById("output-container");
    container.innerHTML = '<div class="loader" style="display: block; margin: 0 auto;"></div>';
    const formData = new FormData(formElem);
    fetch("/", {
      method: "POST",
      body: formData,
      headers: { "X-Requested-With": "XMLHttpRequest" }
    })
    .then(response => {
      if (!response.ok) throw new Error("Network response was not ok");
      return response.text();
    })
    .then(html => {
      container.innerHTML = html;
      executeScripts(container);
    })
    .catch(error => console.error("Error during form submission:", error));
  });

  if (typeof screening_mode !== "undefined" && screening_mode === true) {
    runScreening();
  }
});
