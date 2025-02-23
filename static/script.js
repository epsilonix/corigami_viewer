/*************************************************************
 * Constants, Global Vars
 *************************************************************/
const WINDOW_WIDTH = 2097152;

/*************************************************************
 * Local Storage Helpers
 *************************************************************/

/**
 * Saves relevant form field values to localStorage so they persist.
 */
function storeFormFields() {
  // Save numeric or text inputs
  localStorage.setItem("region_chr", document.getElementById("region_chr").value);
  localStorage.setItem("region_start", document.getElementById("region_start").value);
  localStorage.setItem("del_start", document.getElementById("del_start").value);
  localStorage.setItem("del_width", document.getElementById("del_width").value);
  localStorage.setItem("perturb_width", document.getElementById("perturb_width").value);
  localStorage.setItem("step_size", document.getElementById("step_size").value);
  localStorage.setItem("model_select", document.getElementById("model_select").value);

  // Save ds_option (radio)
  var dsRadios = document.getElementsByName("ds_option");
  for (var i = 0; i < dsRadios.length; i++) {
    if (dsRadios[i].checked) {
      localStorage.setItem("ds_option", dsRadios[i].value);
      break;
    }
  }

  // Save currently selected file from dropdowns
  localStorage.setItem("atac_bw_path", document.getElementById("atac_bw_path").value);
  localStorage.setItem("ctcf_bw_path", document.getElementById("ctcf_bw_path").value);
  localStorage.setItem("peaks_file_path", document.getElementById("peaks_file_path").value);
}

/**
 * Restores saved form field values from localStorage, if available.
 */
function restoreFormFields() {
  // Basic numeric or text inputs
  var regionChr = localStorage.getItem("region_chr");
  if (regionChr) document.getElementById("region_chr").value = regionChr;

  var regionStart = localStorage.getItem("region_start");
  if (regionStart) document.getElementById("region_start").value = regionStart;

  var delStart = localStorage.getItem("del_start");
  if (delStart) document.getElementById("del_start").value = delStart;

  var delWidth = localStorage.getItem("del_width");
  if (delWidth) document.getElementById("del_width").value = delWidth;

  var perturbWidth = localStorage.getItem("perturb_width");
  if (perturbWidth) document.getElementById("perturb_width").value = perturbWidth;

  var stepSize = localStorage.getItem("step_size");
  if (stepSize) document.getElementById("step_size").value = stepSize;

  var modelSelect = localStorage.getItem("model_select");
  if (modelSelect) document.getElementById("model_select").value = modelSelect;

  // ds_option (radio)
  var dsOption = localStorage.getItem("ds_option");
  if (dsOption) {
    var radioElem = document.getElementById("ds_" + dsOption);
    if (radioElem) {
      radioElem.checked = true;
    }
  }

  // File dropdowns – store only the selected value here
  var storedAtac = localStorage.getItem("atac_bw_path");
  if (storedAtac) document.getElementById("atac_bw_path").value = storedAtac;

  var storedCtcf = localStorage.getItem("ctcf_bw_path");
  if (storedCtcf) document.getElementById("ctcf_bw_path").value = storedCtcf;

  var storedPeaks = localStorage.getItem("peaks_file_path");
  if (storedPeaks) document.getElementById("peaks_file_path").value = storedPeaks;
}

/*************************************************************
 * Form Behavior Helpers
 *************************************************************/

/**
 * Update the region_end field based on region_start + WINDOW_WIDTH.
 */
function updateEndPosition() {
  var startField = document.getElementById("region_start");
  var endField = document.getElementById("region_end");
  var startVal = parseInt(startField.value);
  if (!isNaN(startVal)) {
    endField.value = startVal + WINDOW_WIDTH;
  } else {
    endField.value = "";
  }
}

/**
 * Show/hide the deletion/screening parameters based on the ds_option radio.
 */
function toggleOptionalFields() {
  var dsOption = document.querySelector('input[name="ds_option"]:checked');
  var delFields = document.getElementById("deletion-fields");
  var scrFields = document.getElementById("screening-fields-2");
  var peaksContainer = document.getElementById("peaks_file_container");

  if (dsOption) {
    if (dsOption.value === "deletion") {
      if (delFields) delFields.style.display = "block";
      if (scrFields) scrFields.style.display = "none";
      if (peaksContainer) peaksContainer.style.display = "none";
    } else if (dsOption.value === "screening") {
      if (scrFields) scrFields.style.display = "block";
      if (delFields) delFields.style.display = "none";
      if (peaksContainer) peaksContainer.style.display = "block";
    } else {
      if (delFields) delFields.style.display = "none";
      if (scrFields) scrFields.style.display = "none";
      if (peaksContainer) peaksContainer.style.display = "none";
    }
  }
}

/**
 * Validate that the deletion region is within the 2 Mb window.
 */
function validateDeletionArea() {
  var regionStart = parseInt(document.getElementById("region_start").value);
  var regionEnd = regionStart + WINDOW_WIDTH;
  var deletionStart = parseInt(document.getElementById("del_start").value);
  var deletionWidth = parseInt(document.getElementById("del_width").value);
  var deletionEnd = deletionStart + deletionWidth;
  var errorElem = document.getElementById("deletion-error");

  if (deletionStart < regionStart || deletionEnd > regionEnd) {
    errorElem.style.display = "block";
  } else {
    errorElem.style.display = "none";
  }
}

/**
 * Check if required fields (del_start/del_width for deletion,
 * or perturb_width/step_size for screening) are filled before submitting.
 */
function checkFormRequirements() {
  var dsOption = document.querySelector('input[name="ds_option"]:checked');
  if (!dsOption) return true; // If nothing selected, or 'none', no special fields

  var dsVal = dsOption.value;
  if (dsVal === "deletion") {
    var delStartVal = document.getElementById("del_start").value.trim();
    var delWidthVal = document.getElementById("del_width").value.trim();
    if (!delStartVal || !delWidthVal) {
      alert("Please provide both Deletion Start and Deletion Width for Deletion mode.");
      return false;
    }
  }
  else if (dsVal === "screening") {
    var perturbVal = document.getElementById("perturb_width").value.trim();
    var stepVal = document.getElementById("step_size").value.trim();
    if (!perturbVal || !stepVal) {
      alert("Please provide both Perturb Width and Step Size for Screening mode.");
      return false;
    }
  }

  // If dsVal === "none", no extra fields required
  return true;
}

/*************************************************************
 * File Dropdown Handling
 *************************************************************/

/**
 * Call the server to list currently uploaded files for this user session,
 * then populate the <select> with those filenames.
 * Finally, attach a listener to the <input type="file"> for new uploads.
 */
function populateDropdownFromServer(selectId, fileInputId, storageKey) {
  var selectElem = document.getElementById(selectId);

  // 1. Fetch the user’s existing upload list
  var xhr = new XMLHttpRequest();
  xhr.open("GET", "/list_uploads");
  xhr.onload = function() {
    if (xhr.status === 200) {
      var serverFiles = JSON.parse(xhr.responseText); // array of filenames
      // Remove previously added dynamic options (optional)
      for (var i = selectElem.options.length - 1; i >= 0; i--) {
        var val = selectElem.options[i].value;
        // Keep "none" or "preset" or known static options, remove dynamic
        if (
          val !== "none" &&
          !val.startsWith("./corigami_data") &&
          val !== "Preset (hg38/imr90)" &&
          val !== "preset"
        ) {
          selectElem.remove(i);
        }
      }
      // Add the server files to the dropdown
      serverFiles.forEach(function(filename) {
        // Check if it already exists
        var exists = false;
        for (var j = 0; j < selectElem.options.length; j++) {
          if (selectElem.options[j].value === filename) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          var newOption = document.createElement("option");
          newOption.value = filename;
          newOption.text = filename;
          selectElem.appendChild(newOption);
        }
      });

      // Now restore the user’s previously selected file if it still exists
      var storedVal = localStorage.getItem(selectId);
      if (storedVal) {
        // Only re-select if it exists in the current <select> after refresh
        for (var k = 0; k < selectElem.options.length; k++) {
          if (selectElem.options[k].value === storedVal) {
            selectElem.value = storedVal;
            break;
          }
        }
      }
    }
  };
  xhr.send();

  // 2. Attach a listener to the file input so that any new upload
  //    gets stored in localStorage and appended to the dropdown
  updateDropdown(fileInputId, selectId, storageKey);
}

/**
 * Original function to handle user picking a new file in the <input type="file">.
 * This also updates localStorage so next time we see the file name in the dropdown.
 */
function updateDropdown(fileInputId, selectId, storageKey) {
  var fileInput = document.getElementById(fileInputId);
  var selectElem = document.getElementById(selectId);

  fileInput.addEventListener("change", function() {
    if (fileInput.files.length > 0) {
      var filename = fileInput.files[0].name;
      // Add it to the dropdown (client side) - the server will rename it,
      // but for user convenience we show the local name here
      var newOption = document.createElement("option");
      newOption.value = filename;
      newOption.text = filename;

      // Check if it already exists
      var exists = false;
      for (var i = 0; i < selectElem.options.length; i++) {
        if (selectElem.options[i].value === filename) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        selectElem.appendChild(newOption);
      }
      selectElem.value = filename;

      // Save to localStorage for future page loads
      localStorage.setItem(selectId, filename);
    }
  });
}

/*************************************************************
 * Screening
 *************************************************************/

function runScreening() {
  if (parseInt(document.getElementById("region_start").value) < 1048576) {
    alert("Cannot run screening on start position below 1,048,576!");
    return;
  }
  document.getElementById("screening_loader").style.display = "block";
  document.getElementById("screening_image_container").style.display = "none";

  // Use the global screening_params passed from the server (if available)
  var paramsObj = screening_params || {};

  // Fallback: build parameters manually if screening_params not defined
  if (!paramsObj.region_chr) {
    paramsObj = {
      region_chr: document.getElementById("region_chr").value,
      model_select: document.getElementById("model_select").value,
      region_start: document.getElementById("region_start").value,
      region_end: document.getElementById("region_end").value,
      perturb_width: document.getElementById("perturb_width").value,
      step_size: document.getElementById("step_size").value,
      ctcf_bw_path: document.getElementById("ctcf_bw_path").value,
      atac_bw_path: document.getElementById("atac_bw_path").value,
      // The server code in app.py should override or handle if output_dir is session-specific
      output_dir: "./output_new",
      peaks_file: document.getElementById("peaks_file_path").value
    };
  }
  
  // Make an AJAX GET request to /run_screening with these params
  var xhr = new XMLHttpRequest();
  var queryParams = new URLSearchParams(paramsObj).toString();
  
  xhr.onreadystatechange = function() {
    if (xhr.readyState === 4) {
      document.getElementById("screening_loader").style.display = "none";
      if (xhr.status === 200) {
        var response = JSON.parse(xhr.responseText);
        if (response.message) {
          document.getElementById("screening_image_container").innerHTML = "<p>" + response.message + "</p>";
        } else {
          document.getElementById("screening_image_container").innerHTML =
            '<img src="' + response.screening_image + '" alt="Screening Plot" style="max-width:100%;">';
        }
      } else {
        document.getElementById("screening_image_container").innerHTML = "<p>Error generating screening plot.</p>";
      }
      document.getElementById("screening_image_container").style.display = "block";
    }
  };
  xhr.open("GET", "/run_screening?" + queryParams, true);
  xhr.send();
}

/*************************************************************
 * CTCF Normalization Toggle
 *************************************************************/
function updateCtcfNormalization() {
  var ctcfSelect = document.getElementById('ctcf_bw_path');
  var ctcfSwitch = document.getElementById('ctcf-switch');
  var ctcfNoNorm = document.getElementById('ctcf-no-norm');
  var ctcfLog = document.getElementById('ctcf-log');
  var ctcfMinmax = document.getElementById('ctcf-minmax');

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

/*************************************************************
 * Window onload
 *************************************************************/
window.onload = function() {
  // 1. Restore fields from localStorage
  restoreFormFields();

  // 2. Update region_end once region_start is known
  updateEndPosition();

  // 3. Toggle optional fields (deletion/screening) accordingly
  toggleOptionalFields();

  // 4. CTCF normalization toggle
  updateCtcfNormalization();

  // 5. Attach listeners to store new values in localStorage
  document.getElementById("region_start").addEventListener("input", function() {
    storeFormFields();
    updateEndPosition();
  });
  document.getElementById("del_start").addEventListener("input", function() {
    storeFormFields();
    validateDeletionArea();
  });
  document.getElementById("del_width").addEventListener("input", function() {
    storeFormFields();
    validateDeletionArea();
  });
  document.getElementById("perturb_width").addEventListener("input", storeFormFields);
  document.getElementById("step_size").addEventListener("input", storeFormFields);
  document.getElementById("region_chr").addEventListener("change", storeFormFields);
  document.getElementById("model_select").addEventListener("change", storeFormFields);

  var dsRadios = document.getElementsByName("ds_option");
  for (var i = 0; i < dsRadios.length; i++) {
    dsRadios[i].addEventListener("change", function() {
      storeFormFields();
      toggleOptionalFields();
    });
  }

  // 6. When CTCF file selection changes
  var ctcfSelectElem = document.getElementById("ctcf_bw_path");
  ctcfSelectElem.addEventListener("change", function() {
    storeFormFields();
    updateCtcfNormalization();
  });

  // 7. Populate dropdowns from server (actual existing uploaded files)
  populateDropdownFromServer("atac_bw_path", "atac_bw_file", "atac_bw_options");
  populateDropdownFromServer("ctcf_bw_path", "ctcf_bw_file", "ctcf_bw_options");
  populateDropdownFromServer("peaks_file_path", "peaks_file_upload", "peaks_file_options");

  // 8. Modal functionality
  var modal = document.getElementById("instructions-modal");
  var btn = document.getElementById("openModal");
  var span = document.getElementsByClassName("close")[0];
  if (btn) {
    btn.onclick = function(e) {
      e.preventDefault();
      modal.style.display = "block";
    };
  }
  if (span) {
    span.onclick = function() {
      modal.style.display = "none";
    };
  }
  window.onclick = function(event) {
    if (event.target == modal) {
      modal.style.display = "none";
    }
  };

  // 9. Attach a submit event listener on the form to validate required fields dynamically
  var formElem = document.getElementById("corigami-form");
  formElem.addEventListener("submit", function(e) {
    if (!checkFormRequirements()) {
      e.preventDefault(); // Stop submission if validation fails
      return false;
    }
  });

  // 10. If screening_mode was set by the server, automatically run screening
  if (typeof screening_mode !== "undefined" && screening_mode === true) {
    runScreening();
  }
};
