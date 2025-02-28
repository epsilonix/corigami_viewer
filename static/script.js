/*************************************************************
 * Constants, Global Vars
 *************************************************************/
const WINDOW_WIDTH = 2097152;

/*************************************************************
 * Local Storage Helpers
 *************************************************************/
function storeFormFields() {
  localStorage.setItem("region_chr", document.getElementById("region_chr").value);
  localStorage.setItem("region_start", document.getElementById("region_start").value);
  localStorage.setItem("del_start", document.getElementById("del_start").value);
  localStorage.setItem("del_width", document.getElementById("del_width").value);
  localStorage.setItem("perturb_width", document.getElementById("perturb_width").value);
  localStorage.setItem("step_size", document.getElementById("step_size").value);
  localStorage.setItem("model_select", document.getElementById("model_select").value);
  var dsRadios = document.getElementsByName("ds_option");
  for (var i = 0; i < dsRadios.length; i++) {
    if (dsRadios[i].checked) {
      localStorage.setItem("ds_option", dsRadios[i].value);
      break;
    }
  }
  localStorage.setItem("atac_bw_path", document.getElementById("atac_bw_path").value);
  localStorage.setItem("ctcf_bw_path", document.getElementById("ctcf_bw_path").value);
  localStorage.setItem("peaks_file_path", document.getElementById("peaks_file_path").value);
}

function restoreFormFields() {
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
  var dsOption = localStorage.getItem("ds_option");
  if (dsOption) {
    var radioElem = document.getElementById("ds_" + dsOption);
    if (radioElem) radioElem.checked = true;
  }
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

function checkFormRequirements() {
  var dsOption = document.querySelector('input[name="ds_option"]:checked');
  if (!dsOption) return true;
  var dsVal = dsOption.value;
  if (dsVal === "deletion") {
    var delStartVal = document.getElementById("del_start").value.trim();
    var delWidthVal = document.getElementById("del_width").value.trim();
    if (!delStartVal || !delWidthVal) {
      alert("Please provide both Deletion Start and Deletion Width for Deletion mode.");
      return false;
    }
  } else if (dsVal === "screening") {
    var perturbVal = document.getElementById("perturb_width").value.trim();
    var stepVal = document.getElementById("step_size").value.trim();
    if (!perturbVal || !stepVal) {
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
  var selectElem = document.getElementById(selectId);
  var xhr = new XMLHttpRequest();
  xhr.open("GET", "/list_uploads?file_type=" + encodeURIComponent(fileType));
  xhr.onload = function() {
    if (xhr.status === 200) {
      var serverFiles = JSON.parse(xhr.responseText);
      for (var i = selectElem.options.length - 1; i >= 0; i--) {
        var option = selectElem.options[i];
        if (option.value !== "none" && !option.value.startsWith("./corigami_data")) {
          selectElem.remove(i);
        }
      }
      serverFiles.forEach(function(file) {
        var exists = false;
        for (var j = 0; j < selectElem.options.length; j++) {
          if (selectElem.options[j].value === file.value) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          var newOption = document.createElement("option");
          newOption.value = file.value;
          newOption.text = file.name;
          selectElem.appendChild(newOption);
        }
      });
      var storedVal = localStorage.getItem(selectId);
      if (storedVal) {
        for (var k = 0; k < selectElem.options.length; k++) {
          if (selectElem.options[k].value === storedVal) {
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
  var fileInput = document.getElementById(fileInputId);
  if (!fileInput.files || fileInput.files.length === 0) {
    console.log("No file selected for", fileType);
    return;
  }
  var file = fileInput.files[0];
  console.log("Selected file for", fileType + ":", file.name);
  
  var formData = new FormData();
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
      var dropdown = document.getElementById(dropdownId);
      var newOption = document.createElement("option");
      newOption.value = data.saved_path;
      newOption.text = data.display_name;
      dropdown.appendChild(newOption);
      dropdown.value = data.saved_path;
      localStorage.setItem(dropdownId, data.saved_path);
    }
  })
  .catch(error => console.error("Error uploading file for " + fileType + ":", error));
}

/*************************************************************
 * Helper to Execute Inline Scripts After AJAX Update
 *************************************************************/
function executeScripts(container) {
  console.log("Executing scripts from container:", container);
  const scripts = container.querySelectorAll("script");
  scripts.forEach(script => {
    console.log("Found script:", script);
    const newScript = document.createElement("script");
    if (script.src) {
      newScript.src = script.src;
      newScript.onload = function() {
        console.log("Loaded external script:", script.src);
      };
    } else {
      newScript.text = script.textContent;
      console.log("Executing inline script (first 100 chars):", script.textContent.substring(0, 100));
    }
    document.head.appendChild(newScript);
  });
}

/*************************************************************
 * Screening, CTCF Normalization & Training Norm Field
 *************************************************************/
function runScreening() {
  console.log("runScreening() is called");

  const regionStartElem = document.getElementById("region_start");
  if (parseInt(regionStartElem.value) < 1048576) {
    alert("Cannot run screening on start position below 1,048,576!");
    return;
  }

  // Use the screening_chart element as the container.
  const screeningContainerElem = document.getElementById("screening_chart");
  if (screeningContainerElem) {
    screeningContainerElem.innerHTML = "<p>Loading screening plot...</p>";
  }

  var paramsObj = {
    region_chr: document.getElementById("region_chr").value,
    model_select: document.getElementById("model_select").value,
    region_start: document.getElementById("region_start").value,
    region_end: document.getElementById("region_end").value,
    perturb_width: document.getElementById("perturb_width").value,
    step_size: document.getElementById("step_size").value,
    ctcf_bw_path: document.getElementById("ctcf_bw_path").value,
    atac_bw_path: document.getElementById("atac_bw_path").value,
    output_dir: document.getElementById("output_dir").value,
    peaks_file: document.getElementById("peaks_file_path").value
  };

  console.log("Sending screening request with parameters:", paramsObj);
  const queryParams = new URLSearchParams(paramsObj).toString();
  console.log("Query string for screening request:", queryParams);

  const xhr = new XMLHttpRequest();
  xhr.onreadystatechange = function() {
    console.log("XHR readyState:", xhr.readyState, "status:", xhr.status);
    if (xhr.readyState === 4) {
      if (xhr.status === 200) {
        console.log("Screening request successful. Response:", xhr.responseText);
        const response = JSON.parse(xhr.responseText);
        if (response.screening_config) {
          // Parse the screening config JSON
          const screeningConfig = JSON.parse(response.screening_config);
          // Clear the container and render the plot with your D3 drawing function.
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

function updateTrainingNormField() {
  var ctcfFile = document.getElementById("ctcf_bw_path").value;
  var atacState = document.querySelector('input[name="norm_atac"]:checked').value;
  var ctcfState = document.querySelector('input[name="norm_ctcf"]:checked').value;
  var trainingContainer = document.getElementById("training-norm-container");
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

window.onload = function() {
  console.log("Window loaded. Restoring form fields...");
  restoreFormFields();
  updateEndPosition();
  toggleOptionalFields();
  updateCtcfNormalization();
  updateTrainingNormField();

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
  
  var ctcfSelectElem = document.getElementById("ctcf_bw_path");
  ctcfSelectElem.addEventListener("change", function() {
    storeFormFields();
    updateCtcfNormalization();
    updateTrainingNormField();
  });
  
  var normRadios = document.getElementsByName("norm_atac");
  for (var i = 0; i < normRadios.length; i++) {
    normRadios[i].addEventListener("change", updateTrainingNormField);
  }
  
  var ctcfNormRadios = document.getElementsByName("norm_ctcf");
  for (var i = 0; i < ctcfNormRadios.length; i++) {
    ctcfNormRadios[i].addEventListener("change", updateTrainingNormField);
  }

  // Initial population of dropdowns
  populateDropdownFromServer("atac_bw_path", "atac");
  populateDropdownFromServer("ctcf_bw_path", "ctcf");
  populateDropdownFromServer("peaks_file_path", "peaks");

  // Attach change event listeners for immediate AJAX upload of files:
  document.getElementById("atac_bw_file").addEventListener("change", function() {
    ajaxUploadFile("atac_bw_file", "atac", "atac_bw_path");
  });
  document.getElementById("ctcf_bw_file").addEventListener("change", function() {
    ajaxUploadFile("ctcf_bw_file", "ctcf", "ctcf_bw_path");
  });
  document.getElementById("peaks_file_upload").addEventListener("change", function() {
    ajaxUploadFile("peaks_file_upload", "peaks", "peaks_file_path");
  });

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

  var formElem = document.getElementById("corigami-form");
  formElem.addEventListener("submit", function(e) {
    e.preventDefault();
    if (!checkFormRequirements()) {
      return false;
    }
    var formData = new FormData(formElem);
    fetch("/", {
      method: "POST",
      body: formData,
      headers: {
        "X-Requested-With": "XMLHttpRequest"
      }
    })
    .then(response => {
      if (!response.ok) {
        throw new Error("Network response was not ok");
      }
      return response.text();
    })
    .then(html => {
      console.log("Form submission returned HTML (first 200 chars):", html.substring(0, 200));
      var container = document.getElementById("output-container");
      container.innerHTML = html;
      console.log("Executing inline scripts from AJAX response...");
      executeScripts(container);
    })
    .catch(error => {
      console.error("Error during form submission:", error);
    });
  });

  if (typeof screening_mode !== "undefined" && screening_mode === true) {
    runScreening();
  }
};
