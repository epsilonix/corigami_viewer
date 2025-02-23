const WINDOW_WIDTH = 2097152;

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
  // Read the deletion/screening radio group value.
  var dsOption = document.querySelector('input[name="ds_option"]:checked');
  var delFields = document.getElementById("deletion-fields");
  var scrFields = document.getElementById("screening-fields-2");
  var peaksContainer = document.getElementById("peaks_file_container");
  if (dsOption) {
    // Save the current selection so it persists.
    localStorage.setItem("ds_option", dsOption.value);
    if (dsOption.value === "deletion") {
      if (delFields) { delFields.style.display = "block"; }
      if (scrFields) { scrFields.style.display = "none"; }
      if (peaksContainer) { peaksContainer.style.display = "none"; }
    } else if (dsOption.value === "screening") {
      if (scrFields) { scrFields.style.display = "block"; }
      if (delFields) { delFields.style.display = "none"; }
      if (peaksContainer) { peaksContainer.style.display = "block"; }
    } else {
      if (delFields) { delFields.style.display = "none"; }
      if (scrFields) { scrFields.style.display = "none"; }
      if (peaksContainer) { peaksContainer.style.display = "none"; }
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

function updateDropdown(fileInputId, selectId, storageKey) {
  var fileInput = document.getElementById(fileInputId);
  var selectElem = document.getElementById(selectId);
  if (localStorage.getItem(storageKey)) {
    var options = JSON.parse(localStorage.getItem(storageKey));
    options.forEach(function(opt) {
      var exists = false;
      for (var i = 0; i < selectElem.options.length; i++) {
        if (selectElem.options[i].value === opt.value) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        var newOption = document.createElement("option");
        newOption.value = opt.value;
        newOption.text = opt.text;
        selectElem.appendChild(newOption);
      }
    });
  }
  fileInput.addEventListener("change", function() {
    if (fileInput.files.length > 0) {
      var filename = fileInput.files[0].name;
      var newOption = document.createElement("option");
      newOption.value = filename;
      newOption.text = filename;
      var exists = false;
      for (var i = 0; i < selectElem.options.length; i++) {
        if (selectElem.options[i].value === filename) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        selectElem.appendChild(newOption);
        var options = localStorage.getItem(storageKey) ? JSON.parse(localStorage.getItem(storageKey)) : [];
        options.push({ value: filename, text: filename });
        localStorage.setItem(storageKey, JSON.stringify(options));
      }
      selectElem.value = filename;
    }
  });
}

function runScreening() {
  if (parseInt(document.getElementById("region_start").value) < 1048576) {
    alert("Cannot run screening on start position below 1,048,576!");
    return;
  }
  document.getElementById("screening_loader").style.display = "block";
  document.getElementById("screening_image_container").style.display = "none";

  var region_chr = document.getElementById("region_chr").value;
  var model_select = document.getElementById("model_select").value;
  var region_start = document.getElementById("region_start").value;
  var region_end = document.getElementById("region_end").value;
  var perturb_width = document.getElementById("perturb_width").value;
  var step_size = document.getElementById("step_size").value;
  var ctcf_bw_path = document.getElementById("ctcf_bw_path").value;
  var atac_bw_path = document.getElementById("atac_bw_path").value;
  var output_dir = "./output_new";
  var peaks_file = document.getElementById("peaks_file_path").value;
  
  var xhr = new XMLHttpRequest();
  var params = "region_chr=" + encodeURIComponent(region_chr) +
               "&model_select=" + encodeURIComponent(model_select) +
               "&region_start=" + encodeURIComponent(region_start) +
               "&region_end=" + encodeURIComponent(region_end) +
               "&perturb_width=" + encodeURIComponent(perturb_width) +
               "&step_size=" + encodeURIComponent(step_size) +
               "&atac_bw_path=" + encodeURIComponent(atac_bw_path) +
               "&output_dir=" + encodeURIComponent(output_dir);
  params += "&ctcf_bw_path=" + encodeURIComponent(ctcf_bw_path);
  if (peaks_file) {
    params += "&peaks_file=" + encodeURIComponent(peaks_file);
  }
  
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
  xhr.open("GET", "/run_screening?" + params, true);
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

window.onload = function() {
  updateEndPosition();
  document.getElementById("region_start").addEventListener("input", updateEndPosition);
  document.getElementById("del_start").addEventListener("input", validateDeletionArea);
  document.getElementById("del_width").addEventListener("input", validateDeletionArea);
  
  updateDropdown("atac_bw_file", "atac_bw_path", "atac_bw_options");
  updateDropdown("ctcf_bw_file", "ctcf_bw_path", "ctcf_bw_options");
  updateDropdown("peaks_file_upload", "peaks_file_path", "peaks_file_options");
  
  // Modal functionality.
  var modal = document.getElementById("instructions-modal");
  var btn = document.getElementById("openModal");
  var span = document.getElementsByClassName("close")[0];
  btn.onclick = function(e) {
    e.preventDefault();
    modal.style.display = "block";
  }
  span.onclick = function() {
    modal.style.display = "none";
  }
  window.onclick = function(event) {
    if (event.target == modal) {
      modal.style.display = "none";
    }
  }
  
  // Attach event listeners for ds_option radio group.
  var dsRadios = document.getElementsByName("ds_option");
  for (var i = 0; i < dsRadios.length; i++) {
    dsRadios[i].addEventListener("change", toggleOptionalFields);
  }
  
  // On page load, check for saved ds_option.
  var savedOption = localStorage.getItem("ds_option");
  if (savedOption) {
    var radioToCheck = document.querySelector('input[name="ds_option"][value="' + savedOption + '"]');
    if (radioToCheck) {
      radioToCheck.checked = true;
    }
  }
  toggleOptionalFields();
  
  // Attach event listener for CTCF file dropdown.
  var ctcfSelectElem = document.getElementById("ctcf_bw_path");
  ctcfSelectElem.addEventListener("change", updateCtcfNormalization);
  updateCtcfNormalization();
};
