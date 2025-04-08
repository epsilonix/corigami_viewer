// d3_plots.js
/*************************************************************
 * Canvas-Based Plotting for C.Origami
 * (Supports Hi-C heatmap, signal plots with/without x-axis break,
 *  and screening column charts.)
 *************************************************************/

/**
 * SET THIS TO TRUE IF YOU WANT TO HIDE ALL X-AXIS LABELS
 * AND THE ZIG-ZAG BREAK SYMBOL AND THE TICK MARKS.
 */
const HIDE_X_AXIS_LABELS = true;

/**
 * Pixels per Mb in normal (continuous) mode.
 */
const PX_PER_MB = 200;

/**
 * Create a high-DPI (retina) canvas in a container.
 * Sets CSS width/height and scales the drawing context.
 */
function createHiResCanvas(containerSelector, width, height) {
  d3.select(containerSelector).selectAll("*").remove();
  
  const canvasSel = d3.select(containerSelector)
    .append("canvas")
    .style("width", width + "px")
    .style("height", height + "px");
  
  const canvas = canvasSel.node();
  const ratio = window.devicePixelRatio || 1;
  canvas.width = width * ratio;
  canvas.height = height * ratio;
  
  const ctx = canvas.getContext("2d");
  ctx.scale(ratio, ratio);
  return { canvas, ctx };
}


/*************************************************************
 * 1) Hi-C Heatmap (Canvas)
 *************************************************************/
function drawHiCChart(containerSelector, config) {
  // For Hi-C, we plot all the provided matrix data.
  // In deletion mode, we simply hide the x-axis if config.hideXAxis==true.
  const hideTicks = config.hideXAxis === true;
  
  const margin = { top: 5, right: 20, bottom: 5, left: 50 };
  const genomicSpan = config.xAxis.max - config.xAxis.min; // in Mb
  const fullWidth = genomicSpan * PX_PER_MB;               // in pixels
  const height = config.chart.height || 250;
  const canvasWidth  = fullWidth + margin.left + margin.right;
  const canvasHeight = height + margin.top + margin.bottom;
  
  const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
  ctx.translate(margin.left, margin.top);
  
  const n_cols = config.n_cols;
  const maxRow = d3.max(config.series[0].data, d => d[1]);
  const n_rows = maxRow + 1;
  
  // Map matrix column indices to genomic coordinates (Mb)
  const colIndexToMb = d3.scaleLinear()
    .domain([0, n_cols])
    .range([config.xAxis.min, config.xAxis.max]);
  
  // Map genomic coordinates (Mb) to pixels
  const xScale = d3.scaleLinear()
    .domain([config.xAxis.min, config.xAxis.max])
    .range([0, fullWidth]);
  
  // Row scale: row index to vertical pixel
  const rowScale = d3.scaleLinear()
    .domain([0, n_rows])
    .range([0, height]);
  
  // Color scale for heatmap cells
  const colorScale = d3.scaleSequential(d3.interpolateReds)
    .domain([config.colorAxis.min, config.colorAxis.max]);
  
  // Draw all heatmap cells
  config.series[0].data.forEach(d => {
    const col = d[0], row = d[1], val = d[2];
    const mbLeft = colIndexToMb(col);
    const mbRight = colIndexToMb(col + 1);
    const xLeft = xScale(mbLeft);
    const xRight = xScale(mbRight);
    const cellW = xRight - xLeft;
    const yTop = rowScale(row);
    const cellH = rowScale(row + 1) - yTop;
    
    ctx.fillStyle = colorScale(val);
    ctx.fillRect(xLeft, yTop, cellW, cellH);
  });
  
  // If we are not hiding the axis entirely
  if (!hideTicks) {
    ctx.save();
    ctx.translate(0, height);
    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(fullWidth, 0);
    ctx.strokeStyle = "#000";
    ctx.stroke();
    
    if (!HIDE_X_AXIS_LABELS) {
      // Draw the tick marks & labels
      const xTicks = xScale.ticks(10);
      ctx.font = "10px sans-serif";
      ctx.fillStyle = "#000";
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
      xTicks.forEach(tick => {
        const xPos = xScale(tick);
        ctx.beginPath();
        ctx.moveTo(xPos, 0);
        ctx.lineTo(xPos, 6);
        ctx.stroke();
        ctx.fillText(d3.format(".2f")(tick), xPos, 6);
      });
    }
    ctx.restore();
    
    // x-axis title
    if (config.xAxis.title && !HIDE_X_AXIS_LABELS) {
      ctx.save();
      ctx.textAlign = "center";
      ctx.textBaseline = "bottom";
      ctx.font = "12px sans-serif";
      ctx.fillStyle = "#000";
      ctx.fillText(config.xAxis.title, fullWidth / 2, height + margin.bottom - 5);
      ctx.restore();
    }
  }
}


/**
 * Draws a column (bar) chart using the given config.
 * - If config.xAxis.axisBreak is not true, it does the standard single-domain approach.
 * - If config.xAxis.axisBreak === true, and leftDomain/rightDomain are provided,
 *   it uses a two-scale approach to skip the middle region (like a real axis break).
 *
 * The rest of the code is your original logic for binning, screening shading, Y-axis, etc.
 */
function drawColumnChart(containerSelector, config) {
  console.log("=== drawColumnChart START ===");
  console.log("Container selector:", containerSelector);
  console.log("Chart config:", config);

  // -----------------------------------------------------------------
  // Constants & defaults (unchanged)
  // -----------------------------------------------------------------
  const PX_PER_MB = 200;
  const MIN_BARS_PER_MB = 30;
  const HIDE_X_AXIS_LABELS = true;
  const margin = { top: 20, right: 20, bottom: 40, left: 50 };

  // -----------------------------------------------------------------
  // 1) Basic validation
  // -----------------------------------------------------------------
  if (!config.xAxis || typeof config.xAxis.min !== "number" || typeof config.xAxis.max !== "number") {
    console.warn("Missing or invalid xAxis.min/max in config.xAxis:", config.xAxis);
    return;
  }
  const xMin = config.xAxis.min;
  const xMax = config.xAxis.max;
  const domainSpanMb = xMax - xMin;

  // We'll store these in a closure so we can detect axisBreak
  const hasBreak   = !!config.xAxis.axisBreak;           // NEW
  const leftDomain = config.xAxis.leftDomain || [];       // e.g. [xMin, delStartMb]
  const rightDomain= config.xAxis.rightDomain || [];      // e.g. [delEndMb, xMax]

  // -----------------------------------------------------------------
  // 2) Prepare data array (binning/screening) EXACTLY as before
  // -----------------------------------------------------------------
  let originalData = (config.series && config.series[0] && config.series[0].data) || [];
  console.log("Data length:", originalData.length);
  if (!originalData.length) {
    // Show "No data" message
    const fallback = document.querySelector(containerSelector);
    if (fallback) fallback.innerHTML = "<p style='color:red;'>No data available.</p>";
    console.warn("Column chart: data array is empty.");
    return;
  }

  // Sort data by x
  originalData.sort((a, b) => a[0] - b[0]);
  console.log("Sorted data (first few points):", originalData.slice(0, 5));

  // Screening or binning
  let data = [];
  if (config.screening) {
    data = originalData;
  } else {
    const minBars = Math.round(domainSpanMb * MIN_BARS_PER_MB);
    console.log("minBars needed =", minBars, "domainSpanMb:", domainSpanMb);
    if (originalData.length >= minBars) {
      data = originalData;
    } else {
      console.log("Binning to produce", minBars, "bars");
      data = binAndUpsampleData(originalData, xMin, xMax, minBars);
    }
  }

  console.log("Final data length after binning or screening mode:", data.length);

  // -----------------------------------------------------------------
  // 3) Y-domain logic (same as before)
  // -----------------------------------------------------------------
  const rawYMin = d3.min(data, d => d[1]);
  const rawYMax = d3.max(data, d => d[1]);
  let yMin = (typeof rawYMin === "number") ? rawYMin : 0;
  let yMax = (typeof rawYMax === "number") ? rawYMax : 0;

  let domainMin, domainMax, axisValue;
  if (yMax < 0) {
    // entire negative
    domainMin = yMin; 
    domainMax = yMax;
    axisValue = yMin;
  } else if (yMin > 0) {
    // entire positive
    domainMin = 0;
    domainMax = yMax;
    axisValue = 0;
  } else {
    // crosses zero
    domainMin = yMin;
    domainMax = yMax;
    axisValue = 0;
  }

  // We'll define yScale later inside each approach
  // but let's figure out the chart height
  const height = config.chart.height || 150;

  // -----------------------------------------------------------------
  // 4) Check if axisBreak => do the "two-scale" approach
  // -----------------------------------------------------------------
  if (hasBreak && leftDomain.length === 2 && rightDomain.length === 2) {
    // =========== NEW: Two-scale approach =============
    // 4a) Calculate spans
    const leftSpan  = leftDomain[1]  - leftDomain[0];
    const rightSpan = rightDomain[1] - rightDomain[0];
    const totalSpan = leftSpan + rightSpan;
    const gapPx     = 0; // If you want a small gap, e.g. 20
    // Convert to px
    const leftWidthPx  = leftSpan  * PX_PER_MB;
    const rightWidthPx = rightSpan * PX_PER_MB;
    const totalWidthPx = leftWidthPx + rightWidthPx + gapPx;

    // 4b) Create hi-res canvas
    const canvasWidth  = totalWidthPx + margin.left + margin.right;
    const canvasHeight = height       + margin.top  + margin.bottom;
    const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
    if (!ctx) {
      console.error("No 2D context for 2-scale approach");
      return;
    }
    ctx.translate(margin.left, margin.top);

    // 4c) Build 2 x‐scales
    const scaleLeft = d3.scaleLinear()
      .domain([ leftDomain[0],  leftDomain[1] ])
      .range([0, leftWidthPx]);
    const scaleRight = d3.scaleLinear()
      .domain([ rightDomain[0], rightDomain[1] ])
      .range([0, rightWidthPx]);

    // 4d) Filter data into left portion / right portion
    const leftData  = data.filter(d => d[0] >= leftDomain[0]  && d[0] <= leftDomain[1]);
    const rightData = data.filter(d => d[0] >= rightDomain[0] && d[0] <= rightDomain[1]);

    // 4e) Y scale (common for both)
    const yScale = d3.scaleLinear()
      .domain([domainMin, domainMax])
      .nice()
      .range([height, 0]);
    const axisPx = yScale(axisValue);

    // 4f) If screening => shading
    if (config.screening) {
      ctx.save();
      ctx.fillStyle = "rgba(200, 200, 200, 0.4)";
    
      // We'll figure out if/where we shade the left side (start) and the right side (end).
      // 'xMin' and 'xMax' are in Mb, from config.xAxis.min and config.xAxis.max
      // 'chromEndMb' is the entire chromosome length in Mb (if known).
      let chromEndMb = null;
      if (window.chromosomeLengths && config.chromName) {
        let bpLen = window.chromosomeLengths[config.chromName];
        if (bpLen) {
          chromEndMb = bpLen / 1e6; // convert to Mb
        }
      }
    
      // We’ll define two booleans and numeric boundaries:
      let doLeftShading  = false;
      let doRightShading = false;
      let leftLimitMb    = xMin;  // up to where we shade on the left
      let rightLimitMb   = xMax;  // from where we shade on the right
    

      console.log("screening?", config.screening);
      console.log("isChimeric?", config.isChimeric);
      console.log("xMin, xMax:", xMin, xMax);
      console.log("chromName:", config.chromName);
      console.log("chromEndMb computed:", chromEndMb);
      
      // 1) Always shade first & last 1Mb if chimeric:
      if (config.isChimeric) {
        doLeftShading  = true;
        doRightShading = true;
        leftLimitMb    = xMin + 1.0;     // shade from xMin => (xMin + 1Mb)
        rightLimitMb   = xMax - 1.0;     // shade from (xMax - 1Mb) => xMax
      }
      // 2) If NOT chimeric, see if xMin < 1 => that means we’re actually near the chromosome start
      //    so we shade from xMin => 1.0
      //    Similarly, if xMax > chromEndMb - 1 => shade that portion
      else if (chromEndMb) {
        if (xMin < 1.0) {
          doLeftShading = true;
          leftLimitMb   = Math.min(xMax, 1.0); // never exceed xMax
        }
        if (xMax > (chromEndMb - 1.0)) {
          doRightShading = true;
          rightLimitMb   = Math.max(xMin, chromEndMb - 1.0);
        }
      }
    
      // Now do the actual shading, converting Mb => pixels
      if (doLeftShading) {
        const shadeEndPx = xScale(leftLimitMb);
        if (shadeEndPx > 0) {
          ctx.fillRect(0, 0, shadeEndPx, height);
        }
      }
      if (doRightShading) {
        const shadeStartPx = xScale(rightLimitMb);
        if (shadeStartPx < fullWidth) {
          ctx.fillRect(shadeStartPx, 0, fullWidth - shadeStartPx, height);
        }
      }
    
      ctx.restore();
    }
    
    // 4g) Decide bar spacing on left side
    let leftCount = leftData.length;
    let leftSpacing = (leftCount > 1) ? (leftWidthPx / (leftCount - 1)) : leftWidthPx;
    let barMultiplier = config.screening ? 0.02 : 0.05;
    let leftBarWidth = leftSpacing * barMultiplier;
    leftBarWidth = Math.min(Math.max(leftBarWidth, 1), 20);

    // 4h) Draw left bars
    leftData.forEach(d => {
      const xVal = d[0], yVal = d[1];
      const xPx = scaleLeft(xVal) - leftBarWidth/2;
      const yValPx = yScale(yVal);
      const top = Math.min(yValPx, axisPx);
      const bottom = Math.max(yValPx, axisPx);
      const h = bottom - top;
      ctx.fillStyle = (config.series[0].color || "#9FC0DE");
      ctx.fillRect(xPx, top, leftBarWidth, h);
    });

    // Baseline for left chunk
    ctx.beginPath();
    ctx.moveTo(0, axisPx);
    ctx.lineTo(leftWidthPx, axisPx);
    ctx.strokeStyle = "#000";
    ctx.stroke();

    // 4i) Right chunk
    // Move over
    ctx.save();
    ctx.translate(leftWidthPx + gapPx, 0);

    let rightCount = rightData.length;
    let rightSpacing= (rightCount > 1) ? (rightWidthPx / (rightCount - 1)) : rightWidthPx;
    let rightBarWidth = rightSpacing * barMultiplier;
    rightBarWidth = Math.min(Math.max(rightBarWidth,1), 20);

    rightData.forEach(d => {
      const xVal = d[0], yVal = d[1];
      const xPx = scaleRight(xVal) - rightBarWidth/2;
      const yValPx = yScale(yVal);
      const top = Math.min(yValPx, axisPx);
      const bottom = Math.max(yValPx, axisPx);
      const h = bottom - top;
      ctx.fillStyle = (config.series[0].color || "#9FC0DE");
      ctx.fillRect(xPx, top, rightBarWidth, h);
    });

    // Baseline for right chunk
    ctx.beginPath();
    ctx.moveTo(0, axisPx);
    ctx.lineTo(rightWidthPx, axisPx);
    ctx.stroke();

    ctx.restore(); // end right chunk

    // 4j) Y-axis line & ticks (on far left)
    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(0, height);
    ctx.strokeStyle = "#000";
    ctx.stroke();
    // y ticks
    ctx.save();
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    ctx.font = "10px sans-serif";
    yScale.ticks(5).forEach(t => {
      const yPos = yScale(t);
      ctx.beginPath();
      ctx.moveTo(0, yPos);
      ctx.lineTo(-5, yPos);
      ctx.stroke();
      ctx.fillText(t, -7, yPos);
    });
    ctx.restore();

    // 4k) (Optional) Zigzag break line
    if (!HIDE_X_AXIS_LABELS) {
      ctx.save();
      ctx.translate(leftWidthPx, axisPx);
      ctx.beginPath();
      ctx.moveTo(-3, 0);
      ctx.lineTo(0, 6);
      ctx.lineTo(3, 0);
      ctx.strokeStyle = "#000";
      ctx.lineWidth = 1.5;
      ctx.stroke();
      ctx.restore();
    }

    // 4l) X-axis label & chart title if you want them
    if (config.xAxis.title && !HIDE_X_AXIS_LABELS) {
      ctx.textAlign = "center";
      ctx.textBaseline = "bottom";
      ctx.font = "12px sans-serif";
      ctx.fillText(config.xAxis.title, totalWidthPx/2, height + margin.bottom - 5);
    }
    if (config.title && config.title.text) {
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
      ctx.font = "13px sans-serif";
      ctx.fillText(config.title.text, totalWidthPx/2, -margin.top/2);
    }

    console.log("=== drawColumnChart END (two-scale) ===");
    return;  // Done with the axisBreak case
  }

  // -----------------------------------------------------------------
  // 5) Otherwise => Single-scale approach (original code EXACTLY)
  // -----------------------------------------------------------------

  // 5a) Compute fullWidth from domainSpanMb
  let fullWidth = domainSpanMb * PX_PER_MB;
  const canvasWidth  = fullWidth + margin.left + margin.right;
  const canvasHeight = height    + margin.top  + margin.bottom;

  // Create hi-res canvas
  const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
  if (!ctx) {
    console.error("No 2D context returned; createHiResCanvas might have failed?");
    return;
  }
  ctx.translate(margin.left, margin.top);

  // 5b) Build the single xScale
  const xScale = d3.scaleLinear()
    .domain([xMin, xMax])
    .range([0, fullWidth]);

  // 5c) Build the yScale & axisPx
  const yScale = d3.scaleLinear()
    .domain([domainMin, domainMax])
    .nice()
    .range([height, 0]);
  const axisPx = yScale(axisValue);

  // 5d) If screening => shading skip region (unchanged logic)
  if (config.screening) {
    ctx.save();
    ctx.fillStyle = "rgba(200, 200, 200, 0.4)";
    if (config.isChimeric) {
      console.log("COLUMN PLOT SKIP REGION: isChimeric:", config.isChimeric);
      const leftSkipEndX = xScale(Math.min(xMax, xMin + 1.0));
      if (leftSkipEndX > 0) {
        ctx.fillRect(0, 0, leftSkipEndX, height);
      }
      const rightSkipStartX = xScale(Math.max(xMin, xMax - 1.0));
      if (rightSkipStartX < fullWidth) {
        const w = fullWidth - rightSkipStartX;
        ctx.fillRect(rightSkipStartX, 0, w, height);
      }
    } else {
      // normal skip => first/last 1 Mb
      let chromEndMb = null;
      const regionChrElem = document.getElementById("region_chr1");
      if (regionChrElem && window.chromosomeLengths) {
        let regionChr = regionChrElem.value;
        let bpLen = window.chromosomeLengths[regionChr];
        if (bpLen) {
          chromEndMb = bpLen / 1e6;
        }
      }
      if (chromEndMb) {
        const leftSkipEndX = xScale(Math.max(xMin, 1.0));
        if (leftSkipEndX > 0) {
          ctx.fillRect(0, 0, leftSkipEndX, height);
        }
        const rightSkipStartX = xScale(Math.min(xMax, chromEndMb - 1.0));
        if (rightSkipStartX < fullWidth) {
          const w = fullWidth - rightSkipStartX;
          ctx.fillRect(rightSkipStartX, 0, w, height);
        }
      }
    }
    ctx.restore();
  }

  // 5e) Decide bar spacing
  let dataCount = data.length;
  let spacingPixels = (dataCount > 1) ? (fullWidth / (dataCount - 1)) : fullWidth;
  let barMultiplier = config.screening ? 0.02 : 0.05;
  let barWidth = spacingPixels * barMultiplier;
  barWidth = Math.min(barWidth, 20);
  barWidth = Math.max(barWidth, 1);

  // 5f) Draw bars
  data.forEach(d => {
    const xVal = d[0], yVal = d[1];
    const xPx = xScale(xVal) - barWidth/2;
    const yValPx = yScale(yVal);
    const top    = Math.min(yValPx, axisPx);
    const bottom = Math.max(yValPx, axisPx);
    const h      = bottom - top;
    ctx.fillStyle = config.series[0].color || "#9FC0DE";
    ctx.fillRect(xPx, top, barWidth, h);
  });

  // 5g) Horizontal axis line
  ctx.beginPath();
  ctx.moveTo(0, axisPx);
  ctx.lineTo(fullWidth, axisPx);
  ctx.strokeStyle = "#000";
  ctx.stroke();

  // 5h) Y-axis line & ticks
  ctx.save();
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(0, height);
  ctx.strokeStyle = "#000";
  ctx.stroke();
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  ctx.font = "10px sans-serif";
  yScale.ticks(5).forEach(t => {
    const yPos = yScale(t);
    ctx.beginPath();
    ctx.moveTo(0, yPos);
    ctx.lineTo(-5, yPos);
    ctx.stroke();
    ctx.fillText(t, -7, yPos);
  });
  ctx.restore();

  // 5i) X-axis label & chart title
  if (config.xAxis.title && !HIDE_X_AXIS_LABELS) {
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.font = "12px sans-serif";
    const titleText = (typeof config.xAxis.title === "string")
      ? config.xAxis.title
      : (config.xAxis.title.text || "");
    ctx.fillText(titleText, fullWidth / 2, height + margin.bottom - 5);
  }
  if (config.title && config.title.text) {
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.font = "13px sans-serif";
    ctx.fillText(config.title.text, fullWidth / 2, -margin.top / 2);
  }

  console.log("=== drawColumnChart END (single scale) ===");
}



function binAndUpsampleData(originalData, xMin, xMax, desiredBins) {
  let binned = [];
  const binWidth = (xMax - xMin) / desiredBins;

  // We'll iterate from i=0..desiredBins-1
  // For each bin => left = xMin + i*binWidth, right= left+binWidth
  let dataIndex = 0;  // which point in originalData
  const nData = originalData.length;

  for (let i = 0; i < desiredBins; i++) {
    const binLeft  = xMin + i * binWidth;
    const binRight = binLeft + binWidth;

    // collect points that fall in [binLeft..binRight)
    let sumY = 0, countY = 0;
    // we can loop while originalData[dataIndex].x is < binRight
    while (dataIndex < nData) {
      const xVal = originalData[dataIndex][0];
      const yVal = originalData[dataIndex][1];
      if (xVal >= binRight) {
        break; // move to next bin
      }
      if (xVal >= binLeft && xVal < binRight) {
        sumY += yVal;
        countY++;
      }
      dataIndex++;
    }

    let avgY = 0;
    if (countY > 0) {
      avgY = sumY / countY;
    } else {
      // no points => set 0 or last known => we do 0
      // or we could do "carry last bin's average"
      avgY = 0;
    }

    // We'll represent the bin's x as the midpoint:
    const binMid = binLeft + binWidth/2;
    binned.push([binMid, avgY]);
  }

  return binned;
}






const BUFFER_BP = 500000; // used for row assignment

function drawGeneTrackChart(selector, config) {
  console.log("Drawing gene track with config:", config);

  // Parse region boundaries
  const regionStart = parseInt(config.region1.start, 10);
  const regionEnd = parseInt(config.region1.end, 10);
  const regionSpan = regionEnd - regionStart;
  console.log("drawGeneTrackChart: regionStart, regionEnd, regionSpan:", regionStart, regionEnd, regionSpan);

  // Determine deletion mode (if applicable)
  const inDeletionMode = config.deletionStart != null && config.deletionEnd != null;

  // Calculate full width and define xScale.
  let fullWidth, xScale;
  if (config.xAxis && config.xAxis.axisBreak && inDeletionMode) {
    const leftDomainBp = [regionStart, config.deletionStart];
    const rightDomainBp = [config.deletionEnd, regionEnd];
    const leftSpan = config.deletionStart - regionStart;
    const rightSpan = regionEnd - config.deletionEnd;
    const totalSpan = leftSpan + rightSpan;
    fullWidth = (totalSpan / 1e6) * PX_PER_MB;
    const gapPx = 0;
    const leftWidthPx = (leftSpan / totalSpan) * fullWidth;
    const rightWidthPx = (rightSpan / totalSpan) * fullWidth;
    const scaleLeft = d3.scaleLinear().domain(leftDomainBp).range([0, leftWidthPx]);
    const scaleRight = d3.scaleLinear().domain(rightDomainBp).range([0, rightWidthPx]);
    xScale = function(xVal) {
      if (xVal <= config.deletionStart) return scaleLeft(xVal);
      if (xVal >= config.deletionEnd) return leftWidthPx + gapPx + scaleRight(xVal);
      return NaN;
    };
  } else {
    fullWidth = (regionSpan / 1e6) * PX_PER_MB;
    xScale = d3.scaleLinear().domain([regionStart, regionEnd]).range([0, fullWidth]);
  }

  // Set margins and create the SVG container.
  const margin = config.margin || { top: 0, right: 0, bottom: 5, left: 50 };
  const desiredHeight = config.chart.height || 50;
  const container = d3.select(selector);
  container.select("svg").remove();
  const svg = container.append("svg")
                       .attr("width", fullWidth + margin.left + margin.right);
  const g = svg.append("g")
               .attr("transform", `translate(${margin.left},${margin.top})`);

  // Get gene list: use config.genes if provided; otherwise load from annotationFile.
  let genesPromise;
  if (config.genes && config.genes.length) {
    genesPromise = Promise.resolve(config.genes);
  } else if (config.annotationFile) {
    genesPromise = d3.text(config.annotationFile).then(function(text) {
      let lines = text.split("\n").filter(line => line.trim().length > 0);
      return lines.map(function(line) {
        let parts = line.split("\t");
        let coordParts = parts[2].split(/[:\-]/); // e.g., ["chr19", "58345178", "58353492"]
        return {
          gene: parts[0],
          ensembl: parts[1],
          chrom: coordParts[0],
          start: +coordParts[1],
          end: +coordParts[2],
          strand: parts[3],
          type: parts[4]
        };
      });
    });
  } else {
    genesPromise = Promise.resolve([]);
  }

  genesPromise.then(function(genes) {
    // Filter genes for the current region and protein-coding ones
    genes = genes.filter(function(d) {
      return d.chrom === config.region1.chr &&
             d.end >= regionStart &&
             d.start <= regionEnd &&
             d.type === "protein_coding" &&
             !/^\d/.test(d.gene);
    });
    console.log("Number of protein coding genes after filtering:", genes.length);

    if (inDeletionMode) {
      genes = genes.filter(function(d) {
        return (d.end < config.deletionStart || d.start > config.deletionEnd);
      });
    }

    // Sort and assign rows (using BUFFER_BP from your existing code)
    genes.sort((a, b) => a.start - b.start);
    let geneRows = [];
    genes.forEach(function(gene) {
      let placed = false;
      for (let i = 0; i < geneRows.length; i++) {
        if (gene.start > geneRows[i] + BUFFER_BP) {
          gene.row = i;
          geneRows[i] = gene.end;
          placed = true;
          break;
        }
      }
      if (!placed) {
        gene.row = geneRows.length;
        geneRows.push(gene.end);
      }
    });
    const numRows = geneRows.length;
    const MIN_ROW_HEIGHT = 15;
    const adjustedHeight = Math.max(desiredHeight, numRows * MIN_ROW_HEIGHT);
    const rowHeight = adjustedHeight / numRows;
    const canvasHeight = adjustedHeight + margin.top + margin.bottom;
    svg.attr("height", canvasHeight);

    // Draw gene lines
    g.selectAll("line.gene")
     .data(genes)
     .enter()
     .append("line")
     .attr("class", "gene")
     .attr("x1", d => xScale(Math.max(d.start, regionStart)))
     .attr("x2", d => xScale(Math.min(d.end, regionEnd)))
     .attr("y1", d => d.row * rowHeight + rowHeight / 2)
     .attr("y2", d => d.row * rowHeight + rowHeight / 2)
     .attr("stroke", d => d.strand === "+" ? "#9FC0DE" : "#F2C894")
     .attr("stroke-width", Math.min(rowHeight * 0.2, 1.5))
     .attr("stroke-linecap", "round");

    // Draw gene arrow markers
    g.selectAll("path.gene-arrow")
     .data(genes.filter(d => d.strand === "+" || d.strand === "-"))
     .enter()
     .append("path")
     .attr("class", "gene-arrow")
     .attr("d", function(d) {
       const arrowWidth = rowHeight * 0.5;
       const arrowHeight = rowHeight * 0.4;
       const centerY = d.row * rowHeight + rowHeight / 2;
       if (d.strand === "+") {
         const tipX = xScale(Math.min(d.end, regionEnd));
         return `M${tipX},${centerY} L${tipX - arrowWidth},${centerY - arrowHeight/2} L${tipX - arrowWidth},${centerY + arrowHeight/2} Z`;
       } else {
         const tipX = xScale(Math.max(d.start, regionStart));
         return `M${tipX},${centerY} L${tipX + arrowWidth},${centerY - arrowHeight/2} L${tipX + arrowWidth},${centerY + arrowHeight/2} Z`;
       }
     })
     .attr("fill", d => d.strand === "+" ? "#87CEEB" : "#FFA500")
     .attr("stroke", d => d.strand === "+" ? "#87CEEB" : "#FFA500")
     .attr("stroke-width", 1);

    // Draw gene labels
    g.selectAll("text.gene-label")
    .data(genes)
    .enter()
    .append("text")
    .attr("class", "gene-label")
    .attr("x", function(d) {
      const x1 = xScale(Math.max(d.start, regionStart));
      const x2 = xScale(Math.min(d.end, regionEnd));
      let center = (x1 + x2) / 2;

      // Nudge the label if it's too close to the left or right edge
      const buffer = 16;              // how many pixels to pad from edges
      if (center < buffer) {
        center = x1 + buffer;
      } else if (center > fullWidth - buffer) {
        center = x2 - buffer;
      }

      return center;
    })
    .attr("y", d => d.row * rowHeight + rowHeight / 2 + 3)
    .attr("font-size", () => (Math.max(Math.min(rowHeight * 0.8, 16), 12) / 2) + 2)
    .attr("text-anchor", "middle")
    .text(d => d.gene);


    // Collapse view if there are more than 5 rows.
    const maxVisibleRows = 5;
    if (numRows > maxVisibleRows) {
      // Set container height to show only 5 full rows.
      const visibleHeight = maxVisibleRows * rowHeight + margin.top + margin.bottom - 5;
      container.style("height", visibleHeight + "px")
              .style("overflow-y", "hidden");
      // Append a clickable toggle in the upper left corner.
      if (container.select("#gene-track-toggle").empty()) {
        container.append("div")
          .attr("id", "gene-track-toggle")
          .style("position", "absolute")
          .style("top", "0px")
          .style("left", "15px")
          .style("cursor", "pointer")
          .style("width", "16px")
          .style("height", "16px")
          .style("transition", "transform 0.3s ease")
          .html(`<img src="static/caret-down-solid.svg" style="width:16px;height:16px;fill:#333;" />`)
          .on("click", function() {
            const currentHeight = parseFloat(container.style("height"));
            if (currentHeight < canvasHeight) {
              container.style("height", canvasHeight + "px");
              d3.select(this).style("transform", "rotate(180deg)");
            } else {
              container.style("height", visibleHeight + "px");
              d3.select(this).style("transform", "rotate(0deg)");
            }
          });
      }
    }

  }).catch(function(error) {
    console.error("Error loading gene annotation data:", error);
  });
}


/*************************************************************
* Helper: createCanvas
* Creates a canvas, sets the devicePixelRatio scaling,
* and returns { ctx, canvasWidth, canvasHeight, ratio }.
*************************************************************/
function createCanvas(container, totalWidthPx, axisHeight, margin) {
  const canvasWidth  = totalWidthPx + margin.left + margin.right;
  const canvasHeight = axisHeight  + margin.top  + margin.bottom;

  const canvas = document.createElement("canvas");
  canvas.style.width  = canvasWidth  + "px";
  canvas.style.height = canvasHeight + "px";
  container.appendChild(canvas);

  const ratio = window.devicePixelRatio || 1;
  canvas.width  = canvasWidth  * ratio;
  canvas.height = canvasHeight * ratio;
  const ctx = canvas.getContext("2d");
  ctx.scale(ratio, ratio);
  ctx.translate(margin.left, margin.top);

  return { ctx, canvasWidth, canvasHeight, ratio };
}

/*************************************************************
* Helper: drawPastelBar
* Draws a colored rectangle for a region, plus a label on it.
*************************************************************/
function drawPastelBar(ctx, x, widthPx, barHeight, color, chromLabel) {
  ctx.save();
  ctx.fillStyle = color;

  // Example: pick text color. You can refine if needed:
  let textColor = "#000";  // default is black
  if (color.toLowerCase() === "#779ecc" || color.toLowerCase() === "#87ceeb") {
    textColor = "#fff";  // white if background is a darker pastel
  }

  ctx.fillRect(x, -barHeight, widthPx, barHeight);

  // Use the same 10px font size used for ticks:
  ctx.fillStyle = textColor;
  ctx.font = "10px sans-serif";
  ctx.textAlign = "left";
  ctx.textBaseline = "middle";
  ctx.fillText(chromLabel, x + 5, -barHeight / 2);

  ctx.restore();
}

/*************************************************************
* Helper: drawBaseline
* Simple line from xStart to xEnd at y=0
*************************************************************/
function drawBaseline(ctx, xStart, xEnd) {
  ctx.beginPath();
  ctx.moveTo(xStart, 0);
  ctx.lineTo(xEnd, 0);
  ctx.strokeStyle = "#000";
  ctx.stroke();
}

/*************************************************************
* Helper: buildTickVals
* Generates an array of tick positions depending on reversed or not.
*************************************************************/
function buildTickVals(minMb, maxMb, stepMb, reversed=false) {
  const ticks = [];
  if (!reversed) {
    for (let t = minMb; t <= maxMb + 1e-9; t += stepMb) {
      // ensure numerical rounding issues don’t accumulate
      ticks.push(parseFloat(t.toFixed(9)));
    }
  } else {
    for (let t = maxMb; t >= minMb - 1e-9; t -= stepMb) {
      ticks.push(parseFloat(t.toFixed(9)));
    }
  }
  return ticks;
}

/*************************************************************
* Helper: drawTicksBelow
* Draw ticks from the array on the baseline (y=0), 
* with a 6px tick length and labels at y=8 below the line.
*************************************************************/
function drawTicksBelow(ctx, scale, tickVals, labelPrecision=1) {
  ctx.save();
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.font = "10px sans-serif";

  tickVals.forEach(t => {
    const x = scale(t);
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.lineTo(x, 6);
    ctx.strokeStyle = "#000";
    ctx.stroke();

    ctx.fillStyle = "#000";
    ctx.fillText(t.toFixed(labelPrecision), x, 8);
  });

  ctx.restore();
}

/*************************************************************
* drawCustomAxis
* Draws a synthetic axis with up to two segments (for chimeric),
* or a single segment with possible deletion break, *and*
* now correctly handles isReversed = true.
*************************************************************/
function drawCustomAxis(containerSelector, config) {
  // Basic checks
  if (!config?.region1) {
    console.warn("drawCustomAxis: missing config.region1");
    return;
  }

  const PX_PER_MB = 200;
  const margin = { top: 20, right: 30, bottom: 0, left: 50 };
  const axisHeight = 50; 

  // Slightly reduce the pastel bar height
  const barHeight  = 12; 
  const tickStep   = 0.2; 

  const container  = document.querySelector(containerSelector);
  if (!container) {
    console.error("drawCustomAxis: container not found:", containerSelector);
    return;
  }
  // Clear any old axis
  while (container.firstChild) container.removeChild(container.firstChild);

  // Region1
  const r1 = config.region1;
  // Possibly region2
  const r2 = config.region2 || null;

  // Deletion check
  const hasDeletion = (
    !r2 &&
    config.deletionStartMb != null &&
    config.deletionEndMb   != null &&
    config.deletionStartMb < config.deletionEndMb
  );

  /*****************************************************
   * CASE 1) Chimeric => region1 + region2
   *****************************************************/
  if (r2) {
    console.log("custom x axis: chimeric with possibly reversed segments");
    
    // Region1 domain
    const r1Domain = r1.isReversed
      ? [r1.endMb, r1.startMb]
      : [r1.startMb, r1.endMb];
    const c1Mb = Math.abs(r1.endMb - r1.startMb);

    // Region2 domain
    const r2Domain = r2.isReversed
      ? [r2.endMb, r2.startMb]
      : [r2.startMb, r2.endMb];
    const c2Mb = Math.abs(r2.endMb - r2.startMb);

    // total length in Mb
    const totalMb = c1Mb + c2Mb;
    const totalWidthPx = totalMb * PX_PER_MB;

    // create the canvas
    const { ctx } = createCanvas(container, totalWidthPx, axisHeight, margin);

    // chunk1 scale
    const scale1 = d3.scaleLinear()
      .domain(r1Domain)
      .range([0, c1Mb * PX_PER_MB]);

    // chunk2 scale
    const scale2 = d3.scaleLinear()
      .domain(r2Domain)
      .range([0, c2Mb * PX_PER_MB]);

    const chunk1WidthPx = c1Mb * PX_PER_MB;
    const chunk2WidthPx = c2Mb * PX_PER_MB;
    const chunk2Offset  = chunk1WidthPx;  // no gap between the two bars

    // Draw chunk1 bar & baseline
    drawPastelBar(ctx, 0, chunk1WidthPx, barHeight, "#F2C894", r1.chrom);
    drawBaseline(ctx, 0, chunk1WidthPx);

    // Build ticks for chunk1
    const r1Start = Math.min(r1.startMb, r1.endMb);
    const r1End   = Math.max(r1.startMb, r1.endMb);
    let r1Ticks   = buildTickVals(r1Start, r1End, tickStep, r1.isReversed);
    r1Ticks.pop(); // optional: sometimes you may want to pop the last tick

    // Build ticks for chunk2
    const r2Start = Math.min(r2.startMb, r2.endMb);
    const r2End   = Math.max(r2.startMb, r2.endMb);
    let r2Ticks   = buildTickVals(r2Start, r2End, tickStep, r2.isReversed);

    console.log("r1Ticks:", r1Ticks);
    console.log("is r1 reversed?", r1.isReversed);
    drawTicksBelow(ctx, scale1, r1Ticks);

    // chunk2 in separate transform
    ctx.save();
    ctx.translate(chunk2Offset, 0);

    // chunk2 bar & baseline
    drawPastelBar(ctx, 0, chunk2WidthPx, barHeight, "#9FC0DE", r2.chrom);
    drawBaseline(ctx, 0, chunk2WidthPx);

    console.log("r2Ticks:", r2Ticks);
    console.log("is r2 reversed?", r2.isReversed);
    drawTicksBelow(ctx, scale2, r2Ticks);

    ctx.restore();
    return;
  }

  /*****************************************************
   * CASE 2) Deletion => single-chrom with break
   *****************************************************/
  if (hasDeletion) {
    console.log("custom x axis: single-chrom + deletion break");
    const startMb = r1.startMb;
    const endMb   = r1.endMb;
    const delStartMb = config.deletionStartMb;
    const delEndMb   = config.deletionEndMb;
    
    // Is region reversed?
    const reversed = !!r1.isReversed; 
    // We'll treat chunk1 => [start..delStart], chunk2 => [delEnd..end]

    // chunk1 domain
    let c1Min = startMb, c1Max = delStartMb;
    if (reversed) {
      c1Min = delStartMb; 
      c1Max = startMb;
    }
    const c1Mb = Math.abs(delStartMb - startMb);

    // chunk2 domain
    let c2Min = delEndMb, c2Max = endMb;
    if (reversed) {
      c2Min = endMb; 
      c2Max = delEndMb;
    }
    const c2Mb = Math.abs(endMb - delEndMb);

    const chunk1Width = c1Mb * PX_PER_MB;
    const chunk2Width = c2Mb * PX_PER_MB;
    const totalWidthPx = chunk1Width + chunk2Width;

    const { ctx } = createCanvas(container, totalWidthPx, axisHeight, margin);

    // chunk1 scale
    const scale1 = d3.scaleLinear()
      .domain([c1Min, c1Max])
      .range([0, chunk1Width]);

    // chunk2 scale
    const scale2 = d3.scaleLinear()
      .domain([c2Min, c2Max])
      .range([0, chunk2Width]);

    // 1) chunk1 bar, baseline, ticks
    drawPastelBar(ctx, 0, chunk1Width, barHeight, "#F2C894", r1.chrom);
    drawBaseline(ctx, 0, chunk1Width);

    // build ticks for chunk1
    const minMbC1 = Math.min(c1Min, c1Max);
    const maxMbC1 = Math.max(c1Min, c1Max);
    const t1Vals = buildTickVals(minMbC1, maxMbC1, tickStep, reversed);
    drawTicksBelow(ctx, scale1, t1Vals);

    // 2) The red vertical line at boundary
    ctx.save();
    ctx.strokeStyle = "red";
    ctx.beginPath();
    ctx.moveTo(chunk1Width, -barHeight);
    ctx.lineTo(chunk1Width, 0);
    ctx.stroke();
    ctx.restore();

    // 3) chunk2
    ctx.save();
    ctx.translate(chunk1Width, 0);
    drawPastelBar(ctx, 0, chunk2Width, barHeight, "#F2C894", r1.chrom);
    drawBaseline(ctx, 0, chunk2Width);

    // Ticks for chunk2
    const minMbC2 = Math.min(c2Min, c2Max);
    const maxMbC2 = Math.max(c2Min, c2Max);
    const t2Vals = buildTickVals(minMbC2, maxMbC2, tickStep, reversed);
    drawTicksBelow(ctx, scale2, t2Vals);

    ctx.restore();
    return;
  }

  /*****************************************************
   * CASE 3) Single-chrom, no break, no region2
   *****************************************************/
  console.log("custom x axis: single-chrom (no deletion), might be reversed");
  const reversed = !!r1.isReversed;
  const rawStart = r1.startMb;
  const rawEnd   = r1.endMb;
  const domain = reversed
    ? [rawEnd, rawStart]
    : [rawStart, rawEnd];
  const c1Mb = Math.abs(rawEnd - rawStart);
  const totalWidthPx = c1Mb * PX_PER_MB;

  const { ctx } = createCanvas(container, totalWidthPx, axisHeight, margin);

  // scale
  const scale = d3.scaleLinear()
    .domain(domain)
    .range([0, c1Mb * PX_PER_MB]);

  // bar
  drawPastelBar(ctx, 0, c1Mb * PX_PER_MB, barHeight, "#F2C894", r1.chrom);
  // baseline
  drawBaseline(ctx, 0, c1Mb * PX_PER_MB);

  // build ticks
  const minMb = Math.min(rawStart, rawEnd);
  const maxMb = Math.max(rawStart, rawEnd);
  const ticks = buildTickVals(minMb, maxMb, tickStep, reversed);
  drawTicksBelow(ctx, scale, ticks);
}
