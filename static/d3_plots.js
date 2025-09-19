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

// === Hatch helpers ============================================================
const HATCH_BAND_HEIGHT = 20; // px tall band

function hexToRgb(hex) {
  const m = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  if (!m) return { r: 159, g: 192, b: 222 }; // fallback (#9FC0DE)
  return { r: parseInt(m[1], 16), g: parseInt(m[2], 16), b: parseInt(m[3], 16) };
}

function darkenHex(hex, amount = 30) { // amount: 0..255 darker
  const { r, g, b } = hexToRgb(hex);
  const clamp = (v) => Math.max(0, Math.min(255, v));
  const d = (v) => clamp(v - amount);
  const toHex = (v) => v.toString(16).padStart(2, "0");
  return `#${toHex(d(r))}${toHex(d(g))}${toHex(d(b))}`;
}

function makeDiagHatch(colorHex, bgAlpha = 0.12, lineAlpha = 0.5, size = 8, thickness = 1) {
  const off = document.createElement("canvas");
  off.width = off.height = size;
  const c = off.getContext("2d");
  const { r, g, b } = hexToRgb(colorHex);

  // subtle tinted background
  c.fillStyle = `rgba(${r},${g},${b},${bgAlpha})`;
  c.fillRect(0, 0, size, size);

  // diagonal strokes (same angle/feel as your disabled buttons)
  c.strokeStyle = `rgba(${r},${g},${b},${lineAlpha})`;
  c.lineWidth = thickness;
  c.beginPath();
  c.moveTo(-size, size / 2);
  c.lineTo(size / 2, -size);
  c.moveTo(0, size);
  c.lineTo(size, 0);
  c.moveTo(size / 2, size);
  c.lineTo(size * 1.5, size / 2);
  c.stroke();

  return c.createPattern(off, "repeat");
}

function drawHatchedBand(ctx, xStartPx, xEndPx, bandHeightPx, colorHex) {
  const w = Math.max(0, xEndPx - xStartPx);
  if (!w) return;
  const pattern = makeDiagHatch(colorHex);
  ctx.save();
  ctx.fillStyle = pattern;
  ctx.fillRect(xStartPx, 0, w, bandHeightPx);
  ctx.strokeStyle = darkenHex(colorHex, 40); // darker outline for definition
  ctx.lineWidth = 1;
  ctx.strokeRect(Math.floor(xStartPx) + 0.5, 0.5, Math.floor(w) - 1, Math.max(1, bandHeightPx) - 1);
  ctx.restore();
}

/**
 * Draws a column (bar) chart using the given config.
 * - If config.xAxis.axisBreak is not true, it does the standard single-domain approach.
 * - If config.xAxis.axisBreak === true, and leftDomain/rightDomain are provided,
 *   it uses a two-scale approach to skip the middle region (like a real axis break).
 *
 * The rest of the code is your original logic for binning, screening shading, Y-axis, etc.
 */
// === Continuous-series helpers ==============================================

function resampleUniformXY(data, xMin, xMax, stepMb) {
  // data: [[xMb, y], ...], assumed sorted by x
  if (!data.length) return [];
  const out = [];
  let i = 0;

  // clamp bounds to data
  const dMin = data[0][0], dMax = data[data.length - 1][0];
  const X0 = Math.max(xMin, dMin);
  const X1 = Math.min(xMax, dMax);
  if (X1 <= X0) return out;

  for (let x = X0; x <= X1 + 1e-9; x += stepMb) {
    while (i < data.length - 1 && data[i + 1][0] < x) i++;
    // linear interpolate between data[i] and data[i+1]
    if (i >= data.length - 1) {
      out.push([x, data[data.length - 1][1]]);
    } else {
      const [x0, y0] = data[i];
      const [x1, y1] = data[i + 1];
      const t = (x - x0) / Math.max(1e-9, (x1 - x0));
      out.push([x, y0 + t * (y1 - y0)]);
    }
  }
  return out;
}

function movingAverage(vals, winPoints) {
  if (winPoints <= 1) return vals.slice();
  const half = Math.floor(winPoints / 2);
  const out = new Array(vals.length);
  let sum = 0;
  for (let i = 0; i < vals.length; i++) {
    sum += vals[i];
    if (i > winPoints) sum -= vals[i - winPoints - 1];
    const left = Math.max(0, i - half);
    const right = Math.min(vals.length - 1, i + half);
    const len = right - left + 1;
    let s = 0;
    for (let j = left; j <= right; j++) s += vals[j];
    out[i] = s / len;
  }
  return out;
}

function gaussianKernel(winPoints, sigmaPoints) {
  const half = Math.floor(winPoints / 2);
  const k = new Array(winPoints);
  const s2 = 2 * sigmaPoints * sigmaPoints;
  let sum = 0;
  for (let i = -half; i <= half; i++) {
    const v = Math.exp(-(i * i) / s2);
    k[i + half] = v;
    sum += v;
  }
  for (let i = 0; i < k.length; i++) k[i] /= sum;
  return k;
}

function convolve1D(vals, kernel) {
  const n = vals.length, m = kernel.length, half = Math.floor(m / 2);
  const out = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let acc = 0, wsum = 0;
    for (let j = 0; j < m; j++) {
      const idx = i + j - half;
      if (idx >= 0 && idx < n) {
        acc += vals[idx] * kernel[j];
        wsum += kernel[j];
      }
    }
    out[i] = wsum > 0 ? acc / wsum : vals[i];
  }
  return out;
}

function smoothSeriesXY(xy, smoothingCfg, pxPerMb = PX_PER_MB) {
  // smoothingCfg: { method: 'gaussian'|'movingAvg', windowMb: 0.05, sigmaMb?: 0.02, resamplePx?: 1.5 }
  if (!xy.length || !smoothingCfg) return xy;

  const xMin = xy[0][0], xMax = xy[xy.length - 1][0];
  const resamplePx = Math.max(0.5, smoothingCfg.resamplePx ?? 1.5);  // ~1-2 px steps
  const stepMb = resamplePx / pxPerMb;

  // Resample to uniform grid
  const uniformXY = resampleUniformXY(xy, xMin, xMax, stepMb);
  if (!uniformXY.length) return xy;

  const xs = uniformXY.map(d => d[0]);
  let ys = uniformXY.map(d => d[1]);

  // Smoothing
  const windowMb = Math.max(0, smoothingCfg.windowMb ?? 0);
  if (windowMb > 0) {
    const winPoints = Math.max(1, Math.round(windowMb / stepMb));
    if (smoothingCfg.method === 'gaussian') {
      const sigmaMb = Math.max(1e-6, smoothingCfg.sigmaMb ?? (windowMb / 3));
      const sigmaPoints = Math.max(1, Math.round(sigmaMb / stepMb));
      const kernel = gaussianKernel(winPoints | 1, sigmaPoints);
      ys = convolve1D(ys, kernel);
    } else {
      ys = movingAverage(ys, winPoints | 1);
    }
  }

  return xs.map((x, i) => [x, ys[i]]);
}

// replace your current drawLineAndArea with this drop-in
function drawLineAndArea(ctx, pts, xScale, yScale, axisPx, color, doArea, stroke = true, areaAlpha = 0.28) {
  if (!pts.length) return;

  ctx.save();
  ctx.lineWidth = 1.5;
  ctx.lineJoin  = "round";
  ctx.lineCap   = "round";

  const alpha = Number.isFinite(areaAlpha) ? Math.max(0, Math.min(1, areaAlpha)) : 0.28;

  if (doArea) {
    const areaGen = d3.area()
      .x(d => xScale(d[0]))
      .y0(axisPx)
      .y1(d => yScale(d[1]))
      .curve(d3.curveMonotoneX)
      .context(ctx);

    // Convert hex color to RGBA with the specified alpha
    const rgb = hexToRgb(color);
    ctx.fillStyle = `rgba(${rgb.r}, ${rgb.g}, ${rgb.b}, ${alpha})`;
    
    ctx.beginPath();
    areaGen(pts);
    ctx.fill();
  }

  // Draw the stroke if requested (after the area so it's on top)
  if (stroke) {
    const lineGen = d3.line()
      .x(d => xScale(d[0]))
      .y(d => yScale(d[1]))
      .curve(d3.curveMonotoneX)
      .context(ctx);

    ctx.beginPath();
    lineGen(pts);
    ctx.strokeStyle = color;
    ctx.stroke();
  }

  ctx.restore();
}





function drawColumnChart(containerSelector, config) {
  console.log("=== drawColumnChart START ===");
  console.log("config.screening:", config.screening);
  console.log("config.chart.type:", config.chart?.type);

  const PX_PER_MB = 200;
  const MIN_BARS_PER_MB = 30;
  const margin = { top: 20, right: 20, bottom: 40, left: 50 };

  // 1) Validate x-range
  if (!config.xAxis || typeof config.xAxis.min !== "number" || typeof config.xAxis.max !== "number") {
    console.warn("Missing/invalid xAxis.min/max:", config.xAxis);
    return;
  }
  const xMin = config.xAxis.min;
  const xMax = config.xAxis.max;
  const domainSpanMb = xMax - xMin;

  // FIXED: Screening should always be bars, never continuous
  let chartType, isContinuous;
  if (config.screening === true || config.chart?.type === "column") {
    chartType = "column";
    isContinuous = false;
  } else {
    chartType = (config.chart && config.chart.type) || "area";
    isContinuous = (chartType === "line" || chartType === "area" || config.forceContinuous === true);
  }
  console.log("chartType:", chartType);
  console.log("isContinuous:", isContinuous);
  console.log("config.screening (again):", config.screening);
  
  const stroke = (config.stroke !== false);   // default true; false hides the line
  const areaAlpha = config.areaAlpha ?? 0.28; // optional fill opacity
  const hasBreak   = !!config.xAxis.axisBreak;
  const leftDomain = config.xAxis.leftDomain || [];
  const rightDomain= config.xAxis.rightDomain || [];

  // Rest of the function continues as before...
  // (The key change is that isContinuous is now forced to false when config.screening is true)

  // 2) Data prep (keep your binning for bar mode; continuous uses raw for fidelity)
  let originalData = (config.series && config.series[0] && config.series[0].data) || [];
  originalData.sort((a, b) => a[0] - b[0]); // sort by x

  let data = [];
  if (config.screening) {
    data = originalData;
  } else if (!isContinuous) {
    // bars: ensure at least min bar density
    const minBars = Math.round(domainSpanMb * MIN_BARS_PER_MB);
    data = (originalData.length >= minBars) ? originalData
                                            : binAndUpsampleData(originalData, xMin, xMax, minBars);
  } else {
    // continuous: use raw points (we'll resample uniformly later for smoothing)
    data = originalData;
  }

  // 3) Y-domain
  const rawYMin = d3.min(data, d => d[1]);
  const rawYMax = d3.max(data, d => d[1]);
  let yMin = (typeof rawYMin === "number") ? rawYMin : 0;
  let yMax = (typeof rawYMax === "number") ? rawYMax : 0;

  let domainMin, domainMax, axisValue;
  if (yMax < 0) { domainMin = yMin; domainMax = yMax; axisValue = yMin; }
  else if (yMin > 0) { domainMin = 0; domainMax = yMax; axisValue = 0; }
  else { domainMin = yMin; domainMax = yMax; axisValue = 0; }

  const height = config.chart.height || 150;
  const color  = (config.series && config.series[0] && config.series[0].color) || "#9FC0DE";
  const smoothing = config.screening ? null : (config.smoothing || null);


  // ===== Axis-break two-scale =====
  if (hasBreak && leftDomain.length === 2 && rightDomain.length === 2) {
    const leftSpan  = leftDomain[1]  - leftDomain[0];
    const rightSpan = rightDomain[1] - rightDomain[0];
    const leftWidthPx  = leftSpan  * PX_PER_MB;
    const rightWidthPx = rightSpan * PX_PER_MB;
    const totalWidthPx = leftWidthPx + rightWidthPx;

    const canvasWidth  = totalWidthPx + margin.left + margin.right;
    const canvasHeight = height       + margin.top  + margin.bottom;
    const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
    ctx.translate(margin.left, margin.top);

    const scaleLeft  = d3.scaleLinear().domain([leftDomain[0],  leftDomain[1]]).range([0, leftWidthPx]);
    const scaleRight = d3.scaleLinear().domain([rightDomain[0], rightDomain[1]]).range([0, rightWidthPx]);

    const leftData  = data.filter(d => d[0] >= leftDomain[0]  && d[0] <= leftDomain[1]);
    const rightData = data.filter(d => d[0] >= rightDomain[0] && d[0] <= rightDomain[1]);

    const yScale = d3.scaleLinear().domain([domainMin, domainMax]).nice().range([height, 0]);
    const axisPx = yScale(axisValue);
    const smoothing = config.smoothing || { method: 'gaussian', windowMb: 0.06, sigmaMb: 0.02, resamplePx: 0.75 };
    const stroke = (config.stroke !== false);        // default true; set to false to hide stroke
    const areaAlpha = config.areaAlpha ?? 0.28;    
    if (isContinuous) {
      // Left chunk
      const leftSmooth = smoothSeriesXY(leftData, smoothing, PX_PER_MB);
      drawLineAndArea(ctx, leftSmooth, scaleLeft, yScale, axisPx, color, chartType === "area", stroke, areaAlpha);

      // Baseline under left
      ctx.beginPath(); ctx.moveTo(0, axisPx); ctx.lineTo(leftWidthPx, axisPx); ctx.strokeStyle = "#000"; ctx.stroke();

      // Right chunk
      ctx.save();
      ctx.translate(leftWidthPx, 0);
      const rightSmooth = smoothSeriesXY(rightData, smoothing, PX_PER_MB);
      drawLineAndArea(ctx, rightSmooth, scaleRight, yScale, axisPx, color, chartType === "area", stroke, areaAlpha);
      ctx.beginPath(); ctx.moveTo(0, axisPx); ctx.lineTo(rightWidthPx, axisPx); ctx.strokeStyle = "#000"; ctx.stroke();
      ctx.restore();
    } else {
      // === original bars (unchanged) ===
      let leftCount = leftData.length;
      let leftSpacing = (leftCount > 1) ? (leftWidthPx / (leftCount - 1)) : leftWidthPx;
      let barMultiplier = config.screening ? 0.02 : 0.05;
      let leftBarWidth = Math.min(Math.max(leftSpacing * barMultiplier, 1), 20);

      leftData.forEach(d => {
        const xPx = scaleLeft(d[0]) - leftBarWidth/2;
        const yValPx = yScale(d[1]);
        const top = Math.min(yValPx, axisPx), bottom = Math.max(yValPx, axisPx);
        ctx.fillStyle = color;
        ctx.fillRect(xPx, top, leftBarWidth, bottom - top);
      });
      ctx.beginPath(); ctx.moveTo(0, axisPx); ctx.lineTo(leftWidthPx, axisPx); ctx.strokeStyle = "#000"; ctx.stroke();

      ctx.save();
      ctx.translate(leftWidthPx, 0);
      let rightCount = rightData.length;
      let rightSpacing= (rightCount > 1) ? (rightWidthPx / (rightCount - 1)) : rightWidthPx;
      let rightBarWidth = Math.min(Math.max(rightSpacing * barMultiplier, 1), 20);
      rightData.forEach(d => {
        const xPx = scaleRight(d[0]) - rightBarWidth/2;
        const yValPx = yScale(d[1]);
        const top = Math.min(yValPx, axisPx), bottom = Math.max(yValPx, axisPx);
        ctx.fillStyle = color;
        ctx.fillRect(xPx, top, rightBarWidth, bottom - top);
      });
      ctx.beginPath(); ctx.moveTo(0, axisPx); ctx.lineTo(rightWidthPx, axisPx); ctx.strokeStyle = "#000"; ctx.stroke();
      ctx.restore();
    }

    // Y axis
    ctx.beginPath(); ctx.moveTo(0, 0); ctx.lineTo(0, height); ctx.strokeStyle = "#000"; ctx.stroke();
    ctx.save(); ctx.textAlign = "right"; ctx.textBaseline = "middle"; ctx.font = "10px sans-serif";
    yScale.ticks(5).forEach(t => { const y = yScale(t); ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(-5, y); ctx.stroke(); ctx.fillText(t, -7, y); });
    ctx.restore();

    console.log("=== drawColumnChart END (two-scale) ===");
    return;
  }

  // ===== Single-scale =====
  const fullWidth = domainSpanMb * PX_PER_MB;
  const canvasWidth  = fullWidth + margin.left + margin.right;
  const canvasHeight = height    + margin.top  + margin.bottom;
  const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
  ctx.translate(margin.left, margin.top);

  const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, fullWidth]);
  const yScale = d3.scaleLinear().domain([domainMin, domainMax]).nice().range([height, 0]);
  const axisPx = yScale(axisValue);
  if (isContinuous) {
    const smoothXY = smoothSeriesXY(data, smoothing, PX_PER_MB);
    drawLineAndArea(ctx, smoothXY, xScale, yScale, axisPx, color, chartType === "area", stroke, areaAlpha);
  } else {
    // original bars
    let count = data.length;
    let spacing = (count > 1) ? (fullWidth / (count - 1)) : fullWidth;
    let barMultiplier = config.screening ? 0.02 : 0.05;
    let barWidth = Math.min(Math.max(spacing * barMultiplier, 1), 20);
    data.forEach(d => {
      const xPx = xScale(d[0]) - barWidth/2;
      const yValPx = yScale(d[1]);
      const top = Math.min(yValPx, axisPx), bottom = Math.max(yValPx, axisPx);
      ctx.fillStyle = color;
      ctx.fillRect(xPx, top, barWidth, bottom - top);
    });
  }

  // Baseline
  ctx.beginPath(); ctx.moveTo(0, axisPx); ctx.lineTo(fullWidth, axisPx); ctx.strokeStyle = "#000"; ctx.stroke();

  // Y axis
  ctx.beginPath(); ctx.moveTo(0, 0); ctx.lineTo(0, height); ctx.strokeStyle = "#000"; ctx.stroke();
  ctx.save(); ctx.textAlign = "right"; ctx.textBaseline = "middle"; ctx.font = "10px sans-serif";
  yScale.ticks(5).forEach(t => { const y = yScale(t); ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(-5, y); ctx.stroke(); ctx.fillText(t, -7, y); });
  ctx.restore();

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

function drawTicksBelowStacked(
  ctx, scale, tickVals,
  {
    xOffset = 0,        // how much the context is translated for drawing
    worldOffset = 0,    // absolute X offset for collision checks (e.g., chunk1Width)
    labelFormatter = t => t.toFixed(1),
    minGapPx = 6,       // padding between labels
    rowHeight = 12,     // vertical step for stacked rows
    baselineY = 0,
    maxRows = 2,
    font = "10px sans-serif"
  } = {},
  rowRightmost = []     // carry-over state across segments
) {
  ctx.save();
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.font = font;

  // ensure rows
  for (let r = rowRightmost.length; r < maxRows; r++) rowRightmost[r] = -Infinity;

  // sort by on-screen x to get consistent collision behavior with reversed scales
  const ticksSorted = tickVals.slice().sort((a, b) => scale(a) - scale(b));

  // draw ticks + place labels with stacking
  const labels = [];
  ticksSorted.forEach(t => {
    const xLocal = scale(t);
    const xDraw = xLocal + xOffset;          // where we draw on the current ctx
    const xWorld = xLocal + worldOffset;     // where we *measure* collisions globally

    // tick mark
    ctx.beginPath();
    ctx.moveTo(xDraw, baselineY);
    ctx.lineTo(xDraw, baselineY + 6);
    ctx.strokeStyle = "#000";
    ctx.stroke();

    const txt = labelFormatter(t);
    const w = ctx.measureText(txt).width;

    // choose a row that doesn't collide
    let row = 0;
    while (row < maxRows && (xWorld - w / 2) < (rowRightmost[row] + minGapPx)) row++;
    if (row >= maxRows) row = maxRows - 1; // clamp to last row if all crowded

    rowRightmost[row] = Math.max(rowRightmost[row], xWorld + w / 2);
    labels.push({ x: xDraw, y: baselineY + 8 + row * rowHeight, txt });
  });

  // draw labels after placing all
  ctx.fillStyle = "#000";
  labels.forEach(({ x, y, txt }) => ctx.fillText(txt, x, y));

  ctx.restore();
  return rowRightmost;
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
    let rowState = [];
    rowState = drawTicksBelowStacked(ctx, scale1, t1Vals, {
    xOffset: 0,             // we didn't translate yet
    worldOffset: 0          // absolute 0..chunk1Width
    }, rowState);

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
    rowState = drawTicksBelowStacked(ctx, scale2, t2Vals, {
    xOffset: 0,                 // drawing inside the translated ctx
    worldOffset: chunk1Width    // but collisions are checked in world coords
    }, rowState);

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
