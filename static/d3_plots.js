/*************************************************************
 * Canvas-Based Plotting for C.Origami
 * (Supports Hi-C heatmap, signal plots with/without x-axis break,
 *  and screening column charts.)
 *************************************************************/

/**
 * Pixels per Mb in normal (continuous) mode.
 */
const PX_PER_MB = 300;

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
  // In deletion mode, we simply hide the x-axis ticks.
  const hideTicks = config.hideXAxis === true;
  
  const margin = { top: 0, right: 20, bottom: 50, left: 50 };
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
  const colorScale = d3.scaleLinear()
    .domain([config.colorAxis.min, config.colorAxis.max])
    .range([config.colorAxis.minColor, config.colorAxis.maxColor]);
  
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
  
  if (!hideTicks) {
    ctx.save();
    ctx.translate(0, height);
    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(fullWidth, 0);
    ctx.strokeStyle = "#000";
    ctx.stroke();
    
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
    ctx.restore();
    
    if (config.xAxis.title) {
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


/*************************************************************
 * 2) Signal Plots (Canvas) – ATAC/CTCF
 *    Supports continuous mode and deletion-mode (with axis break)
 *************************************************************/
function drawLineChart(containerSelector, config) {
  const xAxisCfg = config.xAxis || {};
  if (xAxisCfg.axisBreak && xAxisCfg.leftDomain && xAxisCfg.rightDomain) {
    drawLineChartWithBreak(containerSelector, config);
  } else {
    drawLineChartContinuous(containerSelector, config);
  }
}

/** Continuous (normal) mode for signal plots. */
function drawLineChartContinuous(containerSelector, config) {
  const margin = { top: 20, right: 20, bottom: 40, left: 50 };
  const xMin = config.xAxis.min, xMax = config.xAxis.max;
  const fullWidth = (xMax - xMin) * PX_PER_MB;
  const height = config.chart.height || 150;
  
  const canvasWidth  = fullWidth + margin.left + margin.right;
  const canvasHeight = height + margin.top + margin.bottom;
  
  const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
  ctx.translate(margin.left, margin.top);
  
  const data = (config.series[0] && config.series[0].data) ? config.series[0].data : [];
  
  const xScale = d3.scaleLinear()
    .domain([xMin, xMax])
    .range([0, fullWidth]);
  
  const yExtent = d3.extent(data, d => d[1]);
  const yScale = d3.scaleLinear()
    .domain(yExtent)
    .nice()
    .range([height, 0]);
  
  // Draw x-axis
  ctx.save();
  ctx.translate(0, height);
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(fullWidth, 0);
  ctx.strokeStyle = "#000";
  ctx.stroke();
  
  ctx.font = "10px sans-serif";
  ctx.fillStyle = "#000";
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  xScale.ticks(10).forEach(tick => {
    const xPos = xScale(tick);
    ctx.beginPath();
    ctx.moveTo(xPos, 0);
    ctx.lineTo(xPos, 6);
    ctx.stroke();
    ctx.fillText(d3.format(".2f")(tick), xPos, 6);
  });
  ctx.restore();
  
  // Draw y-axis (ticks only; no label)
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(0, height);
  ctx.strokeStyle = "#000";
  ctx.stroke();
  
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  ctx.font = "10px sans-serif";
  yScale.ticks(5).forEach(tick => {
    const yPos = yScale(tick);
    ctx.beginPath();
    ctx.moveTo(0, yPos);
    ctx.lineTo(-6, yPos);
    ctx.stroke();
    ctx.fillText(tick, -8, yPos);
  });
  
  // Draw the signal line
  ctx.beginPath();
  data.forEach((d, i) => {
    const x = xScale(d[0]);
    const y = yScale(d[1]);
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  });
  ctx.strokeStyle = config.series[0].color || "steelblue";
  ctx.lineWidth = 1.5;
  ctx.stroke();
  
  // Draw x-axis label if provided
  if (config.xAxis.title) {
    ctx.save();
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.font = "12px sans-serif";
    ctx.fillStyle = "#000";
    ctx.fillText(config.xAxis.title, fullWidth / 2, height + margin.bottom - 5);
    ctx.restore();
  }
}

/** Deletion-mode signal plot with discontinuous (broken) x-axis. */
function drawColumnChart(containerSelector, config) {
  const margin = { top: 20, right: 20, bottom: 40, left: 50 };
  const fullWidth = (config.xAxis.max - config.xAxis.min) * PX_PER_MB;
  const height = config.chart.height || 150;
  
  const canvasWidth = fullWidth + margin.left + margin.right;
  const canvasHeight = height + margin.top + margin.bottom;
  
  // Create high-resolution canvas.
  const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
  ctx.translate(margin.left, margin.top);
  
  const series = config.series && config.series[0];
  if (!series) {
    console.error("No series data found in config for column chart");
    return;
  }
  let data = series.data;
  if (data.length === 0) {
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText("No data available", fullWidth / 2, height / 2);
    return;
  }
  
  // Sort data by x value.
  data.sort((a, b) => a[0] - b[0]);
  
  const xMin = config.xAxis.min;
  const xMax = config.xAxis.max;
  const xScale = d3.scaleLinear()
    .domain([xMin, xMax])
    .range([0, fullWidth]);
  
  // Compute spacing and bar width.
  const spacingPixels = fullWidth / (data.length - 1);
  const barWidth = spacingPixels * 0.02;
  const yScale = d3.scaleLinear()
    .domain([0, d3.max(data, d => d[1])])
    .nice()
    .range([height, 0]);
  
  // Draw bars.
  data.forEach(d => {
    const x = xScale(d[0]) - barWidth / 2;
    const y = yScale(d[1]);
    const h = height - y;
    ctx.fillStyle = series.color || "steelblue";
    ctx.fillRect(x, y, barWidth, h);
  });
  
  // --- Draw x-axis ---
  ctx.save();
  ctx.translate(0, height);
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(fullWidth, 0);
  ctx.strokeStyle = "black";
  ctx.stroke();
  
  const xTicks = xScale.ticks(10);
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.font = "10px sans-serif";
  xTicks.forEach(tick => {
    const xPos = xScale(tick);
    ctx.beginPath();
    ctx.moveTo(xPos, 0);
    ctx.lineTo(xPos, 6);
    ctx.stroke();
    ctx.fillText(d3.format(".2f")(tick), xPos, 6);
  });
  ctx.restore();
  
  // --- Draw y-axis ---
  ctx.save();
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(0, height);
  ctx.strokeStyle = "black";
  ctx.stroke();
  
  const yTicks = yScale.ticks();
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  ctx.font = "10px sans-serif";
  yTicks.forEach(tick => {
    const yPos = yScale(tick);
    ctx.beginPath();
    ctx.moveTo(0, yPos);
    ctx.lineTo(-6, yPos);
    ctx.stroke();
    ctx.fillText(tick, -8, yPos);
  });
  ctx.restore();
  
  // x-axis title.
  ctx.textAlign = "center";
  ctx.textBaseline = "bottom";
  ctx.font = "12px sans-serif";
  ctx.fillText(config.xAxis.title.text, fullWidth / 2, height + margin.bottom - 10);
  
  // y-axis title.
  ctx.save();
  ctx.translate(-margin.left, height / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.font = "12px sans-serif";
  ctx.fillText(config.yAxis.title.text, 0, 0);
  ctx.restore();
  
  // Chart title if provided.
  if (config.title && config.title.text) {
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.font = "13px sans-serif";
    ctx.fillText(config.title.text, fullWidth / 2, -margin.top / 2);
  }
}