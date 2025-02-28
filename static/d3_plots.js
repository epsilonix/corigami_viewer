/* d3_plots.js */

/**
 * Global conversion: pixels per megabase.
 * Adjust this value so that, for example, 1 Mb = 100 pixels.
 */
const PX_PER_MB = 300;

/**
 * Helper function to create a high-resolution canvas.
 * It sets the canvas width/height to the desired size multiplied by devicePixelRatio,
 * then scales the context so that drawing commands use CSS pixels.
 */
function createHiResCanvas(containerSelector, width, height) {
  // Remove any existing elements.
  d3.select(containerSelector).selectAll("*").remove();
  
  // Append canvas.
  const canvasSelection = d3.select(containerSelector)
    .append("canvas")
    .style("width", width + "px")
    .style("height", height + "px");
  
  // Get the raw DOM node.
  const canvas = canvasSelection.node();
  
  // Get device pixel ratio.
  const ratio = window.devicePixelRatio || 1;
  
  // Set the actual canvas size.
  canvas.width = width * ratio;
  canvas.height = height * ratio;
  
  // Return the canvas and the scaled context.
  const ctx = canvas.getContext("2d");
  ctx.scale(ratio, ratio);
  return { canvas, ctx };
}

/**
 * Draws the Hi-C heatmap on a canvas.
 * In deletion mode (if config.deletion exists), the hi-c matrix is stretched horizontally
 * using config.scaleFactor and the x-axis ticks are drawn.
 */
function drawHiCChart(containerSelector, config) {
  // Margins and dimensions.
  const margin = { top: 0, right: 20, bottom: 50, left: 50 };
  const fullWidth = (config.xAxis.max - config.xAxis.min) * PX_PER_MB;
  
  // Compute number of rows from data.
  const maxRow = d3.max(config.series[0].data, d => d[1]);
  const n_rows = maxRow + 1;
  const height = config.chart.height || 250;
  const cellHeight = height / n_rows; // each row gets an equal share of the drawing area
  
  // Total canvas dimensions.
  const canvasWidth = fullWidth + margin.left + margin.right;
  const canvasHeight = height + margin.top + margin.bottom;
  
  // Create a high-resolution canvas.
  const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
  
  // Translate the context by the margins.
  ctx.translate(margin.left, margin.top);
  
  // Define scales.
  const colScale = d3.scaleLinear()
    .domain([0, config.n_cols])
    .range([0, fullWidth]);
  const rowScale = d3.scaleLinear()
    .domain([0, n_rows])
    .range([0, height]);
  
  // Use scaleFactor if deletion mode is active.
  const sf = config.deletion ? config.scaleFactor : 1;
  
  // Color scale.
  const colorScale = d3.scaleLinear()
    .domain([config.colorAxis.min, config.colorAxis.max])
    .range([config.colorAxis.minColor, config.colorAxis.maxColor]);
  
  // Draw each cell.
  config.series[0].data.forEach(d => {
    const x = colScale(d[0]) * sf;
    const y = rowScale(d[1]);
    const w = (colScale(1) - colScale(0)) * sf;
    const h = rowScale(1) - rowScale(0);
    ctx.fillStyle = colorScale(d[2]);
    ctx.fillRect(x, y, w, h);
  });
  
  // --- Draw the x-axis ---
  // Create an x-scale for axis ticks.
  const xScale = d3.scaleLinear()
    .domain([config.xAxis.min, config.xAxis.max])
    .range([0, fullWidth]);
  
  // Save context and move to the bottom of the drawing area.
  ctx.save();
  ctx.translate(0, height);
  
  // Draw axis line.
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(fullWidth, 0);
  ctx.strokeStyle = "black";
  ctx.stroke();
  
  // Generate ticks and draw them.
  const ticks = xScale.ticks(10);
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.font = "10px sans-serif";
  ticks.forEach(tick => {
    const xPos = xScale(tick);
    // Draw tick mark.
    ctx.beginPath();
    ctx.moveTo(xPos, 0);
    ctx.lineTo(xPos, 6);
    ctx.stroke();
    // Draw tick label.
    ctx.fillText(d3.format(".2f")(tick), xPos, 6);
  });
  ctx.restore();
  
  // Draw x-axis title.
  ctx.textAlign = "center";
  ctx.textBaseline = "bottom";
  ctx.font = "12px sans-serif";
  ctx.fillText(config.xAxis.title, fullWidth / 2, height + margin.bottom - 5);
}

/**
 * Draws a line chart for signal plots (CTCF or ATAC) on a canvas.
 * Expects config.series[0].data to be an array of [x, y] pairs.
 */
function drawLineChart(containerSelector, config) {
  console.log('Line chart config:', config);
  
  const margin = { top: 20, right: 20, bottom: 40, left: 50 };
  const fullWidth = (config.xAxis.max - config.xAxis.min) * PX_PER_MB;
  const height = config.chart.height || 150;
  
  const canvasWidth = fullWidth + margin.left + margin.right;
  const canvasHeight = height + margin.top + margin.bottom;
  
  // Create high-resolution canvas.
  const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
  ctx.translate(margin.left, margin.top);
  
  // Create scales.
  const xScale = d3.scaleLinear()
    .domain([config.xAxis.min, config.xAxis.max])
    .range([0, fullWidth]);
  const yExtent = d3.extent(config.series[0].data, d => d[1]);
  const yScale = d3.scaleLinear()
    .domain(yExtent)
    .range([height, 0]);
  
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
  
  const yTicks = yScale.ticks(5);
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
  
  // --- Draw the line ---
  ctx.beginPath();
  config.series[0].data.forEach((d, i) => {
    const x = xScale(d[0]);
    const y = yScale(d[1]);
    if (i === 0) {
      ctx.moveTo(x, y);
    } else {
      ctx.lineTo(x, y);
    }
  });
  ctx.strokeStyle = config.series[0].color || "steelblue";
  ctx.lineWidth = 1.5;
  ctx.stroke();
  
  // x-axis title.
  ctx.textAlign = "center";
  ctx.textBaseline = "bottom";
  ctx.font = "12px sans-serif";
  ctx.fillText(config.xAxis.title, fullWidth / 2, height + margin.bottom - 5);
  
  // Chart title if provided.
  if (config.title && config.title.text) {
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.font = (config.title.style && config.title.style.fontSize ? config.title.style.fontSize : "13px") + " sans-serif";
    ctx.fillText(config.title.text, fullWidth / 2, -margin.top / 2);
  }
}

/**
 * Draws a column (bar) chart for screening on a canvas.
 * Uses a linear scale for the x-axis to sync with other plots.
 * Adjusts the bar width based on the spacing between data points.
 */
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
