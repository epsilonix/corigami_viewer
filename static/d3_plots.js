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
/**
 * Draws the Hi-C heatmap on a canvas, ensuring the x-axis aligns with the other plots.
 * The matrix's column indices are converted to genomic coordinates so that
 * tick marks (and their labels) match the x-axis of the line charts.
 */
/**
 * Draws the Hi-C heatmap on a canvas, ensuring the x-axis aligns with the other plots.
 * The matrix's column indices are converted to genomic coordinates so that
 * tick marks (and their labels) match the x-axis of the line charts.
 */
function drawHiCChart(containerSelector, config) {
  // Margins and dimensions.
  const margin = { top: 0, right: 20, bottom: 50, left: 50 };

  // The total genomic span in megabases.
  const genomicSpan = config.xAxis.max - config.xAxis.min;
  // The width in pixels (shared with other plots).
  const fullWidth = genomicSpan * PX_PER_MB;

  // Number of rows from the data (max row index + 1).
  const maxRow = d3.max(config.series[0].data, d => d[1]);
  const n_rows = maxRow + 1;
  // Number of columns is stored in config.n_cols (from prepare_plot_configs).
  const n_cols = config.n_cols;

  // Canvas height for the heatmap area.
  const height = config.chart.height || 250;
  const canvasWidth = fullWidth + margin.left + margin.right;
  const canvasHeight = height + margin.top + margin.bottom;

  // Create a high-resolution canvas and translate by margins.
  const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
  ctx.translate(margin.left, margin.top);

  // This scale converts matrix column indices (0..n_cols) to genomic coords in Mb.
  // Example: 0 -> config.xAxis.min, n_cols -> config.xAxis.max
  const colIndexToMb = d3.scaleLinear()
    .domain([0, n_cols])
    .range([config.xAxis.min, config.xAxis.max]);

  // This scale is used to map genomic coordinates (in Mb) to final x positions (pixels).
  // It's exactly the same as the line chart's x-scale for consistency.
  const xScale = d3.scaleLinear()
    .domain([config.xAxis.min, config.xAxis.max])
    .range([0, fullWidth]);

  // For rows, we keep a straightforward index-based scale:
  // domain = [0..n_rows], range = [0..height].
  const rowScale = d3.scaleLinear()
    .domain([0, n_rows])
    .range([0, height]);

  // If deletion mode is active, use the provided scale factor for x-stretching.
  // (Alternatively, you can remove or comment this out if you don't need horizontal stretching.)
  const sf = config.deletion ? config.scaleFactor : 1;

  // Define the color scale for cell values.
  const colorScale = d3.scaleLinear()
    .domain([config.colorAxis.min, config.colorAxis.max])
    .range([config.colorAxis.minColor, config.colorAxis.maxColor]);

  // Draw each cell in the heatmap.
  // d = [colIndex, rowIndex, heatmapValue]
  config.series[0].data.forEach(d => {
    const colIdx = d[0];
    const rowIdx = d[1];
    const value = d[2];

    // Convert column index to genomic MB, then to final pixel x.
    const mbLeft = colIndexToMb(colIdx);
    const mbRight = colIndexToMb(colIdx + 1); // next column over
    const xLeft = xScale(mbLeft) * sf;
    const xRight = xScale(mbRight) * sf;
    const w = xRight - xLeft; // cell width in px

    // Convert row index to y coordinate.
    const yTop = rowScale(rowIdx);
    const h = rowScale(rowIdx + 1) - yTop; // cell height

    ctx.fillStyle = colorScale(value);
    ctx.fillRect(xLeft, yTop, w, h);
  });

  // --- Draw the x-axis (same approach as in drawLineChart) ---
  ctx.save();
  ctx.translate(0, height); // move down to the bottom of the heatmap

  // Draw the baseline axis line.
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(fullWidth * sf, 0);
  ctx.strokeStyle = "black";
  ctx.stroke();

  // Generate ticks from the same domain as the line chart.
  // If you want them to compress/expand under deletion, you can apply * sf here.
  const xTicks = xScale.ticks(10);
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.font = "10px sans-serif";
  xTicks.forEach(tick => {
    // xScale(tick) is the pixel position for the coordinate "tick" in Mb
    const xPos = xScale(tick) * sf;

    // Draw the tick mark.
    ctx.beginPath();
    ctx.moveTo(xPos, 0);
    ctx.lineTo(xPos, 6);
    ctx.stroke();

    // Tick label: e.g. "5.00"
    ctx.fillText(d3.format(".2f")(tick), xPos, 6);
  });
  ctx.restore();

  // Optionally, draw an x-axis title if provided in config.
  ctx.textAlign = "center";
  ctx.textBaseline = "bottom";
  ctx.font = "12px sans-serif";
  if (config.xAxis.title) {
    ctx.fillText(config.xAxis.title, (fullWidth * sf) / 2, height + margin.bottom - 5);
  }
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
  if (config.xAxis.title) {
    ctx.fillText(config.xAxis.title, fullWidth / 2, height + margin.bottom - 5);
  }
  
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
  if (config.xAxis.title) {
    ctx.fillText(config.xAxis.title, fullWidth / 2, height + margin.bottom - 5);
  }
  
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
