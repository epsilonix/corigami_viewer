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

const BUFFER_BP = 5000; // used for row assignment

function drawGeneTrackChart(selector, config) {
  console.log("Drawing gene track with config:", config);
  
  // Compute region span in Mb (using region.start/end in base pairs)
  const regionStart = config.region.start;
  const regionEnd = config.region.end;
  const regionSpanMb = (regionEnd - regionStart) / 1e6;
  
  // Determine if we are in deletion mode by checking for deletionStart/End.
  const inDeletionMode = config.deletionStart != null && config.deletionEnd != null;
  
  let fullWidth;
  let xScale;
  
  if (config.xAxis && config.xAxis.axisBreak && inDeletionMode) {
    console.log("gene track: in deletion mode");
    // Deletion mode: use discontinuous x-scale.
    // Define left and right segments in base pairs.
    const leftDomainBp = [regionStart, config.deletionStart];
    const rightDomainBp = [config.deletionEnd, regionEnd];
    
    console.log("leftDomainBp:", leftDomainBp);
    console.log("rightDomainBp:", rightDomainBp);
    
    const leftSpanMb = (config.deletionStart - regionStart) / 1e6;
    const rightSpanMb = (regionEnd - config.deletionEnd) / 1e6;
    const totalSpanMb = leftSpanMb + rightSpanMb;
    
    // Full width in pixels is computed from the total displayed span.
    fullWidth = totalSpanMb * PX_PER_MB;
    
    // For gene track we want the two segments to be right next to each other.
    const gapPx = 0;
    
    // Compute pixel widths for left and right segments.
    const leftWidthPx = (leftSpanMb / totalSpanMb) * fullWidth;
    const rightWidthPx = (rightSpanMb / totalSpanMb) * fullWidth;
    
    // Create discontinuous scales.
    const scaleLeft = d3.scaleLinear()
                        .domain(leftDomainBp)
                        .range([0, leftWidthPx]);
    const scaleRight = d3.scaleLinear()
                         .domain(rightDomainBp)
                         .range([0, rightWidthPx]);
    
    xScale = function(xVal) {
      if (xVal <= config.deletionStart) {
        return scaleLeft(xVal);
      } else if (xVal >= config.deletionEnd) {
        return leftWidthPx + gapPx + scaleRight(xVal);
      }
      // If within deletion, return NaN.
      return NaN;
    };
  } else {
    // Non-deletion mode: use a continuous scale.
    fullWidth = regionSpanMb * PX_PER_MB;
    xScale = d3.scaleLinear()
               .domain([regionStart, regionEnd])
               .range([0, fullWidth]);
  }
  
  // Define margins and overall SVG dimensions.
  const margin = { top: 5, right: 5, bottom: 20, left: 50 };
  const height = config.chart.height || 50;
  const canvasWidth = fullWidth + margin.left + margin.right;
  const canvasHeight = height + margin.top + margin.bottom;
  
  // Create the SVG container.
  const svg = d3.select(selector)
                .append("svg")
                .attr("width", canvasWidth)
                .attr("height", canvasHeight);
  
  const g = svg.append("g")
               .attr("transform", `translate(${margin.left},${margin.top})`);
  
  // Define arrow markers.
  const defs = g.append("defs");
  defs.append("marker")
      .attr("id", "arrowhead-right")
      .attr("viewBox", "0 -5 10 10")
      .attr("refX", 4)
      .attr("refY", 0)
      .attr("markerWidth", 4)
      .attr("markerHeight", 4)
      .attr("orient", "auto")
    .append("path")
      .attr("d", "M0,-5L10,0L0,5")
      .attr("fill", "black");
  
  defs.append("marker")
      .attr("id", "arrowhead-left")
      .attr("viewBox", "0 -5 10 10")
      .attr("refX", 6)
      .attr("refY", 0)
      .attr("markerWidth", 4)
      .attr("markerHeight", 4)
      .attr("orient", "auto")
    .append("path")
      .attr("d", "M10,-5L0,0L10,5")
      .attr("fill", "black");
  
  // Load and parse the gene annotation file.
  d3.text(config.annotationFile).then(function(text) {
    var lines = text.split("\n").filter(line => line.trim().length > 0);
    var genes = lines.map(function(line) {
      var parts = line.split("\t");
      var coordParts = parts[2].split(/[:\-]/); // e.g., ["chr19", "58345178", "58353492"]
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
    
    // Filter genes: same chromosome, overlapping region, protein_coding only.
    genes = genes.filter(function(d) {
      return d.chrom === config.region.chr &&
             d.end >= config.region.start &&
             d.start <= config.region.end &&
             d.type === "protein_coding";
    });
    console.log("Number of protein coding genes after filtering:", genes.length);
    
    // *** If deletion mode, filter out genes that lie entirely within the deletion region ***
    if (inDeletionMode) {
      genes = genes.filter(function(d) {
        // Remove gene if its entire span is inside the deletion.
        return (d.end < config.deletionStart || d.start > config.deletionEnd);
      });
    }
    
    // Sort genes by start coordinate.
    genes.sort((a, b) => a.start - b.start);
    
    // Row assignment: assign genes to rows if they are close.
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
    const rowHeight = numRows > 0 ? height / numRows : height;
    
    // Draw gene lines.
    g.selectAll("line.gene")
     .data(genes)
     .enter()
     .append("line")
     .attr("class", "gene")
     .attr("x1", d => xScale(Math.max(d.start, config.region.start)))
     .attr("x2", d => xScale(Math.min(d.end, config.region.end)))
     .attr("y1", d => d.row * rowHeight + rowHeight / 2)
     .attr("y2", d => d.row * rowHeight + rowHeight / 2)
     .attr("stroke", "black")
     .attr("stroke-width", Math.min(rowHeight * 0.2, 1.5))
     .attr("stroke-linecap", "round")
     .attr("marker-end", d => d.strand === "+" ? "url(#arrowhead-right)" : null)
     .attr("marker-start", d => d.strand === "-" ? "url(#arrowhead-left)" : null);
    
    // Draw gene labels.
    g.selectAll("text.gene-label")
     .data(genes)
     .enter()
     .append("text")
     .attr("class", "gene-label")
     .attr("x", function(d) {
       // Compute center using discontinuous xScale.
       const x1 = xScale(Math.max(d.start, config.region.start));
       const x2 = xScale(Math.min(d.end, config.region.end));
       let center = (x1 + x2) / 2;
       // Clamp center to avoid labels too close to the edges.
       if (center < 5) center = 5;
       if (center > fullWidth - 5) center = fullWidth - 5;
       return center;
     })
     .attr("y", d => d.row * rowHeight + rowHeight / 2 - 4)
     .attr("text-anchor", "middle")
     .attr("font-size", Math.min(rowHeight * 0.4, 12))
     .text(d => d.gene);
    
    // Optionally, add an x-axis if desired.
    
  }).catch(function(error) {
    console.error("Error loading gene annotation data:", error);
  });
}


/** Deletion-mode signal plot with discontinuous (broken) x-axis. */
function drawLineChartWithBreak(containerSelector, config) {
  const margin = { top: 20, right: 20, bottom: 40, left: 50 };
  const { leftDomain, rightDomain } = config.xAxis;
  const leftSpan = leftDomain[1] - leftDomain[0];
  const rightSpan = rightDomain[1] - rightDomain[0];
  const totalSpan = leftSpan + rightSpan;
  // Set gap to 0 so the two halves are adjacent.
  const gapPx = 0;
  
  // The available width for the two segments (in pixels)
  const availableWidth = totalSpan * PX_PER_MB;
  const height = config.chart.height || 150;
  const canvasWidth = availableWidth + margin.left + margin.right + gapPx;
  const canvasHeight = height + margin.top + margin.bottom;
  
  const { ctx } = createHiResCanvas(containerSelector, canvasWidth, canvasHeight);
  ctx.translate(margin.left, margin.top);
  
  // Compute pixel widths for left and right segments
  const leftWidthPx = (leftSpan / totalSpan) * availableWidth;
  const rightWidthPx = (rightSpan / totalSpan) * availableWidth;
  
  // Build x-scales for each segment with a fixed tick step of 0.2 Mb.
  const scaleLeft = d3.scaleLinear()
    .domain(leftDomain)
    .range([0, leftWidthPx]);
  const scaleRight = d3.scaleLinear()
    .domain(rightDomain)
    .range([0, rightWidthPx]);
  
  // Generate ticks for each segment using a fixed step of 0.2 Mb.
  const leftTicks = d3.range(leftDomain[0], leftDomain[1] + 0.001, 0.2);
  let rightTicks = d3.range(rightDomain[0], rightDomain[1] + 0.001, 0.2);
  // Remove the first tick from the right side to avoid overlapping
  if (rightTicks.length && Math.abs(rightTicks[0] - rightDomain[0]) < 1e-6) {
    rightTicks.shift();
  }
  
  // Data: filter to include only points in left or right segments.
  const data = config.series[0].data || [];
  const filtered = data.filter(d => {
    const xVal = d[0];
    return ((xVal >= leftDomain[0] && xVal <= leftDomain[1]) ||
            (xVal >= rightDomain[0] && xVal <= rightDomain[1]));
  });
  
  // Y-scale based on filtered data.
  const yMax = d3.max(filtered, d => d[1]) || 1;
  const yMin = d3.min(filtered, d => d[1]) || 0;
  const yScale = d3.scaleLinear()
    .domain([yMin, yMax])
    .nice()
    .range([height, 0]);
  
  // Helper function: maps x value (Mb) to pixel position.
  function xScale(xVal) {
    if (xVal >= leftDomain[0] && xVal <= leftDomain[1]) {
      return scaleLeft(xVal);
    } else if (xVal >= rightDomain[0] && xVal <= rightDomain[1]) {
      return leftWidthPx + gapPx + scaleRight(xVal);
    }
    return NaN;
  }
  
  // Draw y-axis (ticks only)
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(0, height);
  ctx.strokeStyle = "#000";
  ctx.stroke();
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  ctx.font = "10px sans-serif";
  ctx.fillStyle = "#000";
  yScale.ticks(5).forEach(tick => {
    const yPos = yScale(tick);
    ctx.beginPath();
    ctx.moveTo(0, yPos);
    ctx.lineTo(-6, yPos);
    ctx.stroke();
    ctx.fillText(tick, -8, yPos);
  });
  
  // Draw left x-axis
  ctx.save();
  ctx.translate(0, height);
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(leftWidthPx, 0);
  ctx.stroke();
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  leftTicks.forEach(tick => {
    const xPos = scaleLeft(tick);
    ctx.beginPath();
    ctx.moveTo(xPos, 0);
    ctx.lineTo(xPos, 6);
    ctx.stroke();
    ctx.fillText(d3.format(".2f")(tick), xPos, 6);
  });
  ctx.restore();
  
  // Draw right x-axis
  ctx.save();
  ctx.translate(leftWidthPx + gapPx, height);
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(rightWidthPx, 0);
  ctx.stroke();
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  rightTicks.forEach(tick => {
    const xPos = scaleRight(tick);
    ctx.beginPath();
    ctx.moveTo(xPos, 0);
    ctx.lineTo(xPos, 6);
    ctx.stroke();
    ctx.fillText(d3.format(".2f")(tick), xPos, 6);
  });
  ctx.restore();
  
  // Draw a very small zig-zag break marker (optional)
  ctx.save();
  ctx.translate(leftWidthPx, height);
  ctx.beginPath();
  ctx.moveTo(-3, 0);
  ctx.lineTo(0, 8);
  ctx.lineTo(3, 0);
  ctx.strokeStyle = "#000";
  ctx.lineWidth = 1.5;
  ctx.stroke();
  ctx.restore();
  
  // Draw the signal line
  ctx.beginPath();
  filtered.forEach((d, i) => {
    const px = xScale(d[0]);
    const py = yScale(d[1]);
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
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
    const totalWidth = leftWidthPx + gapPx + rightWidthPx;
    ctx.fillText(config.xAxis.title, totalWidth / 2, height + margin.bottom - 5);
    ctx.restore();
  }
}
