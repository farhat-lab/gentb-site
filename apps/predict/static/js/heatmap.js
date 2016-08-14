
function scatter_plot(data) {
  /*  
   * The idea here is to update an existing plot with new data.
   */
  var svg = d3.select('#scatter svg');

  svg[0][0].__data__ = data;
  var nvChart = svg[0][0].__chart__;

  nvChart.xAxis
      .axisLabel('Genetic Region')
      .tickFormat(function(d) { return data[0].cols[d];})
      .tickValues([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21])
      .rotateLabels(-45);
  svg.data([data]).transition().duration(500).call(nvChart);
  nv.utils.windowResize(nvChart.update);
}

function make_heatmap(data, heatmap_id, scatter_id, no_data) {
  var scat = data["matrix"]["scatter"];
  var options = {
    "no_data": no_data,
    "xaxis_height": 80,
    "yaxis_width": 120,
    "xaxis_font_size": null,
    "yaxis_font_size": null,
    "brush_color": "#0000FF",
    "show_grid": true,
    "anim_duration": 500,
    "on_click": function(elem, x, x_label, y, y_label, datum) {
      $('.datapt.selected').attr('class', 'datapt');
      $(elem).attr('class', 'datapt selected');
      $('#scatter_title').text('Mutation plot for drug='+x_label+', strain=' +y_label);
      scatter_plot(scat.data[y_label][x]);
    },
  }
  heatmap(heatmap_id, data, options);

  nv.addGraph(function() {
    var data = Array();
    var chart = nv.models.scatterChart()
      .yDomain([0, 4]);

    var width = 1000;
    var height = 300;

    chart.noData(options.no_data);
    chart.margin({top: 20, right: 60, bottom: 60, left: 60});
    chart.height(height);
    chart.width(width);
    chart.yAxis.scale(100).orient("left")
  .axisLabel('Number of mutations')
  .tickFormat(d3.format("d"))
  .tickValues([0,1,2,3,4,5]);

    chart.showLegend(true);

    var svg = d3.select(scatter_id + ' svg')
  .attr('perserveAspectRatio', 'xMinYMid')
  .attr('width', width)
  .attr('height', height)
  .datum(data)
  .attr('viewBox', '0 0 ' + width + ' ' + height)
  .transition().duration(1200)
  .call(chart);

    chart.tooltip.contentGenerator(function (data) {
      ret = "<table>";
      $.each(data.point.tip, function(x, value) {
        ret += "<tr><th align=\"right\">" + value + "</th></tr>";
      });
      return ret + "</table>";
    });
    svg[0][0].__chart__ = chart;
  });
}

function heatmap(selector, data, options) {
  var merged = [];
  for (var i = 0; i < data.matrix.data.length; i++) {
    merged.push({
      label: data.matrix.data[i],
      opacity: parseFloat(data.matrix.data[i]),
    })
  }
  data.matrix.merged = merged;

  // ==== BEGIN HELPERS =================================
  
  function htmlEscape(str) {
    return (str+"").replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }
  
  // Given a list of widths/heights and a total width/height, provides
  // easy access to the absolute top/left/width/height of any individual
  // grid cell. Optionally, a single cell can be specified as a "fill"
  // cell, meaning it will take up any remaining width/height.
  // 
  // rows and cols are arrays that contain numeric pixel dimensions,
  // and up to one "*" value.
  function GridSizer(widths, heights, /*optional*/ totalWidth, /*optional*/ totalHeight) {
    this.widths = widths;
    this.heights = heights;
  
    var fillColIndex = null;
    var fillRowIndex = null;
    var usedWidth = 0;
    var usedHeight = 0;
    var i;
    for (i = 0; i < widths.length; i++) {
      if (widths[i] === "*") {
        if (fillColIndex !== null) {
          throw new Error("Only one column can be designated as fill");
        }
        fillColIndex = i;
      } else {
        usedWidth += widths[i];
      }
    }
    if (fillColIndex !== null) {
      widths[fillColIndex] = totalWidth - usedWidth;
    } else {
      if (typeof(totalWidth) === "number" && totalWidth !== usedWidth) {
        throw new Error("Column widths don't add up to total width");
      }
    }
    for (i = 0; i < heights.length; i++) {
      if (heights[i] === "*") {
        if (fillRowIndex !== null) {
          throw new Error("Only one row can be designated as fill");
        }
        fillRowIndex = i;
      } else {
        usedHeight += heights[i];
      }
    }
    if (fillRowIndex !== null) {
      heights[fillRowIndex] = totalHeight - usedHeight;
    } else {
      if (typeof(totalHeight) === "number" && totalHeight !== usedHeight) {
        throw new Error("Column heights don't add up to total height");
      }
    }
  }
  
  GridSizer.prototype.getCellBounds = function(x, y) {
    if (x < 0 || x >= this.widths.length || y < 0 || y >= this.heights.length)
      throw new Error("Invalid cell bounds");
  
    var left = 0;
    for (var i = 0; i < x; i++) {
      left += this.widths[i];
    }
  
    var top = 0;
    for (var j = 0; j < y; j++) {
      top += this.heights[j];
    }
  
    return {
      width: this.widths[x],
      height: this.heights[y],
      top: top,
      left: left
    }
  }
  
  // ==== END HELPERS ===================================


  var el = d3.select(selector);

  var bbox = el.node().getBoundingClientRect();

  var Controller = function() {
    this._events = d3.dispatch("highlight", "datapoint_hover", "transform");
    this._highlight = {x: null, y: null};
    this._datapoint_hover = {x: null, y: null, value: null};
    this._transform = null;
  };
  (function() {
    this.highlight = function(x, y) {
      // Copy for safety
      if (!arguments.length) return {x: this._highlight.x, y: this._highlight.y};

      if (arguments.length == 1) {
        this._highlight = x;
      } else {
        this._highlight = {x: x, y: y};
      }
      this._events.highlight.call(this, this._highlight);
    };

    this.datapoint_hover = function(_) {
      if (!arguments.length) return this._datapoint_hover;
      
      this._datapoint_hover = _;
      this._events.datapoint_hover.call(this, _);
    };

    this.transform = function(_) {
      if (!arguments.length) return this._transform;
      this._transform = _;
      this._events.transform.call(this, _);
    };

    this.on = function(evt, callback) {
      this._events.on(evt, callback);
    };
  }).call(Controller.prototype);

  var controller = new Controller();
  var width = options.width || bbox.width;
  var height = options.height || bbox.height;
  var num_rows = data.matrix.rows.length;
  var num_cols = data.matrix.cols.length;

  // Set option defaults
  var opts = {};
  options = options || {};
  opts.xaxis_height = options.xaxis_height || 80;
  opts.yaxis_width = options.yaxis_width || 120;

  // Make the squares, square
  if(opts.height / num_rows < opts.width / num_cols) {
    opts.height = height
    opts.width = (height / num_rows) * num_cols;
  } else {
    opts.height = (width / num_cols) * num_rows;
    opts.width = width;
  }
  opts.height += opts.xaxis_height;
  opts.width += opts.yaxis_width;

  opts.xclust_height = options.xclust_height || opts.height * 0.12;
  opts.yclust_width = options.yclust_width || opts.width * 0.12;
  opts.link_color = opts.link_color || "#AAA";
  opts.axis_padding = options.axis_padding || 6;
  opts.show_grid = options.show_grid;
  if (typeof(opts.show_grid) === 'undefined') {
    opts.show_grid = true;
  }
  opts.brush_color = options.brush_color || "#0000FF";
  opts.xaxis_font_size = options.xaxis_font_size;
  opts.yaxis_font_size = options.yaxis_font_size;
  opts.anim_duration = options.anim_duration;
  if (typeof(opts.anim_duration) === 'undefined') {
    opts.anim_duration = 500;
  }

  if (!data.rows) {
    opts.yclust_width = 0;
  }
  if (!data.cols) {
    opts.xclust_height = 0;
  }
  
  // Change axis location here, by shifting it's size, then it's grid position.
  var gridSizer = new GridSizer(
    [opts.yclust_width, opts.yaxis_width, "*"],
    [opts.xclust_height, "*", opts.xaxis_height],
    opts.width,
    opts.height
  );

  var cBounds = gridSizer.getCellBounds(2, 1);
  var colDendBounds = gridSizer.getCellBounds(1, 0);
  var rowDendBounds = gridSizer.getCellBounds(0, 1);
  var yBound = gridSizer.getCellBounds(1, 1);
  var xBound = gridSizer.getCellBounds(2, 2);

  var rd = el.select('svg.rowDend');
  var cd = el.select('svg.colDend');
  var rw = rowDendBounds.width, rh = rowDendBounds.height;
  var cw = colDendBounds.width, ch = colDendBounds.height;
  var pad = opts.axis_padding;
  var row = !data.rows ? null : dendrogram(rd, data.rows, false, rw, rh, pad);
  var col = !data.cols ? null : dendrogram(cd, data.cols, true, cw, ch, pad);

  function cssify(styles) {
    return {
      position: "absolute",
      top: styles.top + "px",
      left: styles.left + "px",
      width: styles.width + "px",
      height: styles.height + "px"
    };
  }

  // Create DOM structure
  (function() {
    var inner = el.append("div").classed("inner", true);
    inner.style("width", opts.width + "px");
    inner.style("height", opts.height + "px");
    var info = inner.append("div").classed("info", true);
    var colDend = inner.append("svg").classed("dendrogram colDend", true)
                                     .style(cssify(colDendBounds));
    var rowDend = inner.append("svg").classed("dendrogram rowDend", true)
                                     .style(cssify(rowDendBounds));
    var colmap = inner.append("svg").classed("colormap", true)
                                     .style(cssify(cBounds));
    var xaxis = inner.append("svg").classed("axis xaxis", true)
                                   .style(cssify(xBound));
    var yaxis = inner.append("svg").classed("axis yaxis", true)
                                   .style(cssify(yBound));
    
    // Hack the width of the x-axis to allow x-overflow of rotated labels; the
    // QtWebkit viewer won't allow svg elements to overflow:visible.
    xaxis.style("width", (opts.width - opts.yclust_width) + "px");
    xaxis
      .append("defs")
        .append("clipPath").attr("id", "xaxis-clip")
          .append("polygon")
            .attr("points", "" + [
              [0, 0],
              [xBound.width, 0],
              [xBound.width + yBound.width, xBound.height],
              [0, xBound.height]
            ]);
    xaxis.node(0).setAttribute("clip-path", "url(#xaxis-clip)");

    inner.on("click", function() {
      controller.highlight(null, null);
    });
    controller.on('highlight.inner', function(hl) {
      inner.classed('highlighting',
        typeof(hl.x) === 'number' || typeof(hl.y) === 'number');
    });
  })();
  
  var cm = el.select('svg.colormap');
  var colormap = colormap(cm, data.matrix, cBounds.width, cBounds.height);

  var xa = el.select('svg.xaxis');
  var ya = el.select('svg.yaxis');
  var cols = data.cols || data.matrix.cols;
  var rows = data.rows || data.matrix.rows;

  // This isn't right, shouldn't be here. XXX
  $.each(cols, function(index, item) {
      cols[index] = item.toUpperCase();
  });

  var xax = axisLabels(xa, cols, true, xBound.width, xBound.height, pad);
  var yax = axisLabels(ya, rows, false, yBound.width, yBound.height, pad);
  
  function colormap(svg, data, width, height) {
    // Check for no data
    if (data.length === 0)
      return function() {};

        if (!opts.show_grid) {
          svg.style("shape-rendering", "crispEdges");
        }
 
    var cols = data.dim[1];
    var rows = data.dim[0];
    
    var merged = data.merged;
    
    var x = d3.scale.linear().domain([0, cols]).range([0, width]);
    var y = d3.scale.linear().domain([0, rows]).range([0, height]);
    var tip = d3.tip()
        .attr('class', 'd3heatmap-tip')
        .html(function(d, i) {
          var index = d.row * data.cols.length + d.col;
          // XXX This is specific to one use case, and isn't good here.
          return "<table>" + 
            "<tr><th align=\"right\">Strain</th><td>" + htmlEscape(data.rows[d.row]) + "</td></tr>" +
            "<tr><th align=\"right\">Drug</th><td>" + htmlEscape(data.cols[d.col]) + "</td></tr>" +
            "<tr><th align=\"right\">DR Probability</th><td>" + htmlEscape(d.label) + "</td></tr>" +
            "<tr><th align=\"right\">FP Rate</th><td>" + data.extra[index][0] + "</td></tr>" +
            "<tr><th align=\"right\">FN Rate</th><td>" + data.extra[index][1] + "</td></tr>" +
            "</table>";
        })
        .direction("se")
        .style("position", "absolute");
    
    var brush = d3.svg.brush()
        .x(x)
        .y(y)
        .clamp([true, true])
        .on('brush', function() {
          var extent = brush.extent();
          extent[0][0] = Math.round(extent[0][0]);
          extent[0][1] = Math.round(extent[0][1]);
          extent[1][0] = Math.round(extent[1][0]);
          extent[1][1] = Math.round(extent[1][1]);
          d3.select(this).call(brush.extent(extent));
        })
        .on('brushend', function() {

          if (brush.empty()) {
            controller.transform({
              scale: [1,1],
              translate: [0,0],
              extent: [[0,0],[cols,rows]]
            });
          } else {
            var tf = controller.transform();
            var ex = brush.extent();
            var scale = [
              cols / (ex[1][0] - ex[0][0]),
              rows / (ex[1][1] - ex[0][1])
            ];
            var translate = [
              ex[0][0] * (width / cols) * scale[0] * -1,
              ex[0][1] * (height / rows) * scale[1] * -1
            ];
            controller.transform({scale: scale, translate: translate, extent: ex});
          }
          brush.clear();
          d3.select(this).call(brush).select(".brush .extent")
              .style({fill: opts.brush_color, stroke: opts.brush_color});
        });

    svg = svg
        .attr("width", width)
        .attr("height", height);
    var rect = svg.selectAll("rect").data(merged);
    rect.enter().append("rect").classed("datapt", true)
        .property("colIndex", function(d, i) { return i % cols; })
        .property("rowIndex", function(d, i) { return Math.floor(i / cols); })
        .property("value", function(d, i) { return d.value; })
        .attr("fill-opacity", function(d) {
          if (!d.opacity) {
            return "0.0";
          }
          return d.opacity;
        });
    rect.exit().remove();
    rect.append("title")
        .text(function(d, i) { return d.label; });
    rect.call(tip);

    var spacing;
    if (typeof(opts.show_grid) === 'number') {
      spacing = opts.show_grid;
    } else if (!!opts.show_grid) {
      spacing = 0.25;
    } else {
      spacing = 0;
    }
    function draw(selection) {
      selection
          .attr("x", function(d, i) {
            return x(i % cols);
          })
          .attr("y", function(d, i) {
            return y(Math.floor(i / cols));
          })
          .attr("width", (x(1) - x(0)) - spacing)
          .attr("height", (y(1) - y(0)) - spacing);
    }

    draw(rect);

    controller.on('transform.colormap', function(_) {
      x.range([_.translate[0], width * _.scale[0] + _.translate[0]]);
      y.range([_.translate[1], height * _.scale[1] + _.translate[1]]);
      draw(rect.transition().duration(opts.anim_duration).ease("linear"));
    });
    

    var brushG = svg.append("g");
    brushG.attr('class', 'brush');
    brushG.call(brush);
    if(brush.event) {
      brushG.call(brush.event);
    }
    brushG.select("rect.datapt")

    brushG.select("rect.background")
        .on("mouseenter", function() {
          tip.style("display", "block");
        })
        .on("click", function() {
          if(options['on_click']) {
            var offsetX = d3.event.offsetX;
            var offsetY = d3.event.offsetY;

            var col = Math.floor(x.invert(offsetX));
            var row = Math.floor(y.invert(offsetY));
            var index = row * cols + col;

            var rect = $(svg[0]).children('rect.datapt')[row*cols + col];

            options['on_click'](rect, col, data.cols[col], row, data.rows[row], merged[index]);
          }
        })
        .on("mousemove", function() {
          var e = d3.event;
          var offsetX = d3.event.offsetX;
          var offsetY = d3.event.offsetY;
          if (typeof(offsetX) === "undefined") {
            // Firefox 38 and earlier
            var target = e.target || e.srcElement;
            var rect = target.getBoundingClientRect();
            offsetX = e.clientX - rect.left,
            offsetY = e.clientY - rect.top;
          }
          
          var col = Math.floor(x.invert(offsetX));
          var row = Math.floor(y.invert(offsetY));
          var label = merged[row*cols + col].label;
          tip.show({col: col, row: row, label: label}).style({
            top: d3.event.clientY + 15 + "px",
            left: d3.event.clientX + 15 + "px",
            opacity: 0.9
          });
          controller.datapoint_hover({col:col, row:row, label:label});
        })
        .on("mouseleave", function() {
          tip.hide().style("display", "none");
          controller.datapoint_hover(null);
        });

    controller.on('highlight.datapt', function(hl) {
      rect.classed('highlight', function(d, i) {
        return (this.rowIndex === hl.y) || (this.colIndex === hl.x);
      });
    });
  }

  function axisLabels(svg, data, rotated, width, height, padding) {
    if(!rotated) {
      // Used for axis on the left instead of right.
      svg.attr("transform", "translate(" + (width - (padding * 2)) + ", 0)");
    }
    svg = svg.append('g');

    // The data variable is either cluster info, or a flat list of names.
    // If the former, transform it to simply a list of names.
    var leaves;
    if (data.children) {
      leaves = d3.layout.cluster().nodes(data)
          .filter(function(x) { return !x.children; })
          .map(function(x) { return x.label + ""; });
    } else if (data.length) {
      leaves = data;
    }
    
    // Define scale, axis
    var scale = d3.scale.ordinal()
        .domain(leaves)
        .rangeBands([0, rotated ? width : height]);
    var axis = d3.svg.axis();
    axis.scale(scale);
    axis.orient(rotated ? "bottom" : "left");
    //axis.outerTickSize(0);
    axis.tickPadding(padding);
    axis.tickValues(leaves);

    // Create the actual axis
    var axisNodes = svg.append("g")
        .attr("transform", rotated ? "translate(0," + padding + ")" : "translate(" + padding + ",0)")
        .call(axis);
    var fontSize = opts[(rotated ? 'x' : 'y') + 'axis_font_size']
        || Math.min(18, Math.max(9, scale.rangeBand() - (rotated ? 11: 8))) + "px";
    axisNodes.selectAll("text").style("font-size", fontSize);
    
    var mouseTargets = svg.append("g")
      .selectAll("g").data(leaves);
    mouseTargets
      .enter()
        .append("g").append("rect")
          .attr("transform", rotated ? "rotate(45),translate(0,0)" : "")
          .attr("fill", "transparent")
          .on("click", function(d, i) {
            var dim = rotated ? 'x' : 'y';
            var hl = controller.highlight() || {x:null, y:null};
            if (hl[dim] == i) {
              // If clicked already-highlighted row/col, then unhighlight
              hl[dim] = null;
              controller.highlight(hl);
            } else {
              hl[dim] = i;
              controller.highlight(hl);
            }
            d3.event.stopPropagation();
          });
    function layoutMouseTargets(selection) {
      selection
          .attr("transform", function(d, i) {
            var x = rotated ? scale(d) + scale.rangeBand()/2 : 0;
            var y = rotated ? padding + 6 : scale(d);
            return "translate(" + x + "," + y + ")";
          })
        .selectAll("rect")
          .attr("height", scale.rangeBand() / (rotated ? 1.414 : 1))
          .attr("width", rotated ? height * 1.414 * 1.2 : width);
    }
    layoutMouseTargets(mouseTargets);

    if (rotated) {
      axisNodes.selectAll("text")
        .attr("transform", "rotate(45),translate(6, 0)")
        .style("text-anchor", "start");
    }
    
    controller.on('highlight.axis-' + (rotated ? 'x' : 'y'), function(hl) {
      var ticks = axisNodes.selectAll('.tick');
      var selected = hl[rotated ? 'x' : 'y'];
      if (typeof(selected) !== 'number') {
        ticks.classed('faded', false);
        return;
      }
      ticks.classed('faded', function(d, i) {
        return i !== selected;
      });
    });

    controller.on('transform.axis-' + (rotated ? 'x' : 'y'), function(_) {
      var dim = rotated ? 0 : 1;
      //scale.domain(leaves.slice(_.extent[0][dim], _.extent[1][dim]));
      var rb = [_.translate[dim], (rotated ? width : height) * _.scale[dim] + _.translate[dim]];
      scale.rangeBands(rb);
      var tAxisNodes = axisNodes.transition().duration(opts.anim_duration).ease('linear');
      tAxisNodes.call(axis);
      // Set text-anchor on the non-transitioned node to prevent jumpiness
      // in RStudio Viewer pane
      tAxisNodes.selectAll("g")
          .style("opacity", function(d, i) {
            if (i >= _.extent[0][dim] && i < _.extent[1][dim]) {
              return 1;
            } else {
              return 0;
            }
          });
      mouseTargets.transition().duration(opts.anim_duration).ease('linear')
          .call(layoutMouseTargets)
          .style("opacity", function(d, i) {
            if (i >= _.extent[0][dim] && i < _.extent[1][dim]) {
              return 1;
            } else {
              return 0;
            }
          });
    });

  }
  
  function edgeStrokeWidth(node) {
    if (node.edgePar && node.edgePar.lwd)
      return node.edgePar.lwd;
    else
      return 1;
  }
  
  function maxChildStrokeWidth(node, recursive) {
    var max = 0;
    for (var i = 0; i < node.children.length; i++) {
      if (recursive) {
        max = Math.max(max, maxChildStrokeWidth(node.children[i], true));
      }
      max = Math.max(max, edgeStrokeWidth(node.children[i]));
    }
    return max;
  }
  
  function dendrogram(svg, data, rotated, width, height, padding) {
    var topLineWidth = maxChildStrokeWidth(data, false);
        //maxChildStrokeWidth(data, false);
    
    var x = d3.scale.linear()
        .domain([data.height, 0])
        .range([topLineWidth/2, width-padding]);
    var y = d3.scale.linear()
        .domain([0, height])
        .range([0, height]);
    
    var cluster = d3.layout.cluster()
        .separation(function(a, b) { return 1; })
        .size([rotated ? width : height, NaN]);
    
    var transform = "translate(1,0)";
    if (rotated) {
      // Flip dendrogram vertically
      x.range([topLineWidth/2, -height+padding+2]);
      // Rotate
      transform = "rotate(-90) translate(-2,0)";
    }

    var dendrG = svg
        .attr("width", width)
        .attr("height", height)
      .append("g")
        .attr("transform", transform);
    
    var nodes = cluster.nodes(data),
        links = cluster.links(nodes);

    // I'm not sure why, but after the heatmap loads the "links"
    // array mutates to much smaller values. I can't figure out
    // what's doing it, so instead we just make a deep copy of
    // the parts we want.
    var links1 = links.map(function(link, i) {
      return {
        source: {x: link.source.x, y: link.source.height},
        target: {x: link.target.x, y: link.target.height},
        edgePar: link.target.edgePar
      };
    });
    
    var lines = dendrG.selectAll("polyline").data(links1);
    lines
      .enter().append("polyline")
        .attr("class", "link")
        .attr("stroke", function(d, i) {
          if (!d.edgePar.col) {
            return opts.link_color;
          } else {
            return d.edgePar.col;
          }
        })
        .attr("stroke-width", edgeStrokeWidth)
        .attr("stroke-dasharray", function(d, i) {
          var pattern;
          switch (d.edgePar.lty) {
            case 6:
              pattern = [3,3,5,3];
              break;
            case 5:
              pattern = [15,5];
              break;
            case 4:
              pattern = [2,4,4,4];
              break;
            case 3:
              pattern = [2,4];
              break;
            case 2:
              pattern = [4,4];
              break;
            case 1:
            default:
              pattern = [];
              break;
          }
          for (var i = 0; i < pattern.length; i++) {
            pattern[i] = pattern[i] * (d.edgePar.lwd || 1);
          }
          return pattern.join(",");
        });

    function draw(selection) {
      function elbow(d, i) {
        return x(d.source.y) + "," + y(d.source.x) + " " +
            x(d.source.y) + "," + y(d.target.x) + " " +
            x(d.target.y) + "," + y(d.target.x);
      }
      
      selection
          .attr("points", elbow);
    }

    controller.on('transform.dendr-' + (rotated ? 'x' : 'y'), function(_) {
      var scaleBy = _.scale[rotated ? 0 : 1];
      var translateBy = _.translate[rotated ? 0 : 1];
      y.range([translateBy, height * scaleBy + translateBy]);
      draw(lines.transition().duration(opts.anim_duration).ease("linear"));
    });

    draw(lines);
  }

 
  var dispatcher = d3.dispatch('hover', 'click');
  
  controller.on("datapoint_hover", function(_) {
    dispatcher.hover({data: _});
  });
  
  function on_col_label_mouseenter(e) {
    controller.highlight(+d3.select(this).attr("index"), null);
  }
  function on_col_label_mouseleave(e) {
    controller.highlight(null, null);
  }
  function on_row_label_mouseenter(e) {
    controller.highlight(null, +d3.select(this).attr("index"));
  }
  function on_row_label_mouseleave(e) {
    controller.highlight(null, null);
  }

  return {
    on: function(type, listener) {
      dispatcher.on(type, listener);
      return this;
    }
  };
}


