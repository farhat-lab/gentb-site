/*
 * Copyright 2017, Maha Farhat
 *
 * This file is part of the software inkscape-web, consisting of custom 
 * code for the Inkscape project's django-based website.
 *
 * inkscape-web is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inkscape-web is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with inkscape-web.  If not, see <http://www.gnu.org/licenses/>.
 */


$(document).ready(function() {
  var svg = 'svg.lineages';
  var chart = initialiseLineageChart(svg);
  $('#lineage-store').data('json-signal', function(data) {
    console.log("Lineage ahoy!");
    chartLineageData(svg, chart, data.data);
  }); 
});

function chartLineageData(svg, chart, data) {
    var svg_n = chartData(svg, chart, data);

    svg_n.transition().duration(0);
    // Select the first category by "clicking" on it
    var evt = new MouseEvent("click");
    var node = d3.select(svg).selectAll('.nv-series')
        .filter(function(d, i){return i == 0;}).node().dispatchEvent(evt);

    svg_n.transition().duration(1200);
}

function initialiseLineageChart(svg) {
    var chart = nv.models.multiBarChart()
      .stacked(true)
      .reduceXTicks(false);

    var width = 1000;
    var height = 600;

    chart.margin({top: 20, right: 0, bottom: 120, left: 80});
    chart.height(height);
    chart.width(width);
    chart.yAxis.scale(100).orient("left")
        .axisLabel('Number of strains')
        .tickFormat(d3.format("d"));

    // Single click selection of chart type
    chart.legend.radioButtonMode(true);
    // Disable the stacked/unstacked option
    chart.showControls(false);

    chart.xAxis
        .rotateLabels(-20);

    chart.showLegend(true);

    var svg = d3.select(svg)
          .attr('perserveAspectRatio', 'xMinYMid')
          .attr('width', width)
          .attr('height', height)
          .attr('viewBox', '0 0 ' + width + ' ' + height)
          .datum([])
          .call(chart);

    chart.multibar.dispatch.on("elementClick", function(e) {
        setTabData('lineage', e.data.x, e.data.x, 'ok-circle')
    });

    $('#lineages').parent().click(function(e) {
      unsetTabData('lineage');
    });
    return chart;
}


