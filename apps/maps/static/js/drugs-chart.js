/*
 * Copyright 2017, Maha Farhat
 *
 * This file is part of the software gentb, consisting of custom 
 * code for the GenTB's django-based website.
 *
 * gentb is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gentb is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with gentb.  If not, see <http://www.gnu.org/licenses/>.
 */

$(document).ready(function() {
  var svg = 'svg.drugs';
  var chart = initialiseDrugChart(svg);
  $('#drug-store').data('json-signal', function(data) {
    console.log("Drugging up!");
    chartData(svg, chart, data.data);
  });
});

function initialiseDrugChart(svg) {

    //So we first make the chart
    var chart = nv.models.multiBarChart()
      .stacked(true)
      .reduceXTicks(false);


    var width = 1000;
    var height = 600;

    //Some properties of the chart
    chart.margin({top: 20, right: 0, bottom: 60, left: 80});
    chart.height(height);
    chart.width(width);
    chart.yAxis.scale(100).orient("left")
        .axisLabel('Number of strains')
        .tickFormat(d3.format("d"));

    chart.xAxis
        .axisLabel("Drug Codename")
        .rotateLabels(-45);

    chart.showLegend(true);

    //use d3 (a javascript library) to select the svg tag in the html that you made for the chart and then to put in the chart
    var svg = d3.select(svg)
          .attr('perserveAspectRatio', 'xMinYMid')
          .attr('width', width)
          .attr('height', height)
          .attr('viewBox', '0 0 ' + width + ' ' + height)
          .datum([])
          .transition().duration(1200)
          .call(chart);

    chart.multibar.dispatch.on("elementClick", function(e) {
        setTabData('drug', e.data.x, e.data.x, 'map-marker')
    });

    $('#drugs').parent().click(function(e) {
      unsetTabData('drug');
    });

    return chart;
}
