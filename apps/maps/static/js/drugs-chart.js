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
  var svg = $('svg.drugs');
  $.getJSON(svg.data('drug'), function(data) {
    makeDrugChart(data.data);
  });
});

function makeDrugChart(data) {
    var chart = nv.models.multiBarChart()
      .reduceXTicks(false);

    var width = 1000;
    var height = 600;

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

    chart.tooltip.contentGenerator(function (data) {
      ret = "<table>";
      if(data.data) {
        ret += "<tr><th align=\"right\">" + data.data.x + "</th></tr>";
        ret += "<tr><th align=\"right\">" + data.data.y + "</th></tr>";
      }   
      return ret + "</table>";
    });

    var svg = d3.select('svg.drugs')
          .attr('perserveAspectRatio', 'xMinYMid')
          .attr('width', width)
          .attr('height', height)
          .attr('viewBox', '0 0 ' + width + ' ' + height)
          .datum(data)
          .transition().duration(1200)
          .call(chart);
}


