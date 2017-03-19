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
  var svg = $('svg.mutations');
  makeMutationChart();
});

function makeMutationChart() {
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
        .axisLabel("Mutation Name")
        .rotateLabels(-20);

    chart.showLegend(true);

    var data = [{'values': [], key: 'Resistant'}];

    var svg = d3.select('svg.mutations')
          .attr('perserveAspectRatio', 'xMinYMid')
          .attr('width', width)
          .attr('height', height)
          .attr('viewBox', '0 0 ' + width + ' ' + height)
          .datum(data)
          .transition().duration(500)
          .call(chart);

    var ready = false;
    var x = 0;
    setInterval(function(){
     if(!ready){return}
     data[0].values.push({ x: 'SNP'+x, y: Math.random() * 100});
     if (data[0].values.length > 5) data[0].values.shift();
     console.log("VAL");
     d3.select('svg.mutations').datum(data).call(chart);
     x++;
    }, 3000);

    nv.utils.windowResize(chart.update);
    ready = true;
}


