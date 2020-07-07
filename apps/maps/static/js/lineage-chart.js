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

  $('#lineage-store').data('json-signal', function(data) {
    $("#row1").empty(); $("#row2").empty();
    chartLineageData(data);
  });
});

function chartLineageData(data) {
  var numPlots = data.children.length;
  var height = numPlots <= 2 ? 600 : 300;

  // Splits the charts across two rows, with the bottom row having more elements when `numPlots` is odd
  var counts = {'#row1': Math.floor(numPlots/2), '#row2': Math.ceil(numPlots/2)};
  for (var idx = 0; idx < numPlots; idx++) {
    if (data.children[idx].children.length == 0) {
      data.children[idx].children.push({});
    }
    // Puts two-lineage case in the same row to optimize space
    if (numPlots == 2) {
      var row = '#row1'; var width = 500;
    } else {
      var row = (idx+1)*2 <= numPlots ? '#row1' : '#row2';
      var width = 1000 / counts[row];
    }
    var chart = nv.models.sunburstChart().mode('size');
    var svg = d3.select(row)
          .append('svg');
          svg.datum([data.children[idx]])
          .attr('width', width)
          .attr('height', height)
          .attr('viewBox', '0 0 ' + width + ' ' + height)
          .attr('class', 'lineages')
          .call(chart);
    svg.append('text')
          .text(data.children[idx].name + " " + data.children[idx].value )
          .attr('x', '50%')
          .attr('y', '50%')
          .attr('class', 'lineage_text');
          // Selects/deselects the clicked drug
    chart.sunburst.dispatch.on("elementClick", function(e) {
        console.log("here")
         // toggleTabData('drug', e.data.x, e.data.x, 'map-marker');

          });

  }
}
