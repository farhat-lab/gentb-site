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
  $('div.maps').each(function() {
    var map = $(this);
    $.getJSON(map.data('map'), function(data) {
      makeMap(map, data);
    });
  });
});

function matchKey(datapoint, key_variable) {
  if(key_variable[0]) {
      return (parseFloat(key_variable[0][datapoint]));
  }
  return "black";
}

function makeMap(target, gjson_1) {

  var color = d3.scale.threshold()
    .domain([10, 20, 30, 40])
    .range(['#FFFFCC', '#C7E9B4', '#7FCDBB', '#41B6C4', '#1D91C0']);

  var map = L.map(target[0].id).setView([12, 25], 2);

  L.tileLayer('http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    maxZoom: 18,
    attribution: '' //Map tiles (c) <a href="http://openstreetmap.org">OpenStreetMap</a>, Map data (c) genTB'
  }).addTo(map);

  function style_1(feature) {
    return {
      fillColor: color(feature.properties.values.Total),
      weight: 1,
      opacity: 0.5,
      color: 'black',
      fillOpacity: 0.35
    };
  }

  function onEachFeature(feature, layer) {
    ret = '<h4>' + feature.properties.name + '</h4>';
    $.each(['Sensitive', 'MDR', 'XDR', 'TDR'], function(key, value) {
      if (feature.properties.values[value]) {
        ret += '<div><span>' + value + ':</span><span>' + feature.properties.values[value] + "</span></div>";
      }
    });
    if(Object.keys(feature.properties.values).length > 2) {
      ret += "<hr/><div class='total'><span>Total: </span><span>" + feature.properties.values.Total + "</span></div>";
    }
    ret += "<hr/><button class='btn btn-primary btn-xs disabled'>Select Country</button>"
    layer.bindPopup(ret);
  }
  gJson_layer_1 = L.geoJson(gjson_1, {
    style: style_1,
    onEachFeature: onEachFeature
  }).addTo(map)

  var legend = L.control({
    position: 'topright'
  });

  legend.onAdd = function(map) {
    var div = L.DomUtil.create('div', 'legend');
    return div
  };
  legend.addTo(map);

  var x = d3.scale.linear()
    .domain([0, 40])
    .range([0, 400]);

  var xAxis = d3.svg.axis()
    .scale(x)
    .orient("top")
    .tickSize(1)
    .tickValues(color.domain())

  var svg = d3.select(".legend.leaflet-control").append("svg")
    .attr("id", 'legend')
    .attr("width", 450)
    .attr("height", 40);

  var g = svg.append("g")
    .attr("class", "key")
    .attr("transform", "translate(25,16)");

  g.selectAll("rect")
    .data(color.range().map(function(d, i) {
      return {
        x0: i ? x(color.domain()[i - 1]) : x.range()[0],
        x1: i < color.domain().length ? x(color.domain()[i]) : x.range()[1],
        z: d
      };
    }))
    .enter().append("rect")
    .attr("height", 10)
    .attr("x", function(d) {
      return d.x0;
    })
    .attr("width", function(d) {
      return d.x1 - d.x0;
    })
    .style("fill", function(d) {
      return d.z;
    });

  g.call(xAxis).append("text")
    .attr("class", "caption")
    .attr("y", 21)
    .text('Cases');

};
