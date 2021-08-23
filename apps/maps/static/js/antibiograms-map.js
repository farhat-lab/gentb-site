/*
 * Copyright 2021, Maha Farhat
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
  $('div.maps').each(function(map) {
    var map = L.map(this.id).setView([12, 25], 2);
    var newLegend = function(map, max, color) {
        legend = L.control({
            position: 'topright'
        });

        legend.onAdd = function(map) {
            return L.DomUtil.create('div', 'legend');
        };
        legend.addTo(map);

        var x = d3.scale.linear()
              .domain([0, max])
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
              .text('Number of isolates');

    };

    initialiseStrainMap(map);

    $('#map-store').data('json-signal', function(data) {
      mapAntibiogramData(map, data);
      map.invalidateSize();
    });
  });
});

function matchKey(datapoint, key_variable) {
  if(key_variable[0]) {
      return (parseFloat(key_variable[0][datapoint]));
  }
  return "black";
}

var legend = null;
function initialiseStrainMap(map) {
  L.tileLayer('http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    maxZoom: 18,
    attribution: '' //Map tiles (c) <a href="http://openstreetmap.org">OpenStreetMap</a>, Map data (c) genTB'
  }).addTo(map);
}

var map_layer = null;
function mapAntibiogramData(map, data) {
    let color_max;
    let best_color;

    var existing = getTabData('map');
    if(map_layer) {
        // Remove previous layer
        map.removeLayer(map_layer);
    }
    function get_style(feature) {
        var color_by = feature.properties.row.mean_snp10;
        var color_max = 1.0;
        best_color = d3.scale.threshold()
            .domain([color_max/8, color_max/4, color_max/2, color_max*3/4])
            .range(['#FFFFCC', '#C7E9B4', '#7FCDBB', '#41B6C4', '#1D91C0']);
        return {
            fillColor: best_color(color_by),
            weight: 1,
            opacity: 0.5,
            color: 'black',
            fillOpacity: 0.35
        };
    }
    map_layer = L.geoJson(data, {
      style: get_style,
      onEachFeature: onEachFeature,
    }).addTo(map)
}

function f(v) {
    return v.toLocaleString(undefined, {style: 'percent', minimumFractionDigits: 2})
}

function onEachFeature(feature, layer) {
    var country_code = feature.properties.value;
    var row = feature.properties.row;
    var drug = feature.properties.drug;

    ret = $('<div class="map-hover"></div>');

    ret.append($('<h4>' + feature.properties.name + '</h4>'));
    ret.append($("<div> Drug: <span><strong>" + drug.code + "</strong> (" + drug.name + ")</span></div>"));
    ret.append($("<div> Number of Isolates: <span>" + (row.gentb_snp10_n) + "</span></div>"));
    ret.append($("<div> Marginal Resistance Rate: <span>" + f(row.mean_snp10) + "</span></div>"));
    ret.append($("<div> Lower Marginal Resistance Rate: <span>" + f(row.hi_snp10) + "</span></div>"));
    ret.append($("<div> Upper Marginal Resistance Rate: <span>" + f(row.lo_snp10) + "</span></div>"));

    layer.bindPopup(ret[0]);
}
