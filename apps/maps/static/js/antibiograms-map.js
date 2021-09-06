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
      createFilterDropdowns(data.c_filters);
    });
  });
});

// Create dropdowns which allow the user to filter the data shown in the maps.
function createFilterDropdowns(filters) {

  var template = $('.toolbar .template').clone();
  template.removeClass('template');

  $.each(filters, function(index, f) {
      var existing = $('.toolbar #' + f.key);

      if (existing.length) {
          console.log("Existing!");
      } else {
          console.log("New dropdown!", f);

          var dropdown = template.clone();

          $('.filter-label', dropdown).text(f.label);
          $('.dropdown-menu li', dropdown).remove();
          dropdown[0].id = f.key;

          $.each(f.values, function(index, v) {
              $('.dropdown-menu', dropdown).append($('<li data-value="' + v.value + '"><a href="#">' + v.label + '</a></li>'));
              if (v.selected) {
                  $('.filter-value', dropdown).text(v.value);
              }
          });
          dropdown.show();
          $('.toolbar').append(dropdown);

          // We need to make a new tab store here.
          var store = $('<a href="#" class="list-group-item text-center filter-' + f.key + '" id="' + f.key + '-store" style="display: none;"><h2 class=""></h2><p></p></a>');
          $('#data-store').append(store);

      }
  });
  $('.toolbar .template').remove();

  // Filters toolbar
  $('.toolbar .dropdown-menu li').click(function() {
      var parent = this.parentNode.parentNode;
      var key = parent.id;
      var value = $(this).data('value');
      $('.filter-value', parent).text(value);
      unsetTabData(key);
      addTabData(key, value, value, undefined, undefined);
      reloadMapData();
  });
}

function reloadMapData() {
    var map_tab = $('#map-store');
    map_tab.trigger("data:refresh");
}

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

var details = null;
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
        var conf = data['fill'];
        var color_by = feature.properties.row[conf['column']];
        var color_max = 1.0;
        best_color = d3.scale.threshold()
            .domain(conf['ranges'])
            .range(conf['colors']);
        return {
            fillColor: best_color(color_by),
            weight: 1,
            opacity: 0.5,
            color: 'black',
            fillOpacity: 0.35
        };
    }
    details = data["details"];
    map_layer = L.geoJson(data, {
      style: get_style,
      onEachFeature: onEachFeature,
    }).addTo(map)
}

function type_s(v) {
    return v;
}
function type_f(v) {
    return v.toLocaleString(undefined, {style: 'percent', minimumFractionDigits: 2})
}
function type_i(v) {
    return v; // TODO: ParseInt
}

function onEachFeature(feature, layer) {
    var country_code = feature.properties.value;
    var row = feature.properties.row;
    var drug = feature.properties.drug;

    ret = $('<div class="map-hover"></div>');

    ret.append($('<h4>' + feature.properties.name + '</h4>'));
    ret.append($("<div> Drug: <span><strong>" + drug.code + "</strong> (" + drug.name + ")</span></div>"));
    $.each(details, function(i, e) {
        var type = type_s;
        if (e.type == 'float') { type = type_f; }
        if (e.type == 'int') { type = type_i; }
        var dat = type(row[e.column]);
        ret.append($("<div> " + e.label + ": <span>" + dat + "</span></div>"));
    });

    layer.bindPopup(ret[0]);
}
