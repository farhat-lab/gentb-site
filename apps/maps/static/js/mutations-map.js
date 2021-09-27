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
  $('div.maps').each(function(map) {
//this needs to get the data somehow so we can initialize the legend and scale it dynamically
    var map = L.map(this.id).setView([12, 25], 2);
    var newLegend = function(map, max, color) {

            legend = L.control({
              position: 'topright'
          });

            legend.onAdd = function(map) {
              var div = L.DomUtil.create('div', 'legend');
              return div
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

        }
    initialiseStrainMap(map);
    $('#map-store').data('json-signal', function(data) {
    console.log(data)
      var getMaxAttribute = function (data, attribute) {
          var max = -10000000;
          var maxPlace;
          var roundPlace;
          for (var i=0 ; i<data.length ; i++) {
              if (attribute === "total_strains"){
                  max = Math.max(data[i]["properties"].values.Total, max);
              }
              else if (attribute === "MDR") {
                  if (!(isNaN(data[i]["properties"].values.MDR))){
                      max = Math.max(data[i]["properties"].values.MDR, max);
                  }


              }
              else{
                  max = Math.max(data[i]["properties"][attribute], max);
              }
          }
          maxPlace = Math.ceil(Math.log10(max));
          try {
              return max.toPrecision(maxPlace);
          } catch {
              return -1;
          }
      }
      var maxGDP = getMaxAttribute(data.features, "world_bank_gdp")
      var maxTotal = getMaxAttribute(data.features, "total_strains");
      var maxWealth = getMaxAttribute(data.features, "total_wealth");
      var maxFunding = getMaxAttribute(data.features, "total_funding");
      var maxPop = getMaxAttribute(data.features, "pop_dens");
      var maxTB = getMaxAttribute(data.features, "all_tb_incidence2018");
      var maxMDR = getMaxAttribute(data.features, "MDR");
      var maxWHOMDR = getMaxAttribute(data.features, "who_est_mdr");

      var color_by_dict = {
          // id : ["label", max_value, "key", function() { human readable value }]
          "feature.properties.world_bank_gdp": ["GDP (Trillions $)", maxGDP, "world_bank_gdp", function (val) { return Math.round(val * 100) / 100; } ],
          "feature.properties.values.Total": ["Total isolates", maxTotal, "Total", function (val) { return undefined; }],
          "feature.properties.total_wealth": ["Total Wealth per Capita ($)", maxWealth, "total_wealth", function (val) { return Math.round(val); }],
          "feature.properties.total_funding": ["Total TB Funding (Billions $)", maxFunding, "total_funding", function (val) {
              if (val && val.length != 0) {
                  return Math.round((val/10000000)*10000000)/1000000000;
              }
              return undefined;
          }],
          "feature.properties.pop_dens": ["Population density (people per sq. km of land area)", maxPop, "pop_dens", function (val) { return Math.round(val*10)/10; }],
          "feature.properties.all_tb_incidence2018": ["TB Incidence (all types, cases per 100,000)", maxTB, "all_tb_incidence2018", function (val) { return val; }],
          "feature.properties.values.MDR": ["MDR isolates", maxMDR, "MDR", function (val) { return undefined; }],
          "feature.properties.who_est_mdr": ["WHO Estimated MDR (% of new cases)", maxWHOMDR, "who_est_mdr", function (val) { return (val*100000)/100000; }]
      }

      mapStrainData(map, data, newLegend, maxGDP, maxTotal, maxWealth, color_by_dict);
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
function mapStrainData(map, data, newLegend, maxGDP, maxTotal, maxWealth, color_by_dict) {
  let color_max;
  let best_color;

  var existing = getTabData('map');

  if(map_layer) {
    // Remove previous layer
    map.removeLayer(map_layer);
  }


  function get_style(feature) {
    var color_by;
    if ($(".form-check-input:checked").val() == "feature.properties.values.Total"){
        color_by = feature.properties.values.Total
    }
    else if ($(".form-check-input:checked").val() == "feature.properties.values.MDR"){
        color_by = feature.properties.values.MDR
    }
    else {
        color_by = feature.properties[color_by_dict[$(".form-check-input:checked").val()][2]]
    }

    color_max = color_by_dict[$(".form-check-input:checked").val()][1]
    best_color = d3.scale.threshold()
        .domain([color_max/4, color_max/2, color_max*3/4, color_max])
        .range(['#FFFFCC', '#C7E9B4', '#7FCDBB', '#41B6C4', '#1D91C0']);
    if (!legend) {
        newLegend(map, maxTotal, best_color)
    }


    return {
     fillColor: best_color(color_by),
      weight: 1,
      opacity: 0.5,
      color: 'black',
      fillOpacity: 0.35
    };
  }

  //toggle mapping color
  $('.form-check-input').on('click',function(e){
      if(map_layer) {
        // Remove previous layer
        map.removeLayer(map_layer);
      }
      map_layer = L.geoJson(data,  {
        style: get_style,
        onEachFeature: onEachFeature
      }
    ).addTo(map)

    if(legend) {
      // Remove previous layer
      map.removeControl(legend);
    }

    newLegend(map, color_by_dict[$(".form-check-input:checked").val()][1], best_color);

    d3.select('.caption')
      .text(color_by_dict[$(".form-check-input:checked").val()][0]);


  });

  function onEachFeature(feature, layer) {
    var country_code = feature.properties.value;
    ret = $('<div></div>');
    ret.append($('<h4 style = "font-weight:bold">' + feature.properties.name + '</h4>'));
    $.each(['Sensitive', 'MDR', 'XDR', 'TDR'], function(key, value) {
      if (feature.properties.values[value]) {
        ret.append($('<div style="color:black">Number of <span>' + value + ' isolates: </span><span>' + feature.properties.values[value] + "</span></div>"));
      }
    });
    if(Object.keys(feature.properties.values).length > 2) {
      ret.append($("<hr/><div class='total' style='color:black'><span>Total isolates: </span><span>" + feature.properties.values.Total + "</span></div>"));
    }
    var opt_id = $(".form-check-input:checked").val();
    var opt_col = color_by_dict[opt_id];

    var opt_val = feature.properties.values[opt_col[2]];
    if (opt_val == undefined) {
        opt_val = feature.properties[opt_col[2]];
    }
    opt_val = opt_col[3](opt_val);
    if (opt_val != undefined) {
        ret.append($("<hr/><h6><div style='color:#000080'> " + opt_col[0] + ": <span>" + opt_val + "</span></div></h6>"));
    }

    var button1 = $("<button class='btn btn-primary btn-xs'>Select Country</button>").click(function() {
        // WARNING: jquery class selectors and addClass/removeClass DO NOT work here

        var next = $('.country-'+country_code);
        // Set element id just in case it's useful later
        next[0].id = country_code;

        // Add the countrySelect class to highlight it
        next.attr('class', next.attr('class') + ' countrySelect');
        next.data('deselect', function() { button2.click(); });

        // Set the usable data for other charts
        addTabData('map', country_code, feature.properties.name, 'flag');

        button1.hide();
        button2.show();
        map.closePopup();
    });
    button1.attr('id', 'select-'+country_code);
    var button2 = $("<button class='btn btn-danger btn-xs'>Deselect</button>").click(function() {
        // WARNING: jquery class selectors and addClass/removeClass DO NOT work here
        var previous = $('.country-'+country_code);
        previous.attr('class', previous.attr('class').replace(' countrySelect', ''));
        removeTabData('map', country_code);
        button1.show();
        button2.hide();
        map.closePopup();
    });
    ret.append($("<hr/>"));
    ret.append(button1);
    ret.append(button2);
    layer.bindPopup(ret[0]);

    var id_cls = 'country-' + country_code;
    if(existing.includes(country_code)) {
        id_cls += ' countrySelect';
        button1.hide();
    } else {
        button2.hide();
    }
    layer.setStyle({className: id_cls});

  }

  map_layer = L.geoJson(data,  {
    style: get_style,
    onEachFeature: onEachFeature
  }
).addTo(map)




}
