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
    //can this just go somewhere else?
    var color = d3.scale.threshold()
        .domain([10, 20, 30, 40])
        .range(['#FFFFCC', '#C7E9B4', '#7FCDBB', '#41B6C4', '#1D91C0']);

    // const color = function color(maxValue) {
    //     var scaledColor = d3.scale.threshold()
    //     .domain([maxValue/4, maxValue/2, maxValue*.75, maxValue])
    //     .range(['#FFFFCC', '#C7E9B4', '#7FCDBB', '#41B6C4', '#1D91C0']);
    //     return scaledColor;
    // }


    initialiseStrainMap(map, color);
    $('#map-store').data('json-signal', function(data) {
      mapStrainData(map, color, data);
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
function initialiseStrainMap(map, color) {
    //get data here somehow?
//

  L.tileLayer('http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    maxZoom: 18,
    attribution: '' //Map tiles (c) <a href="http://openstreetmap.org">OpenStreetMap</a>, Map data (c) genTB'
  }).addTo(map);

  legend = L.control({
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

//TODO: look at editing this to change color?
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

// });
}
var map_layer = null;
function mapStrainData(map, color, data) {
    //needs a new color
    //Loop Through GeoJSON Properties
    //https://stackoverflow.com/questions/30133706/loop-through-geojson-properties
    // console.log(data)
    // const maxCount = function(){
    //      let biggest = 0;
    //      data.forEach(function() {
    //          console.log(feature)
    //     });
    // }
    // maxCount()
    //     for (let i=0; i < dataArray.length; i++){
    //         if(d3.max(dataArray[i], (d)=>+d.Count) > biggest){
    //             biggest = d3.max(dataArray[i], (d)=>+d.Count)
    //         }
    //     }
    //     return biggest
    // }



  var existing = getTabData('map');

  if(map_layer) {
    // Remove previous layer
    map.removeLayer(map_layer);
  }

  function get_style(feature) {

      var color_max;
      var color_by;
      if ($(".form-check-input:checked").val() == "feature.properties[1][0].gdp") {
          color_by = feature.properties[1][0].gdp;
          color_max = 400000;
      }
      else {
          color_by = feature.properties[0][0].values.Total;
          color_max = 1000;
      }
  var best_color = d3.scale.threshold()
        .domain([color_max/10, color_max/2, color_max*3/4, color_max])
        .range(['#FFFFCC', '#C7E9B4', '#7FCDBB', '#41B6C4', '#1D91C0']);
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
      var max;
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

    if ($(".form-check-input:checked").val() == "feature.properties[1][0].gdp") {
        var legend_label = "GDP"
        max = 400000;
    }
    else {
        var legend_label = "Number of isolates"
        max = 500;
    }

    newLegend(max);

    d3.select('.caption')
      .text(legend_label);


});

function newLegend(max) {

    var color = d3.scale.threshold()
        .domain([max/4, max/2, max*.75, max])
        .range(['#FFFFCC', '#C7E9B4', '#7FCDBB', '#41B6C4', '#1D91C0']);

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

  //TODO: look at editing this to change color?
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

  function onEachFeature(feature, layer) {
    var country_code = feature.properties[0][0].value;
    ret = $('<div></div>');
    ret.append($('<h4>' + feature.properties[0][0].name + '</h4>'));
    $.each(['Sensitive', 'MDR', 'XDR', 'TDR'], function(key, value) {
      if (feature.properties[0][0].values[value]) {
        ret.append($('<div>Number of <span>' + value + ' isolates: </span><span>' + feature.properties[0][0].values[value] + "</span></div>"));
      }
    });
    if(Object.keys(feature.properties[0][0].values).length > 2) {
      ret.append($("<hr/><div class='total'><span>Total isolates: </span><span>" + feature.properties[0][0].values.Total + "</span></div>"));
    }
    ret.append($('<h6>' + 'Social Determinants of Health:' + '</h6>'));
    ret.append($("<hr/><div> GDP: <span>" + feature.properties[1][0].gdp + "</span></div>"));
    var button1 = $("<button class='btn btn-primary btn-xs'>Select Country</button>").click(function() {
        // WARNING: jquery class selectors and addClass/removeClass DO NOT work here

        var next = $('.country-'+country_code);
        // Set element id just in case it's useful later
        next[0].id = country_code;

        // Add the countrySelect class to highlight it
        next.attr('class', next.attr('class') + ' countrySelect');
        next.data('deselect', function() { button2.click(); });

        // Set the usable data for other charts
        addTabData('map', country_code, feature.properties[0][0].name, 'flag');

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
