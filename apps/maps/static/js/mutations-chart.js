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

function addOption(self, name, value) {
    var opt = $("<option></option>");
    $(self).append(opt);
    opt.attr("value", name.value || value || name);
    opt.text(name.name || name);
    return opt;
}
function replaceOptions(self, options) {
    self.empty();
    if($(self).is('select')) {
      addOption($(self), "---");
    }
    $.each(options, function(index, option) {
      addOption($(self), option).data('children', option.children);
    }); 
}
/*
 * This format function allows us to format our tooltips using a basic
 * templating. Basically each tooltip will be a full text field with html
 * and the values that go into the various parts of the text will be in the
 * format {field_name:format} and will come from the point data source.
 */
String.prototype.format = function (values) {
  return this.replace(/{(([^}:]+):)?([^}]*?)?}/g, function () {
    var name = (arguments[2] || "value").split('.');
    var format = arguments[3] || "s";

    var datum = values;
    for(var i = 0; i < name.length; i++) {
      datum = datum[name[i]];
    }
    var value = datum;
    if(format != 's') {
      value = d3.format(format)(datum);
    }
    if(value == 'NaN' || typeof value == "NaN") {
      value = datum;
    }   
    return typeof value != 'undefined' ? value : ''; 
  }); 
};
function unique(array) {
    return $.grep(array, function(el, index) {
        return index === $.inArray(el, array);
    });
}

$(document).ready(function() {
    $('body').tooltip({
        selector: 'input#snp'
    });

    var svg = 'svg.mutations';
    var chart = initialiseMutationChart(svg);
    $('#mutation-store').data('json-signal', function(data, url, args) {
      if($('#level-0').length == 0) {
        initialiseMutationList(data, url, args, function(mutations) {
        $.getJSON($(svg).data('json-url'), getAllTabData())
          .done(function(json) {
              if(json.data.length > 0) {
                $('#mutations').show();
                $('#mutation_explainer').hide();
              } else {
                $('#mutations').hide();
                $('#mutation_explainer').show();
              }
              chartData(svg, chart, json.data);
          })
          .fail(function(jqxhr, textStatus, error) {
            var err = textStatus + ", " + error;
            console.log("Request Failed: " + err);
          });
        });
      }
      refreshMutation(svg, chart, data);
    }); 
});

function refreshMutation(svg, chart, data) {
  // The goal here is to list all mutations within the
  // selected drug, lineage or country
  $('#mselect').each(function() {
      //replaceOptions($('#level-0'), data.children);
      //$('#level-0').show();
      //$('#level-0').change();
  });
  $('#clear-mutation').click();
}

function initialiseMutationList(data, url, args, refresh_function) {
  $('#mselect').each(function() {
    var container = $(this);

    var select = $('<select id="locus"></select>');
    var input = $('<input type="text" list="mutation-list" id="snp" style="width: 300px;" data-container="body" autocomplete="off" data-toggle="tooltip" title="Select a locus to continue"/>');
    var datalist = $('<datalist id="mutation-list"></datalist>');
    var button_del = $('<a class="btn btn-danger btn-sm pull-right" id="clear-mutation">Clear</button>');
    $('#mutations').hide();

    container.empty();
    container.append(select);
    container.append(input);
    container.append(datalist);
    container.append(button_del);

    // Add this data to the initial locus list
    replaceOptions(select, data.values);

    select.change(function() {
      var val = $(this).val();
      input.val('');
      if(!val || val == '---') {
        input.prop('disabled', true);
      } else {
        input.prop('disabled', false);
        datalist.empty();
        datalist.data('set', '--- ASK AGAIN ---');
        input.keyup();
      }
    }).change();

    input.select(function() {
      var value = $(this).val();
      var mutations = $('#mutation-store').data('value');
      mutations.push(value);
      mutations = unique(mutations);
      $('#mutation-store').data('value', mutations);
      refresh_function(mutations);
      $(this).val('');
      $(this).blur();
      button_del.show();
    });

    input.on('keyup', function(e){
      var selected = $(this).val();
      if($(this).data('previous') == selected) { return; }
      $(this).data('previous', selected);
      args.locus = select.val();
      delete args.snp;
      delete args.ecoli;
      if(selected.toLowerCase().startsWith('e:')) {
          if(isNaN(selected.substr(2)) || selected.length <= 2) {
              return;
          }
          args.ecoli = selected.substr(2);
      } else {
          args.snp = selected;
      }
      $.getJSON(url, args).done(function(json) {
        input.tooltip('hide')
          .attr('data-original-title', json.msg || "Not updated")
          .tooltip('fixTitle')
          .tooltip('show');

        if(json.values) {
          replaceOptions(datalist, json.values);
        } else {
          datalist.empty();
        }
        // Trigger workaround to update datalist.
        input.focus();
      }).fail(function(jqxhr, textStatus, error) {
        var err = textStatus + ", " + error;
        console.log("Request Failed: " + err);
      });
    });

    button_del.click(function() {
      $('#mutation-store').data('value', new Array());
      refresh_function(new Array());
      button_del.hide();
    });

  });
}

function initialiseMutationChart() {
    var chart = nv.models.multiBarChart()
      .stacked(false)
      .reduceXTicks(false);

    var width = 1000;
    var height = 600;

    chart.margin({top: 40, right: 20, bottom: 120, left: 80});
    chart.height(height);
    chart.width(width);
    chart.yAxis.scale(100).orient("left")
        .axisLabel('Percent of strains')
        .tickFormat(d3.format(".2%"));

    chart.tooltip.contentGenerator(function (data) {
	data.extra = {
	   'head': chart.tooltip.headerFormatter()(data.value),
	   'y': chart.tooltip.valueFormatter()(data.data.y),
	};
	template = '<table><thead>' +
	  '<tr><th colspan="3">{extra.head:s}</th></tr>' +
	  '<tr><td class="legend-color-guide">' +
	    '<div style="background-color:{color:s}"></div></td>' +
	    '<td><strong>{data.key:s}</strong></td>' +
	  '<td>{extra.y:s}</td><td>{data.value:d} of {data.total:d}</td></tr></thead></table>';
        return template.format(data);
    });


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

    nv.utils.windowResize(chart.update);
    return chart;
}
