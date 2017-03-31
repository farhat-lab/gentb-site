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
    opt.attr("value", value);
    opt.text(name);
    return opt;
}
function replaceOptions(self, options) {
    self.empty();
    if($(self).is('select')) {
      addOption($(self), "---", '---');
    }
    $.each(options, function(index, option) {
      var opt = addOption($(self), option.name, option.value);
      opt.data('children', option.children);
    }); 
}
function unique(array) {
    return $.grep(array, function(el, index) {
        return index === $.inArray(el, array);
    });
}

$(document).ready(function() {
    var svg = 'svg.mutations';
    var chart = initialiseMutationChart(svg);
    $('#mutation-store').data('json-signal', function(data) {
      if($('#level-0').length == 0) {
        initialiseMutationList(data, function(mutations) {
          
        $.getJSON($(svg).data('json-url'), {'mutations': mutations})
          .done(function(json) {
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
      replaceOptions($('#level-0'), data.children);
      $('#level-0').show();
      $('#level-0').change();
  });
  $('#clear-mutation').click();
}

function initialiseMutationList(data, refreshNow) {
  $('#mselect').each(function() {
    var last_id = null;
    var previous = null;
    var target = $(this);
    var container = $(this).parent();
    target.data('data', data);
    $.each(data.levels, function(index, level) {
        if(index == 2) {
            var input = $('<input type="text"></input>');
            input.attr('id', "level-"+index);
            input.attr('list', "list-"+index);
            input.attr('title', level);
            input.data('next', 'button');
            input.attr('style', 'width: calc(100% - 45px);');
            input.insertBefore(target);

            var select = $('<datalist></datalist>');
            select.attr('id', "list-"+index);
            select.insertBefore(target);

            input.select(function() {
              var value = $(this).val();
              if(value && value != '---') {
                $('#level-button').removeClass('disabled');
              } else {
                $('#level-button').addClass('disabled');
              }
            });
        } else {
            var select = $('<select></select>');
            select.data('previous', index - 1);
            select.data('next', index + 1);
            select.attr('title', level);
            select.attr('id', 'level-' + index);
            select.insertBefore(target);
            replaceOptions(select, []);

            select.change(function() {
              var selected = this.selectedOptions;
              var next_id = $(this).data('next');
              var next = $('#level-' + next_id);
              if(!selected || selected[0].label == '---') {
                next.prop("disabled", true);
                if(next.is('select')) {
                  next.val('---');
                } else {
                  next.val('');
                }
              } else {
                next.prop("disabled", false);
                var children = $(selected[0]).data('children');
                var list = $('#list-' + next_id);
                if(children) {
                  if(list.length == 1) {
                      replaceOptions(list, children);
                  } else {
                      replaceOptions(next, children);
                  }
                }
              }
              next.change();
              next.select();
            });
        }
      });


      var add_button = $('<a class="btn btn-success btn-sm" id="add-mutation">Add</button>');
      add_button.insertBefore(target);
      add_button.click(function() {
        var value = $('#level-1').val();
        var mutations = $('#mutation-store').data('value');
        mutations.push(value);
        mutations = unique(mutations);
        $('#mutation-store').data('value', mutations);
        refreshNow(mutations);
      });
      var clear_button = $('<a class="btn btn-danger btn-sm" id="clear-mutation">Clear</button>');
      clear_button.insertBefore(target);
      clear_button.click(function() {
          $('#mutation-store').data('value', new Array());
          refreshNow(new Array());
      }).click();
  });
}

function initialiseMutationChart() {
    var chart = nv.models.multiBarChart()
      .stacked(true)
      .reduceXTicks(false);

    var width = 1000;
    var height = 600;

    chart.margin({top: 40, right: 20, bottom: 120, left: 80});
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

    nv.utils.windowResize(chart.update);
    return chart;
}
