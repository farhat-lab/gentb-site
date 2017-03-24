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
(function($) {
  $.fn.addOption = function(name, value) {
    var opt = $("<option></option>");
    $(this).append(opt);
    opt.attr("value", value);
    opt.text(name);
    return opt;
  }
  $.fn.replaceOptions = function(options) {
    this.empty();
    var self = this;
    if($(self).is('select')) {
      $(self).addOption("---", '---');
    }
    $.each(options, function(index, option) {
      var opt = $(self).addOption(option.name, option.value);
      opt.data('children', option.children);
    }); 
  }

  $(document).ready(function() {
    var svg = $('svg.mutations');
    makeMutationChart();
    $('#mutation-store').click(mutationRefresh());
  });

function mutationRefresh() {
  // The goal here is to list all mutations within the
  // selected drug, lineage or country
  $('.lister').each(function() {
    var last_id = null;
    var previous = null;
    var target = $(this);
    var container = $(this).parent();
    $.getJSON($(this).data('mutations'), function(data) {
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
            select.attr('style', 'width: 50%;');
            select.attr('id', 'level-' + index);
            select.insertBefore(target);
            select.replaceOptions([]);

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
                      list.replaceOptions(children);
                  } else {
                      next.replaceOptions(children);
                  }
                }
              }
              next.change();
              next.select();
            });
        }
      });

      var add_button = $('<a class="btn btn-success btn-sm">Add</button>');
      add_button.insertBefore(target);
      add_button.attr('id', 'level-button');
      add_button.attr('style', 'width: 45;');
      add_button.click(function() {
        var value = $('#level-2').val();
        if(value && value != '---') {
          $(target).val(value + '\n' + $(target).val());
          $('#level-2').val('').select();
        }
      });
      $('#level-0').replaceOptions(data.children);
      $('#level-0').show();
      $('#level-0').change();
    });
  });
}

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
     d3.select('svg.mutations').datum(data).call(chart);
     x++;
    }, 3000);

    nv.utils.windowResize(chart.update);
    ready = true;
}

})(jQuery);
