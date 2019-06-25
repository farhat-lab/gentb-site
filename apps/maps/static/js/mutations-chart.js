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

var mutations = [];

set_tooltip = function(elem, title) {
    return elem.attr('data-original-title', title)
               .tooltip('fixTitle')
               .tooltip('show');
};

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

var mutations_chart;

$(document).ready(function() {
    $('body').tooltip({
        selector: 'input#snp'
    });

    var svg = 'svg.mutations';
    mutations_chart = initialiseMutationChart(svg);

    $('#mutation-store').data('url-signal', function(url, args) {

        // Setup mutations table first
        var t_table = $('#mutation_map table');
  
        if(t_table.data('loaded')) {
            return
        }   
        t_table.data('loaded', 1); 
  
        var table = t_table.DataTable({
          "processing": true,
          "serverSide": true,
          "ajax": {
              "url": url,
              "data": function ( data ) { 
                  // Sent json, store for future use in selecting
                  return data;
              },  
              "dataSrc": function ( json ) { 
                  // Returned json, filter and etc here.
                  return json.data;
              },  
          },  
          "language": {
            "processing": "Loading...",
          },  
          "columns": [
            {   
              "data": "name",
              "title": "Name",
              "description": "Name of the Gene Locus",
              // "render": $.fn.dataTable.render.number(',', '.', 3, ''),
            },  
            {   
              "data": "gene_locus__name",
              "title": "Gene Locus",
            },  
            {   
              "data": "nucleotide_position",
              "title": "Position",
            },  
            {   
              "data": "nucleotide_reference",
              "title": "Reference",
            },  
            {   
              "data": "nucleotide_varient",
              "title": "Varient",
            },
            {
              "data": "codon_reference",
              "title": "Codon",
            },
            {
              "data": "codon_varient",
              "title": "Varient",
            },
          ],
          'order': [[2, 'asc']],
        })

        table.on('draw', function(a, b, c) {
          $('tbody tr', t_table).click(function() {
            var data = table.row( this.rowIndex - 1 ).data();
            var snp = $('#snp');

            $(this).toggleClass('selected');
            toggleTabData('mutation', data.name, data.name, ' icon-mutation');

            mutations = getTabData('mutation');
            refresh_mutations_graph();

          }).each(function(key, value) {
            // Keep highlights as paging happens.
            var data = table.row( value.rowIndex - 1 ).data();
            if (data && mutations.includes(data.name)) {
              $(value).addClass('selected');
            }
          });
        });

      if($('#level-0').length == 0) {
        initialiseMutationList(url, args);
      }
    });
});

function refresh_mutations_graph() {
    $.getJSON($('svg.mutations').data('json-url'), getAllTabData())
      .done(function(json) {
          if(json.data.length > 0) {
            $('#mutations').show();
            $('#mutation_explainer').hide();
          } else {
            $('#mutations').hide();
            $('#mutation_explainer').show();
          }
          chartData('svg.mutations', mutations_chart, json.data);
      })
      .fail(function(jqxhr, textStatus, error) {
        var err = textStatus + ", " + error;
        console.log("Request Failed: " + err);
      });
}

function refreshMutation(svg, chart, data) {
  // The goal here is to list all mutations within the
  // selected drug, lineage or country
  $('#mselect').each(function() {
      // replaceOptions($('#level-0'), data.children);
      // $('#level-0').show();
      // $('#level-0').change();
  });
  $('#clear-mutation').click();


}

function initialiseMutationList(url, args, refresh_function) {
    var locus = $('#locus');
    var locus_datalist = $('#locus-list');
    var snp = $('#snp');
    var snp_datalist = $('#mutation-list');
    var button_del = $('#clear-mutation');

    var lookup_url = $('#gene_map').data('locus-lookup');
    var range_url = $('#gene_map').data('locus-range');


    var locus_select = function(val) {
      snp.val('');
      if(val && val != '---') {
        snp_datalist.empty();
        snp_datalist.data('set', '--- ASK AGAIN ---');
        snp.keyup();
        set_tooltip(locus, "Gene " + val + " Selected");
        // Refresh range
        reset_args();
        $('#mutation_selector').addClass('loading');
        $('#locus').attr("disabled", "disabled");
        $.getJSON(range_url, args).done(function(json) {
          $('#mutation_selector').removeClass('loading');
          $('#locus').removeAttr("disabled");
          $('#gene_map').show();
          $('#gene_start').text(json['start']);
          $('#gene_end').text(json['end']);
          // XXX $('#gene_label').text(json['title']);
          var max = parseFloat(json['max']);
          for(var i=1;i<=50;i++) {
              var items = json['values'][i - 1];
              if(items == undefined) { items = []; }
              var count = items.length;
              var height = (parseFloat(count) / max) * 50;

              // Back-calc the start of this range (for confirmation)
              var girth = (json['end'] - json['start']) / 50;
              var start = parseInt(((i - 1) * girth) + json['start']);
              var end = parseInt((i * girth) + json['start']);

              $('#ms-'+i).data('items', items)
                  .attr('data-original-title', count + " mutations, " + start + "-" + end)
                  .attr('data-container', 'body')
                  .tooltip();

              $('#ms-'+i+' span').css('line-height', height + 'px');
          }
        });
      }
    }

    $('#mutation_selector > span').click(function() {
        $('#mutation_selector > span.selected').removeClass('selected');
        $(this).addClass('selected');
        replaceOptions(snp_datalist, $(this).data('items'));
        snp.focus();
    });

    locus.select(function(){locus_select($(this).val())});
    locus.data('previous', locus.val());

    function reset_args() {
      args = getAllTabData();
      args.locus = locus.val();
      args.synonymous = $('input[name="synon"]').is(':checked');
      delete args.snp;
      delete args.ecoli;
      delete args.range;
    }

    locus.on('keyup', function(e){
      var selected = $(this).val();
      if($(this).data('previous') == selected) { return; }
      $(this).data('previous', selected);
      reset_args();
      
      $.getJSON(lookup_url, args).done(function(json) {
        set_tooltip(locus, json.msg || "Not updated");

        if(json.values) {
            replaceOptions(locus_datalist, json.values);
            if(json.values.length == 1 && json.values[0].toLowerCase() == selected.toLowerCase()) {
                locus.val(json.values[0]);
                locus_datalist.empty()
                locus_select(json.values[0]);
            } else if(selected.toLowerCase() == 'all') {
                locus.val('All');
                locus_datalist.empty()
                locus_select('All');
            }
        } else {
          locus_datalist.empty();
        }
        // Trigger workaround to update datalist.
        locus.focus();
      }).fail(function(jqxhr, textStatus, error) {
        var err = textStatus + ", " + error;
        console.log("Request Failed: " + err);
      });
    });

    snp.select(function() {
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

    snp.on('keyup', function(e){
      var selected = $(this).val();
      if($(this).data('previous') == selected || selected == '') { return; }
      $(this).data('previous', selected);
      reset_args();
      if(selected.toLowerCase().startsWith('e:')) {
          if(isNaN(selected.substr(2)) || selected.length <= 2) {
              return;
          }
          args.ecoli = selected.substr(2);
      } else {
          args.snp = selected;
      }
      $.getJSON(url, args).done(function(json) {
        snp.attr('data-original-title', json.msg || "Not updated")
          .tooltip('fixTitle')
          .tooltip('show');

        if(json.values) {
          replaceOptions(snp_datalist, json.values);
        } else {
          snp_datalist.empty();
        }
        // Trigger workaround to update datalist.
        snp.focus();
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
