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
  $('#source-store').data('json-signal', function(data) {
    console.log("Sources list!");
    listSources(data.values);
  }); 
});

function listSources(data) {
  var existing = getTabData("map");
  var template = $('#source_template');
  template.hide();
  $("a", "#sources").not(template).remove();

  for(var i in data) {
      var datum = data[i];
      var copy = template.clone(true, true);
      $('#sources').append(copy);
      copy.attr('id', 'source_' + datum[0]);
      $("h3", copy).text(datum[1]);
      $("p", copy).text("(by " + datum[2] + ", " + datum[3] + " records)");
      copy.show();
      copy.data('id', datum[0]);
      copy.data('name', datum[1]);
      copy.click(function() {
          // Select or deselect this source.
          var selected = getTabData('source');
          var selected_element = $('#source_' + selected);
          selected_element.removeClass('btn-primary');
          selected_element.addClass('btn-default');
          if(selected == $(this).data('id')) {
              unsetTabData('source');
          } else {
              $(this).addClass('btn-primary');
              $(this).removeClass('btn-default');
              setTabData('source', $(this).data('id'), $(this).data('name'), 'list')
          }
      });
  }
}
