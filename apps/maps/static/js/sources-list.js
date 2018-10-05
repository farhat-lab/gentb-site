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
  var templates = {
      'source': $('#source_template'),
      'paper': $('#paper_template'),
      'bioproject': $('#bioproject_template'),
  }
  templates.source.hide();
  templates.paper.hide();
  templates.bioproject.hide();
  $("a", "#sources")
        .not(templates.source)
        .not(templates.paper)
        .not(templates.bioproject).remove();

  for(var i in data) {
      var datum = data[i];
      var template = templates[datum.kind];
      var copy = template.clone(true, true);
      $('#sources').append(copy);
      var sel_id = datum.kind + '_' + datum.pk;
      copy.attr('id', sel_id);
      $("h3", copy).text(datum.name);
      if(datum.uploader && datum.uploader != 'None') {
        $("p", copy).text("(by " + datum.uploader + ", " + datum.count + " records)");
      } else {
        $("p", copy).text("(" + datum.count + " records)");
      }
      copy.show();
      copy.data('id', datum.pk);
      copy.data('kind', datum.kind);
      copy.data('name', datum.name);
      copy.click(function() {
          // Select or deselect this source.
          var this_id = $(this).data('kind') + '_' + $(this).data('id');
          var that_id = getTabColumn('source') + '_' + getTabData('source');
          var that = $('#' + that_id);
          that.removeClass('btn-primary');
          that.addClass('btn-default');
          if(this_id == that_id) {
              unsetTabData('source');
          } else {
              $(this).addClass('btn-primary');
              $(this).removeClass('btn-default');
              setTabData('source', $(this).data('id'), $(this).data('name'), 'list', $(this).data('kind'))
          }
      });
  }
}
