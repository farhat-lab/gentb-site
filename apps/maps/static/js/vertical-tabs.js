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

var all_tabs = 'div.vertical-tab-menu > div.list-group > a';

$(document).ready(function() {
  // DEBUG
  $.ajaxSetup({ cache: false });
  $(all_tabs).click(function(e) {
    var tab = $(this);
    e.preventDefault();
    tab.siblings('a.active').removeClass("active");
    tab.addClass("active");
    var index = tab.index();
    $("div.vertical-tab>div.vertical-tab-content").removeClass("active");
    $("div.vertical-tab>div.vertical-tab-content").eq(index).addClass("active");

    var url = tab.data('json-url');
    if(url && !tab.data('done')) {
      var data = getAllTabData(this.id);

      $.getJSON(url, data).done(function(json) {
          tab.data('json-signal')(json, url, data);
          tab.data('done', true);
        })
        .fail(function(jqxhr, textStatus, error) {
          var err = textStatus + ", " + error;
          console.log( "Request Failed: " + err );
        });
    }
  });
  // Activate the existing active tab.
  $(all_tabs+'.active').click();
});

/* Collects all tab values into a dictionary ready for sending to the server */
function getAllTabData(except) {
  var data = {};
  $(all_tabs).each(function() {
    if(this.id != except) {
      if($(this).data('value')) {
        var column = this.id.replace('-store', '');
        if($(this).data('column')) {
          column = $(this).data('column')
        }
        data[column] = $(this).data('value');
      }
    }
  });
  return data;
}

/* Returns data for the given tab key name */
function getTabData(key) {
  var store = $('#'+key+'-store');
  return store.data('value');
}

/*
   key    - This vertical tab that this data is filed under, should match html id.
   value  - The value that should be sent to the server
   text   - The new text name for this tab while selected
   icon   - The new icon for this tab while selected
   column - When the selection can result in different types of values
            the tab can optionally set a column name which replaces
            'key' as the key word argument name sent to the server.
 */
function setTabData(key, value, text, icon, column) {
  var store = $('#'+key+'-store');

  if(!store.data('original-text')) {
    store.data('original', store.data('value'));
    store.data('original-column', store.data('column'));
    store.data('original-text', $('p', store).text());
    store.data('original-icon', $('h2', store).attr('class'));
  }
  if(value != store.data('value')) {
    store.addClass('selected');
    store.data('value', value)
    store.data('column', column)
    $('p', store).text(text);
    $('h2', store).attr('class', 'glyphicon glyphicon-'+icon);

    // Clear all existing graph renderings
    $(all_tabs).not(store).data('done', false);
  }
}

/* Deselect tab, removing it's value and putting all values back to what they
   where when the page loaded */
function unsetTabData(key) {
  var store = $('#'+key+'-store');
  store.removeClass('selected');
  $(all_tabs).not(store).data('done', false);
  store.removeData('value');
  store.removeData('column');

  if(store.data('original-text')) {
    store.data('value', store.data('original'));
    store.data('column', store.data('original-column'));
    $('p', store).text(store.data('original-text'));
    $('h2', store).attr('class', store.data('original-icon'));
  }
}

/* d3 function for adding data to a d3 svg chart */
function chartData(svg, chart, data) {
  return d3.select(svg)
    .datum(data)
    .call(chart);
}

