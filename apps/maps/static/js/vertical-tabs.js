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

//Uses jquery selectors to select all the tabs (created in the beginning of map.html)
var all_tabs = 'div.vertical-tab-menu > div.list-group > a';

$(document).ready(function() {
  // DEBUG

  //Setup the ajax agent that we will use to get data from the server
  $.ajaxSetup({ cache: false });

  //whenever you click on any of the tabs
  $(all_tabs).click(function(e) {
    //we get the current tab that was actually clicked
    var tab = $(this);

    //We don't want to change a tab, we just wanna play with some data so capture the event/prevent default
    e.preventDefault();

    //Lets find the tab that used to be active and make it not active
    tab.siblings('a.active').removeClass("active");

    //Lets add active to the tab that you just clicked
    tab.addClass("active");

    //Gets the index of the tab element (0 indexed of course)
    var index = tab.index();

    //Finds that blue introduction tab and does the same active work with that in accordance with the current tab you picked
    $("div.vertical-tab>div.vertical-tab-content").removeClass("active");
    $("div.vertical-tab>div.vertical-tab-content").eq(index).addClass("active");

    //Access the json-url data currently stored in the tab
    var url = tab.data('json-url');

    // Initializes set of selectors (e.g. drugs, countries)
    var store = $('#'+this.id);
    if (!store.data('values')) {
      store.data('values', []);
      store.data('map', {});
    }


    //If there is actually data and done is false ('meaning that we didn't already do this whole process)
    if(url && !tab.data('done')) {
      //Get the data currently stored in this tab
      var data = getAllTabData(this.id);

      // Remove any querystring from the url
      if (url.indexOf("?") > 0) {
          url = uri.substring(0, uri.indexOf("?"));
      }

      if(tab.data('json-signal')) {
        //sends a json request to the server and sends the data currently stored in the tab and once it is done 
        $.getJSON(url, data).done(function(json) {

          //calls a function using the json data just fetched, the data already stored in the tabs, and the url used to fetch the data
          tab.data('json-signal')(json, url, data);

          //Signal that the tab data has been fetched and 
          tab.data('done', true);
        })
        //for failure
        .fail(function(jqxhr, textStatus, error) {
          var err = textStatus + ", " + error;
          console.log( "Request Failed: " + err );
        });
      } else if(tab.data('url-signal')) {
          tab.data('url-signal')(url, data);
          tab.data('done', true);
      } else {
          console.error("Couldn't find url handler: ", tab.attr('id'));
      }
    }
    var last_tab = localStorage.setItem("last_tab", tab.attr('id'));
  });
  //Okay now that all the data stuff has been dealt with we just move on now and actually move on the new tab that was just clicked

  var last_tab = localStorage.getItem("last_tab");
  if (last_tab) {
      $('#' + last_tab).addClass("active");
  } else {
      $('.defaultactive').addClass("active");
  }
  $('.defaultactive').removeClass("defaultactive");

  // Activate the existing active tab.
  $(all_tabs+'.active').click();
});

/* Collects all tab values into a dictionary ready for sending to the server */
function getAllTabData(except) {
  
  var data = {};
  $(all_tabs).each(function() {
    if(this.id != except) {
      // Check for existing data
      if($(this).data('values')) {
        
        // Change the way the data is indexed in the dictionary if the data comes with a column tagged data
        var column = this.id.replace('-store', '');
        if($(this).data('column')) {
          column = $(this).data('column')
        }
        data[column] = $(this).data('values');
      }
    }
  });
  return data;
}

/* Returns data for the given tab key name */
function getTabData(key) {
  var store = $('#'+key+'-store');
  return store.data('values');
}
/* Returns data for the given tab key name */
function getTabColumn(key) {
  var store = $('#'+key+'-store');
  return store.data('column');
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
function addTabData(key, value, text, icon, column) {
  var store = $('#'+key+'-store');

  if (!store.data('original-text')) {
    store.data('original-column', store.data('column'));
    store.data('original-text', $('p', store).text());
    store.data('original-icon', $('h2', store).attr('class'));
  }

  if (!store.data('values').includes(value)) {
    store.addClass('selected');
    store.data('column', column);
    store.data('values').push(value);

    $('p', store).text(text);
    $('h2', store).attr('class', 'glyphicon glyphicon-'+icon);
  }
  updateVisuals(key);

  // Updates (id -> name) mapping
  if (text) { store.data('map')[value] = text };
}

/* Removes all elements matching `value` from the `key` store */
function removeTabData(key, value) {
  var store = $('#'+key+'-store');
  store.data('values', store.data('values').filter(function(el) { return el != value }));
  updateVisuals(key);
  delete store.data('map')[value];
}

/* Adds `value` to store if not already present, otherwise removes `value` */
function toggleTabData(key, value, text, icon, column) {
  var store = $('#'+key+'-store');
  if (store.data('values').includes(value)) {
    removeTabData(key, value);
  } else {
    addTabData(key, value, text, icon, column);
  }
}

/* Refreshes tab icon and plots to match new data */
function updateVisuals(key) {
  var store = $('#'+key+'-store');
  var num_vals = store.data('values').length;
  switch (num_vals) {
    case 0:
      store.removeClass('selected');
      $('p', store).text(store.data('original-text'));
      $('h2', store).attr('class', store.data('original-icon'));
      break;
    case 1:
      var remaining = store.data('values')[0];
      $('p', store).text(store.data('map')[remaining]);
      break;
    default:
      $('p', store).text(num_vals+' selected');
  }

  // Clears all existing graph renderings
  $(all_tabs).not(store).data('done', false);
}

/* Deselects tab, removing its value and putting all values back to what they
   where when the page loaded */
function unsetTabData(key) {
  var store = $('#'+key+'-store');
  store.removeClass('selected');
  $(all_tabs).not(store).data('done', false);
  store.removeData('column');
  store.data('values', []);
  store.data('map', {});

  if(store.data('original-text')) {
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