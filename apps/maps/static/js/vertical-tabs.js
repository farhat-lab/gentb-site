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

    //If there is actually data and done is false ('meaning that we didn't already do this whole process)
    if(url && !tab.data('done')) {
      //Get the data currently stored in this tab
      var data = getAllTabData(this.id);

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
    }
  });
  //Okay now that all the data stuff has been dealt with we just move on now and actually move on the new tab that was just clicked

  // Activate the existing active tab.
  $(all_tabs+'.active').click();
});

/* Collects all tab values into a dictionary ready for sending to the server */
function getAllTabData(except) {
  
  //Access current tab 
  var data = {};
  $(all_tabs).each(function() {
  	//The ID passed into except is the way we can access the current tab
    if(this.id != except) {
   	  //see that if some 'value' data is storred if it is
      if($(this).data('value')) {
      	
      	//Change the way the data is indexed in the dictionary if the data comes with a column tagged data
        var column = this.id.replace('-store', '');
        if($(this).data('column')) {
          column = $(this).data('column')
        }

        //Storee the value data into the data dictionary
        data[column] = $(this).data('value');
      }
    }
  });

  //Return that data
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

