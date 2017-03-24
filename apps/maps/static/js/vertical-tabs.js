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
var all_tabs = 'div.vertical-tab-menu > div.list-group > a';

$(document).ready(function() {
  $(all_tabs).click(function(e) {
    var tab = $(this);
    e.preventDefault();
    tab.siblings('a.active').removeClass("active");
    tab.addClass("active");
    var index = tab.index();
    $("div.vertical-tab>div.vertical-tab-content").removeClass("active");
    $("div.vertical-tab>div.vertical-tab-content").eq(index).addClass("active");

    if(tab.data('json-url') && !tab.data('done')) {
      $.getJSON(tab.data('json-url'), getAllTabData(this.id))
        .done(function(json) {
          tab.data('json-signal')(json);
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

function getAllTabData(except) {
  var data = {};
  $(all_tabs).each(function() {
    if(this.id != except) {
      if($(this).data('value')) {
        data[this.id.replace('-store', '')] = $(this).data('value');
      }
    }
  });
  return data;
}

function setTabData(key, value, text, icon) {
  var store = $('#'+key+'-store');

  if(!store.data('original-text')) {
    store.data('original', store.data('value'));
    store.data('original-text', $('p', store).text());
    store.data('original-icon', $('h2', store).attr('class'));
  }
  store.addClass('selected');
  store.data('value', value)
  $('p', store).text(text);
  $('h2', store).attr('class', 'glyphicon glyphicon-'+icon);

  // Clear all existing graph renderings
  $(all_tabs).not(store).data('done', false);
}

function unsetTabData(key) {
  var store = $('#'+key+'-store');
  store.removeClass('selected');
  $(all_tabs).not(store).data('done', false);

  if(store.data('original-text')) {
    store.data('value', store.data('original'));
    $('p', store).text(store.data('original-text'));
    $('h2', store).attr('class', store.data('original-icon'));
  }
}

function chartData(svg, chart, data) {
  return d3.select(svg)
    .datum(data)
    .call(chart);
}

