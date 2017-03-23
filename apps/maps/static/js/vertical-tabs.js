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

$(document).ready(function() {
    $("div.vertical-tab-menu>div.list-group>a").click(function(e) {
        e.preventDefault();
        $(this).siblings('a.active').removeClass("active");
        $(this).addClass("active");
        var index = $(this).index();
        $("div.vertical-tab>div.vertical-tab-content").removeClass("active");
        $("div.vertical-tab>div.vertical-tab-content").eq(index).addClass("active");
    });
});

function setTabData(key, value, text, icon) {
  var store = $('#'+key+'-store');

  if(!store.data('original-text')) {
    store.data('original', store.data('value'));
    store.data('original-text', $('p', store).text());
    store.data('original-icon', $('h2', store).attr('class'));
  }
  store.data('value', value)
  $('p', store).text(text);
  $('h2', store).attr('class', 'glyphicon glyphicon-'+icon);
}

function unsetTabData(key) {
  var store = $('#'+key+'-store');

  if(store.data('original-text')) {
    store.data('value', store.data('original'));
    $('p', store).text(store.data('original-text'));
    $('h2', store).attr('class', store.data('original-icon'));
  }
}

