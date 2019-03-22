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
    //$.fn.dataTable.ext.errMode = 'none';

    //var chart = initialiseMutationChart(svg);
    $('#genelocus-store').data('url-signal', function(url, args) {
        var t_table = $('table');

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
            //"render": $.fn.dataTable.render.number(',', '.', 3, ''),
          },
          {
            "data": "start",
            "title": "Start in Reference Genome",
          },
          {
            "data": "length",
            "title": "Size in Reference Genome",
          },
          {
            "data": "mcount",
            "title": "Mutations Count",
          },
          {
            "data": "gene_type",
            "title": "Gene Type",
          },
          {
            "data": "strand",
            "title": "strand",
          },
        ],
        'order': [[1, 'asc']],
      });

      table.on('error.dt', function(e, settings, techNote, message) {
          console.log( 'An error has been reported by DataTables: ', message );
        })
        .on('click', 'tr', function () {
          var data = table.row( this.rowIndex - 1 ).data();
          if ( $(this).hasClass('selected') ) {
            $(this).removeClass('selected');
            unsetTabData('genelocus');
          } else {
            $('tr.selected', t_table).removeClass('selected');
            $(this).addClass('selected');
            setTabData('genelocus', data.name, data.name, ' icon-helix')
          }
      });
    });
});
