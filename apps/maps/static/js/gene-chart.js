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

    $('#genelocus-store').data('url-signal', function(url, args) {
      var t_table = $('#gene_map table');

      if(t_table.data('loaded')) {
          return;
      }
      t_table.data('loaded', 1);

      var table = t_table.DataTable({
        "processing": true,
        "serverSide": true,
        "ajax": {
            "url": url,
            "data": function ( data ) {
                // Remove some of the unneeded parts of the column lookup to conserve http request space.
                delete data['search']['regex'];
                for(var i = 0; i < data.columns.length; i++) {
                    col = data.columns[i];
                    delete col['name'];
                    delete col['orderable'];
                    delete col['searchable'];
                    delete col['search']['regex'];
                }
                // Sent json, store for future use in selecting
                $.extend(data, getAllTabData());
                return data;
            },
            "dataSrc": function ( json ) {
                // Returned json, filter and etc here.
                if(json.error != undefined) {
                    console.error("Error getting json:" + json.error);
                    return [];
                }
                return json.data;
            },
        },
        "rowCallback": function(row, data) {
          if (getTabData('genelocus').includes(data.pk)) { $(row).addClass('selected'); }
        },
        "language": {
          "processing": "Loading...",
        },
        "columns": [
          {
            "data": "pk",
            "title": "ID",
            "visible": false,
            "description": "Locus Primary Key",
          },
          {
            "data": "name",
            "title": "Name",
            "description": "Name of the Gene Locus",
            // "render": $.fn.dataTable.render.number(',', '.', 3, ''),
          },
          {
            "data": "gene_symbol",
            "title": "Gene Symbol",
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
        // commenting this out instead of deleting for now for easy reversal
        //   {
        //     "data": "gene_type",
        //     "title": "Gene Type",
        //   },
        //   {
        //     "data": "strand",
        //     "title": "strand",
        //   },
        ],
        'order': [[2, 'asc']],
      });

      table.on('error.dt', function(e, settings, techNote, message) {
          console.log( 'An error has been reported by DataTables: ', message );
        })
        .on('click', 'tbody tr', function () {
          var data = table.row( this.rowIndex - 1 ).data();
          var name = data.str;
          var pk = data.pk;
          if ($(this).hasClass('selected') ) {
            $(this).removeClass('selected');
            removeTabData('genelocus', pk);
          } else {
            $(this).addClass('selected');
            addTabData('genelocus', pk, name, ' icon-helix')
          }
      });
    });
});
