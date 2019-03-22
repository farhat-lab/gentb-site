/*
 * Enable table uploader for each field.
 */
$(window).on('load', function() {
  $('textarea[type=upload-table]').each(function() {
    var field = this;
    var data = field.dataset;

    var parsers = {};
    var columns = null;
    for(var attr in data) {
      if(attr.startsWith('parser_')) {
        parsers[attr.substring(7)] = new RegExp(data[attr]);
      } else if(attr == 'columns') {
        columns = data[attr].split(',');
      }
    }

    var table = $('<table class="upload-table" id="table-'+field.id+'"><thead><tr class="remote-header"></tr><tr class="local-header"></tr></thead><tbody><tr class="click"></tr></tbody></table>').insertAfter($(field));
    for(var i in columns) {
      $('.remote-header', table).append('<th>'+columns[i]+'</th>');
    }

    /* TASKS TODO:

      * A total line showing number of rows, number of errors (parsers) and number of unmatched rows (based on above unique ID matching)

      * Set column function to tie column index to csv column index
       * Automatically by matching the column names
       * Manually via a user interface addition, maybe a drop down that appears when clicking headers.

      * Have the ability match any column with a unique list from another selection.
       * Get unique ID field/column and record values in a list.
       * Add unique ID listing to the file uploader widget, so files can be matched for unique IDs
       * Add unique ID for site based constants or similar to match them

      * Required columns in order for the form submission to work.

      * Step ONE, load 10 rows from the csv and match up the column names
      * Step TWO, parse the whole csv file converting it to json and passing the filters.
      * Step THREE, save json into the textinput

     */

    $(field).change(function() {
     var row = 0;
     Papa.parse($(this).val(), {
      header: true,
      step: function(results, parser) {
        row += results.data.length;

        if(row == 1) {
          var tr = $('<tr></tr>');
          $('thead', table).append(tr);
          for(var i in results.meta.fields) {
            tr.append($('<th>'+results.meta.fields[i]+'</th>'));
          }
        }

        var tr = $('<tr></tr>');
        $('tbody', table).append(tr);
        for(var j in results.data[0]) {
          tr.append($('<td title="'+j+'">'+results.data[0][j]+'</td>'));
        }
        if(row > 10) {
          parser.pause();
        }
      }
     });
    });

    $('tbody tr.click', table).append($('<td colspan="'+columns.length+'"><input id="file-'+field.id+'" type="file"><label for="file-'+field.id+'" class="resumable-chooser dropbox-dropin-btn"><span class="glyphicon glyphicon-file"></span> Click here to upload</label></td>'));

    var fileInput = $('input', table)[0];
    var readFile = function () {
        var reader = new FileReader();
        reader.onload = function () {
            $(field).val(reader.result);
            $(field).change();
            $('tr.click', table).hide();
        };
        // start reading the file. When it is done, calls the onload event defined above.
        reader.readAsBinaryString(fileInput.files[0]);
    };
    fileInput.addEventListener('change', readFile);

  });
});


