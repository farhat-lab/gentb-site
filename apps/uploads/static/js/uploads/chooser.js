/*
 * Enable multi-uploader for each uploader field.
 */
$(window).load(function() {
  $('input[type=upload-chooser]').each(function() {
    var field = this;

    Dropbox.appKey = $(field).data("app-key");

    var container = $("<ul class='uploads-files'></ul>").insertAfter($(this));
    container.hide();

    var containers = new Array();
    containers.push(container);

    var bucket_inputs = $('input[data-parent=' + field.id + ']');

    if(bucket_inputs.length > 0) {
        bucket_inputs.each(function() {
            var store = $("<ul class='uploads-files' id='"+this.id+"_container'></ul>").insertAfter($(this));
            store.hide();
            containers.push(store);
            if($(this).data('label')) {
              store.append($("<label>"+$(this).data('label')+"</label>"));
              $(this).data('match', new RegExp($(this).data('match')));
            }
        });
        bucket_inputs.each(function() {
          var store = $('#'+this.id+"_container");
          var target = $('[data-bucket="'+$(this).data('link')+'"]');
          store.data('link', $('#' + target[0].id + '_container'));
        });
        container.append($("<label>Unknown Files</label>"));
    }

    var warning = $('<span class="text-danger"><br><label for="id_description">Matched Uploads:</label> Some of the files you are uploading have not been matched, these will be ignored.</span>');
    $(this).parent().append(warning);
    warning.hide();

    // Set the form save to collect all the upload information for the field
    var form = $($(field).closest("form"));
    if(!form.data('uploader-active')) {
      form.on('submit', function() {
        // Save the data to the json field
        var files = new Array();
        $.each(containers, function() {
          $('li.ok', $(this)).each(function() {
            var file = $(this).data('upload');
            file.css = $(this).attr('class');
            files.push(file);
          });
        });
        field.value = JSON.stringify(files);
      });
      form.data('uploader-active', true);
    }

    add_files = function(files) {
      $.each(files, function(index, file) {
          var store = container;
          file.bucket = null;
          file.id = file.name;
          file.css = 'ok';

          bucket_inputs.each(function() {
              var bucket = this;
              var matched = $(bucket).data('match').exec(file.name);
              if(matched && matched[1]) {
                  file.bucket = $(bucket).data('bucket');
                  store = $('#'+bucket.id+"_container");
                  // Cross referencing buckets
                  file.id = matched[1];
                  // 1. match with other buckets if required.
                  var link = store.data('link');
                  if(link) {
                      var linked = $('#' + file.id, link);
                      if(linked.length == 0) {
                          file.css = 'unlinked';
                      } else {
                          linked.attr('class', 'ok');
                      }
                  }
              }
          });

          // Figure out if it's already here.
          var already = false;
          $('li', store).each(function() {
            if($(this).data('upload').name == file.name) {
              already = true;
            }
          });
          if(already) { return false; }

          // Build visual icon
          var li = $("<li><img src='" + file.icon + "'> " + file.name + "</li>");
          var close = $('<span class="glyphicon glyphicon-remove" title="Remove upload"/>');
          li.append(close);
          li.attr('id', file.id);
          li.attr('class', file.css);
          li.data('upload', file);
          store.append(li);
          store.show();
          // Show warning if needed about unlinked files
          var update_linked = function() {
            if($('li.unlinked').length > 0) {
              warning.show();
            } else {
              warning.hide();
            }
          };
          update_linked();
          // Create a close button and connect up the signal
          close.click(function() {
            $('#' + file.id, store.data('link')).attr('class', 'unlinked');
            li.remove();
            if($('li', store).length == 0) { store.hide(); }
            update_linked();
          });
      });
    };

    options = {
      linkType: "direct",
      multiselect: true,

      success: function(files) {
        var event = new CustomEvent('dropboxChooserSuccess', {'files': files});
        field.dispatchEvent(event);
        add_files(files);
      },
      cancel: function() {
        field.dispatchEvent(new CustomEvent('dropboxChooserCancel'));
      }
    };
        
    // add extensions only if specified
    if($(field).data("extensions")) {
      options['extensions'] =  $(field).data("extensions").split(" ");
    }

    var button = Dropbox.createChooseButton(options);
    field.parentNode.insertBefore(button, field);
  });
});
