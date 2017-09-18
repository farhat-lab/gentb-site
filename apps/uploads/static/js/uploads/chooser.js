/*
 * Enable multi-uploader for each uploader field.
 */
$(window).load(function() {
  $('input[type=upload-chooser]').each(function() {
    var field = this;

    // Local uploads
    var resumable = new Resumable({
      target: $(field).data('resumable_url'),
      query: {
        csrfmiddlewaretoken: $('input[name=csrfmiddlewaretoken]').val(),
      }
    });

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
        container.append($("<label>Unknown Files (not uploaded)</label>"));
        container.data('skip', true);
    }

    var warning = $('<span class="text-danger"><br><label for="id_description">Matched Uploads:</label> Some of the files you are uploading have not been matched, these will be ignored.</span>');
    $(this).parent().append(warning);
    warning.hide();

    // Set the form save to collect all the upload information for the field
    var form = $($(field).closest("form"));
    form.on('submit', function(event) {
        // Save the data to the json field
        var files = new Array();
        $.each(containers, function() {
          // Remove any files from skipped buckets
          if($(this).data('skip')) {
            $('li .remove', $(this)).click();
            return;
          }
          // Remove any un-paired files before uploading
          $('li:not(.ok) .remove', $(this)).click();
          // Add each of the uploaded files
          $('li.ok', $(this)).each(function() {
            var file = $(this).data('upload');
            file.css = $(this).attr('class');
            files.push(file);
          });
        });
        field.value = JSON.stringify(files);
        if(resumable.progress() != 1) {
            if(resumable.progress() == 0 && !resumable.isUploading()) {
               $('li.to_upload').addClass('uploading')
                                .removeClass('to_upload');
               resumable.upload();
            }
            event.preventDefault();
            return false;
        }
    });

    add_file = function(file, source) {
        var store = container;
        file.bucket = null;
        file.id = file.name;
        file.id = file.id.replace(/([^A-Za-z0-9[\]{}_.:-])\s?/g, '_');

        // Protect these from re-loading illness
        if(!file.css) { file.css = 'ok'; }
        if(!file.source) { file.source = source; }
        if(!file.name) { return; }

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
                        linked.addClass('ok');
                        linked.removeClass('unlinked');
                    }
                }
            }
            // XXX We could enable drag and drop here if there's multiple buckets.
            // and the file isn't automatically identified.
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
        var close = $('<span class="remove glyphicon glyphicon-remove" title="Remove upload"/>');
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
            var upload_data = li.data('resumable_file');
            if(upload_data && upload_data.source == 'resumable') {
              upload_data.cancel();
            }
            $('#' + file.id, store.data('link')).addClass('unlinked');
            li.remove();
            if($('li', store).length == 0) { store.hide(); }
            update_linked();
        });
        return li;
    };

    options = {
      linkType: "direct",
      multiselect: true,

      success: function(files) {
        var event = new CustomEvent('dropboxChooserSuccess', {'files': files});
        field.dispatchEvent(event);
        $.each(files, function(index, file) {
          add_file(file, 'dropbox');
        });
      },
      cancel: function() {
        field.dispatchEvent(new CustomEvent('dropboxChooserCancel'));
      }
    };
        
    // add extensions only if specified
    if($(field).data("extensions")) {
      options['extensions'] =  $(field).data("extensions").split(" ");
    }

    // Add local uploads button
    var local_button = $('<label for="resumable_upload" class="resumable-chooser dropbox-dropin-btn"><span class="glyphicon glyphicon-upload"></span>Choose from your computer</label>');
    field.parentNode.insertBefore(local_button[0], field);

    // Does the browser support this?
    if(!resumable.support) {
        local_button.addClass('disabled');
        local_button.attr('title', 'Your browser doesn\'t support direct uploads.');
    } else {
      resumable.assignBrowse(local_button[0]);
      var r = resumable;

      r.on('fileSuccess', function(file){
        // Upload is complete for this upload
        if(file.li) {
          file.li.removeClass('uploading');
          file.li.addClass('complete');
        }
      });
      r.on('fileProgress', function(file){
        // File is this much complete
        $('progress', file.li).attr('value', file.progress());
      });
      r.on('fileAdded', function(file, event){
          if(!file || !file.fileName) { return; }
          var data = {
            'id': file.uniqueIdentifier,
            'name': file.fileName,
            'bytes': file.size,
            'uniqueIdentifier': file.uniqueIdentifier,
            'icon': 'https://www.dropbox.com/static/images/icons64/page_white_compressed.png',
          }
          var li = add_file(data, 'resumable');
          li.data('resumable_file', file);
          li.addClass('to_upload');
          li.append($('<progress value="0" max="1"></progress>'));
          file.li = li;
      });
      r.on('fileError', function(file, message){
        // XXX Make UI as error here.
      });
      r.on('complete', function(){
        form.submit();
      });
    }

    // Add dropbox button
    var dropbox_button = Dropbox.createChooseButton(options);
    field.parentNode.insertBefore(dropbox_button, field);

    // Add Manual URL Button
    var url_button = $('<label for="url_upload" class="url-chooser dropbox-dropin-btn"><span class="glyphicon glyphicon-link"></span></label>');
    field.parentNode.insertBefore(url_button[0], field);

    url_button.click(function() {
      bootbox.prompt({
        title: "<span class='glyphicon glyphicon-link'></span> Please enter a URL where your files can be found.<br/><small>This URL can be one of a number of protocols and you should not use it unless you have been instructed to do so.</small>",
        callback: function(result){
          if(result) {
            console.log('EXT');
            console.log($(field).data("extensions"));
            $.ajax({
              type: 'post',
              url: $(field).data('manual_url'),
              data: {
                'url': result,
                'extensions': $(field).data("extensions"),
                'csrfmiddlewaretoken': $('input[name="csrfmiddlewaretoken"]').val(),
              },
              success: function (json) {
                  for(var i = 0; i < json.files.length; i++) {
                    add_file(json.files[i], 'manual');
                  }
                  if(json.files.length == 0) {
                    bootbox.alert({message: 'No matching files found'});
                  }
              },
              error: function (result) {
                  bootbox.alert({message: result.responseText});
              },
            }); 
          }
        }
      });
    });

    // Add any existing files
    if(field.value) {
      var files = JSON.parse(field.value);
      for(var i = 0; i < files.length; i++) {
          add_file(files[i]);
      }
    }
  });
});
