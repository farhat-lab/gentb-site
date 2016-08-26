/*
 * Enable DropBox upload for each uploader field.
 */
$(window).load(function() {
  $('input[type=dropbox-chooser]').each(function() {
    var field = this;
    Dropbox.appKey = $(field).data("app-key");

    var containers = {};
    var container = $("<ul class='dropbox-files'></ul>").insertAfter($(this));
    container.hide();
    var bucket_inputs = $('input[data-parent=' + field.id + ']');

    if(bucket_inputs.length > 0) {
        bucket_inputs.each(function() {
            var cont = $("<ul class='dropbox-files' id='"+this.id+"_container'></ul>").insertAfter($(this));
            cont.hide();
            if($(this).data('label')) {
                cont.append($("<label>"+$(this).data('label')+"</label>"));
            }
        });
        container.append($("<label>Unknown Files</label>"));
    }

    options = {
      linkType: "direct",
      multiselect: true,

      success: function(files) {
        var event = new CustomEvent('dropboxChooserSuccess', {'files': files});
        field.dispatchEvent(event);

        field.value = '[]';
        $('li', container).remove();
        container.hide();
        bucket_inputs.each(function() {
            this.value = '[]';
            var cont = $('#'+this.id+"_container");
            $('li', cont).remove();
            cont.hide();
        });

        $.each(files, function(index, file) {
          var filled = false;
          bucket_inputs.each(function() {
              var bucket = this;
              $.each($(bucket).data('extensions').split(' '), function(x, match) {
                  if(file.name.endsWith(match)) {
                      var cont = $('#'+bucket.id+"_container");
                      cont.append($("<li><img src='" + file.icon + "'> " + file.name + "</li>"));
                      cont.show();
                      var so_far = JSON.parse(bucket.value);
                      so_far.push(file);
                      bucket.value = JSON.stringify(so_far);
                      filled = true;
                  }
              });
              file.name;
          });
          if(!filled) {
              container.append($("<li><img src='" + file.icon + "'> " + file.name + "</li>"));
              container.show();
              var so_far = JSON.parse(field.value);
              so_far.push(file);
              field.value = JSON.stringify(so_far);
          }
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

    var button = Dropbox.createChooseButton(options);
    field.parentNode.insertBefore(button, field);
  });
});
