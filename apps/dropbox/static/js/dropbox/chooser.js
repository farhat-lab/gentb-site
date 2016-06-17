/*
 * Enable DropBox upload for each uploader field.
 */
$(window).load(function() {
  $('input[type=dropbox-chooser]').each(function() {
    var field = this;
    Dropbox.appKey = $(field).data("app-key");

    options = {
      linkType: "direct",
      multiselect: true,

      success: function(files) {
        field.value = JSON.stringify(files);
        var event = new CustomEvent('dropboxChooserSuccess', { 'files': files });
        field.dispatchEvent(event);
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
