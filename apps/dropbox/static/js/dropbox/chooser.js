/*
 * Enable DropBox upload for each uploader field.
 */
$(window).load(function() {
  $('input[type=dropbox-chooser]').each(function() {
    var field = this;
    Dropbox.appKey = $(field).data("app-key");

    options = {
      linkType: "direct",

      success: function(files) {
        field.value = files[0].link;
        var event = new CustomEvent('dropboxChooserSuccess', { 'file': files[0] });
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
