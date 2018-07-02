
$(document).ready(function() {
  $("[data-toggle=tooldesc]").data('container', 'body').tooltip({placement: 'bottom', trigger: 'manual'});
  $("[data-toggle=tooltip], [data-toggle=tooltips] *").data('container', 'body').tooltip();
});
