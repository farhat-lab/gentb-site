
$(document).ready(function() {
  setInterval(function () {
      $('.page-loader:visible').each(function() {
          var elem = this.parentElement;
          $(this).removeClass('page-loader');
          $(this).addClass('page-loading');
          var url = $(this).data('loader-url');
          if (!url) {
              $(this).remove();
              console.error("No url to go to for page loader!");
          }
          $.get(url, function(data) {
              $(elem).html(data);
              $(this).remove(); // Maybe not required
              setup_hover();
          });
      });
  }, 500);
    setup_hover();
});

function setup_hover() {
    $("[data-toggle=tooldesc]").data('container', 'body').tooltip({placement: 'bottom', trigger: 'manual'});
    $("[data-toggle=tooltip], [data-toggle=tooltips] *").data('container', 'body').tooltip();
}
