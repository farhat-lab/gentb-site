
$(document).ready(function() {

  // Disable the click for disabled nav-tabs
  $(".nav-tabs > li.disabled > a").on("click", function(e) {
     e.preventDefault();
     return false;
  });

  // Select the best tab
  if($('.nav-tabs > li.active').length == 0) {
     // This selector could be controlable from a data- attribute
     // if we needed more control later (for now we have one mode).
     var activate;
     $('.nav-tabs > li[class!="disabled"]').each(function() {
       var p = $(this).data('priority');
       if(!activate || $(activate).data('priority') < p) {
         activate = this;
       }
     });
     $(activate).addClass('active');
  }

  // Activate the page of the active nav-tabs
  $('.nav-tabs > li.active > a').each(function() {
    var target = $($(this).attr('href'));
    target.addClass('in active');
  });

  $('#activate-notes').click(function() {
      $('#dataselect').toggleClass('col-lg-12 col-lg-8');
      $('#notes').toggleClass('col-lg-4 col-lg-0 hidden');
      $(this).toggleClass('active');
      $('#heatmap').resize();
  });

  $('#addnote').on('submit', function (event) {
    console.log("Adding note!");
    event.preventDefault();

    console.log($(this).attr('action'));
    $.ajax({
      type: 'post',
      url: $(this).attr('action'),
      data: $(this).serialize(),
      success: function (json) {
        if(json.status == 'OK') {
          $('#note').val('');
          $('#notepad').prepend('<li><hr/><strong class="pull-left primary-font">'+json.title+'</strong></br><span class="ui-state-default">'+json.note+'</span></li>')
        } else {
          $('#note').addClass('has-error');
        }
      },
      error: function (result) {
          $('#note').addClass('has-error');
      },
    });
  });

});
