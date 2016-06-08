(function($) {
  $.fn.addOption = function(name, value) {
    var opt = $("<option></option>");
    $(this).append(opt);
    opt.attr("value", value);
    opt.text(name);
    return opt;
  }
  $.fn.replaceOptions = function(options) {
    this.empty();
    var self = this;
    $(self).addOption("---", '---');

    $.each(options, function(index, option) {
      var opt = $(self).addOption(option.name, option.value);
      opt.data('children', option.children);
    }); 
  }

$(document).ready(function() {
  $('.lister').each(function() {
    var last_id = null;
    var previous = null;
    var target = $(this);
    var container = $(this).parent();
    $.getJSON($(this).data('data-url'), function(data) {
      target.data('data', data);
      $.each(data.levels, function(index, level) {
        var select = $('<select title="'+level+'"></select>');
        if(last_id != null) {
          $('#' + last_id).data('next', select);
          select.data('previous', $('#' + last_id));
        }
        last_id = 'level-'+index;
        select.attr('id', last_id);
        select.insertBefore(target);
        select.replaceOptions([]);
        select.change(function() {
          var selected = this.selectedOptions;
          if(!selected || selected[0].label == '---') {
            $(this).data('next').hide();
            $(this).data('next').val('---');
          } else {
            $(this).data('next').show();
            var children = $(selected[0]).data('children');
            if(children) {
              $(this).data('next').replaceOptions(children);
            }
          }
          $(this).data('next').change();
        })
      });

      var add_button = $('<a class="btn btn-success">Add</button>');
      add_button.insertBefore(target);
      add_button.click(function() {
        var value = $('#'+last_id).val();
        if(value && value != '---') {
          $(target).val(value + '\n' + $(target).val());
        }
      });
      $('#' + last_id).data('next', add_button);

      $('#level-0').replaceOptions(data.children);
      $('#level-0').show();
      $('#level-0').change();
    });
  });
});

})(jQuery);
