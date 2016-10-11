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
    if($(self).is('select')) {
      $(self).addOption("---", '---');
    }
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
        if(index == 2) {
            var input = $('<input type="text"></input>');
            input.attr('id', "level-"+index);
            input.attr('list', "list-"+index);
            input.attr('title', level);
            input.data('next', 'button');
            input.attr('style', 'width: calc(100% - 45px);');
            input.insertBefore(target);

            var select = $('<datalist></datalist>');
            select.attr('id', "list-"+index);
            select.insertBefore(target);

            input.select(function() {
              var value = $(this).val();
              if(value && value != '---') {
                $('#level-button').removeClass('disabled');
              } else {
                $('#level-button').addClass('disabled');
              }
            });
        } else {
            var select = $('<select></select>');
            select.data('previous', index - 1);
            select.data('next', index + 1);
            select.attr('title', level);
            select.attr('style', 'width: 50%;');
            select.attr('id', 'level-' + index);
            select.insertBefore(target);
            select.replaceOptions([]);

            select.change(function() {
              var selected = this.selectedOptions;
              var next_id = $(this).data('next');
              var next = $('#level-' + next_id);
              if(!selected || selected[0].label == '---') {
                next.prop("disabled", true);
                if(next.is('select')) {
                  next.val('---');
                } else {
                  next.val('');
                }
              } else {
                next.prop("disabled", false);
                var children = $(selected[0]).data('children');
                var list = $('#list-' + next_id);
                if(children) {
                  if(list.length == 1) {
                      list.replaceOptions(children);
                  } else {
                      next.replaceOptions(children);
                  }
                }
              }
              next.change();
              next.select();
            });
        }
      });

      var add_button = $('<a class="btn btn-success btn-sm">Add</button>');
      add_button.insertBefore(target);
      add_button.attr('id', 'level-button');
      add_button.attr('style', 'width: 45;');
      add_button.click(function() {
        var value = $('#level-2').val();
        if(value && value != '---') {
          $(target).val(value + '\n' + $(target).val());
          $('#level-2').val('').select();
        }
      });
      $('#level-0').replaceOptions(data.children);
      $('#level-0').show();
      $('#level-0').change();
    });
  });
});

})(jQuery);
