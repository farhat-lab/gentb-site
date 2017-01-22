
$(document).ready(function() {
    $(".nav-tabs > li.disabled > a").on("click", function(e) {
        e.preventDefault();
        return false;
    });
});
