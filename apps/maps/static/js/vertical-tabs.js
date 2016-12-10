$(document).ready(function() {
    $("div.vertical-tab-menu>div.list-group>a").click(function(e) {
        e.preventDefault();
        $(this).siblings('a.active').removeClass("active");
        $(this).addClass("active");
        var index = $(this).index();
        $("div.vertical-tab>div.vertical-tab-content").removeClass("active");
        $("div.vertical-tab>div.vertical-tab-content").eq(index).addClass("active");
    });
});
