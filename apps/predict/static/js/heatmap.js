function scatter_plot(data) {
  /*  
   * The idea here is to update an existing plot with new data.
   */
  var svg = d3.select('#scatter svg');

  svg[0][0].__data__ = data;
  var nvChart = svg[0][0].__chart__;

  nvChart.xAxis
      .axisLabel('Genetic Region')
      .tickFormat(function(d) {
         if (typeof this != 'undefined') {
           var el = d3.select(this);
           // Total width - yaxis_width / number of cols;
           var width = (1000 - 150) / data[0].cols.length;
           var parentNode = d3.select(this.parentNode);

           var p = this.replacement;

           if(p === undefined) {
             p = parentNode.append("foreignObject")
                .attr("width", 200)
                .attr("height", 200)
              .append("xhtml:p")
                .attr('style','word-wrap: break-word; text-align:center;');
             this.replacement = p;
           }

           var vp = d3.select(p[0][0].parentNode);

           vp.attr('x', 0-(width / 2))
             .attr("width", width);

           p.html(data[0].cols[d]);
         }
         // Return blank string to svg text since we've replaced it.
         return '';
      });
  svg.data([data]).transition().duration(500).call(nvChart);
  nv.utils.windowResize(nvChart.update);
}

$(document).ready(function() {

  nv.addGraph(function() {
    var data = Array();
    var chart = nv.models.multiBarChart()
      .yDomain([0, 4])
      .reduceXTicks(false);

    var width = 1000;
    var height = 300;

    chart.noData("click on a box on the heamap to show mutations");
    chart.margin({top: 20, right: 60, bottom: 60, left: 60});
    chart.height(height);
    chart.width(width);
    chart.yAxis.scale(100).orient("left")
         .axisLabel('Number of mutations')
         .tickFormat(d3.format("d"))
         .tickValues([0,1,2,3,4,5]);

    chart.showLegend(true);

    var svg = d3.select('#scatter svg')
          .attr('perserveAspectRatio', 'xMinYMid')
          .attr('width', width)
          .attr('height', height)
          .datum(data)
          .attr('viewBox', '0 0 ' + width + ' ' + height)
          .transition().duration(1200)
          .call(chart);

    chart.tooltip.contentGenerator(function (data) {
      ret = "<table>";
      if(data.data) {
        $.each(data.data.tip, function(x, value) {
          ret += "<tr><th align=\"right\">" + value + "</th></tr>";
        });
      }
      return ret + "</table>";
    });
    svg[0][0].__chart__ = chart;
  });

  $('table.heatmap td div.cell').click(function() {
      $('table.heatmap .selected').removeClass('selected');
      $(this).addClass('selected');
      $('#scatter_title').text('Mutation plot for drug='+$(this).data('col')+', strain=' + $(this).data('row'));
      var data = $(this).data('scatter');
      scatter_plot(data);
  });

});
