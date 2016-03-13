var w = 700;
var h = 800;

var mutnames = Array();
var mutseq = Array();
var mutopacity = Array();


d3.json("/static/local_data/json.txt", function(error, json) {
  if (error) return console.warn(error);
  var outcomedata = json;


  d3.csv("/static/local_data/testdata.csv", function(dataset){
    console.log(dataset[1]);

    var keyArray = Array();
    keyArray = Object.keys(dataset[1]);

    for (var i = 1; i < keyArray.length; i++){   // Should work with shift or slice
      mutnames[i-1] = keyArray[i];
      mutseq[i-1] = i;
      mutopacity[i-1] = 0;
    };

    console.log(mutnames);


  var yshift = 20;
  var xshift = 40;
  var xbase = 110;
  var ybase = 300;
  var ymargin = 0;
  var xmargin = 1;

  var svg =  d3.select("#d3-target-div")
      .append("svg")
      .attr("width", w)
      .attr("height", h);

  var outcome = svg.selectAll("outcome")
      .data(outcomedata.outcome0)
      .enter()
      .append("rect")
      .attr("x", xbase-45)
      .attr("y",  function(d,i){
        return i*yshift + ybase + 2;
      })
      .attr("width", function(d){
        return d*40
      })
      .attr("height", yshift - 2)
      .attr("fill", "teal");

  for (var k = 0; k < mutnames.length; k++){
  var rectname = "rect" +k;

  svg.selectAll(rectname)
      .data(dataset)
      .enter()
      .append("rect")
      .attr("x", xbase + k*xshift)
      .attr("y",  function(d,i){
        return i*yshift + ybase;
      })
      .attr("width", xshift -1)
      .attr("height", yshift -ymargin)
      .attr("fill", function(d){
        return "rgb(" + Math.round(d[mutnames[k]]*250) + ",0,0)"; //"rgb(150,150, " + Math.round(d.aa*250)  + ")";
      });
  };

  svg.selectAll("overlay")
      .data(mutseq)
      .enter()
      .append("rect")
      .attr("x", function(d,i){
        return i*xshift + xbase;
      })
      .attr("y", 0)
      .attr("width", xshift -1)
      .attr("height", dataset.length*yshift + ybase + 8)
      .attr("fill", "pink")
      .attr("opacity", function(d,i){
        return 0.55 * mutopacity[i];
      })
     .on("mouseover", function(d,i){
        var hold = i;
        d3.select(this)
        .attr("fill", "pink") /*function(d){
          if (mutopacity[hold] >0){
            return "pink";
          } else {
            return "none";
          }
        })*/
        .attr("opacity", 0.7) 
        //.attr("stroke","blue");
      })
      .on("mouseout", function(d,i){
        var hold = i;
        d3.select(this)
          .attr("fill","pink")
          .attr("opacity", function(d){
          return 0.55 * mutopacity[hold];
        }) 
        //.attr("stroke","none");
      }) 
      .on("click", function(d,i){
        mutopacity[i] = 1 - mutopacity[i];
        d3.select(this)
        .attr("fill","pink")
        .attr("opacity", function(d){
          return 0.55 * mutopacity[i];
        });
        var myindex = i;
        var mynumber = 0;
        for (i = 0; i < mutopacity.length; i++) { 
          if(mutopacity[i]>0){
            mynumber += Math.pow(2,i);
          }; 
        };
        var myname = "outcome" + mynumber;
        console.log(myname);
        textprob
          .data(outcomedata[myname])
          .transition()
          .duration(300)
          .style("opacity", 0)
          .transition().duration(500)
          .style("opacity", 1)
          .text(function(d) {             
            return Math.round(d*100);
          });
        textvalues
          .data(mutopacity)
          .transition()
          .duration(10)
          .style("opacity", 0)
          .transition().duration(5)
          .style("opacity", 1)
          .text(function(d) {             
            return d;
          });
        outcome
          .data(outcomedata[myname])
          .transition().duration(500)
          .style("opacity", 1)
          .attr("x", xbase-45)
          .attr("y",  function(d,i){
            return i*yshift + ybase + 2;
          })
          .attr("width", function(d){
            return d*40
          })
          .attr("height", yshift - 2)
          .attr("fill", "teal");
      });


  svg.selectAll("text1")
      .data(dataset)
      .enter()
      .append("text")
      .attr("x", 0)
      .attr("y",  function(d,i){
        return (i+0.7)*yshift + ybase;
      })
      .text(function(d){
        return d.id;
      });

  var textprob = svg.selectAll("textprob")
      .data(outcomedata.outcome0)
      .enter()
      .append("text")
      .attr("x", 40)
      .attr("y",  function(d,i){
        return (i+0.7)*yshift + ybase;
      })
      .text(function(d){
        return Math.round(d*100);
      })
      .attr("fill","gray");

 //mutnames = ["aa","bb","cc","dd","ee","ff","gg","hh","ii","jj"];
  svg.selectAll("text2")
      .data(mutnames)
      .enter()
      .append("text")
      .attr("x", 0)
     // .attr("x", function(d,i){
     //   return (i+0.3)*xshift + xbase;
     // })
      .attr("y", 0)
      .text(function(d){
        return d;
      })
      .attr("transform", function(d,i){
        return "translate("+((i+0.6)*xshift+ xbase) +","+(ybase - 30)+") rotate(-90)"  
    })
     // .attr("transform", function(d,i){
     //   return "translate("+(i+0.3)*xshift + xbase+",20)rotate(0)";
     // })
      .attr("fill","black");

  var textvalues = svg.selectAll("textvalues")
      .data(mutopacity)
      .enter()
      .append("text")
      .attr("x", function(d,i){
        return (i+0.4)*xshift + xbase;
      })
      .attr("y", ybase - 10)
      .text(function(d){
        return d;
      })
      .attr("fill","red");

  });

});