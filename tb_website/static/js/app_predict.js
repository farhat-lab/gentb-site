

var hostname="localhost:8080";


// base URL for the R apps: 
var rappURL = "http://0.0.0.0:8000/custom/";


var estimated=false;
var estimateLadda = Ladda.create(document.getElementById("btnEstimate"));

var zparams = { newdata:"" }

function estimate(btn) {
    
    var e = document.getElementById("newseq");
    zparams.newdata = e.options[e.selectedIndex].value;
    //package the zparams object as JSON
    var jsonout = JSON.stringify(zparams);
    var base = rappURL+"rfPredapp?myJSON=";
    var estimated;

    //var test = "{\"x\":[1,2,4,7],\"y\":[3,5,7,9]}";
    //urlcall = base.concat(test);
    urlcall = base.concat(jsonout);
    console.log("urlcall out: ", urlcall);
    
    
    function estimateSuccess(btn,json) {
      estimateLadda.stop();  // stop spinner
      console.log("json in: ", json);
      estimated=true;
      var outmessage = "r: ";
      document.getElementById("output").innerHTML = outmessage.concat(json.r,"<br> s: ",json.s);
    }
    
    function estimateFail(btn) {
      estimateLadda.stop();  // stop spinner
      estimated=false;
    }

    estimateLadda.start();  // start spinner
    makeCorsRequest(urlcall,btn, estimateSuccess, estimateFail);
    
}



// below from http://www.html5rocks.com/en/tutorials/cors/ for cross-origin resource sharing
// Create the XHR object.
function createCORSRequest(method, url, callback) {
    var xhr = new XMLHttpRequest();
    if ("withCredentials" in xhr) {
        // XHR for Chrome/Firefox/Opera/Safari.
        xhr.open(method, url, true);
    } else if (typeof XDomainRequest != "undefined") {
        // XDomainRequest for IE.
        xhr = new XDomainRequest();
        xhr.open(method, url);
    } else {
        // CORS not supported.
        xhr = null;
    }
    return xhr;
    
}


// Make the actual CORS request.
function makeCorsRequest(url,btn,callback, warningcallback) {
    var xhr = createCORSRequest('POST', url);
    if (!xhr) {
        alert('CORS not supported');
        return;
    }
    // Response handlers for asynchronous load
    // onload or onreadystatechange?
    
    xhr.onload = function() {
        
      var text = xhr.responseText;
      console.log("text ", text);
      var json = JSON.parse(text);   // should wrap in try / catch
      var names = Object.keys(json);

      if (names[0] == "warning"){
        warningcallback(btn);
        alert("Warning: " + json.warning);
      }else{
        callback(btn, json);
      }
    }; 
    xhr.onerror = function() {
        // note: xhr.readystate should be 4, and status should be 200.  a status of 0 occurs when the url becomes too large
        if(xhr.status==0) {
            alert('xmlhttprequest status is 0. local server limitation?  url too long?');
        }
        else if(xhr.readyState!=4) {
            alert('xmlhttprequest readystate is not 4.');
        }
        else {
            alert('Woops, there was an error making the request.');
        }
        console.log(xhr);
        estimateLadda.stop();
        selectLadda.stop();
    };
    
    xhr.send();
}



