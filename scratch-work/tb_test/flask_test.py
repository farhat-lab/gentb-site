import sys
from flask import Flask
application = Flask(__name__)

@application.route("/")
def hello():
    return "Hello World!<br /><b>python version:</b> {0}".format(sys.version)

if __name__ == "__main__":
    #127.0.0.1:5000
    application.run(debug=True, port=5000)
