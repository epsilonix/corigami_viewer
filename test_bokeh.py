from flask import Flask, render_template_string
from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.layouts import column

app = Flask(__name__)

@app.route('/test_bokeh')
def test_bokeh():
    # Create a simple Bokeh plot
    p = figure(title="Bokeh Test Plot", width=800, height=400, tools="pan,wheel_zoom,reset")
    p.line([1, 2, 3, 4, 5], [2, 4, 6, 8, 10], line_width=2, line_color="green")
    
    # Embed plot into components
    script, div = components(p)
    
    # Create a simple HTML template inline
    html = render_template_string("""
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="UTF-8">
      <title>Bokeh Test Page</title>
      <link rel="stylesheet" href="https://cdn.bokeh.org/bokeh/release/bokeh-3.6.3.min.css" type="text/css">
        <script src="https://cdn.bokeh.org/bokeh/release/bokeh-3.6.3.min.js"></script>

    </head>
    <body>
      <h1>Bokeh Test Plot</h1>
      <div id="plot-container">
        {{ div|safe }}
      </div>
      {{ script|safe }}
    </body>
    </html>
    """, script=script, div=div)
    
    return html

if __name__ == '__main__':
    app.run(debug=True)
