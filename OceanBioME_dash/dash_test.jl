using Dash
using PlotlyJS

# Initialize the Dash app
app = dash()

# Define the layout of the app
app.layout = html_div([
    html_h1("Interactive Point Plotter"),
    dcc_slider(
        id="x-slider",
        min=0,
        max=10,
        step=0.1,
        value=5,
        marks=Dict(i => string(i) for i in 0:10),
        tooltip=Dict("always_visible" => true)
    ),
    html_label("X Value", style=Dict("margin-top" => "10px")),
    dcc_slider(
        id="y-slider",
        min=0,
        max=10,
        step=0.1,
        value=5,
        marks=Dict(i => string(i) for i in 0:10),
        tooltip=Dict("always_visible" => true)
    ),
    html_label("Y Value", style=Dict("margin-top" => "10px")),
    dcc_graph(id="scatter-plot")
])

# Define the callback to update the graph
callback!(
    app,
    Output("scatter-plot", "figure"),
    Input("x-slider", "value"),
    Input("y-slider", "value")
) do x_value, y_value
    fig = Plot(
        scatter(
            x=[x_value],
            y=[y_value],
            mode="markers",
            marker=attr(size=10, color="blue")
        ),
        Layout(
 #           title="pH: $(round(pH, digits=2)), fCO2: $(round(fCO2, digits=2)) μatm", # Dynamic title
            xaxis_title="pH",
            yaxis_title="fCO2 (μatm)",
            xaxis=attr(range=[7.8, 8.1]),
            yaxis=attr(range=[380, 390])
        )
    )
    return fig
end

# Run the app
run_server(app, "127.0.0.1", 8050)