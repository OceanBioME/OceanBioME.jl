using Dash
using PlotlyJS
using Measurements
using OceanBioME

carbon_chemistry = CarbonChemistry()

# Initialize the Dash app
app = dash()

# Define the layout of the app
app.layout = html_div([
    html_h1("OceanBioME Carbon Chemistry Solver"),
    html_div([
        html_label("DIC (mmol/m³):", style=Dict("margin-right" => "10px", "width" => "150px", "display" => "inline-block")),
        html_div(
            dcc_slider(
                id="x-slider",
                min=2000,
                max=2400,
                step=1,
                value=2150,
                marks=Dict(i => string(i) for i in 2000:100:2400),
                tooltip=Dict("always_visible" => true)
            ),
            style=Dict("display" => "inline-block", "width" => "60%", "margin-top" => "30px") # Added margin-top
        ),
        html_label("error:  ", style=Dict("margin-left" => "10px", "margin-right" => "5px", "display" => "inline-block")),
        dcc_input(
            id="DIC-error",
            type="number",
            value=5.0,
            step=0.1,
            style=Dict("display" => "inline-block", "width" => "10%") # Positioned to the right
        )
    ], style=Dict("display" => "flex", "align-items" => "center", "margin-bottom" => "20px")),
    html_div([
        html_label("Alkalinity (mmol/m³):", style=Dict("margin-right" => "10px", "width" => "150px", "display" => "inline-block")),
        html_div(
            dcc_slider(
                id="y-slider",
                min=2200,
                max=2500,
                step=1,
                value=2300,
                marks=Dict(i => string(i) for i in 2200:100:2500),
                tooltip=Dict("always_visible" => true)
            ),
            style=Dict("display" => "inline-block", "width" => "60%", "margin-top" => "30px") # Added margin-top
        ),
        html_label("error  ", style=Dict("margin-left" => "10px", "margin-right" => "5px", "display" => "inline-block")),
        dcc_input(
            id="Alk-error",
            type="number",
            value=5.0,
            step=0.1,
            style=Dict("display" => "inline-block", "width" => "10%") # Positioned to the right
        )
    ], style=Dict("display" => "flex", "align-items" => "center", "margin-bottom" => "20px")),
    html_div([
        html_label("T (C):", style=Dict("margin-right" => "10px", "width" => "150px", "display" => "inline-block")),
        html_div(
            dcc_slider(
                id="T-slider",
                min=-2,
                max=30,
                step=0.01,
                value=15,
                marks=Dict(i => string(i) for i in -2:5:30),
                tooltip=Dict("always_visible" => true)
            ),
            style=Dict("display" => "inline-block", "width" => "60%", "margin-top" => "30px") # Added margin-top
        ),
        html_label("error:  ", style=Dict("margin-left" => "10px", "margin-right" => "5px", "display" => "inline-block")),
        dcc_input(
            id="T-error",
            type="number",
            value=0.1,
            step=0.01,
            style=Dict("display" => "inline-block", "width" => "10%") # Positioned to the right
        )
    ], style=Dict("display" => "flex", "align-items" => "center", "margin-bottom" => "20px")),
    html_div([
        html_label("S (PSS):", style=Dict("margin-right" => "10px", "width" => "150px", "display" => "inline-block")),
        html_div(
            dcc_slider(
                id="S-slider",
                min=31,
                max=38,
                step=0.01,
                value=33,
                marks=Dict(i => string(i) for i in 31:1:38),
                tooltip=Dict("always_visible" => true)
            ),
            style=Dict("display" => "inline-block", "width" => "60%", "margin-top" => "30px") # Added margin-top
        ),
        html_label("error:  ", style=Dict("margin-left" => "10px", "margin-right" => "5px", "display" => "inline-block")),
        dcc_input(
            id="S-error",
            type="number",
            value=0.1,
            step=0.01,
            style=Dict("display" => "inline-block", "width" => "10%") # Positioned to the right
        )
    ], style=Dict("display" => "flex", "align-items" => "center", "margin-bottom" => "20px")),
    dcc_graph(id="scatter-plot")
])

# Define the callback to update the graph
callback!(
    app,
    Output("scatter-plot", "figure"),
    Input("x-slider", "value"),
    Input("y-slider", "value"),
    Input("T-slider", "value"),
    Input("S-slider", "value"),
    Input("DIC-error", "value"),
    Input("Alk-error", "value"),
    Input("T-error", "value"),
    Input("S-error", "value")
) do DIC_value, Alk_value, T_value, S_value, DIC_error, Alk_error, T_error, S_error
# Define the DIC and Alk measurements with errors
    DIC_measurement = DIC_value ± DIC_error
    Alk_measurement = Alk_value ± Alk_error
    T_measurement = T_value ± T_error
    S_measurement = S_value ± S_error

#    Calculate fCO2 using the carbon_chemistry function
    fCO2 = carbon_chemistry(DIC=DIC_measurement, 
                            Alk=Alk_measurement, 
                            T=T_measurement, 
                            S=S_measurement)

    pH = carbon_chemistry(DIC=DIC_measurement, 
                            Alk=Alk_measurement, 
                            T=T_measurement, 
                            S=S_measurement,
                            return_pH=true)

    # Create the scatter plot
    fig = Plot([
        scatter(
            x=[Measurements.value(pH)],
            y=[Measurements.value(fCO2)],
            mode="markers",
            marker=attr(size=10, color="blue"),
            name="Measurement"
        ),
        scatter(
            x=[Measurements.value(pH) Measurements.value(pH)],
            y=[Measurements.value(fCO2)-Measurements.uncertainty(fCO2) Measurements.value(fCO2)+Measurements.uncertainty(fCO2)],
            mode="lines",
            line=attr(color="blue", width=2),
            showlegend=false
        ),
        scatter(
            x=[Measurements.value(pH)-Measurements.uncertainty(pH) Measurements.value(pH)+Measurements.uncertainty(pH)],
            y=[Measurements.value(fCO2) Measurements.value(fCO2)],
            mode="lines",
            line=attr(color="blue", width=2),
            showlegend=false
        )
    ],
        Layout(
            title=attr(
            text= "pH: $(round(Measurements.value(pH), digits=3)) ± $(round(Measurements.uncertainty(pH), digits=3)), fCO2: $(round(Measurements.value(fCO2), digits=1)) ± $(round(Measurements.uncertainty(fCO2), digits=1)) μatm",
            y=0.9,
            x=0.5,
            xanchor= "center",
            yanchor= "top"),
            xaxis_title="pH",
            yaxis_title="fCO2 (μatm)",
            xaxis=attr(range=[7.8, 8.1]),
            yaxis=attr(range=[300, 800])
        )
    )
    return fig
end

# Run the app

run_server(app, "127.0.0.1", 8050, debug=true)