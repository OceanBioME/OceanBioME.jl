# This script is a simple dashboard for the OceanBioME carbon chemistry solver.
# It allows the user to input DIC, Alkalinity, Temperature, and Salinity, and then
# plots the corresponding pH and fCO2 values.

using Dash
using PlotlyJS
using Measurements
using OceanBioME

# Initialize carbon chemistry solver
carbon_chemistry = CarbonChemistry()

# Initialize the Dash app
app = dash()

# Define styles for consistent appearance
const SLIDER_CONTAINER_STYLE = Dict(
    "display" => "flex",
    "align-items" => "center",
    "margin-bottom" => "20px",
    "padding" => "10px",
    "border-radius" => "8px",
    "background-color" => "#f8f9fa"
)

const LABEL_STYLE = Dict(
    "margin-right" => "15px",
    "width" => "180px",
    "display" => "inline-block",
    "font-weight" => "bold",
    "color" => "#2c3e50",
    "font-family" => "'Segoe UI', Tahoma, Geneva, Verdana, sans-serif"
)

const SLIDER_WRAPPER_STYLE = Dict(
    "display" => "inline-block",
    "width" => "55%",
    "margin-top" => "10px"
)

const ERROR_INPUT_STYLE = Dict(
    "display" => "inline-block",
    "width" => "8%",
    "margin-left" => "15px",
    "padding" => "5px",
    "border" => "1px solid #ddd",
    "border-radius" => "4px",
    "font-family" => "'Segoe UI', Tahoma, Geneva, Verdana, sans-serif"
)

const ERROR_LABEL_STYLE = Dict(
    "margin-left" => "15px",
    "margin-right" => "8px",
    "display" => "inline-block",
    "color" => "#2c3e50",
    "font-weight" => "bold",
    "font-family" => "'Segoe UI', Tahoma, Geneva, Verdana, sans-serif"
)

const HEADER_STYLE = Dict(
    "text-align" => "center",
    "color" => "#2c3e50",
    "margin-bottom" => "30px",
    "font-family" => "'Segoe UI', Tahoma, Geneva, Verdana, sans-serif"
)

const GRAPH_CONTAINER_STYLE = Dict(
    "margin-top" => "20px",
    "padding" => "20px",
    "border-radius" => "10px",
    "background-color" => "white",
    "box-shadow" => "0 2px 10px rgba(0,0,0,0.1)",
    "height" => "600px"
)

const CONTROLS_CONTAINER_STYLE = Dict(
    "width" => "65%",
    "padding" => "20px",
    "background-color" => "white",
    "border-radius" => "10px",
    "box-shadow" => "0 2px 10px rgba(0,0,0,0.1)",
    "margin-right" => "20px"
)

const GRAPH_WRAPPER_STYLE = Dict(
    "width" => "35%",
    "padding" => "20px",
    "background-color" => "white",
    "border-radius" => "10px",
    "box-shadow" => "0 2px 10px rgba(0,0,0,0.1)",
    "height" => "357px",
    "display" => "flex",
    "flex-direction" => "column"
)

const MAIN_CONTAINER_STYLE = Dict(
    "display" => "flex",
    "padding" => "20px",
    "background-color" => "#ecf0f1",
    "min-height" => "100vh",
    "align-items" => "flex-start"
)

# Define the layout of the app
app.layout = html_div([
    html_h1("OceanBioME Carbon Chemistry Solver", style=HEADER_STYLE),
    
    html_div([
        # Left side - Controls
        html_div([
            # DIC control section
            html_div([
                html_label("DIC (mmol/m³):", style=LABEL_STYLE),
                html_div(
                    dcc_slider(
                        id="dic-slider",
                        min=2000,
                        max=2400,
                        step=1,
                        value=2150,
                        marks=Dict(i => string(i) for i in 2000:100:2400),
                        tooltip=Dict("always_visible" => true, "placement" => "bottom")
                    ),
                    style=SLIDER_WRAPPER_STYLE
                ),
                html_label("Error:", style=ERROR_LABEL_STYLE),
                dcc_input(
                    id="dic-error",
                    type="number",
                    value=5.0,
                    step=0.1,
                    style=ERROR_INPUT_STYLE
                )
            ], style=SLIDER_CONTAINER_STYLE),
            
            # Alkalinity control section
            html_div([
                html_label("Alkalinity (mmol/m³):", style=LABEL_STYLE),
                html_div(
                    dcc_slider(
                        id="alk-slider",
                        min=2200,
                        max=2500,
                        step=1,
                        value=2300,
                        marks=Dict(i => string(i) for i in 2200:100:2500),
                        tooltip=Dict("always_visible" => true, "placement" => "bottom")
                    ),
                    style=SLIDER_WRAPPER_STYLE
                ),
                html_label("Error:", style=ERROR_LABEL_STYLE),
                dcc_input(
                    id="alk-error",
                    type="number",
                    value=5.0,
                    step=0.1,
                    style=ERROR_INPUT_STYLE
                )
            ], style=SLIDER_CONTAINER_STYLE),
            
            # Temperature control section
            html_div([
                html_label("Temperature (°C):", style=LABEL_STYLE),
                html_div(
                    dcc_slider(
                        id="temp-slider",
                        min=-2,
                        max=30,
                        step=0.01,
                        value=15,
                        marks=Dict(i => string(i) for i in -2:5:30),
                        tooltip=Dict("always_visible" => true, "placement" => "bottom")
                    ),
                    style=SLIDER_WRAPPER_STYLE
                ),
                html_label("Error:", style=ERROR_LABEL_STYLE),
                dcc_input(
                    id="temp-error",
                    type="number",
                    value=1.0,
                    step=0.01,
                    style=ERROR_INPUT_STYLE
                )
            ], style=SLIDER_CONTAINER_STYLE),
            
            # Salinity control section
            html_div([
                html_label("Salinity (PSS):", style=LABEL_STYLE),
                html_div(
                    dcc_slider(
                        id="salinity-slider",
                        min=31,
                        max=38,
                        step=0.01,
                        value=33,
                        marks=Dict(i => string(i) for i in 31:1:38),
                        tooltip=Dict("always_visible" => true, "placement" => "bottom")
                    ),
                    style=SLIDER_WRAPPER_STYLE
                ),
                html_label("Error:", style=ERROR_LABEL_STYLE),
                dcc_input(
                    id="salinity-error",
                    type="number",
                    value=1.0,
                    step=0.01,
                    style=ERROR_INPUT_STYLE
                )
            ], style=SLIDER_CONTAINER_STYLE)
        ], style=CONTROLS_CONTAINER_STYLE),
        
        # Right side - Graph
        html_div([
            dcc_graph(id="carbon-chemistry-plot", style=Dict("height" => "100%", "width" => "100%"))
        ], style=GRAPH_WRAPPER_STYLE)
    ], style=MAIN_CONTAINER_STYLE)
], style=Dict("background-color" => "#ecf0f1", "min-height" => "100vh"))

# Define the callback to update the graph
callback!(
    app,
    Output("carbon-chemistry-plot", "figure"),
    Input("dic-slider", "value"),
    Input("alk-slider", "value"),
    Input("temp-slider", "value"),
    Input("salinity-slider", "value"),
    Input("dic-error", "value"),
    Input("alk-error", "value"),
    Input("temp-error", "value"),
    Input("salinity-error", "value")
) do dic_value, alk_value, temp_value, salinity_value, dic_error, alk_error, temp_error, salinity_error
    # Convert all values to Float64 to avoid type conversion errors
    dic_value = Float64(dic_value)
    alk_value = Float64(alk_value)
    temp_value = Float64(temp_value)
    salinity_value = Float64(salinity_value)
    dic_error = Float64(dic_error)
    alk_error = Float64(alk_error)
    temp_error = Float64(temp_error)
    salinity_error = Float64(salinity_error)
    
    # Create measurements with uncertainties
    dic_measurement = dic_value ± dic_error
    alk_measurement = alk_value ± alk_error
    temp_measurement = temp_value ± temp_error
    salinity_measurement = salinity_value ± salinity_error

    # Calculate carbon chemistry parameters
    fco2_result = carbon_chemistry(
        DIC=dic_measurement,
        Alk=alk_measurement,
        T=temp_measurement,
        S=salinity_measurement
    )

    ph_result = carbon_chemistry(
        DIC=dic_measurement,
        Alk=alk_measurement,
        T=temp_measurement,
        S=salinity_measurement,
        return_pH=true
    )

    # Extract values and uncertainties
    ph_value = Measurements.value(ph_result)
    ph_uncertainty = Measurements.uncertainty(ph_result)
    fco2_value = Measurements.value(fco2_result)
    fco2_uncertainty = Measurements.uncertainty(fco2_result)

    # Create array of DIC values and calculate corresponding pH and fCO2
    dic_array = range(1800, 2600, length=50)
    ph_contour = []
    fco2_contour = []
    
    for dic in dic_array
        # Calculate pH and fCO2 for this DIC value using current slider settings
        contour_fco2 = carbon_chemistry(DIC=dic, Alk=alk_value, T=temp_value, S=salinity_value)
        contour_ph = carbon_chemistry(DIC=dic, Alk=alk_value, T=temp_value, S=salinity_value, return_pH=true)
        
        push!(ph_contour, Measurements.value(contour_ph))
        push!(fco2_contour, Measurements.value(contour_fco2))
    end

    # Create array of Alkalinity values and calculate corresponding pH and fCO2
    alk_array = range(2000, 2700, length=50)
    ph_alk_contour = []
    fco2_alk_contour = []
    
    for alk in alk_array
        # Calculate pH and fCO2 for this Alkalinity value using current slider settings
        alk_contour_fco2 = carbon_chemistry(DIC=dic_value, Alk=alk, T=temp_value, S=salinity_value)
        alk_contour_ph = carbon_chemistry(DIC=dic_value, Alk=alk, T=temp_value, S=salinity_value, return_pH=true)
        
        push!(ph_alk_contour, Measurements.value(alk_contour_ph))
        push!(fco2_alk_contour, Measurements.value(alk_contour_fco2))
    end

    # Create the plot
    fig = Plot([
        # DIC contour line
        scatter(
            x=fco2_contour,
            y=ph_contour,
            mode="lines",
            line=attr(color="#85c1e9", width=2),
            name="DIC contour",
            showlegend=false,
            hoverinfo="text",
            hovertext=["DIC: $(round(dic, digits=0)) mmol/m³" for dic in dic_array]
        ),
        # Alkalinity contour line
        scatter(
            x=fco2_alk_contour,
            y=ph_alk_contour,
            mode="lines",
            line=attr(color="#82e0aa", width=2),
            name="Alkalinity contour",
            showlegend=false,
            hoverinfo="text",
            hovertext=["Alk: $(round(alk, digits=0)) mmol/m³" for alk in alk_array]
        ),
        # Main measurement point
        scatter(
            x=[fco2_value],
            y=[ph_value],
            mode="markers",
            marker=attr(size=12, color="#3498db", line=attr(color="#2980b9", width=2)),
            name="Measurement",
            showlegend=false
        ),
        # Vertical error bar (fCO2 uncertainty)
        scatter(
            x=[fco2_value - fco2_uncertainty, fco2_value + fco2_uncertainty],
            y=[ph_value, ph_value],
            mode="lines",
            line=attr(color="#3498db", width=2),
            showlegend=false,
            hoverinfo="skip"
        ),
        # Vertical error bar caps
        scatter(
            x=[fco2_value - fco2_uncertainty, fco2_value - fco2_uncertainty, fco2_value - fco2_uncertainty, fco2_value - fco2_uncertainty],
            y=[ph_value - 0.01, ph_value + 0.01, ph_value - 0.01, ph_value + 0.01],
            mode="lines",
            line=attr(color="#3498db", width=2),
            showlegend=false,
            hoverinfo="skip"
        ),
        scatter(
            x=[fco2_value + fco2_uncertainty, fco2_value + fco2_uncertainty, fco2_value + fco2_uncertainty, fco2_value + fco2_uncertainty],
            y=[ph_value - 0.01, ph_value + 0.01, ph_value - 0.01, ph_value + 0.01],
            mode="lines",
            line=attr(color="#3498db", width=2),
            showlegend=false,
            hoverinfo="skip"
        ),
        # Horizontal error bar (pH uncertainty)
        scatter(
            x=[fco2_value, fco2_value],
            y=[ph_value - ph_uncertainty, ph_value + ph_uncertainty],
            mode="lines",
            line=attr(color="#3498db", width=2),
            showlegend=false,
            hoverinfo="skip"
        ),
        # Horizontal error bar caps
        scatter(
            x=[fco2_value - 10, fco2_value + 10, fco2_value - 10, fco2_value + 10],
            y=[ph_value + ph_uncertainty, ph_value + ph_uncertainty, ph_value + ph_uncertainty, ph_value + ph_uncertainty],
            mode="lines",
            line=attr(color="#3498db", width=2),
            showlegend=false,
            hoverinfo="skip"
        ),
        scatter(
            x=[fco2_value - 10, fco2_value + 10, fco2_value - 10, fco2_value + 10],
            y=[ph_value - ph_uncertainty, ph_value - ph_uncertainty, ph_value - ph_uncertainty, ph_value - ph_uncertainty],
            mode="lines",
            line=attr(color="#3498db", width=2),
            showlegend=false,
            hoverinfo="skip"
        )
    ],
        Layout(
            title=attr(
                text="pH: $(round(ph_value, digits=3)) ± $(round(ph_uncertainty, digits=3)) | fCO₂: $(round(fco2_value, digits=1)) ± $(round(fco2_uncertainty, digits=1)) μatm",
                font=attr(size=18, color="#2c3e50", family="'Segoe UI', Tahoma, Geneva, Verdana, sans-serif"),
                y=0.95,
                x=0.5,
                xanchor="center",
                yanchor="top"
            ),
            xaxis=attr(
                title="fCO₂ (μatm)",
                titlefont=attr(size=14, color="#2c3e50"),
                range=[250, 900],
                gridcolor="#ecf0f1",
                zeroline=false
            ),
            yaxis=attr(
                title="pH",
                titlefont=attr(size=14, color="#2c3e50"),
                range=[7.7, 8.2],
                gridcolor="#ecf0f1",
                zeroline=false
            ),
            plot_bgcolor="white",
            paper_bgcolor="white",
            font=attr(family="'Segoe UI', Tahoma, Geneva, Verdana, sans-serif"),
            margin=attr(l=60, r=60, t=80, b=60)
        )
    )
    return fig
end

# Run the app
run_server(app, "127.0.0.1", 8050, debug=true)