# I'm sure there are other libraries that can do this
function parameter_table(labels, values)
    # sanitize a value for use inside an HTML table cell
    clean(v) = let s = "$v"
        s = replace(s, "&" => "&amp;")           # escape first
        s = replace(s, "<" => "&lt;", ">" => "&gt;")
        s = replace(s, r"\r\n|\r|\n" => "<br>")  # then reintroduce breaks as real <br> tags
        s
    end

    rows = join(("<tr><td><code>$label</code></td><td>$(clean(values[idx]))</td></tr>"
             for (idx, label) in enumerate(labels)), "\n")

    return """
```@raw html
    <table><tr><th>Name</th><th>Value</th></tr>
    $rows
    </table>
```
    """
end

function show_parameters(model; exclude = (:optionals, :light_attenuation_model, :particles, :sediment_model, # usually bgc models
                                           :architecture, :x, :y, :z, :A, :N, :C, :nitrate_uptake, :ammonia_uptake, 
                                           :primary_production, :frond_exudation, :nitrogen_erosion, :carbon_erosion, :scalefactor, :custom_dynamics,
                                           :pescribed_temperature, :pescribed_salinity, # usually particles
                                           :field, :fields, :tendencies, # light and sediment
                                           :temperature, :salinity, :pCO₂)) # gas exchange
    names = fieldnames(typeof(model))

    labels = []
    values = []

    for name in names 
        if !(name in exclude)
            if name == :sinking_velocities
                sinking = getproperty(model, name)

                for (name, value) in pairs(sinking)
                    push!(labels, "`$name` sinking speed")
                    push!(values, value)
                end
            else
                push!(labels, "`$name`")
                push!(values, getproperty(model, name))
            end
        end
    end

    return parameter_table(labels, values)
end

function create_parameter_file!(model, name, path)
    output = "# $name default parameters\n\n" * show_parameters(model)
    
    open(path,"w+") do io
        println(io, output)
    end
end
