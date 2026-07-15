# I'm sure there are other libraries that can do this
function parameter_table(labels, values)
    # sanitize a value so it stays inside a single markdown table row
    clean(v) = let s = "$v"
        s = replace(s, r"\r?\n\s*" => "<br>")   # newlines -> HTML break (drop indent)
        s = replace(s, "|" => "\\|")            # escape pipes so they don't split cells
        s
    end

    output = "|Name | Value |\n|---|---|\n"

    for (idx, label) in enumerate(labels)
        output *= "|$label|$(clean(values[idx]))|\n"
    end

    return output 
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
