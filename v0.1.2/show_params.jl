using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years

function show_params(path)
    f = open(path)
    println(read(f, String))
end

