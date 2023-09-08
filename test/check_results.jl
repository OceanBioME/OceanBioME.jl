file_exists = isfile("results.test")

while !(file_exists)
    sleep(10)
    global file_exists = isfile("results.test")
end

result = open("results.test") do file
    parse(Int, read(file, String))
end

output = open("output.test") do file
    read(file, String)
end

println(output)

result == 0 ? error("Tests failed") : nothing
