file_exists = isfile("results.test")

while !(file_exists)
    sleep(10)
    global file_exists = isfile("results.test")
end

result = open("results.test") do file
    read(file, String)
end

result == 0 ? error("Tests failed") : nothing
