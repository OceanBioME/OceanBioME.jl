file_exists = isfile("results.test")

while !(file_exists)
    sleep(60)
    file_exists = isfile("results.test")
end

result = open("results.test") do file
    read(file, String)
end

if result == "true"
    exit()
else
    error("Tests failed")
end