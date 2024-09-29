using StatsBase

include("arithmetic_coding_alt.jl")
# include("arithmetic_coding.jl")

function get_frequencies(text::String)
    freqs = Dict{Char, Int}()
    for c in text
        freqs[c] = get(freqs, c, 0) + 1
    end

    return freqs
end

function test(message, frequencies = nothing)
    if frequencies == nothing
        frequencies = get_frequencies(message)
    end
    endec = make_endec(frequencies, "dc")
    encoded = encode(endec, message)
    decoded = decode(endec, encoded, length(message))

    # encode and decode message
    encoded_string = join(encoded)
    entropy = length(message) * log2(sum(values(frequencies))) - sum(log2(frequencies[c]) for c in message)
    entropy = round(entropy, digits = 2)

    # print info
    println("Entropy of msg = $(entropy) bits")
    println("Encoded length = $(length(encoded)) bits")

    maxlen = 128
    if length(encoded) < maxlen
        println("Message string = \"$message\"")
    else
        println("Message string = \"$(message[1:maxlen])...\"")
    end

    if length(encoded_string) < maxlen
        println("Encoded string = $encoded_string")
    else
        println("Encoded string = $(encoded_string[1:maxlen])...")
    end

    if length(decoded) < maxlen
        println("Decoded string = \"$decoded\"")
    else
        println("Decoded string = \"$(decoded[1:maxlen])...\"")
    end

    @assert message == decoded
end

# test 1
println("\nTest 1\n======")
message = "ABBY CADABBY"
frequencies = Dict('A' => 126, 'B' => 167, 'C' => 116, 'D' => 88, 'Y' => 89, ' ' => 100)
test(message, frequencies)

# test 2
println("\nTest 2\n======")
message = read("test_arithmetic_coding.jl", String)
test(message)

# test 3 - compress random strings
println("\nTest 3\n======")
dc = make_endec(frequencies, "dc", register_size = 16, pad = 1)

letters = collect(keys(frequencies))
weights = collect(values(frequencies))

for trial in 1:10000
    x = join(sample(letters, ProbabilityWeights(weights), 10))
    y = encode(dc, x)
    xh = decode(dc, y, length(x))

    println("\nTrial: $(trial)")
    println(x)
    println(join(y))
    println(xh)

    @assert x == xh
end

# test 4 - generate random strings
println("\nTest 4\n======")
sg = make_endec(frequencies, "sg", register_size = 16, pad = 1)
for trial in 1:10000
    x = rand(0:1, 30)
    y = encode(sg, x)
    xh = decode(sg, y, length(x))

    println("\nTrial: $(trial)")
    println(join(x))
    println(y)
    println(join(xh))

    @assert x == xh
end
