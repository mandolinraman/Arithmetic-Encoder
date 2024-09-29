import math, random
from arithmetic_coding_alt import ArithmeticCoder

# from arithmetic_coding import ArithmeticCoder


def get_frequencies(text):
    "Estimate frequency table from given text."
    freqs = dict()
    for c in text:
        freqs[c] = freqs.get(c, 0) + 1

    return freqs


def test(message, frequencies=None):
    if frequencies is None:
        frequencies = get_frequencies(message)

    endec = ArithmeticCoder(frequencies, mode="dc", register_size=16)
    encoded = endec.encode(message)
    decoded = endec.decode(encoded, len(message))

    # encode and decode message
    encoded_string = "".join(str(bit) for bit in encoded)
    entropy = len(message) * math.log2(sum(frequencies.values())) - sum(
        math.log2(frequencies[c]) for c in message
    )

    # print info
    maxlen = 128
    print(f"Entropy of msg = {entropy:.2f} bits")
    print(f"Encoded length = {len(encoded)} bits")

    if len(encoded) < maxlen:
        print(f'Message string = "{message}"')
    else:
        print(f'Message string = "{message[:maxlen]}..."')

    if len(encoded_string) < maxlen:
        print(f"Encoded string = {encoded_string}")
    else:
        print(f"Encoded string = {encoded_string[:maxlen]}...")

    if len(decoded) < maxlen:
        print(f'Decoded string = "{decoded}"')
    else:
        print(f'Decoded string = "{decoded[:maxlen]}..."')

    # check if decoding was successful
    assert message == decoded


# test 1
print("\nTest 1\n======")
message = "ABBY CADABBY"
frequencies = {"A": 126, "B": 167, "C": 116, "D": 88, "Y": 89, " ": 100}
test(message, frequencies)

# test2
print("\nTest 2\n======")
with open("test_arithmetic_coding.py") as f:
    message = f.read()
test(message)

# test 3 - compress random strings
print("\nTest 3\n======")
dc = ArithmeticCoder(frequencies, mode="dc", register_size=16, pad=1)

letters = list(frequencies.keys())
weights = list(frequencies.values())
for trial in range(10000):
    x = "".join(random.choices(letters, weights=weights, k=10))
    y = dc.encode(x)
    xh = dc.decode(y, len(x))

    print(f"\nTrial: {trial}")
    print(x)
    print("".join(str(i) for i in y))
    print(xh)

    assert x == xh

# test 4 - generate random strings
print("\nTest 4\n======")
sg = ArithmeticCoder(frequencies, mode="sg", register_size=16, pad=1)

for trial in range(10000):
    x = [random.randint(0, 1) for _ in range(30)]
    y = sg.encode(x)
    xh = sg.decode(y, len(x))

    print(f"\nTrial: {trial}")
    print("".join(str(i) for i in x))
    print(y)
    print("".join(str(i) for i in xh))

    assert x == xh
