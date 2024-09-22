import math, random

ALERT = {False: "", True: " ####"}


def divmod_abc(a, b, c):
    """
    Given three unsigned integers a, b, c such that a, b and 2*c - 1 fit in
    L-bit registers, compute:

        q = (a * b) // c and r = (a * b) % c

    without overflowing using only L-bit registers for all arithmetic.
    """

    # It's slightly more efficient if b < a
    if a < b:
        return divmod_abc(b, a, c)

    # Let k > 1 be an integer such that k * c - 1 doesn't overflow. We
    # choose k = 2 which requires that 2 * c - 1 is an L-bit value
    k = 2
    qm, rm = divmod(a, c)
    q, r = 0, 0
    m = b

    # We want to maintain the invariant:
    #     a * b / c = (q + r / c) + (qm + rm / c) * m
    #     where 0 <= r, rm < c
    while m > 0:
        # Let m = k * t + s with 0 <= s < k:
        (t, s) = divmod(m, k)

        # Then,
        #     a * b = (q + r/C) + (qm + rm / c) * (k * t + s)
        #     = (q + qm * s) + (r + rm * s) / c + (qm * k + rm * k / c) * t
        if s > 0:
            r += rm * s  # won't overflow since r + rm * s < c * k
            q += qm * s + r // c  # won't overflow because q <= q_final
            r = r % c

        # At this point q and r have been modified so that we now have
        # a * b / c = (q + r / c) + (qm * k + rm * k / c) * t
        if t > 0:
            # simplify the second term
            rm *= k  # won't overflow since rm * k < c * k
            qm = k * qm + rm // c  # won't overflow because qm * t <= q_final
            rm = rm % c  # won't overflow

        m = t

        # At this point qm, rm and m have been modified so that
        # we maintain the invariant:
        #   a * b / c = (q + r / c) + (qm + rm / c) * m
        #
        # and m is now strictly smaller.

    # # debug
    # assert q == (a * b) // c
    # assert r == (a * b) % c

    return (q, r)


def get_frequencies(text):
    "Estimate frequency table from given text."
    freqs = dict()
    for c in text:
        freqs[c] = freqs.get(c, 0) + 1

    return freqs


class ArithmeticCoder:
    def __init__(self, frequencies, mode="dc", register_size=16, pad=0):
        self.mode = mode.lower()  # "dc" for data compression, "sg" for source generation
        assert self.mode in ["dc", "sg"]

        self.delta = pad  # pad size
        self.ell = register_size  # register size "L"
        self.num = len(frequencies)  # number of symbols

        # compute cmf table
        self.index_of = dict()
        self.symbol = [None] * self.num
        cmf = 0
        for index, (sym, prob) in enumerate(frequencies.items()):
            self.index_of[sym] = (index, prob, cmf)
            self.symbol[index] = (sym, prob, cmf)
            cmf += prob

        # we want 2 * cmf - 1 to fit in a register:
        assert cmf <= (1 << (self.ell - 1))
        self.sum_freq = cmf

    def compress(self, message, length=float("inf")):
        """Convert text to a bit sequence."""

        # local helper function
        def increment():
            nonlocal code_output

            k = len(code_output) - 1
            while code_output[k] == 1:
                code_output[k] = 0
                k -= 1
            code_output[k] = 1

        code_output = []

        full = 1 << self.ell
        half = 1 << (self.ell - 1)

        # low and high are always L-bit values
        low = 0
        high = 2 * half - 1

        n = 0  # keeps track of number of doubling operations
        m = 0  # keeps track of number of symbols generated
        while m < len(message) and n < length:
            sym = message[m]
            # print(f"\nStarting: [{low/full}, {(high+1)/full})" + alert[high < low])
            # print(f"Received: symbol[{m}] = {sym}")
            oldlow = low
            span = (
                (high - low) % full + 1 - self.num * self.delta
            )  # could become (L+1) bit value
            index, prob, cmf = self.index_of[sym]  # extract index, probability and cmf
            high = (
                low
                + divmod_abc(span, cmf + prob, self.sum_freq)[0]
                + (index + 1) * self.delta
                - 1
            ) % full
            low = (
                low + divmod_abc(span, cmf, self.sum_freq)[0] + index * self.delta
            ) % full
            # print(f"       => [{low/full}, {(high+1)/full})")

            if low < oldlow:
                increment()

            m += 1  # increment symbol count

            # check if we can double
            while (high - low) % full < half and n < length:
                if low >= half:
                    code_output.append(1)
                else:
                    code_output.append(0)

                # double the interval size
                low = (2 * low) % full
                high = (2 * high + 1) % full
                n += 1  # count number of doublings

                # print(f"Expanded: [{low/full}, {(high+1)/full})")
            # print("".join(str(x) for x in code_output))

        # print(f"Compressor: n = {n}, mlen = {len(message)}")
        # print(f"low = {low}, high = {high}, half = {half}")

        # needed only for data compression mode:
        if self.mode == "dc":
            # add one more bit to pick a point in the final interval
            code_output.append(1)
            if low >= half:
                increment()

        return code_output

    def expand(self, codeword, length=float("inf")):
        """Convert a bit sequence to text."""

        message_output = []

        full = 1 << self.ell
        half = 1 << (self.ell - 1)

        low = 0
        high = 2 * half - 1
        value = 0  # currrent register contents
        codestream = iter(codeword)
        tail_bit = int(self.mode == "sg")  # = 1 for source generation

        for _ in range(self.ell):
            value = 2 * value + next(codestream, tail_bit)

        n = 0  # keeps track of number of doubling operations
        m = 0  # keeps track of number of symbols generated
        while m < length and n < len(codeword):
            # print(f"\nStarting: [{low/full}, {(high+1)/full})" + alert[high < low] + f" | codeval = {value/full}")

            # decode next symbol
            span = (
                (high - low) % full + 1 - self.num * self.delta
            )  # could become (L+1) bit value

            # Our goal is to find the largest symbol index i such that
            #     value - low >= floor(span * cmf[i] / sum_freq)
            # <=> value - low + 1 - i * delta > span * cmf[i] / sum_freq
            # <=> cmf[i] < (value - low + 1 - i * delta) * sum_freq / span
            # <=> cmf[i] < ceil((value - low + 1 - i * delta) * sum_freq / span)
            #
            # We'll precompute d * sum_freq / span (where d = value - low + 1)
            # and delta * sum_freq / span (when delta != 0). Then we can efficiently
            # compute the RHS above as we increment i

            d = (value - low) % full + 1  # could become (L+1) bit value
            q, r = divmod_abc(d, self.sum_freq, span)
            # Let's convert this fraction into a representation of the form:
            #     (q - r / span) where 0 <= r < span
            # so that q = ceil(d * self.sum_freq / span)
            if r != 0:
                q += 1
                r = span - r

            qi, ri = (
                (0, 0)
                if self.delta == 0
                else divmod_abc(self.sum_freq, self.delta, span)
            )

            for i, (_, _, cmf) in enumerate(self.symbol):
                # print(self.symbol[i], q)
                if cmf < q:
                    index = i  # i could be the next symbol
                else:
                    break

                if self.delta != 0:
                    # update (q, r) -> (q - qi, r + ri) for next iteration
                    r += ri
                    q = q - qi - r // span
                    r = r % span

            sym, prob, cmf = self.symbol[index]  # extract symbol, probability and cmf
            high = (
                low
                + divmod_abc(span, cmf + prob, self.sum_freq)[0]
                + (index + 1) * self.delta
                - 1
            ) % full
            low = (
                low + divmod_abc(span, cmf, self.sum_freq)[0] + index * self.delta
            ) % full

            # print(f"Decoded: symbol[{m}] = {sym}")
            # print(f"       => [{low/full}, {(high+1)/full})")

            message_output.append(sym)
            m += 1  # increment symbol count

            # check if we can double
            while (high - low) % full < half and n < len(codeword):
                # double the interval size
                low = (2 * low) % full
                high = (2 * high + 1) % full
                value = (2 * value + next(codestream, tail_bit)) % full
                n += 1  # count number of doublings
                # print(f"Expanded: [{low/full}, {(high+1)/full}) | codeval = {value/full}")
            # print()

        # print(f"Expander: n = {n}, cwlen = {len(codeword)}")
        # print(f"low = {low}, high = {high}, half = {half}")
        return "".join(message_output)

    def encode(self, seq):
        if self.mode == "dc":
            # seq is text
            return self.compress(seq)
        else:
            # seq is an bit sequence
            return self.expand(seq)

    def decode(self, seq, length):
        if self.mode == "dc":
            # seq is an bit sequence
            return self.expand(seq, length)
        else:
            # seq is text
            return self.compress(seq, length)


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


if __name__ == "__main__":
    # main: example usage

    # test 1
    print("\nTest 1\n======")
    message = "ABBY CADABBY"
    frequencies = {"A": 126, "B": 167, "C": 116, "D": 88, "Y": 89, " ": 100}
    test(message, frequencies)

    # test2
    print("\nTest 2\n======")
    with open("arithmetic_coding_alt.py") as f:
        message = f.read()
    test(message)

    print("\nTest 3\n======")
    endec = ArithmeticCoder(frequencies, mode="sg", register_size=16)

    for trial in range(10000):
        print(f"\nTrial: {trial}")
        x = [random.randint(0, 1) for _ in range(30)]
        y = endec.encode(x)
        xh = endec.decode(y, len(x))

        print("".join(str(i) for i in x))
        print(y)
        print("".join(str(i) for i in xh))
        assert x == xh
