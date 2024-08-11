def divmod_abc(a, b, c):
    """
    Given three L-bit unsigned integers a, b, c compute:
    
        q = (a * b) // c and r = (a * b) % c

    without overflowing using only L-bit registers for all arithmetic.
    The only requirement is that (2 * c - 1) and q do not overflow
    i.e., they fit in L-bit registers.
    """

    k = 2  # Let k be an integer such that k * c - 1 doesn't overflow
    qm, rm, q, r = a // c, a % c, 0, 0
    m = b

    # We want to maintain the invariant:
    #     a * b / c = (q + r / c) + m * (qm + rm / c)
    #     where 0 <= r, rm < c
    while True:
        # print(m, q, r, qm, rm)
        # Let m = k * t + s with 0 <= s < k:
        s = m % k
        t = m // k

        # Then,
        #     a * b = (q + r/C) + (qm + rm / c) * (k * t + s)
        #     = (q + qm * s) + (r + rm * s) / c + (qm * k + rm * k / c) * t
        r += rm * s  # won't overflow since r + rm * s < c * k
        q += qm * s + r // c  # won't overflow because q <= q_final
        r = r % c

        # At this point q and r have been modified so that we now have
        # a * b / c = (q + r / c) + (qm * k + rm * k / c) * t
        if t == 0:
            break  # we will eventually reach t = 0

        rm *= k  # won't overflow since rm * k < c * k
        qm = k * qm + rm // c  # won't overflow because qm * t <= q_final
        rm = rm % c  # won't overflow
        m = t

        # At this point qm, rm and m have been modified so that
        # the invariant is maintained again. Also m is now strictly
        # smaller.

    # # debug
    # assert q == (a * b) // c
    # assert r == (a * b) % c

    return (q, r)


class ArithmeticCoder:
    def __init__(self, frequencies, register_size=16, pad=1):
        # assert pad >= 1
        self.delta = pad  # pad size
        self.ell = register_size  # register size
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

        # encoder specific
        self.num_inversions = 0
        self.code_output = []

        # decoder specific
        self.message_output = []

    def bit_plus_inversions(self, bit):
        "Helper function"
        bitc = 1 - bit
        # print('Output  : ', str(bit) + str(bitc) * self.num_inversions)
        self.code_output.append(bit)
        for _ in range(self.num_inversions):
            self.code_output.append(bitc)

        self.num_inversions = 0

    def encode(self, message):
        """Encode a message."""
        half = 1 << (self.ell - 1)
        quarter = 1 << (self.ell - 2)
        three_quarters = half + quarter

        low = 0
        high = half + (half - 1)
        self.code_output = []

        # process each symbol in a loop
        for j, sym in enumerate(message):
            # print(f'-\nStarting: low = {low}, high = {high}')
            # print(f'Encoding: c[{j}] = {sym}')
            span = high - low - (self.num * self.delta - 1)
            index, prob, cmf = self.index_of[sym]  # extract index, probability and cmf
            high = (
                low
                + divmod_abc(span, cmf + prob, self.sum_freq)[0]
                + (index + 1) * self.delta
                - 1
            )
            low = low + divmod_abc(span, cmf, self.sum_freq)[0] + index * self.delta
            # print(f'       => low = {low}, high = {high}')

            while True:
                if high < half:
                    self.bit_plus_inversions(0)
                elif low >= half:
                    low -= half
                    high -= half
                    self.bit_plus_inversions(1)
                elif low >= quarter and high < three_quarters:
                    low -= quarter
                    high -= quarter
                    self.num_inversions += 1
                else:
                    break

                # double the interval size
                low = 2 * low
                high = 2 * high + 1

            # print(f'Expanded: low = {low}, high = {high}')

        # At this point we have either
        #   low < quarter < half <= high        => contains [quarter, half]
        #   low < half < three_quarters <= high => contains [half, three_quarters]

        # # According to the paper we output two more bits to specify which of the two
        # # quarter width intervals would be contained in our final interval:
        # self.num_inversions += 1
        # if low < quarter:
        #     self.bit_plus_inversions(0)
        # else:
        #     self.bit_plus_inversions(1)

        # But we can output only one more bit "1" to select the point at "half"
        # which is always contained in the final interval:
        self.bit_plus_inversions(1)

        return self.code_output

    def decode(self, codeword, mlen):
        """Decode a message."""
        half = 1 << (self.ell - 1)
        quarter = 1 << (self.ell - 2)
        three_quarters = half + quarter

        low = 0
        high = half + (half - 1)
        self.message_output = []

        value = 0  # currrent register contents
        codestream = iter(codeword)
        for _ in range(self.ell):
            value = 2 * value + next(codestream, 0)

        for _ in range(mlen):
            # decode next symbol
            span = high - low - (self.num * self.delta - 1)
            q, r = divmod_abc(self.sum_freq, value - low + 1, span)  # = (q + r / span)
            # Let's convert this fraction into a representation of the form:
            #     (q - r / span) where 0 <= r < span
            if r != 0:
                q += 1
                r = span - r

            qi, ri = (0, 0) if self.delta == 0 else divmod_abc(
                self.sum_freq, self.delta, span
            )  # = (qi + ri / span)

            # print(low, high)
            # index, new_low = 0, 0
            for i in range(self.num):
                # print(self.symbol[i], q)
                if self.symbol[i][2] < q:
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
            )
            low = low + divmod_abc(span, cmf, self.sum_freq)[0] + index * self.delta

            # to debug:
            # print(sym, value, low, high)

            self.message_output.append(sym)
            while True:
                if high < half:
                    pass  # do nothing
                elif low >= half:
                    value -= half
                    low -= half
                    high -= half
                elif low >= quarter and high < three_quarters:
                    value -= quarter
                    low -= quarter
                    high -= quarter
                else:
                    break

                low = 2 * low
                high = 2 * high + 1
                value = 2 * value + next(codestream, 0)  # default to 0

        return "".join(self.message_output)


# main

frequencies = {"A": 126, "B": 167, "C": 116, "D": 88, "Y": 89}
endec = ArithmeticCoder(frequencies, 12, 0)

message = "ABBYCADABBY"
code = endec.encode(message)
decoded_message = endec.decode(code, len(message))

print(f"message = {message}")
print(f"decoded = {decoded_message}")
