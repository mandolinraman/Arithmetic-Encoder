# Arithmetic Encoder

This is a Python implementation of arithmetic encoding/decoding (mostly) based on this paper:

- Witten, I.H., Neal, R.M. and Cleary, J.G., 1987. Arithmetic coding for data compression. *Communications of the ACM*, *30*(6), pp.520-540. (https://dl.acm.org/doi/abs/10.1145/214762.214771)

The implementation uses the same $L$-bit register size as the code-value for all operations (multiplications are done without overflow).
