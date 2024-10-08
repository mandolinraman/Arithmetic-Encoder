# Arithmetic Encoder

This is a Python implementation of the non-adaptive arithmetic encoding/decoding algorithm based mostly on this paper:

- Witten, I.H., Neal, R.M. and Cleary, J.G., 1987. Arithmetic coding for data compression. *Communications of the ACM*, *30*(6), pp.520-540. (https://dl.acm.org/doi/abs/10.1145/214762.214771)

The algorithm can be for both data compression and source generation. This implementation differs from the paper in the following ways:

- We use an $L$-bit register size to hold the code-value. All arithmetic operations use only $(L+1)$-bit registers unlike the paper that uses $2L$ bit registers to hold intermediate results of multiplications.
- Uses a slightly faster way to find which sub-interval contains the code-value in the decoder.
- A way to handle symbol probabilities that are smaller than $2^{-L}$ by using an optional "padding" around sub-intervals.
- The "alt" versions are alternative implementations where the intervals never shrink to less than half of the maximum interval size.
