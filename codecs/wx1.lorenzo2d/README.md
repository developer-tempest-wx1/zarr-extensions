# `wx1.lorenzo2d` Codec Specification

- **Extension point:** Codec (`Array -> Array`)
- **Name:** `wx1.lorenzo2d`
- **Status:** Experimental / Candidate specification
- **Author:** Jordy van Gils
- **Reference Implementation:** [`wx1-lorenzo`](https://github.com/developer-tempest-wx1/wx1-lorenzo) (Python / `zarr-python >= 3.1`)

## Overview

`wx1.lorenzo2d` is a lossless array-to-array filter codec implementing a 2-D second-order mixed spatial predictor (Lorenzo predictor) combined with ZigZag integer mapping over the last two axes `(axis -2, axis -1)` of fixed-width unsigned integer chunks (`uint8` or `uint16`).

It is intended to precede byte-shuffling and general-purpose entropy compressors (e.g. `bytes` → `numcodecs.shuffle` → `zstd`) on smooth, regular 2-D grids such as latitude/longitude atmospheric fields quantised on an affine integer lattice.

## Applicability & Constraints

- **Input/output rank:** Rank $\ge 2$. Leading dimensions `0 .. ndim - 3` are treated as independent 2-D slices.
- **Supported data types:** `uint8` ($b = 8$) and `uint16` ($b = 16$). Other dtypes MUST be rejected during validation.
- **Chunk shape invariance:** The output chunk shape and dtype are identical to the input chunk shape and dtype.

## Configuration

The codec accepts an empty configuration or no `configuration` key:

```json
{
  "name": "wx1.lorenzo2d"
}
```

## Encoding Algorithm

Given an input chunk $X$ of unsigned $b$-bit integers with shape $(..., N_y, N_x)$:

1. **Horizontal difference** modulo $2^b$:
   - $dx_{..., y, 0} = X_{..., y, 0}$
   - $dx_{..., y, x} = (X_{..., y, x} - X_{..., y, x-1}) \bmod 2^b$ for $1 \le x < N_x$
2. **Vertical difference** modulo $2^b$:
   - $d_{..., 0, x} = dx_{..., 0, x}$
   - $d_{..., y, x} = (dx_{..., y, x} - dx_{..., y-1, x}) \bmod 2^b$ for $1 \le y < N_y$
3. **ZigZag mapping to unsigned $b$-bit integer:**
   Reinterpret $d_{..., y, x}$ as two's-complement signed $b$-bit integer $s \in [-2^{b-1}, 2^{b-1} - 1]$:
   $$z_{..., y, x} = \left((s \ll 1) \oplus (s \gg (b - 1))\right) \bmod 2^b$$
   where $\gg$ denotes arithmetic right shift.

## Decoding Algorithm

Given encoded unsigned $b$-bit integer array $z$ of shape $(..., N_y, N_x)$:

1. **Inverse ZigZag mapping:**
   $$s = (z \gg_{\text{logical}} 1) \oplus -(z \mathbin{\&} 1)$$
   Reinterpret $s \bmod 2^b$ as unsigned $b$-bit integer $d$.
2. **Modular 2-D prefix sums** in unsigned $b$-bit arithmetic (wrapping modulo $2^b$):
   - $dx_{..., y, x} = \sum_{j=0}^{x} d_{..., y, j} \bmod 2^b$
   - $X_{..., y, x} = \sum_{i=0}^{y} dx_{..., i, x} \bmod 2^b$

Reconstruction is deterministic and bit-exact on all valid `uint8` and `uint16` inputs.

## Example `zarr.json` Codec Chain

```json
"codecs": [
  {
    "name": "wx1.lorenzo2d"
  },
  {
    "name": "bytes",
    "configuration": {
      "endian": "little"
    }
  },
  {
    "name": "numcodecs.shuffle",
    "configuration": {
      "elementsize": 2
    }
  },
  {
    "name": "zstd",
    "configuration": {
      "level": 3,
      "checksum": false
    }
  }
]
```
