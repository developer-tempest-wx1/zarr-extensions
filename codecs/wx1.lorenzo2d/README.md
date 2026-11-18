# `wx1.lorenzo2d` codec

Defines an `array -> array` lossless spatial predictor codec (2-D Lorenzo predictor followed by ZigZag mapping) for integer-quantised 2-D or higher-rank grids (e.g. regular latitude/longitude meteorological grids).

- **Extension point:** Codec (`array -> array`)
- **Maintainer:** Jordy van Gils (`developer-tempest-wx1`)
- **JSON Schema:** [`schema.json`](./schema.json)
- **Reference Implementation:** [`wx1-lorenzo`](https://github.com/developer-tempest-wx1/wx1-lorenzo) ([PyPI](https://pypi.org/project/wx1-lorenzo/), Python / `zarr >= 3.1`)

## Codec name

The value of the `name` member in the codec object MUST be `wx1.lorenzo2d`.

## Configuration parameters

The codec takes no required configuration parameters. A `configuration` member MAY be omitted or set to an empty object `{}`.

## Applicability and constraints

- **Rank:** Input arrays MUST have rank ≥ 2. The last two axes (`axis -2` = `y`, `axis -1` = `x`) are treated as the 2-D spatial plane; all leading axes are processed as independent 2-D slices.
- **Data types:** Only unsigned 8-bit (`uint8`, `b = 8`) and unsigned 16-bit (`uint16`, `b = 16`) integers are supported. Other data types MUST be rejected during validation.
- **Chunk invariance:** Output chunk shape and data type are identical to input chunk shape and data type.

## Example

For example, the array metadata below chains `wx1.lorenzo2d` before byte serialisation, 2-byte shuffle, and `zstd`:

```json
{
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
}
```

## Format and algorithm

### Encoding

Given an input chunk `X` of unsigned `b`-bit integers (`uint8` or `uint16`) with shape `(..., Ny, Nx)`:

1. **Horizontal difference** modulo `2^b`:
   - `dx[..., y, 0] = X[..., y, 0]`
   - `dx[..., y, x] = (X[..., y, x] - X[..., y, x - 1]) mod 2^b` for `1 <= x < Nx`
2. **Vertical difference** modulo `2^b`:
   - `d[..., 0, x] = dx[..., 0, x]`
   - `d[..., y, x] = (dx[..., y, x] - dx[..., y - 1, x]) mod 2^b` for `1 <= y < Ny`
3. **ZigZag mapping** to unsigned `b`-bit integer:
   Reinterpret `d[..., y, x]` as a two's-complement signed `b`-bit integer `s` (`int8` or `int16`):
   - `z[..., y, x] = ((s << 1) ^ (s >> (b - 1))) mod 2^b` (using arithmetic right shift `>>`)

### Decoding

Given an encoded unsigned `b`-bit integer array `z` of shape `(..., Ny, Nx)`:

1. **Inverse ZigZag mapping**:
   - `s = (z >> 1) ^ -(z & 1)` (using unsigned/logical right shift `>>`)
   - Reinterpret `s mod 2^b` as unsigned `b`-bit integer `d`.
2. **Modular 2-D prefix sums** in unsigned `b`-bit arithmetic (wrapping modulo `2^b`):
   - `dx[..., y, x] = sum_{j=0..x} d[..., y, j] mod 2^b`
   - `X[..., y, x] = sum_{i=0..y} dx[..., i, x] mod 2^b`

Reconstruction is deterministic and bit-exact across all valid `uint8` and `uint16` inputs.

### Compatibility and interoperability

On smooth regular latitude/longitude meteorological fields stored on an affine GRIB integer lattice, `wx1.lorenzo2d` reduces residual entropy ahead of byte-shuffle and `zstd`. The transform assumes a rectangular 2-D spatial grid and does not directly apply to irregular grids such as reduced Gaussian grids without per-ring handling.
