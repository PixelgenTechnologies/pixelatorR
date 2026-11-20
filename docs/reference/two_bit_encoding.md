# Encode/decode DNA sequences

Encode/decode DNA sequences using 2 bits per base for efficient storage,
comparison and processing.

## Usage

``` r
pack_2bits(seq)

unpack_2bits(packed, k)
```

## Arguments

- seq:

  A character vector of DNA sequences (e.g., "ACGT", "GCAT"). Sequences
  can only contain A, C, G, or T. All sequences must be of the same
  length.

- packed:

  A 64-bit integer (class `integer64`) or character representing the
  packed DNA sequence.

- k:

  Integer. The number of bases in the original DNA sequence.

## Value

- `pack_2bits`: A vector of class `integer64` representing the 2-bit
  packed sequences.

- `unpack_2bits`: A character string of DNA bases (e.g., "ACGT").

## Details

`pack_2bits`:  
Converts DNA strings (containing only A, C, G, T) into a 64-bit integers
using 2 bits per base. This allows compact storage and efficient
comparison of short DNA sequences (kmers).

`unpack_2bits`:  
Given a 64-bit integer encoding a DNA sequence using 2 bits per base (as
from `pack_2bits`), this function reconstructs the original DNA string.

## Examples

``` r
# Pack sequence
packed <- pack_2bits("TCGT") # returns integer64 encoding
packed
#> integer64
#> [1] 231

# Can be vectorized for multiple sequences
packed_seqs <- pack_2bits(rep("TCGTGTCGATCTATGCTGATGTCGTGAT", 1e4))
head(packed_seqs)
#> integer64
#> [1] 57314016814210791 57314016814210791 57314016814210791 57314016814210791
#> [5] 57314016814210791 57314016814210791

# Unpack sequence
unpack_2bits(packed, 4) # returns "TCGT"
#> [1] "TCGT"

# Can be vectorized for multiple sequences
unpacked_seqs <- unpack_2bits(packed_seqs, k = 28)
head(unpacked_seqs)
#> [1] "TCGTGTCGATCTATGCTGATGTCGTGAT" "TCGTGTCGATCTATGCTGATGTCGTGAT"
#> [3] "TCGTGTCGATCTATGCTGATGTCGTGAT" "TCGTGTCGATCTATGCTGATGTCGTGAT"
#> [5] "TCGTGTCGATCTATGCTGATGTCGTGAT" "TCGTGTCGATCTATGCTGATGTCGTGAT"
```
