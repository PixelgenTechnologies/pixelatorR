#' Encode/decode DNA sequences
#'
#' Encode/decode DNA sequences using 2 bits per base for efficient storage, comparison and processing.
#'
#' `pack_2bits`: \cr Converts DNA strings (containing only A, C, G, T) into a 64-bit integers using 2 bits per base.
#' This allows compact storage and efficient comparison of short DNA sequences (kmers).
#'
#' `unpack_2bits`: \cr Given a 64-bit integer encoding a DNA sequence using 2 bits per base (as from \code{pack_2bits}),
#' this function reconstructs the original DNA string.
#'
#' @param seq A character vector of DNA sequences (e.g., "ACGT", "GCAT"). Sequences can only contain
#' A, C, G, or T. All sequences must be of the same length.
#'
#' @rdname two_bit_encoding
#'
#' @return
#' - `pack_2bits`: A vector of class `integer64` representing the 2-bit packed sequences.
#' - `unpack_2bits`: A character string of DNA bases (e.g., "ACGT").
#'
#' @examples
#' # Pack sequence
#' pack_2bits("TCGT") # returns integer64 encoding
#'
#' # Can be vectorized for multiple sequences
#' packed_seqs <- pack_2bits(rep("TCGTGTCGATCTATGCTGATGTCGTGAT", 1e4))
#' head(packed_seqs)
#'
#' @export
pack_2bits <- function(seq) {
  assert_vector(seq, "character", n = 1)
  seq_lens <- sapply(seq, nchar) %>% unique()
  if (!length(seq_lens) == 1) {
    cli::cli_abort(
      c(
        "x" = "Input sequences must be of the same length.",
        "i" = "Lengths: {seq_lens}"
      )
    )
  }

  lookup <- c(A = 0L, C = 1L, G = 2L, T = 3L)

  # Convert sequence string into individual bases
  seq_chars <- strsplit(seq, split = "") %>%
    lapply(rev) %>%
    do.call(rbind, .)

  # Initialize packed as 64-bit integer
  packed <- bit64::as.integer64(0)

  for (i in seq_len(ncol(seq_chars))) {
    val <- lookup[seq_chars[, i]] %>% unname()
    val_is_na <- is.na(val)
    if (any(val_is_na)) {
      cli::cli_abort(
        c(
          "x" = "Invalid base in sequence{?s}: {.str {seq[val_is_na]}}",
          "i" = "Valid bases are A, C, G, T."
        )
      )
    }

    # Shift packed 2 bits to the left (multiply by 4)
    packed <- packed * 4L
    # Add the value corresponding to the current base
    packed <- packed + val
  }

  return(packed)
}


#' @param packed A 64-bit integer (class `integer64`) or character representing the packed DNA sequence.
#' @param k Integer. The number of bases in the original DNA sequence.
#'
#' @rdname two_bit_encoding
#'
#' @examples
#' # Unpack sequence
#' unpack_2bits(packed, 4) # returns "TCGT"
#'
#' # Can be vectorized for multiple sequences
#' unpacked_seqs <- unpack_2bits(packed_seqs, k = 28)
#' head(unpacked_seqs)
#'
#' @export
unpack_2bits <- function(packed, k) {
  assert_class(packed, c("integer64", "character"))

  if (is.character(packed)) {
    packed <- bit64::as.integer64(packed)
  }

  assert_vector(packed, type = "numeric", n = 1)
  assert_single_value(k, type = "numeric")

  lookup <- c("A", "C", "G", "T")

  seq_chars <- matrix(NA_character_, nrow = length(packed), ncol = k)

  for (i in 1:k) {
    # Get last 2 bits
    val <- as.integer(packed %% 4)
    seq_chars[, i] <- lookup[val + 1]

    # Divide by 4 to right shift by 2 bits
    packed <- packed %/% 4
  }

  dna_seq <- apply(seq_chars, 1, paste, collapse = "")

  return(dna_seq)
}
