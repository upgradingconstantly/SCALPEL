"""Quick test of spacer extraction."""

from scalpel.design import SpacerExtractor

# Test sequence with known PAM sites (NGG at position 20-22)
seq = "ATCGATCGATCGATCGATCGAGG"

print(f"Testing spacer extraction on: {seq}")
print(f"Sequence length: {len(seq)}")

extractor = SpacerExtractor()
spacers = extractor.extract_spacers(seq)

print(f"Found {len(spacers)} spacer candidates")
for s in spacers:
    print(f"  Strand: {s.strand.value}, Spacer: {s.spacer_sequence}, PAM: {s.pam_sequence}")
