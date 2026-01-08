"""Test efficiency scoring."""
from scalpel.design.efficiency import score_spacer

# Test the first guide from TP53
spacer = "GCCCAATAAACCACTCTGAC"
pam = "TGG"

print(f"Testing spacer: {spacer}")
print(f"PAM: {pam}")
print()

score, details = score_spacer(spacer, pam)

print(f"Overall Score: {score:.3f}")
print()
print("Rule Details:")
for name, rule_score in details.items():
    print(f"  {name}: {rule_score.score:.2f} - {rule_score.reason}")
