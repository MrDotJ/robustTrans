\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 90 Dualc1 - 80 Dualc2 + 80 DualGreatc3 - 80 DualLessc3 - 25 Dualc4
Subject To
 R0: - 42.4 Dualc1 - 2.3 Dualc2 - Dualc4 + Dualc5 = 3.5
 R1: - 3 Dualc1 - Dualc2 + 33 DualGreatc3 - 33 DualLessc3 - 4 Dualc4
   + Dualc6 = 4
 R2: 2.2 DualGreatc3 - 2.2 DualLessc3 + Dualc7 = 0
Bounds
 -infinity <= Dualc1 <= 0
 -infinity <= Dualc2 <= 0
 -infinity <= DualGreatc3 <= 0
 -infinity <= DualLessc3 <= 0
 -infinity <= Dualc4 <= 0
 -infinity <= Dualc5 <= 0
 -infinity <= Dualc6 <= 0
 -infinity <= Dualc7 <= 0
End
