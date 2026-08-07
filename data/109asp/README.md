# 109Asp data policy

The 109Asp inputs are represented by a provenance manifest only. The complete
two-electron integral file is 1.39 GB and the dense spin-orbital tensor would be
38.9 GB, so neither belongs in this repository or in a laptop-scale workflow.

The older 100,000-line archive is also omitted. It is only 0.26% complete and is
from a different run than the recovered complete opt-KS file; combining the two
would be scientifically invalid. Run `scripts/09_109asp_wall.py` to reproduce the
size and completeness audit without opening either source file.
