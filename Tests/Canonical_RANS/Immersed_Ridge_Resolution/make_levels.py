"""make_levels.py: the layer file of the grid-resolution ridge (README.md).

    python3 make_levels.py DZ TOP > zlevels.txt

Uniform layers of DZ m from 0 to 224 m (the crest at 200 m plus 24 m, or plus 22 m for the immersed mesh, whose
ground is 2 m higher), then a ratio of 1.1, ending exactly at TOP: 1002 for inputs_ridge_immersed, 1000 for
inputs_ridge_body_fitted. The layers above 224 m are then stretched (or squeezed) by one factor so the last one is TOP, so the ratio at
224 m is about 1.2 (DZ = 2) to 1.3 (DZ = 32). DZ must divide 224. The layer count goes to stderr.
"""
import sys

import numpy as np


def levels(dz, top, zu=224.0, r=1.1):
    z = list(np.arange(0.0, zu + 0.5 * dz, dz))
    d = dz
    while z[-1] < top:
        d *= r
        z.append(z[-1] + d)
    z = np.array(z[:-1])
    if top - z[-1] < 0.5 * d:
        z = z[:-1]
    z[z > zu] = zu + (z[z > zu] - zu) * (top - zu) / (z[-1] - zu)
    z[-1] = top
    return z


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    dz, top = float(sys.argv[1]), float(sys.argv[2])
    if not (0 < dz <= 32 and top > 224 and abs(224 / dz - round(224 / dz)) < 1e-9):
        sys.exit("make_levels.py: DZ must be in (0, 32] m and divide 224 m, and TOP must be above 224 m")
    z = levels(dz, top)
    sys.stderr.write(f"{len(z) - 1} layers\n")
    print("\n".join(f"{v:.10g}" for v in z))
