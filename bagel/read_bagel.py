import numpy as np

def read_bagel_dyson(max_rows):

    str_find = "Norms^2 of Dyson orbitals approximately indicate the strength of an inization transitions."
    with open("bagel.out", "r") as f:
        for line in f:
            if str_find in line:    # go to line containing str
                out_array = np.loadtxt(     # numpy loadtxt into an array
                    f,
                    dtype={
                        "names": ("from", "-", "to", "energy", "norm"),
                        "formats": ("i4", "a2", "i4", "f4", "f4"),
                    },
                    skiprows = 4,
                    maxrows = max_rows
                )
    return out_array
