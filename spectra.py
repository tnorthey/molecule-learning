import numpy as np

class Spectra:
    """Manipulate spectra data; apply broadening etc."""

    def __init__(self):
        pass

    def lorenzian_broaden(self, x, y, xmin, xmax, n, fwhm):
        """Apply Lorenzian broadening (FWHM = fwhm) to data y(x),
        outputs new data with length n and min, max = xmin, xmax"""
        x_new = np.linspace(xmin, xmax, n, endpoint=True)
        y_new = np.zeros(n)
        g = (0.5 * fwhm) ** 2  # Factor in Lorentzian function
        for j in range(len(y)):  # loop over original data length
            y_val = y[j]
            x_val = x[j]
            for i in range(n):  # loop over new data size
                lorentz = (
                    y_val * g / ((x_new[i] - x_val) ** 2 + g)
                )  # Lorentzian broadening
                y_new[i] += lorentz
        return x_new, y_new

    def read_bagel_dyson(self, max_rows):
        """ reads dyson norms from bagel output file """
        str_find = "Norms^2 of Dyson orbitals approximately indicate the strength of an inization transitions."
        with open("bagel.out", "r") as f:
            for line in f:
                if str_find in line:  # go to line containing str
                    out_array = np.loadtxt(  # numpy loadtxt into an array
                        f,
                        skiprows=4,  # skip 4 rows of not useful info
                        maxrows=max_rows,  # only read max_rows, errors if this is too large due to eof junk
                        dtype={
                            "names": ("from", "-", "to", "energy", "norm"),
                            "formats": ("i4", "a2", "i4", "f4", "f4"),
                        },
                    )
        return out_array
