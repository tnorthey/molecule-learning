import numpy as np


class Spectra:
    """Manipulate spectra data; apply broadening etc."""
    def __init__(self):
        pass

    def lorenzian_broaden(self, x, y, xmin, xmax, n, fwhm):
        """Apply Lorenzian broadening (w. FWHM = fwhm) to data y(x),
        outputs new data with length n and min, max = xmin, xmax"""
        x_new = np.linspace(xmin, xmax, n, endpoint=True)
        y_broadened = np.zeros(n)

        g = (0.5 * fwhm) ** 2  # Factor in Lorentzian function
        for j in range(len(y)):  # loop over original data length
            y_val = y[j]
            x_val = x[j]
            for i in range(n):  # loop over new data size
                lorentz = (
                    y_val * g / ((x_new[i] - x_val) ** 2 + g)
                )  # Lorentzian broadening
                y_broadened[i] += lorentz
        return x_new, y_broadened
