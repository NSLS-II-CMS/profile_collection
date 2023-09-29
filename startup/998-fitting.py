import numpy as np
from scipy.special import erf
from lmfit import Model


def do_fitting(x, y, *, model_type=None, shift=0.5):
    """
    Do fitting for a peak or an edge. The peak is expected to symmetical. The type
    (peak or edge) is determined automatically based on data or could be explicitly
    set by specifying ``model_type``.

    Parameters
    ----------
    x: iterable
        An array or a list of positions
    y: iterable
        An array or a list of measurements
    model_type: str or None
        None - determine model type based on the number of 'roots',
        'step' - step function, 'peak' - peak.
    shift: float
        Shift applied to the normalized values before finding roots, typically 0.5.

    Returns
    -------
    CEN: float
        Position of the center of the peak or the edge
    FWHM: float
        FWHM of the peak (or similar parameter for the edge). If FWHM is 0 for a peak,
        it means that the scanned range does not contain the full peak. In this case
        the estimate of the center position is selected as 'x' with the largest 'y',
        which is not accurate. In this case, the wider range should be scanned.
    (XMIN, XMAX): type(float)
        The range for positions, which define FWHM of the peak. The range is (None, None)
        for the edge (this could be changed if needed).
    """

    if model_type not in (None, "step", "peak"):
        raise ValueError(f"Unrecognized model type: {model_type!r}")

    x, y = np.array(x), np.array(y)

    if x.ndim != 1:
        raise ValueError(f"Array 'x' must have one dimension: x.ndim={x.ndim}")
    if y.ndim != 1:
        raise ValueError(f"Array 'y' must have one dimension: y.ndim={y.ndim}")
    if x.shape != y.shape:
        raise ValueError(
            f"Arrays 'x' and 'y' have unequal number of elements (x.shape={x.shape}, y.shape={y.shape})"
        )

    # Normalize values first:
    ym = (y - np.min(y)) / (np.max(y) - np.min(y)) - shift  # roots are at Y=0

    CEN, FWHM, XMIN, XMAX = None, None, None, None

    def is_positive(num):
        return True if num > 0 else False

    positive = is_positive(ym[0])
    list_of_roots = []
    for i in range(len(y)):
        current_positive = is_positive(ym[i])
        if current_positive != positive:
            rt = x[i - 1] + (x[i] - x[i - 1]) / (abs(ym[i]) + abs(ym[i - 1])) * abs(ym[i - 1])
            list_of_roots.append(rt)
            positive = not positive

    n_roots = len(list_of_roots)

    if (n_roots >= 2) or (model_type == "peak"):  # Peak
        print(f"Fitting a peak ...")

        nmax = y.argmax()
        xmax = x[nmax]

        root1, root2 = None, None
        for r in list_of_roots:
            if r <= xmax:
                if (root1 is None) or (xmax - r < xmax - root1):
                    root1 = r
            if r > xmax:
                if (root2 is None) or (r - xmax < root2 - xmax):
                    root2 = r

        XMIN = root1 if root1 is not None else x[0]
        XMAX = root2 if root2 is not None else x[-1]

        # Can not find the precise center if the scanned range does not contain the full peak.
        if root1 is None or root2 is None:
            root1, root2 = xmax, xmax

        FWHM = abs(root2 - root1)
        CEN = root1 + 0.5 * (root2 - root1)

    if (n_roots == 1) or (model_type == "step"):  # Step function
        print(f"Fitting a step function ...")
        ym = ym + shift

        def err_func(x, x0, k=2, A=1, base=0):  #### erf fit from Yugang
            return base - A * erf(k * (x - x0))

        mod = Model(err_func)
        x0 = np.mean(x)
        k = 0.1 * (np.max(x) - np.min(x))
        pars = mod.make_params(x0=x0, k=k, A=1.0, base=0.0)
        result = mod.fit(ym, pars, x=x)
        CEN = result.best_values["x0"]
        FWHM = result.best_values["k"]

    return CEN, FWHM, (XMIN, XMAX)
