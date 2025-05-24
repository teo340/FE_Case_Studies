import numpy as np
import pandas as pd
import datetime as dt
from typing import Tuple, List

from scipy.integrate import quad

from utilities.ex1_utilities import (
    year_frac_act_x,
)


def sigmaHJM(a: float, sigma: float, u: float, T: float) -> float:
    """
    Hull-white model volatility function.
    a: mean reversion speed
    sigma: level of volatility
    u: current time
    T: maturity time
    """
    return sigma * (1 - np.exp(-a * (T - u))) / a


def affine_trick(
    t: dt.datetime,
    pricing_grid: List[dt.datetime],
    mean_reversion: float,
    sigma: float,
    discount_factors: pd.Series,
) -> Tuple[pd.Series, pd.Series]:
    """
    Affine trick: Exploits the affine structure of the Hull-White model. Output the functions
    A(t, t_j) and C(t, t_j) pre-computed on the pricing grid s.t.
    B(t, t_j) = A(t, t_j) * exp(-C(t, t_j) * x(t)).

    Parameters:
        t (dt.datetime): Date w.r.t. which the computations are made.
        pricing_grid (List[dt.datetime]): Pricing grid.
        mean_reversion (float): Hull-white mean reversion speed.
        sigma (float): Hull-white interest rate volatility.
        discount_factors (pd.Series): Discount factors curve.

    Returns:
        Tuple[pd.Series, pd.Series]: Tuple with the precomputed functions A(t, t_j) and C(t, t_j).
    """

    A = pd.Series(index=pricing_grid)
    C = pd.Series(index=pricing_grid)

    today = discount_factors.index[0]

    for sim_date in pricing_grid:
        if sim_date > t:
            C[sim_date] = (
                1 - np.exp(-mean_reversion * year_frac_act_x(t, sim_date, 365))
            ) / mean_reversion

            fwd_disc = discount_factors[sim_date] / discount_factors[t]

            sigma_int, _ = quad(
                lambda x: sigmaHJM(
                    mean_reversion, sigma, x, year_frac_act_x(today, sim_date, 365)
                )
                ** 2
                - sigmaHJM(mean_reversion, sigma, x, year_frac_act_x(today, t, 365)) ** 2,
                0,
                year_frac_act_x(today, t, 365),
            )

            A[sim_date] = fwd_disc * np.exp(-0.5 * sigma_int)
        else:
            A[sim_date] = 0
            C[sim_date] = 0

    return A, C


def simulateHW(
    g: np.ndarray,
    a: float,
    sigma: float,
    xprev: np.ndarray,
    tprev,
    tnext,
) -> np.ndarray:
    """
    Simulate the Hull-White process. The matrix of gaussian random variables g
    is generated outside of the function and passed as an argument.
    """
    mean = xprev * np.exp(-a * year_frac_act_x(tprev, tnext, 365))
    var = sigma**2 * (1 - np.exp(-2 * a * year_frac_act_x(tprev, tnext, 365))) / (2 * a)
    std = np.sqrt(var)

    xnext = mean + std * g
    return xnext
