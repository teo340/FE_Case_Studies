"""
Mathematical Engineering - Financial Engineering, FY 2024-2025
Risk Management - Exercise 3: Equity Portfolio VaR/ES and Counterparty Risk
"""

from enum import Enum
from typing import Optional, Tuple, Union

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.stats import norm


class OptionType(Enum):
    """
    Types of options.
    """

    CALL = "call"
    PUT = "put"


def black_scholes_option_pricer(
    S: float,
    K: float,
    ttm: float,
    r: float,
    sigma: float,
    d: float,
    option_type: OptionType = OptionType.PUT,
    return_delta_gamma: bool = False,
) -> Union[float, Tuple[float, float, float]]:
    """
    Return the price (and possibly delta and gamma) of an option according to the Black-Scholes
    formula.

    Parameters:
        S (float): Current stock price.
        K (float): Strike price.
        ttm (float): Time to maturity.
        r (float): Risk-free rate.
        sigma (float): Implied volatility.
        d (float): Dividend yield.
        option_type (OptionType, {'put', 'call'}): Option type, default to put.
        return_delta_gamma (bool): If True the option delta and gamma are returned.

    Returns:
        Union[float, Tuple[float, float, float]]: Option price (and possibly delta and gamma).
    """
    # error r+d
    d1 = (np.log(S / K) + (r - d + sigma**2 / 2) * ttm) / (sigma * np.sqrt(ttm))
    # error d1 -
    d2 = d1 - sigma * np.sqrt(ttm)

    if option_type == OptionType.CALL:
        if return_delta_gamma:
            return (
                S * np.exp(-d * ttm) * norm.cdf(d1) - K * np.exp(-r * ttm) * norm.cdf(d2),
                np.exp(-d * ttm) * norm.cdf(d1),
                np.exp(-d * ttm) * norm.pdf(d1) / (S * sigma * ttm ** (0.5)),
            )
        else:
            return S * np.exp(-d * ttm) * norm.cdf(d1) - K * np.exp(-r * ttm) * norm.cdf(
                d2
            )
    elif option_type == OptionType.PUT:
        if return_delta_gamma:
            return (
                K * np.exp(-r * ttm) * norm.cdf(-d2)
                - S * np.exp(-d * ttm) * norm.cdf(-d1),
                -np.exp(-d * ttm) * norm.cdf(-d1),
                np.exp(-d * ttm) * norm.pdf(d1) / (S * sigma * ttm ** (0.5)),  # error
            )
        else:
            return K * np.exp(-r * ttm) * norm.cdf(-d2) - S * np.exp(-d * ttm) * norm.cdf(
                -d1
            )
    else:
        raise ValueError("Invalid option type.")


def principal_component_analysis(
    matrix: NDArray[np.float64],
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Given a matrix, returns the eigenvalues vector and the eigenvectors matrix.
    """

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    # Sorting from greatest to lowest the eigenvalues and the eigenvectors
    sort_indices = eigenvalues.argsort()[::-1]

    return eigenvalues[sort_indices], eigenvectors[:, sort_indices]


def gaussian_var_es(
    mu: pd.Series,
    sigma: pd.DataFrame,
    alpha: float,
    weights: pd.Series,
    ptf_notional: float = 1e6,
    delta: float = 1,
) -> Tuple[float, float]:
    """
    Return VaR and ES computed via Gaussian parametric approach according to the following formulas:
        VaR_{alpha} = delta * mu + sqrt(delta) * sigma * VaR^{std}_{alpha}, where
            VaR^{std}_{alpha} = N^{-1}(alpha) and N is the standard normal cumulative distribution
            function.
        ES_{alpha} = delta * mu + sqrt(delta) * sigma * ES^{std}_{alpha}, where
            ES^{std}_{alpha} = phi(N^{-1}(alpha)) / (1 - alpha) and phi is the standard normal
            probability density function.

    Parameters:
        mu (pd.Series): Series of mean returns.
        sigma (pd.DataFrame): Returns covariance matrix.
        alpha (float): Confidence level.
        weights (pd.Series): Portfolio weights (considered unchanged).
        ptf_notional (float): Portfolio notional, default to 1MM.
        delta (float): Scaling factor, default to 1 i.e. no adjusment is performed.

    Returns:
        Tuple[float, float]: VaR and ES.
    """
    var_std = norm.ppf(alpha)
    es_std = norm.pdf(norm.ppf(alpha)) / (1 - alpha)

    var: float = (
        delta * np.dot(mu, weights)
        + np.sqrt(delta) * np.sqrt(weights.T @ sigma @ weights) * var_std
    )
    es: float = (
        delta * np.dot(mu, weights)
        + np.sqrt(delta) * np.sqrt(weights.T @ sigma @ weights) * es_std
    )

    return ptf_notional * var, ptf_notional * es


def hs_var_es(
    returns: pd.DataFrame,
    alpha: float,
    weights: pd.Series,
    ptf_notional: float = 1e6,
    delta: float = 1,
    lambda_: Optional[float] = None,
) -> Tuple[float, float]:
    """
    Return VaR and ES computed via possibly weighted historical simulation:

    Parameters:
        returns (pd.DataFrame): Returns.
        alpha (float): Confidence level.
        weights (pd.Series): Portfolio weights (considered unchanged).
        ptf_notional (float): Portfolio notional, default to 1MM.
        delta (float): Scaling factor, default to 1 i.e. no adjusment is performed.
        lambda_ (Optional[float]): Decay factor for weighted historical simulation, default to None
            i.e. standard historical simulation is performed.

    Returns:
        Tuple[float, float]: VaR and ES.
    """

    losses = -ptf_notional * returns @ weights
    sorted_losses = np.sort(losses)[::-1]

    if lambda_ is None:
        m: int = np.floor(len(sorted_losses) * (1 - alpha)).astype(int)
        return sorted_losses[m], sorted_losses[:m].mean()
    else:
        n: int = len(sorted_losses)
        C: float = (1 - lambda_) / (1 - lambda_**n)
        W = C * np.array([lambda_ ** (n - i) for i in range(1, n + 1)])
        W_sorted = np.sort(W)[::-1]

        i_star = 1
        while sum(W_sorted[:i_star]) <= 1 - alpha:
            i_star += 1
        i_star -= 1

        return sorted_losses[i_star], (W_sorted[:i_star+1] @ sorted_losses[:i_star+1]) / sum(
            W_sorted[:i_star+1]
        )


def plausility_check(
    returns: pd.DataFrame,
    weights: pd.Series,
    alpha: float,
    ptf_notional: float = 1e6,
    delta: float = 1,
) -> float:
    """
    Perform plausibility check on a portfolio VaR estimating its order of magnitude.

    Parameters:
        returns (pd.DataFrame): Returns.
        weights (pd.Series): Portfolio weights.
        alpha (float): Confidence level.
        ptf_notional (float): Portfolio notional, default to 1MM.
        delta (float): Scaling factor, default to one, i.e. no scaling is performed.

    Returns:
        float: Portfolio VaR order of magnitude.
    """

    sVaR = (
        -ptf_notional
        * weights
        * returns.quantile(q=[1 - alpha, alpha], axis=0).T.abs().sum(axis=1)
        / 2
    )

    return np.sqrt(delta * np.dot(sVaR, np.dot(returns.corr(), sVaR)))
