"""
Mathematical Engineering - Financial Engineering, FY 2024-2025
Risk Management - Exercise 2: Corporate Bond Portfolio
"""

import datetime as dt
from typing import Union

import numpy as np
import pandas as pd

from bootstrap_py.utils import DayCount, yearfrac
from utilities.ex1_utilities import (
    date_series,
    get_discount_factor_by_zero_rates_linear_interp,
    year_frac_act_x,
)


def bond_cash_flows(
    ref_date: Union[dt.date, pd.Timestamp],
    expiry: Union[dt.date, pd.Timestamp],
    coupon_rate: float,
    coupon_freq: int,
    notional: float = 1.0,
) -> pd.Series:
    """
    Calculate the cash flows of a bond.

    Parameters:
    ref_date (Union[dt.date, pd.Timestamp]): Reference date.
    expiry (Union[dt.date, pd.Timestamp]): Bond's expiry date.
    coupon_rate (float): Coupon rate.
    coupon_freq (int): Coupon frequency in payments per years.
    notional (float): Notional amount.

    Returns:
        pd.Series: Bond cash flows.
    """

    cash_flows_dates = date_series(ref_date, expiry, coupon_freq)[1:]

    yfs = [
        yearfrac(prev, curr, DayCount.EU_30_360)
        for prev, curr in zip(cash_flows_dates[:-1], cash_flows_dates[1:])
    ]
    yfs.insert(0, yearfrac(ref_date, cash_flows_dates[0], DayCount.EU_30_360))
    yfs = np.array(yfs)

    # cash_flows_dates = cash_flows_dates.sort() ERROR!!!!!!!!
    cash_flows_dates.sort()

    # Coupon payments
    cash_flows = pd.Series(
        data=notional * coupon_rate * yfs,
        index=cash_flows_dates,
    )

    # Notional payment
    cash_flows[expiry] += notional

    return cash_flows


def defaultable_bond_dirty_price_from_intensity(
    ref_date: Union[dt.date, pd.Timestamp],
    expiry: Union[dt.date, pd.Timestamp],
    coupon_rate: float,
    coupon_freq: int,
    recovery_rate: float,
    intensity: Union[float, pd.Series],
    discount_factors: pd.Series,
    notional: float = 1.0,
) -> float:
    """
    Calculate the dirty price of a defaultable bond neglecting the recovery of the coupon payments.

    Parameters:
    ref_date (Union[dt.date, pd.Timestamp]): Reference date.
    expiry (Union[dt.date, pd.Timestamp]): Bond's expiry date.
    coupon_rate (float): Coupon rate.
    coupon_freq (int): Coupon frequency in payments a years.
    recovery_rate (float): Recovery rate.
    intensity (Union[float, pd.Series]): Intensity, can be the average intensity (float) or a
        piecewise constant function of time (pd.Series).
    discount_factors (pd.Series): Discount factors.
    notional (float): Notional amount.

    Returns:
        float: Dirty price of the bond.
    """
    fixedDates = date_series(ref_date, expiry, coupon_freq)

    if isinstance(intensity, float):
        prob_survival = np.array(
            [
                np.exp(-intensity * year_frac_act_x(ref_date, date, 365))
                for date in fixedDates
            ]
        )
    elif isinstance(intensity, pd.Series):
        intensity = intensity.reindex(fixedDates[1:]).bfill()

        prob_survival = np.ones(len(intensity.index) + 1)

        for i in range(len(intensity.index)):
            yfs = np.array(
                [
                    year_frac_act_x(fixedDates[j], fixedDates[j + 1], 365)
                    for j in range(i + 1)
                ]
            )

            prob_survival[i + 1] = np.exp(-intensity[: i + 1] @ yfs)
    else:
        raise ValueError("Intensity must be a float or a pd.Series")

    cashFlows = bond_cash_flows(ref_date, expiry, coupon_rate, coupon_freq, notional)

    discount_factors = np.array(
        [
            get_discount_factor_by_zero_rates_linear_interp(  # noqa
                ref_date, date, discount_factors.index, discount_factors.values
            )
            for date in fixedDates[1:]
        ]
    )

    dirtyPrice = sum(
        [
            cf * pdef * dis
            for cf, pdef, dis in zip(cashFlows, prob_survival[1:], discount_factors)
        ]
    )

    dirtyPrice += (
        notional
        * recovery_rate
        * sum(
            [
                discount_factors[i] * (prob_survival[i] - prob_survival[i + 1])
                for i in range(len(discount_factors))
            ]
        )
    )

    return dirtyPrice


def defaultable_bond_dirty_price_from_z_spread(
    ref_date: Union[dt.date, pd.Timestamp],
    expiry: Union[dt.date, pd.Timestamp],
    coupon_rate: float,
    coupon_freq: int,
    z_spread: float,
    discount_factors: pd.Series,
    notional: float = 1.0,
) -> float:
    """
    Calculate the dirty price of a defaultable bond from the Z-spread.

    Parameters:
    ref_date (Union[dt.date, pd.Timestamp]): Reference date.
    expiry (Union[dt.date, pd.Timestamp]): Bond's expiry date.
    coupon_rate (float): Coupon rate.
    coupon_freq (int): Coupon frequency in payments a years.
    z_spread (float): Z-spread.
    discount_factors (pd.Series): Discount factors.
    notional (float): Notional amount.

    Returns:
        float: Dirty price of the bond.
    """

    fixedDates = date_series(ref_date, expiry, coupon_freq)

    cashFlows = bond_cash_flows(ref_date, expiry, coupon_rate, coupon_freq, notional)

    discount_factors = np.array(
        [
            get_discount_factor_by_zero_rates_linear_interp(  # noqa
                ref_date, date, discount_factors.index, discount_factors.values
            )
            for date in fixedDates[1:]
        ]
    )

    B_hat = np.array(
        [
            disc * np.exp(-z_spread * year_frac_act_x(ref_date, date, 365))
            for date, disc in zip(fixedDates[1:], discount_factors)
        ]
    )

    dirtyPrice = sum([cf * disc for cf, disc in zip(cashFlows, B_hat)])
    return dirtyPrice
