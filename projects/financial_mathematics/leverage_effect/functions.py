import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy.stats import norm
from scipy.optimize import curve_fit


# parse the data dir for a given index, return a list of dataframes with stock price data
def load_stocks(index, data_parent_dir):
    data_subdir = os.path.join(data_parent_dir, index)
    dataframes = []
    for filename in os.listdir(data_subdir):
        path = os.path.join(data_subdir, filename)
        name = filename.replace('.csv', '')

        df = pd.read_csv(path, delimiter=',', skiprows=3, header=0, names = ['Date', 'Close','High','Low','Open','Volume'], parse_dates=['Date'], usecols=['Date', 'Close'])
        df.set_index('Date', inplace=True)
        df.rename(columns={"Close": name}, inplace=True)

        dataframes.append(df)
    return dataframes


# concatenate data from all stocks
def concat_and_select(dataframes, min_non_na_fraction_col = 0.85, min_non_na_fraction_row = 0.6):
    df_concat = pd.concat(dataframes, axis=1)
    df_concat.index = df_concat.index.astype(dataframes[0].index.dtype)
    df_concat.sort_index(inplace=True)
    old_shape = df_concat.shape
   
    if min_non_na_fraction_col :
        # only stocks with many sufficiently many time points
        q = min_non_na_fraction_col * df_concat.shape[0]
        df_concat = df_concat.loc[:, df_concat.count() >= q]

    # DATES SELECTION
    if min_non_na_fraction_row:
        # only dates with data from most of the stocks
        threshold = min_non_na_fraction_row * df_concat.shape[1]
        df_concat = df_concat[df_concat.count(axis=1) > threshold] 
    
    new_shape = df_concat.shape

    print(f'% of stocks remaining: {new_shape[1]/old_shape[1]:.2%}')
    print(f'% of dates remaining: {new_shape[0]/old_shape[0]:.2%}')
    print(f'Number of stocks: {new_shape[1]}')
    print(f'Number of dates: {new_shape[0]}')

    return df_concat

# calculate standarized returns form price data, calculate returns of an artificial index 
def calculate_returns(df_concat, type='log', smoothed=False):      #diff, log or pct_change    
    if type=='diff':
        df_diff = df_concat.diff()

    elif type == 'log':
        df_diff = np.log(df_concat).diff()

    elif type == 'pct_change':
        df_diff = df_concat.pct_change(fill_method=None)
  
    # reduce local volatility
    if smoothed:
        rolling_std = df_concat.rolling(window='30D', min_periods=10).std()
        df_diff = df_diff / rolling_std

    # standardize returns
    df_diff = (df_diff - df_diff.mean()) /df_diff.std()
  
    print(f'% of nans: {df_diff.isna().sum().sum()/df_concat.size:.2%}')

    index_returns = df_diff.mean(axis=1)

    return df_diff, index_returns


def LI(df_stocks, index_series, tau_list, gaussianize_I=False):  # influence shifted forward by tau! 
    LI_tau = []

    I = index_series

    if gaussianize_I:
        I = gaussianize(I)

    I2 = I**2

    for tau in tau_list:        

        I2_mean = I2.shift(periods=tau).mean() 

        corr_mean = (I.shift(periods=tau) * I2).mean()

        LI_tau.append(corr_mean/I2_mean)
    return pd.Series(LI_tau, index=tau_list, name='LI_tau')

def Lsigma(df_stocks, index_series, tau_list, gaussianize_I=False):
    # df = df.copy().dropna()
    Lsigma_tau = []

    I = index_series

    if gaussianize_I:
        I = gaussianize(I)
    I2 = I**2

    for tau in tau_list:        

        I2_mean = I2.shift(periods=tau).mean() 

        corr_mean = (I.shift(periods=tau) * (df_stocks**2).mean(axis=1)).mean()

        Lsigma_tau.append(corr_mean/I2_mean)
    return pd.Series(Lsigma_tau, index=tau_list, name='Lsigma_tau')

def rho(df):

    N = df.count(axis=1) # for some dates there are nans
    I = df.mean(axis=1) 

    sigma2 = (df**2).mean(axis=1)  

    numerator = (I**2 * N**2) - (N * sigma2)
    denominator = N * (N - 1) * sigma2

    rho_t = numerator / denominator
    return rho_t


def Lrho(df_stocks, index_series, tau_list, gaussianize_I=False):

    Lrho_tau = []

    I = index_series

    if gaussianize_I:
        I = gaussianize(I)
    I2 = I**2

    rho_vals = rho(df_stocks)

    for tau in tau_list:    
        
        I2_mean = I2.shift(periods=tau).mean() # I2_mean = I2.shift(periods=-tau).mean()

        corr_mean = (I.shift(periods=tau) * rho_vals).mean()

        Lrho_tau.append(corr_mean/I2_mean)
    return pd.Series(Lrho_tau, index=tau_list, name='Lrho_tau')

def calculate_correlation_functions(df_stocks, index_series, gaussianize_I=False):
    tau_list = np.arange(1, 200, 1)

    rho_0 = rho(df_stocks).mean()
    sigma2_0 = (df_stocks**2).mean(axis=1).mean()
    I2_mean = (index_series**2).mean()

    LI_vals = LI(df_stocks, index_series, tau_list, gaussianize_I)
    Lsigma_vals = Lsigma(df_stocks, index_series, tau_list, gaussianize_I)
    Lrho_vals = Lrho(df_stocks, index_series, tau_list, gaussianize_I)

    return tau_list, LI_vals, Lsigma_vals, Lrho_vals, sigma2_0, rho_0, I2_mean




def exp_decay_with_offset(tau, A, tau_c, B):
    return A * np.exp(-tau / tau_c) + B


def plot_correlation_functions_portfolio(data, fig, ax):
    tau_list, LI_vals, Lsigma_vals, Lrho_vals, sigma2_0, rho_0, I2_mean = data
    tau_array = np.array(tau_list)

    # --- Plot L_I ---
    LI_vals.plot(ax=ax[1], label=r'$a(\tau)$', color='tab:green', lw=0.8, alpha=0.6)
    try:
        popt_LI, _ = curve_fit(exp_decay_with_offset, tau_array, LI_vals.values, p0=(-1.0, 1.0, 0.0))
        fit_LI = exp_decay_with_offset(tau_array, *popt_LI)
        ax[1].plot(tau_array, fit_LI, '--', color = 'tab:green', lw=1,)
    except RuntimeError:
        print("Fit failed for LI")

    # --- Plot adjusted L_sigma * rho_0 ---
    adjusted_Lsigma = Lsigma_vals * rho_0
    adjusted_Lsigma.plot(ax=ax[0], label=r'adjusted $b(\tau)$', color='tab:blue', lw=0.7, alpha=0.6)
    try:
        popt_Lsigma, _ = curve_fit(exp_decay_with_offset, tau_array, adjusted_Lsigma.values, p0=(-1.0, 1.0, 0.0))
        fit_Lsigma = exp_decay_with_offset(tau_array, *popt_Lsigma)
        ax[0].plot(tau_array, fit_Lsigma, '--', color = 'tab:blue', lw=1,)
    except RuntimeError:
        print("Fit failed for L_sigma")

    # --- Plot adjusted L_rho * sigma2_0 ---
    adjusted_Lrho = Lrho_vals * sigma2_0
    adjusted_Lrho.plot(ax=ax[0], label=r'adjusted $c(\tau)$', color='tab:orange', lw=0.7, alpha=0.6)
    try:
        popt_Lrho, _ = curve_fit(exp_decay_with_offset, tau_array, adjusted_Lrho.values, p0=(-1.0, 1.0, 0.0))
        fit_Lrho = exp_decay_with_offset(tau_array, *popt_Lrho)
        ax[0].plot(tau_array, fit_Lrho, '--', color='tab:orange', lw=1,)
    except RuntimeError:
        print("Fit failed for L_rho")

    # Axis labels and legends
    ax[0].set_xlabel(r'$\tau$ (days)')
    ax[1].set_xlabel(r'$\tau$ (days)')
    ax[0].set_ylabel('Coefficient value')

    ax[0].legend(loc='lower right')
    ax[1].legend(loc='lower right')


    # Print diagnostic info
    print(f'<I^2> = {I2_mean:.4f}')
    print(f'rho_0*sigma2_0 = {rho_0 * sigma2_0:.4f}')
    print(f'rho_0 = {rho_0:.4f}, sigma2_0 = {sigma2_0:.4f}')


def gaussianize(I):
    """
    Gaussianize an array of index returns I.
    
    Parameters:
    - I: 1D numpy array of index returns.
    
    Returns:
    - I_G: Gaussianized version of I.
    """
    T = len(I)
    
    # Get the ranks (smallest = 1, largest = T)
    ranks = np.argsort(np.argsort(I.values)) + 1  # rank from 1 to T
    
    # Convert ranks to uniform quantiles (in (0,1) range)
    uniform_quantiles = ranks / (T + 1)  # avoid 0 and 1 to keep Φ⁻¹ defined
    
    # Apply the inverse CDF of the standard normal distribution (Gaussianize)
    I_G = norm.ppf(uniform_quantiles)
    
    return pd.Series(I_G, index=I.index)
