from numpy import exp, sqrt
import squmfit

@squmfit.model 
def three_dim_diffusion(lag, tauD, aspect, n, alpha=1):
    """
    Model of correlation function for a three dimensional diffuser
    
    :param lag: Lag time
    :param tauD: Diffusion time
    :param aspect: Aspect ratio
    :param n: Concentration
    :param alpha: Anomalous diffusion exponent
    """
    tau_taud = (lag / tauD)**alpha
    return 1. / sqrt(1 + tau_taud * aspect**-2) / (1 + tau_taud) / n

@squmfit.model
def triplet_correction(lag, tripletFrac, tauF):
    """
    Multiplicative correction to diffusion model accounting for
    triplet transitions of fluorophore.

    :param lag: Lag time
    :param tripletFrac: Fraction of species inhabiting triplet state
    :param tauF: Triplet lifetime (typically on the order of microseconds).
    """
    tau_tauf = lag / tauF
    return (1 - tripletFrac + tripletFrac * exp(-tau_tauf)) / (1 - tripletFrac)

