"""
Post distibution statistics
Chi2 and Jensen-Shanon Divergence caluclators
"""
import numpy as np

def calculateChiCrysol(weightedIns, expIns, expErr):
        """
        Calculates chis same way as it is done in crysol
        """
        #Calculate scaling factor
        chi2_=0.0
        mixed_term_ = 0.0
        square_calc_ = 0.0

        Sindex = 0
        for ins in expIns:
            Iobs=ins
            square_err_ = expErr[Sindex]*expErr[Sindex]
            mixed_term_ += Iobs*weightedIns[Sindex]/square_err_
            square_calc_ += weightedIns[Sindex]*weightedIns[Sindex]/square_err_
            Sindex+=1
        scale_factor = mixed_term_/square_calc_

        Sindex = 0
        for ins in expIns:
            Iobs=ins
            square_err_ = expErr[Sindex]*expErr[Sindex]
            chi2_+=(Iobs-scale_factor*weightedIns[Sindex])*(Iobs-scale_factor*weightedIns[Sindex])/square_err_
            Sindex+=1
        chi2_=chi2_/Sindex
        return chi2_

def calculateChemShiftsChi(weightedIns, expIns, expErr, teoErr):
        """
        Calculates chis for chemical shifts
        """
        #Calculate scaling factor
        chi2_ = 0.0
        square_calc_ = 0.0

        #scale_factor = 1
        Sindex = 0
        for ins in expIns:
            Iobs=ins
            et_err_ = np.sqrt(expErr[Sindex]*expErr[Sindex]+teoErr[Sindex]*teoErr[Sindex])
            square_calc_ += (Iobs-weightedIns[Sindex])*(Iobs-weightedIns[Sindex])
            chi2_+=(Iobs-weightedIns[Sindex])*(Iobs-weightedIns[Sindex])/(et_err_*et_err_)
            Sindex+=1
        return chi2_/Sindex


def InfEntropy(weight):
    S = 0.0
    for wi in weight:
        if wi<0.0000001:
                wi = 0.0000001
        S-=wi*np.log2(wi)
    return S

def JensenShannonDiv(weights_a, weights_b):
    jsd = InfEntropy((weights_a+weights_b)*0.5)-0.5*InfEntropy(weights_a)-0.5*InfEntropy(weights_b)
    return jsd

def mean_for_weights(fit):
    """
    Calculation mean weight from stan fitting object
    :param fit: stan fit model
    :return:
    """
    stan_chain = fit.extract()
    weights = stan_chain["weights"]
    vlen = np.shape(weights)
    mean_weights = []
    for ind in range(vlen[1]):
        mean_weights.append(np.mean(weights[:,ind]))
    return mean_weights

def me_log_lik(fit):
    """
    Calculates model evidence from log likelihood calculations
    :param log_lik:
    :return:
    """
    stan_chain = fit.extract()
    log_lik = stan_chain["loglikes"]
    return np.sum(np.mean(log_lik, axis=1))


def waic(log_lik):
    """
    Calculate the widely available information criterion
    """
    lppd =  np.sum(np.log(np.mean(np.exp(log_lik), axis=0)))

    p_waic = np.sum(np.var(log_lik, axis=0))

    return -2 * lppd + 2 * p_waic