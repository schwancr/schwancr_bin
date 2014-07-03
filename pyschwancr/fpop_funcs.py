
import copy
import numpy as np

RATE_CONSTANTS = {
'cys' : 3.4E10,
'tyr' : 1.3E10,
'his' : 1.3E10,
'met' : 8.3E9,
'phe' : 6.5E9,
'arg' : 3.5E9, # note this may not be a plus n*16 since it often breaks off the nitrogen group on the end
'leu' : 1.7E9,
'ile' : 0, # Not available in paper
'trp' : 1.3E9,
'pro' : 0, # not available in paper
'val' : 7.6E8,
'thr' : 5.1E8,
'ser' : 3.2E8,
'glu' : 2.3E8,
'gln' : 0, # not available in paper
'ala' : 7.7E7,
'asp' : 7.5E7,
'asn' : 4.9E7,
'lys' : 3.5E7,
'gly' : 1.7E7
}

# ^^^^ These rate constants are from Xu et al. (2003) Anal. Chem.
#  (DOI: 10.1021/ac035104h)

RATE_CONSTANTS_INTERP = copy.copy(RATE_CONSTANTS)
RATE_CONSTANTS_INTERP['ile'] = 1E9
RATE_CONSTANTS_INTERP['pro'] = 1E8
RATE_CONSTANTS_INTERP['gln'] = 1E8
# ^^^ Order of magnitude interpolation based on my intuition
# ALL RATE CONSTANTS ARE IN M^-1 s^-1

def _mod_prob(protein_conc, OH_conc, lagtime, rate):
    """
    compute the probability of a residue being modified in lagtime seconds
    
    assuming that OH_conc is constant over the lagtime

    NOTE: This is not necessarily the *right* way to simulate FPOP depending
        on the lifetime of the OH radical
    """

    exp_temp = np.exp((protein_conc - OH_conc) * rate * lagtime)
    
    prob = (1 - exp_temp) / (1 - protein_conc / OH_conc * exp_temp)

    return prob


def cut30(sasa, res_names):

    prob = np.zeros(sasa.shape)
    
    prob[np.where(sasa > 0.50)] = 1

    return prob


def percent_of_ref(sasa, res_names):
    
    ref_sasa_dict = {
    'ala' : 0.67,
    'arg' : 1.96,
    'asn' : 1.13,
    'asp' : 1.06,
    'cys' : 1.04,
    'gln' : 1.44,
    'glu' : 1.38,
    'gly' :   -1,  # label to remove later
    'his' : 1.51,
    'ile' : 1.40,
    'leu' : 1.37, 
    'lys' : 1.67,
    'met' : 1.60,
    'phe' : 1.75,
    'pro' : 1.05,
    'ser' : 0.80,
    'thr' : 1.02,
    'trp' : 2.17,
    'tyr' : 1.87,
    'val' : 1.17
    }

    ref_sasa_per_res = []
    for name in res_names:
        ref_sasa_per_res.append(ref_sasa_dict[name])
    
    ref_sasa_per_res = np.array(ref_sasa_per_res)

    perc_sasa = sasa / np.reshape(ref_sasa_per_res, (1, -1))

    perc_sasa[:, np.where(ref_sasa_per_res < 0)[0]] = 0
    # for glycine, need to set to zero

    return perc_sasa
    

def residue_probabilities(protein_conc, OH_conc, lagtime, pdb, interpolate=False):
    """
    get per residue probabilty of being labeled in some time lagtime.
    
    Parameters:
    -----------
    protein_conc : float
        concentration of protein (in M)
    OH_conc : float
        concentration of protein (in M)
    lagtime : float
        lag time of MSM in seconds
    pdb : msmbuilder.Trajectory
        Trajectory object to get the list of residues

    Output:
    -------
    probs : np.ndarray
        list of probabilities (one for each residue in pdb)
    """

    rates = []

    if interpolate:
        rate_dict = RATE_CONSTANTS_INTERP
    else:
        rate_dict = RATE_CONSTANTS

    for i in np.unique(pdb['ResidueID']):
        name = pdb['ResidueNames'][np.where(pdb['ResidueID'] == i)][0].lower()
        if not name in rate_dict.keys():
            raise Exception("non-standard residue name %s in pdb. rename and try again." % name)

        rates.append(rate_dict[name])

    rates = np.array(rates)

    probs = _mod_prob(protein_conc, OH_conc, lagtime, rates)

    return probs


def products_prob(protein_conc, OH_conc, lagtime, pdb, avg_sasas, sasa_func=cut30, 
        num_products=5, interpolate=False):
    """
    get the probability of multi-label products for an fpop experiment per state
    
    Parameters:
    -----------
    protein_conc : float
        concentration of protein (in M)
    OH_conc : float
        concentration of protein (in M)
    lagtime : float
        lag time of MSM in seconds
    pdb : msmbuilder.Trajectory
        Trajectory object to get the list of residues
    avg_sasas : np.ndarray (2D)
        shape should be (num_states, num_residues) with each row corresponding
        to a state in the msm, and each column corresponding to the average sasa
        of that residue
    sasa_func : function, optional
        function that maps sasa to a number in [0, 1] corresponding to the 
        probability of a residue being modified. default: fpop.cut30 which
        assigns probability=1 for residues with greater than 0.30 nm^2 SASA
    num_products : int, optional
        number of products to calculate the probabilities for (default: 5)
    interpolate : bool, optional
        pass True to use Christian Schwantes' guesses at GLN, PRO, and ILE
        otherwise they have zero probability of being labeled (default: False)

    Outputs:
    --------
    prod_probs : np.ndarray (2d)
        each row (i) corresponds to a state and each column (j) is the probability
        of labeling the protein j times in lagtime seconds while in state i
    """

    num_states = len(avg_sasas)

    res_names = []
    for i in np.unique(pdb['ResidueID']):
        res_names.append(pdb['ResidueNames'][np.where(pdb['ResidueID'] == i)][0].lower())

    per_res_probs = residue_probabilities(protein_conc, OH_conc, lagtime, pdb, interpolate=interpolate)
    print per_res_probs

    state_sasa_probs = sasa_func(avg_sasas, res_names)
    print np.where(state_sasa_probs == 0)[0].shape

    res_probs = state_sasa_probs * np.reshape(per_res_probs, (1, -1))

    if np.any(res_probs >= 1):
        raise Exception("unknown error, a residue probability is greater than or equal to one.")

    # Each row is a state, and each column is a residue's probability of being modified
    # given it is in state i
    #    p[i,j] = sasa_func(avg_sasas[i, j]) * res_prob[j]

    # One can show that P( T = M ) = \sum_i p_i / (1 - p_i) * P( T = M-1 )
    # where T is a random variable corresponding to the number of labelled
    # residues and p_i are the individual probabilities of labeling
    # residue i

    prod_probs = np.zeros((num_states, num_products + 1))  # we also include zero labels

    not_res_probs = 1 - res_probs
    ratio_probs = res_probs / not_res_probs
    ratio_probs_sum = ratio_probs.sum(axis=1)
    print ratio_probs_sum

    prod_probs[:, 0] = np.product(not_res_probs, axis=1)

    for i in xrange(1, num_products + 1):
        prod_probs[:, i] = ratio_probs_sum * prod_probs[:, i - 1]

    return prod_probs
