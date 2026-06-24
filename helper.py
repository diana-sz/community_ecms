def flip_reverse_reactions(model, lb_threshold=-999.0, ub_threshold=0.0, verbose=False):
    """
    Flip reactions in a COBRApy model that only allow reverse flux
    (e.g. lower_bound ~ -1000 and upper_bound = 0) so they become forward reactions.

    Parameters
    ----------
    model : cobra.Model
        The COBRApy model
    lb_threshold : float
        Threshold to detect "very negative" lower bounds (default: -999)
    ub_threshold : float
        Threshold for upper bound (default: 0)
    verbose : bool
        Print flipped reactions

    Returns
    -------
    flipped_rxns : list
        List of reaction IDs that were flipped
    """

    flipped_rxns = []

    for rxn in model.reactions:
        lb = rxn.lower_bound
        ub = rxn.upper_bound

        # Detect reverse-only reactions
        if lb <= lb_threshold and ub == ub_threshold:
            # Flip stoichiometry
            new_stoich = {met: -coeff for met, coeff in rxn.metabolites.items()}
            rxn.add_metabolites(new_stoich, combine=False)

            # Set bounds to forward direction
            rxn.lower_bound = 0.0
            rxn.upper_bound = -lb  # e.g. 1000

            flipped_rxns.append(rxn.id)

            if verbose:
                print(f"Flipped reaction: {rxn.id}")

    return flipped_rxns

