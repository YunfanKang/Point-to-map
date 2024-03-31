def likelyhood_ratio(baseline, observed):
    if observed < 0:
        raise ValueError("The number of observed events is negative")
    if baseline < 0:
        raise ValueError("The number of simulated events is negative")
    if baseline == 0 and observed > 0:
        return 1.0
    if observed > baseline:
        return np.power((observed / baseline), observed) * np.exp(baseline - observed)
    return 1.0
def continuous_homogeneous_possion_likelyhood_test(baseline, observed, total_baseline, total_observed):
    if observed < 0:
        raise ValueError("The number of observed events is negative")
    if baseline < 0:
        raise ValueError("The number of simulated events is negative")
    if observed > baseline:
        return (np.power((observed / baseline), observed) * np.power((total_observed - observed)/ (total_baseline - baseline), total_observed - observed)) / np.power((total_observed/total_baseline), total_observed)
    
    return 1.0
    
    
def expected_counts(study_G_size, G_size, G_counts):
    return G_counts * study_G_size / G_size
def network_scan(G, sub_G):
    #It should be independent of the definition of the search window / subnetwork
    total_Count = 