from MDAnalysis import Universe
import MDAnalysis.analysis.encore as encore

ens1 = Universe("results/full_bayesian/calmodulin/saxs/saxs_ensemble.pdb", multiframe=True)
ens2 = Universe("results/full_bayesian/calmodulin/saxs_CS/saxs_CS_ensemble.pdb", multiframe=True)
ens3 = Universe("results/full_bayesian/calmodulin/saxs_EP/saxs_EP_ensemble.pdb", multiframe=True)
ens4 = Universe("results/full_bayesian/calmodulin/saxs_CS_EP/saxs_EP_CS_ensemble.pdb", multiframe=True)

#rmsd_matrix = encore.get_distance_matrix(encore.utils.merge_universes([ens1, ens2, ens3, ens4]))
#weights = [[0.21,0.05,0.16,0.43,0.16],[0.18,0.19,0.44,0.04,0.18],[0.53,0.42,0.05],[0.26,0.26,0.16,0.32]]
print (encore.ces([ens1, ens2, ens3, ens4]))