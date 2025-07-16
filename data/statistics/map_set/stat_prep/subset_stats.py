
import pickle
from dpluspy import bootstrapping


with open('Behrer_stats.pkl', "rb") as fin:
    stats = pickle.load(fin)
substats = bootstrapping.subset_stats(
    stats, min_r=1e-6, max_r=1e-2, return_dict=True)
with open("BehrerStats16bins.pkl", "wb") as fout:
    pickle.dump(substats, fout)

