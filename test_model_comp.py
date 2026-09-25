import numpy as np
import pandas as pd

topmovers = np.load('datasets/topmovers_IFNG.npy', allow_pickle=True)
essential = pd.read_csv("CEGv2.txt", delimiter='\t')['GENE'].tolist()
topmovers_ne = list(set(topmovers) - set(essential))

sonnet = np.load('sonnet_IFNG/IFNg_test/sampled_genes_5.npy', allow_pickle=True)
haiku = np.load('haiku_IFNG_rerun/test/sampled_genes_5.npy', allow_pickle=True)

for name, pred in [('sonnet', sonnet), ('haiku', haiku)]:
    pred_set = set(pred.tolist())
    hits_all = pred_set & set(topmovers.tolist())
    pred_ne = pred_set - set(essential)
    hits_ne = pred_ne & set(topmovers_ne)
    print(f"{name}: All hits = {len(hits_all)} / {len(topmovers)}, N/E hits = {len(hits_ne)} / {len(topmovers_ne)}")

shared_hits_all = (set(sonnet.tolist()) & set(haiku.tolist())) & set(topmovers.tolist())
print("Shared genes that are hits (All):", len(shared_hits_all))
sonnet_only_hits = (set(sonnet.tolist()) - set(haiku.tolist())) & set(topmovers.tolist())
haiku_only_hits = (set(haiku.tolist()) - set(sonnet.tolist())) & set(topmovers.tolist())
print("Sonnet-only genes that are hits:", len(sonnet_only_hits))
print("Haiku-only genes that are hits:", len(haiku_only_hits))