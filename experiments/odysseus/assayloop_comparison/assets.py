"""Explicit array inputs: no hidden paths, downloads, pickle, or credentials."""
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np

@dataclass
class HistoricalBank:
    genes: list[str]
    residual: np.ndarray
    mean: np.ndarray
    availability: np.ndarray
    scalar_c: float
    query_publication: str
    donor_publications: tuple[str, ...]

    def __post_init__(self):
        self.residual=np.asarray(self.residual,dtype=np.float64)
        self.mean=np.asarray(self.mean,dtype=np.float64)
        self.availability=np.asarray(self.availability,dtype=np.float64)
        if len(set(self.genes))!=len(self.genes) or self.residual.shape!=(len(self.genes),32):
            raise ValueError("Unique vocabulary and 32 historical directions required")
        if self.mean.shape!=(32,) or self.availability.shape!=(len(self.genes),):
            raise ValueError("Mean or availability has wrong shape")
        if not all(np.isfinite(x).all() for x in (self.residual,self.mean,self.availability)):
            raise ValueError("Non-finite bank input")
        if not np.isfinite(self.scalar_c) or self.scalar_c<=0:
            raise ValueError("Prior scale must be finite and positive")
        if np.any((self.availability<0)|(self.availability>1)):
            raise ValueError("Invalid availability")
        if not self.query_publication or not self.donor_publications or self.query_publication in self.donor_publications:
            raise ValueError("Query publication must be excluded from documented donors")

    @classmethod
    def load(cls, path):
        """Load explicit NPZ schema; metadata verifies split-exclusion declaration."""
        with np.load(Path(path),allow_pickle=False) as a:
            return cls(a['genes'].astype(str).tolist(),a['residual'],a['mean'],a['availability'],
                       float(a['scalar_c']),str(a['query_publication'].item()),
                       tuple(a['donor_publications'].astype(str).tolist()))


def build_bank(*, genes, source_weights, donor_mask, donor_hits, donor_publications,
               donor_descriptions, query_description, query_publication):
    """E101/E117 32-D construction for an already-approved, aligned donor matrix.

    Descriptions must use the frozen embedding/normalization convention. Outcome
    matrices contain historical donors only. Query labels have no input channel.
    Exact historical replay prefers the frozen bank (SVD signs/versions matter).
    """
    from scipy.sparse.linalg import svds
    W=np.asarray(source_weights,dtype=np.float64)
    M=np.asarray(donor_mask,dtype=np.float64);Y=np.asarray(donor_hits,dtype=np.float64)
    E=np.asarray(donor_descriptions,dtype=np.float64)
    pubs=np.asarray(donor_publications).astype(str)
    if M.shape!=Y.shape or M.ndim!=2 or M.shape[1]!=len(genes) or W.shape!=(len(genes),10):
        raise ValueError("Donor and source arrays are misaligned")
    if E.ndim!=2 or len(E)!=len(M) or len(pubs)!=len(M):
        raise ValueError("Description/publication rows are misaligned")
    if not np.isin(M,[0,1]).all() or not np.isin(Y,[0,1]).all() or not np.isfinite(E).all():
        raise ValueError("Binary finite historical inputs required")
    keep=pubs!=str(query_publication)
    M,Y,E,pubs=M[keep],Y[keep],E[keep],pubs[keep]
    if min(M.shape)<=32: raise ValueError("More than32 donors and genes required")
    tm=M.sum(0);tmy=(M*Y).sum(0)
    marginal=(tmy+20*tmy.sum()/max(tm.sum(),1))/(tm+20)
    centered=M*(Y-marginal)
    Q,_=np.linalg.qr(W,mode='reduced')
    centered-= (centered@Q)@Q.T
    u,s,vt=svds(centered.T,k=32,random_state=0)
    order=np.argsort(s)[::-1];u,s,vt=u[:,order],s[order],vt[order]
    V=u*np.sqrt(np.maximum(s,0))[None,:]
    norms=np.linalg.norm(V,axis=1);covered=(tm>0)&(norms>1e-8)
    V[covered]/=norms[covered,None];V[~covered]=0
    # Retain the historical float32 storage boundary before downstream float64 use.
    V=V.astype(np.float32).astype(np.float64)
    coords=vt.T*np.sqrt(np.maximum(s,0))[None,:]
    text_map=(E.T@np.linalg.solve(E@E.T+10*np.eye(len(E)),coords)).astype(np.float32)
    desc=np.asarray(query_description,dtype=np.float64)
    if desc.shape!=(E.shape[1],) or not np.isfinite(desc).all():raise ValueError("Invalid query description")
    desc=desc/max(np.linalg.norm(desc),1e-12)
    mean=.1*(desc@text_map.astype(np.float64))
    pub_masks=np.asarray([M[pubs==q].mean(0) for q in sorted(set(pubs))])
    availability=(1+pub_masks.sum(0))/(2+len(pub_masks))
    return HistoricalBank(list(genes),V,mean,availability,float(np.maximum(s,1e-3).mean()),
                          str(query_publication),tuple(sorted(set(pubs))))
