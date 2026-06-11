#!/usr/bin/env python3
"""Recompute WMB Figure-5A concordance + M1-vs-WMB ephys correlation EXACTLY as the
published figure (reference_plot_wmb_main_figure_original.R), but excluding the FULL set
of 35,493 fine-tune+eval cells instead of the 16,000-cell m1_rna_family list.

Prediction source : reviewer_1_comment1/preps_workflow_ctx_m1_family/<...>/all_predictions.csv
Family mapping     : strict assign_wmb_family (numbered WMB subclasses), matching the R script.
Sanity check       : reproduce the published 16,000-exclusion numbers first.
"""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr, spearmanr

ROOT = Path("/mnt/data/serinharmanci/preps_revision_052726")
COMMENT = ROOT / "reviewer_1_comment1"
PWF = COMMENT / "preps_workflow_ctx_m1_family"
PIPE = ROOT / "analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline"
EXCL_DIR = PIPE / "results/wmb_finetune_exclusion"
FIGS = PIPE / "figures"
META = ROOT / "results/preps_ephys_mouse_geneformer__ft_wmb_isocortex_subclass_2000perCellType__hidden+logits+probs/wmb_predictions/wmb_isocortex_combined_predicted_ephys_embs.csv"
M1_MED = ROOT / "reviewer_1_major_point1_clean/figures/DotPlot_M1_patchseq_real_ephys_CTX_M1_family_matched.medians.csv"
M1_N = ROOT / "reviewer_1_major_point1_clean/figures/DotPlot_M1_patchseq_real_ephys_CTX_M1_family_matched.ncells.csv"
EXCL_16K = EXCL_DIR / "wmb_m1_family_finetune_cell_ids.csv"
EXCL_35K = EXCL_DIR / "wmb_subclass_finetune_eval_cell_ids_35493.csv"

FEATURES = ["Input resistance (MOhm)","Latency (ms)","AP amplitude (mV)","Rheobase (pA)",
            "Sag ratio","Membrane time constant (ms)","AP threshold (mV)",
            "Upstroke-to-downstroke ratio","ISI adaptation index","AP width (ms)"]
FAMILY_ORDER = ["IT","ET","CT","Vip","Sncg","Lamp5","Pvalb","Sst"]
DISPLAY = {"IT":"IT","ET":"ET / PT","CT":"CT","Vip":"Vip","Sncg":"Sncg","Lamp5":"Lamp5","Pvalb":"Pvalb","Sst":"Sst"}
FCOL = {"IT":"#4C78A8","ET":"#F58518","CT":"#54A24B","Vip":"#E45756","Sncg":"#B279A2","Lamp5":"#FF9DA6","Pvalb":"#9D755D","Sst":"#72B7B2"}
CONC_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("conc",["#b10026","#f46d43","#ffffbf","#66bd63","#006837"])


def assign_wmb_family(s):
    if not isinstance(s, str): return None
    if re.match(r"^00[4567] L", s) and "IT CTX Glut" in s: return "IT"
    if re.match(r"^022 L5 ET CTX Glut", s): return "ET"
    if re.match(r"^029 L6b CTX Glut", s) or re.match(r"^030 L6 CT CTX Glut", s): return "CT"
    if re.match(r"^046 Vip Gaba", s): return "Vip"
    if re.match(r"^047 Sncg Gaba", s): return "Sncg"
    if re.match(r"^049 Lamp5 Gaba", s) or re.match(r"^050 Lamp5 Lhx6 Gaba", s): return "Lamp5"
    if re.match(r"^051 Pvalb chandelier Gaba", s) or re.match(r"^052 Pvalb Gaba", s): return "Pvalb"
    if re.match(r"^053 Sst Gaba", s): return "Sst"
    return None


def scale01(mat):
    out = mat.copy()
    for f in mat.columns:
        col = mat[f].astype(float); lo,hi = np.nanmin(col),np.nanmax(col)
        out[f] = (col-lo)/(hi-lo) if hi>lo else 0.5
    return out


def load_excl(p): return set(pd.read_csv(p)["cell_id"].astype(str))

_META = None
def meta():
    global _META
    if _META is None:
        m = pd.read_csv(META, usecols=["cell_id","subclass"])
        m["cell_id"] = m["cell_id"].astype(str)
        m["true_family"] = m["subclass"].map(assign_wmb_family)
        _META = m.dropna(subset=["true_family"]).drop_duplicates("cell_id")[["cell_id","true_family"]]
    return _META


def load_preds(subdirs):
    df = pd.concat([pd.read_csv(PWF/s/"all_predictions.csv") for s in subdirs], ignore_index=True)
    df["cell_id"] = df["cell_id"].astype(str)
    return df


def concordance(df):
    rows=[]
    for thr,name in [(0.0,"p >= 0.0\nall cells"),(0.5,"p >= 0.5\nfiltered")]:
        sub=df[df["RNA_family_prob"]>=thr]
        for fam in FAMILY_ORDER:
            f=sub[sub["true_family"]==fam]
            if len(f)==0: continue
            rows.append({"threshold":name,"family":fam,"n_cells":len(f),
                         "n_correct":int((f["RNA_family"]==fam).sum()),
                         "concordance":100*(f["RNA_family"]==fam).mean()})
    return pd.DataFrame(rows)


def run(tag, exclude_ids):
    ct = load_preds(["wmb_isocortex_1_m1_celltype_exons_preds_exact","wmb_isocortex_2_m1_celltype_exons_preds_exact"])
    ct = ct[~ct["cell_id"].isin(exclude_ids)].merge(meta(), on="cell_id", how="inner")
    conc = concordance(ct)
    allc = conc[conc["threshold"].str.contains("0.0")]
    print(f"\n[{tag}] excluded={len(exclude_ids)}  matched-neuronal validation cells={int(allc['n_cells'].sum())}")
    print(allc[["family","n_cells","concordance"]].to_string(index=False))

    ep = load_preds(["wmb_isocortex_1_m1_patchseq_predictions","wmb_isocortex_2_m1_patchseq_predictions"])
    ep = ep[~ep["cell_id"].isin(exclude_ids)].merge(meta(), on="cell_id", how="inner")
    pred_med = ep.groupby("true_family")[FEATURES].median().reindex(FAMILY_ORDER)
    pred_n = ep.groupby("true_family").size().reindex(FAMILY_ORDER)
    m1 = pd.read_csv(M1_MED).rename(columns={"family":"family","RNA family":"family"})
    fcol = [c for c in m1.columns if c.lower() in ("family","rna family","group")][0]
    m1 = m1.rename(columns={fcol:"family"}).set_index("family")[FEATURES].reindex(FAMILY_ORDER)
    m1_s, pred_s = scale01(m1), scale01(pred_med)
    xs,ys=[],[]
    for fam in FAMILY_ORDER:
        for f in FEATURES:
            x,y=m1_s.loc[fam,f],pred_s.loc[fam,f]
            if np.isfinite(x) and np.isfinite(y): xs.append(x); ys.append(y)
    xs,ys=np.array(xs),np.array(ys)
    pr,sr=pearsonr(xs,ys)[0],spearmanr(xs,ys)[0]
    print(f"[{tag}] ephys corr: Pearson {pr:.3f}, Spearman {sr:.3f}, n={len(xs)}")
    return conc,(pr,sr,len(xs)),m1_s,pred_s,pred_n


def make_figure(conc,corr,m1_s,pred_s,pred_n,m1_n,out):
    fig=plt.figure(figsize=(19,13)); gs=GridSpec(2,2,figure=fig,hspace=0.42,wspace=0.28)
    axA,axB,axC,axD=[fig.add_subplot(gs[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
    thr=["p >= 0.0\nall cells","p >= 0.5\nfiltered"]; fams=list(reversed(FAMILY_ORDER))
    grid=np.full((len(fams),len(thr)),np.nan); lab=np.empty_like(grid,dtype=object)
    for _,r in conc.iterrows():
        i=fams.index(r["family"]); j=thr.index(r["threshold"]); grid[i,j]=r["concordance"]; lab[i,j]=f"{r['concordance']:.1f}%\nn={int(r['n_cells']):,}"
    axA.imshow(grid,cmap=CONC_CMAP,vmin=0,vmax=100,aspect="auto")
    axA.set_xticks(range(len(thr))); axA.set_xticklabels(thr,fontsize=10,fontweight="bold")
    axA.set_yticks(range(len(fams))); axA.set_yticklabels([DISPLAY[f] for f in fams],fontsize=11)
    for i in range(len(fams)):
        for j in range(len(thr)):
            if lab[i,j]: axA.text(j,i,lab[i,j],ha="center",va="center",fontsize=9,color="white",linespacing=0.95,fontweight="bold")
    axA.set_title("WMB concordance (PREPS vs subclass)\nall fine-tune+eval cells excluded",fontsize=13,fontweight="bold")
    pr,sr,n=corr; xs,ys,cs=[],[],[]
    for fam in FAMILY_ORDER:
        for f in FEATURES:
            x,y=m1_s.loc[fam,f],pred_s.loc[fam,f]
            if np.isfinite(x) and np.isfinite(y): xs.append(x); ys.append(y); cs.append(FCOL[fam])
    xs,ys=np.array(xs),np.array(ys)
    axB.plot([0,1],[0,1],"--",color="grey",lw=0.8); m,b=np.polyfit(xs,ys,1); xl=np.array([xs.min(),xs.max()]); axB.plot(xl,m*xl+b,color="black",lw=1.0)
    axB.scatter(xs,ys,c=cs,s=42,alpha=0.9)
    axB.text(0.03,0.97,f"Pearson r = {pr:.2f}\nSpearman \u03c1 = {sr:.2f}\nn = {n}",ha="left",va="top",fontsize=11,transform=axB.transAxes,bbox=dict(boxstyle="round",fc="white",ec="grey",alpha=0.9))
    axB.set_xlim(-.02,1.02); axB.set_ylim(-.02,1.02); axB.set_aspect("equal")
    axB.set_xlabel("M1 measured ephys (scaled 0-1)",fontsize=11); axB.set_ylabel("WMB predicted ephys (scaled 0-1)",fontsize=11)
    axB.set_title("M1 vs WMB family-profile correlation",fontsize=13,fontweight="bold")
    axB.legend(handles=[plt.Line2D([0],[0],marker="o",ls="",mfc=FCOL[f],mec="none",label=DISPLAY[f]) for f in FAMILY_ORDER],fontsize=8,loc="lower right")
    for ax,scaled,nc,title in [(axC,m1_s,m1_n,"M1 Patch-seq measured ephys"),(axD,pred_s,pred_n,"WMB PREPS-predicted ephys (FT excluded)")]:
        yorder=list(reversed(FAMILY_ORDER))
        for yi,fam in enumerate(yorder):
            for xi,feat in enumerate(FEATURES):
                v=scaled.loc[fam,feat]
                if np.isfinite(v): ax.scatter(xi,yi,s=40+v*320,c=[v],cmap="magma",vmin=0,vmax=1,edgecolors="0.25",linewidths=0.3)
        ax.set_xticks(range(len(FEATURES))); ax.set_xticklabels(FEATURES,rotation=42,ha="right",fontsize=9)
        ax.set_yticks(range(len(yorder))); ax.set_yticklabels([f"{DISPLAY[f]} (n={int(nc.get(f,0)):,})" for f in yorder],fontsize=10)
        ax.set_title(title,fontsize=13,fontweight="bold"); ax.set_xlim(-.5,len(FEATURES)-.5); ax.set_ylim(-.6,len(yorder)-.4); ax.grid(color="0.9",lw=0.25)
    for ax,t in [(axA,"A"),(axB,"B"),(axC,"C"),(axD,"D")]:
        ax.annotate(t,xy=(-0.08,1.07),xycoords="axes fraction",fontsize=20,fontweight="bold")
    fig.savefig(FIGS/f"{out}.png",dpi=300,bbox_inches="tight",facecolor="white")
    fig.savefig(FIGS/f"{out}.pdf",bbox_inches="tight",facecolor="white"); plt.close(fig)


def main():
    m1n=pd.read_csv(M1_N); m1n=m1n.rename(columns={m1n.columns[0]:"family",m1n.columns[1]:"n"}).set_index("family")["n"]
    run("OLD 16,000-excl (sanity vs published Fig 5A)", load_excl(EXCL_16K))
    conc,corr,m1_s,pred_s,pred_n = run("NEW 35,493-excl (corrected)", load_excl(EXCL_35K))
    conc.to_csv(EXCL_DIR/"WMB_M1_family_concordance_excluding_35493_finetune_cells.csv",index=False)
    pd.DataFrame([{"comparison":"display_scaled01_flattened","n":corr[2],"pearson_r":corr[0],"spearman_rho":corr[1]}]).to_csv(EXCL_DIR/"M1_vs_WMB_excluding_35493_ephys_correlation_summary.csv",index=False)
    make_figure(conc,corr,m1_s,pred_s,pred_n,m1n,"Figure5_WMB_concordance_ephys_excl35493")
    print("\n[saved] corrected tables + Figure5_WMB_concordance_ephys_excl35493.png/pdf")


if __name__=="__main__":
    raise SystemExit(main())
