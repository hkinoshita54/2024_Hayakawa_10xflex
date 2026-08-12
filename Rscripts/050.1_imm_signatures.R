# 1. Sade-Feldman Cell 2018; Guo Nature 2018; MSigDB C7 (CD8 effector up).
Cytotoxic_T = c("Gzma","Gzmb","Gzmk","Prf1","Nkg7","Ifng",
  "Fasl","Klrd1","Ccl5","Il2rb")

# 2. Sade-Feldman Cell 2018; Miller Cell 2019 (TOX-driven exhaustion).
Tcell_Exhaustion = c("Pdcd1","Ctla4","Lag3","Havcr2","Tigit",
  "Tox","Tox2","Nr4a1","Nr4a2","Eomes",
  "Entpd1","Cd160","Cd244")

# 3. Zheng Immunity 2017; MSigDB C7 TREG_UP; FOXP3 regulon signatures.
Treg = c("Foxp3","Il2ra","Ikzf2","Ctla4","Tnfrsf18",
  "Tnfrsf4","Tnfrsf9","Nt5e","Entpd1",
  "Tgfbr1","Tgfbr2","Lag3")

# 4.  Nat Immunol 2019; Im SJ Nat Immunol 2016.
T_Progenitor = c( "Tcf7","Ccr7","Sell","Il7r",
  "Lef1","Bcl2","Cxcr5")

# 5. Victora & Nussenzweig Annu Rev Immunol 2024; MSigDB C7 GC_BCELL_UP.
GCB = c("Aicda","Bcl6","Cd83","Fcer2a","Cxcr4","Cxcr5",
  "Mki67","Mef2b")

# 6. Wu Cell 2020; MSigDB C7 PLASMA_CELL_UP.
Plasma = c("Sdc1","Jchain","Xbp1","Ighm","Igha",
  "Derl3","Prdm1","Irf4")

# 7. Xue Immunity 2014 macrophage spectrum; Zilionis Cell Rep 2019; MSigDB Inflammatory Response.
M1_like = c("Il1b","Tnf","Cxcl9","Cxcl10",
  "Ccl2","Ccl3","Ccl4","Nfkbia",
  "Irf5","Stat1")

# 8. Zilionis Cell Rep 2019; Ponzetta Nat Immunol 2019; MSigDB NEUTROPHIL_UP.
TAM2 = c("Arg1","Mrc1","Retnla","Chil3","Il10",
  "Vegfa","Mmp12","Lgals3","Ccl22")

# 9. Gabrilovich Nat Rev Immunol 2017; broncho-GI tumor scRNA papers; MSigDB MDSC signatures
Mono_MDSC = c("S100a8","S100a9","Ccl2","Lgals3","Il1b",
  "Csf1r","Ccr2","Trem2","Apoe","Lpl")

# 10. Zilionis Cell Rep 2019; Ponzetta Nat Immunol 2019; MSigDB NEUTROPHIL_UP.
TAN = c("S100a8","S100a9","Lcn2","Mmp9","Camp",
  "Cxcr2","Csf3r","Ngp","Elane")

# 11. Brown Cancer Cell 2022; Roberts Nat Immunol 2016.
cDC1 = c("Xcr1","Clec9a","Irf8","Batf3","Cadm1")

# 12. Brown Cancer Cell 2022; MSigDB cDC2_UP.
cDC2 = c("Cd209a","Sirpa","Fcgr2b","Irf4","Clec10a")

# 13. Reizis Nat Rev Immunol 2019; MSigDB PDC_UP.
pDC = c("Siglech","Irf7","Bst2","Gzmb","Tcf4")

# 14. Pauken & Wherry Nat Rev Immunol 2015.
PD1_axis = c("Pdcd1","Cd274","Pdcd1lg2")

# 15. Johnston Nat Immunol 2014.
TIGIT_axis = c("Tigit","Pvr","Nectin2","Nectin3","Nectin4")

# 16. Anderson Immunity 2016
Gal9_TIM3 = c("Lgals9","Havcr2")

# manually created
CXCR2_ligands = c("Cxcl1","Cxcl2","Cxcl3","Cxcl5","Cxcl15")

CCR2_ligands = c("Ccl2","Ccl7","Ccl8")

AntiTumor = c(
  Cytotoxic_T,
  c("Ifng","Cxcl9","Cxcl10"),
  c("Xcr1","Clec9a")
)

ProTumor_Suppression = c(
  TAM2,
  Mono_MDSC,
  Treg,
  c("Tgfb1","Il10","Vegfa")
)

immune_signatures <- list(
  ProTumor_Suppression = ProTumor_Suppression,
  AntiTumor = AntiTumor,
  CCR2_ligands = CCR2_ligands,
  CXCR2_ligands = CXCR2_ligands,
  Gal9_TIM3 = Gal9_TIM3,
  TIGIT_axis = TIGIT_axis,
  PD1_axis = PD1_axis,
  TAN = TAN,
  pDC = pDC,
  cDC2 = cDC2,
  cDC1 = cDC1,
  TAM2 = TAM2,
  M1_like = M1_like,
  Mono_MDSC = Mono_MDSC,
  Plasma = Plasma,
  GCB = GCB,
  Tcell_Exhaustion = Tcell_Exhaustion,
  Cytotoxic_T = Cytotoxic_T,
  Treg = Treg,
  T_Progenitor = T_Progenitor
)

save(immune_signatures, file = "RDSfiles/immune_signatures.RData")

# MDSC manually curated from PMID:32086381---
mdsc_sig <- list(
  MDSC = c("Cd84","Stat1","Stat3","Stat6","Irf1","S100a8","S100a9","Anxa1","Cxcl1","Cxcl2","Cxcr2","Trem1","Arg2"),
  PMN = c("Ly6g","S100a8","S100a9","Csf3r","Cxcr2","Retnlg","Mmp8","Mmp9","Lcn2"),
  Mono = c("Ly6c2","Ccr2","LST1","Fcgr3","Lgals3","Il1b","S100a8","S100a9"),
  Supp = c("Arg1","Arg2","Nos2","Ptgs2","Il10","Tgfb1","Cd274","Cybb","Ncf1","Ncf2")
)
