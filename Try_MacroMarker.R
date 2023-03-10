

FeaturePlot(PBMC.combined, features = c("Top2a", "Ptk2"), min.cutoff = "q9")


## Mac
# Macrophage-Markers
# https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
## M0
FeaturePlot(PBMC.combined, features = c("Cd68", "Adgre1","Cd14","Csf1r","Ly6c1",
                                        "Cx3cr1","Fcgr1a","Itgam","Mertk"), min.cutoff = "q9", coord.fixed = 1) %>% print()

## M1
# M1 http://biocc.hrbmu.edu.cn/CellMarker/
FeaturePlot(PBMC.combined, features = c("Cd16","Cd32","Cd64","Cd68","Cd80","Cd86"), min.cutoff = "q9", coord.fixed = 1) %>% print()
# M1 https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
FeaturePlot(PBMC.combined, features = c("Marco","Nos2","Tlr2","Cd80","Cd86","Csf2",
                                        "Tnf","Il1b","Il6","Tlr4","Cxcl2","Ifng","Il1r1"), min.cutoff = "q9", coord.fixed = 1) %>% print()
FeaturePlot(PBMC.combined, features = c("Il1a","Il1b","Il6","Nos2","Tlr2","Tlr4","Cd80","Cd86"), min.cutoff = "q9", coord.fixed = 1) %>% print()

##*************************************************************************************************##
## M1
# M1 http://biocc.hrbmu.edu.cn/CellMarker/
# TNF-A https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=21926
FeaturePlot(PBMC.combined, features = c("Ccl2", "Ccl5", "Il18", "Il6", "Nos2",
                                        "Cd274", "Tnf","Il1b",
                                        "Cd16", "Cd32","Cd40",
                                        "Cd197"), min.cutoff = "q9", coord.fixed = 1) %>% print()

## M2
# M2 http://biocc.hrbmu.edu.cn/CellMarker/
# TGFbeta1 https://www.ncbi.nlm.nih.gov/gene/21803
FeaturePlot(PBMC.combined, features = c("Chil3","Clec10a", "Il4", "Mrc1",
                                        "Arg1","Tgfb1","Cd163","Il10"), min.cutoff = "q9", coord.fixed = 1) %>% print()


# ## M2c
# # https://www.intechopen.com/chapters/46464
# FeaturePlot(PBMC.combined, features = c("Cd163","Marco", "Cxcl13"), min.cutoff = "q9", coord.fixed = 1) %>% print()


##****
# https://www.bio-rad-antibodies.com/macrophage-m1-m2-tam-tcr-cd169-cd-markers-antibodies.html?JSESSIONID_STERLING=1D9A2780CAC135290345748A64DEAC50.ecommerce1&evCntryLang=TW-en&cntry=TW&thirdPartyCookieEnabled=true
## M2a
# MHCII H2-Ab1 https://www.ncbi.nlm.nih.gov/gene/14961
# SR Srsf2: https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=20382
# CD206 Mrc1: https://www.ncbi.nlm.nih.gov/gene/17533
# IL-1R II Il1r2: https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=16178
FeaturePlot(PBMC.combined, features = c("Cd163", "H2-Ab1","Srsf2","Mrc1","Tgm2",
                                        "Il1r2","Arg1"), min.cutoff = "q9", coord.fixed = 1) %>% print()

# M2b
FeaturePlot(PBMC.combined, features = c("Cd86", "H2-Ab1"), min.cutoff = "q9", coord.fixed = 1) %>% print()

# M2c
# https://www.intechopen.com/chapters/46464
FeaturePlot(PBMC.combined, features = c("Cd163", "Tlr1","Tlr8","Marco", "Cxcl13"), min.cutoff = "q9",  coord.fixed = 1) %>% print()

## TAM
# Cd3 Cd3e https://www.ncbi.nlm.nih.gov/gene/12501
# Cd31 Pecam1 https://www.ncbi.nlm.nih.gov/gene/18613
# iNOS Nos2 https://www.ncbi.nlm.nih.gov/gene/18126
# VEGF Vegfa https://www.ncbi.nlm.nih.gov/gene/22339
FeaturePlot(PBMC.combined, features = c("Ccl2", "Cd3e","Pecam1","Cd68", "Cd163","Il10",
                                        "Nos2", "Pcna","Vegfa"), min.cutoff = "q9",  coord.fixed = 1) %>% print()
##****
##****
# https://www.intechopen.com/chapters/46464
## M1
FeaturePlot(PBMC.combined, features = c("Stat1","Stat4","Socs3",
                                        "Il10","Il12a","Il23a",
                                        "Tnf","Il6","Il1b","Il18"), min.cutoff = "q9", coord.fixed = 1) %>% print()
FeaturePlot(PBMC.combined, features = c("Ccl2","Ccl3","Ccl4","Ccl5","Ccl11","Ccl7","Ccl22",
                                        "Cxcl1","Cxcl2","Cxcl3","Cxcl5","Cxcl8",
                                        "Cxcl9","Cxcl10","Cxcl11","Cxcl16","Nos2"), min.cutoff = "q9", coord.fixed = 1) %>% print()
## M2a
## MR Nr3c2 https://www.ncbi.nlm.nih.gov/gene/110784
FeaturePlot(PBMC.combined, features = c("Stat3","Il10","Il12a","Il23a","Il1ra",
                                        "Tgfb1","Ccl17","Ccl18","Ccl22","Ccl24",
                                        "Srsf2","Nr3c2","Nos2"), min.cutoff = "q9", coord.fixed = 1) %>% print()

## M2b
FeaturePlot(PBMC.combined, features = c("Il10","Il12a","Il23a","Tnf","Il6","Il1b",
                                        "Ccl1","Nos2"), min.cutoff = "q9", coord.fixed = 1) %>% print()

## M2c
FeaturePlot(PBMC.combined, features = c("Stat6","Il10","Il12a","Il23a","Tgfb1",
                                        "Ccl16","Ccl18","Cxcl13",
                                        "Nr3c2","Cd163"), min.cutoff = "q9", coord.fixed = 1) %>% print()


##*************************************************************************************************##


# M2 http://biocc.hrbmu.edu.cn/CellMarker/
FeaturePlot(PBMC.combined, features = c("Chil3","Csf1r","Mrc1","Pparg","Arg1","Cd163","Clec10a","Clec7a",
                                        "Cd206","Cd209","Ccl18","Fizz1"), min.cutoff = "q9", coord.fixed = 1) %>% print()
# https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
FeaturePlot(PBMC.combined, features = c("Cd115", "Cd206", "Pparg", "Arg1", "Cd163", "Cd301",
                                        "Dectin-1", "Pdcd1lg2", "Fizz1"), min.cutoff = "q9", coord.fixed = 1) %>% print()

FeaturePlot(PBMC.combined, features = c("Chil3"), min.cutoff = "q9", coord.fixed = 1) %>% print()


