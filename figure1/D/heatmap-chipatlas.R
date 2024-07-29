#############################
## This script builds different heatmaps based on the chip-atlas results.
##
## The list of sra to filter was determined manually looking at the sra records.
## See the below and the function .buildreplacementvectors() for details.
##
##  resultlistescsramerged[[1]][[1]] - mHG13hist (ok SRX16815311, SRX16815319,
## SRX16815317, SRX19148013, SRX12734384, SRX747752, SRX12735104, SRX5637724,
## SRX1775637, SRX111870, SRX14812477):
## H4 - SRX16815323: Brg1-depleted auxin 20h replicate 2
## H3K4me3 - SRX10304465: treated with WNT3a protein for 12h
## H3K4me2 - SRX4681571:  H3K4me2_dm_RA_d6; Mus musculus; ChIP-Seq
## H3.3 - SRX2245482: H3.3-HA ChIP-seq ES cells 1hr DOX; Mus musculus; ChIP-Seq #nolint
## H3K64ac - SRX1560890: Cell line 46c (sox1_GFP)
## H3K122ac - SRX1560888: Cell line 46c (sox1_GFP)
## H3K79me2 - SRX13332952: ES Bruce 4
## H3K36me3 - SRX12382446: Day 1 SETD8 KO
## H2A.Z - SRX2955747: J1 mouse ESC
## H3K9ac - SRX185848: ES Bruce 4
## H4ac - SRX10445814: SUZ12 KO Day 8
## H4K5ac - SRX14210149: SET1A SET domain deletion
## H4K20ac - SRX13332963: SRM200254_Dnttip1KO1E6_H4K20Ac-AB2
## H4K16ac - SRX298193: Cell line 46c (sox1_GFP)
## H3K9K14ac - ERX245530: Cfp1 KO
## H3K18ac - SRX10424671: Trim66 KO
##
##  resultlistescsramerged[[1]][[2]] - mHG13polTFs (ok SRX8120163, SRX5251389,
## SRX16815350, SRX9446557, SRX6831867, SRX994816, SRX2894853, SRX8556272,
## SRX6827273, SRX14663256,SRX1540623, SRX4580513, SRX1036419, SRX9596422,
## SRX2894851, SRX1670201, SRX2270102, SRX482054, SRX4167136, SRX017058,
## SRX2539827, SRX868619, SRX518333, SRX4167125, SRX5251387, SRX188819,
## SRX2754629, SRX1036420, SRX6879007, SRX1080398, SRX10111410,SRX054580,
## SRX104410, SRX803857, SRX17708154, SRX14134858, SRX6078354, SRX2894852,
## SRX4167134, SRX5247639, SRX3341233, SRX518335, SRX4232696,SRX9637872,
## SRX518334, SRX8330923, SRX1741141, SRX9345864, SRX19935223, SRX6879008,
## SRX823725, SRX1639485, SRX1602644, SRX5436003, SRX759626, SRX14059473,
## SRX236482, SRX6784013, SRX7622740, SRX1829032, SRX142127, SRX104386,
## SRX6033283, SRX039347, SRX14672122, SRX997153, SRX7423745, SRX8722031,
## SRX14210166, SRX6831864, SRX368296, SRX248284, SRX4506772):
## Nanog - SRX9486380:  Nanog ChIP-nexus with CRISPR-knockout of Sox2 binding site #nolint
## Sp3 - SRX5371695: A17Lox_Sp1het1_ESC_Sp3_ChIPseq
## Kmt2d - SRX2762171: WT_RA_ChIP_MLL4; Mus musculus; ChIP-Seq
## Sp1 - SRX5371690: A17Lox_Sp1het2_ESC_Sp1_ChIPseq
## Tbp - SRX11677010: Pol III_degron mESCs IAA 1hr
## Utf1 - SRX404068: KH2 mouse ESC
## Taf12 - SRX11221956: TFIIA degron 0h
## Setd1a - SRX5600782: mES_WT_shNT_Set1aChIP; Mus musculus; ChIP-Seq
## Sox2 - SRX9486383: Sox2 ChIP-nexus with CRISPR-knockout of Nanog binding site #nolint
## RNA polymerase II - SRX4178816: Cell line R1
## Yy1 - SRX5850150: TKO.TSA_SP1_REP3
## Epitope tags - SRX2486173: CFP1 fl/fl cell line, with 3xT7-Strep2 tag knockin at SET1A N-ter #nolint
## Sin3a - SRX9450269: Genotype Ainv15(iCre)
## Smad2  - SRX14663260: PML KO-SMAD2/3 rep2; Mus musculus; ChIP-Seq
## Taf3 - SRX13439109: ChIP-Seq_TAF3_BRD2_degron_3h_rep1
## Ssrp1 - SRX11677018: treatment Pol III_degron mESCs IAA 1hr
## Esrrb - SRX11221882: TAF12 degron 0h
## GFP - SRX385377: E14_Mll2FL_C_GFP_FCS_LIF_2i_double_crosslink
## Taf1 - SRX11221945: TAF12_degron_0h_TAF1_ChIP_Rep2; Mus musculus; ChIP-Seq
## Dendra2 - SRX4178805:Dendra2-RPB1, Halo-MED19, R1 ESC
## Otx2  - SRX9450552: D4_Ctr_OTX2_2; Mus musculus; ChIP-Seq
## Med1 - SRX11221877: TAF12_degron_0h_MED1_ChIP_Rep1
## Banp - SRX10157266: deleted by the GEO staff
## Ints3 - SRX11677013: ChIP_Seq_Pol_III_IAA_anti_INTS3_rep1, treatment Pol III_degron mESCs IAA 1hr #nolint
## Med4 - SRX11221901: MED4_degron_0h_MED6_ChIP_Rep1; Mus musculus; ChIP-Seq
## Kdm6a - SRX13332957: SRM220498_Dnttip1KO1E6_UTX-AB11; Mus musculus; ChIP-Seq
## Brd4 - SRX3898501: Hexanediol treatment
## Kmt2c - SRX2752337: genotype Mll3/4 dCD, R1 ESC
## Mettl3 - SRX2415040: genotype/variation Rosa26:BirA, METTL3 FLAG/Avi
## Kdm4c - SRX213793: Phenotype Agouti
## Med24 - SRX5926394: Treatment 4OHT 2h
## Nelfe - SRX9450256: Genotype Ainv15(iCre)
## Tcf12 - SRX390394: Strain CGR8
## Ctcf - SRX10437500:  F-ALL-INV genotype
## Brd9 - SRX3751909: DMSO treatment
## Rnf2 - SRX482879: RING1A-B_cKO_RING1B_treated_rep2
## Rpa1 - SRX11676997: Pol III_degron mESCs IAA 1hr
## Nipbl - SRX15652647: NIPBL ChIP-M of wt stembryos at 168h
## Rad21 - SRX1670207: Mll3-4 DKO Day1
## Brca2 - SRX15464337: Zfp281 KO
## Med14 - SRX12622651: Treatment tamoxifen 96h
## Ncaph2 - SRX1688321: Transfection shGFP
## Hif1a - SRX13945545: Hypoxia day 6
## DNA-RNA hybrid - SRX13439146:  Rloop_BRD2_BRD3_BRD3_KO_rep1
## Brd3 - SRX13439086: hIP-Seq_BRD3_BRD2_degron_3h_rep1
## Ep400 - SRX9424475: Treatment Auxin 1h
## Max - SRX5850126: TKO TSA Max Rep3
## Ercc3 - SRX11221923: MED4_degron_0h_XPB_ChIP_Rep2
## Tet3 - SRX2545592: Tet3-CXXC
## Dnttip1 - SRX13332921: SRM166780_Bruce4_rep2_DNTTIP1-AB1
## Fgfr1 - SRX868181: Retinoic acid treatment for 48h
## Gtf3c1 - SRX1688323: shGFP_TFIIIC_1
## Klf5 - SRX334975: Cell line GCR8
## Gabpa - SRX5850112: Treatment 5 ng microl trichostatin A for 36h
## Biotin - SRX2388542: Gata3 overexpression
## Pcgf1 - SRX6071847: Treatment tamoxifen 72h
## Maz - SRX13440032: genotype MAZ HoxA delta 5-6
## Runx1 - SRX8754229: R201Q no dox
## Zfp42 - SRX128354: Expressing V5-tagged REX1
## Ctnnb1 - SRX9101843: 3xFlag knock-in to endogeous Kdm2b C-terminus
## Gata4 - SRX191946: Dnmt3a-3b DKO
## Sumo1 - SRX11525492: Cell line R1 Day 8
## 0610010K14Rik - SRX19028417: ChIP_BAP18_dTAG-24h_mRBBP5-dTAG_ES_Rep2
## Zic2 - SRX5817869: 32h N2B27
## Ints5 - SRX19028478: ChIP_INTS5_dTAG-24h_mRBBP5-dTAG_ES_Rep2
## Hoxb1 - SRX4980239: 24h post-RA induction
## Kdm2b - SRX186548: Treatment shRNA
## E2f6 - SRX3581849: V6.5_shGFP_E2f6_StaphSeq_rep2
##
##  resultlistescsramerged[[2]][[1]] - mHG14hist (ok SRX16815311, SRX16815319,
## SRX16815317, SRX19148013, SRX12734384, SRX747752, SRX12735104, SRX175895,
## SRX14812477, SRX1775637)
## H3K4me3 - SRX10304465: treated with WNT3a protein for 12h
## H4 - SRX16815323: Brg1-depleted auxin 20h replicate 2
## H3K79me2 - SRX13332952: ES Bruce 4
## H3K4me2 - SRX4681571:  H3K4me2_dm_RA_d6; Mus musculus; ChIP-Seq
## H3K64ac - SRX1560890: Cell line 46c (sox1_GFP)
## H3K122ac - SRX1560888: Cell line 46c (sox1_GFP)
## H3.3 - SRX2245482: H3.3-HA ChIP-seq ES cells 1hr DOX; Mus musculus; ChIP-Seq #nolint
## H3K9ac - SRX185848: ES Bruce 4
## H2A.Z - SRX2955747: J1 mouse ESC
## H3K36me3 - SRX12382446: Day 1 SETD8 KO
## H4ac - SRX10445814: SUZ12 KO Day 8
## H4K5ac - SRX14210149: SET1A SET domain deletion
## H4K20ac - SRX13332963: SRM200254_Dnttip1KO1E6_H4K20Ac-AB2
## H4K16ac - SRX298193: Cell line 46c (sox1_GFP)
## H3K9K14ac - ERX245530: Cfp1 KO
## H3K18ac - SRX10424671: Trim66 KO
##
## resultlistescsramerged[[2]][[2]] - mHG14polTFs (ok SRX8120163, SRX5251389,
## SRX6831867, SRX16815349, SRX9446557, SRX6827273, SRX994816, SRX8556272,
## SRX11677012, SRX2894853, SRX1540623, SRX4580513, SRX047141, SRX1036419,
## SRX14663256, SRX2270102, SRX9596422, SRX4167136, SRX1670201, SRX017058,
## SRX2539827, SRX868619, SRX868619, SRX2894851, SRX4167125, SRX5251387,
## SRX482054, SRX188819, SRX803857, SRX1080398, SRX518333, SRX2754630,
## SRX054580, SRX10111410, SRX1036420, SRX5247639, SRX6879007, SRX4232696,
## SRX6879008, SRX14134858, SRX4167134, SRX6078354, SRX8330923, SRX2894852,
## SRX17708154, SRX17708154, SRX3341233, SRX823725, SRX9637872, SRX9345864,
## SRX6033283, SRX104410, SRX518335, SRX518334, SRX14672122, SRX039347,
## SRX1517381, SRX236482, SRX1741141, SRX104386, SRX7622740, SRX14059473,
## SRX7622743, SRX1602644, SRX368296, SRX997153, SRX14210166, SRX6784013,
## SRX6831864):
## Nanog - SRX9486380:  Nanog ChIP-nexus with CRISPR-knockout of Sox2 binding site #nolint
## Kmt2d - SRX2762171: WT_RA_ChIP_MLL4; Mus musculus; ChIP-Seq
## Sp3 - SRX5371695: A17Lox_Sp1het1_ESC_Sp3_ChIPseq
## Tbp - SRX11677010: Pol III_degron mESCs IAA 1hr
## RNA polymerase II - SRX19391748: ChIP-Seq_RPB1_RPB3_IAA_0h_rep1;
## Sox2 - SRX9486383: Sox2 ChIP-nexus with CRISPR-knockout of Nanog binding site #nolint
## Setd1a - SRX5600782: mES_WT_shNT_Set1aChIP; Mus musculus; ChIP-Seq
## Sp1 - SRX5371688: cell line A17Lox
## Utf1 - SRX404068: KH2 mouse ESC
## Taf12 - SRX11221956: TFIIA degron 0h
## Sin3a - SRX9450269: Genotype Ainv15(iCre)
## Taf3 - SRX13439109: ChIP-Seq_TAF3_BRD2_degron_3h_rep1
## Ssrp1 - SRX11677018: treatment Pol III_degron mESCs IAA 1hr
## GFP - SRX385377: E14_Mll2FL_C_GFP_FCS_LIF_2i_double_crosslink
## Esrrb - SRX11221882: TAF12 degron 0h
## Smad2  - SRX14663260: PML KO-SMAD2/3 rep2; Mus musculus; ChIP-Seq
## Taf1 - SRX11221945: TAF12_degron_0h_TAF1_ChIP_Rep2; Mus musculus; ChIP-Seq
## Epitope tags - SRX2486173: CFP1 fl/fl cell line, with 3xT7-Strep2 tag knockin at SET1A N-ter #nolint
## Dendra2 - SRX4178805:Dendra2-RPB1, Halo-MED19, R1 ESC
## Banp - SRX10157266: deleted by the GEO staff
## Med4 - SRX11221901: MED4_degron_0h_MED6_ChIP_Rep1; Mus musculus; ChIP-Seq
## Med1 - SRX11221877: TAF12_degron_0h_MED1_ChIP_Rep1
## Yy1 - SRX5850150: TKO.TSA_SP1_REP3
## Brd4 - SRX3898501: Hexanediol treatment
## Otx2  - SRX9450552: D4_Ctr_OTX2_2; Mus musculus; ChIP-Seq
## Kdm6a - SRX13332957: SRM220498_Dnttip1KO1E6_UTX-AB11; Mus musculus; ChIP-Seq
## Med24 - SRX5926395: Treatment 4OHT 2h
## Kmt2c - SRX2752337: genotype Mll3/4 dCD, R1 ESC
## Tcf12 - SRX390394: Strain CGR8
## Nelfe - SRX9450256: Genotype Ainv15(iCre)
## Ctcf - SRX10437500:  F-ALL-INV genotype
## Kdm4c - SRX213793: Phenotype Agouti
## Mettl3 - SRX2415040: genotype/variation Rosa26:BirA, METTL3 FLAG/Avi
## Rpa1 - SRX11676997: Pol III_degron mESCs IAA 1hr
## Rad21 - SRX1670207: Mll3-4 DKO Day1
## Med14 - SRX12622651: Treatment tamoxifen 96h
## Rnf2 - SRX6071873: Treatment tamoxifen 72h
## Nipbl - SRX15652647: NIPBL ChIP-M of wt stembryos at 168h
## Brd9 - SRX3751909: DMSO treatment
## Brd3 - SRX13439086: hIP-Seq_BRD3_BRD2_degron_3h_rep1
## Ncaph2 - SRX1688321: Transfection shGFP
## Dnttip1 - SRX13332921: SRM166780_Bruce4_rep2_DNTTIP1-AB1
## Ep400 - SRX9424476: Dox treatment
## Tet3 - SRX2545592: Tet3-CXXC
## Hif1a - SRX13945545: Hypoxia day 6
## Ercc3 - SRX11221923: MED4_degron_0h_XPB_ChIP_Rep2
## Klf5 - SRX334975: Cell line GCR8
## Fgfr1 - SRX868181: Retinoic acid treatment for 48h
## DNA-RNA hybrid - SRX13439146:  Rloop_BRD2_BRD3_BRD3_KO_rep1
## Max - SRX5850126: TKO TSA Max Rep3
## Biotin - SRX2388542: Gata3 overexpression
## Pcgf1 - SRX6071847: Treatment tamoxifen 72h
## Brca2 - SRX15464337: Zfp281 KO
## Gtf3c1 - SRX1688324: Transfection shGFP
## Gabpa - SRX5850112: Treatment 5 ng microl trichostatin A for 36h
## Sumo1 - SRX11525492: Cell line R1 Day 8
## Maz - SRX13440032: genotype MAZ HoxA delta 5-6
## Myc - SRX17793746: Cleavable Myc
## Ctnnb1 - SRX9101843: 3xFlag knock-in to endogeous Kdm2b C-terminus
## Zic2 - SRX5817869: 32h N2B27
## Zfp42 - SRX128354: Expressing V5-tagged REX1
## 0610010K14Rik - SRX19028417: ChIP_BAP18_dTAG-24h_mRBBP5-dTAG_ES_Rep2
## Ints5 - SRX19028478: ChIP_INTS5_dTAG-24h_mRBBP5-dTAG_ES_Rep2
## Runx1 - SRX8754229: R201Q no dox
##
## The list of replacement data are as follow:
##
## completeresultlistescsramerged[[1]][[1]] - mHG13hist:
## H4 - 78.776291 - SRX16815331
## H3K4me3 - 74.569790 - SRX5382140
## H3K4me2 - 55.066922 - SRX8818318
## NA (H3.3 - SRX2245482: SRX15799297, SRX2245484, SRX2245483)
## H3K64ac - 20.65010 - SRX4090625
## H3K122ac - 23.326960 - SRX4090627
## NA (H3K79me2 - SRX13332952: SRX747753)
## H3K36me3  - 34.416826 - SRX12382379 (SRX12382438, SRX12382488, SRX12382439, SRX12382441, SRX12382444, SRX12382440, SRX12382443, SRX12382445, SRX12382447) # nolint
## H2A.Z - 41.682600 - SRX7568534
## H3K9ac - 38.814532 - SRX13341956
## NA (H4ac - SRX10445814: SRX4384421, SRX4384469, SRX10445817, SRX10445828, SRX4384429, SRX4384445, SRX4384461, SRX10445815, SRX10445816, SRX10445830) # nolint
## NA (H4K5ac - SRX14210149)
## NA (H4K20ac - SRX13332963: SRX13332993)
## NA (H4K16ac - SRX298193)
## NA (H3K9K14ac - ERX245530: ERX245522, SRX396249, ERX245529, ERX245520, SRX396247) # nolint
## NA (H3K18ac - SRX10424671)
##
## completeresultlistescsramerged[[1]][[2]] - mHG13polTFs:
## Nanog - 85.468451 - SRX9486378 (Nanog - SRX9486380: SRX9486382)
## NA (Sp3 - SRX5371695: SRX5371703, SRX5371705, SRX5371702, SRX5371704, SRX5371694, SRX5371693, SRX5371696, SRX5371695) # nolint
## Kmt2d - 76.099426 - SRX2762151
## NA (Sp1 - SRX5371690: SRX5850140","SRX5850149", "SRX5850143","SRX5850146","SRX5850147","SRX5850141","SRX5371701","SRX5850148","SRX5371700","SRX5850145","SRX5850142","SRX5850144","SRX5850139","SRX5371689","SRX5371698","SRX5371688","SRX5371690")# nolint
## Tbp - 68.068834 - SRX9195301 (Tbp - SRX11677010: SRX6003871, SRX18520802, SRX6003867, SRX6003873, SRX1089849, SRX11677007, SRX11677008, SRX11677009) # nolint
## NA (Utf1 - SRX404068)
## Taf12 - 56.40535 - SRX11221932 (Taf12 - SRX11221956: SRX11221953, SRX11221948, SRX11221957, SRX11221940, SRX11221941, SRX11221954, SRX11221955) # nolint
## Setd1a       64.24474  SRX2762152 (Setd1a - SRX5600782: SRX5600783, SRX2762172) # nolint
## Sox2      68.068834  SRX9486379: (Sox2 - SRX9486383: SRX9486381) #nolint
## RNA polymerase II - 71.892925 - SRX8556273 (RNA polymerase II - SRX4178816: SRX19391750, SRX19391751, SRX6427836, SRX6427835, SRX19391749, SRX15799260, SRX15799262, SRX15799259, SRX19391748, SRX4178816) # nolint
## Yy1 - 47.227533 - SRX2873907
## Epitope tags - 72.275335 - SRX15888978 (Should be named Ascl1 - Epitope tags - SRX2486173) #nolint
## Epitope tags - 40.726577 -  SRX147771 (Should be named FBXL10 - SRX17045820, SRX9187061, SRX17045815, SRX9295987, SRX17045816, SRX2498439, SRX2498434, SRX9195279, SRX2498437, SRX470058, SRX15888976, SRX9295984, SRX9295985, SRX15888975, SRX15888979, SRX6813450) # nolint
## Epitope tags - 39.579350 - SRX15888985 (Should be named Ngn2 - SRX2498435, SRX9295986) # nolint
## Epitope tags - 36.520076 - SRX8832797 (Should be named Chd8 - SRX15888984, SRX470060, SRX19028394, SRX9195286, SRX9195292) # nolint
## Epitope tags - 33.078394 - SRX191070 (Should be name Kdm2b, see below if other dataset has more overlap - SRX9187041, SRX9195289, SRX14134853, SRX9195293, SRX2486172, SRX17045819) # nolint
## Epitope tags - 30.019120 - SRX9195280 (Should be named Taf1, see below if other dataset has more overlap - SRX8003646, SRX9195283, SRX1845479, SRX2646179, SRX9195290) # nolint
## Epitope tags  - 29.827916 - SRX5247634 (Should be named Zbtb11)
## Sin3a      56.787763 SRX4403317 (Sin3a - SRX9450269: SRX9450273, SRX9450270, SRX9450271, SRX9450268) # nolint
## Smad2      22.944551 SRX12702093 (Smad2  - SRX14663260: SRX5251386, SRX12702095, SRX12702096, SRX4578585, SRX14663257, SRX14663258, SRX14663259) # nolint
## NA (Taf3 - SRX13439109: SRX19028410, SRX19028412,SRX19028411, SRX19028413, SRX13439133, SRX13439132, SRX13439131) # nolint
## NA (Ssrp1 - SRX11677018: SRX11677016, SRX11677017, SRX11677015)
## NA (Esrrb - SRX11221882: SRX16038359 - SRX5023711, SRX6032065, SRX6032054, SRX6032063, SRX5023708, SRX6032053, SRX5023707, SRX6032051, SRX16038361, SRX6032064, SRX16038363, SRX6032052, SRX11221881, SRX16038364, SRX16038362, SRX16038360, SRX11221883, SRX11221884) # nolint
## GFP - 66.347992 - SRX373167 (Must be named Ash2L - GFP - SRX385377: SRX16844910, SRX16844908, SRX16844917) # nolint
## GFP - 63.097514 - SRX11417001 (Must be named RPB3 - SRX16844885, SRX16844916, SRX16844918, SRX16844909) # nolint
## GFP - 58.508604 - SRX11417015 (Must be named RPB10 - SRX16844890, SRX16844911, SRX16844924, SRX16844886, SRX16844919, SRX16844888, SRX16844893, SRX16844892, SRX16844891, SRX16844904, SRX16844896, SRX16844902, SRX16844889, SRX16844920, SRX16844921, SRX16844925) # nolint
## GFP - 53.346080 - SRX11416997 (Must be named RPB1 - SRX11417002, SRX390111, SRX16844905, SRX16844894, SRX16844897, SRX16844884, SRX16844903) # nolint
## GFP - 51.051625 - SRX11416999 (Must be named RPB2 - SRX16844887, SRX390112, SRX16844895, SRX16844923) # nolint
## GFP - 47.036329 - SRX11417014 (Must be named RPB9 - SRX11416998, SRX16844900)
## GFP - 35.181644 - SRX11417011 (Must be named RPB8 - SRX11417016, SRX11417013, SRX14663241, SRX16844912, SRX16844906, SRX175893, SRX16844901) # nolint
## GFP - 33.652008 - SRX11417005 (Must be named RPB5 - SRX16844898, SRX16844926)
## GFP - 26.003824 - SRX11417007 (Must be named RPB6 - SRX16844914, SRX14663243, SRX175894, SRX12702091, SRX16844907, SRX16844913, SRX16844928, SRX16844899, SRX9410191, SRX16844922, SRX14663240, SRX14663242, SRX14663239) # nolint
## NA (Taf1 - SRX11221945: SRX11221912, SRX11221909, SRX11221911, SRX11221910, SRX11221947, SRX19028407, SRX11221946, SRX19028408, SRX19028406, SRX19028409, SRX11221944) # nolint
## NA (Dendra2 - SRX4178805: SRX7568551, SRX7568556, SRX8003650, SRX7568545, SRX8003649, SRX8003652, SRX8003651, SRX7568552, SRX7568544, SRX4178813) # nolint
## NA (Otx2  - SRX9450552: SRX499141, SRX499117, SRX499143, SRX9450558, SRX9450557, SRX9450551) # nolint
## Med1 - 58.126195 - SRX9195310 (Med1 - SRX11221877: SRX11221878)
## Banp - 21.032505 - SRX8873444 (Banp - SRX10157266: SRX10978091, SRX10978090, SRX10978089, SRX10157269, SRX10157267, SRX10157268) # nolint
## NA (Ints3 - SRX11677013: SRX11677014, SRX11677011, SRX11677012)
## NA (Med4 - SRX11221901: SRX11221903, SRX11221904, SRX11221902)
## NA (Kdm6a - SRX13332957: SRX3141818, SRX13332951)
## Brd4 - 49.521989 - SRX4506780 (Brd4 - SRX3898501)
## NA (Kmt2c - SRX2752337: SRX2752334, SRX2752335, SRX2752333)
## NA (Mettl3 - SRX2415040: SRX2415029, SRX2415026)
## NA (Kdm4c - SRX213793: SRX424007)
## Med24 - 56.59656 - SRX5926394 (Med24 - SRX5926394)
## Nelfe - 41.873805 - SRX14134855 (Nelfe - SRX9450256: SRX17045818, SRX9450261, SRX9450257) # nolint
## NA (Tcf12 - SRX390394)
## Ctcf - 20.650096 - SRX8994331 (Ctcf - SRX10437500:  SRX8994333, SRX13228162, SRX13228165,  SRX5292078,  SRX2402721, SRX11237493,  SRX5125732, SRX12124750, SRX10437501, SRX11237523, SRX13439115, SRX10437504,  SRX9418410, SRX11237535, SRX11237517, SRX10437499, SRX10437506, SRX12204132,  SRX5125625, SRX11237511, SRX10560012, SRX10560014, SRX11237529, SRX11237505,  ERX2021911, SRX13439112, SRX10560013, SRX13439111, SRX10560015, SRX10437497,SRX10437505,SRX10437502, SRX10437496, SRX10437503, SRX13439114, SRX10437498) # nolint
## NA (Brd9 - SRX3751909:  SRX3751912, SRX4580516, SRX3751911)
## NA (Rnf2 - SRX482879: SRX186550, SRX6071871, SRX12124753, SRX381488, SRX7205310, SRX482877, SRX6071872, SRX7205309, SRX5242841, SRX1158299, SRX482878, SRX381487, SRX6071873) # nolint
## NA (Rpa1 - SRX11676997: SRX11676995, SRX11676994, SRX11676996)
## Nipbl - 23.900574 - SRX022697 (Nipbl - SRX15652647: SRX3753618, SRX15652646, SRX15652597) # nolint
## Rad21 - 26.959847 - SRX9107990 (Rad21 - SRX1670207: SRX6831871, SRX1670191, SRX7685626, SRX6831870, SRX1670203, SRX8642890) # nolint
## Brca2 - 27.91587 - SRX15464336 (Brca2 - SRX15464337: SRX15464338)
## NA (Med14 - SRX12622651: SRX12622647,SRX12622652, SRX12622648,SRX12622649,SRX12622654,SRX12622650, SRX12622658) # nolint
## NA (Ncaph2 - SRX1688321)
## NA (Hif1a - SRX13945545)
## NA (DNA-RNA hybrid - SRX13439146:  to remove anyway)
## NA (Brd3 - SRX13439086: SRX13439088, SRX13439089, SRX13439087)
## NA (Ep400 - SRX9424475: SRX9424473, SRX9424477, SRX9424476)
## NA (Max - SRX5850126: SRX5850121, SRX5850117, SRX5850119, SRX5850123, SRX5850125, SRX5850122) # nolint
## NA (Ercc3 - SRX11221923)
## NA (Tet3 - SRX2545592: SRX2545593, SRX2545591, SRX2545596, SRX2545595)
## NA (Dnttip1 - SRX13332921: SRX13332923, SRX13332925, SRX13332919)
## NA (Fgfr1 - SRX868181)
## NA (Gtf3c1 - SRX1688323: SRX1688324, SRX1688327)
## NA (Klf5 - SRX334975: SRX468970)
## NA (Gabpa - SRX5850112: SRX5850110, SRX5850105, SRX5850114)
## NA (Biotin - SRX2388542)
## NA (Pcgf1 - SRX6071847)
## NA (Maz - SRX13440032: SRX13440028)
## NA (Runx1 - SRX8754229)
## NA (Zfp42 - SRX128354)
## NA (Ctnnb1 - SRX9101843: SRX9101842)
## NA (Gata4 - SRX191946)
## NA (Sumo1 - SRX11525492: SRX11525495, SRX11525491, SRX11525496)
## NA (0610010K14Rik - SRX19028417: SRX19028415)
## NA (Zic2 - SRX5817869: SRX5817870, SRX5817882, SRX5817867)
## NA (Ints5 - SRX19028478)
## NA (Hoxb1 - SRX4980239)
## NA (Kdm2b - SRX186548)
## NA (E2f6 - SRX3581849)
##
## completeresultlistescsramerged[[2]][[1]] - mHG14hist:
## H3K4me3 - 73.556797 - SRX5382140
## H4 - 75.698324 - SRX16815331
## NA (H3K79me2 - SRX13332952: SRX747753)
## H3K4me2 - 60.893855 - SRX8818318
## NA (H3K64ac - SRX1560890)
## H3K122ac - 22.998138 - SRX4090627
## NA (H3.3 - SRX2245482: SRX15799297, SRX2245484, SRX2245483)
## H3K9ac - 34.078212 - SRX13341956
## H2A.Z - 36.312849 - SRX7568534
## H3K36me3 - 26.163873 - SRX12382379
## NA (H4ac - SRX10445814)
## NA (H4K5ac - SRX14210149)
## NA (H4K20ac - SRX13332963: SRX13332993)
## NA (H4K16ac - SRX298193)
## NA (H3K9K14ac - ERX245530)
## NA (H3K18ac - SRX10424671)
##
## completeresultlistescsramerged[[2]][[2]] - mHG14polTFs:
## Nanog - 84.823091 - SRX9486378
## Kmt2d - 76.536313 - SRX2762151
## NA (Sp3 - SRX5371695)
## Tbp - 66.945996 - SRX9195301
## RNA polymerase II - 75.884544 - SRX8556273
## Sox2 - 68.156425 - SRX9486379
## Setd1a - 59.86965 - SRX2762152
## NA (Sp1 - SRX5371688)
## NA (Utf1 - SRX404068)
## Taf12 - 57.35568 - SRX11221932
## Sin3a - 63.314711 - SRX4403317
## NA (Taf3 - SRX13439109)
## NA (Ssrp1 - SRX11677018)
## GFP - 65.921788 - SRX11417001 (Must be named RPB3)
## GFP - 64.990689 -  SRX373167 (Must be named Ash2L)
## GFP - 59.124767 - SRX11416997 (Must be named RPB1)
## GFP - 57.728119 - SRX11417015 (Must be named RPB10)
## GFP - 54.748603 - SRX11416999 (Must be named RPB2)
## GFP - 46.461825 - SRX11417014 (Must be named RPB9)
## GFP - 34.264432 - SRX11417005 (Must be named RPB5)
## GFP - 33.240223 - SRX11417011 (Must be named RPB8)
## GFP - 25.884544 - SRX11417007 (Must be named RPB6)
## NA (Esrrb - SRX11221882)
## NA (Smad2  - SRX14663260)
## NA (Taf1 - SRX11221945)
## Epitope tags - 68.528864 - SRX15888978 (Should be named Ascl1)
## Epitope tags - 44.040968 - SRX147771 (Should be named FBXL10)
## Epitope tags - 36.778399 - SRX15888985 (Should be named Ngn2)
## Epitope tags - 33.705773 - SRX8832797 (Should be named Chd8)
## Epitope tags - 29.329609 - SRX9195280 (Should be named Taf1)
## Epitope tags - 28.305400 - SRX5247634 (Should be named Zbtb11)
## Epitope tags - 23.74302 - SRX191070 (Should be name Kdm2b)
## NA (Dendra2 - SRX4178805)
## Banp - 17.783985 - SRX8873444
## NA (Med4 - SRX11221901)
## Med1 - 58.286778 - SRX9195310
## Yy1 - 47.299814 - SRX2873907
## Brd4 - 48.789572 - SRX4506780
## NA (Otx2  - SRX9450552)
## NA (Kdm6a - SRX13332957)
## Med24 - 58.47300 - SRX5926394
## NA (Kmt2c - SRX2752337)
## NA (Tcf12 - SRX390394)
## Nelfe - 40.875233 - SRX14134855
## Ctcf - 19.553073 - SRX8994331
## NA (Kdm4c - SRX213793)
## NA (Mettl3 - SRX2415040)
## NA (Rpa1 - SRX11676997)
## Rad21 - 29.981378 - SRX9107990
## NA (Med14 - SRX12622651)
## NA (Rnf2 - SRX6071873)
## Nipbl - 21.694600 -  SRX022697
## NA (Brd9 - SRX3751909)
## NA (Brd3 - SRX13439086)
## NA (Ncaph2 - SRX1688321)
## NA (Dnttip1 - SRX13332921)
## NA (Ep400 - SRX9424476)
## NA (Tet3 - SRX2545592)
## NA (Hif1a - SRX13945545)
## NA (Ercc3 - SRX11221923)
## NA (Klf5 - SRX334975)
## NA (Fgfr1 - SRX868181)
## NA (DNA-RNA hybrid - to remove anyway)
## NA (Max - SRX5850126)
## NA (Biotin - SRX2388542)
## NA (Pcgf1 - SRX6071847)
## Brca2 - 17.13222 - SRX15464336
## NA (Gtf3c1 - SRX1688323: SRX1688324, SRX1688327)
## NA (Gabpa - SRX5850112)
## NA (Sumo1 - SRX11525492)
## NA (Maz - SRX13440032)
## Myc - 21.042831 - SRX5276471 (Myc - SRX17793746: SRX17793736)
## NA (Ctnnb1 - SRX9101843)
## NA (Zic2 - SRX5817869)
## NA (Zfp42 - SRX128354)
## NA (0610010K14Rik - SRX19028417)
## NA (Ints5 - SRX19028478)
## NA (Runx1 - SRX8754229)
##
## Run on R 4.3.2
## Descostes June 2020 - modified dec 2023
#############################

library("pheatmap")


#############
## PARAMS
#############

resultpathvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG13_hist_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG13_pol_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG13_TF_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG14_hist_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG14_pol_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG14_TF_PST.txt") #nolint

resultnamevec <- c("mHG13hist", "mHG13pol", "mHG13TF", "mHG14hist", "mHG14pol",
    "mHG14TF")
repprefixvec <- c("mHG13", "mHG14")
suffixmerged <- c("hist", "polTFs")
experimentname <- c("peaksHG")
outputfolder <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/heatmaps") #nolint
percentthreshold <- c(20)
ignoreqval <- FALSE


#############
## FUNCTIONS
#############

.checkparams <- function(outputfolder, resultpathvec, resultnamevec,
        experimentname, repprefixvec) {

    if (!file.exists(file.path(outputfolder, experimentname)))
        dir.create(file.path(outputfolder, experimentname), recursive = TRUE)

    if (!isTRUE(all.equal(length(resultpathvec), length(resultnamevec))))
        stop("One name should be given for each result file.")

    if (!isTRUE(all.equal(length(experimentname), 1)))
        stop("Experiment name should be unique.")

    if (!isTRUE(all.equal(length(outputfolder), 1)))
        stop("Output folder shoulld be unique.")

    if (!isTRUE(all.equal(length(repprefixvec), 2)))
        stop("This script was designed to handle two replicates.")
}

.readandfilter <- function(resultpathvec, resultnamevec, chipatlascolnames,
    ignoreqval) {

        resultlist <- mapply(function(currentpath, currentname, colnamevec) {

            message("\t Processing ", currentname)
            fi <- read.delim(currentpath, stringsAsFactors = FALSE,
                header = FALSE)
            colnames(fi) <- colnamevec
            initialnb <- nrow(fi)
            nbquery <- as.numeric(unlist(lapply(strsplit(fi$OverlapQuery, "/"),
              "[", 1)))
            nbcontrol <- as.numeric(unlist(lapply(strsplit(fi$OverlapControl,
              "/"), "[", 1)))
            fi$QueryHigher <- nbquery > nbcontrol

            ## Filtering results having more overlap on random control
            if (isTRUE(all.equal(length(which(fi$QueryHigher)), 0)))
                stop("All random controls have more peaks than query.")
            fi <- fi[fi$QueryHigher, ]

            ## Filtering on the Qvalue
            if (!ignoreqval) {
                if (isTRUE(all.equal(length(which(fi$LogQval < -1.3)), 0)))
                    stop("All q-value are not significant. See line 73.")
                fi <- fi[fi$LogQval < -1.3, ]
            }

            ## Filtering the input expriments
            idxinput <- fi$AntigenClass != "Input control"
            if (isTRUE(all.equal(length(which(idxinput)), 0)))
                stop("All results are mock IP. See line 92.")
            fi <- fi[idxinput, ]

            percentkept <- round((nrow(fi) * 100) / initialnb)
            message("\t Keeping ", nrow(fi), "/", initialnb, "(", percentkept,
                "%)")
            if (!isTRUE(all.equal(percentkept, 0)))
                return(fi)
            else
                return(NA)
        }, resultpathvec, resultnamevec, MoreArgs = list(chipatlascolnames),
        SIMPLIFY = FALSE)

    return(resultlist)
}

.removeelements <- function(resultlist, idxremove) {

    if (!isTRUE(all.equal(length(idxremove), 0))) {
        message("\t Removing the following categories because of lack of ",
                "significant results: ", paste(resultnamevec[idxremove],
                        collapse = "-"))
        resultlist <- resultlist[-idxremove]
    }else {
        message("\t Nothing to remove.")
    }
    return(resultlist)
}

.filteroncells <- function(resultlist, lookesc = TRUE) {

    result <- mapply(function(currentdf, currentname, lookesc) {

                if (lookesc)
                    res <- currentdf[which(currentdf$Cell == "ES cells"), ]
                else
                    res <- currentdf[which(currentdf$Cell != "ES cells"), ]

                if (isTRUE(all.equal(nrow(res), 0))) {
                    if (lookesc)
                        stop("\t No lookesc found in ", currentname)
                    else
                        message("\t No other tissues found in ", currentname)
                    return(NA)
                }else {
                    message("\t Number of results for ", currentname, ": ",
                            nrow(res))
                    return(res)
                }
            }, resultlist, names(resultlist), MoreArgs = list(lookesc),
            SIMPLIFY = FALSE)
    return(result)
}

.calculatepercentoverlap <- function(dfname, df, pthres) {

    message("\t Processing ", dfname)
    ##Calculate and order by percentage of Overlap
    nbquery <- as.numeric(unlist(lapply(strsplit(df$OverlapQuery,"/"), "[", 1)))
    nbpeaks <- unique(as.numeric(unlist(lapply(strsplit(df$OverlapQuery, "/"),
        "[", 2))))
    if (!isTRUE(all.equal(length(nbpeaks), 1)))
        stop("Pb with the number of query peaks.")

    percentoverlap <- (nbquery * 100) / nbpeaks
    df <- cbind(df, percentoverlap)

    idxover <- which(df$percentoverlap >= pthres)
    if (isTRUE(all.equal(length(idxover), 0)))
        stop("No lines were kept, decrease the percentage value.\n")
    else
        message("\t\t Keeping ", length(idxover), "/", nrow(df))
    df <- df[idxover, ]
    df <- df[order(df$percentoverlap), ]
    return(df)
}

.buildantigenlist <- function(resultlist, thres = 1, keepmaxonly = FALSE) {

    result <- mapply(function(df, dfname, thres) {

        df <- .calculatepercentoverlap(dfname, df, thres)
        nbquery <- as.numeric(unlist(lapply(strsplit(df$OverlapQuery, "/"),
                "[", 1)))

        if (keepmaxonly) {

            ## Keeping the line with the highest percentage of overlap by
            ## antigen
            message("\t\t\t Keeping max overlap by antigens")
            antigenlist <- split(df, df$Antigen)
            antigenlist <- lapply(antigenlist, function(currentantigen) {
                if (isTRUE(all.equal(nrow(currentantigen), 1)))
                    return(currentantigen)
                idxmax <- which.max(currentantigen$percentoverlap)
                return(currentantigen[idxmax, ])})
            antigendf <- do.call(rbind, antigenlist)
        } else {
            message("\t\t\t Computing percentages only")
            antigendf <- df
        }

        ## Keeping only the name of the antigen and the percentage
        message("\t\t Returning ", nrow(antigendf), "/", length(nbquery))
        result <- data.frame(Antigen = antigendf$Antigen,
                            PercentOverlap = antigendf$percentoverlap,
                            SRA = antigendf$ID)
        return(result)
    }, resultlist, names(resultlist), MoreArgs = list(thres),
            SIMPLIFY = FALSE)
    return(result)
}

.mergepolandtf <- function(reslist, prefvec) {

    res <- lapply(prefvec, function(currentpref, reslist) {

        message("Processing ", currentpref)
        suffvec <- c("hist", "pol", "TF")
        reslistnames <- names(reslist)

        idxres <- sapply(suffvec,
            function(currentsuff, currentpref, reslistnames){
                idx <- grep(paste0(currentpref, currentsuff), reslistnames)
                if (isTRUE(all.equal(length(idx), 0)))
                    stop("The element was not found in the list when merging",
                        " pol and TF")
                return(idx)
            }, currentpref, reslistnames)
        resnames <- c(paste0(currentpref,"hist"), paste0(currentpref, "polTFs"))
        poltfdf <- rbind(reslist[[idxres[2]]], reslist[[idxres[3]]])
        res <- list(reslist[[idxres[1]]], poltfdf)
        names(res) <- resnames
        return(res)
    }, reslist)
    return(res)
}

.buildreplacementvectors <- function() {

    sra_to_filter <- list(
        "mHG13" = list(
            "mHG13hist" = c("SRX16815323", "SRX10304465", "SRX4681571",
                "SRX2245482", "SRX1560890", "SRX1560888", "SRX13332952",
                "SRX12382446", "SRX2955747", "SRX185848", "SRX10445814",
                "SRX14210149", "SRX13332963", "SRX298193", "ERX245530",
                "SRX10424671"),
            "mHG13polTFs" = c("SRX9486380", "SRX5371695", "SRX2762171",
            "SRX5371690", "SRX11677010", "SRX404068", "SRX11221956",
            "SRX5600782", "SRX9486383", "SRX4178816", "SRX5850150",
            "SRX2486173", "SRX2486173", "SRX2486173", "SRX2486173",
            "SRX2486173", "SRX2486173", "SRX2486173", "SRX9450269",
            "SRX14663260", "SRX13439109", "SRX11677018", "SRX11221882",
            "SRX385377", "SRX385377", "SRX385377", "SRX385377", "SRX385377",
            "SRX385377", "SRX385377", "SRX385377", "SRX385377", "SRX11221945",
            "SRX4178805", "SRX9450552", "SRX11221877", "SRX10157266",
            "SRX11677013", "SRX11221901", "SRX13332957", "SRX3898501",
            "SRX2752337", "SRX2415040", "SRX213793", "SRX5926394", "SRX9450256",
            "SRX390394", "SRX10437500", "SRX3751909", "SRX482879",
            "SRX11676997", "SRX15652647", "SRX1670207", "SRX15464337",
            "SRX12622651", "SRX1688321", "SRX13945545", "SRX13439146",
            "SRX13439086", "SRX9424475", "SRX5850126", "SRX11221923",
            "SRX2545592", "SRX13332921", "SRX868181", "SRX1688323", "SRX334975",
            "SRX5850112", "SRX2388542", "SRX6071847", "SRX13440032",
            "SRX8754229", "SRX128354", "SRX9101843", "SRX191946", "SRX11525492",
            "SRX19028417", "SRX5817869", "SRX19028478", "SRX4980239",
            "SRX186548", "SRX3581849")),
        "mHG14" = list(
            "mHG14hist" = c("SRX10304465", "SRX16815323", "SRX13332952",
                "SRX4681571", "SRX1560890", "SRX1560888", "SRX2245482",
                "SRX185848", "SRX2955747", "SRX12382446", "SRX10445814",
                "SRX14210149", "SRX13332963", "SRX298193", "ERX245530",
                "SRX10424671"),
            "mHG14polTFs" = c("SRX9486380", "SRX2762171", "SRX5371695",
                "SRX11677010", "SRX19391748", "SRX9486383", "SRX5600782",
                "SRX5371688", "SRX404068", "SRX11221956", "SRX9450269",
                "SRX13439109", "SRX11677018", "SRX385377", "SRX385377",
                "SRX385377", "SRX385377", "SRX385377", "SRX385377",
                "SRX385377", "SRX385377", "SRX385377", "SRX11221882",
                "SRX14663260", "SRX11221945", "SRX2486173", "SRX2486173",
                "SRX2486173", "SRX2486173", "SRX2486173", "SRX2486173",
                "SRX2486173", "SRX4178805", "SRX10157266", "SRX11221901",
                "SRX11221877", "SRX5850150", "SRX3898501", "SRX9450552",
                "SRX13332957", "SRX5926395", "SRX2752337", "SRX390394",
                "SRX9450256", "SRX10437500", "SRX213793", "SRX2415040",
                "SRX11676997", "SRX1670207", "SRX12622651", "SRX6071873",
                "SRX15652647", "SRX3751909", "SRX13439086", "SRX1688321",
                "SRX13332921", "SRX9424476", "SRX2545592", "SRX13945545",
                "SRX11221923", "SRX334975", "SRX868181", "SRX13439146",
                "SRX5850126", "SRX2388542", "SRX6071847", "SRX15464337",
                "SRX1688324", "SRX5850112", "SRX11525492", "SRX13440032",
                "SRX17793746", "SRX9101843", "SRX5817869", "SRX128354",
                "SRX19028417", "SRX19028478", "SRX8754229")))

    perc_replace <- list(
        "mHG13" = list(
            "mHG13hist" = c(78.776291, 74.569790, 55.066922, NA, 20.65010,
                23.326960, NA, 34.416826, 41.682600, 38.814532, NA, NA, NA, NA,
                NA, NA),
            "mHG13polTFs" = c(85.468451, NA, 76.099426, NA, 68.068834, NA,
                56.40535, 64.24474, 68.068834, 71.892925, 47.227533, 72.275335,
                40.726577, 39.579350, 36.520076, 33.078394, 30.019120,
                29.827916, 56.787763, 22.944551, NA, NA, NA, 66.347992,
                63.097514, 58.508604, 53.346080, 51.051625, 47.036329,
                35.181644, 33.652008, 26.003824, NA, NA, NA, 58.126195,
                21.032505, NA, NA, NA, 49.521989, NA, NA, NA, 56.59656,
                41.873805, NA, 20.650096, NA, NA, NA, 23.900574, 26.959847,
                27.91587, rep(NA, 28))),
        "mHG14" = list(
            "mHG14hist" = c(73.556797, 75.698324, NA, 60.893855, NA, 22.998138,
                NA, 34.078212, 36.312849, 26.163873, NA, NA, NA, NA, NA, NA),
            "mHG14polTFs" = c(84.823091, 76.536313, NA, 66.945996, 75.884544,
                68.156425, 59.86965, NA, NA, 57.35568, 63.314711, NA, NA,
                65.921788, 64.990689, 59.124767, 57.728119, 54.748603,
                46.461825, 34.264432, 33.240223, 25.884544, NA, NA, NA,
                68.528864, 44.040968, 36.778399, 33.705773, 29.329609,
                28.305400, 23.74302, NA, 17.783985, NA, 58.286778, 47.299814,
                48.789572, NA, NA, 58.47300, NA, NA, 40.875233, 19.553073, NA,
                NA, NA, 29.981378, NA, NA, 21.694600, rep(NA, 14), 17.13222,
                NA, NA, NA, NA, 21.042831, NA, NA, NA, NA, NA, NA)))

    sra_replace <- list(
        "mHG13" = list(
            "mHG13hist" = c("SRX16815331", "SRX5382140", "SRX8818318", NA,
                "SRX4090625", "SRX4090627", NA, "SRX12382379", "SRX7568534",
                "SRX13341956", NA, NA, NA, NA, NA, NA),
            "mHG13polTFs" = c("SRX9486378", NA, "SRX2762151", NA, "SRX9195301",
                NA, "SRX11221932", "SRX2762152", "SRX9486379", "SRX8556273",
                "SRX2873907", "SRX15888978", "SRX147771", "SRX15888985",
                "SRX8832797", "SRX191070", "SRX9195280", "SRX5247634",
                "SRX9450269", "SRX12702093", NA, NA, NA, "SRX373167",
                "SRX11417001", "SRX11417015", "SRX11416997", "SRX11416999",
                "SRX11417014", "SRX11417011", "SRX11417005", "SRX11417007", NA,
                NA, NA, "SRX9195310", "SRX8873444", NA, NA, NA, "SRX4506780",
                NA, NA, NA, "SRX5926394", "SRX14134855", NA, "SRX8994331", NA,
                NA, NA, "SRX022697", "SRX9107990", "SRX15464336", rep(NA, 28))),
        "mHG14" = list(
            "mHG14hist" = c("SRX5382140", "SRX16815331", NA, "SRX8818318", NA,
                "SRX4090627", NA, "SRX13341956", "SRX7568534", "SRX12382379",
                NA, NA, NA, NA, NA, NA),
            "mHG14polTFs" = c("SRX9486378", "SRX2762151", NA, "SRX9195301",
                "SRX8556273", "SRX9486379", "SRX2762152", NA, NA, "SRX11221932",
                "SRX4403317", NA, NA, "SRX11417001", "SRX373167", "SRX11416997",
                "SRX11417015", "SRX11416999", "SRX11417014", "SRX11417005",
                "SRX11417011", "SRX11417007", NA, NA, NA, "SRX15888978",
                "SRX147771", "SRX15888985", "SRX8832797", "SRX9195280",
                "SRX5247634", "SRX191070", NA, "SRX8873444", NA, "SRX9195310",
                "SRX2873907", "SRX4506780", NA, NA, "SRX5926394", NA, NA,
                "SRX14134855", "SRX8994331", NA, NA, NA, "SRX9107990", NA, NA,
                "SRX022697", rep(NA, 14), "SRX15464336", NA, NA, NA, NA,
                "SRX5276471", NA, NA, NA, NA, NA, NA)))

    antigen_replace <- list(
        "mHG13" = list(
            "mHG13hist" = c("H4", "H3K4me3", "H3K4me2", NA, "H3K64ac",
                "H3K122ac", NA, "H3K36me3", "H2A.Z", "H3K9ac", NA, NA, NA, NA,
                NA, NA),
            "mHG13polTFs" = c("Nanog", NA, "Kmt2d", NA, "Tbp", NA, "Taf12",
                "Setd1a", "Sox2", "RNA polymerase II", "Yy1", "Ascl1", "FBXL10",
                "Ngn2", "Chd8", "Kdm2b", "Taf1", "Zbtb11", "Sin3a", "Smad2", NA,
                NA, NA, "Ash2L", "RPB3", "RPB10", "RPB1", "RPB2", "RPB9",
                "RPB8", "RPB5", "RPB6", NA, NA, NA, "Med1", "Banp", NA, NA, NA,
                "Brd4", NA, NA, NA, "Med24", "Nelfe", NA, "Ctcf", NA, NA, NA,
                "Nipbl", "Rad21", "Brca2", rep(NA, 28))),
        "mHG14" = list(
            "mHG14hist" = c("H3K4me3", "H4", NA, "H3K4me2", NA, "H3K122ac", NA,
                "H3K9ac", "H2A.Z", "H3K36me3", NA, NA, NA, NA, NA, NA),
            "mHG14polTFs" = c("Nanog", "Kmt2d", NA, "Tbp", "RNApolII", "Sox2",
                "Setd1a", NA, NA, "Taf12", "Sin3a", NA, NA, "RPB3", "Ash2L",
                "RPB1", "RPB10", "RPB2", "RPB9", "RPB5", "RPB8", "RPB6", NA,
                NA, NA, "Ascl1", "FBXL10", "Ngn2", "Chd8", "Taf1", "Zbtb11",
                "Kdm2b", NA, "Banp", NA, "Med1", "Yy1", "Brd4", NA, NA, "Med24",
                NA, NA, "Nelfe", "Ctcf", NA, NA, NA, "Rad21", NA, NA,
                "Nipbl", rep(NA, 14), "Brca2", NA, NA, NA, NA, "Myc", NA, NA,
                NA, NA, NA, NA)))

    return(list(sra_to_filter, perc_replace, sra_replace, antigen_replace))
}

.testrepnames <- function(namesrepsra, namesreppercreplace, namesrepsrareplace,
    namesantigenreplace) {

        if (!isTRUE(all.equal(namesrepsra, namesreppercreplace)) ||
            !isTRUE(all.equal(namesrepsra, namesrepsrareplace)) ||
            !isTRUE(all.equal(namesrepsra, namesantigenreplace)))
                stop("Problem in names of replicates")
}

.testna <- function(percreplace, srareplace, antireplace) {

    idxnaperc <- which(is.na(percreplace))
    idxnasra <- which(is.na(srareplace))
    idxnaanti <- which(is.na(antireplace))

    if (!isTRUE(all.equal(idxnaperc, idxnasra)) ||
        !isTRUE(all.equal(idxnaperc, idxnaanti)))
            stop("There is an error in NA values. They are not at",
                " the same indexes between vectors")
}

.testelementvec <- function(namessra, namespercreplace, namessrareplace,
    namesantigenreplace, sra, percreplace, srareplace, antireplace) {

                    if (!isTRUE(all.equal(namessra, namespercreplace)) ||
                    !isTRUE(all.equal(namessra, namessrareplace)) ||
                    !isTRUE(all.equal(length(namessra),
                        length(namesantigenreplace))))
                        stop("Problem in names of elements")

                if (!isTRUE(all.equal(length(sra), length(percreplace))) ||
                    !isTRUE(all.equal(length(sra), length(srareplace))) ||
                    !isTRUE(all.equal(length(sra), length(antireplace)))) {
                        message("sra: ", length(sra))
                        message("percreplace: ", length(percreplace))
                        message("srareplace: ", length(srareplace))
                        message("antireplace: ", length(antireplace))
                        stop("Vectors differs in length")
                }

                .testna(percreplace, srareplace, antireplace)
}

.testvaluesvec <- function(sra_to_filter, perc_replace, sra_replace,
    antigen_replace) {

        # repsra <- sra_to_filter[[1]]
        # reppercreplace <- perc_replace[[1]]
        # repsrareplace <- sra_replace[[1]]
        # repantigenreplace <- antigen_replace[[1]]
        # namesrepsra <- names(sra_to_filter)[1]
        # namesreppercreplace <- names(perc_replace)[1]
        # namesrepsrareplace <- names(sra_replace)[1]
        # namesantigenreplace <- names(antigen_replace)[1]
    invisible(mapply(function(repsra, reppercreplace, repsrareplace,
        repantigenreplace, namesrepsra, namesreppercreplace, namesrepsrareplace,
        namesantigenreplace) {

            .testrepnames(namesrepsra, namesreppercreplace, namesrepsrareplace,
                namesantigenreplace)

            # sra <- repsra[[1]]
            # percreplace <- reppercreplace[[1]]
            # srareplace <- repsrareplace[[1]]
            # antireplace <- repantigenreplace[[1]]
            # namessra <- names(repsra)[1]
            # namespercreplace <- names(reppercreplace)[1]
            # namessrareplace <- names(repsrareplace)[1]
            # namesantireplace <- names(repantigenreplace)[1]
            mapply(function(sra, percreplace, srareplace, antireplace, namessra,
            namespercreplace, namessrareplace, namesantireplace) {

                message("Testing ", namespercreplace)

                .testelementvec(namessra, namespercreplace, namessrareplace,
                    namesantigenreplace, sra, percreplace, srareplace,
                    antireplace)

            }, repsra, reppercreplace, repsrareplace, repantigenreplace,
                names(repsra), names(reppercreplace), names(repsrareplace),
                names(repantigenreplace))
    }, sra_to_filter, perc_replace, sra_replace, antigen_replace,
        names(sra_to_filter), names(perc_replace), names(sra_replace),
        names(antigen_replace)))
}

.replaceandremovena <- function(repsrafilter, reppercreplace, repsrareplace,
    repantireplace, resrep) {

                    reslist <- mapply(function(currentname, repsrafilter,
                        reppercreplace, repsrareplace, repantireplace, resrep) {

                            message("\t\t Elements: ", currentname)
                            elsratofilter <- repsrafilter[[currentname]]
                            newperc <- reppercreplace[[currentname]]
                            newsra <- repsrareplace[[currentname]]
                            newanti <- repantireplace[[currentname]]
                            restomodif <- resrep[[currentname]]

                            idx <- match(elsratofilter, restomodif$SRA)
                            if (any(is.na(idx))) {
                                idxna <- which(is.na(idx))
                                stop("Problem in retrieving sra to filter: ",
                                    elsratofilter[idxna])
                            }
                            restomodif$Antigen[idx] <- newanti
                            restomodif$PercentOverlap[idx] <- newperc
                            restomodif$SRA[idx] <- newsra

                            ## Removing lines containing NA
                            idxremove <- which(is.na(restomodif$PercentOverlap))
                            if (!isTRUE(all.equal(length(idxremove), 0)))
                                restomodif <- restomodif[-idxremove, ]
                            return(restomodif)
                    }, names(resrep), MoreArgs = list(repsrafilter,
                        reppercreplace, repsrareplace, repantireplace, resrep),
                        SIMPLIFY = FALSE)
                    return(reslist)
}

.replaceelementschipatlas <- function(sra_to_filter, perc_replace, sra_replace,
    antigen_replace, resultlistescsramerged) {

        res <- mapply(function(currentrepname, sratofilter, percreplace,
            srareplace, antigenreplace, res) {

        ## Retrieving the lists of histones and polTFs for current rep
        message("\t Processing ", currentrepname)
        repsrafilter <- sratofilter[[currentrepname]]
        reppercreplace <- percreplace[[currentrepname]]
        repsrareplace <- srareplace[[currentrepname]]
        repantireplace <- antigenreplace[[currentrepname]]
        resrep <- res[[currentrepname]]

        reslist <- .replaceandremovena(repsrafilter, reppercreplace,
                repsrareplace, repantireplace, resrep)
        return(reslist)
    }, names(sra_to_filter),
        MoreArgs = list(sra_to_filter, perc_replace, sra_replace, 
            antigen_replace, resultlistescsramerged), SIMPLIFY = FALSE)
    return(res)
}



.retrieveall <- function(res, sectionvec) {

    message("\t Retrieving all ")
    resdflist <- mapply(function(currentsection, currentrep, repname) {

        message("\t\t Extracting ", currentsection, " from ", repname)
        sravec <- currentrep[[currentsection]]$SRA
        antigenvec <- currentrep[[currentsection]]$Antigen
        return(data.frame(SRA = sravec,
            Antigen = antigenvec))
    }, sectionvec, res, names(res), SIMPLIFY = FALSE)
    resdf <- do.call("rbind", resdflist)

    ## Only missing antigen will be used for completing data, following line is
    ## ok. See function .completemissingdata
    ## The sra of rep1 will be used for missing datas in rep2 and the other
    ## way around
    dupantigen <- which(duplicated(resdf$Antigen))
    if (!isTRUE(all.equal(length(dupantigen), 0)))
        resdf <- resdf[-dupantigen, ]

    return(resdf)
}

.retrievepercentvec <- function(missingsra, repcomplete) {

    return(unlist(lapply(missingsra, function(currentsra, repcomplete) {
                    idxsracomplete <- which(repcomplete$SRA == currentsra)
                    if (!isTRUE(all.equal(length(idxsracomplete), 0))) {
                        message("\t\t\t Complementary sra found")
                        return(repcomplete$PercentOverlap[idxsracomplete])
                    }else {
                        return(0)
                    }
                }, repcomplete)))
}

.completemissingdata <- function(sectionvec, res, rescomplete, allsra,
    allantigen) {

    rescompletedlist <- mapply(function(currentsection, currentrep, repname,
        currentrepcomplete, allsra, allantigen) {

            rep <- currentrep[[currentsection]]
            repcomplete <- currentrepcomplete[[currentsection]]
            idx <-  match(allantigen, rep$Antigen)
            idxna <- which(is.na(idx))

            if (!isTRUE(all.equal(length(idxna), 0))) {
                message("\t Completing missing data for ", currentsection,
                    " from ", repname)
                missingantigen <- allantigen[idxna]
                missingsra <- allsra[idxna]
                message("\t\t The missing data are ",
                    paste(missingantigen, collapse = "-"), " having sra ",
                    paste(missingsra, collapse = "-"))
                percentvec <- .retrievepercentvec(missingsra, repcomplete)
                missingdf <- data.frame(Antigen = missingantigen,
                                        PercentOverlap = percentvec,
                                        SRA = missingsra)
                rep <- rbind(rep, missingdf)
            }
            return(rep)
    }, sectionvec, res, names(res), rescomplete,
        MoreArgs = list(allsra, allantigen), SIMPLIFY = FALSE)
    return(rescompletedlist)
}

.plotheatmap <- function(mat, scalemethod, ignoreqval, experimentname,
    percentthreshold, outputfolder, sectiontype) {

        main_title <- if (!ignoreqval) experimentname else paste(experimentname,
                    "No Qval", sep = "-")
        outfile <- file.path(outputfolder, experimentname, sectiontype,
                paste0(experimentname, "-", percentthreshold,
                        if (ignoreqval) "-noQval", "-", scalemethod, ".pdf"))

        message("\t Plotting to ", outfile)

        pheatmap::pheatmap(mat,
        scale = scalemethod,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        clustering_method = "complete",
        annotation_names_row = TRUE,
        annotation_names_col = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        main = main_title,
        fontsize_row = 4,
        filename = outfile)
}

.callheatmapgeneration <- function(res, rescomplete, sectionvec, sectiontype,
    outputfolder, experimentname, percentthreshold, ignoreqval) {

    message("Generating heatmap for ", sectiontype)
    alldf <- .retrieveall(res, sectionvec)
    allsra <- alldf$SRA
    allantigen <- alldf$Antigen

    ## Completing missing data if necessary
    rescompletedlist <- .completemissingdata(sectionvec, res, rescomplete,
        allsra, allantigen)

    ## Merge results by antigen (assuming two replicates, verified when checking
    ## parameters)
    message("Merging replicates and create result table")
    resmerged <- merge(rescompletedlist[[1]], rescompletedlist[[2]],
        by = "Antigen")

    ## Writing the result table
    message("Writing output table result")
    outfold <- file.path(outputfolder, experimentname, sectiontype)
    if(!file.exists(outfold))
        dir.create(outfold, recursive = TRUE)
    filename <- paste0(experimentname, "-", percentthreshold,
        if (ignoreqval) "-noQval", ".txt")
    outfile <- file.path(outfold, filename)
    write.table(resmerged, file = outfile, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = TRUE)

    ## Restricting table for plotting
    mat <- cbind(resmerged$PercentOverlap.x, resmerged$PercentOverlap.y)
    rownames(mat) <- resmerged$Antigen
    colnames(mat) <- sectionvec

    ## Plotting heatmap
    .plotheatmap(mat, "none", ignoreqval, experimentname, percentthreshold,
        outputfolder, sectiontype)
    .plotheatmap(mat, "row", ignoreqval, experimentname, percentthreshold,
        outputfolder, sectiontype)
}



#############
## MAIN
#############

####
## PART 1: Preparing the data
####

## Checking parameters
.checkparams(outputfolder, resultpathvec, resultnamevec, experimentname,
    repprefixvec)

if (ignoreqval)
    warning("!!!!!!!! #############\n Q-val is ignored. Only matches are ",
            "reported without guarantee that they are significant.\n",
            "############# !!!!!!!!")

## Reading and filtering results files
message("Reading ChIP-Atlas results")
chipatlascolnames <- c("ID", "AntigenClass", "Antigen", "CellClass", "Cell",
        "NumOfPeaks", "OverlapQuery", "OverlapControl", "LogPVal", "LogQval",
        "FoldEnrichment")
resultlist <- .readandfilter(resultpathvec, resultnamevec, chipatlascolnames,
    ignoreqval)
names(resultlist) <- resultnamevec
idxremove <- which(is.na(resultlist))
resultlist <- .removeelements(resultlist, idxremove)
nbpeaksvec <- unlist(lapply(resultlist, function(currentcat){
                    nbpeaks <- unique(as.numeric(unlist(lapply(
                        strsplit(currentcat$OverlapQuery, "/"), "[", 2))))
                    if (!isTRUE(all.equal(length(nbpeaks), 1)))
                        stop("nbpeaks should be unique")
                    return(nbpeaks)
                }))


## Retrieving embryonic stem cells and other cells
message("\n\n Retrieving embryonic stem cells")
resultlistesc <- .filteroncells(resultlist)
idxremove <- which(is.na(resultlistesc))
resultlistesc <- .removeelements(resultlistesc, idxremove)
if (!isTRUE(all.equal(length(idxremove), 0)))
    nbpeaksvec <- nbpeaksvec[-idxremove]

## Building antigen unique list
message("Building antigen unique lists for Embryonic Stem Cells:")
resultlistescsra <- .buildantigenlist(resultlistesc, thres = percentthreshold,
    keepmaxonly = TRUE)
completeresultlistescsra <- .buildantigenlist(resultlistesc)

if (!isTRUE(all.equal(names(resultlistescsra),
    names(completeresultlistescsra))))
    stop("The two lists do not contain the same elements.")

## Merging pol and TFs for each replicates
resultlistescsramerged <- .mergepolandtf(resultlistescsra, repprefixvec)
completeresultlistescsramerged <- .mergepolandtf(completeresultlistescsra,
    repprefixvec)
names(resultlistescsramerged) <- repprefixvec
names(completeresultlistescsramerged) <- repprefixvec


####
## PART 2: Filtering the data
####

## The list of sra to filter was determined manually looking at the sra records.
##
## See the top of the script and the function .buildreplacementvectors() for
## details
resultsreplacement <- .buildreplacementvectors()
sra_to_filter <- resultsreplacement[[1]]
perc_replace <- resultsreplacement[[2]]
sra_replace <- resultsreplacement[[3]]
antigen_replace <- resultsreplacement[[4]]

## Verifying elements of each vector that must be of same length
.testvaluesvec(sra_to_filter, perc_replace, sra_replace, antigen_replace)

## Replacing values for each vector
message("Replacing values")
resultlistescsramerged <- .replaceelementschipatlas(sra_to_filter, perc_replace,
    sra_replace, antigen_replace, resultlistescsramerged)


####
## PART 3: Generating the heatmap
####

## Completing each element of the list with missing sra
## First retrieving all sra of the list
histsections <- paste0(repprefixvec, suffixmerged[1])
poltfsections <- paste0(repprefixvec, suffixmerged[2])
.callheatmapgeneration(resultlistescsramerged, completeresultlistescsramerged,
    histsections, "Histones", outputfolder, experimentname, percentthreshold,
    ignoreqval)
.callheatmapgeneration(resultlistescsramerged, completeresultlistescsramerged,
    poltfsections, "PolII-TFs", outputfolder, experimentname,
    percentthreshold, ignoreqval)
