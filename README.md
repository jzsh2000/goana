goana
=====

GO analysis, use `goana` function in
[limma](http://bioconductor.org/packages/release/bioc/html/limma.html) package
under the hook.

usage
-----

> Note: use `goana.mouse.R` for mouse gene list.

```bash
Rscript goana.R <geneset.txt> [output.csv]
```

Here `geneset.txt` is a gene list file with one gene name per line, see
[geneset.txt](shiny/geneset.txt) as an example. If the parameter `output.csv`
isn't set, the output csv file will be saved in the same directory as
`geneset.txt` file.

output
------

* __GO\_id__ - Gene ontology term ID (e.g. GO:0006915)
* __Term__ - Gene ontology term (e.g. apoptotic process)
* __Ont__ - Gene ontology category (here always BP)
* __N__ - Number of genes in this term (e.g. 1804)
* __DE__ - Number of user-provided genes in the term (e.g. 109)
* __P.DE__ - Enrichment P-value (e.g. 1.83e-73)
* __FDR__ - Adjusted P-value (e.g. 1.18e-69)
* __gene\_list__ - Names of user-provided genes in this term (e.g. "AIFM3,
  ANXA1, APP, ATF3, BAX, BCAP31, BCL10, BCL2L1, BCL2L10, BCL2L11, BCL2L2, BID,
  BIK, BIRC3, BMF, BMP2, BNIP3L, BRCA1, BTG2, CASP1, CASP2, CASP3, CASP4,
  CASP6, CASP7, CASP8, CASP9, CAV1, CCND2, CD14, CD2, CD38, CD44, CDKN1A,
  CDKN1B, CFLAR, CLU, CREBBP, CTH, CTNNB1, CYLD, DAP, DAP3, DDIT3, DFFA,
  DIABLO, DNAJA1, DNAJC3, DNM1L, EGR3, ERBB3, F2R, FAS, FASLG, GADD45A,
  GADD45B, GPX1, GSN, H1F0, HGF, HMGB2, HMOX1, HSPB1, IER3, IFNB1, IGF2R, IL1A,
  IL1B, IL6, IRF1, JUN, KRT18, LEF1, LGALS3, LMNA, MADD, MCL1, MGMT, PAK1,
  PDCD4, PDGFRB, PEA15, PMAIP1, PPP3R1, PPT1, PRF1, PSEN1, PSEN2, PTK2, RARA,
  RELA, RHOB, RHOT2, ROCK1, SATB1, SOD1, SOD2, SQSTM1, TGFB2, TIMP1, TIMP3,
  TNF, TNFRSF12A, TNFSF10, TOP2A, TSPO, TXNIP, VDAC2, XIAP")


shiny application
-----------------

Interactive web application is available at the [shiny](shiny) directory, and
you should have [Rstudio](https://www.rstudio.com/) to start this application.
