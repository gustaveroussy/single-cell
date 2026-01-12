<p align="center">
  <img src="https://github.com/gustaveroussy/single-cell/blob/master/images/wiki_header.png" title="wiki_header">
</p>


# Single-cell RNA-seq analysis pipeline

---
### Pipeline goal:  
Perform single-cell RNA-seq analysis from FastQ files to cerebro file for 10XGenomics technology data.

---

### Available Steps:
<b>Individual Analysis</b>
* Alignments (Alignment_countTable_GE, Alignment_countTable_ADT, Alignment_annotations_TCR_BCR),
* QC of Droplets and filetring (Droplets_QC_GE, Filtering_GE),
* Normalization and dimension Reduction (Norm_DimRed_Eval_GE),
* Clustering, Identification of Marker Genes and Annotation of clusters (Clust_Markers_Annot_GE),
* Integration of several additional "omics" (Adding_ADT, Adding_TCR, Adding_BCR),
* Creation of a Cerebro object to help vizualisation of results (Cerebro).

<p align="center">
<img src="https://github.com/gustaveroussy/single-cell/blob/master/images/individual_analysis_pipeline.png" width="600" title="individual_analysis_pipeline">
</p>

<b>Integrated Analysis of several samples</b>
* Integration, Normalization and dimension Reduction (Int_Norm_DimRed_Eval_GE)
* Clustering, Identification of Marker Genes and Annotation of clusters (Int_Clust_Markers_Annot_GE),
* Integration of several additional "omics" (Int_Adding_ADT, Int_Adding_TCR, Int_Adding_BCR),
* Creation of a Cerebro object to help vizualisation of results (Cerebro).

<b>Grouped Analysis (no integration) of several samples</b>
* Merger, Normalization and dimension Reduction (Grp_Norm_DimRed_Eval_GE)
* Clustering, Identification of Marker Genes and Annotation of clusters (Grp_Clust_Markers_Annot_GE),
* Integration of several additional "omics" (Grp_Adding_ADT, Grp_Adding_TCR, Grp_Adding_BCR),
* Creation of a Cerebro object to help vizualisation of results (Cerebro). 

### Future Developments:
* Improve cell annotation,
* Update to CellRanger 6.0.0,
* Analysis of scATAC-seq,
* ...

<br>

See complete documentation on the [wiki](https://github.com/gustaveroussy/single-cell/wiki)

<p align="center">
  <img src="https://github.com/gustaveroussy/single-cell/blob/master/images/wiki_footer.jpg" title="wiki_footer">
</p>
