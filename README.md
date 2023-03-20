# SpaLine
this repository is for spatial transcriptome analysis

## Supported spatial transcriptome technologies
- 10x visium 
- stereoseq 

## Dependencies 
ggplot2
Seurat
SpaTalk
CellphoneDB
celldex


## Usage
### preprocess 
##### read file and generate quality control figure
```R
x<-Standard_Seu(x,project,QC_dir=NULL,spatial=NULL) 
### x could be a seurat object, a count matrix or a parent directory path to cellranger filtered_feature_bc_matrix
### project should be your sample name for identification
### spatial_dir should be the cellranger out spatial directory. Due to lack of spatial HE images, if your data is from stereoseq, this parameter should set to NULL.
### QC_dir, if not null, should be the QC outdir path you want save your QC figures. Make sure you have access to the QC directory 
```

##### marker nuclei  (combined infromation in HE images, only used in cellranger out data with images)
```R
x<-Mark_Nuclei(path,project,QC_dir=NULL)
### path, a parent directory path to cellranger filtered_feature_bc_matrix
```
the spot with low nuclei number(less than 1) can be removed by subset
```
real_spots<-colnames(x@meta.data)[which(x@meta.data$n_nuclei > 1)]
x<-subset(x,cells=real_spots)
```

### cell annotation
cell annotaion can be done by single cell RNAseq annotation method like singleR using a public label reference or spatial deconvolution method with a sc reference
```R
#### using singleR with public database
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
ref=hpca.se
label=hpca.se$label.main
x<-Cell_Annotation(x,method="singleR",ref=ref,label=label,cluster=TRUE)
#### deconvolution by spatalk
x<-Standard_Spatalk(sp=x,sc,sc_celltype,sp_meta=NULL,geneinfo,lrpairs,pathways)
```
### cell communication 
the output file of Standard_Spatalk already containing the LR and target path information, but we can also do cell communication mannually using CellphoneDB
```R
### first, select your interest area in seurat object and save the neccesary file for CellphoneDB
select_x<-SpArea_Select(x,ident="Tumor_Epi_Boundary",group="predicted.id",interestsID=c("Tumor","Epi"),spatial10x=FALSE,saveCPDB_dir=NULL)
### second, run CellphoneDB in shell
cellphonedb method statistical_analysis yourmetafile.txt yourcountsfile.txt --iterations=10 --threads=2
### extract Ligand-receptor interaction 
mat<-Cell_Communication(pathlist,rm_complex=FALSE,top=20)
########pathlist could be a list of path to the out of cellphonedb of different selected area 
```
### visualization
#### Celltype distribution --bar
```R
CellRatio_bar(x)  ##x should be a seurat object or a SpaTalk object
```
#### Celltype distribution -- spatial
```R

``` 
#### Function analysis 
```
GO_KEGG<-Fun_Plot(x,species,showCategory=20,celltype_sender=NULL,celltype_receiver=NULL,top=20)
###### celltype_sender and celltype_receiver and top only needed when x is a SpaTalk object
```

#### Ligand-Receptor interaction 
```R
LR_plots(mat)  #mat is the output from Cell_Communication
```
#### 

