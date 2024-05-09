# scRNAseq_project
Self-learning project for scRNAseq

For the porject, I used the NSCLC project data (https://www.10xgenomics.com/datasets/20-k-mixture-of-nsclc-dt-cs-from-7-donors-3-v-3-1-3-1-standard-6-1-0).
I followed the tutorial of Dr. Khusbu Patel, available via the link
https://youtu.be/5HBzgsz8qyk?si=ifN56YeU0pylqrfG

Following several online tutorials, I completed the analysis in R using respective packages.

As a learning exercise, I used different dimensions (10 and 15) for clustering, and tried to 
interpret the UMAP plots.

Some key indings were:
1. UMAP images appear different for different dimensions. One possibility is that
   cahnging dimensions can change the patterns that can be finally observed as
   data may be lost in lower dimensions.
2. Using PCA, I identified similar number of clusters although number of dimensions
   were different. However, as the resolution increased, the number of clusters increased.
   A fine balance has to be obtained between the number of possible clusters and resolution.

Acknowledgements: I thank Dr. Patel for providing these tutorials via opensource platforms.
