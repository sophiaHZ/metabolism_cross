# metabolism_cross pipeline

##usage
Rscript /home/swf/bin/metabolism/new_pipeline.R -d /home/swf/Project/S9540CA/new -q /home/swf/Project/S9540CA/MS_identified_information.xlsx -m BvsA_meta.txt -s /home/swf/Project/S9540CA/sample.txt -g BvsA -r 1.5 -p 0.05 -l 1 -f 0.05 -v 1 -k /home/swf/Project/S9540CA/new/annotation/KEGG_pathway_annotation.xls -c /home/swf/Project/S9540CA/new/annotation/KEGG_compound_annotation.xls -o hsa

-d: working dir (results dir),default MS file dir
-q: MS information xlsx
-m: metabolism file including sample quantification information and Fold.Change, FDR or log2FC columns
-s: sample.txt including two columns named Sample and Type
-g: analysis compare group, such as BvsA
-r: Diff protein ratio, default 1.5
-p: Diff protein p value, default 0.05
-l: Diff metabolite log2fc, default 1
-f: Diff metabolite FDR, default 0.05
-v: Diff metabolite VIP value, default 1
-k: Protein KEGG_pathway_annotation.xls
-c: Metabolite Compound_pathway_annotation.xls
-o: Organism abbreviation, default hsa

KEGG_pathway_annotation.xls dir has to have two files: all_pathway.txt and all_list.txt
