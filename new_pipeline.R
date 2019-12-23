library(openxlsx)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(RColorBrewer) 
library(pheatmap)
library(data.table)
library(getopt)

commands=matrix(c('help','h',0,'logical','help document',
		'dir','d',2,'character','Analysis path',
		'Msfile','q',1,'character','MS_identified_information.xlsx',
		'Metafile','m',1,'character','Metabolism_information.txt',
		'samplefile','s',1,'character','sample.txt with column Sample and Type',
		'group','g',1,'character','compare group',
		'protein_fold','r',2,'double','protien diff ratio',
		'protein_pvalue','p',2,'double','protein diff pvalue',
		'meta_log2fc','l',2,'double','metabolite diff ratio',
		'meta_fdr','f',2,'double','metabolite diff FDR',
		'meta_vip','v',2,'double','metabolite VIP threshold',
		'protein_annotation','k',1,'character','Protein Annotation file',
		'meta_annotation','c',1,'character','Metabolite Annotation file',
		'organism','o',2,'character','KEGG organism abbreviation'),byrow=T,ncol=5)
#commands
args=getopt(commands)
annotation_dir = dirname(args$meta_annotation)
#args
if(!is.null(args$help)){
	cat(getopt(commands, usage=TRUE))
	q(status=1)
}
if(is.null(args$dir)){args$dir = dirname(args$Msfile)}
if(is.null(args$protein_fold)){args$protein_fold=1.5}
if(is.null(args$protein_pvalue)){args$protein_pvalue=0.05}
if(is.null(args$meta_log2fc)){args$meta_log2fc=1}
if(is.null(args$meta_fdr)){args$meta_fdr=0.05}
if(is.null(args$meta_vip)){args$meta_vip=1}
if(is.null(args$organism)){args$organism="hsa"}
if(is.null(args$Msfile) || is.null(args$Metafile) || is.null(args$samplefile) || is.null(args$group) || is.null(args$protein_annotation) ||is.null(args$meta_annotation)){
	cat(getopt(commands, usage=TRUE))
	q(status=1)
}

###read data
ms = read.xlsx(args$Msfile,sheet=2,rowNames=T,colNames=T)
meta = read.table(args$Metafile,header=T,row.names=1,check.names=F,sep="\t",quote="",dec=".",stringsAsFactors=F)
samples = read.table(args$samplefile,header=T,row.names=1,check.names=F,sep="\t",quote="",dec=".")
colnames(ms) = sub("/","vs",colnames(ms),)

##plot figure 2
gpath=paste(args$dir,args$group,sep="/")
if(!dir.exists(gpath)){
	dir.create(gpath)	
}
table_path = paste0(gpath,"/table")
if(!dir.exists(table_path)){
	dir.create(table_path)
}
figure_path = paste0(gpath,"/Figure")
if(!dir.exists(figure_path)){
	dir.create(figure_path)
}
kegg_path = paste0(gpath,"/pathway_image")
if(!dir.exists(kegg_path)){
        dir.create(kegg_path)
}

print("Finish preparations at")
t0 = Sys.time()
t0

setwd(gpath)
sample = rownames(samples)[samples$Type %in% unlist(strsplit(args$group,"vs"))]
d1 = ms[,c(sample,paste0(args$group,".Ratio"),paste0(args$group,".P.value"))]
colnames(d1)=c(sample,"Ratio","P.value")
d1$log2FC=log2(d1$Ratio)
d1$log10p=-log10(d1$P.value)
d1$Regulation = 'unchange'
d1$Regulation[d1$Ratio>args$protein_fold&d1$P.value<args$protein_pvalue]='up'
d1$Regulation[d1$Ratio<1/args$protein_fold&d1$P.value<args$protein_pvalue]='down'
ymax=max(d1$log10p,na.rm=T)+1
xmax=max(abs(d1$log2FC),na.rm=T)+1
if(ymax > 20){
        ymax <- 20
}
if(xmax >6){
        xmax=6
}
p1=ggplot(data=d1,aes(x=log2FC,y=log10p,colour=Regulation))+geom_point(size=2)+
  scale_colour_manual(values=c('unchange'="gray",'down'="#4575B4",'up'="#D73027"))+
  labs(x='Log2FC proteins',y="-log10 P value")+theme_bw()+ylim(0,ymax)+xlim(-xmax,xmax)

if ("log2FC" %in% colnames(meta)){
	d2 = meta[,c(sample,"VIP","FDR","Fold.Change","log2FC")]
	colnames(d2)=c(sample,"VIP","P.value","Ratio","log2FC")
}else{
	d2 = meta[,c(sample,"VIP","FDR","Fold.Change")]
	colnames(d2)=c(sample,"VIP","P.value","Ratio")
	d2$log2FC = log2(d1$Ratio)
}

d2$log10p=-log10(d2$P.value)
d2$Regulation='unchange'
max(d2$log10p,na.rm=T)
max(abs(d2$log2FC),na.rm=T)
ymax=max(d2$log10p,na.rm=T)+1
xmax=max(abs(d2$log2FC),na.rm=T)+1
if(ymax > 10){
        ymax <- 10
}
if(xmax >6){
        xmax=6
}
d2$Regulation[d2$VIP>args$meta_vip&d2$P.value<args$meta_fdr&d2$log2FC>args$meta_log2fc]='up'
d2$Regulation[d2$VIP>args$meta_vip&d2$P.value<args$meta_fdr&d2$log2FC<(-args$meta_log2fc)]='down'
p2=ggplot(data=d2,aes(x=log2FC,y=log10p,colour=Regulation))+geom_point(aes(size=VIP))+
  scale_colour_manual(values=c('unchange'="gray",'down'="#4575B4",'up'="#D73027"))+
  labs(x='Log2FC metabolites',y="-log10 FDR")+theme_bw()+ylim(0,ymax)+xlim(-xmax,xmax)
ggarrange(p1,p2,labels='AUTO',nrow=1)
ggsave(paste0(figure_path,"/Figure 2. volcano.pdf"),width=12,height=8)
ggsave(paste0(figure_path,"/Figure 2. volcano.png"),width=12,height=8)

print("Finish figure2")
Sys.time()-t0
t0 = Sys.time()

####plot figure 1
statistics = reshape2::melt(cbind(Proteome = table(d1$Regulation),Metabolism = table(d2$Regulation))[c("up","down"),],varnames=c("Type","Sample"),value.name = "Number")
maxn=max(statistics$Number)
maxn=ceiling(maxn*1.1)
ggplot(statistics,aes(x=Sample,y=Number,fill=Type))+geom_bar(stat="identity",width=0.7,position=position_dodge(width =1))+theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1),axis.text=element_text(size=10,color='black',face='bold'),axis.title=element_text(size=12,color='black',face='bold'),strip.text=element_text(size=10,color='black',face='bold'),axis.title.x=element_text(vjust=-5),axis.title.y=element_text(vjust=5),plot.margin=unit(c(0,0,2,1),units='cm'),legend.title=element_text(size=12),legend.text=element_text(size=10))+labs(x='',y='Number of proteins and metabolites')+geom_text(aes(label=Number),position=position_dodge(1),color='black',vjust=-0.2,size=4.5)+scale_fill_brewer(palette="Set1")+theme(legend.position='right',legend.direction='vertical')+ylim(0,maxn)
ggsave(paste0(figure_path,"/Figure 1.Summary.pdf"),width=6,height=6)
ggsave(paste0(figure_path,"/Figure 1.Summary.png"),width=6,height=6)

print("Finish figure1")
Sys.time()-t0
t0 = Sys.time()

###plot figure 3
protein.annotation = read.table(args$protein_annotation,header=T,row.names=1,check.names=F,sep="\t",quote="",dec=".",fill=T)
metabolite.annotation = read.table(args$meta_annotation,header=T,row.names=1,check.names=F,sep="\t",quote="",dec=".",fill=T)
diff_data = rbind(protein.annotation[rownames(d1)[d1$Regulation != "unchange"],],metabolite.annotation[rownames(d2)[d2$Regulation != "unchange"],])
write.table(diff_data,"combine_protein_metabolites_KEGG_annotation.txt",sep="\t",quote=F)
system("perl /home/swf/bin/metabolism/2.pl combine_protein_metabolites_KEGG_annotation.txt >metabolite_and_protein_pathway_classify.txt")

data = read.table("metabolite_and_protein_pathway_classify.txt",header=T,row.names=1,check.names=F,sep="\t",quote="",dec=".",stringsAsFactors=F)
data1 =t(data[,c(1,2)])
w = max(colSums(data1))
h=ncol(data1)
n = max(sapply(colnames(data1),nchar))
pdf(paste0(figure_path,"/Figure 3.Differential Metabolite and Protein pathway classify.pdf"),width=w/8+5,height=h/5+5)
par(mar=c(5,n/3,2,5),las=2)
b = barplot(as.matrix(data1),horiz = T,col=c("lightblue","orange"),cex.names = 0.75,border = F,xlim = c(0,w*1.1))
legend("topright",legend = rownames(data1),pch=19,col = c("lightblue","orange"))
text(data1[1,]/2,b,labels =data1[1,])
text((data1[2,]/2+data1[1,]),b,labels = data1[2,])
dev.off()

print("Finish figure3")
Sys.time()-t0
t0 = Sys.time()

#####write Table 1
wb1=createWorkbook()
addWorksheet(wb1,"Statistics")
writeData(wb1,"Statistics",x=cbind(Proteome = table(d1$Regulation),Metabolism = table(d2$Regulation))[c("up","down"),],colNames = T,rowNames = T)
addWorksheet(wb1,"Protein")
hs1 = createStyle(fgFill = "grey", halign = "CENTER", valign = "center",textDecoration = "bold",border = "TopBottomLeftRight",wrapText = T)
writeData(wb1,"Protein",data.table(d1,keep.rownames = T),rowNames = F,colNames = T,headerStyle = hs1,keepNA = F)
addWorksheet(wb1,"Metabolite")
writeData(wb1,"Metabolite",data.table(d2,keep.rownames=T),rowNames = F,colNames = T,headerStyle = hs1,keepNA = F)
addWorksheet(wb1,"Diff Protein")
writeData(wb1,"Diff Protein",x=data.table(d1[d1$Regulation != "unchange",],keep.rownames = T),rowNames = F,colNames = T,headerStyle = hs1,keepNA = F)
addWorksheet(wb1,"Diff Meta")
writeData(wb1,"Diff Meta",x=data.table(d2[d2$Regulation != "unchange",],keep.rownames = T),rowNames = F,colNames = T,headerStyle = hs1,keepNA = F)
addWorksheet(wb1,"KEGG annotation")
writeData(wb1,"KEGG annotation",x=data.table(diff_data,keep.rownames=T),rowNames = F,colNames = T,headerStyle = hs1,keepNA = F)
saveWorkbook(wb1,paste0(table_path,"/Table 1. Data summary.xlsx"),overwrite=T)

print("Finish table1")
Sys.time()-t0
t0 = Sys.time()

####write table 2
data=data.table(data,keep.rownames = T)
names(data)=c("Pathway","Number of proteins","Number of metabolites","Protein IDs","Metabolites IDs")
wb=createWorkbook()
hs2 = createStyle(fgFill = "lightblue", halign = "CENTER", valign = "center",textDecoration = "bold",border = "TopBottomLeftRight",wrapText = T)
bs2 = createStyle(halign="left",valign="center",fgFill = "#8FBC8F",border="TopBottomLeftRight")
addWorksheet(wb,"pathway")
writeData(wb,"pathway",x = data[,c(1,2,3)],keepNA = F,rowNames = F,headerStyle = hs2)
addStyle(wb,"pathway",bs2,rows=2:(nrow(data)+1),cols=1:3,gridExpand = TRUE)
setColWidths(wb,"pathway",1:3,widths = "auto")
setRowHeights(wb, "pathway", rows = 1, heights = 40)

addWorksheet(wb,"pathway detail")
writeData(wb,"pathway detail",x=data,keepNA = F,rowNames = F,headerStyle = hs2)
addStyle(wb,"pathway detail",bs2,rows=2:(nrow(data)+1),cols=1:5,gridExpand = TRUE)
setColWidths(wb,"pathway detail",1:5,widths = "auto")
setRowHeights(wb, "pathway detail", rows = 1, heights = 40)
saveWorkbook(wb,file = paste0(table_path,"/Table 2. Metabolite and protein pathway classify.xlsx"),overwrite = T)

print("Finish table2")
Sys.time()-t0
t0 = Sys.time()

setwd(kegg_path)
all_diff = rbind(d1[d1$Regulation != "unchange",c("Ratio","P.value","Regulation")],d2[d2$Regulation !="unchange",c("Ratio","P.value","Regulation")])
write.table(all_diff,"protein_metabolites_diff_quant.txt",sep="\t",quote=F)
system(paste0("perl /home/swf/bin/metabolism/KEGG_pathway_detail_itraq_protein_new.pl ",annotation_dir,"/all_pathway.txt protein_metabolites_diff_quant.txt ",annotation_dir,"/all_list.txt ",paste0(gpath,"/metabolite_and_protein_pathway_classify.txt "),args$organism))

print("Finish kegg pathway plot")
Sys.time()-t0
t0 = Sys.time()

####enrich bubble plot figure5 and table3

setwd(gpath)
d1$Regulation[d1$Regulation=="unchange"]=""
d2$Regulation[d2$Regulation=="unchange"]=""
write.table(data.table(d1[,c("Ratio","P.value","Regulation")],keep.rownames = T),"protein_diff.list",sep="\t",quote=F,row.names=F)
write.table(data.table(d2[,c("Ratio","P.value","Regulation")],keep.rownames = T),"metabolite_diff.list",sep="\t",quote=F,row.names=F)

system(paste0("perl /home/swf/bin/metabolism/enrich_kegg_new.pl ",args$protein_annotation," protein_diff.list 1 1.5 Protein"))
system(paste0("perl /home/swf/bin/metabolism/enrich_kegg_new.pl ",args$meta_annotation," metabolite_diff.list 1 1.5 Metabolite"))
system("cat Protein_pathway_enrichment.xls Metabolite_pathway_enrichment.xls >pathway_enrich.txt")

enrich_data = read.table("pathway_enrich.txt",header=F,sep="\t",comment.char = "",quote="")
colnames(enrich_data)=c("Type","KEGG.pathway","Mapping","Background","All.Mapping","All.Background","Fold.enrichment","p.value","Related.proteins")
getPalette = colorRampPalette(brewer.pal(3, "Set1"))
enrich_data=subset(enrich_data,enrich_data$p.value<0.05)
h = nrow(enrich_data)
enrich_data$Fold.enrichment=log2(enrich_data$Fold.enrichment)
ggplot(enrich_data,aes(x=Fold.enrichment,y=KEGG.pathway,color=p.value))+
  geom_point(aes(size=Mapping,shape=Type))+
  scale_color_gradientn(colours = getPalette(3))+
  scale_size_area(trans='log2')+theme_bw()
ggsave(file=paste0(figure_path,"/Figure 5. Metabolite and protein pathway enrichment.pdf"),width=8,height=5+h/10)
ggsave(file=paste0(figure_path,"/Figure 5. Metabolite and protein pathway enrichment.png"),width=8,height=5+h/10)

wb3 = createWorkbook()
addWorksheet(wb3,"Pathway enrich")
hs3 = createStyle(halign = "CENTER", valign = "center",textDecoration = "bold",wrapText = T)
writeData(wb3,"Pathway enrich",enrich_data,rowNames = F,colNames = T,headerStyle = hs3)
setColWidths(wb3,"Pathway enrich",cols = 1:ncol(enrich_data),widths = "auto")
saveWorkbook(wb3,paste0(table_path,"/Table 3. Metabolite and protein pathway enrichment.xlsx"),overwrite=T)

print("Finish figure5 and table3 kegg enrich")
Sys.time()-t0
t0 = Sys.time()

######spearman heatmap figure 6
protein=d1[d1$Regulation != "",sample]
submeta = d2[d2$Regulation != "",sample]
subprotein=protein[rowSums(is.na(protein))<length(sample)/2,]
cor=matrix(0,nrow(subprotein),nrow(submeta))
cor=cor(t(subprotein),t(submeta),method = 'spearman',use='pairwise.complete.obs')

print("finish cor test")
Sys.time()-t0
t0 = Sys.time()

h = nrow(cor)
w = ncol(cor)
p3=pheatmap(cor,color = colorRampPalette(c("#024e68","#06799f","white","#ffa700","#fb000d"))(100),silent = T,fontsize_row = h/(100+h)+4,fontsize_col = w/(100+w)+4)
ggsave(plot=p3,file=paste0(figure_path,"/Figure 6. Metabolite and protein correaltion heatmap.pdf"),width=w/40+5,height = h/40+5)
ggsave(plot=p3,file=paste0(figure_path,"/Figure 6. Metabolite and protein correaltion heatmap.png"),width=w/40+5,height = h/40+5)

print("Finish figure6")
Sys.time()-t0
t0 = Sys.time()

########network for figure7
node1=array()
node2=array()
score=array()
t=0
for(i in 1:nrow(cor)){
  for(j in 1:ncol(cor)){
    if(abs(cor[i,j])>0.8){
      t=t+1
      node1[t]=rownames(cor)[i]
      node2[t]=colnames(cor)[j]
      score[t]=cor[i,j]
    }
  }
}
nodedata = data.table(name=c(unique(node1),unique(node2)),Type=c(rep("Protein",length(unique(node1))),rep("Metabolism",length(unique(node2)))))
write.table(file="node.txt",nodedata,sep="\t",quote=F,row.names=F)
network=data.table(node1=node1,node2=node2,Pearson.correlation=score)
network[,type := ceiling(Pearson.correlation)]
write.table(file="network.txt",network,sep="\t",quote=F,row.names = F)
wb4=createWorkbook()
addWorksheet(wb4,"network edge")
hs4 = createStyle(halign = "CENTER", valign = "center",textDecoration = "bold",wrapText = T)
writeData(wb4,"network edge",network,rowNames = F,colNames = T,headerStyle = hs4)
saveWorkbook(wb4,file=paste0(table_path,"/Table 4. Metabolite and protein co-expression network.xlsx"),overwrite = T)

print("Finish table4 network")
Sys.time()-t0

