library(ggplot2)
library(ggtranscript)
library(dplyr)
library(ggtree)
library(cowplot)

gtf = read.delim("QueSub.ggtranscript.gff", header = F, quote = "")
colnames(gtf) = c("Chr","source","type","start","end","score","strand","phase","transcript_name","gene_name","domain")
gtf$Chr = gsub("Chromosome", "Chr",gtf$Chr)
gtf$transcript_name=as.factor(gtf$transcript_name)

tree = read.tree("QueSub.pep.trim.aln.treefile")
tree_plot = ggtree(TCP_tree) + theme_tree2()

gtf_draw =
  ggplot(filter(gtf, type=="exon"),aes(xstart=start,xend=end,y=transcript_name))+
  ggtranscript::geom_range(fill = "white",height=0.25)+
  ggtranscript::geom_range(data = filter(gtf, type=="CDS"))+
  ggtranscript::geom_range(data = filter(gtf, type=="domain"&domain=="TCP"),
             fill = "orange")+
  ggtranscript::geom_range(data = filter(gtf, type=="domain"&domain!="TCP"),
             fill = "purple", alpha = 0.5)+
  geom_intron(data = to_intron(filter(gtf,type=="exon"),"transcript_name"),
              aes(strand=strand),
              arrow.min.intron.length=100)+
  scale_y_discrete(limits=rev(na.omit(tree_plot$data)$label))+
  theme_bw()+
  theme(axis.title = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,color = "black"))
ggsave("test.jpeg",width = 40,height = 40,dpi = 300)

