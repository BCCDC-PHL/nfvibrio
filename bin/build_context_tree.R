library(ggtree)
library(ggplot2)
library(ape)
library(stringr)
library(grid)

args = commandArgs(trailingOnly = T)

if (length(args) != 2 && length(args) != 3){
    print("USAGE: Rscript build_context_tree.R INPUT_TREE OUTPUT_PDF [HEADER_DELIM]")
    quit()
}

if (length(args) == 2){
  delim = '|'
}else{
  delim = args[3]
}

delim = '|'
# 1: input tree
# 2: metadata file indicating target / background 
# 3: output PDF plot 

# setwd("/Users/jpalmer/work/norovirus/references/G9-test")

tree = read.tree(args[1])
# tree = read.tree('~/work/norovirus/results-jan-31-9/phylo/capsid/tree/jan-31-9-gtype.nwk')

print("subsetting...")
get_mask <- function(labels){
  return(grepl("^C.+$", labels))
}

target_mask = get_mask(tree$tip.label)
target_labels = tree$tip.label[target_mask]
background_labels = tree$tip.label[!target_mask]

branches <- list(Clinical=target_labels, Environmental=background_labels)
tree <- groupOTU(tree, branches)
print(tree)
print("ggtree")
p <- ggtree(tree, aes(color=group)) +
  geom_tiplab(align=T, linesize=0.25) + xlim(NA,0.5) +   
  scale_color_manual(name="Sequence Type", values=c(Clinical = "darkblue",Environmental= "red")) #+ 
  #guides(color = guide_legend(override.aes = list(size = 4))) 

print("ggsave")
# pdf(file = args[2], width = 12, height = 8)

ggsave(args[2], device="pdf",   width = 10, height = 65, limitsize=F)
# dev.off()


#tree$tip.label = sapply(tree$tip.label, function(x){ if (grepl("\\|", x)){ paste(strsplit(x, '\\|')[[1]][2:3], collapse="_")} else {x}})

