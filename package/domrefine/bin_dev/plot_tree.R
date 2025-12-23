#!/home/chiba/bin/rscript -n
#!/home/chiba/bin/rscript -ng
USAGE =
"Usage: 
-o NAME: file name prefix
-c CEX: character size
-l: show node label
-u: unroted tree
-h: help

-r NODES: red
-b NODES: bule
-v NODES: violet
-V NODES: dark violet
-g NODES: green
-G NODES: dark green
-y NODES: yellow
-Y NODES: dark yellow
-t NODES: tomato
-a NODES: aqua
-O NODES: dark orange
-p NODES: pink
-1 NODES: gray
-2 NODES: dark gray
"
### Get command line argumants ###
getopts("o:c:luhr:b:v:V:g:G:y:Y:t:a:O:p:P:1:2:")
if (opt.h) {
  cat(USAGE)
  q()
}
##################################

library(ape)

tree = read.tree(STDIN)

### Settings ###
cex = 1
if (length(tree$tip.label) > 200) {
  cex = 0.1;
} else if (length(tree$tip.label) > 70) {
  cex = 0.3;
} else if (length(tree$tip.label) > 30) {
  cex = 0.5;
}
if (opt.c != F) {
  cex = as.numeric(opt.c)
}

if (! is.na(ARGV[1])) {
  pdf(paste(ARGV[1], ".pdf", sep=""))
}

# type="cladogram",
PLOT_TYPE = "phylogram"
if (opt.u != F) {
  PLOT_TYPE = "unrooted"
}


LEAF_COLOR = rep("black", length(tree$tip.label))
color.leaves=function(leaf_color, leaves, color){
  if (leaves != F) {
    for (label in unlist(strsplit(leaves, ","))) {
      for (i in 1:length(tree$tip.label)) {
        if (tree$tip.label[i] == label) {
          leaf_color[i] = color
        }
      }
    }
  }
  return(leaf_color)
}
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.r, "red");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.b, "blue");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.v, "darkorchid1");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.V, "darkviolet");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.g, "green");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.G, "forestgreen");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.y, "yellow");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.Y, "yellow3");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.t, "tomato");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.a, "skyblue");
# LEAF_COLOR = color.leaves(LEAF_COLOR, opt.a, "cornflowerblue");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.O, "darkorange1");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.p, "hotpink");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.P, "deeppink");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.1, "gray60");
LEAF_COLOR = color.leaves(LEAF_COLOR, opt.2, "gray45");

EDGE_COLOR = rep("black", dim(tree$edge)[1])

# tree$node.label
# tree$edge

plot(tree,
     type=PLOT_TYPE,
     show.node.label=opt.l,
     edge.color = EDGE_COLOR,
     tip.color = LEAF_COLOR,
     cex=cex)

# nodelabels(tree$node.label,
#            frame="none",
# #            adj=c(0,1.2),
# #            adj=c(1,0),
#            bg="none"
#            )

# edgelabels(text=tree$edge.length
#            ,frame="none"
#            ,bg="none"
#            ,col="darkgray"
#            )
