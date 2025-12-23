#!/home/chiba/bin/rscript -g
USAGE =
"Usage: 
-l: lowess
-f LOWESS_PARAMETER (default: 0.05)
-h: help
"
### Get command line argumants ###
getopts("hlf:zm")
if (opt.h) {
  cat(USAGE)
  q()
}
##################################

### input ###
score.change = t$V1
cog.change = t$V2

if (opt.z) {
} else {
  score.change = score.change[cog.change!=0]
  cog.change = cog.change[cog.change!=0]
}
if (opt.m) {
  min.score.change = min(score.change)
  cog.change = cog.change[score.change!=min.score.change]
  score.change = score.change[score.change!=min.score.change]
}

### output ###
length(score.change)

cor.test(score.change, cog.change)

plot(score.change, cog.change
     # , xlim=c(-1,1)
     # , xlim=c(-0.45,0.3)
     )
abline(v=
       c(0
         # ,-0.05
         )
       ,lt=2)
abline(h=0
       # , col="red"
       )

lm = lm(cog.change~score.change)
lm
# abline(lm, col="red")

if (opt.l) {
  lowess.f = 0.05
  if (opt.f != F) {
    lowess.f = as.numeric(opt.f)
  }
  lines(lowess(score.change, cog.change, f=lowess.f), col="red")
}
