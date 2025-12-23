#!/home/chiba/etc/rscript/rscript -ng
USAGE =
"Usage: 
-h: help
"
### Get command line argumants ###
getopts("h123")
if (opt.h) {
  cat(USAGE)
  q()
}
##################################

f1 = read.table(ARGV[1])
if (length(ARGV) >= 2) {
  f2 = read.table(ARGV[2])
}
# f3 = read.table(ARGV[3])


if (length(ARGV) == 1) {
# plot(log10(f1$V1), log10(f1$V2), xaxt="n", yaxt="n", xlim=c(0,5.2), ylim=c(0,5.4))
  plot(log10(f1$V1), log10(f1$V2), xaxt="n", yaxt="n", xlim=c(0,max(log10(f1$V1))), ylim=c(0,max(log10(f1$V2))))
} else if (length(ARGV) >= 2) {
  plot(log10(f1$V1), log10(f1$V2), col="blue", xaxt="n", yaxt="n", xlim=c(0,max(log10(f1$V1), log10(f2$V1))), ylim=c(0,max(log10(f1$V2), log10(f2$V2))))
  points(log10(f2$V1), log10(f2$V2), col="red")
}
# points(log10(f3$V1), log10(f3$V2), col="green")

log.axis(s=1, max=100000)
log.axis(s=2, max=100000)

if (opt.1 || opt.2 || opt.3) {
  if (opt.1) {
    x = f1$V1
    y = f1$V2
  } else if (opt.2) {
    x = f2$V1
    y = f2$V2
  } else if (opt.3) {
    x = f3$V1
    y = f3$V2
  }
  lm = lm(log10(y)~log10(x))
  print(lm)
  print(summary(lm))
  abline(lm)
}
