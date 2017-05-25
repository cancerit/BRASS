tumour_file = commandArgs(T)[1]
normal_file = commandArgs(T)[2]
ploidy = as.numeric(commandArgs(T)[3])
purity = as.numeric(commandArgs(T)[4])
output_body = commandArgs(T)[5]
cent_file = commandArgs(T)[6]
genomebuild = commandArgs(T)[7]
cytoband_file = commandArgs(T)[8]

# set some defaults
if (purity > 1) { purity = purity / 100 }
# load the chr list from the cent_file
library(data.table)
chrs <- as.vector(fread(cent_file, select=c("chr"), colClasses=c("character"))$chr)
cat(paste0("Reading in file ", tumour_file, "...\n"), file = stderr())
d.t = read.table(gzfile(tumour_file), header = F, colClasses = c("character", rep("numeric", 6)))
d.t = d.t[d.t[,1] %in% chrs, ]
cat(paste0("Reading in file ", normal_file, "...\n"), file = stderr())
d.n = read.table(gzfile(normal_file), header = F, colClasses = c("character", rep("numeric", 6)))
d.n = d.n[d.n[,1] %in% chrs, ]

# Sanity check
d.t = d.t[d.t[,4] == d.t[,3] - d.t[,2], ]
d.n = d.n[d.n[,4] == d.n[,3] - d.n[,2], ]
if (!all(d.t[,2] == d.n[,2] & d.t[,3] == d.n[,3])) {
    stop()
}

# Normalise read counts
normalisation_factor = sum(d.n[,6]) / sum(d.t[,6])
# pseudocount = 1 / min(sum(d.n[,6]), sum(d.t[,6])) * 1e8
pseudocount = 0.1
if (ploidy >= 2) {
    d.t[,6] = d.t[,6] * normalisation_factor
    d.t[,7] = d.t[,7] * normalisation_factor
} else {
    d.n[,6] = d.n[,6] / normalisation_factor
    d.n[,7] = d.n[,7] / normalisation_factor
}

logratio = log2((d.t[,6] + pseudocount) / (d.n[,6] + pseudocount))
gc = d.t[,5]
avg_t_fb_count = mean(d.t[,7] / (d.t[,6] + d.t[,7]), na.rm = T)
avg_n_fb_count = mean(d.n[,7] / (d.n[,6] + d.n[,7]), na.rm = T)
tumour_fb_ratio = (d.t[,7] + avg_t_fb_count) / (d.t[,6] + d.t[,7] + 1)
tumour_fb_ratio = pmin(tumour_fb_ratio, quantile(tumour_fb_ratio, 0.99))
normal_fb_ratio = (d.n[,7] + avg_n_fb_count) / (d.n[,6] + d.n[,7] + 1)
normal_fb_ratio = pmin(normal_fb_ratio, quantile(normal_fb_ratio, 0.99))

cat("Normalising copy number using GAM...\n", file = stderr())
library(gam)
gam_m = gam(logratio ~ s(gc) + s(tumour_fb_ratio) + s(normal_fb_ratio), data = data.frame(logratio, gc, tumour_fb_ratio, normal_fb_ratio))
normalised_logratio = logratio - predict(gam_m)
normalised_logratio = normalised_logratio - median(normalised_logratio) + log2((purity * ploidy + (1 - purity) * 2) / 2)

# Smooth and segment data per chromosome
cat("Runmed-smoothing data\n", file = stderr())
for (chr in unique(d.t[,1])) {
    idx = d.t[,1] == chr
    normalised_logratio[idx] = runmed(normalised_logratio[idx], 5)
}

cat("Segmenting copy number...\n", file = stderr())
library(copynumber)
options(scipen=999)
data_for_pcf = data.frame(chr = d.t[,1], pos = rowMeans(d.t[,2:3]), cn = normalised_logratio)
segs = pcf(data_for_pcf, gamma = 200, assembly = genomebuild , cytoband_file = cytoband_file)
out_table = data.frame(chr = segs[,2], start = segs[,4], end = segs[,5], cn = segs[,7], nwin = segs[,6])
out_table = out_table[order(out_table[,1], out_table[,2]), ]
i = 1
while (i < nrow(out_table)) {
    if (out_table[i, 1] == out_table[i+1, 1]) {
        bkpt = floor( (out_table[i, 3] + out_table[i+1, 2])/2 )
        out_table[i, 3] = bkpt
        out_table[i+1, 2] = bkpt
    }
    i = i + 1
}

# skipping diagnostinc plots for chromosome 9

if(!exists do_diagnostic) {
cat("Creating diagnostic plots...\n", file = stderr())
# Some diagnostic plots
pdf(paste0(output_body, ".diagnostic_plots.pdf"))
smoothScatter(gc, logratio, pch = ".", xlab = "GC-content", ylab = "Log-ratio")
smoothScatter(tumour_fb_ratio, logratio, pch = ".", xlab = "Tumour FB ratio", ylab = "Log-ratio")
smoothScatter(normal_fb_ratio, logratio, pch = ".", xlab = "Normal FB ratio", ylab = "Log-ratio")
plot.gam(gam_m, rugplot = F)
idx = d.t[,1] == "9"
plot(d.t[idx, 3], d.t[idx, 6], xlab = "Chr 9 position", ylab = "Read depth", pch = ".")
plot(d.t[idx, 3], normalised_logratio[idx], xlab = "Chr 9 position", ylab = "Normalised log-ratio", pch = ".")
dev.off()
}

cat("Outputting data...\n", file = stderr())

# Outputting data
write.table(
    data.frame(
        d.t[,1:3],
        ( 2^(1+normalised_logratio) - 2*(1-purity) ) / purity,  # Absolute copy number
        normalised_logratio  # Normalised copy number ratio
    ),
    paste0(output_body, ".abs_cn.bg"),
    row.names = F,
    col.names = F,
    sep = "\t",
    quote = F
)

write.table(
    data.frame(
        out_table[, 1:3],
        ( 2^(1+out_table[,4]) - 2*(1-purity) ) / purity,
        out_table[,5],
        out_table[,4]
    ),
    paste0(output_body, ".segments.abs_cn.bg"),
    row.names = F,
    col.names = F,
    sep = "\t",
    quote = F
)
