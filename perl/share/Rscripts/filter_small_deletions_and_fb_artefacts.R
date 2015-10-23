rg_file = commandArgs(T)[1]
full_bam_insert_size_distribution = commandArgs(T)[2]
foldback_artefact_classification_file = commandArgs(T)[3]
output_file_name = commandArgs(T)[4]

d = read.table(rg_file, header = F, sep = "\t", stringsAsFactors = F, comment = "#")

isize = read.table(full_bam_insert_size_distribution, header = F, sep = " ")
quantile.from.freq = function(vals, freq, q) {
    freq = ifelse(is.na(freq), 0, freq)
    o = order(vals)
    vals = vals[o]
    freq = freq[o]
    cs <- cumsum(freq)
    total = sum(freq)
    sapply(
        q,
        function(qq) vals[max(which(cs < qq * total)) + 1]
    )
}
insert_size_threshold = quantile.from.freq(isize[,1], isize[,2], 0.99)
is_small_deletion_artefact = d[,1] == d[,4] & d[,9] == "+" & d[,10] == "-" & d[,5] - d[,3] <= insert_size_threshold

foldback_classifications = read.table(foldback_artefact_classification_file, header = F, sep = "\t")
is_foldback_artefact = d[,7] %in% foldback_classifications[!foldback_classifications[,3], 1]

idx = !is_small_deletion_artefact & !is_foldback_artefact
write.table(
    d[idx, ],
    output_file_name,
    row.names = F,
    col.names = F,
    sep = "\t",
    quote = F
)
