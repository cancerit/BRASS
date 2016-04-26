library(stringr)

# MAX_FOLDBACK_DIST = 5000
MIN_DIST_OF_CN_SEG_BKPT_TO_RG = 20000
MAX_GET_READS_EXTEND_DIST = 10000
args = commandArgs(trailingOnly = T)
options(scipen = 200)
MIN_WINDOW_BIN_COUNT = 10

if (length(args) < 6) {
    cat("Usage: \n")
    cat("  Rscript get_rg_cns.R rgs_file cn_file segments_file bam_file acf cent_file\n")
    stop()
} else {
    cat(sprintf("Using following settings:\nMIN_DIST_OF_CN_SEG_BKPT_TO_RG = %d\nMAX_GET_READS_EXTEND_DIST = %d\nMIN_WINDOW_BIN_COUNT = %d\n", MIN_DIST_OF_CN_SEG_BKPT_TO_RG, MAX_GET_READS_EXTEND_DIST, MIN_WINDOW_BIN_COUNT), file = stderr())
}

cat("Reading in data...\n")
rgs_file = args[1]
cn_file = args[2]
segs_file = args[3]
bam_file = args[4]
acf = as.numeric(args[5])
cent_file = args[6]
gender_chr = args[7]
gender_present = args[8]


# Telomeric and centromeric regions are determined based on standard coordinates
# + surrounding unmappable regions. The matrix has coordinates for where
# p-tel ends, regions of centromeres and where q-tel starts.
centromere_telomere_coords<-as.matrix(
                              read.table( cent_file,
                                          sep="\t",
                                          header=TRUE,
                                          stringsAsFactor=FALSE,
                                          row.names = 1,
                                          fill = TRUE)[,1:4] # this drops the comment column
                                       )

coord_within_tel_or_cent = function(chr, pos) {
    mapply(
        function(c, p) {
            temp = centromere_telomere_coords[c, ]
            olap_ptel = p <= temp[1]
            olap_cen  = p >= temp[2] && p <= temp[3]
            olap_qtel = p >= temp[4]
            return(olap_ptel || olap_cen || olap_qtel)
        },
        chr,
        pos
    )
}


if (acf > 1) acf = acf * 0.01

if (file.info(rgs_file)$size > 0) {
    rgs = read.table(rgs_file, header = F, sep = "\t", stringsAsFactors = F, comment = "")
} else {
    rgs = data.frame(
        chr_l = c(),
        l5 = c(),
        l3 = c(),
        chr_h = c(),
        h5 = c(),
        h3 = c(),
        id = c(),
        score = c(),
        stl = c(),
        sth = c(),
        pos_l = c(),
        pos_h = c()
    )
}
if (!("cn" %in% ls())) {
    cn = read.table(cn_file, header = F, sep = "\t", colClasses = c("character", rep("numeric", 3)))
}
segs = read.table(segs_file, header = F, sep = "\t", colClasses = c("character", rep("numeric", 3)))


cat("Organizing and sorting data...\n")
bkpts_of_chr  = list()
rg_end_of_chr = list()
rg_idx_of_chr = list()

chrs<-rownames(centromere_telomere_coords);

for (c in chrs) {
    bkpts_of_chr[[c]] = numeric()
    rg_end_of_chr[[c]] = numeric()
    rg_idx_of_chr[[c]] = numeric()
}

if (nrow(rgs) > 0 ) {
    for (i in 1:nrow(rgs)) {
        # skip male gender chr if it is not indicated as present by ASCAT summary
        if(gender_present == "N") {
          if (rgs[i, 1] == gender_chr || rgs[i, 4] == gender_chr) next()
        }

        # Skip rearrangements where one end falls onto centromere or telomere
        if (any(coord_within_tel_or_cent(c(rgs[i, c(1,4)]), c(rgs[i, 11:12])))) {
            next()
        }

        l = rgs[i, 11]
        h = rgs[i, 12]

        bkpts_of_chr[[rgs[i,1]]] = append(bkpts_of_chr[[rgs[i,1]]], l)  # Breakpoint positions
        bkpts_of_chr[[rgs[i,4]]] = append(bkpts_of_chr[[rgs[i,4]]], h)
        rg_end_of_chr[[rgs[i,1]]] = append(rg_end_of_chr[[rgs[i,1]]], 1)  # Low/high end
        rg_end_of_chr[[rgs[i,4]]] = append(rg_end_of_chr[[rgs[i,4]]], 2)
        rg_idx_of_chr[[rgs[i,1]]] = append(rg_idx_of_chr[[rgs[i,1]]], i)  # Generic index for this rearrangement
        rg_idx_of_chr[[rgs[i,4]]] = append(rg_idx_of_chr[[rgs[i,4]]], i)
    }
}


# Order the rearrangement breakpoints falling into each chromosome.
for (chr in names(bkpts_of_chr)) {
    o = order(bkpts_of_chr[[chr]])
    bkpts_of_chr[[chr]] = bkpts_of_chr[[chr]][o]
    rg_end_of_chr[[chr]] = rg_end_of_chr[[chr]][o]
    rg_idx_of_chr[[chr]] = rg_idx_of_chr[[chr]][o]
}

# Go through every chromosome and perform segmentation
cat("Segmenting data...\n")
seg_chr = character()
seg_start_coord = numeric()
seg_end_coord = numeric()
seg_start_bkpt = character()
seg_end_bkpt = character()
seg_cn = numeric()
seg_bin_count = numeric()
for (c in chrs) {
  cat(paste0("  Chr ", c, "...\n"))

    cn_idx = cn[,1] == c
    pos = cn[cn_idx, 2]/2 + cn[cn_idx, 3]/2

    segs.f = segs[segs[,1] == c, , drop = F]

    if (nrow(segs.f) > 1) {
        cn_seg_bkpts = segs.f[-1, 2]
        if (length(bkpts_of_chr[[c]]) > 0) {
            far_from_rg_bkpts = sapply(
                cn_seg_bkpts,
                function(x) min(abs(x - bkpts_of_chr[[c]])) > MIN_DIST_OF_CN_SEG_BKPT_TO_RG
            )
        }
        else {
            far_from_rg_bkpts = T
        }
        not_bordering_small_seg = segs.f[-1, 5] >= MIN_WINDOW_BIN_COUNT & segs.f[-nrow(segs.f), 5] >= MIN_WINDOW_BIN_COUNT
        has_big_cn_change = abs(segs.f[-1, 4] - segs.f[-nrow(segs.f), 4]) > 0.3
        cn_seg_bkpts = cn_seg_bkpts[far_from_rg_bkpts & not_bordering_small_seg & has_big_cn_change]

        # Remove copy number breakpoints that are within centromeres or telomeres
        bad_bkpt_idx = coord_within_tel_or_cent(
            rep(c, length(cn_seg_bkpts)),
            cn_seg_bkpts
        )
        cn_seg_bkpts = cn_seg_bkpts[!bad_bkpt_idx]


        # Add the centromeric coordinates to the boundaries
        boundaries = c(
            bkpts_of_chr[[c]],  # Rearrangement breakpoint positions
            cn_seg_bkpts,       # (Remaining) copy number breakpoint positions
            centromere_telomere_coords[c, 2:3][centromere_telomere_coords[c, 2:3] != 0]  # Centromere boundaries
        )
        bkpts = c(
            bkpts_of_chr[[c]],
            rep("cn_bkpt", length(cn_seg_bkpts)),
            rep("centromere", sum(!(centromere_telomere_coords[c, 2:3] == 0)))
        )
        rg_idx = c(
            rgs[rg_idx_of_chr[[c]], 7],
            rep("cn_bkpt", length(cn_seg_bkpts)),
            rep("centromere", sum(!(centromere_telomere_coords[c, 2:3] == 0)))
        )
        rg_end = c(
            rg_end_of_chr[[c]],
            rep("cn_bkpt", length(cn_seg_bkpts)),
            rep("centromere", sum(!(centromere_telomere_coords[c, 2:3] == 0)))
        )
    } else {
        boundaries = c(
            bkpts_of_chr[[c]],
            centromere_telomere_coords[c, 2:3][centromere_telomere_coords[c, 2:3] != 0]
        )
        bkpts = c(
            bkpts_of_chr[[c]],
            rep("centromere", sum(!(centromere_telomere_coords[c, 2:3] == 0)))
        )
        rg_idx = c(
            rgs[rg_idx_of_chr[[c]], 7],
            rep("centromere", sum(!(centromere_telomere_coords[c, 2:3] == 0)))
        )
        rg_end = c(
            rg_end_of_chr[[c]],
            rep("centromere", sum(!(centromere_telomere_coords[c, 2:3] == 0)))
        )
    }

    o = order(boundaries)
    boundaries = boundaries[o]

    bkpts = bkpts[o]
    rg_idx = rg_idx[o]
    rg_end = rg_end[o]
    boundaries = c(
        centromere_telomere_coords[c, "ptel"],  # Coordinate where mappable region for p-telomere starts
        boundaries,
        centromere_telomere_coords[c, "qtel"]   # Coordinate where mappable region for q-telomere starts
    )

    seg_chr = append(
        seg_chr,
        rep(c, length(boundaries) - 1)
    )
    seg_start_coord = append(
        seg_start_coord,
        boundaries[-length(boundaries)]
    )
    seg_end_coord = append(
        seg_end_coord,
        boundaries[-1]
    )
    seg_start_bkpt = append(
        seg_start_bkpt,
        c(
            paste0(c, "_ptel"),
            ifelse (
                rg_idx %in% c("cn_bkpt", "centromere"),
                rg_idx,
                paste0(rg_idx, ":", rg_end)
            )
        )
    )
    seg_end_bkpt = append(
        seg_end_bkpt,
        c(
            ifelse (
                rg_idx %in% c("cn_bkpt", "centromere"),
                rg_idx,
                paste0(rg_idx, ":", rg_end)
            ),
            paste0(c, "_qtel")
        )
    )
    temp_seg_cn = sapply(
        by(
            cn[cn_idx, 4],
            findInterval(pos, boundaries),
            function(x) c(x)
        ),
        function(x) {
            if (length(x) < MIN_WINDOW_BIN_COUNT) NA else median(x)
        }
    )[as.character(1:(length(boundaries)-1))]
    if (centromere_telomere_coords[c, 2] == 0) {
        temp_seg_cn[1] = NA  # Remove CN of first segment of an acrocentric chromosome?
    }
    seg_cn = append(
        seg_cn,
        temp_seg_cn
    )
    seg_bin_count = append(
        seg_bin_count,
        sapply(
            by(
                cn[cn_idx, 4],
                findInterval(pos, boundaries),
                function(x) c(x)
            ),
            length
        )[as.character(1:(length(boundaries)-1))]
    )

}

# Go through all segments and attempt to get copy number for small segments
# directly from the BAM file.
i = 1
cat("Estimating copy number in short segments...\n")
prev_chr = "0"
while (i <= length(seg_chr)) {
    if (seg_chr[i] == "Y") {
        i = i + 1
        next()
    }

    # For debugging
    # if (seg_chr[i] != "5") { i = i + 1; next() }

    if (seg_chr[i] != prev_chr) {
        cat(paste0("  Chr ", seg_chr[i], "...\n"))
        prev_chr = seg_chr[i]
    }

    if (!is.na(seg_cn[i])) {
        i = i + 1
        next
    }

    # Now a stretch of CN = NA starts. Find end point of the NA stretch
    j = i
    while (j < length(seg_cn) && is.na(seg_cn[j+1]) && seg_chr[i] == seg_chr[j+1]) { j = j + 1 }

    # Check if there's a non-NA CN segment on either side of the NA stretch
    left_side_not_NA = i > 1 && seg_chr[i-1] == seg_chr[i] && !is.na(seg_cn[i-1])
    rght_side_not_NA = seg_chr[j+1] == seg_chr[i] && !is.na(seg_cn[j+1]) && j < length(seg_cn)

    if (!left_side_not_NA && !rght_side_not_NA) {
        i = j+1
        next
    }

    # Check if any rearrangement breakpoints involved
    rgs_involved = F
    for (x in i:j) {
        if ( !(seg_start_bkpt[x] %in% c("cn_bkpt", "centromere")) || !(seg_end_bkpt[x] %in% c("cn_bkpt", "centromere"))) {
            rgs_involved = T
        }
    }

    # Exception: the first (empty) segment of an acrocentric chromosome
    if (seg_start_bkpt[i] == paste0(seg_chr[i], "_ptel") && seg_end_bkpt[i] == "centromere" && centromere_telomere_coords[seg_chr[i], 1] == 0) {
        rgs_involved = F
    }

    if (!rgs_involved) {
        i = j+1
        next
    }

    cat(paste0("    ", seg_start_coord[i], " ", seg_end_coord[j], "\n"))

    # Get the BAM file subset out first
    loc = paste0(
        seg_chr[i], ":",
        (if (left_side_not_NA) pmax(seg_start_coord[i-1], seg_start_coord[i]-MAX_GET_READS_EXTEND_DIST) else seg_start_coord[i]), "-",
        (if (rght_side_not_NA) pmin(seg_end_coord[j+1], seg_end_coord[j]+MAX_GET_READS_EXTEND_DIST) else seg_end_coord[j])
    )

    cmd = paste0(
        "samtools view -hb -q 1 -F 3852 -f 2 ",
        bam_file, " ",
        loc,
        " > ",
        bam_file, ".subset"
    )
    system(cmd)

    # Create a BED file for each segment
    k = j - i + 3
    if (!left_side_not_NA || !rght_side_not_NA) { k = k - 1 }
    bed_chrs = rep(seg_chr[i], k)

    bed_starts = c(
        (if (left_side_not_NA) pmax(seg_start_coord[i-1], seg_start_coord[i]-MAX_GET_READS_EXTEND_DIST) else c()),
        (if (rght_side_not_NA) seg_start_coord[i:(j+1)] else seg_start_coord[i:j])
    )
    bed_ends = c(
        (if (left_side_not_NA) seg_end_coord[(i-1):j] else seg_end_coord[i:j]),
        (if (rght_side_not_NA) pmin(seg_end_coord[j+1], seg_end_coord[j]+MAX_GET_READS_EXTEND_DIST) else c())
    )
    write.table(data.frame(bed_chrs, bed_starts, bed_ends), paste0(bam_file, ".subset.segs"), col.names = F, row.names = F, quote = F, sep = "\t")
    seg_coords = paste(bed_starts, bed_ends, sep = "-")

    # Compute average coverage for each segment
    cmd = paste0(
        "bedtools coverage -d -abam ",
        bam_file, ".subset ",
        "-b ", bam_file, ".subset.segs ",
        "| bedtools groupby -g 1,2,3 -c 5 -o mean | sort -k2,2n "
    )

    res = system(cmd, intern = T)
    res = t(as.data.frame(strsplit(res, "\t")))
    coverage_of_coord = as.numeric(res[,4]) + 0.01
    names(coverage_of_coord) = paste(res[,2], res[,3], sep = "-")
    coverages = coverage_of_coord[seg_coords]
    coverages[c(-1, -length(coverages))] = ifelse(
        bed_starts[-c(1, length(bed_starts))] - bed_ends[-c(1, length(bed_ends))] >= -1,
        NA,
        coverages[-c(1, length(coverages))]
    )

    # Estimate copy numbers for NA regions
    if (!left_side_not_NA) {
        o2 = acf * seg_cn[j+1] + (1 - acf) * 2
        cn_estimates = ( (coverages[-length(coverages)] * o2 / coverages[length(coverages)]) - (1-acf)*2 )/acf
    } else if (!rght_side_not_NA) {
        o1 = acf * seg_cn[i-1] + (1 - acf) * 2
        cn_estimates = ( (coverages[-1] * o1 / coverages[1]) - (1-acf)*2 )/acf
    } else {
        o1 = acf * seg_cn[i-1] + (1 - acf) * 2
        o2 = acf * seg_cn[j+1] + (1 - acf) * 2
        cn_estimates = sqrt(
            pmax(0, ( (coverages[c(-1, -length(coverages))] * o1 / coverages[1]) - (1-acf)*2 )/acf ) *
            pmax(0, ( (coverages[c(-1, -length(coverages))] * o2 / coverages[length(coverages)]) - (1-acf)*2 )/acf )
        )
    }

    if (length(i:j) != length(cn_estimates)) stop("Length of i:j != length of cn_estimates")
    seg_cn[i:j] = cn_estimates

    i = i + 1
}

seg_bin_count[is.na(seg_bin_count)] = 1

write.table(
    data.frame(
        seg_chr,
        seg_start_coord,
        seg_end_coord,
        seg_cn,
        seg_start_bkpt,
        seg_end_bkpt,
        seg_bin_count
    ),
    paste0(cn_file, ".rg_cns"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
)


# Cleanup
file.remove(paste0(bam_file, ".subset"))
file.remove(paste0(bam_file, ".subset.segs"))
