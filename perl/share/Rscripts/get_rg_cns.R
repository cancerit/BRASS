library(stringr)

# MAX_FOLDBACK_DIST = 5000
MIN_DIST_OF_CN_SEG_BKPT_TO_RG = 20000
MAX_GET_READS_EXTEND_DIST = 10000
args = commandArgs(trailingOnly = T)
options(scipen = 200)
MIN_WINDOW_BIN_COUNT = 10

if (length(args) < 4) {
    cat("Usage: \n")
    cat("  Rscript get_rg_cns.R rgs_file cn_file segments_file bam_file acf\n")
    stop()
} else {
    cat(sprintf("Using following settings:\nMIN_DIST_OF_CN_SEG_BKPT_TO_RG = %d\nMAX_GET_READS_EXTEND_DIST = %d\nMIN_WINDOW_BIN_COUNT = %d\n", MIN_DIST_OF_CN_SEG_BKPT_TO_RG, MAX_GET_READS_EXTEND_DIST, MIN_WINDOW_BIN_COUNT), file = stderr())
}


chr_lens = c(
    "1" = 249250621,
    "2" = 243199373,
    "3" = 198022430,
    "4" = 191154276,
    "5" = 180915260,
    "6" = 171115067,
    "7" = 159138663,
    "8" = 146364022,
    "9" = 141213431,
    "10" = 135534747,
    "11" = 135006516,
    "12" = 133851895,
    "13" = 115169878,
    "14" = 107349540,
    "15" = 102531392,
    "16" = 90354753,
    "17" = 81195210,
    "18" = 78077248,
    "19" = 59128983,
    "20" = 63025520,
    "21" = 48129895,
    "22" = 51304566,
    "X" = 155270560,
    "Y" = 59373566
)

# Telomeric and centromeric regions are determined based on standard coordinates
# + surrounding unmappable regions. The matrix has coordinates for where
# p-tel ends, regions of centromeres and where q-tel starts.

# centromere_coords = matrix(
#     c(
#         121535434, 142535434,  # chr1, centromere + heterochromatin
#          90545103,  91545103,  # chr2, centromere
#          90504854,  93504854,  # chr3, centromere
#          49660117,  52660117,  # chr4, centromere
#          46405641,  49405641,  # chr5, centromere
#          58830166,  61830166,  # chr6, centromere
#          58054331,  61054331,  # chr7, centromere
#          43838887,  46838887,  # chr8, centromere
#          47367679,  65367679,  # chr9, centromere + heterochromatin
#          39254935,  42254935,  # chr10, centromere
#          51644205,  54644205,  # chr11, centromere
#          34856694,  37856694,  # chr12, centromere
#                NA,  19020000,  # chr13, short arm + centromere + heterochromatin
#                NA,  19000000,  # chr14. short arm + centromere
#                NA,  20000000,  # chr15, short arm + centromere
#          35335801,  38335801,  # chr16, centromere + heterochromatin
#          22263006,  25263006,  # chr17, centromere
#          15460898,  18460898,  # chr18, centromere
#          24631782,  27731782,  # chr19, heterochromatin + centromere + heterochromatin
#          26369569,  29369569,  # chr20, centromere
#          11238129,  14288129,  # chr21, heterochromatin + centromere
#                NA,  16000000,  # chr22, short arm + centromere
#          58632012,  61632012,  # chrX, centromere
#          10104553,  13104553   # chrY, centromere
#      ),
#      ncol = 2,
#      byrow = T
# )
centromere_telomere_coords = matrix(
    c(
        # 750000, 121270001, 144840000, 249220001, # chr 1: ...
        750000, 121270001, 150000000, 249220001, # ... q-side of centromere manually set
         10000,  89330001,  95390000, 242950001, # chr 2
         60000,  90500001,  93510000, 197820001, # chr 3
         40000,  49090001,  52680001, 190910001, # chr 4
         10000,  46400001,  49440000, 180720001, # chr 5
        200000,  58770001,  61880000, 170920001, # chr 6
         80000,  58050001,  61980000, 159130001, # chr 7
        160000,  43790001,  46880000, 146300001, # chr 8
        200000,  38770001,  70990000, 141090001, # chr 9
        100000,  39150001,  42400000, 135230001, # chr 10
        190000,  51580001,  54800000, 134940001, # chr 11
        180000,  34850001,  37860000, 133840001, # chr 12
             0,         0,  19360000, 115110001, # chr 13
             0,         0,  20190000, 107290001, # chr 14
             0,         0,  20030000, 102280001, # chr 15
         80000,  35240001,  46490000,  90160001, # chr 16
             0,  22240001,  25270000,  81110001, # chr 17
        130000,  15410001,  18540000,  78010001, # chr 18
        250000,  24600001,  27740000,  59100001, # chr 19
        120000,  26290001,  29420000,  62920001, # chr 20
             0,         0,  14340000,  48100001, # chr 21 - don't use the p-arm data
             0,         0,  16850000,  51200001, # chr 22
        310000,  58500001,  61730000, 155240001  # chr X
   ),
   ncol = 4,
   byrow = T
)

rownames(centromere_telomere_coords) = c(1:22, "X")
colnames(centromere_telomere_coords) = c("ptel", "cen_start", "cen_end", "qtel")

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


cat("Reading in data...\n")
rgs_file = args[1]
cn_file = args[2]
segs_file = args[3]
bam_file = args[4]
bam_body = sub(".+/", "", bam_file)
acf = as.numeric(args[5])
if (acf > 1) acf = acf * 0.01
rgs = read.table(rgs_file, header = F, sep = "\t", stringsAsFactors = F, comment = "")
if (!("cn" %in% ls())) {
    cn = read.table(cn_file, header = F, sep = "\t", colClasses = c("character", rep("numeric", 3)))
}
segs = read.table(segs_file, header = F, sep = "\t", colClasses = c("character", rep("numeric", 3)))


cat("Organizing and sorting data...\n")
bkpts_of_chr  = list()
rg_end_of_chr = list()
rg_idx_of_chr = list()
chrs = c(1:22, "X")
for (c in chrs) {
    bkpts_of_chr[[c]] = numeric()
    rg_end_of_chr[[c]] = numeric()
    rg_idx_of_chr[[c]] = numeric()
}

for (i in 1:nrow(rgs)) {
    # Skip rearrangements where one end falls onto centromere or telomere
    if (rgs[i, 1] == "Y" || rgs[i, 4] == "Y") next()
    # if (any(coord_within_tel_or_cent(c(rgs[i, c(1,4)]), c(rgs[i, 13:14])))) {
    if (any(coord_within_tel_or_cent(c(rgs[i, c(1,4)]), c(rgs[i, 11:12])))) {
        next()
    }

    # l = rgs[i, 13]
    # h = rgs[i, 14]
    l = rgs[i, 11]
    h = rgs[i, 12]

    bkpts_of_chr[[rgs[i,1]]] = append(bkpts_of_chr[[rgs[i,1]]], l)  # Breakpoint positions
    bkpts_of_chr[[rgs[i,4]]] = append(bkpts_of_chr[[rgs[i,4]]], h)
    rg_end_of_chr[[rgs[i,1]]] = append(rg_end_of_chr[[rgs[i,1]]], 1)  # Low/high end
    rg_end_of_chr[[rgs[i,4]]] = append(rg_end_of_chr[[rgs[i,4]]], 2)
    rg_idx_of_chr[[rgs[i,1]]] = append(rg_idx_of_chr[[rgs[i,1]]], i)  # Generic index for this rearrangement
    rg_idx_of_chr[[rgs[i,4]]] = append(rg_idx_of_chr[[rgs[i,4]]], i)
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
        bam_body, ".subset"
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
    write.table(data.frame(bed_chrs, bed_starts, bed_ends), paste0(bam_body, ".subset.segs"), col.names = F, row.names = F, quote = F, sep = "\t")
    seg_coords = paste(bed_starts, bed_ends, sep = "-")

    # Compute average coverage for each segment
    cmd = paste0(
        "bedtools coverage -d -abam ",
        bam_body, ".subset ",
        "-b ", bam_body, ".subset.segs ",
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
file.remove(paste0(bam_body, ".subset"))
file.remove(paste0(bam_body, ".subset.segs"))
