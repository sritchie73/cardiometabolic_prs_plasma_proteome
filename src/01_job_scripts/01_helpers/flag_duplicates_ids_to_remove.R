library(data.table)

args <- commandArgs(trailingOnly = TRUE)
info_file <- args[1]
out_file <- args[2]

info <- fread(info_file, header=FALSE)
setnames(info, c("id", "chr", "pos", "A1", "A2", "info"))

# First handle cases where variant has the same chr, pos, and alleles:
info[, row := .I]
info[, sorted_alleles := paste(sort(c(A1, A2)), collapse=":"), by=row]
info[, real_id := paste(chr, pos, sorted_alleles, sep=":")]

max_info <- info[,.(info=max(info)), by=real_id]
to_keep <- info[max_info, on = .(real_id, info), mult="first"]
to_remove <- info[!to_keep, on = .(id)]

# Also handle cases where the rsID is a duplicate, but has a slightly different position
to_keep[, rsid := gsub(":[0-9]*$", "", id)]
to_keep[, rsid_alleles := paste(rsid, sorted_alleles, sep=":")]
dup_bad_pos <- to_keep[rsid != ".", .N, by=rsid_alleles][N > 1, .(rsid_alleles)]
dup_bad_pos <- to_keep[dup_bad_pos, on = .(rsid_alleles)]
dup_to_keep <- dup_bad_pos[, .SD[which.max(info)], by=rsid_alleles]
dup_to_remove <- dup_bad_pos[!dup_to_keep, on = .(id)]

to_remove <- rbind(to_remove[, .(id)], dup_to_remove[,.(id)])
fwrite(to_remove[,.(id)], file=out_file, quote=FALSE, col.names=FALSE)

