library(data.table)

poly <- fread("analyses/pub/cardiometabolic_proteins/review2/leaveout_assocs_with_aptamer.txt")
poly <- poly[LD_blocks_removed != 0]

rm <- poly[P > 0.05, .SD[which.min(LD_blocks_removed)], by=.(PRS, Gene)]
poly <- poly[rm[,.(PRS, Gene, LD_blocks_removed)], on = .(PRS, Gene, LD_blocks_removed <= LD_blocks_removed),
             .(PRS, Gene, x.LD_blocks_removed, BP_removed, pct_removed, Beta, L95, U95, P,
               Aptamer, Beta.aptamer, L95.aptamer, U95.aptamer, P.aptamer)]
poly[Gene == "C5orf38", Gene := "CEI"]
poly[Gene == "PDE4D", Gene := "PDE4D/A"]
setnames(poly, c("Gene", "x.LD_blocks_removed"), c("Protein", "LD_blocks_removed"))

fwrite(poly, sep="\t", quote=FALSE, file="analyses/pub/cardiometabolic_proteins/review3/leaveout_chunk_assocs.txt")
