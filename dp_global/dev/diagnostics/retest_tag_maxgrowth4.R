library(data.table)
source("dp_global/dev/diagnostics/utils.R")
compile_transition_cost()
Dt_all <- load_dataset("data_simulation/data/simulated_data_1.csv")
Dt2747 <- Dt_all[Tag == 2747]
# set hard max growth to 4.0 (blocks 2.7->22.2)
Dt2747[, Bio_Max_Growth := 4]
Dt2747[, Bio_Max_Growth_Soft := 0]
Dt2747[, Bio_K_Growth := 0]
# ensure other Bio_* present (copied from previous script)
Dt2747[, Bio_Mu_Growth := 0.350960664887803]
Dt2747[, Bio_Gamma_Growth := 0.314805038516982]
Dt2747[, Bio_Sigma0_Growth := 0.269453148265253]
Dt2747[, Bio_Sigma1_Growth := 0.0105915836625149]
Dt2747[, Bio_H0_Mortality := 0.00873479569505417]
Dt2747[, Bio_Beta_Mortality := 0.00834401015479029]
Dt2747[, Bio_Recruit_Meanlog := -0.250726454828908]
Dt2747[, Bio_Recruit_Sdlog := 0.838668739452498]
Dt2747[, Bio_Recruit_MaxDBH_unit := 38.4999]
Dt2747[, Bio_Recruitment_lambda := 0.0552806029785688]
Dt2747[, Bio_Max_Shrink := -0.5]
Dt2747[, Bio_K_Shrink := 0]

res <- run_dp_on_dt(Dt2747, anchor_start = 7, min_growth = -2.5, max_growth = 7.5, posterior_top_k = 1, verbose = TRUE)
cat("\nReconstruction with Bio_Max_Growth=4.0 (hard cap)\n")
print(res[, .(CensusID, DBH, OriginalStemID, ReconstructedStemID, ReconstructionMethod, ConstraintViolation)])

with(res, {
    stems <- unique(ReconstructedStemID)
    cols <- rainbow(length(stems))
    plot(NULL,
        xlim = c(0.5, 4.5),
        ylim = c(0, 25),
        main = "Tag 2747 DP Reconstruction with Bio_Max_Growth=4.0",
        xlab = "CensusID",
        ylab = "DBH"
    )
    for (i in seq_along(stems)) {
        idx <- ReconstructedStemID == stems[i]
        lines(CensusID[idx], DBH[idx],
            type = "b",
            pch = 16,
            col = cols[i]
        )
    }
    legend("topleft", legend = paste("Stem", stems), col = cols, pch = 16, lty = 1)
})

