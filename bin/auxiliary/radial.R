#' run radial IVW
#' 
#' @param dat
#' The harmonised dataset, as produced by 
#' \code{\link[TwoSampleMR]{harmonise_data}}.
#' 
#' @return 
#' The results of the IVW analysis, as a list.

get_ivw_radial_results <- function(
    dat
) {
    return(try(ivw_radial_modified(RadialMR::format_radial(
        dat$beta.exposure,
        dat$beta.outcome,
        seBXG = dat$se.exposure,
        seBYG = dat$se.outcome,
        RSID = dat$SNP),
        alpha = 0.05/nrow(dat)
    )))
}

#' remove outliers using radial IVW
#' 
#' @param dat
#' The harmonised dataset, as produced by 
#' \code{\link[TwoSampleMR]{harmonise_data}}.
#' 
#' @return 
#' Returns the harmonised dataset, as produced by 
#' \code{\link[TwoSampleMR]{harmonise_data}}, with outlying variants, as 
#' defined by radial IVW, removed.

remove_outliers_radial <- function(
    dat
) {
    ivw_radial_results <- get_ivw_radial_results(dat)
    
    if (ivw_radial_results$outliers != "No significant outliers") {
        dat <- dat[!(dat$SNP %in% ivw_radial_results$outliers$SNP),]
    }
    
    return(dat)
}

#' MR IVW radial method
#'
#' runs a modified version of \code{\link[RadialMR]{ivw_radial}} that won't fail
#'
#' @references
#' https://github.com/WSpiller/RadialMR/

ivw_radial_modified <- function(
    r_input, 
    alpha, 
    weights, 
    tol
) {
    Ratios <- r_input[, 3] / r_input[, 2]
    F <- r_input[, 2]^2 / r_input[, 4]^2
    mf <- mean(F)
    cat()
    
    if (missing(alpha)) {
        alpha <- 0.05
        warning("Significance threshold for outlier detection not specified: Adopting a 95% threshold")
    }
    
    if (missing(weights)) {
        weights <- 3
        warning("Weights not specified: Adopting modified second-order weights")
    }
    
    if (missing(tol)) {
        tol <- 1e-05
    }
    
    summary <- TRUE
    
    if (weights == 1) {
        W <- ((r_input[, 2]^2) / (r_input[, 5]^2))
    }
    
    if (weights == 2) {
        W <- ((r_input[, 5]^2 / r_input[, 2]^2) + ((r_input[, 3]^2 *
            r_input[, 4]^2) / r_input[, 2]^4))^-1
    }
    
    if (weights == 3) {
        W <- ((r_input[, 2]^2) / (r_input[, 5]^2))
    }
    
    Wj <- sqrt(W)
    BetaWj <- Ratios * Wj
    IVW.Model <- lm(BetaWj ~ -1 + Wj)
    EstimatesIVW <- summary(lm(IVW.Model))
    IVW.Slope <- EstimatesIVW$coefficients[1]
    IVW.SE <- EstimatesIVW$coefficients[2]
    IVW_CI <- confint(IVW.Model)
    DF <- length(r_input[, 1]) - 1
    Qj <- W * (Ratios - IVW.Slope)^2
    Total_Q <- sum(Qj)
    Total_Q_chi <- pchisq(Total_Q, length(r_input[, 2]) - 1,
        lower.tail = FALSE
    )
    
    if (weights == 3) {
        
        W <- ((r_input[, 5]^2 + (IVW.Slope^2 * r_input[, 4]^2)) / r_input[
            ,
            2
        ]^2)^-1
        Wj <- sqrt(W)
        BetaWj <- Ratios * Wj
        IVW.Model <- lm(BetaWj ~ -1 + Wj)
        EstimatesIVW <- summary(lm(BetaWj ~ -1 + Wj))
        IVW.Slope <- EstimatesIVW$coefficients[1]
        IVW.SE <- EstimatesIVW$coefficients[2]
        IVW_CI <- confint(IVW.Model)
        Qj <- W * (Ratios - IVW.Slope)^2
        Total_Q <- sum(Qj)
        Total_Q_chi <- pchisq(Total_Q, length(r_input[, 2]) -
            1, lower.tail = FALSE)
    }
    
    Iterative_ivw <- function(int.tol) {
        
        Diff <- 1
        Bhat1.Iterative <- 0
        count <- 0
        
        while (Diff >= tol) {
            
            W <- 1 / (r_input[, 5]^2 / r_input[, 2]^2 + (Bhat1.Iterative^2) *
                r_input[, 4]^2 / r_input[, 2]^2)
            Wj <- sqrt(W)
            BetaWj <- Ratios * Wj
            new.IVW.Model <- lm(BetaWj ~ -1 + Wj)
            new.EstimatesIVW <- summary(lm(BetaWj ~ -1 + Wj))
            new.IVW.Slope <- new.EstimatesIVW$coefficients[1]
            new.IVW.SE <- new.EstimatesIVW$coefficients[2]
            new.IVW_CI <- confint(new.IVW.Model)
            new.Qj <- W * (Ratios - new.IVW.Slope)^2
            new.Total_Q <- sum(new.Qj)
            new.Total_Q_chi <- pchisq(new.Total_Q, length(r_input[
                ,
                2
            ]) - 1, lower.tail = FALSE)
            Diff <- abs(Bhat1.Iterative - new.IVW.Slope)
            Bhat1.Iterative <- new.IVW.Slope
            Bhat1.SE <- new.IVW.SE
            Bhat1.t <- summary(new.IVW.Model)$coefficients[
                1,
                3
            ]
            Bhat1.p <- summary(new.IVW.Model)$coefficients[
                1,
                4
            ]
            count <- count + 1
        }
        It.Dat <- data.frame(
            Bhat1.Iterative, Bhat1.SE, Bhat1.t,
            Bhat1.p
        )
        
        multi_return2 <- function() {
            Out_list2 <- list(
                It.Res = It.Dat, count = count,
                It.CI = new.IVW_CI
            )
            return(Out_list2)
        }
        
        OUT2 <- multi_return2()
    }
    
    Bhat1.Iterative <- Iterative_ivw(tol)
    
    PL2 <- function(a) {
        b <- a[1]
        w <- 1 / ((phi) * r_input[, 5]^2 / r_input[, 2]^2 + (b^2) *
            r_input[, 4]^2 / r_input[, 2]^2)
        q <- sum(w * (Ratios - b)^2)
    }
    
    PLfunc <- function(a) {
        
        phi <- a[1]
        
        PL2 <- function(a) {
            beta <- a[1]
            w <- 1 / (phi * r_input[, 5]^2 / r_input[, 2]^2 + (beta^2) *
                r_input[, 4]^2 / r_input[, 2]^2)
            q <- (sum(w * (Ratios - beta)^2))
        }
        b <- optimize(PL2, interval = c(lb, ub))$minimum
        w <- 1 / (phi * r_input[, 5]^2 / r_input[, 2]^2 + (b^2) *
            r_input[, 4]^2 / r_input[, 2]^2)
        q <- (sum(w * (Ratios - b)^2) - DF)^2
    }
    
    BootVar <- function(sims = 1000) {
        
        B <- NULL
        pp <- NULL
        
        for (hh in 1:sims) {
            L <- length(r_input[, 2])
            choice <- sample(seq(1, L), L, replace = TRUE)
            while (length(unique(choice)) == 1) {
                choice <- sample(seq(1, L), L, replace = TRUE)
            }
            bxg <- r_input[, 2][choice]
            seX <- r_input[, 4][choice]
            byg <- r_input[, 3][choice]
            seY <- r_input[, 5][choice]
            Ratios <- byg / bxg
            W1 <- 1 / (seY^2 / bxg^2)
            BIVw1 <- Ratios * sqrt(W1)
            sW1 <- sqrt(W1)
            IVWfitR1 <- summary(lm(BIVw1 ~ -1 + sW1))
            phi_IVW1 <- IVWfitR1$sigma^2
            W2 <- 1 / (seY^2 / bxg^2 + (byg^2) * seX^2 / bxg^4)
            BIVw2 <- Ratios * sqrt(W2)
            sW2 <- sqrt(W2)
            IVWfitR2 <- summary(lm(BIVw2 ~ -1 + sW2))
            phi_IVW2 <- IVWfitR2$sigma^2
            phi_IVW2 <- max(1, phi_IVW2)
            phi_IVW1 <- max(1, phi_IVW1)
            lb <- IVWfitR1$coef[1] - 10 * IVWfitR1$coef[2]
            ub <- IVWfitR1$coef[1] + 10 * IVWfitR1$coef[2]
            
            PL2 <- function(a) {
                b <- a[1]
                w <- 1 / ((phi) * seY^2 / bxg^2 + (b^2) * seX^2 / bxg^2)
                q <- sum(w * (Ratios - b)^2)
            }
            
            PLfunc <- function(a) {
                phi <- a[1]
                PL2 <- function(a) {
                    beta <- a[1]
                    w <- 1 / (phi * seY^2 / bxg^2 + (beta^2) * seX^2 / bxg^2)
                    q <- (sum(w * (Ratios - beta)^2))
                }
                b <- optimize(PL2, interval = c(-lb, ub))$minimum
                w <- 1 / (phi * seY^2 / bxg^2 + (b^2) * seX^2 / bxg^2)
                q <- (sum(w * (Ratios - b)^2) - DF)^2
            }
            
            phi <- optimize(PLfunc, interval = c(phi_IVW2, phi_IVW1 +
                0.001))$minimum
            B[hh] <- optimize(PL2, interval = c(lb, ub))$minimum
        }
        
        se <- sd(B)
        mB <- mean(B)
        return(list(mB = mB, se = se))
    }
    
    CIfunc <- function() {
        
        z <- qt(df = DF, 0.975)
        z2 <- 2 * (1 - pnorm(z))
        
        PL3 <- function(a) {
            b <- a[1]
            w <- 1 / (r_input[, 5]^2 / r_input[, 2]^2 + (b^2) * r_input[
                ,
                4
            ]^2 / r_input[, 2]^2)
            q <- (sum(w * (Ratios - b)^2) - qchisq(1 - z2, DF))^2
        }
        lb <- Bhat - 10 * SE
        ub <- Bhat + 10 * SE
        low <- optimize(PL3, interval = c(lb, Bhat))$minimum
        high <- optimize(PL3, interval = c(Bhat, ub))$minimum
        CI <- c(low, high)
        return(list(CI = CI))
    }
    
    phi <- 1
    Bhat <- optimize(PL2, interval = c(-2, 2))$minimum
    W <- 1 / (r_input[, 5]^2 / r_input[, 2]^2 + (Bhat^2) * r_input[
        ,
        4
    ]^2 / r_input[, 2]^2)
    SE <- sqrt(1 / sum(W))
    FCI <- CIfunc()
    QIVW <- sum(W * (Ratios - Bhat)^2)
    Qp <- 1 - pchisq(QIVW, DF)
    Qind <- W * (Ratios - Bhat)^2
    ExactQ <- c(QIVW, Qp)
    ExactQind <- Qind
    FE_EXACT <- t(c(Bhat, SE, Bhat / SE, 2 * (1 - pt(
        abs(Bhat / SE),
        DF
    ))))
    
    FE_EXACT <- data.frame(FE_EXACT)
    names(FE_EXACT) <- c(
        "Estimate", "Std.Error", "t value",
        "Pr(>|t|)"
    )
    
    BIVW1 <- Ratios * sqrt(1 / (r_input[, 5]^2 / r_input[, 2]^2))
    IVWfit1 <- summary(lm(BIVW1 ~ -1 + sqrt(1 / (r_input[, 5]^2 / r_input[
        ,
        2
    ]^2))))
    
    phi_IVW1 <- IVWfit1$sigma^2
    BIVW2 <- Ratios * sqrt(1 / (r_input[, 5]^2 / r_input[, 2]^2 +
        (r_input[, 3]^2) * r_input[, 4]^2 / r_input[, 2]^4))
    IVWfit2 <- summary(lm(BIVW2 ~ -1 + sqrt(1 / (r_input[, 5]^2 / r_input[
        ,
        2
    ]^2 + (r_input[, 3]^2) * r_input[, 4]^2 / r_input[, 2]^4))))
    phi_IVW2 <- IVWfit2$sigma^2
    phi_IVW2 <- max(1, phi_IVW2)
    phi_IVW1 <- max(1, phi_IVW1) + 0.001
    lb <- Bhat - 10 * SE
    ub <- Bhat + 10 * SE

    phi <- optimize(PLfunc, interval = c(phi_IVW2, phi_IVW1))$minimum
    Bhat <- optimize(PL2, interval = c(lb, ub))$minimum

    Boot <- BootVar()
    SE <- Boot$se
    RCI <- Bhat + c(-1, 1) * qt(df = DF, 0.975) * SE
    RE_EXACT <- t(c(Bhat, SE, Bhat / SE, 2 * (1 - pt(
        abs(Bhat / SE),
        DF
    ))))
    
    RE_EXACT <- data.frame(RE_EXACT)
    names(RE_EXACT) <- c(
        "Estimate", "Std.Error", "t value",
        "Pr(>|t|)"
    )
    
    Qj_Chi <- 0
    
    for (i in 1:length(Qj)) {
        Qj_Chi[i] <- pchisq(Qj[i], 1, lower.tail = FALSE)
    }
    
    r_input$Qj <- Qj
    r_input$Qj_Chi <- Qj_Chi
    Out_Indicator <- rep(0, length(r_input[, 2]))
    
    for (i in 1:length(r_input[, 2])) {
        if (Qj_Chi[i] < alpha) {
            Out_Indicator[i] <- 1
        }
    }
    
    r_input$Outliers <- factor(Out_Indicator)
    levels(r_input$Outliers)[levels(r_input$Outliers) == "0"] <- "Variant"
    levels(r_input$Outliers)[levels(r_input$Outliers) == "1"] <- "Outlier"
    
    if (sum(Out_Indicator == 0)) {
        outlier_status <- "No significant outliers"
        outtab <- "No significant outliers"
    }
    
    if (sum(Out_Indicator > 0)) {
        outlier_status <- "Outliers detected"
        Out_Dat <- subset(r_input, Outliers == "Outlier")
        outtab <- data.frame(Out_Dat[, 1], Out_Dat$Qj, Out_Dat$Qj_Chi)
        colnames(outtab) <- c("SNP", "Q_statistic", "p.value")
    }
    
    if (summary == TRUE) {
        cat("\n")
        cat("Radial IVW\n")
        cat("\n")
        Sum.Dat <- data.frame(coef(EstimatesIVW))
        names(Sum.Dat) <- c(
            "Estimate", "Std.Error", "t value",
            "Pr(>|t|)"
        )
        names(Bhat1.Iterative$It.Res) <- names(Sum.Dat)
        combined.dat <- (rbind(Sum.Dat, Bhat1.Iterative$It.Res))
        combined.dat <- rbind(combined.dat, FE_EXACT)
        combined.dat <- rbind(combined.dat, RE_EXACT)
        row.names(combined.dat) <- c(
            "Effect", "Iterative", "Exact (FE)",
            "Exact (RE)"
        )
        if (weights == 1) {
            row.names(combined.dat)[1] <- "Effect (1st)"
        }
        if (weights == 2) {
            row.names(combined.dat)[1] <- "Effect (2nd)"
        }
        if (weights == 3) {
            row.names(combined.dat)[1] <- "Effect (Mod.2nd)"
        }
        print(combined.dat)
        cat("\n")
        cat("\nResidual standard error:", round(
            EstimatesIVW$sigma,
            3
        ), "on", EstimatesIVW$df[2], "degrees of freedom")
        cat("\n")
        cat(paste(c("\nF-statistic:", " on", " and"), round(
            EstimatesIVW$fstatistic,
            2
        ), collapse = ""), "DF, p-value:", format.pval(pf(EstimatesIVW$fstatistic[1L],
            EstimatesIVW$fstatistic[2L], EstimatesIVW$fstatistic[3L],
            lower.tail = FALSE
        ), digits = 3))
        cat("\n")
        cat(
            "Q-Statistic for heterogeneity:", Total_Q, "on",
            length(r_input[, 2]) - 1, "DF", ",", "p-value:",
            Total_Q_chi
        )
        cat("\n")
        cat("\n", outlier_status, "\n")
        cat("Number of iterations =", Bhat1.Iterative$count)
        cat("\n")
    }
    
    out_data <- data.frame(r_input[, 1], r_input[, 6], r_input[
        ,
        7
    ], r_input[, 8])
    out_data$Wj <- Wj
    out_data$BetaWj <- BetaWj
    out_data <- out_data[c(1, 5, 6, 2, 3, 4)]
    names(out_data) <- c(
        "SNP", "Wj", "BetaWj", "Qj", "Qj_Chi",
        "Outliers"
    )
    
    multi_return <- function() {
        Out_list <- list(
            coef = EstimatesIVW$coef, qstatistic = Total_Q,
            df = length(r_input[, 2]) - 1, outliers = outtab,
            data = out_data, confint = confint(IVW.Model), it.coef = combined.dat[2, ],
            fe.coef = combined.dat[3, ], re.coef = combined.dat[4, ], it.confint = Bhat1.Iterative$It.CI, fe.confint = FCI$CI,
            re.confint = RCI, meanF = mf, Total_Q_chi = Total_Q_chi
        )
        class(Out_list) <- "IVW"
        return(Out_list)
    }
    
    OUT <- multi_return()
}
