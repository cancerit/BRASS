artefact_read_pair_count_prob = function(rp_counts, alphas) {
    xmin = 4
    res = sapply(
        alphas,
        function(a) {
            dpldis(rp_counts, xmin, a)
        }
    )
    apply(
        res,
        2,
        function(x) {
            if (any(x == Inf) | any(x <= 0)) {
                rep(-Inf, length(x))
            } else {
                log(x)
            }
        }
    )
}

sample_new_alpha = function(rp_counts, old_alpha, new_proposals = 50, max_alpha = 30, sample_size = 50000) {
    # proposals = runif(new_proposals, 1, max_alpha)
    # probability_matrix = artefact_read_pair_count_prob(rp_counts, proposals)
    # max_log_prob = max(colSums(probability_matrix))
    # sample(
    #     proposals[idx],
    #     1,
    #     prob = exp(colSums(probability_matrix[, idx]) - max_log_prob)
    # )

    # Sample using Metropolis-Hastings-type method instead
    proposal_sd = 0.02
    if (length(rp_counts) > sample_size) {
        rp_counts = sample(rp_counts, sample_size, replace = F)
    }
    old_log_prob = sum(artefact_read_pair_count_prob(rp_counts, old_alpha))
    new_alpha = old_alpha + rnorm(1, 0, proposal_sd)
    while (new_alpha < 1) {
        new_alpha = old_alpha + rnorm(1, 0, proposal_sd)
    }
    new_log_prob = sum(artefact_read_pair_count_prob(rp_counts, new_alpha))
    new_log_prob = new_log_prob + pnorm(0, old_alpha, proposal_sd, lower.tail = F, log = T) - pnorm(0, new_alpha, proposal_sd, lower.tail = F, log = T)
    if (new_log_prob > old_log_prob) {
        new_alpha
    } else {
        prob = exp(new_log_prob - old_log_prob)
        sample(c(new_alpha, old_alpha), 1, prob = c(prob, 1 - prob))
    }
}

#
# Use gamma instead of sn
#
sample_new_gamma_params = function(isize, artefact_categories, gamma_params, proposal_sd = 0.1, sample_size = 50000) {
    if (length(isize) > sample_size) {
        idx = sample(1:length(isize), sample_size, replace = F)
        isize = isize[idx]
        artefact_categories = artefact_categories[idx]
    }

    old_mean = gamma_params[1]
    old_shape_1 = gamma_params[2]
    old_shape_2 = gamma_params[3]
    old_rate_1 = old_shape_1 / old_mean
    old_rate_2 = old_shape_2 / old_mean

    old_log_prob_1 = sum(dgamma(log10(isize[artefact_categories == 1]), shape = old_shape_1, rate = old_rate_1, log = T))
    old_log_prob_2 = sum(dgamma(log10(isize[artefact_categories == 2]), shape = old_shape_2, rate = old_rate_2, log = T))
    old_log_prob = old_log_prob_1 + old_log_prob_2

    new_params = c(
        mean = exp(rnorm(1, log(old_mean), proposal_sd)),
        shape_1 = exp(rnorm(1, log(old_shape_1), proposal_sd)),
        shape_2 = exp(rnorm(1, log(old_shape_2), proposal_sd))
    )
    new_rate_1 = new_params["shape_1"] / new_params["mean"]
    new_rate_2 = new_params["shape_2"] / new_params["mean"]

    new_log_prob_1 = sum(dgamma(log10(isize[artefact_categories == 1]), shape = new_params["shape_1"], rate = new_rate_1, log = T))
    new_log_prob_2 = sum(dgamma(log10(isize[artefact_categories == 2]), shape = new_params["shape_2"], rate = new_rate_2, log = T))
    new_log_prob = new_log_prob_1 + new_log_prob_2
    if (new_log_prob > old_log_prob) {
        new_params
    } else if (runif(1) < exp(new_log_prob - old_log_prob)) {
        new_params
    } else {
        c(mean = old_mean, shape_1 = old_shape_1, shape_2 = old_shape_2)
    }
}

# Class assignment function
is_artefact_class = function(size, rp_count, alpha, gamma_params, artefacts_prior_probs) {
    mean = gamma_params[1]
    shape_1 = gamma_params[2]
    shape_2 = gamma_params[3]

    # Compute likelihood for non-artefact
    non_artefact_prob = artefacts_prior_probs[1] / log10(1e8) / max(rp_count)

    # Likelihood for first and second artefacts class
    rp_count_prob = exp(artefact_read_pair_count_prob(rp_count, alpha))[,1]

    size_prob = dgamma(log10(size), shape = shape_1, rate = shape_1 / mean)
    combined_artefact_prob_1 = size_prob * rp_count_prob * artefacts_prior_probs[2]

    size_prob = dgamma(log10(size), shape = shape_2, rate = shape_2 / mean)
    combined_artefact_prob_2 = size_prob * rp_count_prob * artefacts_prior_probs[3]

    apply(
        matrix(c(rep(non_artefact_prob, length(size)), combined_artefact_prob_1, combined_artefact_prob_2), ncol = 3),
        1,
        function(probs) {
            sample(0:2, 1, prob = probs)
        }
    )
}

run_mcmc_for_artefacts = function(rp_count, size, prev_alpha, prev_gamma_params, prev_artefact_category, n_iter, burn_in) {
    sampled_alpha = c()
    sampled_gamma_params = matrix(NA, nrow = 0, ncol = 3)
    sampled_artefacts_1_sum = rep(0, length(rp_count))
    sampled_artefacts_2_sum = rep(0, length(rp_count))

    cat("Iterating...\n", file = stderr())
    for (i in 1:n_iter) {
        if (i %% 100 == 0) {
            cat(i, "\n", sep = "", file = stderr())
            cat(prev_gamma_params, sum(prev_artefact_category == 0), sum(prev_artefact_category == 1), sum(prev_artefact_category == 2), "\n", file = stderr())
        }
        alpha = sample_new_alpha(rp_count[prev_artefact_category > 0], prev_alpha)
        gamma_params = sample_new_gamma_params(size, prev_artefact_category, prev_gamma_params)
        artefact_category = is_artefact_class(
            size,
            rp_count,
            prev_alpha,
            prev_gamma_params,
            c(
                sum(prev_artefact_category == 0) + 1,
                sum(prev_artefact_category == 1) + 1,
                sum(prev_artefact_category == 2) + 1
            ) / (length(size) + 3)
        )

        if (i > burn_in) {
            sampled_alpha = append(sampled_alpha, alpha)
            sampled_gamma_params = rbind(sampled_gamma_params, gamma_params)
            sampled_artefacts_1_sum = sampled_artefacts_1_sum + (artefact_category == 1)
            sampled_artefacts_2_sum = sampled_artefacts_2_sum + (artefact_category == 2)
        }

        prev_alpha = alpha
        prev_gamma_params = gamma_params
        prev_artefact_category = artefact_category
    }

    list(
        alpha = mean(sampled_alpha),
        gamma_params = colMeans(sampled_gamma_params),
        artefact_prob = data.frame(
            true = 1 - (sampled_artefacts_1_sum + sampled_artefacts_2_sum) / (n_iter - burn_in),
            artefact_1 = sampled_artefacts_1_sum / (n_iter - burn_in),
            artefact_2 = sampled_artefacts_2_sum / (n_iter - burn_in)
        )
    )
}


# Read in files
f = commandArgs(T)[1]
s = sub("_vs_.+", "", f)
s = sub(".+/", "", s)
d = read.table(f, header = F, sep = "\t", stringsAsFactors = F, comment = "")
idx = d[,1] == d[,4] & d[,9] == d[,10]
d = d[idx,]
rp_count = d[,8]
size = rowMeans(d[,5:6]) - rowMeans(d[,2:3])

# Initiate parameters
library(MASS)
library(poweRlaw)
library(RColorBrewer)
idx = size <= 1e5

# Make decision on what to plot
pdf_file = sub("\\..+", ".inversions.pdf", f)
pdf(pdf_file, h = 10)
layout(matrix(c(4, 4, 1, 2, 4, 4, 3, 3), ncol = 2))
bad_groups_count = sum(size <= 1e5 & rp_count == 4)
if (sum(idx) == 0) {
    plot(c(), axes = F, xlab = "", ylab = "", main = "", xlim = 0:1, ylim = 0:1)
    plot(c(), axes = F, xlab = "", ylab = "", main = "", xlim = 0:1, ylim = 0:1)
    plot(c(), axes = F, xlab = "", ylab = "", main = "", xlim = 0:1, ylim = 0:1)
} else {
    y = table(rp_count[idx])
    plot(y, type = "h", xlab = "Read pair count", ylab = "Frequency", col = "indianred", xlim = c(1, 50))

    # Perform power law modelling?
    title(main = paste0("bad groups count: ", bad_groups_count))
    if (bad_groups_count >= 50) {
        # The MCMC
        m = displ$new(rp_count[idx])
        m$setXmin(4)
        alpha = estimate_pars(m)$pars
        gamma_params = coef(fitdistr(log10(size[idx]), "gamma"))  # Gives shape and rate
        gamma_params = c(gamma_params[1] / gamma_params[2], gamma_params[1])  # Change into mean and shape_1
        gamma_params = c(gamma_params, gamma_params[2] / 10)  # Add initial shape_2
        mcmc_res = run_mcmc_for_artefacts(
            rp_count,
            size,
            alpha,
            gamma_params,
            ifelse(
                idx,
                1,
                0
            ),
            5000,
            1000
        )

        alpha = mcmc_res[["alpha"]]
        lines(
            4:50,
            dpldis(4:50, 4, alpha) * sum(idx),
            col = "blue",
            lwd = 0.5
        )
        par(usr = c(0, 1, 0, 1))
        text(0.8, 0.8, paste0("alpha: ", alpha), pos = 1)

        plot(
            4:max(as.numeric(names(y))),
            ppldis(4:(max(as.numeric(names(y)))), 4, alpha),
            type = "l",
            xlab = "Read pair count",
            ylab = "Cumulative density"
        )
        temp = rp_count[idx]
        lines(
            as.numeric(names(y)),
            cumsum(ifelse(
                as.numeric(names(y)) < 4,
                0,
                y
            )) / sum(y),
            col = "blue"
        )
        legend("bottomright", bty = "n", col = c("black", "blue"), lty = 1, legend = c("Predicted", "Observed"))
    }
    else {
        plot(c(), axes = F, xlab = "", ylab = "", main = "", xlim = 0:1, ylim = 0:1)
    }

    x = seq(0, 5, 0.02)
    hist(log10(size[idx]), ylab = "", col = "indianred", border = NA, breaks = x, main = "Insert size distribution", xlab = "log10(insert size)")

    if (bad_groups_count >= 50) {
        gamma_params = mcmc_res[["gamma_params"]]
        shape_1 = gamma_params[2]
        rate_1 = gamma_params[2] / gamma_params[1]
        gamma1_prob = sum(mcmc_res[["artefact_prob"]][,2]) / sum(idx)
        temp = gamma1_prob * pgamma(x, shape = shape_1, rate = rate_1)
        lines(x[-1] - (x[2]-x[1])/2, sum(idx) * (temp[-1] - temp[-length(temp)]), col = "green")

        gamma2_prob = sum(mcmc_res[["artefact_prob"]][,3]) / sum(idx)
        shape_2 = gamma_params[3]
        rate_2 = gamma_params[3] / gamma_params[1]
        temp = gamma2_prob * pgamma(x, shape = shape_2, rate = rate_2)
        lines(x[-1] - (x[2]-x[1])/2, sum(idx) * (temp[-1] - temp[-length(temp)]), col = "blue")

        par(usr = c(0, 1, 0, 1))
        text(
            0.1,
            0.9,
            adj = c(0, 1),
            bty = "n",
            labels = paste0(
                "mean = ", round(gamma_params[1], 4), "\n",
                "shape_1 = ", round(gamma_params[2], 4), "\n",
                "shape_2 = ", round(gamma_params[3], 4), "\n",
                "count1 = ", round(sum(mcmc_res[["artefact_prob"]][,2])), "\n",
                "count2 = ", round(sum(mcmc_res[["artefact_prob"]][,3]))
            )
        )
    }
}

cols = brewer.pal(11, "RdYlBu")
if (sum(idx) > 0 && bad_groups_count >= 50) {
    col = cols[findInterval(mcmc_res[["artefact_prob"]][,1], seq(0, 1, length.out = 12), all.inside = T)]
    threshold = 5
    threshold_idx = which.max(cumsum(sort(1 - mcmc_res[["artefact_prob"]][,1])) > threshold)
    idx = rank(1 - mcmc_res[["artefact_prob"]][,1]) < threshold_idx
} else {
    col = cols[1]
    threshold_idx = length(size)
    idx = rep(T, length(size))
}

if(length(rp_count) > 0) {
  x = jitter(rp_count)
  plot(
      x,
      size / 1e3,
      xlab = "Read pair count",
      ylab = "Inner size (kb)",
      main = paste(s, " inversion calls\n", threshold_idx, "/", length(size), " kept", sep = ""),
      log = "xy",
      pch = 16,
      col = col,
      ylim = c(0.01, 1000000)
  )
  points(x[idx], size[idx]/1e3, pch = ".", col = "white")
  abline(h = 10^(-2:6), col = "grey", lty = 3)
  abline(h = c(300, 500)/1e3, col = "grey", lty = 3)
}
dev.off()


# Write out results
output_file = sub(".inversions.pdf", ".is_fb_artefact.txt", pdf_file)
if (bad_groups_count >= 50) {
    write.table(
        data.frame(
            d[,7],  # ID
            mcmc_res[["artefact_prob"]][,1],  # Probability of being true
            rank(1 - mcmc_res[["artefact_prob"]][,1]) < threshold_idx  # Whether the rearrangement is to be kept
        ),
        output_file,
        row.names = F,
        col.names = F,
        sep = "\t",
        quote = F
    )
} else {
    write.table(
        data.frame(
            d[,7],  # ID
            rep(1, nrow(d)),
            rep(TRUE, nrow(d))
        ),
        output_file,
        row.names = F,
        col.names = F,
        sep = "\t",
        quote = F
    )
}

