## sliding-window analysis of state transition rates through time (migration)
##
## inputs (hopefully you have these):
## - tree: ultrametric, time-scaled phylo (branch lengths in Myr)
## - dat: data.frame(species, migration) as seen in corHMM()
## - fit: corHMM() fit used for marginal reconstruction
##
## workflow:
## 1. draw stochastic character maps from the fitted corHMM model, each map
##    gives transition times across branches
## 2. for each time window, count transitions and divide by the time
##    lineages spent in the window, giving a per-lineage rate, per-direction
##    rates also calculated
## 3. repeat step 2 on histories simulated under the constant-rate model, 
##    windows above null have elevated rates
## 4. (to do) likelihood tests, fit model with a separate Q
##    inside window vs. outside, compare with constant-Q model


library(ape)
library(phytools)
library(corHMM)


## Helpers

## convert corHMM solution (NA diagonal, may contain NA for disallowed
## transitions) into proper Q matrix (rows = from, cols = to, rows sum to 0)
## added fix for makeSimmap parsing states differently (0,1,2 vs 1,2,3)
corhmm_Q <- function(fit, dat = NULL) {
  Q <- fit$solution
  Q[is.na(Q)] <- 0
  diag(Q) <- -rowSums(Q)
  ## "(1,R1)" -> "1"
  rownames(Q) <- colnames(Q) <- sub("^\\((.*),R1\\)$", "\\1", rownames(Q))
  if (!is.null(dat)) {
    leg <- getStateMat4Dat(dat)$legend # names = corHMM index, values = data codes
    rownames(Q) <- colnames(Q) <- unname(leg[rownames(Q)])
  }
  Q
}

## flatten simmap tree into state segments and transition events,
## with time expressed as age (Myr before present, more = older)
map_events <- function(st) {
  H     <- nodeHeights(st)
  T0    <- max(H)
  n     <- lengths(st$maps)
  e     <- rep(seq_along(n), n) # edge index of each segment
  dur   <- unlist(st$maps, use.names = FALSE)
  state <- unlist(lapply(st$maps, names), use.names = FALSE)
  first <- !duplicated(e) # first segment on each edge
  cs    <- cumsum(dur)
  off   <- (cs - dur)[first] # cumu. sum offset per edge
  end_h <- H[e, 1] + cs - off[e] # height at segment end
  seg   <- data.frame(state = state,
                      old   = T0 - (end_h - dur), # age at segment start
                      young = T0 - end_h) # age at segment end
  k  <- which(!first) # change happens at start of every non-first segment
  ev <- data.frame(from = state[k - 1], to = state[k], age = seg$old[k])
  list(seg = seg, ev = ev, T0 = T0)
}

## window stats for one flattened map
## lo/hi are window age bounds (lo = younger edge)
window_stats <- function(me, lo, hi, states) {
  ## time each segment spends inside each window (segments x windows)
  ov <- pmax(outer(me$seg$old, hi, pmin) - outer(me$seg$young, lo, pmax), 0)
  Ti <- t(matrix(vapply(states, function(s) colSums(ov[me$seg$state == s, , drop = FALSE]),
                        numeric(length(lo))), nrow = length(lo))) # states x windows
  rownames(Ti) <- states
  L  <- colSums(Ti) # total lineage-Myr per window

  ## events inside each window by direction
  inw <- outer(me$ev$age, lo, ">=") & outer(me$ev$age, hi, "<")
  pairs <- expand.grid(from = states, to = states, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$from != pairs$to, ]
  Nij <- t(mapply(function(f, t) colSums(inw[me$ev$from == f & me$ev$to == t, , drop = FALSE]),
                  pairs$from, pairs$to))
  if (length(lo) == 1) Nij <- matrix(Nij, ncol = 1)
  rownames(Nij) <- paste0(pairs$from, "->", pairs$to)

  qij <- Nij / Ti[pairs$from, , drop = FALSE] # per-direction rate (per Myr in state i)
  rownames(qij) <- paste0("q", pairs$from, pairs$to)

  data.frame(lo = lo, hi = hi, mid = (lo + hi) / 2,
             lineage_time = L,
             n_trans = colSums(Nij),
             rate = colSums(Nij) / L, # transitions per lineage-Myr
             t(Nij), t(qij), t(Ti / rep(L, each = length(states))),
             check.names = FALSE)
}

## run window_stats over a list of simmap trees
## tree-wide value of every statistic (one window covering the whole tree)
## is stored per map in attr(, "overall"), used in high_rate_windows()
windows_over_maps <- function(simmaps, width, step, max_age = NULL, states) {
  flat <- lapply(simmaps, map_events)
  if (is.null(max_age)) max_age <- flat[[1]]$T0
  lo <- seq(0, max_age - width, by = step)
  hi <- lo + width
  out <- do.call(rbind, lapply(seq_along(flat), function(i)
    cbind(rep = i, window_stats(flat[[i]], lo, hi, states))))
  attr(out, "overall") <- do.call(rbind, lapply(seq_along(flat), function(i)
    cbind(rep = i, window_stats(flat[[i]], 0, flat[[i]]$T0 * (1 + 1e-9), states))))
  out
}

## summarize replicate maps per window, median and 95% interval of columns
summarise_windows <- function(df, col) {
  s <- aggregate(df[[col]], by = list(mid = df$mid),
                 FUN = function(x) quantile(x, c(.025, .5, .975), na.rm = TRUE))
  data.frame(mid = s$mid, lwr = s$x[, 1], med = s$x[, 2], upr = s$x[, 3])
}


## stochastic maps from the fitted model (observed) and null histories
## use_null = FALSE skips the constant-rate simulations and you get only the
## observed rates, and high_rate_windows() picks out the fast windows
run_sliding_window <- function(tree, dat, fit, width, step,
                               nsim = 100, max_age = NULL, seed = 1,
                               use_null = TRUE) {
  set.seed(seed)
  Q      <- corhmm_Q(fit)
  states <- rownames(Q)

  ## observed: maps conditioned on the tip data under the fitted model
  obs_maps <- makeSimmap(tree = tree, data = dat, model = fit$solution,
                         rate.cat = 1, root.p = fit$root.p, nSim = nsim)

  obs  <- windows_over_maps(obs_maps, width, step, max_age, states)
  rm(obs_maps)

  null <- NULL
  if (use_null) {
    ## null: unconditional histories under the same constant-rate Q
    ## root state is drawn from marginal root probabilities
    root_p <- fit$states[1, ] # row 1 = root node
    names(root_p) <- states
    null_maps <- lapply(seq_len(nsim), function(i)
      sim.history(tree, Q, anc = sample(states, 1, prob = root_p),
                  nsim = 1, direction = "row_to_column", message = FALSE))
    null <- windows_over_maps(null_maps, width, step, max_age, states)
  }
  list(obs = obs, null = null, Q = Q, states = states,
       width = width, step = step)
}

## qualitative screen for fast windows (no null model/likelihood testing)
## on each map the window value of `col` is divided by the same statistic
## for the whole tree, e.g. rate / (all transitions / total tree length),
## or q23 / tree-wide q23. these "base rates" lead to:
##
## - rel_med: median relative rate across maps (1 = tree-wide average)
## - prop_above: fraction of maps where the relative rate >= ratio,
##                measures robustness to mapping uncertainty
## - n_trans_med: median number of transitions in the window
## - rank: 1 = fastest window (among windows with >= min_trans transitions)
##
## windows are considered "high" when rel_med >= ratio + prop_above >= support,
## ones with few transitions (n_trans_med < min_trans) are never called
## high because rates there are mostly noise (typically near the root)
high_rate_windows <- function(res, col = "rate", ratio = 1.5, support = 0.5,
                              min_trans = 5) {
  o   <- res$obs
  ov  <- attr(o, "overall")
  ref <- ov[[col]][match(o$rep, ov$rep)]
  o$rel <- o[[col]] / ref
  agg <- function(x, f) tapply(x, o$mid, f)
  out <- data.frame(
    mid         = sort(unique(o$mid)),
    lo          = agg(o$lo, min),
    hi          = agg(o$hi, max),
    value_med   = agg(o[[col]], function(x) median(x, na.rm = TRUE)),
    tree_wide   = median(ov[[col]], na.rm = TRUE),
    rel_med     = agg(o$rel, function(x) median(x, na.rm = TRUE)),
    rel_lwr     = agg(o$rel, function(x) quantile(x, .025, na.rm = TRUE)),
    rel_upr     = agg(o$rel, function(x) quantile(x, .975, na.rm = TRUE)),
    prop_above  = agg(o$rel, function(x) mean(x >= ratio, na.rm = TRUE)),
    n_trans_med = agg(o$n_trans, median),
    lineage_time = agg(o$lineage_time, median),
    row.names = NULL)
  ok <- !is.na(out$rel_med) & out$n_trans_med >= min_trans
  out$rank <- NA_integer_
  out$rank[ok] <- rank(-out$rel_med[ok], ties.method = "min")
  out$high <- ok & out$rel_med >= ratio & out$prop_above >= support
  out[order(-out$mid), ] # oldest first
}

## flag windows where observed median rate exceeds null 97.5% quantile
flag_elevated <- function(res, col = "rate") {
  if (is.null(res$null))
    stop("no null simulations (use_null = FALSE); use high_rate_windows() instead")
  o <- summarise_windows(res$obs, col)
  n <- summarise_windows(res$null, col)
  data.frame(mid = o$mid, obs_med = o$med, obs_lwr = o$lwr, obs_upr = o$upr,
             null_med = n$med, null_upr = n$upr,
             elevated = o$med > n$upr)
}



