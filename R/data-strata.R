#' Combine levels of a stratum variable
#'
#' Adds an aggregation layer by combining levels of one stratum variable
#' into groups: causes of death into cause groups, deprivation deciles into
#' quintiles, districts into regions, single nationalities into larger
#' groups. Death counts are summed. Exposures and populations are summed
#' when the variable divides the population (regions, deciles,
#' `additive = "deaths"`) and kept when every level shares the same
#' population (causes of death, `additive = "rates"`). A group with one
#' member renames the level.
#'
#' A residual group ("all other causes") is defined in one of two ways.
#' With `other` alone it is the sum of the levels not listed in `groups`.
#' With `other` and `total` it is the total minus all other groups, so it
#' also holds deaths that no level covers (for example when the source lists
#' selected causes and all causes). Negative remainders (the listed levels
#' exceed the total) are set to zero with a warning.
#'
#' To add totals and groups as extra strata while keeping the original
#' levels, use [mt_aggregate_strata()] instead.
#'
#' @param data An [mt_data] object.
#' @param var Stratum variable whose levels are combined.
#' @param groups Named list (the levels of `var` in each group), or a
#'   function that maps a vector of levels to group names (for example
#'   `function(d) paste0("Q", ceiling(as.numeric(d) / 2))` for deciles to
#'   quintiles).
#' @param other Name of a residual group, or `NULL` to keep unlisted levels
#'   as they are.
#' @param total Level of `var` that is the total over all levels (for
#'   example all causes), or `NULL`. It is kept, and with `other` the
#'   residual group is the total minus the other groups.
#' @param keep Levels kept as they are when `other` is given.
#' @param additive `"deaths"` (exposures are summed) or `"rates"`
#'   (exposures are shared by the levels and kept).
#' @param name New name of the stratum variable (for example
#'   `"imd_quintile"`), or `NULL` to keep `var`.
#' @return The `mt_data` object with the combined levels. Cells are summed
#'   where every member level is observed; others are dropped with a
#'   warning.
#' @examples
#' x <- mt_example_causes(years = 2015:2019)
#' y <- mt_group_levels(x, "cause", list(noncancer = c("cvd", "other")), additive = "rates")
#' mt_strata(y)
#' # deciles to quintiles, with a new name for the aggregation layer
#' d <- expand.grid(year = 2020, age = "90+", decile = 1:10)
#' d$deaths <- 10 + d$decile
#' d$exposure <- 1000
#' z <- mt_data(deaths = d, exposure = d, by = "decile")
#' z <- mt_group_levels(z, "decile", function(l) paste0("Q", ceiling(as.numeric(l) / 2)), name = "quintile")
#' mt_strata(z)
#' @export
mt_group_levels <- function(data, var, groups, other = NULL, total = NULL, keep = character(),
                            additive = c("deaths", "rates"), name = NULL) {
  check_class(data, "mt_data", "data")
  additive <- match.arg(additive)
  if (!is.character(var) || length(var) != 1 || !var %in% data$by) {
    mt_abort(sprintf("`var` must be one of the stratum variables (%s).", paste(data$by, collapse = ", ")),
             class = "argument")
  }
  levels_all <- unique(data$strata[[var]])
  if (is.function(groups)) {
    g <- as.character(groups(levels_all))
    if (length(g) != length(levels_all)) mt_abort("`groups` must return one group per level.", class = "argument")
    groups <- split(levels_all[!is.na(g)], factor(g[!is.na(g)], levels = unique(g[!is.na(g)])))
  }
  if (!is.list(groups) || (length(groups) && (is.null(names(groups)) || any(!nzchar(names(groups)))))) {
    mt_abort("`groups` must be a named list of levels or a function.", class = "argument")
  }
  groups <- lapply(groups, as.character)
  listed <- unlist(groups, use.names = FALSE)
  unknown <- setdiff(c(listed, total, keep), levels_all)
  if (length(unknown)) {
    mt_abort(sprintf("Unknown level(s) of `%s`: %s.", var, paste(unknown, collapse = ", ")), class = "argument")
  }
  if (anyDuplicated(listed)) mt_abort("A level is in more than one group.", class = "argument")
  if (!is.null(total) && total %in% listed) mt_abort("The `total` level cannot be in a group.", class = "argument")
  unlisted <- setdiff(levels_all, c(listed, total))
  if (!is.null(other)) {
    if (other %in% c(names(groups), keep, total)) mt_abort("`other` must be a new group name.", class = "argument")
    rest <- setdiff(unlisted, keep)
    unlisted <- intersect(unlisted, keep)
    if (is.null(total) && length(rest)) groups[[other]] <- rest
  }
  for (l in unlisted) groups[[l]] <- l
  if (!is.null(total)) groups[[total]] <- total
  out <- regroup_strata(data, var, groups, additive)
  if (!is.null(total) && !is.null(other)) {
    out <- add_remainder(out, var, total, setdiff(names(groups), total), other, additive)
  }
  if (!is.null(name) && !identical(name, var)) {
    if (name %in% data$by) mt_abort(sprintf("`%s` is already a stratum variable.", name), class = "argument")
    setnames(out$strata, var, name)
    out$by[out$by == var] <- name
  }
  out
}

# Residual level: the total minus the other levels of `var`, per stratum of
# the other variables and cell
add_remainder <- function(x, var, total, parts, other, additive) {
  st <- x$strata
  others <- setdiff(x$by, var)
  tot <- st[st[[var]] == total]
  new <- copy(tot)
  new[, (var) := other]
  new[, stratum := max(st$stratum) + seq_len(.N)]
  new[, label := do.call(paste, c(lapply(x$by, function(b) new[[b]]), sep = "."))]
  map <- data.table(tot = tot$stratum, new = new$stratum)
  cells <- c("age_lo", "age_hi", "period_lo", "period_hi")
  st_key <- st[, c("stratum", others, var), with = FALSE]
  negative <- 0L
  for (m in names(x$obs)) {
    o <- x$obs[[m]]
    if (m != "deaths" && additive == "rates") {
      r <- o[stratum %in% map$tot][, stratum := map$new[match(stratum, map$tot)]]
    } else {
      o2 <- merge(o, st_key, by = "stratum")
      grp <- c(others, cells)
      tt <- o2[o2[[var]] == total, c("stratum", grp, "value"), with = FALSE]
      pp <- o2[o2[[var]] %in% parts, list(s = sum(value), n = .N), by = grp]
      r <- merge(tt, pp, by = grp)
      r <- r[n == length(parts)]
      r[, value := value - s]
      negative <- negative + sum(r$value < -1e-8)
      r[, value := pmax(value, 0)]
      r <- r[, c("stratum", cells, "value"), with = FALSE]
      r[, stratum := map$new[match(stratum, map$tot)]]
    }
    r[, provenance := "aggregated"]
    setcolorder(r, c("stratum", cells, "value", "provenance"))
    x$obs[[m]] <- rbind(o, r)
    setorderv(x$obs[[m]], c("stratum", "period_lo", "age_lo"))
  }
  if (negative) {
    mt_warn(sprintf("The other levels exceed the total in %d cell(s); the remainder `%s` is set to zero there.",
                    negative, other), class = "group")
  }
  x$strata <- rbind(st, new)
  x
}

#' Rename levels of a stratum variable
#'
#' Gives levels of a stratum variable new labels, for example sex codes
#' `"1"` and `"2"` as `"male"` and `"female"`. Levels mapped to the same new
#' label are combined (see [mt_group_levels()]).
#'
#' @param data An [mt_data] object.
#' @param var Stratum variable.
#' @param labels Named character vector, `c(old = "new")`.
#' @inheritParams mt_group_levels
#' @return The `mt_data` object with the new labels.
#' @examples
#' d <- data.frame(Sex = rep(1:2, each = 2), year = 2020:2021, age = "90+", deaths = 5:8, exposure = 100)
#' x <- mt_data(deaths = d, exposure = d, by = "Sex")
#' mt_strata(mt_relabel(x, "Sex", c("1" = "male", "2" = "female")))
#' @export
mt_relabel <- function(data, var, labels, additive = c("deaths", "rates")) {
  check_class(data, "mt_data", "data")
  if (!is.character(labels) || is.null(names(labels))) {
    mt_abort("`labels` must be a named character vector, c(old = \"new\").", class = "argument")
  }
  groups <- split(names(labels), factor(unname(labels), levels = unique(labels)))
  mt_group_levels(data, var, groups, additive = match.arg(additive))
}

#' Restrict data to ages, periods and strata
#'
#' @param data An [mt_data] object.
#' @param ... Conditions on the stratum variables, for example
#'   `region != "Germany"` or `sex == "female"`.
#' @param ages Age range `c(lower, upper)`: cells with lower bound at or
#'   above `lower` and upper bound at or below `upper` (`Inf` for an open
#'   interval) are kept.
#' @param periods Period range `c(first, last)` (calendar years, inclusive;
#'   `NA` for an open end, for example `c(NA, 2019)`).
#' @param drop Remove stratum variables that have a single level left (for
#'   example `cause` after `cause == "all"`).
#' @return The restricted `mt_data` object; strata are renumbered and the
#'   grid, if set, is restricted to the same ranges.
#' @examples
#' x <- mt_example(years = 2001:2019)
#' y <- mt_subset(x, sex == "female", ages = c(60, Inf), periods = c(2010, 2019))
#' y
#' @export
mt_subset <- function(data, ..., ages = NULL, periods = NULL, drop = FALSE) {
  check_class(data, "mt_data", "data")
  cond <- substitute(list(...))
  st <- data$strata
  keep <- rep(TRUE, nrow(st))
  if (length(cond) > 1) {
    res <- eval(cond, as.list(st), parent.frame())
    for (r in res) {
      if (!is.logical(r) || length(r) != nrow(st)) {
        mt_abort("Conditions must be logical expressions on stratum variables.", class = "argument")
      }
      keep <- keep & !is.na(r) & r
    }
  }
  if (!any(keep)) mt_abort("No stratum satisfies the conditions.", class = "argument")
  data <- strata_subset(data, keep)
  if (!is.null(ages)) {
    if (length(ages) != 2) mt_abort("`ages` must be c(lower, upper).", class = "argument")
    data$obs <- lapply(data$obs, function(o) o[age_lo >= ages[1] - 1e-9 & age_hi <= ages[2] + 1e-9])
    if (!is.null(data$grid)) {
      a <- data$grid$age
      k <- a$lower >= ages[1] - 1e-9 & a$upper <= ages[2] + 1e-9
      if (!any(k)) mt_abort("No grid age is inside `ages`.", class = "grid")
      data$grid$age <- mt_axis(a$lower[k], a$upper[k], a$labels[k], kind = "age")
    }
  }
  if (!is.null(periods)) {
    if (length(periods) != 2) mt_abort("`periods` must be c(first, last).", class = "argument")
    data <- mt_window(data, if (is.na(periods[1])) -Inf else periods[1],
                      if (is.na(periods[2])) Inf else periods[2])
  }
  data$obs <- data$obs[vapply(data$obs, nrow, 1L) > 0L]
  if (!length(data$obs)) mt_abort("No observations left.", class = "argument")
  if (drop) {
    single <- data$by[vapply(data$by, function(b) length(unique(data$strata[[b]])) == 1L, logical(1))]
    if (length(single)) {
      data$by <- setdiff(data$by, single)
      data$strata[, (single) := NULL]
      lab <- if (length(data$by)) do.call(paste, c(lapply(data$by, function(b) data$strata[[b]]), sep = ".")) else "all"
      data$strata[, label := lab]
    }
  }
  data
}

# Combine levels of `var` into groups (named list of member levels)
regroup_strata <- function(x, var, groups, additive = "deaths") {
  st <- copy(x$strata)
  others <- setdiff(x$by, var)
  map <- data.table(level = unlist(groups, use.names = FALSE),
                    group = rep(names(groups), lengths(groups)))
  st[, group_ := map$group[match(st[[var]], map$level)]]
  st <- st[!is.na(group_)]
  st[, group_ := factor(group_, levels = names(groups))]
  # keep the order of the strata, with groups in the order given
  for (b in x$by) {
    set(st, j = paste0("o_", b), value = if (b == var) as.integer(st$group_) else match(st[[b]], unique(st[[b]])))
  }
  setorderv(st, c(paste0("o_", x$by), "stratum"))
  st[, new_ := .GRP, by = c("group_", others)]
  st[, n_ := .N, by = new_]
  st[, first_ := stratum == stratum[1], by = new_]
  dropped <- 0L
  for (m in names(x$obs)) {
    o <- merge(x$obs[[m]], st[, list(stratum, new_, n_, first_)], by = "stratum")
    cells <- c("new_", "age_lo", "age_hi", "period_lo", "period_hi")
    if (m != "deaths" && additive == "rates") {
      a <- o[first_ == TRUE, c(cells, "value", "provenance"), with = FALSE]
    } else {
      a <- o[, list(value = sum(value), n = .N, members = n_[1],
                    provenance = if (.N > 1) "aggregated" else provenance[1]), by = cells]
      dropped <- dropped + sum(a$n < a$members)
      a <- a[n == members][, c("n", "members") := NULL]
    }
    setnames(a, "new_", "stratum")
    setcolorder(a, c("stratum", "age_lo", "age_hi", "period_lo", "period_hi", "value", "provenance"))
    setorderv(a, c("stratum", "period_lo", "age_lo"))
    x$obs[[m]] <- a
  }
  if (dropped) {
    mt_warn(sprintf("%d cell(s) dropped because not every member level was observed.", dropped), class = "group")
  }
  new <- unique(st[, c("new_", "group_", others), with = FALSE], by = "new_")
  setorderv(new, "new_")
  new[, (var) := as.character(group_)]
  new[, group_ := NULL]
  setnames(new, "new_", "stratum")
  new[, label := do.call(paste, c(lapply(x$by, function(b) new[[b]]), sep = "."))]
  setcolorder(new, c("stratum", x$by, "label"))
  x$strata <- new[]
  x
}

strata_subset <- function(x, keep) {
  st <- x$strata[keep]
  map <- stats::setNames(seq_len(nrow(st)), st$stratum)
  for (m in names(x$obs)) {
    o <- x$obs[[m]][stratum %in% st$stratum]
    o[, stratum := as.integer(map[as.character(stratum)])]
    x$obs[[m]] <- o
  }
  st[, stratum := seq_len(.N)]
  x$strata <- st
  x
}
