# EcoTuneR 2.0
# Compare food web model (Ecopath with Ecosim) trophic levels with Bayesian stable
# isotope (SIA) trophic levels.
#
# Quick start:
#   source('https://raw.githubusercontent.com/LewisLab-URI/EcoTuneR/main/EcoTuneR.R')
#   ecotuneR()                      # opens the app; upload a CSV or load the example
#   ecotuneR(combined_output)       # opens the app with your data frame
#   ecotuneR(combined_output, export = TRUE)   # writes a report and plots to a folder
#
# Input: a data frame with 5 columns in this order, plus an optional 6th:
#   1 Consumer (functional group)   2 Lower bound of the 95% credibility interval
#   3 Upper bound                   4 SIA mode (trophic position)
#   5 Ecopath trophic level (NA if none)   6 SIA sample size (optional)

local({
  need <- c("shiny", "bslib", "ggplot2")
  missing_pkgs <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("EcoTuneR needs the following R package(s): ", paste(missing_pkgs, collapse = ", "),
         ".\nInstall with:\n  install.packages(c(",
         paste0('"', missing_pkgs, '"', collapse = ", "), "))", call. = FALSE)
  }
})
library(shiny)
library(bslib)
library(ggplot2)

et_version <- "2.0"
et_repo_raw <- "https://raw.githubusercontent.com/LewisLab-URI/EcoTuneR/main"
et_citation <- "Lewis et al. (in review). Bayesian stable isotope analysis as a validation approach for marine food web models: a coastal Louisiana, USA case study. Ecosphere."

et_status_levels <- c("Within 95% CI", "Overestimated", "Underestimated", "No model trophic level")
et_status_colors <- c("Within 95% CI" = "#2B7BBA", "Overestimated" = "#E08A00",
                      "Underestimated" = "#7B5EA7", "No model trophic level" = "#9AA3AD")
et_roles <- c("Consumer (functional group)", "Lower bound of 95% CI", "Upper bound of 95% CI",
              "SIA mode (trophic position)", "Ecopath trophic level", "SIA sample size (optional)")

# ---------------------------------------------------------------- data checks

et_as_number <- function(x) {
  # converts a column to numbers; returns the values and which entries were not numbers
  if (is.numeric(x)) return(list(value = as.numeric(x), bad = rep(FALSE, length(x))))
  txt <- trimws(as.character(x))
  blank <- is.na(txt) | txt == "" | toupper(txt) %in% c("NA", "N/A", "NAN", "NULL", "-")
  value <- suppressWarnings(as.numeric(txt))
  list(value = value, bad = is.na(value) & !blank)
}

et_clean_names <- function(x) trimws(gsub("_", " ", as.character(x)))

et_validate <- function(df, cols = NULL) {
  # checks an input table and returns list(ok, errors, warnings, data)
  # cols: which columns of df hold the six roles, in order (6th may be NA)
  errors <- character(0); warnings <- character(0)
  fail <- function(msg) list(ok = FALSE, errors = c(errors, msg), warnings = warnings, data = NULL)
  if (!is.data.frame(df)) return(fail("The input is not a data frame."))
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (ncol(df) < 5) {
    return(fail(paste0("The input has ", ncol(df), " column(s). EcoTuneR needs at least 5: ",
                       "consumer, lower bound, upper bound, SIA mode, Ecopath trophic level.")))
  }
  if (is.null(cols)) cols <- c(1:5, if (ncol(df) >= 6) 6 else NA)
  cols <- suppressWarnings(as.integer(cols))
  if (length(cols) < 6) cols <- c(cols, rep(NA, 6 - length(cols)))
  if (any(is.na(cols[1:5])) || any(cols[1:5] < 1 | cols[1:5] > ncol(df))) {
    return(fail("Each of the five required columns must be assigned."))
  }
  used <- cols[!is.na(cols)]
  if (anyDuplicated(used)) return(fail("The same column is assigned to more than one role."))

  keep <- rowSums(!is.na(df) & trimws(as.matrix(df)) != "") > 0   # drop fully empty rows
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) return(fail("The input has no rows of data."))

  consumer <- et_clean_names(df[[cols[1]]])
  rows <- function(i) paste0(utils::head(consumer[i], 6), collapse = ", ")
  more <- function(i) if (sum(i) > 6) paste0(" and ", sum(i) - 6, " more") else ""

  num <- lapply(cols[2:5], function(j) et_as_number(df[[j]]))
  labels <- c("lower bound", "upper bound", "SIA mode", "Ecopath trophic level")
  for (k in 1:4) {
    if (any(num[[k]]$bad)) {
      errors <- c(errors, paste0("The ", labels[k], " column has values that are not numbers (",
                                 rows(num[[k]]$bad), more(num[[k]]$bad), ")."))
    }
  }
  n <- rep(NA_real_, nrow(df))
  if (!is.na(cols[6])) {
    nn <- et_as_number(df[[cols[6]]])
    if (any(nn$bad)) warnings <- c(warnings, paste0("Some sample sizes are not numbers and were ignored (",
                                                    rows(nn$bad), more(nn$bad), ")."))
    n <- nn$value
  }
  if (length(errors) > 0) return(list(ok = FALSE, errors = errors, warnings = warnings, data = NULL))

  d <- data.frame(Consumer = consumer, Lower = num[[1]]$value, Upper = num[[2]]$value,
                  Mode = num[[3]]$value, TL = num[[4]]$value, N = n, stringsAsFactors = FALSE)

  no_name <- is.na(d$Consumer) | d$Consumer == ""
  if (any(no_name)) {
    warnings <- c(warnings, paste0(sum(no_name), " row(s) with no consumer name were removed."))
    d <- d[!no_name, , drop = FALSE]
  }
  incomplete <- is.na(d$Lower) | is.na(d$Upper) | is.na(d$Mode)
  if (any(incomplete)) {
    warnings <- c(warnings, paste0("Removed ", sum(incomplete), " group(s) missing a lower bound, upper bound or SIA mode (",
                                   paste0(utils::head(d$Consumer[incomplete], 6), collapse = ", "),
                                   if (sum(incomplete) > 6) paste0(" and ", sum(incomplete) - 6, " more") else "", ")."))
    d <- d[!incomplete, , drop = FALSE]
  }
  if (nrow(d) == 0) return(list(ok = FALSE, errors = "No complete rows of data were found.",
                                warnings = warnings, data = NULL))
  flipped <- d$Lower > d$Upper
  if (any(flipped)) {
    errors <- c(errors, paste0("The lower bound is greater than the upper bound for: ",
                               paste0(utils::head(d$Consumer[flipped], 6), collapse = ", "),
                               if (sum(flipped) > 6) paste0(" and ", sum(flipped) - 6, " more") else "",
                               ". Check that the columns are in the right order."))
    return(list(ok = FALSE, errors = errors, warnings = warnings, data = NULL))
  }
  outside <- d$Mode < d$Lower | d$Mode > d$Upper
  if (any(outside)) {
    warnings <- c(warnings, paste0("The SIA mode is outside its own credibility interval for: ",
                                   paste0(utils::head(d$Consumer[outside], 6), collapse = ", "),
                                   if (sum(outside) > 6) paste0(" and ", sum(outside) - 6, " more") else "", "."))
  }
  if (anyDuplicated(d$Consumer)) {
    dup <- unique(d$Consumer[duplicated(d$Consumer)])
    warnings <- c(warnings, paste0("Repeated consumer names were numbered to keep them apart (",
                                   paste0(utils::head(dup, 6), collapse = ", "), ")."))
    d$Consumer <- make.unique(d$Consumer, sep = " ")
  }
  if (all(is.na(d$TL))) warnings <- c(warnings, "No Ecopath trophic levels were found, so no comparison can be made.")
  if (all(is.na(d$N))) d$N <- NA_real_
  rownames(d) <- NULL
  list(ok = TRUE, errors = character(0), warnings = warnings, data = d)
}

et_classify <- function(d) {
  # adds the validation outcome for each functional group
  status <- rep("No model trophic level", nrow(d))
  has <- !is.na(d$TL)
  status[has & d$TL <= d$Upper & d$TL >= d$Lower] <- "Within 95% CI"
  status[has & d$TL > d$Upper] <- "Overestimated"
  status[has & d$TL < d$Lower] <- "Underestimated"
  d$Status <- factor(status, levels = et_status_levels)
  d$Difference <- d$TL - d$Mode
  d
}

et_counts <- function(d) {
  tab <- table(d$Status)
  with_tl <- sum(!is.na(d$TL))
  list(total = nrow(d), with_tl = with_tl, within = tab[["Within 95% CI"]],
       over = tab[["Overestimated"]], under = tab[["Underestimated"]],
       none = tab[["No model trophic level"]],
       pct = if (with_tl > 0) round(100 * tab[["Within 95% CI"]] / with_tl) else NA)
}

et_validate_posterior <- function(post, groups) {
  # posterior draws: column 1 = consumer, column 2 = one trophic position draw per row
  if (is.null(post)) return(list(data = NULL, message = NULL, ok = TRUE))
  if (!is.data.frame(post) || ncol(post) < 2) {
    return(list(data = NULL, ok = FALSE,
                message = "The posterior file needs two columns: consumer, then one trophic position draw per row."))
  }
  tp <- et_as_number(post[[2]])
  p <- data.frame(Consumer = et_clean_names(post[[1]]), TP = tp$value, stringsAsFactors = FALSE)
  p <- p[!is.na(p$TP) & p$Consumer %in% groups, , drop = FALSE]
  if (nrow(p) == 0) {
    return(list(data = NULL, ok = FALSE,
                message = "No posterior draws matched the consumer names in the main table."))
  }
  matched <- length(unique(p$Consumer))
  list(data = p, ok = TRUE,
       message = paste0("Posterior draws read for ", matched, " of ", length(groups),
                        " groups. Other groups use the approximate curve."))
}

et_example_data <- function() {
  local_file <- file.path("Example Datasets", "Example_combined_output.csv")
  if (file.exists(local_file)) return(utils::read.csv(local_file, check.names = FALSE))
  if (file.exists("Example_combined_output.csv")) return(utils::read.csv("Example_combined_output.csv", check.names = FALSE))
  utils::read.csv(paste0(et_repo_raw, "/Example%20Datasets/Example_combined_output.csv"), check.names = FALSE)
}

# ---------------------------------------------------------------------- plots

et_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#E6E8EE"),
          axis.title = element_text(color = "#33364D"),
          axis.text = element_text(color = "#33364D"),
          plot.title = element_text(face = "bold", color = "#12112A"),
          plot.subtitle = element_text(color = "#5A5E73"),
          plot.title.position = "plot",
          legend.position = "top", legend.justification = "left",
          legend.title = element_blank())
}

et_plot_summary <- function(d) {
  # one stacked bar: how many groups fall in each outcome
  tab <- as.data.frame(table(Status = d$Status), stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, , drop = FALSE]
  tab$Status <- factor(tab$Status, levels = rev(et_status_levels))
  tab$label <- ifelse(tab$Freq / sum(tab$Freq) < 0.08, "",
                      paste0(tab$Freq, " (", round(100 * tab$Freq / sum(tab$Freq)), "%)"))
  ggplot(tab, aes(x = Freq, y = "", fill = Status)) +
    geom_col(width = 0.6, color = "white", linewidth = 1) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), color = "white",
              fontface = "bold", size = 4.6) +
    scale_fill_manual(values = et_status_colors, breaks = et_status_levels) +
    labs(x = "Number of functional groups", y = NULL) +
    et_theme() +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank())
}

et_plot_overview <- function(d, order_by = "SIA trophic level") {
  # every group on one chart: SIA interval and mode, with the model trophic level on top
  ord <- switch(order_by,
                "Name" = order(d$Consumer, decreasing = TRUE),
                "Outcome" = order(-as.integer(d$Status), d$Mode),
                "Difference (model minus SIA)" = order(d$Difference, na.last = FALSE),
                order(d$Mode))
  d$Consumer <- factor(d$Consumer, levels = d$Consumer[ord])
  with_tl <- d[!is.na(d$TL), , drop = FALSE]
  ggplot(d, aes(y = Consumer)) +
    geom_linerange(aes(xmin = Lower, xmax = Upper), color = "#B9BFCC", linewidth = 2.2, lineend = "round") +
    geom_point(aes(x = Mode, shape = "SIA mode (bar = 95% CI)"), color = "#12112A", size = 2.4) +
    (if (nrow(with_tl) > 0) geom_point(data = with_tl, aes(x = TL, fill = Status, shape = "Ecopath trophic level"),
                                       size = 3.4, color = "white", stroke = 0.6)) +
    scale_shape_manual(values = c("SIA mode (bar = 95% CI)" = 16, "Ecopath trophic level" = 23),
                       breaks = c("SIA mode (bar = 95% CI)", "Ecopath trophic level")) +
    scale_fill_manual(values = et_status_colors, breaks = et_status_levels, drop = TRUE) +
    guides(fill = guide_legend(override.aes = list(shape = 23, size = 4, color = "white")),
           shape = guide_legend(override.aes = list(fill = "#5A5E73", color = "#12112A", size = 3))) +
    labs(x = "Trophic level", y = NULL) +
    et_theme() +
    theme(panel.grid.major.y = element_blank(), legend.box = "vertical",
          legend.box.just = "left", legend.margin = margin(0, 0, 0, 0))
}

et_plot_detail <- function(d, group, posterior = NULL) {
  # distribution plot for one group. With posterior draws it is the real density;
  # without them it is an approximate curve drawn from the lower bound, upper bound and mode.
  g <- d[d$Consumer == group, , drop = FALSE][1, ]
  draws <- if (!is.null(posterior)) posterior$TP[posterior$Consumer == group] else numeric(0)
  real <- length(draws) >= 10
  values <- if (real) draws else c(g$Lower, g$Upper, g$Mode)
  marks <- c(g$Lower, g$Upper, g$Mode, g$TL)
  span <- diff(range(marks, na.rm = TRUE)); if (span == 0) span <- 0.5
  lines <- data.frame(value = c(g$Lower, g$Upper, g$Mode),
                      name = c("SIA 95% CI", "SIA 95% CI", "SIA mode (trophic position)"),
                      stringsAsFactors = FALSE)
  if (!is.na(g$TL)) lines <- rbind(lines, data.frame(value = g$TL, name = "Ecopath trophic level"))
  keys <- c("SIA 95% CI", "SIA mode (trophic position)", "Ecopath trophic level")
  lines$name <- factor(lines$name, levels = keys)
  tl_color <- unname(et_status_colors[as.character(g$Status)])
  subtitle <- if (real) {
    paste0("Posterior density from ", format(length(draws), big.mark = ","), " draws")
  } else {
    "Approximate curve drawn from the lower bound, upper bound and mode"
  }
  n_text <- if (!is.na(g$N)) paste0(" (n = ", g$N, ")") else ""
  p <- ggplot(data.frame(value = values), aes(value)) +
    geom_density(fill = "#7189D9", color = "#4A5FC1", alpha = 0.35, linewidth = 0.6) +
    geom_vline(data = lines, aes(xintercept = value, color = name, linetype = name), linewidth = 1) +
    scale_color_manual(values = c("SIA 95% CI" = "#5A5E73", "SIA mode (trophic position)" = "#12112A",
                                  "Ecopath trophic level" = tl_color), breaks = keys, drop = TRUE) +
    scale_linetype_manual(values = c("SIA 95% CI" = "dashed", "SIA mode (trophic position)" = "solid",
                                     "Ecopath trophic level" = "longdash"), breaks = keys, drop = TRUE) +
    labs(title = paste0(group, n_text), subtitle = subtitle, x = "Trophic level", y = "Density") +
    et_theme() +
    theme(legend.key.width = grid::unit(1.6, "cm"))
  if (!real) {
    p <- p + expand_limits(x = c(min(marks, na.rm = TRUE) - 0.8 * span, max(marks, na.rm = TRUE) + 0.8 * span))
  } else {
    p <- p + expand_limits(x = range(marks, na.rm = TRUE))
  }
  p
}

# ---------------------------------------------------------------- sample size

samplesizeFit <- function(combined_output, flag_cutoff = 2) {
  # fits a power function decay curve (CI range = a * n^b) using the optional
  # 6th column (SIA sample size) and flags groups that deviate from the curve.
  # returns NULL if sample size data are not available.
  if (ncol(combined_output) < 6) return(NULL)
  ss <- combined_output[, c(1, 2, 3, 6)]
  colnames(ss) <- c("Consumer", "PD.L", "PD.U", "N")
  ss$N <- suppressWarnings(as.numeric(ss$N))
  ss$Range <- ss$PD.U - ss$PD.L
  ss <- ss[!is.na(ss$N) & !is.na(ss$Range) & ss$N > 0 & ss$Range > 0, ]
  if (nrow(ss) < 5 || length(unique(ss$N)) < 3) return(NULL)

  fit <- lm(log(Range) ~ log(N), data = ss) # power function fit on the log-log scale
  a <- exp(coef(fit)[[1]])
  b <- coef(fit)[[2]]
  ss$Expected <- a * ss$N^b
  ss$Deviation <- log(ss$Range / ss$Expected)
  cutoff <- flag_cutoff * mad(ss$Deviation)
  ss$Flag <- "Near curve"
  ss$Flag[ss$Deviation > cutoff] <- "Wider than expected"
  ss$Flag[ss$Deviation < -cutoff] <- "Narrower than expected"
  list(data = ss, a = a, b = b, flag_cutoff = flag_cutoff)
}

samplesizePlot <- function(ss) {
  # plots CI range against sample size with the fitted curve from samplesizeFit
  curve_df <- data.frame(N = exp(seq(log(min(ss$data$N)), log(max(ss$data$N)), length.out = 200)))
  curve_df$Range <- ss$a * curve_df$N^ss$b
  flagged <- ss$data[ss$data$Flag != "Near curve", ]
  ggplot(ss$data, aes(x = N, y = Range)) +
    geom_line(data = curve_df, linetype = "dotted", linewidth = 1, color = "#12112A") +
    geom_point(aes(color = Flag), size = 3) +
    geom_text(data = flagged, aes(label = gsub("_", " ", Consumer)), vjust = -1, size = 4) +
    scale_color_manual(values = c("Near curve" = "#5A5E73",
                                  "Wider than expected" = "#D1495B",
                                  "Narrower than expected" = "#2B7BBA")) +
    scale_x_log10() +
    labs(x = "SIA sample size (n, log scale)", y = "95% CI range (upper minus lower)") +
    et_theme()
}

samplesizeTable <- function(ss) {
  # table of the groups flagged by samplesizeFit
  flagged <- ss$data[ss$data$Flag != "Near curve", ]
  data.frame("Functional Group" = gsub("_", " ", flagged$Consumer),
             "n" = as.integer(flagged$N),
             "CI Range" = round(flagged$Range, 2),
             "Expected CI Range" = round(flagged$Expected, 2),
             "Flag" = flagged$Flag,
             check.names = FALSE)
}

# --------------------------------------------------------------------- report

et_results_table <- function(d) {
  out <- data.frame("Functional group" = d$Consumer, "Lower" = d$Lower, "Upper" = d$Upper,
                    "SIA mode" = d$Mode, "Ecopath TL" = d$TL,
                    "Difference (Ecopath minus SIA)" = round(d$Difference, 2),
                    "Outcome" = as.character(d$Status), check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(is.na(d$N))) out[["SIA sample size"]] <- d$N
  out
}

et_png_uri <- function(plot, width, height, dpi = 120) {
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f))
  ggsave(f, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  base64enc::dataURI(file = f, mime = "image/png")
}

et_headline <- function(d) {
  k <- et_counts(d)
  if (is.na(k$pct)) return("No Ecopath trophic levels were supplied, so no comparison was made.")
  paste0(k$within, " of ", k$with_tl, " functional groups (", k$pct,
         "%) have an Ecopath trophic level within the 95% credibility interval of the SIA trophic level.")
}

et_report_html <- function(d, posterior = NULL, title = "EcoTuneR results", sample = FALSE,
                           detail_plots = TRUE, progress = NULL) {
  # one self-contained HTML file: summary, overview, results table, sample size, group plots
  esc <- htmltools::htmlEscape
  k <- et_counts(d)
  img <- function(uri, alt) paste0('<img src="', uri, '" alt="', esc(alt), '" style="max-width:100%;height:auto">')
  tbl <- et_results_table(d)
  row_html <- vapply(seq_len(nrow(tbl)), function(i) {
    col <- unname(et_status_colors[tbl$Outcome[i]])
    cells <- vapply(names(tbl), function(nm) {
      v <- tbl[[nm]][i]; v <- if (is.na(v)) "" else as.character(v)
      if (nm == "Outcome") paste0('<td><span class="dot" style="background:', col, '"></span>', esc(v), "</td>")
      else paste0("<td>", esc(v), "</td>")
    }, character(1))
    paste0("<tr>", paste0(cells, collapse = ""), "</tr>")
  }, character(1))
  parts <- c(
    "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><title>", esc(title), "</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;color:#12112A;max-width:980px;margin:32px auto;padding:0 20px;line-height:1.5}",
    "h1{font-size:28px;margin-bottom:4px}h2{font-size:20px;margin-top:36px;border-bottom:2px solid #7189D9;padding-bottom:4px}",
    "table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:6px 10px;border-bottom:1px solid #E6E8EE;text-align:left}",
    "th{background:#F3F4F9}.dot{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:6px}",
    ".note{background:#FFF6DD;border-left:4px solid #E0A800;padding:10px 14px;margin:16px 0}.muted{color:#5A5E73;font-size:14px}",
    ".headline{font-size:19px;font-weight:600;margin:18px 0}</style></head><body>",
    "<h1>", esc(title), "</h1><div class='muted'>Made with EcoTuneR ", et_version, " on ", format(Sys.Date(), "%d %B %Y"), "</div>",
    if (sample) "<div class='note'><strong>This is a sample.</strong> It uses example data from Barataria Bay, Louisiana, to show what EcoTuneR produces.</div>",
    "<div class='headline'>", esc(et_headline(d)), "</div>",
    "<p>Overestimated: ", k$over, ". Underestimated: ", k$under,
    if (k$none > 0) paste0(". No Ecopath trophic level: ", k$none), ".</p>",
    img(et_png_uri(et_plot_summary(d), 9, 2.4), "Summary of outcomes"),
    "<h2>All functional groups</h2>",
    img(et_png_uri(et_plot_overview(d), 9, max(4, 0.24 * nrow(d) + 1.6)), "Overview of all groups"),
    "<h2>Results table</h2><table><thead><tr>",
    paste0("<th>", esc(names(tbl)), "</th>", collapse = ""), "</tr></thead><tbody>",
    paste0(row_html, collapse = ""), "</tbody></table>")
  ss <- samplesizeFit(d[, c("Consumer", "Lower", "Upper", "Mode", "TL", "N")])
  if (!is.null(ss)) {
    flags <- samplesizeTable(ss)
    parts <- c(parts, "<h2>Sample size</h2><p class='muted'>Power function decay curve: CI range = ",
               round(ss$a, 2), " &times; n<sup>", round(ss$b, 2), "</sup>. Groups are flagged when their credibility interval range deviates from the curve by more than ",
               ss$flag_cutoff, " median absolute deviations on the log scale.</p>",
               img(et_png_uri(samplesizePlot(ss), 9, 5.2), "Sample size plot"),
               if (nrow(flags) > 0) c("<table><thead><tr>", paste0("<th>", esc(names(flags)), "</th>", collapse = ""),
                                      "</tr></thead><tbody>",
                                      paste0(apply(flags, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>")), collapse = ""),
                                      "</tbody></table>"))
  }
  if (detail_plots) {
    parts <- c(parts, "<h2>Functional group plots</h2>",
               "<p class='muted'>Curves are approximations drawn from the lower bound, upper bound and mode unless posterior draws were supplied.</p>")
    for (i in seq_len(nrow(d))) {
      if (!is.null(progress)) progress(i, nrow(d))
      parts <- c(parts, img(et_png_uri(et_plot_detail(d, d$Consumer[i], posterior), 8, 4.6, dpi = 100), d$Consumer[i]))
    }
  }
  paste0(c(parts, "<p class='muted'>Cite: ", esc(et_citation), "</p></body></html>"), collapse = "")
}

et_export_folder <- function(d, posterior = NULL, data_name = "combined_output") {
  # writes the report, results table and one PNG per functional group to a folder
  folder <- paste0("Exported Results of ", gsub("[^[:alnum:]._ -]", ".", data_name))
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  message("Exporting summary and plots to '", folder, "'")
  ggsave(file.path(folder, "Summary Plot.png"), et_plot_summary(d), width = 9, height = 2.4, dpi = 200, bg = "white")
  ggsave(file.path(folder, "Overview Plot.png"), et_plot_overview(d), width = 9,
         height = max(4, 0.24 * nrow(d) + 1.6), dpi = 200, bg = "white", limitsize = FALSE)
  utils::write.csv(et_results_table(d), file.path(folder, "Results Table.csv"), row.names = FALSE)
  ss <- samplesizeFit(d[, c("Consumer", "Lower", "Upper", "Mode", "TL", "N")])
  if (!is.null(ss)) {
    ggsave(file.path(folder, "Sample Size Plot.png"), samplesizePlot(ss), width = 9, height = 5.8, dpi = 200, bg = "white")
    utils::write.csv(samplesizeTable(ss), file.path(folder, "Sample Size Flags.csv"), row.names = FALSE)
  }
  for (w in seq_len(nrow(d))) {
    file_name <- gsub("[^[:alnum:]._[:space:]]", ".", paste0(d$Consumer[w], ".png"))
    ggsave(file.path(folder, file_name), et_plot_detail(d, d$Consumer[w], posterior),
           width = 9, height = 5.8, dpi = 200, bg = "white")
    cat(".")
  }
  writeLines(et_report_html(d, posterior, title = paste0("A Summary of ", data_name), detail_plots = FALSE),
             file.path(folder, paste0("A Summary of ", gsub("[^[:alnum:]._ -]", ".", data_name), ".html")))
  message("\nSummary and plots exported")
  if (interactive()) utils::browseURL(folder)
  invisible(folder)
}

# ------------------------------------------------------------------------- UI

et_css <- "
:root { --et-navy: #12112A; --et-blue: #4A5FC1; --et-periwinkle: #7189D9; --et-soft: #F3F4F9; }
body { color: var(--et-navy); background: #FAFAFC; }
.navbar { background: #ffffff !important; border-bottom: 1px solid #E6E8EE; padding-top: 6px; padding-bottom: 6px; }
.navbar-brand img { height: 44px; width: auto; }
.navbar .nav-link { font-weight: 600; color: #5A5E73 !important; padding: 8px 14px !important; border-radius: 8px; }
.navbar .nav-link.active { color: #ffffff !important; background: var(--et-blue); }
.btn-primary, .btn-default.et-primary { background: var(--et-blue); border-color: var(--et-blue); color: #fff; }
.btn-primary:hover, .btn-default.et-primary:hover { background: #3B4DA3; border-color: #3B4DA3; color: #fff; }
.card { border: 1px solid #E6E8EE; border-radius: 12px; box-shadow: 0 1px 2px rgba(18,17,42,.04); }
.card-header { background: #ffffff; font-weight: 700; border-bottom: 1px solid #E6E8EE; }
h1.et-page { font-size: 30px; font-weight: 700; margin: 6px 0 2px; }
p.et-lead { color: #5A5E73; max-width: 80ch; margin-bottom: 18px; }
.et-banner { background: #FFF6DD; border: 1px solid #E9CF7A; border-left: 6px solid #E0A800; border-radius: 10px;
  padding: 12px 16px; margin: 14px 0 6px; font-size: 16px; }
.et-headline { font-size: 22px; font-weight: 700; margin: 4px 0 14px; }
.et-stat { background: #fff; border: 1px solid #E6E8EE; border-radius: 12px; padding: 14px 16px; border-top: 5px solid var(--c); height: 100%; }
.et-stat .v { font-size: 34px; font-weight: 700; line-height: 1.1; }
.et-stat .t { color: #5A5E73; font-size: 14px; }
.et-table { font-size: 15px; margin-bottom: 0; }
.et-table tbody tr { cursor: pointer; }
.et-table tbody tr:hover td { background: #EEF1FB; }
.et-dot { display: inline-block; width: 11px; height: 11px; border-radius: 50%; margin-right: 7px; }
.et-scroll { max-height: 520px; overflow-y: auto; }
.et-msg { border-radius: 10px; padding: 10px 14px; margin-bottom: 10px; }
.et-msg.err { background: #FDECEC; border: 1px solid #F1B5B5; }
.et-msg.warn { background: #FFF6DD; border: 1px solid #E9CF7A; }
.et-msg.ok { background: #E8F5EE; border: 1px solid #A9D8BE; }
.et-footer { margin-top: 40px; padding: 22px 0 28px; border-top: 1px solid #E6E8EE; background: #fff; }
.et-logos { display: flex; flex-wrap: wrap; align-items: center; gap: 28px; margin-bottom: 12px; }
.et-logos img { height: 40px; width: auto; }
.et-logos img.tall { height: 84px; }
.et-footer .t { color: #5A5E73; font-size: 14px; }
.et-empty { color: #5A5E73; padding: 30px 0; font-size: 17px; }
.et-q { font-weight: 600; margin-bottom: 4px; }
"

et_www <- function() {
  # logos come from a local www folder when there is one, otherwise from GitHub
  if (file.exists(file.path("www", "LSU.png"))) {
    addResourcePath("etwww", normalizePath("www"))
    "etwww"
  } else paste0(et_repo_raw, "/www")
}

et_ui <- function(www) {
  logo <- function(file, href, alt, class = NULL) tags$a(href = href, target = "_blank", rel = "noopener",
                                                         tags$img(src = paste0(www, "/", file), alt = alt, class = class))
  footer <- tags$div(class = "et-footer", tags$div(class = "container-fluid",
    tags$div(class = "et-logos",
      logo("LewisLab.png", "https://web.uri.edu/lewis-lab/", "Lewis Lab", class = "tall"),
      logo("LSU.png", "https://www.lsu.edu/", "Louisiana State University"),
      logo("USM.png", "https://www.usm.edu/", "The University of Southern Mississippi"),
      logo("UCSC.png", "https://www.ucsc.edu/", "University of California Santa Cruz"),
      logo("NASEM.png", "https://www.nationalacademies.org/gulf/gulf-research-program", "National Academies Gulf Research Program")),
    tags$div(class = "t", "This research was funded by the National Academy of Sciences Gulf Research Program. EcoTuneR ",
             et_version, " is developed by the ",
             tags$a(href = "https://web.uri.edu/lewis-lab/", target = "_blank", rel = "noopener", "Lewis Lab", .noWS = "after"),
             ", University of Rhode Island Graduate School of Oceanography. ",
             tags$a(href = "https://github.com/LewisLab-URI/EcoTuneR", target = "_blank", rel = "noopener", "Code and instructions on GitHub", .noWS = "after"), ".")))

  page_navbar(
    id = "nav",
    title = tags$img(src = paste0(www, "/ecotuner0.png"), alt = "EcoTuneR"),
    window_title = "EcoTuneR",
    theme = bs_theme(version = 5),
    fillable = FALSE,
    header = tagList(tags$head(tags$style(HTML(et_css)), tags$script(HTML(
      "$(document).on('shiny:connected', function() {
         var q = ''; try { q = window.top.location.search; } catch (e) { q = window.location.search; }
         if (/[?&]example=1/.test(q)) Shiny.setInputValue('url_example', 1);
       });"))), uiOutput("banner")),
    footer = footer,

    nav_panel("Data", value = "Data",
      tags$h1(class = "et-page", "Data"),
      tags$p(class = "et-lead", "EcoTuneR compares the trophic levels from a food web model (Ecopath with Ecosim) with the trophic levels estimated from stable isotope data using a Bayesian model. Upload your combined table, or load the example to see how it works."),
      layout_columns(col_widths = c(4, 8),
        tagList(
          card(card_header("1. Load a table"),
               fileInput("file", "Upload a CSV file", accept = c(".csv", "text/csv")),
               actionButton("load_example", "Load the Barataria Bay example", class = "et-primary"),
               tags$div(style = "margin-top:10px", downloadLink("dl_example", "Download the example file to use as a template"))),
          card(card_header("Optional: posterior draws"),
               tags$p(class = "t", style = "color:#5A5E73;font-size:14px",
                      "Without this file, each group's curve is an approximation drawn from its lower bound, upper bound and mode. To show the real posterior distribution, upload a CSV with two columns: consumer, then one trophic position draw per row."),
               fileInput("post_file", NULL, accept = c(".csv", "text/csv")),
               uiOutput("post_msg"))),
        tagList(
          card(card_header("2. Check the columns"),
               tags$p(style = "color:#5A5E73;font-size:14px", "Columns are read by position. Column names do not need to match. Change an assignment below if your columns are in a different order."),
               uiOutput("mapping"),
               uiOutput("messages"),
               uiOutput("go_results")),
          card(card_header("First rows of your table"), tags$div(style = "overflow-x:auto", tableOutput("preview")))))),

    nav_panel("Summary", value = "Summary",
      tags$h1(class = "et-page", "Summary"),
      uiOutput("summary_body")),

    nav_panel("Details", value = "Details",
      tags$h1(class = "et-page", "Details"),
      uiOutput("details_body")),

    nav_panel("Sample Size", value = "Sample Size",
      tags$h1(class = "et-page", "Sample Size"),
      uiOutput("samplesize_body")),

    nav_panel("Interpretation", value = "Interpretation",
      tags$h1(class = "et-page", "Considerations for interpreting results"),
      tags$p(class = "et-lead", "EcoTuneR is a visualization tool and certain trends within the data may not be apparent when using this tool. Understanding the limitations of using Trophic Level data to validate an ecosystem model is important for a balanced interpretation of results. Please take time to thoroughly analyze results, and refer to the suggestions found in the table below. Species probability distributions are approximations due to the constraints of model outputs."),
      layout_columns(col_widths = c(6, 6),
        card(card_header("Trophic level is within bounds of 95% CI"),
             tags$div(class = "et-q", "Are any of the credibility interval ranges larger than 5?"),
             tags$p("If yes, large credibility intervals increase the chance of a trophic level falling within the bounds. This may not be a result of agreement between the two methods."),
             tags$div(class = "et-q", "Do all groups have an adequate sample size within the stable isotope data?"),
             tags$p("If no, small sample sizes of any relevant groups may lead to large credibility intervals (see above).")),
        card(card_header("Trophic level is outside bounds of 95% CI"),
             tags$div(class = "et-q", "If the system is detritus based, do any species rely heavily on detritus as a part of their diet? (>50%)"),
             tags$p("If yes, the detritus trophic level of 1 (EwE standard) may be offsetting species trophic levels to be lower than they are in the natural system."),
             tags$div(class = "et-q", "Were the stable isotope samples muscle or bone? How much diet information was integrated into the samples?"),
             tags$p("If muscle, these samples only integrate about 3 months of diet information and may not agree with EwE data collected over a longer period of time. Conversely, bone sample integrate up to 1 year of diet information"),
             tags$div(class = "et-q", "Were the SIA and EwE data collected in generally the same time frame (i.e., year, decade)?"),
             tags$p("If no, there may be lack of agreement in the data due to changes over time in the system of interest (i.e., fishing, construction, climate change)"),
             tags$div(class = "et-q", "What season were the SIA and EwE samples collected in? Is there any mismatch (i.e., summer vs fall)?"),
             tags$p("If yes, there may be lack of agreement in the data if data were collected in different seasons because of ecological changes in systems over the course of a year."))),
      card(card_header("How to cite"),
           tags$p(et_citation),
           tags$p(style = "color:#5A5E73;font-size:14px;margin-bottom:0",
                  "Built with shiny (Chang et al.), bslib (Sievert et al.) and ggplot2 (Wickham 2016). Input tables are typically produced with tRophicPosition (Quezada-Romegialli et al. 2018).")))
  )
}

# --------------------------------------------------------------------- server

et_server <- function(init = NULL, init_posterior = NULL, start_example = FALSE) {
  function(input, output, session) {
    rv <- reactiveValues(raw = init, source = if (is.null(init)) "none" else "r",
                         post_raw = init_posterior, read_error = NULL)

    load_example <- function() {
      ex <- tryCatch(et_example_data(), error = function(e) NULL)
      if (is.null(ex)) { rv$read_error <- "The example file could not be loaded. Check your internet connection."; return(FALSE) }
      rv$read_error <- NULL; rv$raw <- ex; rv$source <- "example"; rv$post_raw <- NULL
      TRUE
    }
    read_csv <- function(path) {
      d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
      if (ncol(d) == 1 && any(grepl(";", c(names(d), d[[1]])))) d <- utils::read.csv2(path, check.names = FALSE, stringsAsFactors = FALSE)
      d
    }

    observeEvent(input$load_example, { if (load_example()) nav_select("nav", "Summary") })
    observeEvent(input$url_example, { if (rv$source == "none" && load_example()) nav_select("nav", "Summary") })
    observeEvent(input$file, {
      d <- tryCatch(read_csv(input$file$datapath), error = function(e) NULL)
      if (is.null(d)) { rv$read_error <- "That file could not be read as a CSV table."; return() }
      rv$read_error <- NULL; rv$raw <- d; rv$source <- "user"; rv$post_raw <- NULL
    })
    observeEvent(input$post_file, {
      rv$post_raw <- tryCatch(read_csv(input$post_file$datapath), error = function(e) data.frame())
    })
    observe({
      q <- parseQueryString(isolate(session$clientData$url_search))
      if ((identical(q$example, "1") || start_example) && isolate(rv$source) == "none") {
        if (load_example()) nav_select("nav", "Summary")
      } else if (isolate(rv$source) == "r") nav_select("nav", "Summary")
    })

    cols <- reactive({
      req(rv$raw)
      default <- c(1:5, if (ncol(rv$raw) >= 6) 6 else NA)
      picked <- suppressWarnings(as.integer(vapply(1:6, function(i) {
        v <- input[[paste0("map", i)]]; if (is.null(v) || v == "") NA_character_ else v
      }, character(1))))
      if (all(is.na(picked)) || any(picked > ncol(rv$raw), na.rm = TRUE) || any(is.na(picked[1:5]))) default else picked
    })
    validated <- reactive({ req(rv$raw); et_validate(rv$raw, cols()) })
    dat <- reactive({ v <- validated(); req(v$ok); et_classify(v$data) })
    has_data <- reactive(!is.null(rv$raw) && isTRUE(validated()$ok))
    post <- reactive({
      if (is.null(rv$post_raw) || !has_data()) return(list(data = NULL, message = NULL, ok = TRUE))
      et_validate_posterior(rv$post_raw, dat()$Consumer)
    })
    ss <- reactive({ d <- dat(); samplesizeFit(d[, c("Consumer", "Lower", "Upper", "Mode", "TL", "N")]) })

    observe({
      if (has_data() && !all(is.na(dat()$N))) nav_show("nav", "Sample Size") else nav_hide("nav", "Sample Size")
    })

    output$banner <- renderUI({
      if (rv$source != "example") return(NULL)
      tags$div(class = "container-fluid", tags$div(class = "et-banner",
        tags$strong("This is a sample. "),
        "You are looking at example data from Barataria Bay, Louisiana, to show what EcoTuneR produces. Upload your own data on the Data tab to see your results."))
    })

    # ---- Data tab
    output$mapping <- renderUI({
      if (is.null(rv$raw)) return(tags$div(class = "et-empty", "No table loaded yet."))
      choices <- stats::setNames(seq_along(rv$raw), paste0(seq_along(rv$raw), ": ", names(rv$raw)))
      default <- c(1:5, if (ncol(rv$raw) >= 6) 6 else NA)
      layout_columns(col_widths = c(4, 4, 4), !!!lapply(1:6, function(i) {
        ch <- if (i == 6) c("(none)" = "", choices) else choices
        selectInput(paste0("map", i), et_roles[i], choices = ch,
                    selected = if (is.na(default[i])) "" else default[i])
      }))
    })
    output$messages <- renderUI({
      if (!is.null(rv$read_error)) return(tags$div(class = "et-msg err", rv$read_error))
      if (is.null(rv$raw)) return(NULL)
      v <- validated()
      tagList(
        lapply(v$errors, function(m) tags$div(class = "et-msg err", tags$strong("Problem: "), m)),
        lapply(v$warnings, function(m) tags$div(class = "et-msg warn", tags$strong("Note: "), m)),
        if (v$ok) {
          k <- et_counts(et_classify(v$data))
          tags$div(class = "et-msg ok", tags$strong("Table read. "),
                   paste0(k$total, " functional groups, ", k$with_tl, " with an Ecopath trophic level",
                          if (!all(is.na(v$data$N))) ", sample sizes included" else ", no sample sizes", "."))
        })
    })
    output$go_results <- renderUI({
      if (!has_data()) return(NULL)
      actionButton("go_summary", "View results", class = "et-primary")
    })
    observeEvent(input$go_summary, nav_select("nav", "Summary"))
    output$preview <- renderTable({
      validate(need(!is.null(rv$raw), "No table loaded yet."))
      utils::head(rv$raw, 8)
    }, na = "")
    output$post_msg <- renderUI({
      p <- post()
      if (is.null(p$message)) return(NULL)
      tags$div(class = paste("et-msg", if (p$ok) "ok" else "err"), p$message)
    })
    output$dl_example <- downloadHandler(
      filename = "Example_combined_output.csv",
      content = function(file) utils::write.csv(et_example_data(), file, row.names = FALSE))

    empty <- function() tags$div(class = "et-empty", "Load a table on the Data tab to see results here.")

    # ---- Summary tab
    output$summary_body <- renderUI({
      if (!has_data()) return(empty())
      d <- dat(); k <- et_counts(d)
      stat <- function(value, label, color) tags$div(class = "et-stat", style = paste0("--c:", color),
                                                     tags$div(class = "v", value), tags$div(class = "t", label))
      tagList(
        tags$div(class = "et-headline", et_headline(d)),
        layout_columns(col_widths = if (k$none > 0) c(2, 3, 2, 2, 3) else c(3, 3, 3, 3),
          stat(k$total, "functional groups", "#12112A"),
          stat(if (is.na(k$pct)) "-" else paste0(k$pct, "%"), "of model trophic levels within the SIA 95% CI", et_status_colors[[1]]),
          stat(k$over, "overestimated by the model", et_status_colors[[2]]),
          stat(k$under, "underestimated by the model", et_status_colors[[3]]),
          if (k$none > 0) stat(k$none, "with no model trophic level", et_status_colors[[4]])),
        card(card_header("Outcome of the comparison"), plotOutput("summary_plot", height = "190px")),
        card(card_header("All functional groups"),
             tags$p(style = "color:#5A5E73;font-size:14px;margin-bottom:6px",
                    "Each gray bar is the 95% credibility interval of the SIA trophic level and the dark point is its mode. The diamond is the Ecopath trophic level, colored by outcome."),
             selectInput("order_by", "Order by", width = "280px",
                         c("SIA trophic level", "Difference (model minus SIA)", "Outcome", "Name")),
             plotOutput("overview_plot", height = "auto")),
        card(card_header("Results by functional group"),
             tags$p(style = "color:#5A5E73;font-size:14px;margin-bottom:6px", "Click a row to open that group on the Details tab."),
             radioButtons("status_filter", NULL, inline = TRUE,
                          choices = c("All", et_status_levels[et_status_levels %in% as.character(d$Status)])),
             tags$div(class = "et-scroll", uiOutput("group_table"))),
        card(card_header("Download"),
             tags$div(downloadButton("dl_report", "Download report (HTML)", class = "et-primary"),
                      downloadButton("dl_table", "Download results table (CSV)", style = "margin-left:8px")),
             tags$p(style = "color:#5A5E73;font-size:14px;margin:8px 0 0", "The report is a single file with the summary, every group's plot and the results table. It can take a minute to build.")))
    })
    output$summary_plot <- renderPlot(et_plot_summary(dat()), res = 96)
    output$overview_plot <- renderPlot(
      et_plot_overview(dat(), if (is.null(input$order_by)) "SIA trophic level" else input$order_by),
      height = function() if (has_data()) max(420, 24 * nrow(dat()) + 170) else 420, res = 96)
    output$group_table <- renderUI({
      d <- dat()
      f <- if (is.null(input$status_filter)) "All" else input$status_filter
      if (f != "All") d <- d[as.character(d$Status) == f, , drop = FALSE]
      show_n <- !all(is.na(dat()$N))
      num <- function(x) if (is.na(x)) "" else format(round(x, 2), nsmall = 2)
      tags$table(class = "table table-sm et-table",
        tags$thead(tags$tr(tags$th("Functional group"), tags$th("Outcome"), tags$th("SIA 95% CI"),
                           tags$th("SIA mode"), tags$th("Ecopath TL"), tags$th("Difference"), if (show_n) tags$th("n"))),
        tags$tbody(lapply(seq_len(nrow(d)), function(i) {
          tags$tr(onclick = sprintf("Shiny.setInputValue('goto_group', %s, {priority: 'event'})",
                                    jsonlite::toJSON(d$Consumer[i], auto_unbox = TRUE)),
                  tags$td(d$Consumer[i]),
                  tags$td(tags$span(class = "et-dot", style = paste0("background:", et_status_colors[[as.character(d$Status[i])]])),
                          as.character(d$Status[i])),
                  tags$td(paste0(num(d$Lower[i]), " to ", num(d$Upper[i]))),
                  tags$td(num(d$Mode[i])), tags$td(num(d$TL[i])),
                  tags$td(if (is.na(d$Difference[i])) "" else sprintf("%+.2f", d$Difference[i])),
                  if (show_n) tags$td(if (is.na(d$N[i])) "" else d$N[i]))
        })))
    })
    observeEvent(input$goto_group, {
      rv$goto <- input$goto_group
      updateSelectInput(session, "consumer", selected = input$goto_group)
      nav_select("nav", "Details")
    })
    output$dl_table <- downloadHandler(
      filename = "EcoTuneR_results.csv",
      content = function(file) utils::write.csv(et_results_table(dat()), file, row.names = FALSE))
    output$dl_report <- downloadHandler(
      filename = "EcoTuneR_report.html",
      content = function(file) {
        withProgress(message = "Building report", value = 0, {
          html <- et_report_html(dat(), post()$data, sample = rv$source == "example",
                                 progress = function(i, n) setProgress(i / n, detail = paste("plot", i, "of", n)))
          writeLines(html, file)
        })
      })

    # ---- Details tab
    output$details_body <- renderUI({
      if (!has_data()) return(empty())
      groups <- dat()$Consumer
      selected <- isolate(if (!is.null(rv$goto) && rv$goto %in% groups) rv$goto else groups[1])
      tagList(
        card(
          layout_columns(col_widths = c(6, 6),
            selectInput("consumer", "Functional group", groups, selected = selected, width = "100%"),
            tags$div(style = "padding-top:32px",
                     actionButton("back", "Previous"), actionButton("forward", "Next", style = "margin-left:6px"),
                     downloadButton("dl_plot", "Download plot (PNG)", style = "margin-left:6px"))),
          plotOutput("detail_plot", height = "460px"),
          uiOutput("detail_note")),
        card(card_header("Numbers for this group"), tableOutput("detail_table")))
    })
    current <- reactive({
      d <- dat(); g <- input$consumer
      if (is.null(g) || !(g %in% d$Consumer)) d$Consumer[1] else g
    })
    step <- function(by) {
      groups <- dat()$Consumer
      i <- min(max(match(current(), groups) + by, 1), length(groups))
      updateSelectInput(session, "consumer", selected = groups[i])
    }
    observeEvent(input$forward, step(1))
    observeEvent(input$back, step(-1))
    output$detail_plot <- renderPlot(et_plot_detail(dat(), current(), post()$data), res = 96)
    output$detail_note <- renderUI({
      draws <- if (is.null(post()$data)) 0 else sum(post()$data$Consumer == current())
      if (draws >= 10) return(NULL)
      tags$p(style = "color:#5A5E73;font-size:14px;margin:6px 0 0",
             "This curve is an approximation drawn from three numbers (lower bound, upper bound and mode), not the full posterior distribution. Upload posterior draws on the Data tab to show the real distribution.")
    })
    output$detail_table <- renderTable({
      g <- dat()[dat()$Consumer == current(), , drop = FALSE][1, ]
      out <- data.frame("Lower bound" = g$Lower, "Upper bound" = g$Upper, "SIA mode" = g$Mode,
                        "Ecopath trophic level" = g$TL, "Difference (Ecopath minus SIA)" = g$Difference,
                        check.names = FALSE)
      if (!is.na(g$N)) out[["SIA sample size"]] <- as.integer(g$N)
      out[["Outcome"]] <- as.character(g$Status)
      out
    }, na = "none", digits = 2)
    output$dl_plot <- downloadHandler(
      filename = function() paste0(gsub("[^[:alnum:]._ -]", ".", current()), ".png"),
      content = function(file) ggsave(file, et_plot_detail(dat(), current(), post()$data),
                                      width = 9, height = 5.8, dpi = 200, bg = "white"))

    # ---- Sample Size tab
    output$samplesize_body <- renderUI({
      if (!has_data()) return(empty())
      s <- ss()
      if (is.null(s)) {
        return(tags$div(class = "et-empty", "Sample size data were not provided, or there are too few groups to fit a curve. Add SIA sample size as a 6th column of the input data to use this tab."))
      }
      tagList(
        tags$p(class = "et-lead", paste0("Each point is a functional group. The dotted line is a power function decay curve (CI range = ",
                                         round(s$a, 2), " x n^", round(s$b, 2),
                                         ") fitted to all groups. Groups are flagged when their credibility interval range deviates from the curve by more than ",
                                         s$flag_cutoff, " median absolute deviations on the log scale.")),
        card(plotOutput("samplesize_plot", height = "480px")),
        card(card_header("Functional groups that deviate from the curve"), tableOutput("samplesize_table")))
    })
    output$samplesize_plot <- renderPlot({ req(ss()); samplesizePlot(ss()) }, res = 96)
    output$samplesize_table <- renderTable({ req(ss()); samplesizeTable(ss()) })
  }
}

# -------------------------------------------------------------- main function

ecotuneR <- function(combined_output = NULL, export = FALSE, bypass_check = FALSE,
                     posterior = NULL, example = FALSE) {
  # combined_output: data frame in the 5 (or 6) column format; leave empty to upload in the app
  # export:          TRUE writes a report and plots to a folder instead of opening the app
  # bypass_check:    kept so older scripts still run; the checks no longer need a dialog
  # posterior:       optional data frame of posterior draws (consumer, trophic position)
  # example:         TRUE opens the app with the Barataria Bay example loaded
  data_name <- paste(deparse(substitute(combined_output)), collapse = "")
  if (is.null(combined_output)) {
    if (isTRUE(export)) stop("Give ecotuneR() a data frame to export results.", call. = FALSE)
    return(shinyApp(ui = et_ui(et_www()), server = et_server(start_example = isTRUE(example))))
  }
  v <- et_validate(combined_output)
  for (w in v$warnings) message("Note: ", w)
  if (!v$ok) {
    stop("'", data_name, "' is not ready for ecotuneR:\n  ", paste(v$errors, collapse = "\n  "),
         "\nColumns must be, in order: consumer, lower bound, upper bound, SIA mode, Ecopath trophic level, and optionally SIA sample size.",
         call. = FALSE)
  }
  d <- et_classify(v$data)
  p <- et_validate_posterior(posterior, d$Consumer)
  if (!is.null(p$message)) message(p$message)
  if (isTRUE(export)) return(et_export_folder(d, p$data, data_name))
  shinyApp(ui = et_ui(et_www()), server = et_server(init = combined_output, init_posterior = posterior))
}
