# Race Course Analyzer ---------------------------------------------------------
# Upload a GPX file to get:
#   1. An elevation profile (distance vs. elevation, y-axis not starting at 0)
#   2. Elevation gain/loss for each mile (final partial mile labeled, e.g. "3.1")
#   3. Where each mile marker falls on the course (1/2, 1, 2, ...)
#   4. Flat-equivalent distance from Minetti's energy cost of running vs. grade
#
# Packages: install.packages(c("shiny", "xml2", "ggplot2", "leaflet"))
# Run:      shiny::runApp("app.R")

library(shiny)
library(xml2)
library(ggplot2)
library(leaflet)

`%||%` <- function(a, b) if (is.null(a)) b else a

M_TO_FT  <- 3.28084
M_PER_MI <- 1609.344
FT_PER_MI <- 5280

# ---- Helpers -----------------------------------------------------------------

haversine_m <- function(lat1, lon1, lat2, lon2) {
  r <- 6371008.8
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}

parse_gpx <- function(path) {
  doc <- read_xml(path)
  xml_ns_strip(doc)
  pts <- xml_find_all(doc, "//trkpt")
  if (length(pts) == 0) pts <- xml_find_all(doc, "//rtept")
  if (length(pts) < 2) stop("No track or route points found in this GPX file.")

  lat <- as.numeric(xml_attr(pts, "lat"))
  lon <- as.numeric(xml_attr(pts, "lon"))
  ele <- suppressWarnings(as.numeric(xml_text(xml_find_first(pts, "./ele"))))

  if (all(is.na(ele))) {
    stop("This GPX file has no elevation (<ele>) data. ",
         "Try re-exporting it from Strava/Garmin with elevation included.")
  }
  # Fill any gaps in elevation by linear interpolation
  if (anyNA(ele)) {
    ok <- !is.na(ele)
    ele <- approx(which(ok), ele[ok], xout = seq_along(ele), rule = 2)$y
  }

  step_m <- c(0, haversine_m(head(lat, -1), head(lon, -1), lat[-1], lon[-1]))
  data.frame(
    lat    = lat,
    lon    = lon,
    dist   = cumsum(step_m) / M_PER_MI,  # miles
    ele_ft = ele * M_TO_FT               # GPX elevation is in meters
  )
}

# Moving average over a distance window (not a point count), so smoothing is
# consistent no matter how densely the GPS recorded points.
smooth_by_distance <- function(d, y, window_mi) {
  if (window_mi <= 0) return(y)
  h  <- window_mi / 2
  cs <- c(0, cumsum(y))
  lo <- findInterval(d - h, d, left.open = TRUE) + 1
  hi <- findInterval(d + h, d)
  (cs[hi + 1] - cs[lo]) / (hi - lo + 1)
}

interp_at <- function(x, y, xout) {
  approx(x, y, xout = xout, ties = mean, rule = 2)$y
}

gain_loss <- function(e) {
  de <- diff(e)
  c(gain = sum(de[de > 0]), loss = abs(sum(de[de < 0])))
}

# Minetti et al. (2002) energy cost of running, J/kg/m, for grade i (rise/run,
# decimal). The polynomial was fit over -0.45..+0.45, so grades are clamped.
minetti_cr <- function(i) {
  i <- pmin(pmax(i, -0.45), 0.45)
  155.4 * i^5 - 30.4 * i^4 - 43.3 * i^3 + 46.3 * i^2 + 19.5 * i + 3.6
}
MINETTI_CR0 <- minetti_cr(0)  # 3.6 J/kg/m on the flat

# Split the course into equal chunks of about `seg_mi`, take each chunk's grade
# from its endpoint elevations, and weight its length by Cr(i) / Cr(0).
# Returns one row per chunk with a running total of flat-equivalent miles.
flat_equivalent <- function(d, e, seg_mi) {
  total  <- max(d)
  n      <- max(1, round(total / seg_mi))
  br     <- seq(0, total, length.out = n + 1)
  z      <- interp_at(d, e, br)                  # ft
  run_ft <- diff(br) * FT_PER_MI
  rise   <- diff(z)
  grade  <- rise / run_ft
  len_mi <- sqrt(run_ft^2 + rise^2) / FT_PER_MI  # distance along the slope
  equiv  <- len_mi * minetti_cr(grade) / MINETTI_CR0
  data.frame(start = head(br, -1), end = br[-1], grade = grade,
             equiv = equiv, cum = cumsum(equiv))
}

# Flat-equivalent miles from the start to each distance in `x`
flat_equiv_at <- function(fe, x) approx(c(0, fe$end), c(0, fe$cum), xout = x, rule = 2)$y

mile_splits <- function(d, e, fe = NULL) {
  total <- max(d)
  whole <- floor(total)
  breaks <- c(0, seq_len(whole))
  # Add the final partial mile unless it's a negligible sliver (< 0.02 mi)
  if (total - whole >= 0.02 || whole == 0) breaks <- c(breaks, total)
  else breaks[length(breaks)] <- total

  rows <- lapply(seq_len(length(breaks) - 1), function(i) {
    a <- breaks[i]; b <- breaks[i + 1]
    inside <- d > a & d < b
    ys <- c(interp_at(d, e, a), e[inside], interp_at(d, e, b))
    gl <- gain_loss(ys)
    is_last_partial <- i == length(breaks) - 1 && abs(b - round(b)) > 1e-9
    label <- if (is_last_partial) sprintf("%.1f", b) else sprintf("%d", round(b))
    net <- ys[length(ys)] - ys[1]
    data.frame(
      Mile         = label,
      Segment      = sprintf("%.2f \u2013 %.2f mi", a, b),
      `Gain (ft)`  = sprintf("+%.0f", gl[["gain"]]),
      `Loss (ft)`  = sprintf("\u2212%.0f", gl[["loss"]]),
      `Net (ft)`   = sprintf("%+.0f", net),
      `Avg grade`  = sprintf("%+.1f%%", 100 * net / ((b - a) * FT_PER_MI)),
      `Flat-eq (mi)` = if (is.null(fe)) NA else
        sprintf("%.3f", diff(flat_equiv_at(fe, c(a, b)))),
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(fe)) out$`Flat-eq (mi)` <- NULL
  out
}

KM_TO_MI <- 0.621371

fmt_num <- function(x) sub("\\.?0+$", "", sprintf("%.2f", x))

# Marks are kept as a small table of (value, unit); unit is "mi" or "km".
empty_marks <- function() data.frame(value = numeric(), unit = character())

marks_to_mi <- function(m) m$value * ifelse(m$unit == "km", KM_TO_MI, 1)

mark_label <- function(value, unit) {
  ifelse(unit == "mi" & abs(value - 0.5) < 1e-9, "\u00bd mi",
         paste(sapply(value, fmt_num), unit))
}

mark_key <- function(value, unit) paste0(unit, "_", sapply(value, fmt_num))

preset_marks <- function(kind, total) {
  switch(kind,
    mile = data.frame(value = c(0.5, seq_len(floor(total))), unit = "mi"),
    half = data.frame(value = seq(0.5, total, by = 0.5), unit = "mi"),
    km   = data.frame(value = seq_len(floor(total / KM_TO_MI)), unit = "km"),
    empty_marks())
}

locate_marks <- function(trk, e, marks) {
  if (nrow(marks) == 0) return(NULL)
  marks$lat    <- interp_at(trk$dist, trk$lat, marks$mile)
  marks$lon    <- interp_at(trk$dist, trk$lon, marks$mile)
  marks$ele_ft <- interp_at(trk$dist, e, marks$mile)
  marks$change <- marks$ele_ft - c(e[1], head(marks$ele_ft, -1))
  marks
}

# Nearest point on the track to a clicked lat/lng: returns its distance (mi),
# or NA if the click was more than `max_m` meters from the course.
nearest_dist <- function(trk, lat, lng, max_m = 150) {
  m <- haversine_m(lat, lng, trk$lat, trk$lon)
  i <- which.min(m)
  if (m[i] > max_m) NA_real_ else trk$dist[i]
}

# Keep only the part of the track between `lo` and `hi` miles, with exact
# interpolated endpoints, and re-zero distance so the new start is mile 0.
trim_track <- function(trk, lo, hi) {
  ends   <- as.data.frame(lapply(trk, function(col) interp_at(trk$dist, col, c(lo, hi))))
  inside <- trk$dist > lo & trk$dist < hi
  out    <- rbind(ends[1, ], trk[inside, , drop = FALSE], ends[2, ])
  out$dist <- out$dist - lo
  rownames(out) <- NULL
  out
}

# "+12", "-8", "0" style labels for heights relative to the start
signed_ft <- function(x) ifelse(is.na(x), NA,
                                ifelse(round(x) == 0, "0", sprintf("%+.0f", x)))

# Elevation profile chart, shared by the app and the PNG export.
# relative = TRUE plots feet above/below the start instead of elevation.
build_profile <- function(trk, e, mk, show_raw = FALSE, relative = FALSE,
                          title = NULL, subtitle = NULL) {
  base <- if (relative) e[1] else 0
  df <- data.frame(dist = trk$dist, ele = e - base, raw = trk$ele_ft - base)
  if (!is.null(mk)) mk$ele_ft <- mk$ele_ft - base

  yr  <- range(if (show_raw) c(df$ele, df$raw) else df$ele)
  pad <- max(10, diff(yr) * 0.12)
  ylo <- yr[1] - pad
  yhi <- yr[2] + pad

  p <- ggplot(df, aes(dist)) +
    geom_ribbon(aes(ymin = ylo, ymax = ele), fill = "#3b82f6", alpha = 0.18) +
    geom_line(aes(y = ele), colour = "#1d4ed8", linewidth = 0.9)

  if (show_raw) {
    p <- p + geom_line(aes(y = raw), colour = "grey50", linewidth = 0.3, alpha = 0.7)
  }
  if (relative) {
    p <- p + geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
      annotate("text", x = 0, y = 0, label = " start", hjust = 0, vjust = -0.5,
               size = 3.2, colour = "grey30")
  }
  if (!is.null(mk)) {
    p <- p +
      geom_vline(xintercept = mk$mile, linetype = "dotted", colour = "grey55") +
      geom_point(data = mk, aes(mile, ele_ft), colour = mk$colour, size = 2.2) +
      geom_text(data = mk, aes(x = mile, y = yhi, label = label,
                               hjust = ifelse(mile > max(df$dist) * 0.97, 1.05, 0.5)),
                size = 3.2, vjust = 1.2, check_overlap = TRUE)
  }

  p +
    scale_x_continuous(breaks = c(0, 0.5, seq_len(floor(max(df$dist)))),
                       labels = function(x) ifelse(x == 0.5, "\u00bd", x),
                       expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = if (relative) signed_ft else waiver()) +
    coord_cartesian(ylim = c(ylo, yhi), expand = FALSE) +
    labs(x = "Distance (miles)",
         y = if (relative) "Elevation vs. start (ft)" else "Elevation (ft)",
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey35"))
}

# Write the course line, start/finish and marks as KML (Google Earth, Google
# My Maps, CalTopo, Gaia, etc.). KML colors are aabbggrr hex.
write_course_kml <- function(file, trk, e, mk, name, trimmed = FALSE) {
  hex_to_kml <- function(hex, alpha = "ff") {
    h <- sub("^#", "", hex)
    paste0(alpha, substr(h, 5, 6), substr(h, 3, 4), substr(h, 1, 2))
  }
  icon_href <- "http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png"

  doc  <- xml_new_root("kml", xmlns = "http://www.opengis.net/kml/2.2")
  d    <- xml_add_child(doc, "Document")
  xml_add_child(d, "name", name)
  gl <- gain_loss(e)
  xml_add_child(d, "description", sprintf(
    "%.2f mi, +%.0f ft / -%.0f ft%s. Exported from Race Course Analyzer.",
    max(trk$dist), gl[["gain"]], gl[["loss"]], if (trimmed) " (trimmed)" else ""))

  add_point_style <- function(id, hex, scale = 1) {
    st <- xml_add_child(d, "Style", id = id)
    is <- xml_add_child(st, "IconStyle")
    xml_add_child(is, "color", hex_to_kml(hex))
    xml_add_child(is, "scale", scale)
    xml_add_child(xml_add_child(is, "Icon"), "href", icon_href)
    xml_add_child(xml_add_child(st, "LabelStyle"), "scale", 0.9)
  }
  ls <- xml_add_child(xml_add_child(d, "Style", id = "course"), "LineStyle")
  xml_add_child(ls, "color", hex_to_kml("#1d4ed8"))
  xml_add_child(ls, "width", 4)
  add_point_style("start",  "#16a34a", 1.2)
  add_point_style("finish", "#111827", 1.2)
  add_point_style("mark_mi", "#dc2626")
  add_point_style("mark_km", "#4338ca")

  coord <- function(lon, lat, ele_ft) sprintf("%.7f,%.7f,%.1f", lon, lat, ele_ft / M_TO_FT)

  add_point <- function(parent, nm, style, lon, lat, ele_ft, desc = NULL) {
    pm <- xml_add_child(parent, "Placemark")
    xml_add_child(pm, "name", nm)
    if (!is.null(desc)) xml_add_child(pm, "description", desc)
    xml_add_child(pm, "styleUrl", paste0("#", style))
    pt <- xml_add_child(pm, "Point")
    xml_add_child(pt, "altitudeMode", "clampToGround")
    xml_add_child(pt, "coordinates", coord(lon, lat, ele_ft))
  }

  # Course
  f  <- xml_add_child(d, "Folder"); xml_add_child(f, "name", "Course")
  pm <- xml_add_child(f, "Placemark")
  xml_add_child(pm, "name", name)
  xml_add_child(pm, "styleUrl", "#course")
  line <- xml_add_child(pm, "LineString")
  xml_add_child(line, "tessellate", 1)
  xml_add_child(line, "altitudeMode", "clampToGround")
  xml_add_child(line, "coordinates",
                paste(coord(trk$lon, trk$lat, e), collapse = " "))
  n <- nrow(trk)
  add_point(f, "Start",  "start",  trk$lon[1], trk$lat[1], e[1],
            sprintf("Elevation %.0f ft", e[1]))
  add_point(f, "Finish", "finish", trk$lon[n], trk$lat[n], e[n],
            sprintf("%.2f mi, elevation %.0f ft", max(trk$dist), e[n]))

  # Marks
  if (!is.null(mk) && nrow(mk) > 0) {
    f <- xml_add_child(d, "Folder"); xml_add_child(f, "name", "Marks")
    for (i in seq_len(nrow(mk))) {
      add_point(f, mk$label[i],
                if (grepl("km$", mk$label[i])) "mark_km" else "mark_mi",
                mk$lon[i], mk$lat[i], mk$ele_ft[i],
                sprintf("%.2f mi (%.2f km). Elevation %.0f ft (%s ft vs. start, %+.0f ft since previous mark).",
                        mk$mile[i], mk$mile[i] / KM_TO_MI, mk$ele_ft[i],
                        signed_ft(mk$ele_ft[i] - e[1]), mk$change[i]))
    }
  }
  write_xml(doc, file, options = c("format", "as_xml"))
}

stat_card <- function(label, value, note = NULL) {
  div(class = "stat",
      div(class = "stat-value", value),
      div(class = "stat-label", label),
      if (!is.null(note)) div(class = "stat-note", note))
}

# ---- UI ----------------------------------------------------------------------

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background: #f7f7f5; }
    .stats { display: flex; flex-wrap: wrap; gap: 12px; margin: 8px 0 16px; }
    .stat { background: #fff; border: 1px solid #e3e3e0; border-radius: 8px;
            padding: 10px 16px; min-width: 130px; }
    .stat-value { font-size: 22px; font-weight: 600; color: #1f2d3d; }
    .stat-label { font-size: 12px; color: #6b7280; text-transform: uppercase;
                  letter-spacing: .04em; }
    .stat-note { font-size: 12px; color: #4b5563; margin-top: 2px; }
    .panel-box { background: #fff; border: 1px solid #e3e3e0; border-radius: 8px;
                 padding: 12px 16px; margin-bottom: 16px; }
    .panel-box h4 { margin-top: 4px; }
    .panel-head { display: flex; justify-content: space-between; align-items: center;
      gap: 12px; }
    .panel-head .btn { padding: 4px 10px; font-size: 13px; }
    .head-tools { display: flex; align-items: center; gap: 10px; }
    .seg .form-group, .seg .shiny-options-group { margin: 0; }
    .seg .shiny-options-group { display: inline-flex; border: 1px solid #ccc;
      border-radius: 4px; overflow: hidden; }
    .seg .radio-inline { margin: 0 !important; padding: 0; }
    .seg input[type=radio] { position: absolute; opacity: 0; pointer-events: none; }
    .seg .radio-inline span { display: block; padding: 4px 10px; font-size: 13px;
      cursor: pointer; color: #374151; background: #fff; }
    .seg .radio-inline + .radio-inline span { border-left: 1px solid #ccc; }
    .seg input[type=radio]:checked + span { background: #1d4ed8; color: #fff; }
    .seg input[type=radio]:focus-visible + span { outline: 2px solid #93c5fd; }
    .marks-label { font-weight: 700; margin-bottom: 6px; display: block; }
    .add-row { display: flex; gap: 6px; align-items: stretch; }
    .add-row input { flex: 1; min-width: 0; height: 34px; padding: 4px 8px;
      border: 1px solid #ccc; border-radius: 4px; font-size: 14px; }
    .unit-toggle { display: inline-flex; border: 1px solid #ccc; border-radius: 4px;
      overflow: hidden; }
    .unit-toggle button { border: 0; background: #fff; padding: 0 10px;
      font-size: 13px; color: #374151; cursor: pointer; }
    .unit-toggle button.active { background: #1d4ed8; color: #fff; }
    .add-btn { height: 34px; padding: 0 12px; border: 0; border-radius: 4px;
      background: #1d4ed8; color: #fff; font-weight: 600; cursor: pointer; }
    .add-btn:hover { background: #1e40af; }
    .mark-msg { min-height: 1.3em; font-size: 13px; color: #b45309; margin: 4px 0; }
    .chips { display: flex; flex-wrap: wrap; gap: 6px; margin: 4px 0 8px; }
    .chip { display: inline-flex; align-items: center; gap: 2px; padding: 2px 4px 2px 10px;
      border-radius: 999px; background: #fee2e2; color: #991b1b; font-size: 13px;
      font-variant-numeric: tabular-nums; }
    .chip.km { background: #e0e7ff; color: #3730a3; }
    .chip.out { background: #f3f4f6; color: #9ca3af; text-decoration: line-through; }
    .chip button { border: 0; background: transparent; color: inherit; cursor: pointer;
      font-size: 15px; line-height: 1; padding: 0 5px; border-radius: 999px; }
    .chip button:hover { background: rgba(0,0,0,.08); }
    .chips-empty { color: #9ca3af; font-size: 13px; }
    .quick { font-size: 13px; color: #6b7280; }
    .quick a { margin-right: 8px; white-space: nowrap; }
    .profile-wrap { position: relative; }
    .profile-wrap .cursor-line { position: absolute; width: 0; display: none;
      border-left: 1.5px solid #f59e0b; pointer-events: none; z-index: 5; }
    .profile-wrap .cursor-dot { position: absolute; width: 10px; height: 10px;
      margin: -5px 0 0 -5px; border-radius: 50%; background: #f59e0b;
      border: 2px solid #fff; display: none; pointer-events: none; z-index: 6; }
    .hover-info { color: #4b5563; min-height: 1.4em; font-variant-numeric: tabular-nums; }
    .hint { color: #6b7280; font-size: 13px; }
    .trim-box .form-group { margin-bottom: 0; }
    .trim-note { font-size: 13px; color: #4b5563; margin: -4px 0 4px;
      font-variant-numeric: tabular-nums; }
    .trim-note.on { color: #b45309; }
    table { width: 100%; }
  "))),

  titlePanel("Race Course Analyzer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("gpx", "Upload a GPX file", accept = c(".gpx")),
      uiOutput("trim_ui"),

      tags$label(class = "marks-label", `for` = "new-mark", "Marks"),
      div(class = "add-row",
          tags$input(id = "new-mark", type = "number", min = "0", step = "0.25",
                     placeholder = "e.g. 1.5", inputmode = "decimal"),
          div(class = "unit-toggle",
              tags$button(type = "button", class = "active", `data-unit` = "mi", "mi"),
              tags$button(type = "button", `data-unit` = "km", "km")),
          tags$button(type = "button", class = "add-btn", id = "add-mark-btn", "Add")),
      uiOutput("mark_msg"),
      uiOutput("mark_chips"),
      div(class = "quick", "Quick add: ",
          actionLink("preset_mile", "\u00bd + every mile"),
          actionLink("preset_half", "every \u00bd mi"),
          actionLink("preset_km",   "every km"),
          actionLink("clear_marks", "clear all")),
      helpText("Or click the course on the map, or the elevation chart,",
               "to add a mark at that spot."),
      hr(),
      sliderInput("smooth_ft", "Elevation smoothing window (ft)",
                  min = 0, max = 500, value = 100, step = 25),
      helpText("GPS elevation is noisy; raw data overstates climbing.",
               "A ~100 ft window usually lands close to official course",
               "numbers. Set to 0 to use raw elevations."),
      checkboxInput("show_raw", "Show raw elevation on chart", FALSE),
      hr(),
      sliderInput("grade_seg_ft", "Flat-equivalent segment length (ft)",
                  min = 25, max = 1320, value = 200, step = 25),
      helpText("Flat-equivalent distance uses Minetti's energy cost of running",
               "vs. grade (the physiology behind Strava's GAP). The course is",
               "cut into equal segments of about this length; each counts as",
               "its length \u00d7 Cr(grade) / Cr(0). Shorter segments catch",
               "short, steep hills but are more sensitive to GPS noise.")
    ),

    mainPanel(
      width = 9,
      uiOutput("summary"),

      div(class = "panel-box",
          div(class = "panel-head",
              h4("Elevation profile"),
              div(class = "head-tools",
                  div(class = "seg",
                      radioButtons("y_mode", NULL, inline = TRUE,
                                   choices = c("Absolute" = "abs",
                                               "Relative" = "rel"))),
                  uiOutput("download_ui"))),
          div(id = "profile-wrap", class = "profile-wrap",
              plotOutput("profile", height = "360px",
                         hover = hoverOpts("profile_hover", delay = 40,
                                           delayType = "throttle", clip = TRUE),
                         click = "profile_click"),
              div(class = "cursor-line"),
              div(class = "cursor-dot")),
          div(class = "hover-info", textOutput("hover_info"))),

      div(class = "panel-box",
          div(class = "panel-head",
              h4("Course map"),
              uiOutput("download_kml_ui")),
          div(class = "hint", "Click anywhere on the course to drop a mark.",
              "Hover the elevation chart to see that spot here."),
          leafletOutput("map", height = "500px")),

      fluidRow(
        column(6, div(class = "panel-box",
                      h4("Elevation by mile"),
                      tableOutput("splits"))),
        column(6, div(class = "panel-box",
                      h4("Marks"),
                      tableOutput("markers")))
      )
    )
  ),

  # Cursor line on the elevation chart. It follows the mouse in the browser;
  # the server only supplies the plot panel's edges and the curve's height.
  tags$script(HTML("
    (function() {
      var wrap = document.getElementById('profile-wrap');
      var line = wrap.querySelector('.cursor-line');
      var dot  = wrap.querySelector('.cursor-dot');
      var b = null, off = 0;
      function hide() { line.style.display = 'none'; dot.style.display = 'none'; }
      Shiny.addCustomMessageHandler('profile-cursor', function(m) {
        b = m;
        if (m.dot == null) { dot.style.display = 'none'; }
        else if (line.style.display === 'block') {
          dot.style.top = (off + m.dot) + 'px';
          dot.style.display = 'block';
        }
      });
      wrap.addEventListener('mousemove', function(e) {
        var img = wrap.querySelector('img');
        if (!b || !img) return;
        var r = img.getBoundingClientRect(), w = wrap.getBoundingClientRect();
        var x = e.clientX - r.left;
        if (x < b.left || x > b.right) { hide(); return; }
        off = r.top - w.top;
        line.style.left   = (r.left - w.left + x) + 'px';
        line.style.top    = (off + b.top) + 'px';
        line.style.height = (b.bottom - b.top) + 'px';
        line.style.display = 'block';
        dot.style.left = line.style.left;
      });
      wrap.addEventListener('mouseleave', hide);
    })();

    // Mark editor: send the typed value with the event so it is never stale
    (function() {
      function unit() {
        var b = document.querySelector('.unit-toggle button.active');
        return b ? b.getAttribute('data-unit') : 'mi';
      }
      function send() {
        var el = document.getElementById('new-mark');
        var v = parseFloat(el.value);
        Shiny.setInputValue('add_mark_req',
          { value: isNaN(v) ? null : v, unit: unit(), raw: el.value },
          { priority: 'event' });
        el.value = ''; el.focus();
      }
      document.addEventListener('click', function(e) {
        var t = e.target.closest('.unit-toggle button');
        if (t) {
          t.parentNode.querySelectorAll('button').forEach(function(b) {
            b.classList.toggle('active', b === t);
          });
          document.getElementById('new-mark').focus();
          return;
        }
        if (e.target.closest('#add-mark-btn')) { send(); return; }
        var x = e.target.closest('.chip button');
        if (x) Shiny.setInputValue('remove_mark', x.getAttribute('data-key'),
                                   { priority: 'event' });
      });
      document.addEventListener('keydown', function(e) {
        if (e.target.id === 'new-mark' && e.key === 'Enter') { e.preventDefault(); send(); }
      });
    })();
  "))
)

# ---- Server ------------------------------------------------------------------

server <- function(input, output, session) {

  # The whole recording, exactly as uploaded
  raw_track <- reactive({
    validate(need(input$gpx, "Upload a GPX file to get started."))
    tryCatch(parse_gpx(input$gpx$datapath),
             error = function(e) validate(need(FALSE, conditionMessage(e))))
  })
  raw_total <- reactive(max(raw_track()$dist))

  # ---- Trim: drop a warm-up before the start or a late watch stop ----
  # The chosen range is stored with the file it belongs to, so a newly
  # uploaded course always starts untrimmed.
  trim_state <- reactiveVal(NULL)
  trim_in    <- debounce(reactive(input$trim), 250)
  observeEvent(trim_in(), {
    trim_state(list(path = isolate(input$gpx$datapath), range = trim_in()))
  })

  trim_range <- reactive({
    total <- raw_total()
    t <- trim_state()
    if (is.null(t) || !identical(t$path, input$gpx$datapath)) return(c(0, total))
    c(max(0, t$range[1]), min(total, t$range[2]))
  })
  is_trimmed <- reactive({
    r <- trim_range()
    r[1] > 0.005 || r[2] < raw_total() - 0.005
  })

  slider_max <- reactive(ceiling(raw_total() * 100) / 100)

  output$trim_ui <- renderUI({
    req(input$gpx)
    req(tryCatch(raw_track(), error = function(e) NULL))
    div(class = "trim-box",
        sliderInput("trim", "Trim course (mi)", min = 0, max = slider_max(),
                    value = c(0, slider_max()), step = 0.01, ticks = FALSE,
                    width = "100%"),
        uiOutput("trim_note"))
  })

  output$trim_note <- renderUI({
    r <- trim_range()
    if (!is_trimmed()) {
      return(div(class = "trim-note",
                 "Drag either end to cut off a warm-up or a late watch stop."))
    }
    div(class = "trim-note on",
        sprintf("Using %.2f\u2013%.2f mi (%.2f of %.2f mi). ",
                r[1], r[2], diff(r), raw_total()),
        actionLink("trim_reset", "Reset"))
  })

  observeEvent(input$trim_reset, {
    updateSliderInput(session, "trim", value = c(0, slider_max()))
  })

  # Smooth the full recording first, so the trimmed ends don't suffer edge
  # effects, then cut. Everything downstream works on the trimmed course,
  # with distance re-zeroed at the new start.
  track <- reactive({
    trk <- raw_track()
    trk$ele_s <- smooth_by_distance(trk$dist, trk$ele_ft, input$smooth_ft / FT_PER_MI)
    r <- trim_range()
    validate(need(diff(r) >= 0.05, "The trimmed course is too short; widen the trim range."))
    trim_track(trk, r[1], r[2])
  })

  ele <- reactive(track()$ele_s)

  flat_eq <- reactive({
    trk <- track()
    flat_equivalent(trk$dist, ele(), input$grade_seg_ft / FT_PER_MI)
  })

  # ---- Marks: one list, edited by the add row, presets, chips and clicks ----

  marks      <- reactiveVal(preset_marks("mile", 3))
  undo_marks <- reactiveVal(NULL)
  mark_msg   <- reactiveVal(NULL)

  course_total <- reactive({
    if (is.null(input$gpx)) return(NA_real_)
    tryCatch(max(track()$dist), error = function(e) NA_real_)
  })

  # Add marks, skipping duplicates and anything past the finish.
  # Returns a message for a single rejected mark, or NULL.
  add_marks <- function(new, quiet = FALSE) {
    cur   <- marks()
    total <- course_total()
    mi    <- marks_to_mi(new)
    dup   <- round(mi, 3) %in% round(marks_to_mi(cur), 3) | duplicated(round(mi, 3))
    past  <- !is.na(total) & mi > total + 1e-6
    ok    <- !dup & !past
    marks(rbind(cur, new[ok, , drop = FALSE]))
    if (quiet || nrow(new) != 1 || ok) return(NULL)
    lab <- mark_label(new$value, new$unit)
    if (past) sprintf("%s is past the finish (%.2f mi).", lab, total)
    else sprintf("%s is already on the list.", lab)
  }

  observeEvent(input$add_mark_req, {
    r <- input$add_mark_req
    if (is.null(r$value) || !is.finite(r$value) || r$value <= 0) {
      mark_msg(if (nzchar(r$raw %||% "")) "Enter a distance greater than 0."
               else "Type a distance, then press Enter or Add.")
      return()
    }
    mark_msg(add_marks(data.frame(value = r$value, unit = r$unit)))
  })

  observeEvent(input$remove_mark, {
    m <- marks()
    marks(m[mark_key(m$value, m$unit) != input$remove_mark, , drop = FALSE])
    mark_msg(NULL)
  })

  # Clicks on the map or chart add a mile mark at that spot
  add_click <- function(d) {
    if (is.na(d) || d <= 0) return()
    mark_msg(add_marks(data.frame(value = round(d, 2), unit = "mi")))
  }
  observeEvent(input$map_click, {
    add_click(nearest_dist(track(), input$map_click$lat, input$map_click$lng))
  })
  observeEvent(input$profile_click, {
    x <- input$profile_click$x
    if (!is.null(x)) add_click(min(max(x, 0), max(track()$dist)))
  })

  quick_add <- function(kind) {
    total <- course_total()
    if (is.na(total)) total <- 3.1   # before a course is loaded, assume a 5K
    add_marks(preset_marks(kind, total), quiet = TRUE)
    mark_msg(NULL)
  }
  observeEvent(input$preset_mile, quick_add("mile"))
  observeEvent(input$preset_half, quick_add("half"))
  observeEvent(input$preset_km,   quick_add("km"))

  observeEvent(input$clear_marks, {
    if (nrow(marks()) == 0) return()
    undo_marks(marks())
    marks(empty_marks())
    mark_msg("cleared")
  })
  observeEvent(input$undo_clear, {
    if (!is.null(undo_marks())) marks(undo_marks())
    undo_marks(NULL); mark_msg(NULL)
  })

  output$mark_msg <- renderUI({
    msg <- mark_msg()
    if (identical(msg, "cleared")) {
      div(class = "mark-msg", "All marks cleared. ", actionLink("undo_clear", "Undo"))
    } else div(class = "mark-msg", msg)
  })

  output$mark_chips <- renderUI({
    m <- marks()
    if (nrow(m) == 0) return(div(class = "chips chips-empty", "No marks yet."))
    m$mi <- marks_to_mi(m)
    m <- m[order(m$mi), , drop = FALSE]
    total <- course_total()
    div(class = "chips", lapply(seq_len(nrow(m)), function(i) {
      out <- !is.na(total) && m$mi[i] > total + 1e-6
      span(class = paste("chip", if (m$unit[i] == "km") "km", if (out) "out"),
           title = if (out) "Past the finish of this course" else NULL,
           mark_label(m$value[i], m$unit[i]),
           tags$button(type = "button", `data-key` = mark_key(m$value[i], m$unit[i]),
                       title = "Remove", "\u00d7"))
    }))
  })

  markers <- reactive({
    trk <- track()
    m <- marks()
    df <- data.frame(mile = marks_to_mi(m), label = mark_label(m$value, m$unit),
                     colour = ifelse(m$unit == "km", "#4338ca", "#dc2626"))
    df <- df[df$mile <= max(trk$dist) + 1e-6, , drop = FALSE]
    locate_marks(trk, ele(), df[order(df$mile), , drop = FALSE])
  })

  output$summary <- renderUI({
    trk <- track(); e <- ele(); gl <- gain_loss(e)
    div(class = "stats",
        stat_card("Distance",       sprintf("%.2f mi", max(trk$dist))),
        stat_card("Elevation gain", sprintf("%.0f ft", gl[["gain"]])),
        stat_card("Elevation loss", sprintf("%.0f ft", gl[["loss"]])),
        {
          fe <- sum(flat_eq()$equiv); actual <- max(trk$dist)
          stat_card("Flat-equivalent", sprintf("%.2f mi", fe),
                    sprintf("%+.1f%% vs. actual", 100 * (fe / actual - 1)))
        },
        if (relative()) {
          tagList(stat_card("Low vs. start",  paste(signed_ft(min(e) - e[1]), "ft")),
                  stat_card("High vs. start", paste(signed_ft(max(e) - e[1]), "ft")))
        } else {
          tagList(stat_card("Low point",  sprintf("%.0f ft", min(e))),
                  stat_card("High point", sprintf("%.0f ft", max(e))))
        })
  })

  relative <- reactive(identical(input$y_mode, "rel"))
  base_ft  <- reactive(if (relative()) ele()[1] else 0)

  output$profile <- renderPlot({
    build_profile(track(), ele(), markers(), input$show_raw, relative())
  })

  course_name <- reactive({
    req(input$gpx)
    tools::file_path_sans_ext(input$gpx$name)
  })

  output$download_ui <- renderUI({
    req(input$gpx)
    downloadButton("download_profile", "Download PNG", class = "btn-default btn-sm")
  })

  output$download_kml_ui <- renderUI({
    req(input$gpx)
    downloadButton("download_kml", "Download KML", class = "btn-default btn-sm")
  })

  output$download_kml <- downloadHandler(
    filename = function() paste0(course_name(), "-course.kml"),
    content = function(file) {
      write_course_kml(file, track(), ele(), markers(), course_name(),
                       trimmed = is_trimmed())
    },
    contentType = "application/vnd.google-earth.kml+xml"
  )

  output$download_profile <- downloadHandler(
    filename = function() paste0(course_name(), "-elevation.png"),
    content = function(file) {
      trk <- track(); e <- ele(); gl <- gain_loss(e)
      range_txt <- if (relative()) {
        sprintf("%s to %s ft vs. start", signed_ft(min(e) - e[1]), signed_ft(max(e) - e[1]))
      } else sprintf("%.0f\u2013%.0f ft", min(e), max(e))
      p <- build_profile(
        trk, e, markers(), input$show_raw, relative(),
        title    = course_name(),
        subtitle = sprintf("%.2f mi (%.2f flat-eq)  \u00b7  +%.0f ft / \u2212%.0f ft  \u00b7  %s",
                           max(trk$dist), sum(flat_eq()$equiv),
                           gl[["gain"]], gl[["loss"]], range_txt)
      )
      ggsave(file, p, width = 10, height = 4.5, dpi = 200, bg = "white")
    }
  )

  hover_dist <- reactive({
    x <- input$profile_hover$x
    if (is.null(x)) NULL else min(max(x, 0), max(track()$dist))
  })

  # Send the plot panel's edges (in CSS pixels) and the profile's height at the
  # cursor, so the browser can draw the cursor line and dot.
  observe({
    h <- input$profile_hover
    req(h, h$range, h$domain)
    rx <- h$img_css_ratio$x; ry <- h$img_css_ratio$y
    if (is.null(rx)) rx <- h$coords_img$x / h$coords_css$x
    if (is.null(ry)) ry <- h$coords_img$y / h$coords_css$y
    y  <- interp_at(track()$dist, ele(), hover_dist()) - base_ft()
    frac <- (y - h$domain$bottom) / (h$domain$top - h$domain$bottom)
    dot_px <- h$range$bottom - frac * (h$range$bottom - h$range$top)
    session$sendCustomMessage("profile-cursor", list(
      left   = h$range$left / rx,  right  = h$range$right / rx,
      top    = h$range$top / ry,   bottom = h$range$bottom / ry,
      dot    = if (is.finite(dot_px)) dot_px / ry else NULL
    ))
  })

  output$hover_info <- renderText({
    d <- hover_dist()
    if (is.null(d)) return("Hover the chart to read distance and elevation; click to add a mark.")
    y <- interp_at(track()$dist, ele(), d)
    sprintf("%.2f mi (%.2f km)  \u00b7  %s", d, d / KM_TO_MI,
            if (relative()) paste(signed_ft(y - base_ft()), "ft vs. start")
            else sprintf("%.0f ft", y))
  })

  output$splits <- renderTable({
    trk <- track()
    mile_splits(trk$dist, ele(), flat_eq())
  }, striped = TRUE, align = "llrrrrr")

  output$markers <- renderTable({
    mk <- markers()
    validate(need(!is.null(mk), "No marks yet. Type distances or click the course."))
    out <- data.frame(
      Mark          = mk$label,
      `Dist (mi)`   = sprintf("%.2f", mk$mile),
      Elev          = if (relative()) signed_ft(mk$ele_ft - base_ft())
                      else sprintf("%.0f", mk$ele_ft),
      `Change (ft)` = sprintf("%+.0f", mk$change),
      Lat           = sprintf("%.5f", mk$lat),
      Lon           = sprintf("%.5f", mk$lon),
      check.names = FALSE
    )
    names(out)[3] <- if (relative()) "vs. start (ft)" else "Elev (ft)"
    out
  }, striped = TRUE, align = "lrrrrr")

  # ---- Map: drawn once per course; marks and hover dot update in place ----

  add_marks_layer <- function(m, mk) {
    m <- clearGroup(m, "marks")
    if (is.null(mk)) return(m)
    addCircleMarkers(m, lng = mk$lon, lat = mk$lat, group = "marks",
                     radius = 6, color = mk$colour, fillColor = "#fff",
                     weight = 3, fillOpacity = 1, label = mk$label,
                     labelOptions = labelOptions(noHide = TRUE, direction = "right",
                                                 textsize = "12px"))
  }

  # The (trimmed) course line with its start and finish dots
  add_course_layer <- function(m, trk) {
    m <- clearGroup(m, "course")
    n <- nrow(trk)
    m |>
      addPolylines(lng = trk$lon, lat = trk$lat, color = "#1d4ed8", group = "course",
                   weight = 5, opacity = 0.85) |>
      addCircleMarkers(lng = trk$lon[1], lat = trk$lat[1], radius = 7, group = "course",
                       color = "#16a34a", fillOpacity = 1, stroke = FALSE,
                       label = "Start") |>
      addCircleMarkers(lng = trk$lon[n], lat = trk$lat[n], radius = 7, group = "course",
                       color = "#111827", fillOpacity = 1, stroke = FALSE,
                       label = "Finish")
  }

  # Drawn once per upload, so trimming doesn't reset the zoom. The full
  # recording sits underneath as a grey dashed line, which is what shows up
  # in the trimmed-off ends.
  output$map <- renderLeaflet({
    raw <- raw_track()
    # Both basemaps are free and need no API key
    m <- leaflet() |>
      addTiles(group = "Map") |>
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") |>
      addLayersControl(baseGroups = c("Map", "Satellite"),
                       options = layersControlOptions(collapsed = FALSE)) |>
      addPolylines(lng = raw$lon, lat = raw$lat, color = "#6b7280", weight = 3,
                   opacity = 0.8, dashArray = "4 6", label = "Trimmed off",
                   options = pathOptions(interactive = FALSE))
    m <- add_course_layer(m, isolate(track()))
    add_marks_layer(m, isolate(markers()))
  })

  observe({
    m <- add_course_layer(leafletProxy("map"), track())
    add_marks_layer(m, markers())
  })

  observe({
    d <- hover_dist()
    proxy <- clearGroup(leafletProxy("map"), "hover")
    if (is.null(d)) return()
    trk <- isolate(track())
    addCircleMarkers(proxy, group = "hover",
                     lng = interp_at(trk$dist, trk$lon, d),
                     lat = interp_at(trk$dist, trk$lat, d),
                     radius = 8, color = "#f59e0b", fillColor = "#f59e0b",
                     fillOpacity = 0.9, weight = 2,
                     options = pathOptions(interactive = FALSE))
  })
}

shinyApp(ui, server)
