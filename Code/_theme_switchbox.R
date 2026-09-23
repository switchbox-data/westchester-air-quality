# Switchbox ggplot2 house style for the EJLPC air-quality project.
# Source from a script that produces figures: source("Code/_theme_switchbox.R")
#
# Ported/adapted from reports2/lib/ggplot/switchbox_theme.R. Fonts and logo are
# loaded from bundled local files (Code/assets/) — NO network download, so this
# works in the offline container as well as a bind-mounted Code/.
#
# Public API (see DELIVERABLE notes at the bottom of this file):
#   sb_colors                          named brand palette vector
#   theme_set(theme_sb())              activate the house theme
#   scale_fill_sb_map(name = ...)      sequential map fill ramp (white->sky->midnight)
#   scale_color_sb() / scale_fill_sb() categorical palette for multi-series charts
#   sb_basemap(bbox_sf, ...)           greyscale basemap raster data frame (hex col)
#   add_sb_logo(plot, ...)             overlay the Switchbox logo on a figure
#   save_sb_fig(plot, path, ...)       ggsave wrapper @ 300 dpi, white bg

suppressPackageStartupMessages({
  library(sysfonts)
  library(showtext)
  library(ggplot2)
})

# ---- Asset path resolution --------------------------------------------------
# Resolve the assets dir robustly: works when run as `Rscript Code/01_*.R` from
# /project (path "Code/assets") and from a few other plausible working dirs.
.sb_find_asset <- function(filename) {
  candidates <- c(
    file.path("Code", "assets", filename),
    file.path("assets", filename),
    file.path(getwd(), "Code", "assets", filename)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop(sprintf("Switchbox theme: asset '%s' not found. Looked in: %s",
                 filename, paste(candidates, collapse = ", ")))
  }
  hit[1]
}

# ---- Font registration (local files, no download) --------------------------
font_add(
  family  = "IBM-Plex-Sans",
  regular = .sb_find_asset("IBMPlexSans-Regular.otf"),
  bold    = .sb_find_asset("IBMPlexSans-Bold.otf")
)
showtext_auto()
showtext::showtext_opts(dpi = 300)

# ---- Brand palette ----------------------------------------------------------
sb_colors <- c(
  "sky"           = "#68BED8", # primary
  "midnight"      = "#023047", # primary
  "carrot"        = "#FC9706", # primary
  "saffron"       = "#FFC729", # secondary
  "pistachio"     = "#A0AF12", # secondary
  "black"         = "#000000", # utilitarian
  "white"         = "#FFFFFF", # utilitarian
  "midnight_text" = "#0B6082", # lighter, text only
  "pistachio_text" = "#546800" # darker, text only
)

# ---- Base theme -------------------------------------------------------------
# Returns a theme object so callers can `theme_set(theme_sb())` or `p + theme_sb()`.
theme_sb <- function() {
  theme_minimal() %+replace%
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      legend.title = element_text(family = "IBM-Plex-Sans", size = 12, hjust = 0.5),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      text = element_text(family = "IBM-Plex-Sans", size = 12),
      axis.text = element_text(family = "IBM-Plex-Sans", size = 12),
      axis.title = element_text(family = "IBM-Plex-Sans", size = 12),
      strip.text = element_text(family = "IBM-Plex-Sans", size = 12),
      axis.title.x = element_text(margin = margin(t = 3)),
      axis.title.y = element_text(margin = margin(r = 3), angle = 90)
    )
}

# ---- Map fill ramp (sequential: white -> sky -> midnight) -------------------
scale_fill_sb_map <- function(name = waiver(), ...) {
  scale_fill_gradientn(
    name    = name,
    colours = c("#FFFFFF", "#68BED8", "#023047"),
    ...
  )
}

# ---- Categorical palette (multi-series) -------------------------------------
.sb_cat_pal <- unname(sb_colors[c("midnight", "carrot", "sky", "pistachio")])

scale_color_sb <- function(...) {
  scale_color_manual(values = .sb_cat_pal, ...)
}
scale_colour_sb <- scale_color_sb

scale_fill_sb <- function(...) {
  scale_fill_manual(values = .sb_cat_pal, ...)
}

# ---- Greyscale basemap ------------------------------------------------------
# Mirrors the bm_df pattern in 02_Health_risk.R: returns a data frame with x, y
# and an `hex` column for geom_raster(... fill = hex) + scale_fill_identity().
# Uses a greyscale tile provider so colored tracts dominate. If tiles can't be
# fetched (no network in the container), returns NULL and warns; callers should
# fall back to a plain white panel.
sb_basemap <- function(bbox_sf, provider = "CartoDB.PositronNoLabels", zoom = 14) {
  if (!requireNamespace("maptiles", quietly = TRUE)) {
    warning("Switchbox theme: maptiles not available; no basemap.")
    return(NULL)
  }
  bm <- tryCatch(
    maptiles::get_tiles(bbox_sf, provider = provider, zoom = zoom),
    error = function(e) {
      warning(sprintf("Switchbox theme: basemap tiles unavailable (%s); ",
                      conditionMessage(e)),
              "fall back to white panel.")
      NULL
    }
  )
  if (is.null(bm)) return(NULL)
  bm_df <- as.data.frame(bm, xy = TRUE)
  bm_df$hex <- rgb(bm_df[[3]], bm_df[[4]], bm_df[[5]], maxColorValue = 255)
  bm_df
}

# ---- Logo overlay -----------------------------------------------------------
# Overlays the bundled Switchbox logo in a corner of a finished plot using
# annotation_custom + rasterGrob. Returns the modified ggplot. Position is
# given as normalized-plot coordinates (0-1); defaults to bottom-right.
# Works under both coord_sf (maps) and the default Cartesian coord — it adds
# only an annotation layer and does NOT touch the plot's coordinate system.
add_sb_logo <- function(plot, x = 0.92, y = 0.06, width = 0.08) {
  logo <- png::readPNG(.sb_find_asset("sb_logo.png"))
  grob <- grid::rasterGrob(
    logo,
    x = grid::unit(x, "npc"), y = grid::unit(y, "npc"),
    width = grid::unit(width, "npc"),
    just = c("center", "center"),
    interpolate = TRUE
  )
  plot + annotation_custom(grob)
}

# ---- Export -----------------------------------------------------------------
# PNG @ 300 dpi, white background. Maps: width=10, height=8. Line charts: 7x4.5.
save_sb_fig <- function(plot, path, width = 7, height = 4.5, dpi = 300, ...) {
  ggsave(filename = path, plot = plot,
         width = width, height = height, dpi = dpi, bg = "white", ...)
}
