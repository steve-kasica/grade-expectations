# Course Maps

Shiny app for profiling cross country courses from GPX files.

## Run

Open `course-maps.Rproj` in RStudio, then:

```r
install.packages(c("shiny", "xml2", "ggplot2", "leaflet"))  # first time only
shiny::runApp()
```

Upload a `.gpx` file to see:

- **Elevation profile:** distance (mi) vs. elevation (ft), y-axis fit to the course, with mile markers
- **Elevation by mile:** gain, loss, net and average grade per mile (final partial mile labeled by its distance, e.g. 3.1)
- **Mile markers:** location and elevation of the ½, 1, 2, … marks, plus a course map

- **Flat-equivalent distance:** total and per-mile distance adjusted for grade (see below)

## Flat-equivalent distance

Uses Minetti et al. (2002), the energy cost of running (J/kg/m) at grade *i*, the same
physiology behind Strava's Grade Adjusted Pace:

Cr(i) = 155.4i⁵ − 30.4i⁴ − 43.3i³ + 46.3i² + 19.5i + 3.6

The course is cut into equal segments (slider, default 200 ft). Each segment counts as
its length × Cr(i) / Cr(0), using the smoothed elevation. Grades are clamped to ±45%,
the range the curve was fit on. Shorter segments pick up short, steep hills but also
more GPS noise.

## Smoothing

GPS elevation is noisy, and raw data overstates climbing. The smoothing slider
(default 100 ft) averages elevation over a distance window. Set it to 0 for raw data.

## Course files

Each course has its own folder of GPX tracks (e.g. `all-city/`, `ev-berg/`).
Different watches measure slightly different distances (3.10 to 3.24 mi for the same 5K),
so the last split label follows the track: 3.1, 3.2 and so on.
