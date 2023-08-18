require(dplyr)
require(tidyr)
require(raster)
require(chron)
require(detect)

GIS.data <- "Rdata/GIS.rda"
if (!file.exists(GIS.data))
  download.file(url="https://figshare.com/ndownloader/files/37547995",destfile=GIS.data)
load(GIS.data)

event.data <- "Rdata/all-events.csv"
if (!file.exists(event.data))
  download.file(url="https://figshare.com/ndownloader/files/42055824",destfile=event.data)
eventos_completos <- read.csv2(event.data)
coordinates(eventos_completos) <- c("long","lat")
crs(eventos_completos) <- grd@proj4string
qry <- over(eventos_completos,grd)
eventos_completos$grid <- qry$OID_

## differences between both event data frames:
table(droplevels(subset(eventos,!species %in% eventos_completos$species)$species))

# Walks per grid

tps_xy <- spTransform(tps, CRS("+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +type=crs"))
xys <- coordinates(tps_xy)


walk <- tps@data %>% 
  group_by(grid) %>% 
  summarise(walk = n(),
            dates = n_distinct(time))

walk_ms <- data.frame()
for (k in unique(tps_xy@data$grid)) {
  ss <- tps_xy@data$grid == k
  mtz <- as.matrix(dist(xys[ss,]))
  mtz[upper.tri(mtz,diag = TRUE)] <- NA
  walk_ms <- rbind(walk_ms, data.frame(grid=k, distance=sum(apply(mtz[-1,],1, min, na.rm=TRUE))))
}

walk <- walk %>% left_join(walk_ms, by = "grid")

plot(distance~walk, walk)

cam <- camaras %>% 
  group_by(grid) %>% 
  summarise(cam = sum(dias.de.trabajo),
            caz = max(caza.celda)) %>% 
  transmute(grid, caz, cam)


event_summary <- eventos_completos@data %>%
  group_by(grid) %>%
  summarise(
    total_spp=n_distinct(species))

on_camera_event_summary <- eventos_completos@data %>%
  filter(!camara %in% "RAS") %>% 
  group_by(grid) %>%
  summarise(
    on_camera_spp=n_distinct(species),
    on_camera_events=n())

off_camera_event_summary <- eventos_completos@data %>%
  filter(camara %in% "RAS") %>% 
  group_by(grid) %>%
  summarise(
    off_camera_spp=n_distinct(species),
    off_camera_events=n())

table(eventos_completos@data$bloque)

event_summary <- event_summary %>% 
  left_join(on_camera_event_summary, by = "grid") %>% 
  left_join(off_camera_event_summary, by = "grid")

effort_data <- walk %>%
  left_join(cam, by = "grid") %>%
  left_join(event_summary, by = "grid") %>%
  filter(!is.na(grid) & (!is.na(walk) | !is.na(cam))) %>%
  transmute(grid, 
            walk=coalesce(walk,0L),
            distance=coalesce(distance,0L),
            cam=coalesce(cam,0L),
            total_spp=coalesce(total_spp,0L),
            on_camera_spp=coalesce(on_camera_spp,0L),
            off_camera_spp=coalesce(off_camera_spp,0L),
            on_camera_events=coalesce(on_camera_events,0L),
            off_camera_events=coalesce(off_camera_events,0L)) 


effort_data$bloque <- factor(grd@data$cuadrado[match(effort_data$grid,grd@data$OID_)])

effort_data %>% 
  group_by(bloque) %>% 
  summarise(`nr. of s. u.`=n_distinct(grid),
            `gps log points` = sum(walk),
            `distance in m` = sum(distance),
            `camera days` = sum(cam),
            `events on camera` = sum(on_camera_events),
            `events off camera` = sum(off_camera_events))
