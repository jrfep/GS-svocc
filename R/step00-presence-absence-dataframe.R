#! R --vanilla
args <- commandArgs(TRUE)
sp <- args[1]

out.file <- sprintf("Rdata/padata/%s.rda",sp)
if (!file.exists(out.file)) {
  require(dplyr)
  require(tidyr)
  require(raster)
  require(chron)
  require(detect)

  GIS.data <- "Rdata/GIS.rda"
  options(timeout = max(300, getOption("timeout")))
  
  if (!file.exists(GIS.data))
    download.file(url="https://figshare.com/ndownloader/files/37547995",destfile=GIS.data, extra="--continue")
  event.data <- "Rdata/all-events.csv"
  if (!file.exists(event.data))
    download.file(url="https://figshare.com/ndownloader/files/42055824",destfile=event.data)
  
  load(GIS.data)
  eventos_actualizados <- read.csv2(event.data)
  coordinates(eventos_actualizados) <- c("long","lat")
  crs(eventos_actualizados) <- grd@proj4string
  qry <- over(eventos_actualizados,grd)
  eventos_actualizados$grid <- qry$OID_
  
  ## differences between both event data frames:
  eventos_adicionales <- eventos %>% 
    filter(!species %in% eventos_actualizados$species) %>%
    mutate(fotos=as.character(fotos))
  
  eventos <- eventos_actualizados@data %>% 
    bind_rows(eventos_adicionales)
  
  
  #tps@data %>% group_by(grid) %>% tally() %>% transmute(grid,walk=(n-mean(n))/sd(n))-> walk

  tps_xy <- spTransform(tps, CRS("+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +type=crs"))
  xys <- coordinates(tps_xy)
  
  
  walk <- tps@data %>% 
    group_by(grid) %>% 
    summarise(gps_log = n(),
              dates = n_distinct(time))
  
  walk_ms <- data.frame()
  for (k in unique(tps_xy@data$grid)) {
    ss <- tps_xy@data$grid == k
    mtz <- as.matrix(dist(xys[ss,]))
    mtz[upper.tri(mtz,diag = TRUE)] <- NA
    walk_ms <- rbind(walk_ms, data.frame(grid=k, distance=sum(apply(mtz[-1,],1, min, na.rm=TRUE))))
  }
  
  walk <- walk %>% left_join(walk_ms, by = "grid") %>% 
    transmute(grid, walk=(distance-mean(distance))/sd(distance))
  
  
  camaras %>% group_by(grid) %>% summarise(cam=sum(dias.de.trabajo),caz=max(caza.celda)) %>% transmute(grid,caz,cam=(cam-mean(cam))/sd(cam)) -> cams

  dts <- data.frame(
      bsq=values(vbsq),
      dcon=values(dist.conucos),
      dcom=values(dist.comunidades),
      frs=values(dist.frs),
      dbsq=values(dist.dbsq),
      dcaz1=values(dist.caza1),
      dcaz2=values(dist.caza2),
      grid=values(rgrd)) %>%
    group_by(grid) %>%
    summarise(
      bsq=mean(bsq),
      dbsq=mean(dbsq),
      dcon=mean(dcon),
      dcom=mean(dcom),
      dpob=mean(min(dcon,dcom)),
      dhum=mean(min(dcon,dcom,dcaz1,dcaz2)),
      dcaz1=mean(dcaz1),
      dcaz2=mean(dcaz2),
      dcaz=mean(min(dcaz1,dcaz2)),
      frs=mean(frs))

  ##for (sp in levels(eventos$species)) {
  eventos %>%
    mutate(
      target=species %in% sp,
      in_camera=if_else(camara %in% "RAS",FALSE,species %in% sp),
      in_walk=if_else(camara %in% "RAS",species %in% sp,FALSE)) %>%
    group_by(grid) %>%
    summarise(
      pa=max(target),
      pa_cam=max(in_camera),
      pa_ras=max(in_walk)) -> species


  dts  %>%
    left_join(walk) %>% full_join(cams) %>% full_join(species) %>%
    filter(!is.na(grid) & (!is.na(walk) | !is.na(cam))) %>%
    transmute(grid, walk, cam,
      caz=as.factor(!caz %in% 0),
      pa=coalesce(pa,0L),
      pa_cam=coalesce(pa_cam,0L),
      pa_ras=coalesce(pa_ras,0L),
      bsq=(bsq-mean(bsq))/sd(bsq),
      dcon=(dcon-mean(dcon,na.rm=T))/sd(dcon,na.rm=T),
      dcom=(dcom-mean(dcom,na.rm=T))/sd(dcom,na.rm=T),
      dhum=(dhum-mean(dhum,na.rm=T))/sd(dhum,na.rm=T), ## all human activities together (com, con, caz)
      dpob=(dpob-mean(dpob,na.rm=T))/sd(dpob,na.rm=T),
      dcaz=(dcaz-mean(dcaz,na.rm=T))/sd(dcaz,na.rm=T),
      frs=(frs-mean(frs,na.rm=T))/sd(frs,na.rm=T),
      dbsq=(dbsq-mean(dbsq,na.rm=T))/sd(dbsq,na.rm=T)) -> pa.data

  ##colSums(is.na(pa.data))
  pa.data$walk <- coalesce(pa.data$walk,min(walk$walk,na.rm=T))
  pa.data$cam <- coalesce(pa.data$cam,min(cams$cam,na.rm=T))
  pa.data$bloque <- factor(grd@data$cuadrado[match(pa.data$grid,grd@data$OID_)])

  ## sampling area (Warapata vs. Kavanayen) and method (Cam+walk vs Only walk)
  pa.data$muestreo <- pa.data$bloque %in% 1:6
  pa.data$metodo <- pa.data$grid %in% cams$grid


  save(file=out.file,pa.data)
}
