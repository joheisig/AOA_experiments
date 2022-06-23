s = stack(sen_ms, prediction, AOA)
d = as.data.frame(s)
d = tidyr::pivot_longer(d, -layer_value)
str(d)

ggplot(d, aes(x=layer_value, y=value, fill=layer_value)) +
  geom_boxplot(outlier.alpha = .5, outlier.color = "grey70") +
  facet_wrap(~name, scales = "free") +
  theme_minimal()

##############
library(stars)
library(sf)
library(raster)
library(mapview)
setwd("/home/jheisig/PhD_home/repos/OpenGeoHub_2020/practice/")
sen_ms <- stack("data/Sen_Muenster.grd") %>% scale()
sen_ms = aggregate(sen_ms, 2)
#sen_ms_df = as.data.frame(sen_ms)
sen_ms_st = st_as_stars(sen_ms)
sen_bb = sen_ms_st %>%  
  st_bbox() %>% 
  st_as_sfc()

reg1 = senST %>% st_make_grid(2000, what = "centers") %>% st_extract(sen_ms_st, .) %>% st_as_sf()
reg2 = senST %>% st_make_grid(1000, what = "centers") %>% st_extract(sen_ms_st, .) %>% st_as_sf()
reg3 = senST %>% st_make_grid(500, what = "centers") %>% st_extract(sen_ms_st, .) %>% st_as_sf()

set.seed(12)
rand1 = sen_bb %>% st_sample(20) %>% st_extract(sen_ms_st, .) %>% st_as_sf()
rand2 = sen_bb %>% st_sample(60) %>% st_extract(sen_ms_st, .) %>% st_as_sf()
rand3 = sen_bb %>% st_sample(120) %>% st_extract(sen_ms_st, .) %>% st_as_sf()


#m = mapview(list(rand1,rand3, strat1, rand2), col.regions = list("lightgreen", "orange", "red", "darkgreen"))
m = mapview(list(rand1,rand2, rand3), col.regions = list("red4", "darkorange","gold2"))
m
###################
sen_df = sen_ms %>% as.data.frame()
sen_reg10 = sampleRegular(sen_ms, ncell(sen_ms)/10)

distfun <- function(x) {
  if (any(is.na(x))) {
    return(NA)
  }
  else {
    tmp <- FNN::knnx.dist(t(matrix(x)), sen_reg10, k = 1)
    return(sd(tmp))
  }
}

library(parallel)
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)
clusterExport(cl=cl, list("sen_reg10"))
system.time(sd_dist <- parApply(cl = cl, X = sen_df[], MARGIN = 1, FUN=distfun))
beepr::beep(3)

sd_ras = sen_ms[[1]] %>% setValues(sd_dist) 



plot(log(sd_ras))
hist(log(sd_ras))



n = 6  # number of equally sized distance classes
q = quantile(log(sd_ras), probs = seq(0,1, length.out = n+1))
rcl_mat = matrix(c(q[1:n], q[2:(n+1)], 1:n), nrow=n)

sd_class = reclassify(log(sd_ras), rcl_mat, include.lowest=T)
plot(sd_class, col=terrain.colors(n))
hist(sd_class)

spplot(log(sd_ras), col.regions = rev(grey.colors(100)), main="sd() of within-image pixel distances")
spplot(sd_class, col.regions=rainbow(n+1), main = "reclassified distance image",
       colorkey = list(space = "right", at = c(0:n)+.5))

#----

set.seed(11)
set.seed(123)
strat1 = sampleStratified(sd_class, 3, xy=T, sp=T)[,2:3] %>% st_as_sf() %>% st_extract(sen_ms_st, .) %>% st_as_sf()
strat2 = sampleStratified(sd_class, 10, xy=T, sp=T)[,2:3] %>% st_as_sf() %>% st_extract(sen_ms_st, .) %>% st_as_sf()
strat3 = sampleStratified(sd_class, 20, xy=T, sp=T)[,2:3] %>% st_as_sf() %>% st_extract(sen_ms_st, .) %>% st_as_sf()

mapview(sd_class, col.regions=rev(terrain.colors(n-1))) +
mapview(list(strat1, strat2, strat3), col.regions=list("darkgreen", "olivedrab", "lightgreen")) 

m + mapview(list(strat1, strat2, strat3), col.regions=list("darkgreen", "olivedrab", "lightgreen")) 


mapview(log(sd_ras), col.region = rev(grey.colors(100)), alpha.regions=1, legend=F) +
mapview(list(strat3, rand3), color=NA, alpha=0, col.regions = list("darkgreen", "red4"),
        layer.name = list("stratified-random", "random"))

#===============
train1=st_drop_geometry(rand1)
train2=st_drop_geometry(rand2)
train3=st_drop_geometry(rand3)
train4=st_drop_geometry(strat1)
train5=st_drop_geometry(strat2)
train6=st_drop_geometry(strat3)

library(CAST)
aoa1 = aoa(sen_ms, train = train1, cl = cl) %>% setNames(paste0(names(.), "_n_",nrow(rand1),"_rand."))
aoa2 = aoa(sen_ms, train = train2, cl = cl) %>% setNames(paste0(names(.), "_n_",nrow(rand2),"_rand."))
aoa3 = aoa(sen_ms, train = train3, cl = cl) %>% setNames(paste0(names(.), "_n_",nrow(rand3),"_rand."))
aoa4 = aoa(sen_ms, train = train4, cl = cl) %>% setNames(paste0(names(.), "_n_",nrow(strat1),"_strat."))
aoa5 = aoa(sen_ms, train = train5, cl = cl) %>% setNames(paste0(names(.), "_n_",nrow(strat2),"_strat."))
aoa6 = aoa(sen_ms, train = train6, cl = cl) %>% setNames(paste0(names(.), "_n_",nrow(strat3),"_strat."))

di = stack(aoa1[[1]], aoa4[[1]],  aoa2[[1]], aoa5[[1]], aoa3[[1]],  aoa6[[1]]) 
dilog = di %>% log() %>% 
  spplot(layout = c(2,3), colorkey = list(space = "bottom"), main="log( Dissimilarity Index )")
di = di %>% spplot(layout = c(2,3))

ao = stack(aoa1[[2]], aoa4[[2]],  aoa2[[2]], aoa5[[2]], aoa3[[2]],  aoa6[[2]]) %>% 
  spplot(col.regions = c("tomato1","lightseagreen"), layout = c(2,3), main = "AOA",
       colorkey = list(space = "bottom", at = c(0,.5,1),
                       labels = list(labels =c("outside", "inside"),
                                     at = c(0.25, 0.75))))

gridExtra::grid.arrange(dilog, ao, ncol=2)

di
dilog
ao
di+ao


stack(aoa1[[1]], aoa2[[1]],  aoa3[[1]]) %>% plotRGB(stretch='lin')
stack(aoa1[[2]], aoa2[[2]],  aoa3[[2]]) %>% plotRGB(stretch='lin')

sum(values(aoa1$AOA_rand1) == 0) / ncell(aoa1)
sum(values(aoa2$AOA_rand2) == 0) / ncell(aoa1)
sum(values(aoa3$AOA_rand3) == 0) / ncell(aoa1)
sum(values(aoa4$AOA_strat1) == 0) / ncell(aoa1)



#------------------
senS = scale(sen) %>% as.data.frame()
reg = sampleRegular(scale(sen), ncell(sen)/10)

distfun <- function(x) {
  if (any(is.na(x))) {
    return(NA)
  }
  else {
    tmp <- FNN::knnx.dist(t(matrix(x)), reg, k = 1)
    return(mean(tmp))
  }
}


cl <- makeCluster(6)
registerDoParallel(cl)
clusterExport(cl=cl, list("reg"))
meandist <- parApply(cl = cl, X = senS[], MARGIN = 1, FUN=distfun)

mdras = sen[[1]]
values(mdras) = meandist
plot(log(mdras))
hist(log(mdras))

q=quantile(log(mdras))
rclmat = matrix(c(q[1:length(q)-1], q[2:length(q)], 2:length(q)-1), nrow=length(q)-1)

mdclass = reclassify(log(mdras), rclmat)
plot(mdclass)

strat2 = sampleStratified(mdclass, 10, xy=T, sp=T)[,2:3] %>% st_as_sf() %>% st_extract(senST, .) %>% st_as_sf()

mapview(strat2, col.regions="black") + m

train5=st_drop_geometry(strat2)
aoa5 = aoa(sen, train = train5, cl = cl) %>% setNames(paste0(names(.), "_strat2"))
spplot(aoa5)
#------------------------





distfun <- function(x) {
  if (any(is.na(x))) {
    return(NA)
  }
  else {
    tmp <- FNN::knnx.dist(t(matrix(x)), train, k = 1)
    return(min(tmp))
  }
}

rowdistfun <- function(x){
  tmp = FNN::knnx.dist(train, t(matrix(x)), k = nrow(train))
}

train=st_drop_geometry(rand1)
mindist1 <- apply(senD, 1, FUN = distfun)
tDist_avrg1 <- apply(senD, 1, FUN = function(x) mean(rowdistfun(x)[-1], na.rm = TRUE))
trainDist_avrgmean1 <- mean(trainDist_avrg, na.rm = TRUE)
trainDist_min1 <- c()
for (i in 1:nrow(train)) {
  trainDist <- FNN::knnx.dist(t(matrix(train[i, ])), train, k = 1)
  trainDist[i] <- NA
  trainDist_min1 <- c(trainDist_min, min(trainDist, na.rm = T))
}
DI_out <- mindist1/trainDist_avrgmean1

train=st_drop_geometry(rand2)
mindist2 <- apply(senD, 1, FUN = distfun)
tDist_avrg2 <- apply(sen, 1, FUN = function(x) mean(rowdistfun(x)[-1], na.rm = TRUE))

train=st_drop_geometry(rand3)
mindist3 <- apply(senD, 1, FUN = distfun)
tDist_avrg3 <- apply(sen, 1, FUN = function(x) mean(rowdistfun(x)[-1], na.rm = TRUE))

#####

k=kmeans(as.matrix(sen_ms), 6)

sen_ms$cluster = as.factor(k$cluster)
plot(sen_ms$center, col=sf.colors(6, categorical=T))
