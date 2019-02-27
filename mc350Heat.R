library(data.table)
library(ggplot2)
library(geoR)
library(MASS)
library(plotrix)
library(maptools)
library(sp)
library(rgeos)
library(rgdal)

data = as.data.table(read.table("~/Downloads/mc350.csv",sep='\t',dec = ',',head = T))

values = data[,lapply(.SD, mean), by = regiao]

border = data.table(x = c(0,0.920,0.920,2.502,2.502,9.562,9.562,0,0), y = c(5.695,5.695,0.988,0.988,0,0,9.739,9.7390,5.695))

regions = list()
regions$Mesas1 = data.table(x = c(5.875,8.69,8.69,5.875,5.875), y = c(0,0,0.658,0.658,0))
regions$Mesas2 = data.table(x = c(5.475,9.562,9.562,5.475,5.475), y = c(1.851,1.851,3.086,3.086,1.851))
regions$Mesas3 = data.table(x = c(5.475,9.562,9.562,5.475,5.475), y = c(4.436,4.436,5.671,5.671,4.436))
regions$CDO = data.table(x = c(5.7,9.562,9.562,5.7,5.7), y = c(6.797,6.797,9.739,9.739,6.797))
regions$Reuniao = data.table(x = c(2.377,5.7,5.7,2.377,2.377), y = c(6.797,6.797,9.739,9.739,6.797))
regions$Mesas4 = data.table(x = c(0,0.768,0.768,0,0), y = c(8.48,8.48,5.695,5.695,8.48))
regions$Copa = data.table(x = c(0.92,3.62,3.62,0.92,0.92), y = c(4.651,4.651,0.988,0.988,4.651))

coords = data.table(regiao = c("mesaTV","mesa_marilia","mesa_kenji","mesa_chris","mesa_gabriel","cdo","reuniao","lounge","copa","cozinha"), x = c(7.6797404, 6.6670895, 8.8113676, 6.4774595, 8.7967807, 7.4256097, 4.0414429, 0.3363637, 2.4368810, 1.6637740), y = c(0.3388008, 2.5414266, 2.5706005, 5.1816602, 5.1378994, 8.6095878, 8.4928924, 7.5739160, 5.5463334, 3.5770984))

regions_df = rbindlist(regions,idcol = T)

plot(border,asp=1,type='l')
lines(regions_df[.id=='Mesas1',.(x,y)])
lines(regions_df[.id=='Mesas2',.(x,y)])
lines(regions_df[.id=='Mesas3',.(x,y)])
lines(regions_df[.id=='CDO',.(x,y)])
lines(regions_df[.id=='Reunião',.(x,y)])
lines(regions_df[.id=='Mesas4',.(x,y)])
lines(regions_df[.id=='Copa',.(x,y)])
points(coords[,.(x,y)], pch=4)

values = merge(values,coords,all.x = T, by="regiao")

#Cria geodata
geodata = list(coords = values[,.(x,y)], data = values$temperature_C, borders = border)
class(geodata) = 'geodata'

points(geodata)
points(geodata, trend="1st")
points(geodata, trend="2nd")

plot(geodata,low=T)
plot(geodata, low=T, trend="1st")
plot(geodata, low=T, trend="2nd")

## Visualizando os dados originais e permutados de posicao para avaliar tendencia
par(mfrow=c(2,2))
points(geodata, pt.div="quint", cex.min=1, cex.max=1)
points(geodata, pt.div="quint", cex.min=1, cex.max=1, permute=T)

points(geodata, trend="1st", pt.div="quint")
points(geodata, trend="1st", pt.div="quint", perm=T)
par(mfrow=c(1,1))

# Verificando necessidade de transformacao e normalidade
boxcox(geodata,trend='1st',lambda=seq(-10,10))
shapiro.test(geodata$data)

qqnorm(geodata$data)
qqline(geodata$data)

grid_length = c(30,30)

# Define grid to interpolate
gr0 = expand.grid(seq(min(geodata$borders[,1]),max(geodata$borders[,1]), len=grid_length[1]), 
                  seq(min(geodata$borders[,2]),max(geodata$borders[,2]), len=grid_length[2]))

# Crop using shape
gr = locations.inside(gr0, geodata$borders)

# Verify spatial dependency via semi-variograms
get_variog = function(geodata){
  variogram = list()
  variogram.env = list()
  
  variogram$v1 = tryCatch(variog(geodata), error = function(cond){return(NULL)})
  variogram.env$v1 = tryCatch(variog.mc.env(geodata, obj=variogram$v1), error = function(cond){return(NULL)})
  
  variogram$v2 = tryCatch(variog(geodata, fix.lam=F), error = function(cond){return(NULL)})
  variogram.env$v2 = tryCatch(variog.mc.env(geodata, obj=variogram$v2), error = function(cond){return(NULL)})
  
  variogram$v3 = tryCatch(variog(geodata, trend="1st"), error = function(cond){return(NULL)})
  variogram.env$v3 = tryCatch(variog.mc.env(geodata, obj=variogram$v3), error = function(cond){return(NULL)})
  
  variogram$v3a = tryCatch(variog(geodata, trend="1st", fix.lam=F), error = function(cond){return(NULL)})
  variogram.env$v3a = tryCatch(variog.mc.env(geodata, obj=variogram$v3a), error = function(cond){return(NULL)})
  
  variogram$v4 = tryCatch(variog(geodata, trend="2nd"), error = function(cond){return(NULL)})
  variogram.env$v4 = tryCatch(variog.mc.env(geodata, obj=variogram$v4), error = function(cond){return(NULL)})
  
  variogram$v4a = tryCatch(variog(geodata, trend="2nd", fix.lam=F), error = function(cond){return(NULL)})
  variogram.env$v4a = tryCatch(variog.mc.env(geodata, obj=variogram$v4a), error = function(cond){return(NULL)})
  
  variogram$v5 = tryCatch(variog(geodata,max.dist=(max(abs(geodata$coords[,1]))-min(abs(geodata$coords[,1])))/2), error = function(cond){return(NULL)})
  variogram.env$v5 = tryCatch(variog.mc.env(geodata, obj=variogram$v5), error = function(cond){return(NULL)})
  
  variogram$v5a = tryCatch(variog(geodata,fix.lam = F,max.dist=(max(abs(geodata$coords[,1]))-min(abs(geodata$coords[,1])))/2), error = function(cond){return(NULL)})
  variogram.env$v5a = tryCatch(variog.mc.env(geodata, obj=variogram$v5a), error = function(cond){return(NULL)})
  
  variogram$v5b = tryCatch(variog(geodata,trend = '1st',max.dist=(max(abs(geodata$coords[,1]))-min(abs(geodata$coords[,1])))/2), error = function(cond){return(NULL)})
  variogram.env$v5b = tryCatch(variog.mc.env(geodata, obj=variogram$v5b), error = function(cond){return(NULL)})
  
  variogram$v5c = tryCatch(variog(geodata,trend = '2nd',max.dist=(max(abs(geodata$coords[,1]))-min(abs(geodata$coords[,1])))/2), error = function(cond){return(NULL)})
  variogram.env$v5c = tryCatch(variog.mc.env(geodata, obj=variogram$v5c), error = function(cond){return(NULL)})
  
  variogram$v5d = tryCatch(variog(geodata,fix.lam=F,trend = '1st',max.dist=(max(abs(geodata$coords[,1]))-min(abs(geodata$coords[,1])))/2), error = function(cond){return(NULL)})
  variogram.env$v5d = tryCatch(variog.mc.env(geodata, obj=variogram$v5d), error = function(cond){return(NULL)})
  
  variogram$v5e = tryCatch(variog(geodata,fix.lam = F,trend = '2nd',max.dist=(max(abs(geodata$coords[,1]))-min(abs(geodata$coords[,1])))/2), error = function(cond){return(NULL)})
  variogram.env$v5e = tryCatch(variog.mc.env(geodata, obj=variogram$v5e), error = function(cond){return(NULL)})
  return(list(variogram = variogram, variogram.env = variogram.env))
}

modeling = function(dependency,variogram,geodata,type){
  # Choosed semi-variogram to define model initial parameters
  if(any(dependency)){
    semivariogram = eval(parse(text = paste0('variogram$',names(dependency[dependency][1]))))
    
    # Define initial parameters
    partial_still = seq(min(semivariogram$v),max(semivariogram$v),length.out = 20)
    range = seq(min(semivariogram$u),max(semivariogram$u),length.out = 20)
    initial = as.data.table(expand.grid(partial_still = partial_still,range = range))
    nugget = semivariogram$v[1]
    
    # List of models to test
    # trends = c("cte","1st",'2nd')
    # distros = c("matern", "exponential", "gaussian", "spherical", "circular", "cubic", "wave", "power", "powered.exponential", "cauchy", "gencauchy", "gneiting", "gneiting.matern", "pure.nugget")
    # 
    # parameters = as.data.table(expand.grid(trends,distros),keep.rownames = T)
    
    mod = list()
    mod$m1 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='matern'), error = function(cond){return(NULL)})
    mod$m2 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='matern'), error = function(cond){return(NULL)})
    mod$m3 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='matern'), error = function(cond){return(NULL)})
    mod$m4 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='exponential'), error = function(cond){return(NULL)})
    mod$m5 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='exponential'), error = function(cond){return(NULL)})
    mod$m6 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='exponential'), error = function(cond){return(NULL)})
    mod$m7 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='gaussian'), error = function(cond){return(NULL)})
    mod$m8 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='gaussian'), error = function(cond){return(NULL)})
    mod$m9 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='gaussian'), error = function(cond){return(NULL)})
    mod$m10 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='spherical'), error = function(cond){return(NULL)})
    mod$m11 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='spherical'), error = function(cond){return(NULL)})
    mod$m12 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='spherical'), error = function(cond){return(NULL)})
    mod$m13 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='circular'), error = function(cond){return(NULL)})
    mod$m14 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='circular'), error = function(cond){return(NULL)})
    mod$m15 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='circular'), error = function(cond){return(NULL)})
    mod$m16 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='cubic'), error = function(cond){return(NULL)})
    mod$m17 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='cubic'), error = function(cond){return(NULL)})
    mod$m18 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='cubic'), error = function(cond){return(NULL)})
    mod$m19 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='wave'), error = function(cond){return(NULL)})
    mod$m20 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='wave'), error = function(cond){return(NULL)})
    mod$m21 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='wave'), error = function(cond){return(NULL)})
    mod$m22 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='power'), error = function(cond){return(NULL)})
    mod$m23 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='power'), error = function(cond){return(NULL)})
    mod$m24 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='power'), error = function(cond){return(NULL)})
    mod$m25 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='powered.exponential'), error = function(cond){return(NULL)})
    mod$m26 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='powered.exponential'), error = function(cond){return(NULL)})
    mod$m27 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='powered.exponential'), error = function(cond){return(NULL)})
    mod$m28 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='cauchy'), error = function(cond){return(NULL)})
    mod$m29 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='cauchy'), error = function(cond){return(NULL)})
    mod$m30 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='cauchy'), error = function(cond){return(NULL)})
    mod$m31 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='gencauchy'), error = function(cond){return(NULL)})
    mod$m32 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='gencauchy'), error = function(cond){return(NULL)})
    mod$m33 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='gencauchy'), error = function(cond){return(NULL)})
    mod$m34 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='gneiting'), error = function(cond){return(NULL)})
    mod$m35 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='gneiting'), error = function(cond){return(NULL)})
    mod$m36 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='gneiting'), error = function(cond){return(NULL)})
    mod$m37 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='gneiting.matern'), error = function(cond){return(NULL)})
    mod$m38 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='gneiting.matern'), error = function(cond){return(NULL)})
    mod$m39 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='gneiting.matern'), error = function(cond){return(NULL)})
    mod$m40 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='cte', cov.model='pure.nugget'), error = function(cond){return(NULL)})
    mod$m41 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='1st', cov.model='pure.nugget'), error = function(cond){return(NULL)})
    mod$m42 = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='2nd', cov.model='pure.nugget'), error = function(cond){return(NULL)})
    
    # models = do.call(function(x,y,z)paste0("mod$m",x ," = tryCatch(likfit(geodata, ini=initial, fix.nug=F, fix.kappa=F, trend='",y,"', cov.model='",z,"'), error = function(cond){return(NULL)})"),list(parameters$rn,parameters$Var1,parameters$Var2))
    # 
    # mod = lapply(as.list(models),function(x) eval(parse(text=x)))
    # print('mod')
    # print(mod)
    print('mod.bic')
    # Choose best model using BIC
    mod.bic = sapply(mod,function(x)x$BIC)
    print(mod.bic)
    print('best_model_names')
    best_model_name = names(sort(abs(mod.bic))[1])
    print(best_model_name)
    print('best_model')
    best_model = eval(parse(text = paste0('mod$',best_model_name)))
    # print('save')
    # save(best_model,file = paste0('../models_',type,'/model_',today(),'.RData'))
    print('kc')
    
    kc = lapply(mod,function(x)tryCatch(krige.conv(geodata, loc=gr0,krige=krige.control(obj.model=x),output=output.control(n.predictive=500)), error = function(cond){return(NULL)}))
    return(list(kc = kc, model = mod))
  }else{
    mod_files = list.files(paste0('../models_',type))
    last_model = sort(mod_files,decreasing = T)[1]
    load(paste0('../models_',type,'/',last_model))
    
    last_model = eval(parse(text = as.character(as.data.table(ll())[data.class == 'likGRF']$member)))
    kc = krige.conv(geodata, loc=gr0,krige=krige.control(obj.model=last_model),output=output.control(n.predictive=500))
    return(list(kc = kc, model = last_model))
  }
}

variogram = get_variog(geodata)$variogram
variogram.env = get_variog(geodata)$variogram.env

## Verificando dependencia espacial por variogramas
plot(variogram$v3,env=variogram.env$v3,xlab="u", ylab=expression(gamma(u)),main='Semi-Variograma')

# Checks for spatial dependency on semi-variogram
dependency = mapply(function(x,y){
  any(x$v < y$v.lower | x$v > y$v.upper)
},variogram,variogram.env)

model = modeling(dependency,variogram,geodata,'density')
kc = model$kc
mod = model$model

my.color2 = colorRampPalette(c('cyan','green', 'yellow','orange',"red"))

# Areas as SpatialPolygons object
# Transform object to SpatialPolygon
pos_x = c(-3,-3)
pos_y = c(10,0)
  
  # Split coords by area
  regions_df2 <- split(regions_df,regions_df$.id)
  # Remove column area
  regions_df2 <- lapply(regions_df2, function(x) { x[,.id:=NULL]})
  # Creates a Polygon object
  regions_poly = sapply(regions_df2, Polygon)
  # add id variable
  regions_poly = lapply(seq_along(regions_poly), function(i) Polygons(list(regions_poly[[i]]), ID = names(regions_poly)[i]))
  # SpatialPolygons object
  regions_poly <- SpatialPolygons(regions_poly)
  # Centroids coords
  labels <- as.data.table(gCentroid(regions_poly,byid=TRUE))
  # Add Labels
  labels[,labels:=row.names(gCentroid(regions_poly,byid=TRUE))]
  # Creates de_para
  labels[,de_para:=1:.N][,`:=`(text=paste(de_para,labels,sep = '. '), leg.x=seq(pos_x[1],pos_x[2],length.out = .N), leg.y=seq(pos_y[1],pos_y[2],length.out = .N))]
  

# par(mar = rep(0, 4))
#10,13,16,34
image(model$kc$m10,col=my.color2(100),xlab = '',ylab='',main='',axes=T,ylim=c(-5,20))
color.legend(xl = 13, xr = 15, yb = 0, yt = 5, gradient =  'y', rect.col = my.color2(100), legend = round(seq(min(kc$m10$predict, na.rm = T),max(kc$m10$predict),length.out = 5),2),cex = 1, srt=0,pos = 2,offset=0)

plot(regions_poly, add = T)
with(labels,text(x, y, de_para, cex = 1,srt=0,family = 'Roboto'))
with(labels,text(leg.x, leg.y, text, cex = 1,srt=0,family = 'Roboto',adj=0))
# dev.off()

load('~/Downloads/mc350.RData')

grid_length = c(100,100)
gr0 = expand.grid(seq(min(geodata$borders[,1]),max(geodata$borders[,1]), len=grid_length[1]), 
                  seq(min(geodata$borders[,2]),max(geodata$borders[,2]), len=grid_length[2]))
gr = locations.inside(gr0, geodata$borders)

kc_density = krige.conv(geodata, loc=gr0,krige=krige.control(obj.model=model$model$m10),output=output.control(n.predictive=500))

save.image("~/Downloads/mc350Completo.RData")

png(paste0('~/Downloads/mc350.png'),units = 'px', width = 1000, height = 1000,res=96)
par(mar = rep(0, 4))
#10,13,16,34
image(kc_density,col=my.color2(100),xlab = '',ylab='',main='',axes=T,ylim=c(-5,15))
color.legend(xl = 12, xr = 14, yb = 0, yt = 5, gradient =  'y', rect.col = my.color2(100), legend = round(seq(min(kc_density$predict, na.rm = T),max(kc_density$predict),length.out = 5),2),cex = 1, srt=0,pos = 2,offset=0.5)

plot(regions_poly, add = T)
with(labels,text(x, y, de_para, cex = 1,srt=0,family = 'Roboto'))
with(labels,text(leg.x, leg.y, text, cex = 1,srt=0,family = 'Roboto',adj=0))
dev.off()
