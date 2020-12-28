library(ecopower)
library(ecoCopula)
library(mvabund)

data(spider)
spiddat = mvabund(spider$abund)

X = data.frame(spider$x)
X$Treatment2 = rep(c("A","B"),each=14)
X$Treatment4 = rep(c("A","B","C","D"),each=7)
X$month = rep(c("Jan","Feb","Mar","Apr"),each=7)

Xnew = X
Xnew$Treatment2[6:7] = c("B","B")

Xvec = Xvec2 = data.frame(Treatment2 = X$Treatment2)
Xvec2$Treatment2[6:7] = c("B","B")

term_cont = "bare.sand"

increasers = c("Alopacce", "Arctlute", "Arctperi", "Pardnigr", "Pardpull")
decreasers = c("Alopcune", "Alopfabr", "Zoraspin")

fit0.glm = manyglm(spiddat~1, family="negative.binomial")
fit0.cord = cord(fit0.glm)

fit0_X.glm = manyglm(spiddat~1, family="negative.binomial", data=X)
fit0_X.cord = cord(fit0_X.glm)

fit1.glm = manyglm(spiddat~bare.sand, family="negative.binomial", data=X)
fit1.cord = cord(fit1.glm)

fit2.glm = manyglm(spiddat~soil.dry+bare.sand, family="negative.binomial", data=X)
fit2.cord = cord(fit2.glm)

fit_mix.glm = manyglm(spiddat~month+bare.sand, family="negative.binomial", data=X)
fit_mix.cord = cord(fit_mix.glm)

fit_fac2.glm = manyglm(spiddat~Treatment2, family="negative.binomial", data=X)
fit_fac2.cord = cord(fit_fac2.glm)

fit_vec.glm = manyglm(spiddat~Treatment2, family="negative.binomial", data=Xvec)
fit_vec.cord = cord(fit_vec.glm)

fit_fac4.glm = manyglm(spiddat~Treatment4, family="negative.binomial", data=X)
fit_fac4.cord = cord(fit_fac4.glm)

fit_mth.glm = manyglm(spiddat~month+Treatment2, family="negative.binomial", data=X)
fit_mth.cord = cord(fit_mth.glm)

effect_size = 1.5
effect_mat = effect_alt(fit1.glm, effect_size, increasers, decreasers, term_cont)
