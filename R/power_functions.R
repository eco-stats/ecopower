utils::globalVariables(c(
  "object",
  "N",
  "coeffs",
  "newdata",
  "n_replicate",
  "do.fit",
  "test",
  "term",
  "object0",
  "n.samp",
  "nlv",
  "nsim"
))

MVApowerstat = function(stats, coeffs) {
  anova(
    extend(
      object=object,
      N=N,
      coeffs=coeffs,
      newdata=newdata,
      n_replicate=n_replicate,
      do.fit=do.fit
    ),
    nBoot=1,
    test=test,
    show.time = "none"
  )$table[term,3]
}



MVApowerstat_long_alt = function(stats, coeffs) {
  # alt_mod = tryCatch(cord(extend(
  #   object=object,
  #   N=N,
  #   coeffs=coeffs,
  #   newdata=newdata,
  #   n_replicate=n_replicate,
  #   do.fit=TRUE
  # ),n.samp=n.samp,nlv = nlv),
  # error = message("Error in factanal, try increasing n.samp and decreasing nlv"),
  # finally=cord(extend(
  #   object=object,
  #   N=N,
  #   coeffs=coeffs,
  #   newdata=newdata,
  #   n_replicate=n_replicate,
  #   do.fit=TRUE
  # ),n.samp=n.samp,nlv = nlv))

  boolFalse<-F
  while(boolFalse==F)
  {
    tryCatch({
      alt_mod = cord(extend(
        object=object,
        N=N,
        coeffs=coeffs,
        newdata=newdata,
        n_replicate=n_replicate,
        do.fit=TRUE
      ),n.samp=n.samp,nlv = nlv)
      boolFalse<-T
    },error=function(e){warning("Error in factanal, try increasing n.samp and decreasing nlv if this code takes too long to run.")
    },finally={})
  }

  # extended_data <<- data.frame(alt_mod$obj$x)
  assign("extended_data", data.frame(alt_mod$obj$x), inherits=TRUE)
  coeffs0_l = effect_null(alt_mod$obj, term=term)

  stat = anova(alt_mod$obj,
               nBoot=1,
               test=test,
               show.time = "none"
  )$table[term,3]
  results = list(stat,alt_mod,coeffs0_l)
}



MVApowerstat_long_null = function(stats,alt_mod,coeffs) {
  anova(
    extend(
      object=alt_mod,
      N=N,
      coeffs=coeffs,
      newdata=newdata,
      n_replicate=n_replicate,
      do.fit=do.fit
    ),
    nBoot=1,
    test=test,
    show.time = "none"
  )$table[term,3]
}

MVApowerstat_long_null_2 = function(stats,alt_mods,coeffs) {
  mod = sample(1:nsim,1)
  object = alt_mods[[mod]]
  coeffs = coeffs[[mod]]
  res = anova(
    extend(
      object=object,
      N=N,
      coeffs=coeffs,
      newdata=newdata,
      n_replicate=n_replicate,
      do.fit=do.fit
    ),
    nBoot=1,
    test=test,
    show.time = "none"
  )$table[term,3]
}



MVApowerstatObj = function(stats) {
  object_sim = extend(
    object=object,
    N=N,
    coeffs=coeffs,
    newdata=newdata,
    n_replicate=n_replicate,
    do.fit=do.fit
  )

  object0_sim = update(
    object_sim,
    formula=object0$formula[-2],
    data=object_sim$data
  )

  anova(
    object0_sim, object_sim, nBoot=1, test=test, show.time="none"
  )$table[2,3]
}
