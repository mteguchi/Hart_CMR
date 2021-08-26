
capture.recapture.stats <- function(loc, dat.1, Ei.inputs, save.fig, fig.height, fig.width){
  
  CJS.data <- dat2CJS(dat.1, save.file = FALSE)
  n.occ <- ncol(CJS.data$data)
  n.caught <- colSums(CJS.data$data)
  
  p.captures <- ggplot(dat.1) + 
    geom_point(aes(x = DATE,  y = ID)) +
    geom_path(aes(x = DATE, y = ID)) +
    theme(axis.text.y = element_blank()) +
    labs(title = loc)
  
  if (save.fig)
    ggsave(filename = paste0("figures/Ei_capture_history_", loc, ".png"),
           plot = p.captures, device = "png", dpi = 600,
           height = fig.height, width = fig.width)
  
  n.caps <- rowSums(CJS.data$data) %>% 
    data.frame() %>% 
    rownames_to_column(var = "ID") 
  
  colnames(n.caps) <- c("ID", "n")   # there has to be a better way but... 
  
  p.n.captures <- ggplot(data = n.caps) + 
    geom_histogram(aes(x = n), binwidth = 1) +
    labs(title = loc)
  
  if (save.fig)
    ggsave(filename = paste0("figures/Ei_capture_histogram_", loc, ".png"),
           plot = p.n.captures, device = "png", dpi = 600,
           height = fig.height, width = fig.width)
  
  tmp.counts <- t(apply(CJS.data$data, MARGIN = 1, FUN = n.captures))
  tmp.counts.df <- data.frame(count.captures(tmp.counts)) 
  rownames(tmp.counts.df) <- 1:max(tmp.counts)
  colnames(tmp.counts.df) <- colnames(CJS.data$data)
  tmp.counts.df <- rownames_to_column(tmp.counts.df, var = "n.recaptures")
  tmp.counts.df %>% pivot_longer(!n.recaptures, names_to = "Date", values_to = "Freq") -> counts.freq
  
  p.n.recap <- ggplot(counts.freq) + 
    geom_col(aes(x = Date, y = Freq, fill = n.recaptures)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "top") +
    ylab("# individuals") + labs(fill = "# recaptures", title = loc)
  
  
  if (save.fig)
    ggsave(filename = paste0("figures/Ei_recaptures_", loc, ".png"),
           plot = p.n.recap, device = "png", dpi = 600,
           height = fig.height, width = fig.width)
  
  loc.stats <- list(CJS.data = CJS.data,
                    n.occ = n.occ,
                    n.caught = n.caught,
                    n.caps = n.caps)
  
  return(loc.stats)
}

# Counts the number of individuals that are caught x times at each capture occasion. Use this with output from n.captures.
count.captures <- function(x){
  max.n <- max(x)
  out.table <- matrix(nrow = max.n, ncol = ncol(x))
  for (k in 1:ncol(x)){
    for (k1 in 1:max.n){
      out.table[k1, k] <- sum(x[,k] == k1)
    }
  }
  return(out.table)
}


# this function enumerate the capture sequence for each individual. Rather than having 0s and 1s, they are turned into 0s and 1s, 2s, 3s, etc. Cumulative sum but keeping zeros as zeros. Use it with apply to a capture history dataframe. When used with apply, the output is individuals in columns and capture occasions in rows. So, transform the output to make it comparable to the capture history matrix.
#t(apply(X, MARGIN = 1, FUN = n.captures)), where X is the capture history matrix/dataframe
n.captures <- function(x){
  #x.1 <- vector(mode = "numeric", length = ncol(x))
  x.1 <- x[1]
  for (k in 2:length(x)){
    if (x[k] > 0) { 
      x[k] <- x.1 + x[k]
      x.1 <- x[k]
    }
  }
  
  return(x)
}



extract.Nhats <- function(Cm.inputs, real.estimates){
  data.0 <- Cm.inputs$CJS.data$data
  
  phats <- real.estimates[grep("p", real.estimates$parameter),]
  
  phats[which(phats[,"estimate"] < 0.001), c("estimate", "se", "lcl", "ucl")] <- NA
  phats$season <- colnames(data.0)[1:(ncol(data.0)-1)]
  
  n.caught <- colSums(data.0)
  
  model.averaged.Phi <- model.average(Cm.results, parameter = "Phi")
  
  Nhats.df <- data.frame(season = colnames(data.0)[1:(ncol(data.0)-1)],
                         Nhat = (n.caught[1:(length(n.caught) - 1)]/phats$estimate) ) %>%
    mutate(SE_Nhat = (n.caught[1:(length(n.caught) - 1)]/phats$estimate) * phats$se/phats$estimate,
           #lcl  = (n.caught[2:length(n.caught)]/phats$lcl) * p.residents,
           #ucl = (n.caught[2:length(n.caught)]/phats$ucl) * p.residents,
           lcl = (n.caught[1:(length(n.caught)-1)]/phats$estimate)  - 1.96 * SE_Nhat,
           ucl = (n.caught[1:(length(n.caught)-1)]/phats$estimate)  + 1.96 * SE_Nhat,
           lcl2 = ifelse(lcl < 0, 0, lcl))
  
  return(Nhats.df)
}



do_analysis <- function(dp, ddl)
{
  # create formulas for Phi
  # tsm is time-since-marking; check for transient effects
  Phi.dot <-  list(formula = ~ 1)  
  #Phi.weight <- list(formula= ~ min_weight)   # many missing data 
  #Phi.t <- list(formula = ~ time)             # we never have this model worked for turtles... 
  #Phi.season <- list(formula = ~ sum_win)      # this also is unlikely... 
  #Phi.transience <- list(formula = ~ Transient)
  Phi.tsm <- list(formula = ~ tsm)
  
  #create formulas for p
  p.dot <- list(formula = ~ 1)
  p.t <- list(formula = ~ time)
  #p.tsm <- list(formula = ~ tsm)
  #p.transience <- list(formula = ~ Transient)
  #p.tsm.transience <- list(formula = ~ tsm + Transient)
  #p.t.transience <- list(formula = ~ time + Transient)
  p.effort <- list(formula = ~ effort)
  p.season <- list(formula = ~ sum_win)
  p.tsm.season <- list(formula = ~ tsm + sum_win)
  p.tsm.effort <- list(formula = ~ tsm + effort)
  
  # create all combinations 
  cml <- create.model.list("CJS")
  
  # run all all models and return as a list with class marklist
  results <- mark.wrapper(cml,
                          data=dp,
                          ddl=ddl,
                          output=FALSE,
                          silent=TRUE)
  return(results)
}


Ei.vonBert.jags.data <- function(dat.1){
  # remove rows with is.na(CCL) is true:
  dat.1 %>% filter(!is.na(CCL)) -> dat.2
  
  n.cap.ID <- table(dat.2$ID)
  recap.ID <- data.frame(n.cap.ID[n.cap.ID > 2])
  colnames(recap.ID) <- c("ID", "Freq")
  
  recap.ID %>% left_join(dat.2, by = "ID") -> recap.data
  
  # Make length and capture date matrices
  unique.ID <- recap.ID$ID
  size.mat <- date.mat <- matrix(nrow = length(unique.ID),
                                 ncol = max(recap.data$Freq))
  
  date.1 <- structure(numeric(length(unique.ID)), class = "Date")
  n.vec <- vector(mode = "numeric", length = length(unique.ID))
  
  k <- 1
  for (k in 1:length(unique.ID)){
    tmp.ID <- filter(recap.data, ID == as.character(unique.ID[k]))
    size.mat[k, 1:nrow(tmp.ID)] <- tmp.ID$CCL
    date.mat[k, 1:nrow(tmp.ID)] <- tmp.ID$DATE - min(tmp.ID$DATE)
    date.1[k] <- min(tmp.ID$DATE)
    n.vec[k] <- nrow(tmp.ID)
  }
  
  date.mat <- date.mat[, 2:ncol(date.mat)]/365
  
  jags.data <- list(nIndiv = length(unique.ID),
                    n = n.vec,
                    L = size.mat,
                    t = date.mat)
  
  return(list(jags.data = jags.data,
              ID = unique.ID))
}



vonBert.jags.data <- function(dat.1){
  # remove rows with is.na(CCL) is true:
  dat.1 %>% 
    filter(!is.na(CCL)) -> dat.2
  
  n.cap.ID <- table(dat.2$ID)
  recap.ID <- data.frame(n.cap.ID[n.cap.ID > 2])
  colnames(recap.ID) <- c("ID", "Freq")
  
  recap.ID %>% left_join(dat.2, by = "ID") -> recap.data
  
  # Make length and capture date matrices
  unique.ID <- recap.ID$ID
  size.mat <- date.mat <- matrix(nrow = length(unique.ID),
                                 ncol = max(recap.data$Freq))
  
  date.1 <- structure(numeric(length(unique.ID)), class = "Date")
  n.vec <- vector(mode = "numeric", length = length(unique.ID))
  
  k <- 1
  for (k in 1:length(unique.ID)){
    tmp.ID <- filter(recap.data, ID == as.character(unique.ID[k]))
    size.mat[k, 1:nrow(tmp.ID)] <- tmp.ID$CCL
    date.mat[k, 1:nrow(tmp.ID)] <- tmp.ID$DATE - min(tmp.ID$DATE)
    date.1[k] <- min(tmp.ID$DATE)
    n.vec[k] <- nrow(tmp.ID)
  }
  
  date.mat <- date.mat[, 2:ncol(date.mat)]/365
  
  jags.data <- list(nIndiv = length(unique.ID),
                    n = n.vec,
                    L = size.mat,
                    t = date.mat)
  
  return(list(jags.data = jags.data,
              ID = unique.ID))
}

# type_site_spec, site_name were deleted because there are commas in the field. 2021-08-25
# Also, Data_owner, Contact and Use were deleted 2021-08-25
get.data.Ei <- function(filename){
  col.def <- cols(monitoring_event = col_character(),
                  value = col_integer(),
                  season = col_character(),
                  year = col_integer(),
                  month = col_integer(),
                  day = col_integer(),
                  species = col_factor(),
                  turtle_code = col_factor(),
                  recapture = col_character(),
                  community = col_factor(),
                  start_date = col_date(format = "%m/%d/%Y"),
                  tot_hours = col_double(),
                  type_monitoring = col_factor(),
                  methodology = col_character(),
                  longitude_net = col_double(),
                  type_site_gen = col_factor(),
                  latitude = col_double(),
                  longitude = col_double(),
                  SCL_CCL = col_double(),
                  capture_date = col_date(format = "%m/%d/%Y"),
                  name = col_character(),
                  SCL_min = col_double(),
                  SCL_max = col_double(),
                  SCW = col_double(),
                  CCL_min = col_double(),
                  CCL_max = col_double(),
                  CCW = col_double(),
                  BD = col_double(),
                  PL = col_double(),
                  TTL = col_double(),
                  Weight_kg = col_double(),
                  sex = col_factor(),
                  Marca_n_dx = col_character())
  
  dat.1 <- read_csv(file = filename, col_types = col.def)
  
  dat.1 %>% mutate(ID = turtle_code,
                   CDATE = capture_date) %>%
    transmute(ID = ID,
              detect = 1,
              DATE = CDATE,
              season = as.factor(season),
              SCL = SCL_max,
              CCL = CCL_max,
              weight_kg = Weight_kg,
              sex = sex,
              community = community)  -> dat.2
  
  return(dat.2)
}



dat2CJS <- function(dat.1, save.file = FALSE, filename = "not.saved"){
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.1, 
             id.vars = c("ID", "season"), 
             measure.vars = "detect")
  
  # make a table with ID by season
  dat.01 <- cast(tmp, 
                 formula = ID ~ season,
                 fun.aggregate = length)
  
  # replace > 1 with ones
  dat.01 <- as.data.frame(dat.01) %>%
    remove_rownames() %>%
    column_to_rownames(var = "ID")
  
  dat.01[(dat.01 > 1)] <- 1
  
  # save file for later
  if (save.file){
    out.name <- filename
    write.csv(dat.01, 
              file = out.name, 
              row.names = T,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.01)
  
  return(out)
  
}


known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

dat2dat01.year <- function(dat.1, year, save.file = FALSE){
  dat.year <- filter(dat.1, YEAR == year)
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.year, 
             id.var = c("ID", "DATE"), 
             measure.var = "detect")
  
  # make a table with ID by Date
  dat.01.year <- cast(tmp, ID ~ DATE)
  
  # replace NAs with zeros
  dat.01.year[is.na(dat.01.year)] <- 0
  dat.01.year <- as.data.frame(dat.01.year)
  
  # save file for later
  if (save.file){
    out.name = paste0("data/Cm_01_", year, ".csv")
    write.csv(dat.01.year, 
              file = out.name, 
              row.names = F,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.01.year)
  
  return(out)
  
}


dat2dat01 <- function(dat.1, save.file = FALSE){
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.1, 
             id.var = c("ID", "season"), 
             measure.var = "detect")
  
  # make a table with ID by year
  dat.01 <- cast(tmp, ID ~ season)
  
  # replace > 1 with ones
  dat.01 <- as.data.frame(dat.01) %>%
    remove_rownames() %>%
    column_to_rownames(var = "ID")
  
  dat.01[(dat.01 > 1)] <- 1
  
  # save file for later
  if (save.file){
    out.name <- "data/Ei_01_all.csv"
    write.csv(dat.01, 
              file = out.name, 
              row.names = T,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.01)
  
  return(out)
  
}


dat2CJS_covCCL <- function(dat.1){
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.1, 
             id.vars = c("ID", "season"), 
             measure.vars = "CCL")
  
  # make a table with ID by season
  dat.CCL <- reshape2::dcast(tmp, 
                             formula = ID ~ season,
                             value.var = "value",
                             fun.aggregate = mean)
  
  out <- dat.CCL
  return(out)
  
}


compute.LOOIC <- function(loglik, data.vector, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  loglik.vec <- as.vector(loglik)
  loglik.mat <- matrix(loglik.vec[!is.na(data.vector)], 
                        nrow = MCMC.params$n.chains * n.per.chain)
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = 4)
  
  loo.out <- rstanarm::loo(loglik.mat, 
                           r_eff = Reff, 
                           cores = 4, k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}


# Extracting posterior samples of deviance or any other variable from jags output:
extract.samples <- function(varname, zm){
  dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
  return(dev)
}

