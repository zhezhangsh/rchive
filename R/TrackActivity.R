TrackActivity <- function(session, action = '', data = NULL, dir = '.') {
  cdata <- session$clientData; 
  ddata <- lapply(names(cdata), function(name) cdata[[name]]);
  names(ddata) <- names(cdata);
  
  fn <- paste(dir, '/', as.character(Sys.Date()), '_', as.integer(Sys.time()), '_', session$token, '.rds', sep='');
  
  saveRDS(list(clientData = ddata, action = action, data = data), fn); 
}