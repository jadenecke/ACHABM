# Custom functions
rotateCoordinatesByY <- function(v, d){
  matrix(c(cospi(d/180), 0, sinpi(d/180), 0, 1, 0, -sinpi(d/180), 0, cospi(d/180)), ncol = 3, nrow = 3, byrow = TRUE) %*% v
}
rotateCoordinatesByZ <- function(v, d){
  matrix(c(cospi(d/180), -sinpi(d/180), 0, sinpi(d/180), cospi(d/180), 0, 0, 0, 1), ncol = 3, nrow = 3, byrow = TRUE) %*% v
}
rotateCoordinatesByX <- function(v, d){
  matrix(c(1, 0, 0, 0, cospi(d/180), -sinpi(d/180), 0, sinpi(d/180), cospi(d/180)), ncol = 3, nrow = 3, byrow = TRUE) %*% v
}
resetAMonomers <- function(...){
  for(i in seq_along(neurons)){
    neurons[[i]][["aMonomer"]] <<- matrix(rep(0, length(neurons[[i]][["aMonomer"]])), nrow = nrow(neurons[[i]][["aMonomer"]]), ncol = ncol(neurons[[i]][["aMonomer"]]))
  }
}