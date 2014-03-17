calc_relative_humidity_in_soil <- function(soil_moisture, layers){
  
  len_soil_moisture <- length(soil_moisture)/layers
  rel_humidity = array(NA, dim=c(layers,len_soil_moisture))
  fract_wilting
  for (L in 1:layers){
    for (LM in 1:len_soil_moisture){
      evap_stop <- min(0.1, max(soil_moisture[L,]))
	  moisture_limit =  max(soil_moisture[L,]) - evap_stop
	
	  if (soil_moisture[L,LM] > moisture_limit & soil_moisture[L,LM] > fract_wilting*max(soil_moisture[L,])) {
	     rel_humidity[L,LM] = 0.5*(1-cos((soil_moisture[L,LM]-moisture_limit)*pi/evap_stop))   
	  }
	  else {
	    rel_humidity[L,LM] = 0   
	  }
	}  
  }
}  
