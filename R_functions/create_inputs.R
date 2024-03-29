create_inputs <- function(param, val, lh_dat){
	
	if(any(param==FALSE) & any(val==FALSE)){
		dat_input <- lh_dat
		dat_input$CV_c <- lh_dat$CVcatch
		dat_input$CV_l <- lh_dat$CVlen
	}

	## adjust DataList
	if(any(param!=FALSE) & any(val!=FALSE)){
		dat_input <- lh_dat
		dat_input$CV_c <- lh_dat$CVcatch
		dat_input$CV_l <- lh_dat$CVlen
		
		for(i in 1:length(param)){
			dat_input[[param[i]]] <- val[i]
		}
	}

	return(dat_input)
}