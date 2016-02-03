create_inputs <- function(simdir, param, val, lh_dat){
	
	if(simdir != "sensitivity"){
		dat_input <- lh_dat
		dat_input$SigmaR <- 0.6
		dat_input$SigmaF <- 0.2
		dat_input$CV_c <- lh_dat$CVcatch
		dat_input$CV_l <- lh_dat$CVlen
	}

	## adjust DataList
	if(grepl("sensitivity", simdir)){
		dat_input <- lh_dat
		dat_input$SigmaR <- 0.6
		dat_input$SigmaF <- 0.2
		dat_input$CV_c <- lh_dat$CVcatch
		dat_input$CV_l <- lh_dat$CVlen
		dat_input[[param]] <- val
	}

	return(dat_input)
}