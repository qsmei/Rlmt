#jobs-xml 
xml_lmt_jobs<-function(jobs   #  lmt_jobs object  
				        ){

	
if(jobs$type=="mcemreml"&is.null(jobs$samplers$type)){	
	jobs$samplers$type="singlepass"	
}

if(jobs$type=="pevsolve"&is.null(jobs$pevsolve_factor)){	
	jobs$pevsolve_factor="g"	
}
 
jobs_xml="  <jobs>" 
jobs_xml=rbind(jobs_xml,paste0("    jobs: ",paste(jobs$type,collapse=",")))		

for(i_job in jobs$type){

	if(i_job=="default"){
	
		jobs_xml=rbind(jobs_xml,"    <default>")	
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$conv),paste0("      conv:",jobs$conv), NULL)) #ifelse can't return NULL value		
		jobs_xml=rbind(jobs_xml,"    </default>")			
		
	}else if(i_job=="solve"){
	
		jobs_xml=rbind(jobs_xml,"    <solve>")	
		jobs_xml=rbind(jobs_xml,paste0("      solver: ",paste("mysolver",1:length(jobs$solvers$type),sep="",collapse =",")))			
		jobs_xml=rbind(jobs_xml,"    </solve>")	
		
	}else if(i_job=="sample"){
	
		jobs_xml=rbind(jobs_xml,"    <sample>")	
		jobs_xml=rbind(jobs_xml,paste0("      sampler: ",paste("mysampler",1:length(jobs$samplers$type),sep="",collapse =",")))			
		jobs_xml=rbind(jobs_xml,"    </sample>")			

	}else if(i_job=="pevsample"){
	
		jobs_xml=rbind(jobs_xml,"    <pevsample>")	
		jobs_xml=rbind(jobs_xml,paste0("      sampler: ",paste("mysampler",1:length(jobs$samplers$type),sep="",collapse =",")))
		jobs_xml=rbind(jobs_xml,"    </pevsample>")	

	}else if(i_job=="pevsolve"){
	
		jobs_xml=rbind(jobs_xml,"    <pevsolve>")	
		jobs_xml=rbind(jobs_xml,paste0("      solver: ",paste("mysolver",1:length(jobs$solvers$type),sep="",collapse =",")))		
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$pevsolve_factor),paste0("      factor:",jobs$pevsolve_factor), NULL)) 
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$pevsolve_levels),paste0("      levels:",jobs$pevsolve_levels), NULL)) 
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$pevsolve_levelsfile),paste0("      levelfile:",jobs$pevsolve_levelsfile), NULL)) 		
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$pevsolve_nrhs),paste0("      nrhs:",jobs$pevsolve_nrhs), NULL)) 
		jobs_xml=rbind(jobs_xml,"    </pevsolve>")			

	}else if(i_job=="aireml"){
	
		jobs_xml=rbind(jobs_xml,"    <aireml>")	
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$rounds),paste0("      rounds:",jobs$rounds), NULL)) 
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$convtype),paste0("      convtype:",jobs$convtype), NULL)) 		
		jobs_xml=rbind(jobs_xml,switch(sum(jobs$convtype%in%c("ll","ng","cd"))>=1,paste0("      convtype",jobs$convtype,":",jobs$conv), NULL)) 
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$nscale),paste0("      nscale:",jobs$nscale), NULL)) 		
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$switch),paste0("      switch:",jobs$switch), NULL)) 

		jobs_xml=rbind(jobs_xml,"    </aireml>")		

	}else if(i_job=="airemlc"){
	
		jobs_xml=rbind(jobs_xml,"    <airemlc>")	
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$rounds),paste0("      rounds:",jobs$rounds), NULL)) 
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$convtype),paste0("      convtype:",jobs$convtype), NULL)) 		
		jobs_xml=rbind(jobs_xml,switch(sum(jobs$convtype%in%c("ll","ng","cd"))>=1,paste0("      convtype",jobs$convtype,":",jobs$conv), NULL)) 
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$nscale),paste0("      nscale:",jobs$nscale), NULL)) 		
		jobs_xml=rbind(jobs_xml,ifelse(!is.null(jobs$switch),paste0("      switch:",jobs$switch), "      switch:writeai,residuals,solutions")) 

		jobs_xml=rbind(jobs_xml,"    </airemlc>")		

	}else if(i_job=="mcemreml"){
	
		jobs_xml=rbind(jobs_xml,"    <mcemreml>")	
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$rounds),paste0("      rounds:",jobs$rounds), NULL)) 
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$conv),paste0("      conv:",jobs$conv), NULL)) #ifelse can't return NULL value			
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$rounds),paste0("      sampler:mysampler1"), NULL)) 
		jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$rounds),paste0("      solver:mysolver1"), NULL)) 		
		jobs_xml=rbind(jobs_xml,"    </mcemreml>")				

	}else if(i_job=="yhat"){
	
		jobs_xml=rbind(jobs_xml,"    <yhat>")	
		jobs_xml=rbind(jobs_xml,"    </yhat>")	
	}

}
jobs_xml=rbind(jobs_xml,"  </jobs>")	
#for samplers
if(!is.null(jobs$samplers$type)){
	
		jobs_xml=rbind(jobs_xml,"  <samplers>")			
		jobs_xml=rbind(jobs_xml,paste0("   samplers: ",paste("mysampler",1:length(jobs$samplers$type),sep="",collapse =",")))			

		for(i_sampler in 1:length(jobs$samplers$type)){
			
			jobs_xml=rbind(jobs_xml,paste0("     <",paste0("mysampler",i_sampler,">")))							
			jobs_xml=rbind(jobs_xml,paste0("       <",jobs$samplers$type[i_sampler],">"))				
			jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$samplers$samples),paste0("    samples:",jobs$samplers$samples), NULL)) 			
			jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$samplers$burnin),paste0("    rounds:",jobs$samplers$burnin), NULL)) 		
			jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$samplers$chains),paste0("    rounds:",jobs$samplers$chains), NULL)) 			
			jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$samplers$switch),paste0("    rounds:",jobs$samplers$switch), NULL)) 
			jobs_xml=rbind(jobs_xml,paste0("       </",jobs$samplers$type[i_sampler],">"))	
			jobs_xml=rbind(jobs_xml,paste0("     </",paste0("mysampler",i_sampler,">")))			
		}
		
		jobs_xml=rbind(jobs_xml,"  </samplers>")	
}


#for solvers
if(!is.null(jobs$solvers$type)){
	
		jobs_xml=rbind(jobs_xml,"  <solvers>")			
		jobs_xml=rbind(jobs_xml,paste0("   solvers: ",paste("mysolver",1:length(jobs$solvers$type),sep="",collapse =",")))			

		for(i_solver in 1:length(jobs$solvers$type)){
			
			jobs_xml=rbind(jobs_xml,paste0("     <",paste0("mysolver",i_solver,">")))							
			jobs_xml=rbind(jobs_xml,paste0("       <",jobs$solvers$type[i_solver],">"))
			
			if(jobs$solvers$type[i_solver]=="pcgiod"){
			
				jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$solvers$rounds),paste0("         rounds:",jobs$solvers$rounds), NULL)) 			
				jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$solvers$conv),paste0("        conv:",jobs$solvers$conv), NULL)) 		
				jobs_xml=rbind(jobs_xml,switch(!is.null(jobs$solvers$convtype),paste0("        convtype:",jobs$solvers$convtype), NULL)) 
			}
			jobs_xml=rbind(jobs_xml,paste0("       </",jobs$solvers$type[i_solver],">"))	
			jobs_xml=rbind(jobs_xml,paste0("     </",paste0("mysolver",i_solver,">")))			
		}
		
		jobs_xml=rbind(jobs_xml,"  </solvers>")	
}	


return(jobs_xml)
}
		
