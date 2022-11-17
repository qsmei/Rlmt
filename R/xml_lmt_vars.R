xml_lmt_vars<-function(pars,  #include model_type and model_file 
                           data
						){

if(!is.null(pars$res_gamma$file)|!is.null(pars$res_gamma$value)){

	cat("Using user-prvoded  Gamma matrix to construct Residual variance structures! \n")
		# under development, this function was achieved by using res_block.csv 
		# vars_xml=rbind(vars_xml,"      <gamma>")
		# vars_xml=rbind(vars_xml,"        <D>")
		# vars_xml=rbind(vars_xml,paste0("          file: ",pars$res_gamma$file)) #only one pedigree exists in this term			
		# vars_xml=rbind(vars_xml,"        </D>")				
		# vars_xml=rbind(vars_xml,"      </gamma>")

	     #write  res.block.csv  for  handling inhomogeneous residual error 
		 phe_data=read.csv(data$phe$file,header=T)
		 if(!is.null(pars$res_gamma$value)){
		 phe_data$diag_res=pars$res_gamma$value[,1]
		 }else{
		 diag_res=read.csv(pars$res_gamma$file,header=T)
		 phe_data$diag_res=diag_res[,1]
		 }
		 if(!is.null(pars$vars$r$value)){
		 res_value=pars$vars$r$value
		 }else{
		 res_value=read.csv(pars$vars$r$file,header=F)[,]
		 }		 
		 
		 #code provided by Vinzent for writing 
		 CP <- list();RP <- list();SI <- list()
		 for(i in 1:nrow(phe_data)){
		  RP[[i]] <- as.integer(i)
		  CP[[i]] <- as.integer(1)
		  SI[[i]] <- matrix(1/res_value*1/phe_data$diag_res[i],1,1)
		 }
		 names(RP) <- paste("rp_",seq(1,length(RP)),sep="")
		 names(CP) <- paste("cp_",seq(1,length(CP)),sep="")
		 names(SI) <- paste("sigma_",seq(1,length(SI)),sep="")
		 O <- c(SI,CP,RP)
		 writelmtblockfile(O,"res.blkcsv","txt")
		 pars$vars$r$file=paste0(getwd(),"/res.blkcsv")	
}

vars=pars$vars 
vars_name=sapply(vars,function(x)x$name) #get the name of vars
vars=vars[order(match(vars_name,c("r","g",paste0("g",1:100),"p",setdiff(vars_name,c("r","g",paste0("g",1:100),"p")))))] #order the postion of vars 
vars_name=sapply(vars,function(x)x$name) #get the name of vars after sorting 
vars_xml="  <vars>" 

#res

if(!is.na(match("r",vars_name))){

vars_xml=rbind(vars_xml,"    <res>")
vars_xml=rbind(vars_xml,"      <sigma>")

vars_res=vars[[match("r",vars_name)]] 
	
if(is.null(pars$res_gamma$file)&is.null(pars$res_gamma$value)){ #consider the situation that inhomogeneous residual error 
	if(!is.null(vars_res$file)){
	vars_xml=rbind(vars_xml,paste0("        file:",vars_res$file))
	}else if(!is.null(vars_res$value)){
	vars_xml=rbind(vars_xml,"        <matrix attributes=\"array\">")	
	value_matrix=as.matrix(vars_res$value)
	vars_xml=rbind(vars_xml,matrix(apply(value_matrix,1,function(x)paste0("          ",paste(x,collapse = ","))),nrow=nrow(value_matrix)))	
	vars_xml=rbind(vars_xml,"        </matrix>")	
	}
	
}else{ 

	vars_xml=rbind(vars_xml,paste0("        type: ","block"))
	vars_xml=rbind(vars_xml,paste0("        file: ",vars_res$file))
}	
vars_xml=rbind(vars_xml,"      </sigma>")


vars_xml=rbind(vars_xml,"    </res>")
}

i=0
#vars, g,p......
if(length(setdiff(vars_name,"r"))!=0){

vars_xml=rbind(vars_xml,paste0("    vars: ",paste(setdiff(vars_name,"r"),collapse=",")))

	for(vars_i in vars[-1]){
		i=i+1
		if(gsub('[[:digit:]]+', '', vars_i$name)=="g"){  #animal effect 
		
			vars_xml=rbind(vars_xml,paste0("    <",vars_i$name,">"))				
			vars_xml=rbind(vars_xml,"      <sigma>")		

			if(!is.null(vars_i$value)){ #user provided value in R environment
			
			vars_xml=rbind(vars_xml,"        <matrix attributes=\"array\">")
			value_matrix=as.matrix(vars_i$value)
			vars_xml=rbind(vars_xml,matrix(apply(value_matrix,1,function(x)paste0("          ",paste(x,collapse = ","))),nrow=nrow(value_matrix)))				
			vars_xml=rbind(vars_xml,"        </matrix>")		


			}else if(!is.null(vars_i$file)){
				vars_xml=rbind(vars_xml,paste0("        file:",vars_i$file))
			}

			vars_xml=rbind(vars_xml,"      </sigma>")
		
			vars_xml=rbind(vars_xml,"      <gamma>")

			if(pars$blup_type=="PBLUP"){			
				vars_xml=rbind(vars_xml,"        <A>")			
				vars_xml=rbind(vars_xml,"          pedigree: myped1") #only one pedigree exists in this term
				vars_xml=rbind(vars_xml,"        </A>")					
			}else if(pars$blup_type%in%c("SS_TBLUP","SS_GBLUP","SS_SNPBLUP")){
			
				vars_xml=rbind(vars_xml,switch(pars$blup_type%in%c("SS_GBLUP","SS_TBLUP"),"        <H>",NULL))	#there is no need of <H> tag  for SS_SNPBLUP			
				vars_xml=rbind(vars_xml,na.omit(ifelse(pars$blup_type=="SS_TBLUP","          type: tblup",ifelse(pars$blup_type=="SS_SNPBLUP", "          type: snpblup1",NA))))
				vars_xml=rbind(vars_xml,switch(pars$blup_type%in%c("SS_GBLUP","SS_TBLUP"),"          pedigree: myped1",NULL)) #there is no need of pedigree in here for SS_SNPBLUP				
				
				vars_xml=rbind(vars_xml,ifelse(pars$blup_type%in%"SS_GBLUP",paste0("          grm: mygrm",i), paste0("          genotype: mygeno",i)))  #there is only need of grm  for SS_GBLUP
						                                                #not yet for multiple grms or multple genotypes in H	......					
				vars_xml=rbind(vars_xml,switch(!is.null(pars$aweight),paste0("          aweight:",pars$aweight), NULL)) 
				vars_xml=rbind(vars_xml,switch(!is.null(pars$switch),paste0("          switch:",pars$switch), NULL))

				vars_xml=rbind(vars_xml,switch(pars$blup_type%in%c("SS_GBLUP","SS_TBLUP"),"        </H>",NULL)) #there is no need of <H> tag  for SS_SNPBLUP											
			
			}else if(pars$blup_type%in%c("External_BLUP")){
			
				vars_xml=rbind(vars_xml,"        <E>")
				
				if(is.null(data$grm$file)){stop("User must provide grm_file!")}
				vars_xml=rbind(vars_xml,paste0("          file: ",data$grm$file))	
				
				vars_xml=rbind(vars_xml,"        </E>")		
				
			}else if(pars$blup_type%in%c("GBLUP")){ #
			
				vars_xml=rbind(vars_xml,"        <G>")
				
				if(is.null(data$grm$file)){stop("User must provide grm_file!")}
				#vars_xml=rbind(vars_xml,paste0("          file: ",data$grm$file))	
				vars_xml=rbind(vars_xml,paste0("      grm: ","mygrm1"))  		
				
				vars_xml=rbind(vars_xml,"        </G>")			
			}
			
			
			vars_xml=rbind(vars_xml,"      </gamma>")
			vars_xml=rbind(vars_xml,paste0("    </",vars_i$name,">"))
	}else{   # for permanent effect and non-genetic effect 
	
			vars_xml=rbind(vars_xml,paste0("    <",vars_i$name,">"))				
			vars_xml=rbind(vars_xml,"      <sigma>")		

			if(!is.null(vars_i$value)){ #user provided value in R environment

			vars_xml=rbind(vars_xml,"        <matrix attributes=\"matrix\">")
			value_matrix=as.matrix(vars_i$value)
			vars_xml=rbind(vars_xml,matrix(apply(value_matrix,1,function(x)paste0("          ",paste(x,collapse = ","))),nrow=nrow(value_matrix)))				
			vars_xml=rbind(vars_xml,"        </matrix>")
			
			}else if(!is.null(vars_i$file)){
				vars_xml=rbind(vars_xml,paste0("        file:",vars_i$file))
			}

			vars_xml=rbind(vars_xml,"      </sigma>")
			vars_xml=rbind(vars_xml,paste0("    </",vars_i$name,">"))			
	
	}

	}
vars_xml=rbind(vars_xml,"  </vars>")		

}
return(vars_xml)
}


#assign initial vars 
initial_vars<-function(random_type  #list of lmt_models object
					   ){

types=data.frame(table(random_type))
vars=NULL
for(i in 1:nrow(types)){
vars_name=as.character(types[i,1])
vars_value=matrix(0.5,nrow=types[i,2],ncol=types[i,2])
diag(vars_value)=1
vars=c(vars,lmt_vars$new(name=vars_name,value=vars_value))
}
names(vars)=sapply(vars,function(x)x$name)
return(vars)
}



