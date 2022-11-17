xml_lmt_data<-function(data, #lmt_data object
					   pars   #lmt_par object
                           ){

 
data$write_file()  #write file and get the file name

#for phenotype data
data_xml="  <data>" 
data_xml=rbind(data_xml,paste0("    datafile:",data$phe$file))
data_xml=rbind(data_xml,paste0("    	missingthreshold:",data$phe$missing_threshold))
data_xml=rbind(data_xml,"  </data>") #not yet for missing phenotype file 

#for pedigree 
if(!is.null(data$ped$file)){  #consider the situation that pedigree is not need, e.g. GBLUP......

	for(i in 1:length(data$ped$file)){

		data_xml=rbind(data_xml,"  <pedigrees>") 
		data_xml=rbind(data_xml,paste0("    pedigrees: myped",i)) 
		data_xml=rbind(data_xml,paste0("    <myped",i,">"))
		data_xml=rbind(data_xml,paste0("      file:",data$ped$file))
		
		if(!is.null(pars$meta_gamma)){

			data_xml=rbind(data_xml,"        <metamatrix attributes=\"array\">")
			value_matrix=as.matrix(pars$meta_gamma)
			data_xml=rbind(data_xml,matrix(apply(value_matrix,1,function(x)paste0("          ",paste(x,collapse = ","))),nrow=nrow(value_matrix)))				
			data_xml=rbind(data_xml,"        </metamatrix>")		
		
			if(is.null(data$ped$switch)){data$ped$switch="selfing"}
		}
		
		data_xml=rbind(data_xml,switch(!is.null(data$ped$switch),paste0("      switch:",data$ped$switch), NULL)) 		
		data_xml=rbind(data_xml,paste0("    </myped",i,">"))
		data_xml=rbind(data_xml,"  </pedigrees>") 
	
	}

}

if(pars$blup_type!="Defined_BLUP"|!is.null(data$grm$file)){
#for genotype  
if(!is.null(data$geno$file)){  #consider the situation that pedigree is not need, e.g. GBLUP......

	data_xml=rbind(data_xml,"  <genotypes>") 
	data_xml=rbind(data_xml,paste0("      genotypes: ",paste("mygeno",1:length(data$geno$file),sep="",collapse =",")))
	for(i in 1:length(data$geno$file)){
		data_xml=rbind(data_xml,paste0("    <mygeno",i,">"))
		data_xml=rbind(data_xml,paste0("      file:",data$geno$file[i]))  
		data_xml=rbind(data_xml,paste0("      cross:",data$geno$id_file[i]))
		data_xml=rbind(data_xml,paste0("      pedigree: myped",1))		
		data_xml=rbind(data_xml,paste0("    </mygeno",i,">"))	
	}
	data_xml=rbind(data_xml,"  </genotypes>") 
}
}
#for grm  
if(!is.null(data$grm$file)){  #using provided grm for external analyasis

	data_xml=rbind(data_xml,"  <grms>") 
	data_xml=rbind(data_xml,paste0("      grms: ",paste("mygrm",1:length(data$grm$file),sep="",collapse =",")))
	for(i in 1:length(data$grm$file)){  #not yet ......
	
			data_xml=rbind(data_xml,paste0("    <mygrm",i,">"))
			data_xml=rbind(data_xml,paste0("      file: ",data$grm$file[i]))  	
			data_xml=rbind(data_xml,switch(!is.null(data$grm$id_file),paste0("      cross:",data$grm$id_file[i]),NULL)) 	
					
			data_xml=rbind(data_xml,paste0("    </mygrm",i,">"))
	
	}
	data_xml=rbind(data_xml,"  </grms>")

}else if(!is.null(data$geno$file)){  #seems only for SS-GBLUP situation

	if(pars$blup_type%in%"SS_GBLUP"){
	
		data_xml=rbind(data_xml,"  <grms>") 
		data_xml=rbind(data_xml,paste0("      grms: ",paste("mygrm",1:length(data$geno$file),sep="",collapse =",")))
	
		for(i in 1:length(data$geno$file)){

			data_xml=rbind(data_xml,paste0("    <mygrm",i,">"))
			data_xml=rbind(data_xml,paste0("      genotype: mygeno",i))  	
			#not solved yet
			#data_xml=rbind(data_xml,switch(!is.null(pars$grm_method),paste0("      method:",pars$grm_method), NULL)) 		
			data_xml=rbind(data_xml,paste0("    </mygrm",i,">"))
			 
	
		}
		data_xml=rbind(data_xml,"  </grms>")
	}

}

return(data_xml)
}