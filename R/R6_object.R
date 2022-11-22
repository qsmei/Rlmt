library(R6) 
lmt_vars <- R6Class("lmt_vars",
  public = list(
	value=NULL,
	name="g",
	file=NULL,

initialize=function(value=NULL,
					 name=NULL,
					 file=NULL
					 ){
		
	self$value=value
	self$name=name 
	self$file=file

		  }
)
)

lmt_formulas <- R6Class("lmt_formulas",
  public = list(
	fixed=NULL,
	covariate=NULL,
	random=NULL,
	polyno=NULL,
initialize=function(fixed=NULL,
					covariate=NULL,
					random=NULL,
					polyno=NULL){
		
	self$fixed=fixed
	self$covariate=covariate 
	self$random=random
	self$polyno=polyno	

		  }),
active=list(

	trait=function(fixed=self$fixed)all.vars(fixed)[1]  #dynamic changing with the value of self$fixed
														#user can't modify trait directly
		   )		  
)



lmt_jobs <- R6Class("lmt_jobs",
  public = list(
	type=NULL,  #default,solve,sample,pevsample,pevsolve,airemlc,mcemreml,yhat
	conv=NULL,
	pevsolve_factor=NULL,
	pevsolve_levels=NULL,
	pevsolve_levelsfile=NULL,
	pevsolve_nrhs=NULL,
	rounds=NULL,
	convtype=NULL,
	nscale=NULL,
	switch=NULL,
	samplers=NULL, #list
	solvers=NULL,
initialize=function(type="default",
					conv=NULL,
					pevsolve_factor=NULL,
					pevsolve_levels=NULL,
					pevsolve_levelsfile=NULL,
					pevsolve_nrhs=NULL,
					rounds=NULL,
					convtype=NULL,
					nscale=NULL,
					switch=NULL,
					samplers_type=NULL,  # singlepass, blocked, pev
					samplers_samples=NULL,
					samplers_burnin=NULL,
					samplers_chains=NULL,
					samplers_switch=NULL,
					samplers_conv=NULL,
					samplers_convtype=NULL,
					solvers_type=c("pcgiod"),  #pcgiod, direct 
					solvers_conv=NULL,
					solvers_convtype=NULL,
					solvers_rounds=NULL
					){
		
	self$type=type
	self$conv=conv	
	self$convtype=convtype	
	self$rounds=rounds	
	self$samplers$type=samplers_type  # singlepass, blocked, pev
	self$samplers$samples=samplers_samples
	self$samplers$burnin=samplers_burnin
	self$samplers$chains=samplers_chains
	self$samplers$switch=samplers_switch
	self$samplers$conv=samplers_conv
	self$samplers$convtype=samplers_convtype	
	self$solvers$type=solvers_type	
	self$solvers$conv=solvers_conv
	self$solvers$convtype=solvers_convtype
	self$solvers$rounds=solvers_rounds

	}
)
)



lmt_pars <- R6Class("lmt_pars",
  public = list(
	blup_type=NULL,  #PBLUP,GBLUP,SNPBLUP,SS-SNPBLUP,SS-GBLUP,SS-TBLUP,Defined-BLUP(sparase and dense)
					  #for Defined-BLUP,see /usr/home/qgg/vinzent/lmt/examples/blup/multivariate_special/gamma_external_dense
	switch=NULL, #adjustg2a
	aweight=NULL,
	meta_gamma=NULL,
	grm_method=NULL,
	gg_type=NULL, # absorb, extra
	jobs=NULL,
	vars=NULL,
	res_gamma=NULL, #
initialize=function(blup_type="PBLUP",
					aweight=0.05,
					switch="adjustg2a",
					meta_gamma=NULL,
					gg_type=NULL,
					grm_method="VR",  # VR, YA
					jobs=lmt_jobs$new(),
					vars=lmt_vars$new(),
					res_gamma=NULL,
					res_gamma_file=NULL
					){
		
	self$blup_type=blup_type
	self$aweight=aweight 
	self$switch=switch
	self$gg_type=gg_type
	self$jobs=jobs
	self$grm_method=grm_method
	self$res_gamma$value=res_gamma
	self$res_gamma$file=res_gamma_file		
	if("lmt_vars"%in%class(vars)&is.null(vars$name)&is.null(vars$name)&is.null(vars$name)){  #in the situation that vars was null object	
		self$vars=vars 	
	}else{
		self$vars=c(NULL,vars)
		names(self$vars)=sapply(self$vars, function(x)x$name)
	}
		  },

	 add_vars=function(value=NULL,
					   name=NULL,
					   file=NULL){

	 if("lmt_vars"%in%class(value)){
		i_vars=value         # R6 didin't support function polymorphism, thus we used the first parameter as the default parameter
	 }else{
		i_vars=lmt_vars$new(value,name,file)
	 }
	 
	 if("lmt_vars"%in%class(self$vars)&is.null(self$vars$name)){ #consider that vars is null lmt_vars object
	 
		self$vars=c(NULL,i_vars)
	 
	 }else{

		vars_name=sapply(self$vars,function(x)x$name) 	
		
	 
		 if(i_vars$name%in%vars_name){ 
			cat(paste0("vars:",i_vars$name," was already exists, the old value will be coverd by the new value!","\n"))
			self$vars[[match(i_vars$name,vars_name)]]=i_vars
		 }else{
			self$vars=c(self$vars,i_vars)
		 }
	}	 
	     names(self$vars)=sapply(self$vars, function(x)x$name)	
		
		 invisible(self)
	 },
	 
	 rm_vars=function(name=NULL){
	 
		vars_name=sapply(self$vars,function(x)x$name)	
	 
		 if(name%in%traits_name){ 
			self$vars=self$formulas[-c(match(name,traits_name))]
		 }else{
		 	cat(paste0("vars:",name," is not exists in vars,please check your input!","\n"))
		 }	
		names(self$vars)=sapply(self$vars, function(x)x$name)
		invisible(self)
	 }		  
		  
)
)

# how to handle two genetic effect (g1,g2) by providing two genotype data, e.g. dominance effect
# ssSNPBLUP model with a separate polygenic factor

lmt_models<-R6Class("lmt_models",
    public = list(
    formulas=NULL,
    pars=NULL,
    id_name=NULL,
	dam_name=NULL,
	pe_name=NULL,
initialize=function(fixed=NULL,
					covariate=NULL,
					random=NULL,
					polyno=NULL,
					blup_model=NULL,
					id_name="id",
					dam_name="dam",
					pe_name=NULL,
					formulas=lmt_formulas$new(), 
					pars=lmt_pars$new()){
					

	self$pars=pars
	self$id_name=id_name
	self$dam_name=dam_name	

	if("lmt_formulas"%in%class(formulas)&is.null(formulas$fixed)&is.null(fixed)){ # for null lmt_formulas object  
	
		self$formulas=formulas
	
	}else{ 

		if(!is.null(fixed)){
			formulas$fixed=fixed
			formulas$covariate=covariate 
			formulas$random=random
			formulas$polyno=polyno
		}
		self$formulas=c(NULL,formulas)	#make sure formulas is a list 

		names(self$formulas)=sapply(self$formulas, function(x)x$trait)
	}		
	},
		  
	 add_formulas=function(fixed=NULL,
						   covariate=NULL,
					        random=NULL,
					        polyno=NULL){
	 
	 if("lmt_formulas"%in%class(fixed)){
		i_formula=fixed         # R6 didin't support function polymorphism, thus we used the first parameter as the default parameter
	 }else{
		i_formula=lmt_formulas$new(fixed,covariate,random,polyno)
	 }

	if("lmt_formulas"%in%class(self$formulas)&is.null(self$formulas$fixed)){ #consider that vars is null lmt_vars object
	 
		self$formulas=c(NULL,i_formula)
	 
	 }else{

		traits_name=sapply(self$formulas,function(x)x$trait) 	
	 
		 if(i_formula$trait%in%traits_name){ 
			cat(paste0("Trait:",i_formula$trait," was already exists in formulas, the old one will be coverd by the new one!","\n"))
			self$formulas[[match(i_formula$trait,traits_name)]]=i_formula
		 }else{
			self$formulas=c(self$formulas,i_formula)
		 }	
	}	
		names(self$formulas)=sapply(self$formulas, function(x)x$trait)
	
	    invisible(self)
	} ,

	 rm_formulas=function(trait=NULL){
	 
		traits_name=sapply(self$formulas,function(x)x$trait) 	
	 
		 if(trait%in%traits_name){ 
			self$formulas=self$formulas[-c(match(trait,traits_name))]
		 }else{
		 	cat(paste0("Trait:",trait," is not exists in formulas,please check your input!","\n"))
		 }	
		names(self$formulas)=sapply(self$formulas, function(x)x$trait)
		invisible(self)
	 }	
	 
)
)

lmt_data<-R6Class("lmt_data",
    public = list(
    phe=NULL,
    ped=NULL,
    geno=NULL,
    output_path=NULL,
	grm=NULL,
initialize=function(phe=NULL,
					phe_file=NULL,
					phe_missing_indicator=NULL,
					phe_missing_indicator_file=NULL,
					ped_file=NULL,  #for multiple files, it should keep normal pedigree in the first position
					phe_missing_threshold=-9999,
					ped=NULL,
					ped_switch=NULL,
					ped_phanto=NULL,
					ped_qfile=NULL,
					geno=NULL,
					geno_file=NULL,
					geno_id=NULL,
					geno_id_file=NULL,
					geno_switch=NULL,
					geno_pq=NULL,
					geno_pq_file=NULL,
					grm=NULL,
					grm_file=NULL,
					grm_id=NULL,
					grm_id_file=NULL,
					#grm_method=NULL,
					grm_outfile=NULL,
					output_path=NULL){
	self$phe$value=phe
	self$phe$file=phe_file
	self$phe$missing_indicator=phe_missing_indicator
	self$phe$missing_indicator_file=phe_missing_indicator_file	
	self$phe$missing_threshold=phe_missing_threshold
	self$ped$value=ped
	self$ped$file=ped_file
	self$ped$switch=ped_switch
	self$ped$phanto=ped_phanto	
	self$ped$qfile=ped_qfile		
	self$geno$value=geno
	self$geno$file=geno_file
	self$geno$id=geno_id
	self$geno$id_file=geno_id_file
	self$geno$pq=geno_pq 
	self$geno$pq_file=geno_pq_file
	self$geno$switch=geno_switch
	self$grm$value=grm
	self$grm$file=grm_file
	self$grm$id=grm_id 
	self$grm$id_file=grm_id_file
	#self$grm$method=grm_method						
	self$output_path=output_path

		  }, 
	 #writing data as local file
	 
     write_file=function(){
	 
	 if(is.null(self$phe$file)&!is.null(self$phe$value)){	 #phenotype
	 
		if(substr(colnames(self$phe$value)[1],1,1)!="#"){
			colnames(self$phe$value)[1]=paste0("#",colnames(local_phe)[1],collapse="") # the requiremnt of the first column name of phe data
		}
	 
		write.csv(self$phe$value,paste0(self$output_path,"/lmt_pheno.csv"),row.names=F)
		self$phe$file=paste0(self$output_path,"/lmt_pheno.csv")
	 }
	 
	 if(is.null(self$ped$file)&!is.null(self$ped$value)){	 #pedigree	
	 
		for(i in 1:length(self$ped$value)){  #many pedigree file may exist		
			write.csv(self$ped$value[i],paste0(self$output_path,"/lmt_ped",i,".csv"),row.names=F)
			self$ped$file=c(self$ped$file,paste0(self$output_path,"/lmt_ped",i,".csv"))			
		}	 
	 }
	 
	 if(is.null(self$geno$file)&!is.null(self$geno$value)){  #genotype	
	 
		for(i in 1:length(self$geno$value)){  #many pedigree file may exist		
			write.csv(self$geno$value[i],paste0(self$output_path,"/lmt_geno",i,".csv"),row.names=F)
			self$geno$file=c(self$geno$file,paste0(self$output_path,"/lmt_geno",i,".csv"))
		}
	 }

	 if(is.null(self$geno$id_file)&!is.null(self$geno$id)){	 #id of genotype
	 
		for(i in 1:length(self$geno$id)){  #many pedigree file may exist		
			write.csv(self$geno$id[i],paste0(self$output_path,"/lmt_genoid",i,".csv"),row.names=F)  #id[[i]], maybe be a possible bug in the future version of R
			self$geno$id_file=c(self$geno$id_file,paste0(self$output_path,"/lmt_genoid",i,".csv"))
		}
	 }	 

	#not yet for grm ......
	 
	 }	 
)
)

lmt<-R6Class("lmt",
    #inherit = lmt_models,
    public = list(
	models=NULL,
	data=NULL,
	xml_file=NULL,
	software_path=NULL,
	software_name=NULL,
	output_path=NULL,
	show_info=TRUE,
initialize=function(models=lmt_models$new(),
					data=lmt_data$new(),
					xml_file=NULL,
					software_path=NULL,
					software_name="lmt",
					output_path=NULL,
					show_info=TRUE
					){
		
	self$data=data
	self$models=models
	self$xml_file=xml_file
	self$software_path=software_path
	self$software_name=software_name
	self$output_path=ifelse(is.null(output_path),getwd(),output_path)
		  },
   
   #running lmt 
   run_lmt=function(output_path=NULL){
   
    trait_name=sapply(self$models$formulas,function(x)x$trait)
   
	cat("R package:Rlmt is only the wrapper of LMT in the field of academic research! \n")
	cat("For commercial use of the LMT, please contact Vinzent Boerner of QGG. Email: vinzent.boerner@gmx.de,vinzent.boerner@mbg.au.dk  !\n")

	if(length(grep(" ",output_path))>0){stop("output_path shouldn't contain blank, please choose another path!")}
	
	if(is.null(output_path)){output_path=self$output_path}
	
	output_path=paste0(output_path,"/",self$models$pars$blup_type,"_",paste(trait_name,collapse = "_"))	
	if(!file.exists(output_path)){dir.create(output_path,recursive=TRUE)}
	setwd(output_path)

	if(!is.null(self$xml_file)){
	
		system(paste0("cp -r ",self$xml_file," ",output_path,"/",paste(trait_name,collapse = "_"),"_lmt_par.xml"))
	
	}else{
	
		write.table(mylmt$xml,paste0(output_path,"/",paste(trait_name,collapse = "_"),"_lmt_par.xml"),quote=F,row.names=F,col.names=F)
		
	}
	
	cat(paste0("Results of LMT are saved in path: ",output_path,"\n"))
	cat(paste0("Start the ",self$models$pars$blup_type," analyse of ",length(trait_name),
	    ifelse(length(trait_name)==1," trait"," traits")," model:",paste(trait_name,collapse = " & ")," \n"))

	system(paste0("cp -r ",self$software_path,"/",self$software_name," ",output_path))	
	system(paste0("cp -r ",self$software_path,"/",self$software_name,".key ",output_path))
	# system("ulimit -s unlimited")
	# system("export KMP_AFFINITY=granularity=core,scatter")
	# system("export OMP_NUM_THREADS=32")
	# system("export OMP_STACKSIZE=2000M")
	# system("export OMP_DYNAMIC=false")
	# system("export OMP_PLACES=cores")
	# system("export OMP_PROC_BIND=true")
	# system("export OMP_MAX_ACTIVE_LEVELS=2147483647")
	if(self$show_info==TRUE){
		system(paste0(self$software_path,"/",self$software_name," -f ",paste0(output_path,"/",paste(trait_name,collapse = "_"),"_lmt_par.xml")),ignore.stdout=F)					
	}else{
		system(paste0(self$software_path,"/",self$software_name," -f ",paste0(output_path,"/",paste(trait_name,collapse = "_"),"_lmt_par.xml")),ignore.stdout=T)
	}
	
	cat(paste0("Complete the ",self$models$pars$blup_type," analyse of ",length(trait_name),
	    ifelse(length(trait_name)==1," trait"," traits")," model:",paste(trait_name,collapse = " & ")," \n"))						
					   
    
   },
   
    #show in the console
    print = function(...) {
      cat("<Rlmt> of ",2, " elements\n", sep = "")
    }
   ),
   
   active=list(
   
   xml=function(models=self$models,data=self$data){
   
		if(!is.null(data$phe$file)&!is.null(c(NULL,models$formulas)[[1]]$fixed)){
			
			#for data 
			data_xml=xml_lmt_data(data,models$pars)
		     
			#for models
			tmp_xmls=xml_lmt_models(models)
			random_type=tmp_xmls$type
			models_xml=tmp_xmls$xml

			#initial vars 
			lmt_initial_vars=initial_vars(random_type)
			
			if("lmt_vars"%in%class(models$pars$vars)&is.null(models$pars$vars$name)){ #for null lmt_vars object
				
				cat("User doesn't provide intial values for variance components, software will assign initial values for all random effects in formulas automaticaly!\n")
				models$pars$vars=lmt_initial_vars
			
			}else {
				if(!is.null(models$pars$vars))names(models$pars$vars)=sapply(models$pars$vars,function(x)x$name)
				for(initial_i_vars in lmt_initial_vars){
					
					name=initial_i_vars$name
					
					i_vars=models$pars$vars[[match(name,names(models$pars$vars))]]
					if(is.null(i_vars)|!identical(nrow(matrix(i_vars$value)),nrow(matrix(initial_i_vars$value)))){					
						cat("User-provided variance components are not compatiable with the effects in formulas, software will assign initial values for all effects in formulas automaticaly!\n")
						models$pars$vars=lmt_initial_vars
						break;
					}
				
				}
				
			}
			#for vars 
			vars_xml=xml_lmt_vars(models$pars,data)
			
			#for jobs
			jobs_xml=xml_lmt_jobs(models$pars$jobs)
			
			xml="<root>"
			xml=rbind(xml,data_xml,vars_xml,models_xml,jobs_xml,"</root>")
			return(xml)
		}
   
   
         }
   )	  
)

