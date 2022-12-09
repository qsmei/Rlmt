get_models<-function(models
					){

 
fixed=models$fixed      #include trait and fixed(class effect)
covariate=models$covariate
random=models$random
polyno=models$polyno   		#days~{log(sqrt(x))}+{3*exp(x)}+{x^2}
						    #assume all traits have the same polynominal expression

#for fixed term 
if(is.null(fixed)){stop("Please check your input of the fixed term,it should not be empty!")}
trait=all.vars(fixed)[1]
tmp_fixed=get_nested(fixed)
fixed_effect=tmp_fixed$effect_name
in_nested_fixed_effect=tmp_fixed$in_nested_effect
out_nested_fixed_effect=tmp_fixed$out_nested_effect
#fixed_effect=setdiff(fixed_effect,in_nested_fixed_effect) #without nested random effect

if(!is.null(in_nested_fixed_effect)){stop("Please check your input of the fixed term,there should be no nested structure!")}

#for covariate
covariate_effect=NULL
in_nested_covariate_effect=NULL
out_nested_covariate_effect=NULL
if(!is.null(covariate)){
tmp_covariate=get_nested(covariate)
covariate_effect=tmp_covariate$effect_name
in_nested_covariate_effect=tmp_covariate$in_nested_effect
out_nested_covariate_effect=tmp_covariate$out_nested_effect
covariate_effect=setdiff(covariate_effect,in_nested_covariate_effect) #without nested random effect
}
#for random 
random_effect=NULL
in_nested_random_effect=NULL
out_nested_random_effect=NULL
if(!is.null(random)){
tmp_random=get_nested(random)
random_effect=tmp_random$effect_name
in_nested_random_effect=tmp_random$in_nested_effect
out_nested_random_effect=tmp_random$out_nested_effect
random_effect=setdiff(random_effect,in_nested_random_effect) #without nested random effect
}
#if(is.null(tmp_random$effect_name)){stop("Please check your input of the random term,there should be random effect exists!")}

#for polynominal
poly_effect=NULL 
poly_expression=NULL
if(!is.null(polyno)){
poly_effect=all.vars(polyno)[1]
poly_expression=as.character(polyno)[-c(1:2)]	
poly_expression=gsub("\n","",poly_expression)
poly_expression=gsub(" ","",poly_expression)
poly_expression=unlist(regmatches(poly_expression, gregexpr("(?<=\\{).+?(?=\\})", poly_expression,perl=T))) #using regrex to distinguish multple polynominal
}

return(list(trait=trait,
			fixed_effect=fixed_effect,
			in_nested_fixed_effect=in_nested_fixed_effect,
			out_nested_fixed_effect=out_nested_fixed_effect,
			covariate_effect=covariate_effect,
			in_nested_covariate_effect=in_nested_covariate_effect,
			out_nested_covariate_effect=out_nested_covariate_effect,
			random_effect=random_effect,
			in_nested_random_effect=in_nested_random_effect,
			out_nested_random_effect=out_nested_random_effect,			
			poly_effect=poly_effect,
			poly_expression=poly_expression
	   )
	   )
	
	}
#example m1=lmt_model(t1~mu+sex+herd,covariate=~age+days(sex)+days(herd),random=~days(id)+days(dam)+pe,polyno=days~{log(sqrt(x))}+{exp(x)}) 
#example m2=lmt_model(t1~mu+sex+herd,covariate=~age+days(sex)+weight(herd),random=~days(id)+days(dam)+pe,polyno=days~{log(sqrt(x))}+{exp(x)}) 
#example m3=lmt_model(t1~mu+sex+herd,covariate=~days+age(sex)+weight(herd),random=~days(id)+days(dam)+pe,polyno=days~{log(sqrt(x))}+{exp(x)}) 

#function for nested situation
get_nested<-function(effect){

effect_expression=attr(terms(effect), which = "term.labels") # get the name of effect (include expression format)	

is_response=(attr(terms(effect),which = "response")==1) #if response variable is exist

if(is_response){ # there is reponse variable in the left of ~
effect_name=all.vars(effect)[-1] # get the name of effect(not include outside of  bracket)
}else{
effect_name=all.vars(effect) # get the name of effect(not include outside of  bracket)
}
	
in_nested_effect=NULL    #inside the bracket
out_nested_effect=NULL   #outside the bracket

if(!identical(effect_name,effect_expression)){ #  nested effect exists

nested_pos=which(effect_name!=effect_expression)

if(is_response){
tmp=attr(terms(effect), which = "variables")[-c(1:2)]  # get all variable outside and inside the bracket
}else{
tmp=attr(terms(effect), which = "variables")[-1] # get all variable outside and inside the bracket
}

	in_nested_effect=effect_name[nested_pos]
	
	for(pos in nested_pos){
	
	out_nested_effect=c(out_nested_effect,as.character(tmp[[pos]])[1]) #correspond to the nested_random_effect
	
	}

}

return(list(effect_name=effect_name,
		    in_nested_effect=in_nested_effect,
			out_nested_effect=out_nested_effect))
}	

#get the xml format of lmt_model
xml_lmt_models<-function(models #list of multiple lmt_models 
						){
id_name=models$id_name
dam_name=models$dam_name
pe_name=models$pe_name

is_gg=ifelse(identical(models$pars$gg_type,"extra"),TRUE,FALSE)

final_formula=NULL	
g_value=0 #the order of genetic effect 
gg_value=0 #the order of genetic group
ng_value=0 #the order of non-genetic effect
p_value=0 #the order of permanent effect
ma_value=0 # the order of maternal effect, for multple traits model, maternal effects only have one component
m_random_type=NULL	# the type of random effect	
formulas=models$formulas	
t_random=NULL #trait speficied random with number, for speficing var-name in AI-inv matrix 
trait_name=NULL		
for(i in 1:length(formulas)){ # c(NULL,models) make models becomes list 

m=get_models(formulas[[i]])
trait=m$trait
trait_name=c(trait_name,trait)
#for fixed effect, it doesn't have nested  and polynominal structure
m_fixed=NULL
if(length(m$fixed_effect)!=0){m_fixed=c(m_fixed,paste0(m$fixed_effect,"*T",i,m$fixed_effect))}

#for covariate effect
m_covariate=NULL
if(length(m$covariate_effect)!=0){   #for non-nested situation

	if(length(setdiff(m$covariate_effect,m$poly_effect))!=0){ #poly_expression doesn't  exists
	
		m_covariate=c(m_covariate,paste0(setdiff(m$covariate_effect,m$poly_effect),"(t(co",
									      "))*T",i,setdiff(m$covariate_effect,m$poly_effect)))
	}

	if(sum(m$poly_effect%in%m$covariate_effect)!=0){  #poly_expression exists
		
		m_covariate=c(m_covariate,paste0(m$poly_effect,"(t(co(",
										 paste0("p(",paste(1:length(m$poly_expression),collapse = ","),")"),
									      ")))*T",i, m$poly_effect))		
	}
}

if(length(m$out_nested_covariate_effect)!=0){ #for nested situation

	if(length(setdiff(m$out_nested_covariate_effect,m$poly_effect))!=0){ #poly_expression doesn't  exists
	     tmp_pos=!m$out_nested_covariate_effect%in%m$poly_effect
		m_covariate=c(m_covariate,paste0(m$out_nested_covariate_effect[tmp_pos],"(t(co(n(",
										 m$in_nested_covariate_effect[tmp_pos],
									      "))))*T",i,m$in_nested_covariate_effect[tmp_pos]))
	}

	if(sum(m$poly_effect%in%m$out_nested_covariate_effect)!=0){  #poly_expression  exists
		tmp_pos=m$out_nested_covariate_effect%in%m$poly_effect	
		m_covariate=c(m_covariate,paste0(m$out_nested_covariate_effect[tmp_pos],"(t(co(",
										 paste0("p(",paste(1:length(m$poly_expression),collapse = ","),");"),
										 paste0("n(",m$in_nested_covariate_effect[tmp_pos],")"),
									      ")))*T",i,m$in_nested_covariate_effect[tmp_pos],
										  1:sum(tmp_pos)# consider multple polynominal-nested structure
										  ))		
	}


}


#for random effect
m_random=NULL
m_random_level=NULL #for sorting the postion of different random effect in the final formula
tmp_t_random=NULL  #trait speficied random with number, for speficing var-name in AI-inv matrix 
		
if(length(m$random_effect)!=0){   #for non-nested situation

	
	for(r in m$random_effect){   
	
		if(!identical(r,m$poly_effect)){               #poly_expression doesn't  exists
		                                        # animal effect doesn't have 
		
			if(gsub('[[:digit:]]+', '', r)==id_name){     #animal effect
										                    #gsub('[[:digit:]]+', '', r) can be used for removing number from r
														 #unlist(regmatches(r, gregexpr("[[:digit:]]+", r))) can be used for getting number from r
		    g_number=unlist(regmatches(r, gregexpr("[[:digit:]]+", r)))
			g_type=paste0("g",g_number)
			m_random_type=c(m_random_type,g_type)			 
			g_value=g_value+1	
			m_random=c(m_random,paste0(id_name,"*T",i,"u",g_number,"(v(",g_type,"(",sum(m_random_type%in%g_type),")))"))	 #id*u1(v(g(1)))	
			tmp_random_name=paste0(g_type,sum(m_random_type%in%g_type)) # eg. g1,g2			
			
			m_random_level=c(m_random_level,1+g_value)

					if(is_gg){   #exist genetic group effect
					gg_value=gg_value+1
					m_random_type=c(m_random_type,"gg")
					m_random=c(m_random,paste0(id_name,"(t(gg(gg_ped)))*",id_name,"gg",i,"(v(gg(",sum(m_random_type%in%"gg"),")))"))	 #id*p1(v(p(1)))	
					tmp_random_name=paste0("gg",sum(m_random_type%in%"gg")) # eg. gg1,gg2
					m_random_level=c(m_random_level,50+gg_value)		
					}
					
			}else if(gsub('[[:digit:]]+', '', r)==dam_name){   #maternal effect	

		     g_number=unlist(regmatches(r, gregexpr("[[:digit:]]+", r)))
			g_type=paste0("g",g_number)
			g_value=g_value+1

			# if(ma_value!=0&i>=2){ #for multple traits 
				# m_random_type=c(m_random_type,NULL)
			# }else {
				# m_random_type=c(m_random_type,g_type)
				# ma_value=sum(m_random_type%in%g_type)
			# }
	
				m_random_type=c(m_random_type,g_type)
				ma_value=sum(m_random_type%in%g_type)
			
			m_random=c(m_random,paste0(dam_name,"*T",i,"m",g_number,"(v(",g_type,"(",ma_value,")))"))	 #id*m1(v(g(1)))	
			tmp_random_name=paste0(g_type,ma_value) # eg. g1,g2
			m_random_level=c(m_random_level,1+g_value)

					if(is_gg){   #exist genetic group effect
					gg_value=gg_value+1
					m_random_type=c(m_random_type,"gg")
					m_random=c(m_random,paste0(dam_name,"(t(gg(gg_ped)))*",dam_name,"gg",i,"(v(gg(",sum(m_random_type%in%"gg"),")))"))	 #id*p1(v(p(1)))	
					tmp_random_name=paste0("gg",sum(m_random_type%in%"gg")) # eg. g1,g2
					m_random_level=c(m_random_level,50+gg_value)		
					}

			}else if(r=="pe"){   #permanent effect
			
			if(is.null(pe_name)&ma_value!=0&p_value==0){
			pe_name=dam_name;cat("Without user-speficied, software will set dam_name as pe_name! \n")
			}else if(is.null(pe_name)&ma_value==0&p_value==0){
			pe_name=id_name;cat("Without user-speficied, software will set id_name as pe_name! \n")
			}
			p_value=p_value+1
			
			m_random=c(m_random,paste0(pe_name,"*T",i,"u","(v(pe(",p_value,")))"))	 #id*p1(v(p(1)))
			tmp_random_name=paste0("pe",p_value) # eg. pe1,pe2
			m_random_level=c(m_random_level,100+p_value)			
			m_random_type=c(m_random_type,"pe")
			
			}else{                #non-(genetic effect,maternal effect, permanent effect) 
							    #for this situation, the symbol of random effect is themself
			ng_value=ng_value+1
			m_random_type=c(m_random_type,r)			
			m_random=c(m_random,paste0(r,"*T",i,r,"(v(",r,"(",sum(m_random_type%in%r),")))"))	 #id*ng1(v(q(1)))
			tmp_random_name=paste0(r,sum(m_random_type%in%r)) # eg. pe1,pe2			
			m_random_level=c(m_random_level,1000+ng_value)

			}		
	
	     }else{                                 #poly_expression  exists

		       #it seems that there are no poly_expression for non-nested random effect
			
			}
		tmp_t_random=c(tmp_t_random,tmp_random_name)
	}
}


if(length(m$out_nested_random_effect)!=0){ #for nested situation


	for(r in m$in_nested_random_effect){   
		
		tmp_pos=match(r,m$in_nested_random_effect)
		
		if(!identical(m$out_nested_random_effect[tmp_pos],m$poly_effect)){               #poly_expression doesn't  exists
		
			if(gsub('[[:digit:]]+', '', r)==id_name){     #animal effect
			
		     g_number=unlist(regmatches(r, gregexpr("[[:digit:]]+", r)))
			g_type=paste0("g",g_number)			
			m_random_type=c(m_random_type,g_type)
			g_value=g_value+1	
			m_random=c(m_random,paste0(m$out_nested_random_effect[tmp_pos],"(t(co(n(",id_name,
			                                "))))","*T",i,"u",g_number,"(v(",g_type,"(",sum(m_random_type%in%g_type),")))"))	 #id*u1(v(g(1)))	
			tmp_random_name=paste0(g_type,sum(m_random_type%in%g_type)) # eg. g									
			m_random_level=c(m_random_level,1+g_value)			

					if(is_gg){   #exist genetic group effect
					gg_value=gg_value+1
					m_random_type=c(m_random_type,"gg")
					m_random=c(m_random,paste0(id_name,"(t(gg(gg_ped)))*",id_name,"gg",i,"(v(gg(",sum(m_random_type%in%"gg"),")))"))	 #id*p1(v(p(1)))	
					tmp_random_name=paste0("gg",sum(m_random_type%in%"gg")) # eg. g	
					m_random_level=c(m_random_level,50+gg_value)		
					}

			}else if(gsub('[[:digit:]]+', '', r)==dam_name){   #maternal effect
			
		    g_number=unlist(regmatches(r, gregexpr("[[:digit:]]+", r)))
			g_type=paste0("g",g_number)	

			g_value=g_value+1
			

				m_random_type=c(m_random_type,g_type)
				ma_value=sum(m_random_type%in%g_type)
			
			m_random=c(m_random,paste0(m$out_nested_random_effect[tmp_pos],"(t(co(n(",dam_name,
			                                "))))","*T",i,"m",g_number,"(v(",g_type,"(",ma_value,")))"))	 #id*m1(v(g(1)))	
			tmp_random_name=paste0(g_type,ma_value) # eg. g1									
			m_random_level=c(m_random_level,1+g_value)

					if(is_gg){   #exist genetic group effect
					gg_value=gg_value+1
					m_random_type=c(m_random_type,"gg")
					m_random=c(m_random,paste0(dam_name,"(t(gg(gg_ped)))*",dam_name,"gg",i,"(v(gg(",sum(m_random_type%in%"gg"),")))"))	 #id*p1(v(p(1)))
					tmp_random_name=paste0("gg",sum(m_random_type%in%"gg")) # eg. gg1		
					m_random_level=c(m_random_level,50+gg_value)		
					}

			}else if(r=="pe"){   #permanent effect
			
			if(is.null(pe_name)&ma_value!=0&p_value==0){
			pe_name=dam_name;cat("Without user-speficied, software will set dam_name as pe_name! \n")
			}else if(is.null(pe_name)&ma_value==0&p_value==0){
			pe_name=id_name;cat("Without user-speficied, software will set id_name as pe_name! \n")
			}
			p_value=p_value+1

			m_random=c(m_random,paste0(m$out_nested_random_effect[tmp_pos],"(t(co(n(",pe_name,
			                                "))))","*T",i,"pe","(v(pe(",p_value,")))"))	 #id*m1(v(g(1)))	
			tmp_random_name=paste0("pe",p_value) # eg. g								
			m_random_level=c(m_random_level,100+p_value)						
			m_random_type=c(m_random_type,"pe")
			}else{
			ng_value=ng_value+1
			m_random_type=c(m_random_type,r)			
			m_random=c(m_random,paste0(m$out_nested_random_effect[tmp_pos],"(t(co(n(",r,
			                                "))))","*T",i,r,"(v(",r,"(",sum(m_random_type%in%r),")))"))	 #id*ng1(v(q(1)))
			tmp_random_name=paste0(r,sum(m_random_type%in%r)) # eg. little1,little2									
			m_random_level=c(m_random_level,1000+ng_value)				
			}		
	
	     }else{                                 #poly_expression  exists

			
			if(gsub('[[:digit:]]+', '', r)==id_name){     #animal effect
		     g_number=unlist(regmatches(r, gregexpr("[[:digit:]]+", r)))
			g_type=paste0("g",g_number)					
			m_random_type=c(m_random_type,rep(g_type,length(m$poly_expression)))
			g_value=g_value+1	

			m_random=c(m_random,paste0(m$out_nested_random_effect[tmp_pos],"(t(co(",
									   paste0("p(",paste(1:length(m$poly_expression),collapse = ","),");n("),id_name,
			                                "))))","*T",i,"u",g_number,"(v(",g_type,"(",
									   paste((sum(m_random_type%in%g_type)-length(m$poly_expression)+1):(sum(m_random_type%in%g_type)),collapse=","),")))"))	 
			tmp_random_name=paste0(g_type,(sum(m_random_type%in%g_type)-length(m$poly_expression)+1):(sum(m_random_type%in%g_type))) # eg. g1,g2,g3...

			g_value=g_value+length(m$poly_expression)-1		
			m_random_level=c(m_random_level,1+g_value)				
			
					if(is_gg){   #exist genetic group effect
					gg_value=gg_value+1
					m_random_type=c(m_random_type,"gg")
					m_random=c(m_random,paste0(id_name,"(t(gg(gg_ped)))*",id_name,"gg",i,"(v(gg(",sum(m_random_type%in%"gg"),")))"))	 #id*p1(v(p(1)))	
					tmp_random_name=paste0(gg,sum(m_random_type%in%"gg")) # eg. gg1...
					m_random_level=c(m_random_level,50+gg_value)		
					}
					
			}else if(gsub('[[:digit:]]+', '', r)==dam_name){   #maternal effect
		     g_number=unlist(regmatches(r, gregexpr("[[:digit:]]+", r)))
			g_type=paste0("g",g_number)					
			g_value=g_value+1


			#m_random_type=c(m_random_type,g_type)
			m_random_type=c(m_random_type,rep(g_type,length(m$poly_expression)))
				ma_value=sum(m_random_type%in%g_type)
			
			
			
			m_random=c(m_random,paste0(m$out_nested_random_effect[tmp_pos],"(t(co(",
									   paste0("p(",paste(1:length(m$poly_expression),collapse = ","),");n("),dam_name,
			                                "))))","*T",i,"m",g_number,"(v(",g_type,"(",
									     paste((ma_value-length(m$poly_expression)+1):(ma_value),collapse=","),")))"))					 
			tmp_random_name=paste0(g_type,(ma_value-length(m$poly_expression)+1):(ma_value)) # eg. g1,g2,g3...							 									 
			g_value=g_value+length(m$poly_expression)-1
			m_random_level=c(m_random_level,1+g_value)			

					if(is_gg){   #exist genetic group effect
					gg_value=gg_value+1
					m_random_type=c(m_random_type,"gg")
					m_random=c(m_random,paste0(dam_name,"(t(gg(gg_ped)))*",dam_name,"gg",i,"(v(gg(",sum(m_random_type%in%"gg"),")))"))	 #id*p1(v(p(1)))	
					m_random_level=c(m_random_level,50+gg_value)		
					}
					
			}else if(r=="pe"){   #permanent effect

			if(is.null(pe_name)&ma_value!=0&p_value==0){
			pe_name=dam_name;cat("Without user-speficied, software will set dam_name as pe_name! \n")
			}else if(is.null(pe_name)&ma_value==0&p_value==0){
			pe_name=id_name;cat("Without user-speficied, software will set id_name as pe_name! \n")
			}
			p_value=p_value+1

			m_random=c(m_random,paste0(m$out_nested_random_effect[tmp_pos],"(t(co(",
									   paste0("p(",paste(1:length(m$poly_expression),collapse = ","),");n("),pe_name,
			                                "))))","*T",i,"pe","(v(pe(",
										 paste(p_value:(p_value+length(m$poly_expression)-1),collapse=","),")))"))
			tmp_random_name=paste0("pe",p_value:(p_value+length(m$poly_expression)-1)) # eg. pe1,pe2,pe3...								
			p_value=p_value+length(m$poly_expression)-1			
			m_random_level=c(m_random_level,100+p_value)
			m_random_type=c(m_random_type,rep("pe",length(m$poly_expression)))
			#m_random_type=c(m_random_type,"pe")
			}else{       #non-genetic effect
			ng_value=ng_value+1
			m_random_type=c(m_random_type,r)			
			m_random=c(m_random,paste0(m$out_nested_random_effect[tmp_pos],"(t(co(",
									   paste0("p(",paste(1:length(m$poly_expression),collapse = ","),");n("),r,
			                                "))))","*T",i,r,"(v(",r,"(", 
										  paste(sum(m_random_type%in%r):(sum(m_random_type%in%r)+length(m$poly_expression)-1),collapse=","),")))"))
			tmp_random_name=paste0("r",sum(m_random_type%in%r):(sum(m_random_type%in%r)+length(m$poly_expression)-1)) # eg. little1...	
										  
			ng_value=ng_value+length(m$poly_expression)-1		
			m_random_level=c(m_random_level,1000+ng_value)				

			}		
			
		}
	tmp_t_random=c(tmp_t_random,tmp_random_name)
	}
}


m_random=m_random[order(m_random_level)] #make sure the formula easliy readable in the final formula
m_formula=paste0(trait,"=",paste(c(m_fixed,m_covariate,m_random),collapse="+"))
final_formula=c(final_formula,m_formula)
t_random=c(t_random,list(c(tmp_t_random,paste0("r",i)))) #trait_specified_random effect name with number
#print(m)
# print(m$trait)
# print(m_random)
# print(m_formula)
}
names(t_random)=trait_name
final_random_type=m_random_type # for assign intial variance components
final_random_type=c(final_random_type,rep("r",length(formulas))) #including residual variance

m_random_type=unique(m_random_type)
m_random_type=m_random_type[order(match(m_random_type,c("g","pe",setdiff(m_random_type,c("g","pe")))))]

#model xml 
model_xml="  <models>" 
if(length(final_formula)!=0){
	
	model_xml=rbind(model_xml,"    <eqn attributes=\"strings\">")
	
	for(formula in final_formula){
	
		model_xml=rbind(model_xml,paste0("      ",formula))	
	
	}
	model_xml=rbind(model_xml,"    </eqn>")	
}

if(length(m$poly_expression)!=0){
	
	model_xml=rbind(model_xml,"    <poly attributes=\"strings\">")
	
	for(poly_exp in m$poly_expression){
	
		model_xml=rbind(model_xml,paste0("      ",poly_exp))	
	
	}
	model_xml=rbind(model_xml,"    </poly>")	

}
	model_xml=rbind(model_xml,"  </models>" )	
	
return(list(xml=model_xml,type=final_random_type,t_random=t_random))
}


 #add EBV
 getEBV = function(trait_name) {
		sol=read.csv("results.csv",stringsAsFactors=F,header=F)
		colnames(sol)=c("effect","subname","level","value")
		sol$trait=substr(sol$subname,1,2)
		sol$sub_effect=substr(sol$subname,3,50)

		EBV=sol[sol$effect%in%"g",]

		tEBV=NULL
		for(i in 1:length(trait_name)){

			ebv=sol[sol$effect%in%"g"&sol$trait%in%paste0("T",i),c("sub_effect","level","value")]
			colnames(ebv)=c("genetic effect","id","value")
			i_sol=sol[sol$trait%in%paste0("T",i),]
			write.csv(ebv,paste0("Rlmt.",trait_name[i],".ebv.csv"),row.names=F)
			#write.csv(i_sol,paste0(trait_name[i],".Rlmt.solution.csv"),row.names=F) #maybe no need to write 
			tEBV=c(tEBV,list(list(ebv=ebv,sol=i_sol)))
		}
		names(tEBV)=trait_name
		return(tEBV)
 }	