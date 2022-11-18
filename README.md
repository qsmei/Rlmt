
### Rlmt: R wrapper of the Linear Mixed Models Toolbox (LMT)

## Contents

-   [OVERVIEW](#overview)

-   [GETTING STARTED](#getting-started)

-   [USAGE](#usage)

------------------------------------------------------------------------


### OVERVIEW 

​	The **L**inear mixed **M**odels **T**oolbox (**lmt**) is a stand-alone single executable software for for large scale linear mixed model analysis. **lmt** has been used successfully for genetic evaluation data sets with >>200k genotyped animals, >>15m animals, >>500m equations. More details about lmt can be found in website: [lmt website](https://dmu.ghpc.au.dk/lmt/wiki/index.php?title=The_Linear_Mixed_Models_Toolbox)

​	In order to further assistance the use of lmt, I developed an R package: Rlmt for interfacing with lmt. This R package is developed via R6 package which can provides an implementation of encapsulated object-oriented programming for R.  And this makes the grammar of Rlmt looks like easy and elegant (i hope :smile:). 

​	About the usage of lmt in academic and commercial,  user should follow the rules  as mentioned in lmt website. 



## GETTING STARTED

### 🙊Installation

`Rlmt` links to R packages `R6`. This package should be installed before installing `Rlmt`.   Then, user can install `Rlmt` from github directly. 

```R
install.packages("R6")
devtools::install_github("qsmei/Rlmt")
```

If you have any questions about Rlmt, please contact the *author* (quanshun1994@gmail.com) or (olef.christensen@qgg.au.dk).

### 🙊Examples

-   Ex 1. PBLUP
-   Ex 2. SSGBLUP
-   Ex 3. Multiple traits analysis
-   Ex 4. Several special models 
-   Ex 5. Modify jobs

## Usage

Due to the requirement of LMT, user need to set up the following parameters in the Linux terminal:

```shell
ulimit -s unlimited
export OMP_NUM_THREADS=32
export OMP_STACKSIZE=2000M
export OMP_DYNAMIC=false
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_MAX_ACTIVE_LEVELS=2147483647
```

`Rlmt` provides several files which are saved in `~/Rlmt/extdata`. We can get the path of these files by typing

``` {.r}
system.file("extdata", package = "Rlmt") # path of example files
```

#### Ex 1. PBLUP

``` R
library(Rlmt)
# path of example data  
example_path=system.file("extdata", package = "Rlmt") 
#construct lmt_data object
mydata=lmt_data$new(phe_file=paste0(example_path,"/data.csv"),					
                    ped_file=paste0(example_path,"/pedigree_unsorted.csv")) 

#construct  lmt_models object
mymodels=lmt_models$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id)

#construct Rlmt object
mylmt=lmt$new(models=mymodels,data=mydata,software_path="/usr/home/qgg/vinzent/lmt")

#default of Rlmt: only solved mixed model with provided variance components
mylmt$run_lmt("/usr/home/qgg/qumei/lmt/test_Result") #output path

#modify variacne components directly, Rlmt also allows user to provide the file of variance components
mymodels$pars$add_vars(value=49,name="g")
mymodels$pars$add_vars(value=15,name="r")

#estimate variance components by airemlc
m_jobs=mymodels$pars$jobs
m_jobs$type="airemlc"

mylmt$run_lmt("/usr/home/qgg/qumei/lmt/test_Result1") #output path
```

#### Ex 2. SSGBLUP

``` R
library(Rlmt)
# path of example data  
example_path=system.file("extdata", package = "Rlmt") 
#construct lmt_data object
mydata=lmt_data$new(phe_file=paste0(example_path,"/data.csv"),				    
                    ped_file=paste0(example_path,"/pedigree_unsorted.csv"),					
                    geno_file=paste0(example_path,"/genotypes.txt"),					
                    geno_id_file=paste0(example_path,"/genoids.csv"))

#construct  lmt_models object
mymodels=lmt_models$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id)
mymodels$pars$blup_type="SS_GBLUP"

#modify variance components directly, Rlmt also allows user to provide the file of variance components
mymodels$pars$add_vars(value=49,name="g")
mymodels$pars$add_vars(value=15,name="r")

#estimate variance components by airemlc
m_jobs=mymodels$pars$jobs
m_jobs$type="airemlc"

#construct Rlmt object
mylmt=lmt$new(models=mymodels,data=mydata,software_path="/usr/home/qgg/vinzent/lmt")

#do analysis
mylmt$run_lmt("/usr/home/qgg/qumei/lmt/test_Result2") #output path


#################Metafounder###################
mydata=lmt_data$new(phe_file=paste0(example_path,"/data.csv"),				    
                    ped_file=paste0(example_path,"/pedigreeGG_unsorted.csv"),					
                    geno_file=paste0(example_path,"/genotypes.txt"),					
                    geno_id_file=paste0(example_path,"/genoids.csv"))
mymodels=lmt_models$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id)
mylmt=lmt$new(models=mymodels,data=mydata,software_path="/usr/home/qgg/vinzent/lmt")
mymodels$pars$blup_type="SS_GBLUP"
mylmt$models$pars$meta_gamma=0.5
mylmt$run_lmt("/usr/home/qgg/qumei/lmt/test_Result") #output path


#################ssGBLUP with two genomic factors###################
mydata=lmt_data$new(phe_file=paste0(example_path,"/data.csv"),
					ped_file=paste0(example_path,"/pedigree_unsorted.csv"),
					geno_file=c(paste0(example_path,"/genotypes1.txt"),paste0(example_path,"/genotypes2.txt")),
					geno_id_file=c(paste0(example_path,"/genoids1.csv"),paste0(example_path,"/genoids2.csv")))


mymodels=lmt_models$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id1+id2)
mymodels$pars$blup_type="SS_GBLUP"
mylmt=lmt$new(models=mymodels,data=mydata,software_path="/usr/home/qgg/vinzent/lmt")
mylmt$models$pars$jobs$type="airemlc"
mylmt$run_lmt("/usr/home/qgg/qumei/lmt/test_Result") #output path

##################SS_TBLUP###################
mydata=lmt_data$new(phe_file=paste0(example_path,"/data.csv"),				    
                    ped_file=paste0(example_path,"/pedigree_unsorted.csv"),					
                    geno_file=paste0(example_path,"/genotypes.txt"),					
                    geno_id_file=paste0(example_path,"/genoids.csv"))
mymodels=lmt_models$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id)
mymodels$pars$blup_type="SS_TBLUP"
mylmt=lmt$new(models=mymodels,data=mydata,software_path="/usr/home/qgg/vinzent/lmt")
mylmt$run_lmt("/usr/home/qgg/qumei/lmt/test_Result") #output path

##################SS_SNPBLUP################### 
#need to fix some errors in the next step

#################ssGBLUP model with GRM supplied externally###################
#need to fix some errors in the next step


```

#### Ex 3. Multiple traits analysis

``` R
library(Rlmt)
# path of example data  
example_path=system.file("extdata", package = "Rlmt") 
#construct lmt_data object
mydata=lmt_data$new(phe_file=paste0(example_path,"/data.csv"),					
                    ped_file=paste0(example_path,"/pedigree_unsorted.csv")) 
m1=lmt_models$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id) #model of trait1
mylmt=lmt$new(models=m1,data=mydata,software_path="/usr/home/qgg/vinzent/lmt")

#construct another lmt_models object
m2=lmt_formulas$new(fixed=tr2~f11+f12+f13,covariate=~c1c1,random=~id) #model of trait2

mylmt$models$add_formulas(m2) #add new trait in analysis
mylmt$models$pars$jobs$type="airemlc" #estimate variance components by airemlc
mylmt$run_lmt("/usr/home/qgg/qumei/lmt/test_Result2") #multple traits analysis
```

#### Ex 4.  Several special models 

``` R
library(Rlmt)
example_path=system.file("extdata", package = "Rlmt") 
mydata=lmt_data$new(phe_file=paste0(example_path,"/data.csv"),					
                    ped_file=paste0(example_path,"/pedigree_unsorted.csv")) 
# maternal effect
m1=lmt_models$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id+dam) 

# permanent effect
m1=lmt_models$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id+dam+pe) 

#nested effect
m1=lmt_models$new(fixed=tr1~f11+f12,covariate=~f13(c1c1),random=~id+dam+pe)

#user defined polynomial expression 
m1=lmt_models$new(fixed=tr1~f11+f12,covariate=~c1c1,random=~id+dam+pe,,polyno=c1c1~{x^2}+{exp(x)}) 

#3-order Legendre polynomials
m1=lmt_models$new(fixed=tr1~f11+f12,covariate=~c1c1,random=~id+pe,,polyno=c1c1~{l1}+{l2}+{l3}) 

mylmt=lmt$new(models=m1,data=mydata,software_path="/usr/home/qgg/vinzent/lmt")
mylmt$models$pars$jobs$type="airemlc"
mylmt$run_lmt("/usr/home/qgg/qumei/lmt/test_Result3") #do analysis
```

#### Ex  5.  Modify jobs

``` R
#sampler 
m_jobs=lmt_jobs$new(type="sample")
mymodels$pars$jobs=m_jobs
m_jobs$samplers$type=c("singlepass","blocked")
mylmt$run_lmt("....../test_Result") #output path

		
#solver 
m_jobs=mymodels$pars$jobs
m_jobs$type="mcemreml"
m_jobs$conv=-9.21034
m_jobs$rounds=300
mylmt$run_lmt("....../test_Result") #output path

#pevsolve
m_jobs=lmt_jobs$new(type="pevsolve")
mymodels$pars$jobs=m_jobs
mylmt$run_lmt("....../test_Result") #output path
```
