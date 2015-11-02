library(MASS)
###########################################
############### Train Model ###############
###########################################
# Find rank one matrix. Called by eigen.model
rankOne.calc <- function( Qmat, blocks ) {
   M <- length(blocks)
   N <- 0
   for ( i in 1:(M-1) ) {
      for ( j in (i+1):M ) {
         N <- N + length(blocks[[i]])*length(blocks[[j]])
      }
   }
   L <- dim(Qmat)[1]

   # Set up system of equations
   X <- matrix( rep(0,L*N) ,ncol=L )
   q.vec <- rep(0,N)
   ind <- 1
   for ( h in 1:(M-1) ) {
      for ( i in (h+1):M ) {
         for ( j in blocks[[h]] ) {
            for ( k in blocks[[i]] ) {
               q.vec[ind] <- log(abs(Qmat[j,k]))
               X[ind,j] <- 1
               X[ind,k] <- 1
               ind <- ind + 1
            }
         }
      }
   }

   # Least squares solution
   T <- c()
   if ( M == 2 ) {
      X.inv <- ginv(X[,2:L])
      T <- X.inv %*% q.vec
      T <- c(0,T[,1])
   } else {
      X.inv <- ginv(X)
      T <- X.inv %*% q.vec
   }

   # Calculate Rmat
   Rmat <- exp(T) %*% t(exp(T))
   Rmat
}

######################### Read in training data and fit model ###########################
# annotations: File with annotations
# blocks: List of arrays giving block structure ( e.g. blocks[[1]]=c(1,2,3,4), blocks[[2]]=c(5,6,7,8,9,10,11,12) )
# columns: Array listing which columns to read in from annotations ( e.g. 5:20 )
eigen.model <- function( annotations, blocks, columns, type=NULL, save.file=NULL, invert=NULL, header=TRUE, sep="\t" ) {
	    # Find number of columns in annotations
	    anno1 <- read.table(annotations,nrows=1,sep=sep)
	    dim.anno <- dim(anno1)

	    # Read in full dataset
	    cols <- rep("NULL",dim.anno[2])
	    cols[columns] <- "character"
	    anno <- read.table(annotations,colClasses=cols,header=header,sep=sep)
	    d.anno <- dim(anno)
	    print(d.anno)

	    for ( i in 1:d.anno[2] ) { 
	    	anno[,i] <- as.numeric(anno[,i])
	    }

	    keep <- 1:d.anno[1]
	    if ( type == "coding" ) {
	      keep <- c()
	      for ( i in 1:4 ) {
	       	  keep <- union(keep,which(!is.na(anno[,i])))
	      }
	    } else if ( type == "noncoding" ) {
	      for ( i in 1:4 ) {
	       	  keep <- intersect(keep,which(is.na(anno[,i])))
	      }
	    }

	    anno <- anno[keep,]

	    if ( !is.null(invert) ) {
	       for ( i in invert ) {
	       	   anno[,i] <- 1 - anno[,i]
	       }
	    }

	    means <- rep(0,d.anno[2])
	    sds <- rep(0,d.anno[2])
	    for ( i in 1:(d.anno[2]) ) {
	    	means[i] <- mean(anno[,i],na.rm=TRUE)
		sds[i] <- sd(anno[,i],na.rm=TRUE)
		anno[,i] <- (anno[,i]-means[i])/sds[i]
	    }
	    print(dim(anno))

	    include <- c()
	    blocks2 <- list()
	    start <- 1
	    for ( i in 1:length(blocks) ) {
	    	include <- c(include,blocks[[i]])
		blocks2[[i]] <- seq(start,(start+length(blocks[[i]])-1),1)
		start <- start + length(blocks[[i]])
	    }
	    print(include)
	    means <- means[include]
	    sds <- sds[include]
	    anno <- anno[,include]

	    # Calculate covariance
	    cov1 <- cor(anno,use="pairwise.complete.obs")
	    print("cov1")
	    print(dim(cov1))
	    
	    # Calculate R matrix and weights
	    R.mat <- rankOne.calc( cov1, blocks2 )
	    print("R.mat")
	    print(dim(R.mat))
	    eigen.decomp <- eigen(R.mat)
	    weights <- eigen.decomp$vectors[,1]

	    # Record and save/return results
	    model <- list()
	    model$cov.matrix <- cov1
	    model$rank.one.matrix <- R.mat
	    model$weights <- weights
	    model$means <- means
	    model$sds <- sds

	    if ( !is.null(save.file) ) {
	       save(model,file=save.file)
	    }
	    model
}

######################### Eigen-PC ###########################
eigen.model.pc <-
function( annotations, columns, save.file=NULL, invert=NULL, header=TRUE, sep="\t" ) {
    # Find number of columns in annotations
    anno1 <- read.table(annotations,nrows=1,sep=sep)
    dim.anno <- dim(anno1)
    
    # Read in full dataset
    cols <- rep("NULL",dim.anno[2])
    cols[columns] <- "character"
    anno <- read.table(annotations,colClasses=cols,header=header,sep=sep)
    d.anno <- dim(anno)
    
    for ( i in 1:d.anno[2] ) {
        anno[,i] <- as.numeric(anno[,i])
    }
    
    if ( !is.null(invert) ) {
        for ( i in invert ) {
            anno[,i] <- 1 - anno[,i]
        }
    }
    
    means <- rep(0,d.anno[2])
    sds <- rep(0,d.anno[2])
    for ( i in 1:(d.anno[2]) ) {
        means[i] <- mean(anno[,i],na.rm=TRUE)
		 sds[i] <- sd(anno[,i],na.rm=TRUE)
		 	anno[,i] <- (anno[,i]-means[i])/sds[i]
    }
    
    # Calculate covariance
    cov1 <- cor(anno,use="pairwise.complete.obs")
    
    # Calculate R matrix and weights
    eigen.decomp <- eigen(cov1)
    weights <- eigen.decomp$vectors[,1]
    
    # Record and save/return results
    model <- list()
    model$cov.matrix <- cov1
    model$weights <- weights
    model$means <- means
    model$sds <- sds
    
    if ( !is.null(save.file) ) {
        save(model,file=save.file)
    }
    model
}


###########################################
############# Calculate Scores ############
###########################################
### Called by eigen.score
calc.scores <- function( x, columns, means, sds, weights, consequence, invert, impute.val, impute.which, impute.cons, outfile ) {
	    # Adjust annotations (rescale, normalize, impute)
	    adj <- as.numeric(as.character(x[columns])) # Make datatype numeric
	    if ( !is.null(invert) ) {
	    	    adj[invert] <- 1 - adj[invert]  # switch scale (e.g. SIFT, allele frequency)
	    }
	    missing <- which(is.na(adj))  # Find missing values
	    adj[missing] <- means[missing] # Impute missing values
	    adj <- (adj-means)/sds  # Standardize annotations
	    
	    # Impute special values for select consequences
	    if ( !is.null(impute.cons) ) {
	       impute <- 0
	       cons <- as.character(x[consequence])
	       for ( c in impute.cons ) {
	       	   impute <- impute + grepl(c,cons)
	       }
	       if ( impute > 0 ) {
	       	  adj[impute.which] <- (impute.val-means[impute.which])/sds[impute.which]
	       }
	    }

	    score.N <- t(adj) %*% weights
	    out <- paste(c(as.character(x),score.N),collapse="\t")
	    cat(out,file=outfile,sep="\n")
}
	    

############# Read in model file and setup score calculation ############
# annotations: File with annotations
# model: Comma separated list of R data files containing model information
# use: Indicate which annotations to use in the score (e.g. 1:12)
# columns: Columns containing annotations specified by use (e.g. 5:16)
# header: Does the annotation file have a header line? (TRUE/FALSE)
# consequence: Column containing functional consequence (e.g. 21 )
# invert: Annotations which should be rescaled as 1-x (e.g. c(1,17,18,19,20) )
# impute.val: Vector of values to use for imputation of special values (e.g. c(1,1,1,5.37) )
# impute.which: Vector of annotations to apply imputation to (e.g. c(1,2,3,4) )
# impute.cons: Use non-mean imputation for consequences in impute.cons (e.g. c("stop_lost","splice_site") )
eigen.score <- function( annotations, model1, model2, use, columns, save.file, header=NULL, consequence=NULL, invert=NULL, impute.val=NULL, impute.which=NULL, impute.cons=NULL ) {
	    anno.input <- gzfile(annotations,open="r")	   # Open connection to read in annotations

	    # Process model files
	    means.list <- c()
	    sds.list <- c()
	    weights.list <- c()
	    model.files <- strsplit(model1,",")[[1]]
	    #print(use)
	    for ( i in 1:length(model.files) ) {
	    	try.result <- try(load(model.files[i]),silent=TRUE)
		if ( class(try.result) != "try-error" ) {
		   weights.list <- rbind(weights.list,abs(model$weights[use]/sum(model$weights[use])))
		} else { 
		  print("not found")
		}
	    }

	    model.files <- strsplit(model2,",")[[1]]
	    #print(use)
	    for ( i in 1:length(model.files) ) {
	    	try.result <- try(load(model.files[i]),silent=TRUE)
		if ( class(try.result) != "try-error" ) {
		   means.list <- rbind(means.list,model$means[use])
		   sds.list <- rbind(sds.list,model$sds[use])
		} else { 
		  print("not found")
		}
	    }


	    means <- rep(0,length(use))
	    sds <- rep(0,length(use))
	    weights <- rep(0,length(use))
	    for ( i in 1:length(use) ) {
	    	means[i] <- median(means.list[,i],na.rm=TRUE)
		sds[i] <- median(sds.list[,i],na.rm=TRUE)
		weights[i] <- median(weights.list[,i],na.rm=TRUE)
	    }

	    output <- gzfile(save.file,open="w")
	    if ( !is.null(header) ) {
	       out <- paste(header,collapse="\t")
	       cat(out,file=output,sep="\n")
	    }
	    calc_score <- function(x) { calc.scores(x,columns,means,sds,weights,consequence,invert,impute.val,impute.which,impute.cons,output) }
	    while ( length( line <- readLines(anno.input, n = 1, warn = FALSE ) ) > 0 ) {
	    	  parts <- strsplit(as.character(line),"\t")
		  lapply(parts,calc_score)
	    }
	    close(anno.input)
	    close(output)
}




#####################################