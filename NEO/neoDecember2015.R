# neo.txt: Network Edge Orienting. Learn a Bayesian network
#          for gene expression traits and clinical phenotype trait
#          data, using SNPs as causal anchors and compute
#          structural confidence scores for each directed
#          edge between traits.
# 
# Uses SNPs or genetic variation to orient edges in
# an (mRNA trait,clinical trait) correlation network. 
#
# Methods implemented (but not necessary active by default)
# include ZEO (Z-score Edge Orienting), LEO (Local-structure 
# Edge Orienting), RVV (Recursive V-structures with Verification),
# and refinements of Cohen et al's FTC filtering approach for
# learning SEMs.
# 
# As currently configured, NEO is set up to automatically
# identify SNPs (columns starting with 'SNP', snd a single
# Trait column in the data x from the column names, and 
# process these data accordingly; searching for a causal
# directed but not-necessarily-acyclic graph within the 
# SNP-gene-trait variation data.


# Scott Ritchie (scottr@student.unimelb.edu.au)
# made the code more robust June 2014.
# Peter Langfelder updated the code in the past.

# Update December2015: comment out code around require(rcom). 
# The package rcom is not available and the error the code has been
# working around may not exist anymore.
#
# Copyright (c) 2007, 2008, Jason E. Aten  <j.e.aten@gmail.com>
# Intial implementation: 14 January 2007
#
# Licensing: 
#   Licensed non-exclusively under the GNU GPL version >= 3.
#   http://www.gnu.org/licenses/gpl-3.0.txt
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, 
# INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND 
# FITNESS FOR A PARTICULAR PURPOSE.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
# THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
# OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
# OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# neover(): Automatically updated CVS source code control information; 
#     the R user can call neover() in R to get the version of the
#     NEO code base they are running.
#
# CVS Version = Id: neo.txt,v 1.339 2008/02/12 08:26:41 jaten Exp 
# Became Subversion = Id: neo.txt 23 2008-02-12 08:38:20Z jaten 
#
# Currently: $Id: neo.txt 78 2008-03-14 00:48:54Z jaten $

if(exists("neover") ) rm(neover);
neover=function(quiet=FALSE) {a=paste(sep="","neo version: ","$Id: neo.txt 78 2008-03-14 00:48:54Z jaten $"); if (!quiet) { cat(a); cat("\n") }; invisible(a)}

#
#
# neo(): Network Edge Orienting, top level wrapper function
#
# Author: Jason E. Aten <j.e.aten@gmail.com>
# Jan-June 2007
#
# Arguments:
#
# datCombined = combined data frame containing SNPs and traits. SNPs 
#             are also called exogenous variables in Structural 
#             Equation Modeling. This means they have no parents in the graph.
#
#             In contrast to SNPs, traits (either clinical traits, mRNA 
#             levels, cluster centers, or Princple components from a PCA
#             or cluster analysis) are endogenous variables and so are
#             expected to have parents in the graph, namely one or more
#             of the SNPs. 
#
#             Every column of datCombined that is specified in snpcols
#             (or has a name that begins with "SNP" or "snp"), 
#             is considered an exogenous SNP column. For instance, it
#             may be useful to specify sex (male/female) as an
#             exogenous variable in snpcols to see when sex is
#             playing an important role in gene expression causality.
#
#             Every other column (not in snpcols) is considered an 
#             endogenous trait column.
#
#             Typically there are two kinds of traits used in this analysis:
#             the mRNA or gene expression traits, and the clinical traits.
#
#             The clinical traits are typically the focus of the
#             analysis, and we usually want to find if genes (mRNA traits)
#             are causal or reactive to the (clinical) traits. Typically 
#             a single trait at a time would be considered (set traitcols= to 
#             just the single column containing the clinical trait of interest),
#             but multiple trait analysis is allowed and performed if more
#             that one traitcols is specified.
#
#             The only distinction that setting traitcols offers over
#             all the given non-snp columns, is that NEO will always try 
#             to oreint edges into these special variables, even if they are 
#             in a supernode and so edge orienting scores show
#             higher variation and deserve a more cautious interpretation.
#
# snpcols  =  columns of datCombined that contain SNPs. The default is to 
#             call noe.guess.snp(datCombined) which guesses that the snpcols
#             are those that start with "SNP" (upper or lowercase)
#
# traitcols = Traits are special in that we always attempt to orient
#             edges impinging on a trait, even if the trait is in a
#             supernode with other variables it is highly correlated with.
#
# skip.LEO  = if TRUE, we don't even bother to do the LEO computations at all.
#             Default is false. The LEO methods are much more sensitive, but 
#             are slower. For a quick look at the ZEO results, set 
#             skip.LEO=TRUE.
#
# quiet     = if TRUE, no log output is displayed on the R console. A 
#             logfile is still generated, given a current timestamp and 
#             pm$run.title, and written to disk. If quiet is FALSE, then
#             everything to the screen is also sink()-ed to a log file
#             whose name ends with '.txt'
#
# pm        = a vector of parameters that control the computation. To 
#             modify these, it is probably easiest to get some good 
#             defaults by calling neo.get.param(), which is also the default 
#             action, and then modify just the members of the named list
#             that you wish to change. See the neo.get.param() function
#             just below neo() for documentation of the parameters.
#
# use.ranks = if TRUE, replace the data in the datCombined frame by their
#             ranks before doing any other processing, so as to get a 
#             non-parametric computation (e.g. Spearman correlations and
#             partials are then automatically computed rather than Pearson 
#             correlations and partials). FALSE by default.
#
# repos.nodes = if the final graph summary picture is drawn and repos.nodes
#             is TRUE, then the program will pause to allow the user to
#             click on the picture to rearrange the nodes in the drawn graph.
#             The node nearest to the click is moved to the point of the
#             click.
#
# pm$ignorable.edge.list = a list of edge pairs (the indices of columns of 
#             datCombined) that will be ignored. Useful for speeding up
#             processing when you only care about relationships between
#             classes of traits, such as between genes and clinical traits,
#             rather than relationships between the clinical traits 
#             themselves.
#
# pm$run.title = string giving name of experiment, appended to log directory
#             name (should have no spaces) to help track different experiments.
#
# pm$no.log = defaults to FALSE. If TRUE, don't produce the (somewhat 
#             extensive) logging of results. Useful for doing lots of 
#             simulation calls to neo() where you just what the programatic
#             output in the returned value.
#
if(exists("neo") ) rm(neo);
neo=function(datCombined, snpcols=neo.guess.snps(datCombined,pm), 
             traitcols=neo.guess.traits(datCombined,pm), pm=neo.get.param(),quiet=FALSE,
             repos.nodes=FALSE, skip.LEO=FALSE, use.ranks=FALSE) {

   # we crash if starting directory is too long.
   if (nchar(getwd()) > 150) {
      stop("Must run NEO from a directory whose path is short: nchar(getwd()) < 150. Otherwise the file system calls can fail when logging. Use Sysinternals' 'junction' if necessary to creation symlinks.");
   }

   pm.print(pm,paste("All starting...",date(), "with version",neover()))
   sanity.column.names.check = sapply(colnames(datCombined),starts.with.zero.thru.nine)
   if(any(sanity.column.names.check)) {
        stop(paste("column names of datCombined cannot start with 0-9 as digits. Aborting. See column", colnames(datCombined)[sanity.column.names.check]))
   }

    # screen out factors in snpcols, because we want alleles numeric and ordered for F2 crosses
    # and factors as default read in may not be ordered
    for (i in snpcols) {
       if (is.factor(datCombined[,i])) stop("The specified snpcols must be numeric rather than factors so that correlations are computed correctly.")
    }
   
   # sanity check the supplied data frame -- it should not have any NA or non-varying columns
   check.rather.than.impute=TRUE
   if (pm$rough.and.ready.NA.imputation) {
      check.rather.than.impute=FALSE;
   }

   # get names of any forced SNP assignments before anything changes as far as datCombined
   # columns goes.

   if (!is.null(pm$forced.MA.colnum)) {
      pm$forced.MA.cn = colnames(datCombined)[pm$forced.MA.colnum] # cn stands for names
   } else {
      pm$forced.MA.cn = NULL
   }

   if (!is.null(pm$forced.MB.colnum)) {
      pm$forced.MB.cn = colnames(datCombined)[pm$forced.MB.colnum]
   } else {
      pm$forced.MB.cn = NULL
   }


   datCombined=check.or.impute.missing.data(datCombined=datCombined,
                 snpcols=snpcols,pm=pm,check=check.rather.than.impute)

#   co.lin=check.multi.collinearity(datCombined,pm=pm)
#   if (co.lin$has.redundant) {
#      stop(co.lin$msg)
#   }

   if (skip.LEO) pm$use.LEO = FALSE; # force use of just ZEO if we aren't 
                                     # computing LEO

   if (use.ranks) {
      for (j in 1:ncol(datCombined)) {
         datCombined[,j]=rank(datCombined[,j])
      }
   }

   # setup for graphing at the end.
   zGRAPH=list() # need to define: zGRAPH$numsnps, zGRAPH$nc, zGRAPH$num.non.snps
   # defaults:
   zGRAPH$numsnps = length(snpcols)
   zGRAPH$num.non.snps = length(setdiff(1:ncol(datCombined),snpcols))
   zGRAPH$nc = zGRAPH$numsnps + zGRAPH$num.non.snps

   # set defaults for permutation block, if needed.
   if (is.null(pm$block.to.permute) | length(pm$block.to.permute) ==0) {
      pm$block.to.permute=snpcols
   }
   if (any(!(pm$block.to.permute %in% (1:ncol(datCombined))))) {
      stop("Bad specification of pm$block.to.permute: out of range")
   }
   # preserve names in case columns get re-arranged.
   pm$block.to.permute.cn = colnames(datCombined)[pm$block.to.permute]


   if (!is.null(pm$A) & !is.null(pm$B)) {
      pm$just.AB = TRUE

      # we have A and B and so we need only consider these columns in 
      # addition to the snpcols
      if (length(intersect(pm$A,pm$B)) > 0) { 
            stop("pm$A and pm$B sets must be disjoint. Halting.")
      }

      if (any(!(pm$A %in% 1:ncol(datCombined)))) stop ("pm$A values out of bounds.")
      if (any(!(pm$B %in% 1:ncol(datCombined)))) stop ("pm$B values out of bounds.")

      if (any((pm$A %in% snpcols))) stop ("pm$A values out of bounds: cannot be in snpcols.")
      if (any((pm$B %in% snpcols))) stop ("pm$B values out of bounds: cannot be in snpcols.")

      pm$A.cn = colnames(datCombined)[pm$A] # so we can re-match them easily 
             # later, even after shrinking the matrix after SNP selection.
      pm$B.cn = colnames(datCombined)[pm$B]

      keepers=sort(union(c(snpcols,traitcols,pm$A),pm$B))
      snpcols=match(snpcols,keepers) # if rearranged, get back proper snpcols
      traitcols = match(union(pm$A,pm$B),keepers)
      pm$A   =match(pm$A, keepers)   # ditto
      pm$B   =match(pm$B, keepers)   # ditto
      datCombined=datCombined[,keepers]

      # adjust graph drawing stuff.
      zGRAPH$numsnps = length(snpcols)
      zGRAPH$num.non.snps = length(traitcols)
      zGRAPH$nc = length(keepers)

      # enumerate the edges to be checked. Save in pm$enumerated.A.cn, pm$enumerated.B.cn
      pm=get.enumab.cn.set(pm)

   } else {
      # pm$A/B not called for, but still use the pm$enumerated.A.cn, pm$enumerated.B.cn
      # lists to organize execution into batch execution.
      
      # Can't use get.enumab.cn.set when A and B lists overlap, so do it manually:
      cn=colnames(datCombined)

      non.snp.cols=setdiff(1:ncol(datCombined),snpcols)
      lenT = length(non.snp.cols)
      n.rows = lenT * (lenT-1) /2
      pm$enumerated.A.cn = rep(NA,n.rows) # initialize
      pm$enumerated.B.cn = rep(NA,n.rows)
      k=1
      for (i in 1:(lenT-1)) {
         for (j in (i+1):lenT) {
             pm$enumerated.A.cn[k] = cn[non.snp.cols[i]]
             pm$enumerated.B.cn[k] = cn[non.snp.cols[j]]
             k=k+1
         }
      }
   }

   if (quiet) pm$quiet=quiet;

   if (pm$no.log) {
      neo.log.file="no.neo.log"
   } else {
     # setup logging - put all log files in a new, timestamped directory
     neo.log.file = neo.mkdir.logfile(pm$run.title);
     pm$neo.log.file=neo.log.file

     pm$orig.neo.log.file=neo.log.file
     sink(file = paste(sep="",neo.log.file,".txt"), append=TRUE, type="output",split = !pm$quiet);
     on.exit(while (sink.number()) sink()); # turn off log in case of crash, and upon normal exit.
     pm.print(pm,paste("Logfile:",neo.log.file));
   }

   # do we open excel at the end?
   pm$open.excel.at.end = FALSE # default
   if (!pm$no.log & !pm$quiet) pm$open.excel.at.end = TRUE
   pm$open.excel.now = FALSE # set to true after last permutation

     # setup control structure--parameters in pm, for the run
     if(is.null(pm)) {
       pm = neo.get.param(logpath=neo.log.file);
     }

    # not necessary, these are the defaults
    #   pm$use.leo=TRUE; # LEO if true, ZEO if false.
    #   pm$eo.type="forward.step.selection"; # additional options: "all.snp.permutations.equal.votes","max.vs.max"

    #necessary
    pm$no.obs.Z = nrow(datCombined); # If you want, you can over-ride no.obs for Z-score computations here/and below.

   # INVAR: we have pm, set true no.obs rows
   pm$no.obs = nrow(datCombined);

   if (nrow(datCombined) < 5) {
      stop("datCombined must have 5 or more rows (individuals) for neo() to work.")
   }

   # before any columns get shuffled, translate the pm$only.score.edges.into into column names
   # so we can re-match them later during filtering.
   if (!is.null(pm$only.score.edges.into)) {
      pm$only.score.edges.into = colnames(datCombined)[pm$only.score.edges.into]
   }

   # auto-farm out permutations to SGE if requested
   if (pm$run.perms.on.sge & pm$number.BLOCK.permutations > 0) {
      # save the current environment
      # after changing pm$run.perms.on.sge to avoid infinite loop
      pm$run.perms.on.sge = FALSE

      # and give sge jobs their own directory, inside the main log directory
      pm$sge.run.dir = create.logpath.enclosed.sge.directory(pm)
      pm.print(pm,paste("Created directory for Sun Grid Engine job output: ",pm$sge.run.dir))      

      logfile.to.restore = pm$neo.log.file
      pm$neo.log.file = "sge_neo."

      # temporarily change directory into sge run directory
      cwd=getwd()
      setwd(pm$sge.run.dir)
      
      save(list=ls(), file="master.sge.image.rdat")

      # now that sge jobs have their own logfile directory, get back ours!
      pm$neo.log.file = logfile.to.restore

      # function neo.sge.perm loads the saved file "master.sge.image.rdat" and runs one permutation (could be more).
      #n.perm.job=ceiling(pm$number.BLOCK.permutations/10)

      n.perm.job=pm$number.BLOCK.permutations
      #job.list=rep("neo.sge.perm",n.perm.job)
      job.list=rep("neo.sge.perm",n.perm.job)

      pm.print(pm,paste("Submitting ",n.perm.job," permutation jobs to Sun Grid Engine."))
      submit.to.sge(job.list=job.list,
                    neo.install.path=pm$sge.needs.neo.install.path,
                    r.program.install.path=pm$sge.needs.r.program.install.path)

      setwd(cwd) # restore directory, now that jobs are submitted.
      
      pm$number.BLOCK.permutations = 0
    }
   pm$doing.permutation = FALSE; # first time through

   # Implement permutations loop

   # handle permutations in a loop
   pm$perm.count = 1;
   while(1) {
      # permutation loop will end with a break;
      # at end of this while loop, put: if (pm$perm.count > pm$number.BLOCK.permutations) { break; # we are done} else { # permute data and repeat}
   
   # skip initial data analysis on unpermuted data if only doing a permutation
   if (pm$perm.count == 1 & pm$just.do.permutation.no.unpermuted) {
       # then skip the first unpermuted data run, and go straight down to permuting the data
   } else {  # begin skip.first.time.test

   # Implement batch processing for all processing...
   pm$neo.batch.run = TRUE
   pm$first.batch.pass = TRUE
   pm$last.batch.pass = FALSE

   # IMPORTANT: Inside the batch loop, we modify pm$A and pm$B, as well
   # as pm$A.cn and pm$B.cn,  but we don't change
   # the fully enumerated lists, pm$enumerated.A.cn and pm$enumerated.B.cn

   batches = length(pm$enumerated.A.cn)
   for (batch.num in 1:batches) {
   tryCatch({ # Prevent NEO from aborting mid run
        k=batch.num

        print(paste(sep="","*** Edge ",batch.num," of ",batches,
         ", running: ",pm$enumerated.A.cn[k]," - ",pm$enumerated.B.cn[k]," ...",date()));
   if (batch.num > 1) pm$first.batch.pass = FALSE
   if (batch.num == batches) pm$last.batch.pass = TRUE

   traitA.cn = pm$enumerated.A.cn[k]
   traitB.cn = pm$enumerated.B.cn[k]
   traitA = match(traitA.cn, colnames(datCombined))
   traitB = match(traitB.cn, colnames(datCombined))

   dc = datCombined[,c(snpcols,traitA,traitB)]
   my.snpcols = 1:length(snpcols)

   pm$A = trA.loc = ncol(dc)-1
   pm$B = trB.loc = ncol(dc)

   pm$A.cn = colnames(dc)[pm$A] # keep consistency
   pm$B.cn = colnames(dc)[pm$B] # keep consistency
   
   # main four calls

   #
   # 1) initial SNP primary is not the final SNP election
   #

   # modified for batch, was:  r1 = neo.snp.primary(datCombined, snpcols=snpcols, traitcols=traitcols, 
   r1 = neo.snp.primary(datCombined=dc, snpcols=my.snpcols, traitcols=c(trA.loc,trB.loc), 
                        pm=pm, 
                        skip.snp.select=pm$skip.snp.select,
                        neo.log.file=neo.log.file);
   if (!pm$no.log && pm$save.intermediate.files) save(r1,file=paste(sep="",neo.log.file,".batch_",batch.num,".r1.rdat"));

   # if requested by setting pm$number.BLOCK.permutations > 0, we evaluate 
   # the choice of SNPs by shuffling (generating a permutation of) the 
   # assignment of SNPs to cases (individuals, mice), and then 
   # re-selecting SNPs.
   #
   # For each SNP chosen, we can then assign a permtation frequency: how often
   # the k-th ranked SNP exceeded the value observed in the k-th ranked SNP in
   # the actual data.

   #
   # 2) variables go through two more processes...supernoding and filtering
   #
   r2 = neo.create.super.nodes(r1);
   if (!pm$no.log && pm$save.intermediate.files) save(r2,file=paste(sep="",neo.log.file,"r2.rdat"));

   #
   # 3) the final "election" happens after applying various filters
   #
   r3 = neo.edge.filters(r2,ignorable.edge.list);
   if (!pm$no.log && pm$save.intermediate.files)   save(r3,file=paste(sep="",neo.log.file,"r3.rdat"));

   #
   # 4) edge orient. Walk two steps down from each SNP, and integrate 
   # SNPs signals by the various edge orienting strategies.
   #
   r4 = neo.edge.orient(r3,
                        fit.full.model=FALSE,skip.LEO=skip.LEO)
   if (!pm$no.log && pm$save.intermediate.files)  save(r4,file=paste(sep="",neo.log.file,"r4.rdat"));

   #    Excel sheet should be opened automatically 
   #    at end of summary, if not, set pm$excel.path correctly.
   z = summary(r4);
   }, error = function(e) {
     warning("NEO crashed fatally when attempting to orient edges between ",
             pm$enumerated.A.cn[k], " and ", pm$enumerated.B.cn[k],
             " with error:\n", e, "\n",
             " skipping this edge.")
     if (k == 1) {
       stop("NEO crashed when orienting first edge. No point continuing because",
            " the final CSV file will be malformed.\n",
            "Try shuffling your non-SNP columns so the offending edge does not",
            "occur first.")
     }
   })

 } # end batch.num loop

 } # end skip.first.time.test, implemented in case we are only doing a permutation...

  # process permutations, if any.
  if (pm$perm.count > pm$number.BLOCK.permutations) { 
       if (pm$number.BLOCK.permutations) {print(paste("done with all ",pm$number.BLOCK.permutations," permutations, at ",date()))}
       break; # we are done (thus the default of perm.count=1 when pm$number.BLOCK.permutations=0
     } else {

       # permute data and repeat
       pm$perm.count = pm$perm.count+1

       pm$no.log = TRUE
       pm$quiet = TRUE
       pm$save.intermediate.files=FALSE
       pm$doing.permutation = TRUE

       print(paste("Starting on ",pm$perm.count-1," out of ",pm$number.BLOCK.permutations," permutations, at ",date()))

       # save permutations to their own file
       pm$perm.file = paste(sep="",pm$neo.log.file,".perm.csv")


       # re-set the permute colnums in case of re-arrangement of datCombined
       pm$block.to.permute = match(pm$block.to.permute.cn,colnames(datCombined))
       datCombined=permute.snps.traits(pm, datCombined, pm$block.to.permute)
       print(paste(".............permutation of SNP:gene data done....starting NEO evalutaion",date()))

       if (pm$perm.count > pm$number.BLOCK.permutations) { pm$open.excel.now = TRUE }

  } # conclusion/goto top of while(1) loop.

  } # end permutation loop


 #
 # possible halting point: if just doing permutations, skip the graphing and the messing with the excel file output.
 #
 if (pm$just.do.permutation.no.unpermuted) {
    pm.print(pm,paste("Done with permutation and pm$just.do.permutation.no.unpermuted=TRUE. Halting at ",date()))
    return(invisible(c()))
 }
   
   # read in permutations and summarize in .csv file for Excel. (Assuming we also computed the observed stats)
   if (pm$doing.permutation & !pm$just.do.permutation.no.unpermuted) create.permutation.report(pm,colnames(datCombined),snpcols)

   # read back in the final matrix to get the whole thing
   # for graph drawing.

   # only go past here if we have a valid excel file...
   excel.file=get.excel.file(pm)

   if(!file.exists(excel.file)) {
      pm.print(pm,paste(sep="","Could not load excel file: '",excel.file,"'. (Probably no pm$neo.log.file set. Inspect pm$neo.log.file='",pm$neo.log.file,"') Finishing before graphing or (trying to) open spreadsheet."))
      return(invisible(z))
    }

   # set a default so we know if we got it by checking !is.null(fullcsv), even if we don't do graph drawing...
   fullcsv=NULL
   
   if (1) { # to comment in/out graph drawing.

      # don't insit on rownames=1 b/c snp edges may get duplicated
      fullcsv = read.csv(file=excel.file,header=T,sep=",",quote="\"",dec=".",fill=T,comment.char="")
      w.have.leo=which(!is.na(fullcsv$LEO.NB.OCA))

      zGRAPH$pm = pm
      non.snp.cols=setdiff(1:ncol(datCombined),snpcols) # not sure why, but this variable not in scope now. So re-define.
      zGRAPH$non.snp.cols=non.snp.cols
      zGRAPH$coor = do.graph.placement(snpcols,zGRAPH$numsnps, zGRAPH$nc, zGRAPH$num.non.snps, non.snp.cols)

      traits=unique(c(pm$enumerated.A.cn,pm$enumerated.B.cn))

      # the matrix to graph; get this from the final (accumulated after batches) excel file, now in fullcsv
      zGRAPH$M = matrix(nrow=zGRAPH$nc, ncol=zGRAPH$nc,data=0); 
      colnames(zGRAPH$M)=rownames(zGRAPH$M) = colnames(datCombined)

      zGRAPH$final.eo = zGRAPH$M # also make a score matrix, for the labelling of the graph
      zGRAPH$final.eo[] = NA # start unlablled though.

      non.snp.cols=setdiff(1:ncol(datCombined),snpcols)


      rn=as.character(fullcsv[,1]) # rownames now in first column
      splt=strsplit(rn,split=" -> ",fixed=TRUE)
      FROM.A=sapply(splt,function(x) x[1])
      TO.B  =sapply(splt,function(x) x[2])
      
      # create the graph
      for (i in w.have.leo) {
          if (fullcsv$LEO.NB.OCA[i] > 0) {
                zGRAPH$final.eo[FROM.A[i],TO.B[i]] = fullcsv$LEO.NB.OCA[i]
                zGRAPH$M[FROM.A[i],TO.B[i]] = 1
                eval(parse(text=as.character(fullcsv$Final.SNPs.LEO.NB.OCA[i]))) # sets fsnp.cpa,fsnp.oca
                if (!is.null(fsnp.cpa)) zGRAPH$M[fsnp.cpa,FROM.A[i]]=1
                if (!is.null(fsnp.oca)) zGRAPH$M[fsnp.oca,TO.B[i]]=1
          }
      }

      zGRAPH$my.title=pm$run.title
      
      # draw the graph; use adjust=TRUE to click with mouse and 
      # reposition the nearest node.
      # Skip drawing on unix b/c x11 frequently doesn't have the fonts.
      #
      if (!pm$quiet & !is.unix()) { z$final.coor = neo.graph(zGRAPH,adjust=repos.nodes);}
   } # end comment in/out graph drawing


   # Since z now just has the stats for last edge evaluated, and some of the earlier
   # scripts (before we implemented batch processing for scalability)
   # expect there to be multiple edges evaluated in z$stats, we'll try
   # to keep back-compatibility by reading in the full set of edges
   # that were evaluated, and assigning them to z$stats
   #

   if (!is.null(fullcsv)) {

        excel.file=get.excel.file(pm)

        # try to sort first, if possible. Use check.names=FALSE to prevent ___ -> X___
        #
        excel.file.to.sort = read.csv(file=excel.file,header=T,sep=",",quote="\"",dec=".",fill=T,comment.char="",check.names=FALSE)


        # old: sorted.df=excel.file.to.sort[order(excel.file.to.sort$LEO.NB.OCA,decreasing=TRUE),]
        # new:
        sorted.df=excel.file.to.sort[order(excel.file.to.sort$P.weighted.LEO.NB.OCA,decreasing=TRUE),]

        colnames(sorted.df)[1]="edge"

        z$stats = sorted.df

        write.csv(sorted.df,file=excel.file,na ="")

       if (pm$open.excel.at.end) {
            open.excel.file(excel.file,pm)
       }
    }

      # re-save, with the coordinates, and the fullcvs spreadsheet in z$stats
   if (!pm$no.log && pm$save.intermediate.files)  save(z,file=z$post.save.fname);

   invisible(z)

} # end neo() top level wrapper function


# PARAMETER CHOICES DOCUMENTATION FOLLOWS!
#
# neo.get.param(): Parameter defaults for neo functions, neo.get.param()
#
# give the user sane default parameters and let them avoid having to type all 
# the defaults out each time
#
# logpath gets pre-generated and passed in each time neo() calls.
#
if(exists("neo.get.param") ) rm(neo.get.param);
neo.get.param=function(logpath="neo.logfile") {
   list(
        graph.arrow.edge.th=.3, # an edge score must surpass this threshold to be graphed
                      # by neo.graph,
                      # (and kept if using ZEO to orient edges)

        # a pair with correlation below this are considered independent.
        # This parameter is only used for LCD, RVV, and the max.max special
        # case handling; not used in OCA and CPA computations.
        cor.ind.th=0.2, 

        # a pair with correlation above this are considered dependent.
        # The absolute correlation between a marker and a trait has
        # to surpass this threshold to be considered at all (including in OCA and CPA).
        cor.dep.th=0.1, 

        # if the partial correlation falls below this upon conditioning,
        # the edge is eliminated by a filter.
        #pcor.th=0.05, 
        # Not relevant for OCA and CPA, and in particular is turned off
        # by pm$skip.pcor.filters.currently, which is TRUE by default.
        pcor.th=0,  # turn off this filter.

        # threshold for Cohen's omega filter. We keep small omega edges 
        # (below omega.th), because this implies small indirect influences
        # between a pair of variables.
        #omega.th = 0.85, 
        omega.th = 1, # 1 will effectively turn off this pcor threshold filters.

        k.pcor.th=-0.05, # negative value for k.pcor.th will turn off the 
                         # VCdep loops-- which take a long time, and are 
                         # only used if no SNPs are availabe; the negative 
        # default makes for faster runs. In general, k.pcor.th
        # is the point at which: if the partial correlation between a
        # pair when conditioned upon a single variable k falls below
        # the k.pcor.th value, then the edge between the pair is filtered out.

        induced.dep.th = .2, # if checking for conditional dependence 
        # (v-structure location),
        # --that forms between independent parents of a common child,
        # we require that the induced correlation (upon conditioning 
        # on the prospective child) reach at least this level.

        snp.below=.3,   # if SNPs are correlated at 1-snp.below or higher, 
        # then they are clustered
        # together and a single representative from the cluster is chosen to
        # represent the Super(SNP)node. Only relevant if pm$do.snp.pick.hcluster = TRUE,
        # but by default currently it is FALSE.

        cor.above.supernode.def=.7, # if PCs are correlated at this value or 
        # higher, then
        # they are clumped into a supernode. Variables within
        # a supernode are not allowed to filter each other's edges.

        force.cutree.count=25, # when picking SNPs from the hcluster-ing of 
        # all SNPs, cut the
        # SNP tree into this many branches sets. Within each set (cluster)
        # the SNP with the highest correlation with any non.snp trait or PC
        # is chosen as a representative. This is an attempt to get some 
        # orthogonal SNPs into the playing field.

        top.N.snps.per.trait=4, # During greedy SNP selection, pick this 
        # many top SNPs in terms of
        # the SNPs that are correlated with each pc or trait.

        snp.redundancy.cor.level=.9, # retract redundant SNPs at this level 
        # of abs correlation

        skip.greedy.snps.if.redundant = TRUE, # if TRUE,
        # don't even bother selecting subsequent greedy snps that don't
        # differ from previous snps by having a correlation below the
        # pm$snp.redundancy.cor.level set above.

        # At or above this correlation level, columns in datCombined
        # are flagged as too multi-collinear to continue. This check
        # is implemented in check.multi.collinearity().
        # (Must be positive; negative values will have abs() applied first.)
        datCombined.column.redundancy.reject.cor.level = .95, 

        warn.if.zeo.and.leo.disagree = FALSE, # if TRUE,
        # then issue a warning as we notice that the ZEO and LEO scores
        # differ in sign.

        show.snp.selection.log = FALSE, # if TRUE,
        # show the details of why the Greedy SNPs were chosen the way they were.
        

        skip.snp.retraction=TRUE, # if TRUE, we don't retract redundant SNPs
        # correlated with
        # previously chosen SNPs at or above snp.redundancy.cor.level

        turn.NA.scores.to.zero.and.zero.to.666 =FALSE, # when SNPs aren't
        # correlated with traits, we ususally just
        # report NA for the edge scores without enough information to be 
        # scored. However, if this option is TRUE, we report 0 instead.
        # And formerly zero scores are turned to 666e-8. Implemented per
        # Steve's request.

        no.obs=NULL, # The true value of nrow(datCombined), for when it is
        # needed.

        no.obs.Z=NULL, # if set, over-rides nrow(datCombined) in the 
        # Z-score calculations

        # The user should supply data with NA already imputed.
        # We hacked in some solutions, but probably these are not optimal.
        # They are only used if this is not FALSE.
        rough.and.ready.NA.imputation = FALSE,

        # when fitting multiple marker models in local.sem.four.var.m1m2(), include
        # covariance terms within the M_A set and within the M_B set?
        # FALSE is traditional; TRUE may give better or more accurate model fits
        # because we aren't penalized for some of the markers being correlated
        # with others in their group. Now confirmed, we default to TRUE now.
        add.MA.MB.covar.in.local.sem.four.var.m1m2 = TRUE,

        use.leo=TRUE,
        # Use the LEO (Local SEM Edge Orienting) method if true, otherwise 
        # use the ZEO to draw the final graph and set the "keep" flag in
        # the output spreadsheet/dataframe. Note that LEO is still calculated,
        # even if it won't be used (for comparison). If you wish to circumvent
        # the LEO computations, use instead the flag skip.LEO=TRUE
        # when called neo() or neo.edge.orient().

        # Edge Orienting method to use (to draw graph from)
        #
        # Use the following for eo.type to specify the edge orienting used
        # in drawing the final graph in neo.graph from z$M.
        # Options: "all.snp.permutations.equal.votes", "forward.step.selection" or # "max.vs.max"
        eo.type = "forward.step.selection", # "all.snp.permutations.equal.votes", # default

        # Where Excel is installed.
        excel.path="C:\\Program Files\\Microsoft Office\\Office10\\EXCEL.EXE",

        # parameters for using the fisher transform of the correlation 
        # coefficient to (semi)normality.
        fisher.dof.cor=3, # 3 usually
        fisher.dof.pcor1=4, # 4 usually

        # the fit function sem (from R package sem) crashes less often with
        # correlation matrices instead of covariance matrices.
        sem.fit.correlation.instead.of.covariance = TRUE,

        # Control how deep the forward step-wise regression to select SNPs goes.
        # This is the number of steps beyond the first greedy step. For instance,
        # a value of 1 here means that 2 SNPs will be chosen (assuming they are
        # available) into A, and 2 SNPs chosen into B, for each A-B edge oriented.
        # A value of 0 is allowed; in which case if do.snp.pick.fstep==TRUE then
        # the same SNP as the first greedy SNP will be selected.
        orthog.search.depth = 5,

        # Over-ride orthog.search.depth with a smaller number if this is not NULL.
        # Only applied in OCA model building at present (losem.walk.two.steps.weighted.mean.mlreg).
        # Also, this gets applied BEFORE the final forward stepwise regression,
        # so this will likely be the strongest marginal set...not necessarily
        # what you always want, but good for emulating manual selection...
        # Applied after marker assignment consistency is obtained.
        max.snps.in.OCA.CPA.sets = 3,

        #logfile dir/path.prefix
        neo.log.file=logpath,

        # may wish (with multiple runs under SGE especially) to re-use existing
        # log directory....in which case neo.log.file is kept, and is assumed
        # to be a file path prefix that points into an already existing directory.
        # Otherwise (if TRUE) we generate such a path prefix by combining pm$run.title
        # and a time stamp.
        keep.existing.neo.log.file=FALSE,

        # for SunGridEngine runs, these need to be specified...
        sge.needs.neo.install.path="/home/jaten/dev/peculiar/neo.txt",
        sge.needs.r.program.install.path="/usr/local/bin/R",

        # if quiet is TRUE, we endeavor not to draw graphs and open excel 
        # spreadsheets and display progress.
        quiet=FALSE,

        # number.BLOCK.permutations is > 0, then we re-run the analysis (much 
        # slower) the given number of times, each time shuffling the 
        # relationship between cases and this block of columns (defaults to SNPs--see next parameter).
        number.BLOCK.permutations = 0,

        # Specify the block of column numbers to permute.
        # Typically this would be the snpcols, or the traitcol or some genecol,
        # depending on which hypothesis of independence you wish to evaluate.
        # DEFAULT: use snpcols if left to NULL.
        block.to.permute = NULL,

        # optional additional blocks to permute (must be non-overlapping
        # with each other and the first block).
        block2.to.permute = NULL,
        block3.to.permute = NULL,

        # options to permutation runs...

        # TRUE means ignore the original data and just permute data and run neo.
        just.do.permutation.no.unpermuted = FALSE,

        # if a run on a sun grid engine based cluster, we
        # can parcel out the permutation jobs to different cluster nodes
        # for faster execution, then collect at the end.
        run.perms.on.sge = FALSE,


        # impute.nn() will use as nearby SNPs to impute, only those SNPs that are correlated
        # at this minimum level or higher
        impute.nn.minimum.correlation.th = .8,

        # Do final collider computation at each v-structure inferred in the 
        # network. Or not.
        do.final.colliders = FALSE,

        # skip the SNP selection and just use all provided SNPs?
        skip.snp.select = FALSE,

        # should we save after each stage for diagnostics/recovery during
        # a long run?
        save.intermediate.files = FALSE, # make batch processing go faster
        #save.intermediate.files = TRUE, 

        # edges to ignore; ex: pm$ignorable.edge.list=list(c(1,2),c(2,3)) 
        # to ignore edges between columns 1 & 2 (and 2 & 1)
        # and also ignore any edges between columns 2 and 3.
        ignorable.edge.list = NULL,

        # if specified, only score edges into this column or columns of 
        # datCombined
        only.score.edges.into = NULL,

        minimum.minor.allele.freq = NULL,

        run.title = "", # name of the run/experiment

        no.log = FALSE,   # whether to log or not

        # set this to TRUE if you want to mix single marker scores with 
        # double marker model
        # scores, even though they are on very different scales and not 
        # really comparable.
        work.hard.at.EO.contingencies.if.not.both.M.B.and.M.A.markers = FALSE, 

        # shall we do Minority Reporting? If so, we compare max,for,all scores
        # and issue a warning if they differ in sign.
        do.minority.report.max.for.all = FALSE,

        # m1m2.average: This is a very, very slow score, the average over all the m1m2
        # models, using each of the M.A markers paired with each of the M.B
        # markers at a time.  Typically off (FALSE) by default because it's slow.
        do.m1m2.average = FALSE,
#        do.m1m2.average = TRUE,

        do.leo.max = TRUE, # fit and compare MA->A->B<-MB models for best single marker MA and MB
        do.leo.all = TRUE, # fit and compare the average of MA->A->B models

        use.traitcols.guess = FALSE, # Turn back on if you want neo
              # to use trait names that start with T as clinical traits by defaulat

        zeo.proxy.for.leo.th = NA, # if NA, compute all LEO scores. If not NA,
        # use this threshold above which
        # we compute LEO scores (slow), below which the ZEO (fast) score is
        # saying the LEO is probably not worth computing. Since ZEO is more
        # than 100 times faster than LEO, using this proxy yields massive 
        # speedups. A 2.0 threshold is more realistic and conservative, but 
        # we want to not miss subtle signals too.
        # So we would default to a somewhat low proxy score of 0.

        do.snp.pick.greedy = TRUE, # if TRUE, SNP selection includes the 
        # greedy SNPs
        do.snp.pick.fstep  = TRUE, # if TRUE, pick SNPs by (Gram-Schmidt) 
        # forward step-wise regression
        do.snp.pick.hcluster = FALSE, # if TRUE, also cluster the SNPs and 
        # include the single SNP from each cluster with highest trait 
        # correlation
        
        A = NULL, # if both A and B are non-NULL, we only check edges between
        # the A set and B set of column numbers.
        B = NULL,

        just.AB = FALSE, # if true (which is useful for big runs using A, B specification)
                        # then we don't try to create the big square matrices M, af,
                        # etc, and in this way we can avoid an out-of-memory crash.
                        # just.AB will be set to true automatically by the program
                        # if A and B are non-null.

        skip.pcor.filters.currently = TRUE, # turn off the rarely used partial correlation
                                           # filters, especially because they will crash
                                           # on a big batch multi-marker pm$A, pm$B set run.

        neo.batch.run = TRUE,   # flag set to true when running a large job
                                # whose results won't all stay in memory, and so must be written out to
                                # disk as we go. 
        first.batch.pass = TRUE, # used to track batch processing...do we need to write the headers?
                                # Should be TRUE only the first batch of edges processed.
        last.batch.pass = FALSE, # set to false on the last pass to signal to open the Excel file.

        # chromsome.wide.strong.marker.consisteny:
        #
        # If TRUE, when enforcing marker consistency (all M.A markers have stronger correlation
        # with A than with B, and vice-versa for M.B), let the strongest trait-marker
        # correlation within each chromosome win that entire chromosome for the trait (A or B),
        # and subsequently don't let the other trait (B, if for instance A wins) have any markers
        # on that whole chromosome.
        #
        # This prevents semi-overlapping peak eQTL from having M.A and M.B both
        # end up with markers nearby from the same chromosome. This happened in
        # the BxH.Apoe-null data for females with Insig1 and B4301..Riken on 
        # chr8:100Mb or so.  It allows automatic SNP selection to model what
        # would likely happen under manual SNP curation more closely.
        #
        # If TRUE, the marker labels must include the chromsome number after '.chr' as
        # in gene.mmt23434.chr1.Morgrans0.01,  or geneB.mmt23223.chrX.Morgrans0.232.
        #
        chromsome.wide.strong.marker.consistency = TRUE,


        # Forcing SNP choices for M.A or for M.B. Specify the column numbers of the snps within datCombined.
        # If set, then automatic SNP selection will be ignored.
        forced.MA.colnum = NULL,
        forced.MB.colnum = NULL
        
        )
}

#
# install required libraries, if not already...


repos=getOption("repos"); repos["CRAN"] = "http://cran.stat.ucla.edu";
options(repos=repos);
# to choose a different mirror: 
# chooseCRANmirror()

# to auto-install libraries, if missing, when first sourced.
if (!require(GeneCycle,quietly=TRUE)) install.packages("GeneCycle");
if (!require(GeneNet,quietly=TRUE)) install.packages("GeneNet");
if (!require(ggm,quietly=TRUE)) install.packages("ggm");
if (!require(e1071,quietly=TRUE)) install.packages("e1071");
if (!require(sem,quietly=TRUE)) install.packages("sem");
if (!require(MASS,quietly=TRUE)) install.packages("MASS");
if (!require(matrixStats,quietly=TRUE)) { install.packages("matrixStats"); library(matrixStats);}

# try to pre-empt Excel opening bug...
#if (.Platform$OS.type == "windows") {
#   if (!require(rcom,quietly=TRUE)) {
#      install.packages("rcom"); # try to pro-actively prevent Excel-VBA error in R-2.5 and R(D)COM-2.5
#      library(rcom)
#      comRegisterRegistry()
#      comRegisterServer()
#   }
#}

#if (!require(Matrix,quietly=TRUE)) install.packages("Matrix"); # for sparse Matrix


#install.packages(c("GeneTS","ggm","e1071","sem")); 
# Could also: use ignoreWarnings() to ignore warnings if we already have the packages.

#library(GeneCycle) #for pcor.shrink to compute conditioning on all remaining vars
#library(GeneNet) #for pcor.shrink to compute conditioning on all remaining vars
#library(ggm) # for pcor to compute conditioning on individual variables k
#library(e1071) # for impute: to handle NA in the data by filling in the median
#library(sem) # for RAM model format 'mod' and Maximum likelihood indices of fit.
#library(MASS) # for stepAIC
#library(Matrix)

library(rpart) # CART models, comes with R.


# function to run a bunch of (multimaker/neo) permutations under sun.grid.engine
#
# some other function should have save datCombined,snpcols,traitcols,pm 
# with something like: save(list=ls(name=parent.frame()), file="master.sge.image.rdat")
# before calling this one.
#
if(exists("neo.sge.perm") ) rm(neo.sge.perm);
neo.sge.perm=function() {

  # load our workspace
  load("master.sge.image.rdat")

  # We do use our own dedicated directory -- but this should be set before qsub-ing.

  # replace old variable with new value
  sge.task = get.sge.task.or.zero()

  # only permutations, no observed stats
  pm$just.do.permutation.no.unpermuted = TRUE

  pm$number.BLOCK.permutations = 10 # do 10 permutations per sge task (can be more if desired)
#  pm$number.BLOCK.permutations = 1

  pm$quiet=TRUE
  pm$no.log=TRUE

  # just in case we do turn logging on; don't trample the original.
  if (sge.task > 0) {
     pm$neo.log.file = paste(sep="",pm$neo.log.file,".sgetask.",sge.task)
     pm$keep.existing.neo.log.file = TRUE
  }

  print(paste("Starting neo inside sge run...pm$number.BLOCK.permutations=",pm$number.BLOCK.permutations," at ",date()))
  z=neo(datCombined=datCombined,snpcols=snpcols,traitcols=traitcols,pm=pm)
  invisible(z)
}




# sma.neo.sge.perm: function to run a bunch of single.marker.analysis permutations under sun.grid.engine
#
# some other function should have save datCombined,snpcols,traitcols,pm 
# with something like: save(list=ls(name=parent.frame()), file="permute.sma.master.sge.image.rdat"))
# before calling this one.
#
# This function can't take any arguments. All data is loaded from "permute.sma.master.sge.image.rdat"
#
if(exists("sma.neo.sge.perm") ) rm(sma.neo.sge.perm);
sma.neo.sge.perm=function() {

  # unique timestamp for produced files.
  one.date=paste(sep="",gsub(":",".",gsub(" ","_",paste(sep="",date()))))

  # load our workspace; including pm
  load("permute.sma.master.sge.image.rdat")

  # We do use our own dedicated directory, but this should be set before qsub-ing

  # replace old variable with new value
  sge.task = get.sge.task.or.zero()

  # only permutations, no observed stats
  pm$just.do.permutation.no.unpermuted = TRUE

  # just in-case, to avoid infinite loop
  pm$run.perms.on.sge = FALSE

  # run just one permutation per SGE job
  pm$number.BLOCK.permutations = 1
  num.perm = pm$number.BLOCK.permutations # reuse code below

  pm$quiet=TRUE
  pm$no.log=TRUE

  # just in case we do turn logging on; don't trample the original.
  if (sge.task > 0) {
     pm$neo.log.file = paste(sep="",pm$neo.log.file,".sgetask.",sge.task)
     pm$keep.existing.neo.log.file = TRUE
     sma.perm.file = paste(sep="","sgetask.",sge.task,".",sma.perm.file)
  }

  print(paste("sma.neo.sge.perm(): Starting single.marker.analysis inside sge run...pm$number.BLOCK.permutations=",pm$number.BLOCK.permutations," at ",one.date,"; output to dir: ",getwd(),"  file: ",sma.perm.file))


   # Now start the permutations.

   # execute SMA
   for (i in 1:num.perm) {
      print(paste("single.marker.analysis permutations, doing : ",i," out of ",num.perm," at ",date()))

      datCombined=permute.snps.traits(pm=pm, datCombined, block.to.permute, save.perm.to.file=FALSE)

      sma.pre = single.marker.analysis(datCombined=datCombined,snpcols=snpcols,genecols=genecols,traitcols=traitcols,leo.i.th=leo.i.th, leo.o.th=leo.o.th, leo.nb.th=leo.nb.th,pm=pm,use.ranks=use.ranks,build.multi.marker.to.gene.model=build.multi.marker.to.gene.model, impute.na=impute.na)

      # write out just statistics of interest...actually keep most everything now, adding :ncol
      sma = sma.pre[,c(1,1,3,5,6,7,8,9,11,12:ncol(sma.pre))] 
      colnames(sma)[1]="model"
      colnames(sma)[2]="M"
      colnames(sma)[3]="A"
      colnames(sma)[4]="B"

      sma[,2] = as.character(sma[,2])
      sma[,3] = as.character(sma[,3])
      sma[,4] = as.character(sma[,4])

      # create unique row identifier
      sma[,1]=paste(sep="->",sma[,2],sma[,3],sma[,4])

      if (i==1) {
         write.table(sma, file=sma.perm.file, col.names = NA, sep=",",dec=".",qmethod="double")
      } else {
         write.table(sma, file=sma.perm.file, append=TRUE, col.names = FALSE, sep=",",dec=".",qmethod="double")
      }
   } # end i loop

  invisible(c())
}


# seed.root.leaf.neo: function to go
# over the whole processing of seeding and then 2 or 3 marker (pm$max.snps.in.OCA.CPA.sets
# parameter controls this) rooting,
# followed by leafing and neo multimarker analysis.
#
# Basically an automation of the process described in
# Methods_Insig1.Supplement.11Sept2007.doc
# so that we can also do permutation testing on
# this proceedure.
#
# INVAR: known.downstream and known.upstream are single genes
# (might relax this later).
# 
if(exists("seed.root.leaf.neo") ) rm(seed.root.leaf.neo);
seed.root.leaf.neo=function(pm, datCombined, snpcols, known.upstream, known.downstream, genecols, require.sma.model.prob.over=.89, require.sma.leo.nb.over=.75, num.seed.snps=2, require.cor.with.known.up.over=.5) {

   cn=colnames(datCombined)
   w.up = known.upstream # match(known.upstream,cn)
   w.dn = known.downstream # match(known.downstream,cn)

   if (w.up > ncol(datCombined) | any(is.na(w.up))) stop("bad known.upstream specification")
   if (w.dn > ncol(datCombined) | any(is.na(w.dn))) stop("bad known.downstream specification")

   # seed: find top loci on different chromosomes that look like
   #       causal drives of w.up -> w.dn
   # 
   sma.seed = single.marker.analysis(datCombined=datCombined,snpcols=snpcols,genecols=w.up,traitcols=w.dn,pm=pm)
   ord.seed = sma.seed[order(sma.seed$leo.nb.AtoB ,decreasing=TRUE),]

   # take the top distinct markers...into uniq.mark
   markers = as.character(ord.seed[,1])
   uniq.mark = unique(markers)

   # take out the M: annotation in front...so "M:snp" becomes "snp"
   uniq.mark.cn = sapply(strsplit(uniq.mark,"M:",fixed=TRUE),function(x) x[2])
   
   # find the chromosomes, so we can pick markers on distinct chromosomes...
   chrm1=sapply(strsplit(uniq.mark,".chr",fixed=TRUE),function(x) x[2])
   chrm=sapply(strsplit(chrm1,".",fixed=TRUE),function(x) x[1])
   if (all(!is.na(chrm)) & all(chrm != "")) { got.chrom = TRUE; }

   # get the set of unique chromosomes, in order.
   uniq.chrm = unique(chrm)
   
   chosen.seed.snp = rep(NA,num.seed.snps)
   chosen.seed.snp[1]=uniq.mark.cn[1]
   
   if (num.seed.snps > 1) {
      for (i in 2:num.seed.snps) {
         chosen.seed.snp[i]=uniq.mark.cn[min(which(chrm == uniq.chrm[i]))]
      }
    }
   
   uniq.mark.w = match(chosen.seed.snp ,cn)
   uniq.mark.cn = cn[uniq.mark.w]
   
   if (any(is.na(uniq.mark.w))) { stop(paste("problem locating uniq.mark.cn genes in datCombined, these being:",uniq.mark.cn))}

   
   # pick out genes highly correlated with known.upstream to begin with
   use.genecols = setdiff(genecols, w.up)
   cx2=cor(datCombined[,w.up],datCombined[,use.genecols],use="p")
   known.up.cor.genes=colnames(cx2)[which(abs(cx2)[1,]>= require.cor.with.known.up.over)]
   w.known.up.cor.genes = match(known.up.cor.genes, cn)

   # now, any good SMA results for these highly correlated genes using the chosen markers...
   require.mlogp.below = -log10(require.sma.model.prob.over)
   #sma.res=list()
   high.leo.and.model.prob=c()
   for (i in 1:length(uniq.mark.w)) {
      sma=single.marker.analysis(datCombined=datCombined,snpcols=uniq.mark.w[i],genecols=w.up,traitcols=w.known.up.cor.genes,pm=pm)
      ord.sma=sma[order(sma$leo.nb.AtoB,decreasing=TRUE),]
      # sma.res[[i]]=ord.sma
      high.leo = which(ord.sma$leo.nb.AtoB >= require.sma.leo.nb.over)

      high.model.prob = which(ord.sma$mlogp.M.AtoB <= require.mlogp.below)

      both.criteria = intersect(high.leo, high.model.prob)
      
      pass.criteria.genes = sapply(strsplit(as.character(ord.sma[both.criteria,5]),"B:",fixed=TRUE),function(x) x[2])
      w.pass.criteria.genes = match(pass.criteria.genes, cn)

      # only have to look good in response to one of the markers.
      high.leo.and.model.prob = union(high.leo.and.model.prob, w.pass.criteria.genes)

      if (!pm$no.log) {
         save(ord.sma, pass.criteria.genes, file=gsub(" ","_",paste(sep="",pm$neo.log.file,".seed.sma.for.",uniq.mark.cn[i],".to.",cn[known.upstream],".rdat")))

      }
    }

   high.leo.and.model.prob.cn = cn[high.leo.and.model.prob]

   if (!pm$no.log) {
       save(list = c("ord.seed","high.leo.and.model.prob.cn"), file=gsub(" ","_",paste(sep="",pm$neo.log.file,"post.sma.ord.seed.for.known.upstream.",cn[known.upstream],".rdat")))


      my.file = gsub(" ","_",paste(sep="",pm$neo.log.file,".total.sma.leafings_num.pass.genes.is_",length(high.leo.and.model.prob),"_.for.known.",cn[known.upstream],".pass.criteria.genes.rdat"))
      save(high.leo.and.model.prob.cn, file=my.file)
      pm.print(pm,paste("Saved results to ",my.file))
   }
   pm.print(pm,paste("Number of genes passing criteria: ",length(high.leo.and.model.prob)))

   # and run a neo analysis on these
   high.leo.and.model.prob = setdiff(high.leo.and.model.prob,w.up) # just in case, so we don't overlap
   pm$A=w.up
   pm$B=high.leo.and.model.prob
   z=neo(datCombined=datCombined,snpcols=snpcols,pm=pm)
   invisible(z)
}


#
# reproduce the Insig1 analysis (previously done by hand)
# with an automated version 
#

if (exists("neo.seed.insig1.test")) { rm(neo.seed.insig1.test) }
neo.seed.insig1.test=function() { 

  a=load("/home/jaten/dev/peculiar/bxh/liver.1146snps.23388mrna.21clinical.bxh.apoe.null.rdat")
  liver.bxh.male.female=liver.bxh.male.female.cm
  x=liver.bxh.male.female.cm

  female=liver.bxh.male.female$sex==2
  male=liver.bxh.male.female$sex==1
  datCombined=liver.bxh.male.female[female,-str.me$sex]
#  datCombined=liver.bxh.male.female[male,-str.me$sex]

  pm=neo.get.param()  

  which.Insig1=pmatch("Insig1",colnames(liver.bxh.male.female))
  which.Fdft1=pmatch("Fdft1.MMT00082285",colnames(liver.bxh.male.female))
  which.Dhcr7=pmatch("Dhcr7",colnames(liver.bxh.male.female))
  which.B4301=pmatch("B430110G05Rik",colnames(liver.bxh.male.female))
  
  pm$run.title="known.up.Insig1.known.dn.Fdft1"

   one.date=date() 

  my.new.dir=gsub(":",".",gsub(" ","_",paste(sep="","seed.insig1.",one.date)))
  dir.create(my.new.dir)
  setwd(my.new.dir)
  
  seed.root.leaf.neo(pm, datCombined, snpcols=str.me$snpcols, known.upstream=which.Insig1, known.downstream=which.Fdft1, genecols=str.me$genecols, require.sma.model.prob.over=.5, require.sma.leo.nb.over=.75, num.seed.snps=2, require.cor.with.known.up.over=.5) 

  
}

# MALE VERSION
if (exists("neo.seed.insig1.test.MALE")) { rm(neo.seed.insig1.test.MALE) }
neo.seed.insig1.test.MALE=function() { 

  a=load("/home/jaten/dev/peculiar/bxh/liver.1146snps.23388mrna.21clinical.bxh.apoe.null.rdat")
  liver.bxh.male.female=liver.bxh.male.female.cm
  x=liver.bxh.male.female.cm

  female=liver.bxh.male.female$sex==2
  male=liver.bxh.male.female$sex==1
#  datCombined=liver.bxh.male.female[female,-str.me$sex]
  datCombined=liver.bxh.male.female[male,-str.me$sex]

  pm=neo.get.param()  

  which.Insig1=pmatch("Insig1",colnames(liver.bxh.male.female))
  which.Fdft1=pmatch("Fdft1.MMT00082285",colnames(liver.bxh.male.female))
  which.Dhcr7=pmatch("Dhcr7",colnames(liver.bxh.male.female))
  which.B4301=pmatch("B430110G05Rik",colnames(liver.bxh.male.female))
  
  pm$run.title="known.up.Insig1.known.dn.Fdft1"

   one.date=date() 

  my.new.dir=gsub(":",".",gsub(" ","_",paste(sep="","seed.insig1.",one.date)))
  dir.create(my.new.dir)
  setwd(my.new.dir)
  
  seed.root.leaf.neo(pm, datCombined, snpcols=str.me$snpcols, known.upstream=which.Insig1, known.downstream=which.Fdft1, genecols=str.me$genecols, require.sma.model.prob.over=.5, require.sma.leo.nb.over=.75, num.seed.snps=2, require.cor.with.known.up.over=.5) 

  
}




#
# permutation version of the above, for FEMALE data. See next for MALE.
#
# started by doing, on calypso (the SGE cluster headnode):
#
###   job.list=rep("neo.seed.insig1.perm",20)
###   getwd()
###   #[1] "/home/jaten/dev/peculiar/insig1.perms"
###   source("../neo.txt")
###   submit.to.sge(job.list)
###   #your job-array 361256.1-20:1 ("neo.sge.2.qsub.me.sh") has been submitted

if (exists("neo.seed.insig1.perm")) { rm(neo.seed.insig1.perm) }
neo.seed.insig1.perm=function() { 

  a=load("/home/jaten/dev/peculiar/bxh/liver.1146snps.23388mrna.21clinical.bxh.apoe.null.rdat")
  liver.bxh.male.female=liver.bxh.male.female.cm
  x=liver.bxh.male.female.cm

  female=liver.bxh.male.female$sex==2
  datCombined=liver.bxh.male.female[female,-str.me$sex]

  pm=neo.get.param()

  which.Insig1=pmatch("Insig1",colnames(liver.bxh.male.female))
  which.Fdft1=pmatch("Fdft1.MMT00082285",colnames(liver.bxh.male.female))
  which.Dhcr7=pmatch("Dhcr7",colnames(liver.bxh.male.female))
  which.B4301=pmatch("B430110G05Rik",colnames(liver.bxh.male.female))
  
  pm$run.title="known.up.Insig1.known.dn.Fdft1"

  one.date=date()
  wd=getwd()
  sge.task = get.sge.task.or.zero()


  # do 10 permutations at once
  for (i in 1:10) {  

     my.new.dir=gsub(":",".",gsub(" ","_",paste(sep="","seed.insig1.perm.",one.date,".task.",sge.task,".perm.",i)))

     # start from where we started
     setwd(wd) 
     dir.create(my.new.dir)
     setwd(my.new.dir)

     dc=permute.snps.traits(pm, datCombined, block.of.cols.to.permute.labels.on=str.me$snpcols, save.perm.to.file=FALSE)

     seed.root.leaf.neo(pm, datCombined=dc, snpcols=str.me$snpcols, known.upstream=which.Insig1, known.downstream=which.Fdft1, genecols=str.me$genecols, require.sma.model.prob.over=.5, require.sma.leo.nb.over=.75, num.seed.snps=2, require.cor.with.known.up.over=.5) 

     

  } 
  
} # end permutation version/ neo.seed.insig1.perm()

if (exists("analyze.female.insig1.perms")) { rm(analyze.female.insig1.perms) }
analyze.female.insig1.perms=function() {

    my.dir="/home/jaten/dev/peculiar/insig1.FEMALE3.perms"
#   setwd("/home/jaten/dev/peculiar/insig1.perms")
    setwd(my.dir)
    system(paste(sep="","find . -name 'neo.logs.*.csv' -print > ",my.dir,"/dirs.with.logs.txt"))

   a=read.table(paste(sep="",my.dir,"/dirs.with.logs.txt"))
   a=as.character(a$V1)

   oca.over.th = rep(NA,length(a))
   oca.over.th.weighted = rep(NA,length(a))
   cpa = rep(NA,length(a))
   cpa.weighted = rep(NA,length(a))

   th = 0.05
    
   for (i in 1:length(a)) {
      my=read.csv(a[i])

      oca.over.th[i]=sum(my$LEO.NB.OCA > th,na.rm=T)
      oca.over.th.weighted[i]=sum(my$P.weighted.LEO.NB.OCA > th,na.rm=T)
      cpa[i]=sum(my$LEO.NB.CPA > th,na.rm=T)
      cpa.weighted[i]=sum(my$LEO.NB.CPA * my$Model.P.value.AtoB > th,na.rm=T)
      
    }

    summary(oca.over.th)
    summary(oca.over.th.weighted)
    summary(cpa)
    summary(cpa.weighted)

    # how does this compare to the actual observed female insig1 run?

    # obs.female.insig1 = read.csv("/home/jaten/dev/peculiar/insig1.FEMALE/seed.insig1.Mon_Feb_18_15.56.29_2008/neo.logs.Mon_Feb_18_16.28.21_2008.known.up.Insig1.known.dn.Fdft1/neo.logs.Mon_Feb_18_16.28.21_2008.known.up.Insig1.known.dn.Fdft1.csv")

    obs.female.insig1 = read.csv("/home/jaten/dev/peculiar/insig1.FEMALE3/seed.insig1.Wed_Feb_20_00.59.24_2008/neo.logs.Wed_Feb_20_01.24.38_2008.known.up.Insig1.known.dn.Fdft1/neo.logs.Wed_Feb_20_01.24.38_2008.known.up.Insig1.known.dn.Fdft1.csv")

     obs.oca.over.th=sum(obs.female.insig1$LEO.NB.OCA > th,na.rm=T)
     obs.oca.over.th.weighted=sum(obs.female.insig1$P.weighted.LEO.NB.OCA > th,na.rm=T)
     obs.cpa=sum(obs.female.insig1$LEO.NB.CPA > th,na.rm=T)
     obs.cpa.weighted=sum(obs.female.insig1$LEO.NB.CPA * obs.female.insig1$Model.P.value.AtoB > th,na.rm=T)

     fdr.oca =  mean(oca.over.th) / obs.oca.over.th 
     fdr.weighted.oca = mean(oca.over.th.weighted) / obs.oca.over.th.weighted    
     fdr.cpa = mean(cpa) / obs.cpa
     fdr.weighted.cpa = mean(cpa.weighted) /obs.cpa.weighted

     print(paste("At th=",th,"   we see FDR for          OCA of ",signif(fdr.oca,3),"  (",signif(mean(oca.over.th),3),"/",obs.oca.over.th,")" ))
     print(paste("At th=",th,"   we see FDR for weighted OCA of ",signif(fdr.weighted.oca,3),"  (",signif(mean(oca.over.th.weighted),3),"/",obs.oca.over.th.weighted,")" ))
     print(paste("At th=",th,"   we see FDR of           CPA of ",signif(fdr.cpa,3),"  (",signif(mean(cpa),3),"/",obs.cpa,")" ))
     print(paste("At th=",th,"   we see FDR of  weighted CPA of ",signif(fdr.weighted.cpa,3),"  (",signif(mean(cpa.weighted),3),"/",obs.cpa.weighted,")" ))
    

# at th=.2, for females, 100 permutations done.
## > fdr.oca
## [1] 0.4174603
## > fdr.weighted.oca
## [1] 0.2948649
## > fdr.cpa
## [1] 0.5705217
## > fdr.weighted.cpa
## [1] 0.4858537
## > 
    

}



###   job.list=rep("MALE.neo.seed.insig1.perm",20)
###   getwd()
###   #[1] "/home/jaten/dev/peculiar/insig1.MALE.perms"
###   source("../neo.txt")
###   submit.to.sge(job.list)
###   #your job-array 361253.1-20:1 ("neo.sge.2.qsub.me.sh") has been submitted

if (exists("MALE.neo.seed.insig1.perm")) { rm(MALE.neo.seed.insig1.perm) }
MALE.neo.seed.insig1.perm=function() { 

  a=load("/home/jaten/dev/peculiar/bxh/liver.1146snps.23388mrna.21clinical.bxh.apoe.null.rdat")
  liver.bxh.male.female=liver.bxh.male.female.cm
  x=liver.bxh.male.female.cm

  female=liver.bxh.male.female$sex==2
  male=liver.bxh.male.female$sex==1
#  datCombined=liver.bxh.male.female[female,-str.me$sex]
  datCombined=liver.bxh.male.female[male,-str.me$sex]

  pm=neo.get.param()

  which.Insig1=pmatch("Insig1",colnames(liver.bxh.male.female))
  which.Fdft1=pmatch("Fdft1.MMT00082285",colnames(liver.bxh.male.female))
  which.Dhcr7=pmatch("Dhcr7",colnames(liver.bxh.male.female))
  which.B4301=pmatch("B430110G05Rik",colnames(liver.bxh.male.female))
  
  pm$run.title="known.up.Insig1.known.dn.Fdft1"

  one.date=date()
  wd=getwd()
  sge.task = get.sge.task.or.zero()


  # do 10 permutations at once
  for (i in 1:10) {  

     my.new.dir=gsub(":",".",gsub(" ","_",paste(sep="","seed.insig1.perm.",one.date,".task.",sge.task,".perm.",i)))

     # start from where we started
     setwd(wd) 
     dir.create(my.new.dir)
     setwd(my.new.dir)

     dc=permute.snps.traits(pm, datCombined, block.of.cols.to.permute.labels.on=str.me$snpcols, save.perm.to.file=FALSE)

     seed.root.leaf.neo(pm, datCombined=dc, snpcols=str.me$snpcols, known.upstream=which.Insig1, known.downstream=which.Fdft1, genecols=str.me$genecols, require.sma.model.prob.over=.5, require.sma.leo.nb.over=.75, num.seed.snps=2, require.cor.with.known.up.over=.5) 

     

  } 
  
} # end permutation version/ MALE.neo.seed.insig1.perm()



if (exists("confirm.insig1.in.male")) { rm(confirm.insig1.in.male) }
confirm.insig1.in.male=function() { 

  # .3 weighted OCA or higher

  my.female.genes =c("Fasn.MMT00018029.chr11.bp120359985","Cyp2r1.MMT00058486.chr7.bp12048523","Cyp51.MMT00071535.chr5.bp4095740","Sc4mol.MMT00059097.chr8.bp64574528","Sc5d.MMT00032020.chr9.bp43970434","Paox.MMT00036605.chr7.bp132018721","Khk.MMT00011470.chr5.bp29258979","Slc25a1.MMT00054225.chr16.bp17530063","Wdr5.MMT00035187.chr2.bp27759428","Pfkfb1.MMT00078546.chrX.bp138292585","MGC18837.MMT00022933.chr9.bp22526985","Slc23a1.MMT00022094.chr18.bp35778092","Zdhhc6.MMT00047560.chr19.bp54943601","Acac.MMT00048501.chr11.bp83907075","Gale.MMT00077659.chr4.bp134743172","Scd1.MMT00022878.chr19.bp44497177","MMT00012966.MMT00012966.chr5.bp85212161","Qdpr.MMT00063965.chr5.bp44662614","x5730427C23Rik.MMT00054020.chr7.bp133270810","Thrsp.MMT00053283.chr7.bp89494761","AV047578.MMT00024954.chr15.bp32068789")

#  my.female.genes = c("Supt3h.MMT00063933.chr17.bp44317383","Pmvk.MMT00030469.chr3.bp91702875","Acly.MMT00022101.chr11.bp100059972","Sqle.MMT00000743.chr15.bp60651347","Fdft1.MMT00080352.chr14.bp54296665","x0610007P14Rik.MMT00054614.chr12.bp81807416","Slc25a1.MMT00054225.chr16.bp17530063","Leng1.MMT00012778.chr7.bp3727812","BC036563.MMT00012566.chr4.bp124076539","x6030440G05Rik.MMT00031386.chr6.bp98942613")

  a=load("/home/jaten/dev/peculiar/bxh/liver.1146snps.23388mrna.21clinical.bxh.apoe.null.rdat")
  liver.bxh.male.female=liver.bxh.male.female.cm
  x=liver.bxh.male.female.cm

  female=liver.bxh.male.female$sex==2
  male=liver.bxh.male.female$sex==1
#  datCombined=liver.bxh.male.female[female,-str.me$sex]
  datCombined=liver.bxh.male.female[male,-str.me$sex]

  pm=neo.get.param()

  which.Insig1=pmatch("Insig1",colnames(liver.bxh.male.female))
  which.Fdft1=pmatch("Fdft1.MMT00082285",colnames(liver.bxh.male.female))
  which.Dhcr7=pmatch("Dhcr7",colnames(liver.bxh.male.female))
  which.B4301=pmatch("B430110G05Rik",colnames(liver.bxh.male.female))
  
  pm$run.title="confirm.in.male.Insig1.top10"

  cn=colnames(datCombined)

  pm$A=which.Insig1
  pm$B=match(my.female.genes,cn)
  
  neo(datCombined,snpcols=str.me$snpcols,pm=pm)
  
  }




# run the above neo.seed.insig1.perm on SGE



if (exists("is.na.or.nan")) { rm(is.na.or.nan) }
is.na.or.nan=function(x) { is.na(x) | is.nan(x)}


# helper function to determine if a named list has a particular variable
if (exists("has.name")) { rm(has.name) }
has.name=function(name,my.list) !is.na(match(name,names(my.list)))

if(exists("is.unix")) rm(is.unix)
is.unix = function() {
  if (.Platform$OS.type == "windows") return(FALSE);
  if (has.name("pkgType",.Platform) & .Platform$pkgType == "mac.binary") return(FALSE);
  return(TRUE);
}

# open a window for plotting, on windows, mac, or unix
if (exists("cross.platform.windows")) { rm(cross.platform.windows) }
cross.platform.windows=function() {
  if (.Platform$OS.type == "windows") {
     return(windows());
  }
  if (has.name("pkgType",.Platform)) {
     if (.Platform$pkgType == "mac.binary") {
        return(quartz());
     }
  }
  return(X11());
}

# check.or.impute.missing.data()
#
# if CHECK == TRUE, then stop if any missing data.
# otherwise, impute NA values the best we can. Since this
# is hackish, we really want to insist that the user take
# care of missing data on their own first.
#
if(exists("check.or.impute.missing.data")) rm(check.or.impute.missing.data);
check.or.impute.missing.data=function(datCombined,snpcols,pm,check=TRUE) {

   if (check) {
      # check and insist upon good data in datCombined
      if (any(is.na(datCombined))) {
         stop("Missing data found in the supplied datCombined data frame. Please impute the NA and missing values first. The e1071::impute() function and the NEO::impute.nn() function may be of help here.")
      }

      # check for non-varying columns
      s=colSds(as.matrix(datCombined),na.rm=TRUE)
      if (any(s==0)) {
         stop(paste(sep="","Some datCombined columns had no variation. Please remove these columns from datCombined: ",colnames(datCombined)[s==0]))
      }

      # data looks okay, use it
      return(datCombined)
   }

   # not checking, do rough and ready imputation instead.
   datCombined.snpfilled=datCombined

   # impute missing data by predicting using the other columns...in the non-ignored columns.
   if (any(is.na(datCombined[,snpcols]))) {
       pm.print(pm,paste("imputing missing values in SNPs...",date()));
       datCombined.snpfilled=impute.nn(datCombined,snpcols,pm);
       pm.print(pm,paste("done with imputation of NA for SNPs...",date()));
   }

   # impute any remainder traits using their median
   non.snp.cols=setdiff(1:ncol(datCombined),snpcols)
   x = as.data.frame(impute(datCombined.snpfilled[,non.snp.cols],what="median"))
   datCombined.snpfilled[,non.snp.cols]=x
   datCombined.snpfilled
}

# Implementation of Steve's request to have NA turned to 0, 
# and 0 turned to small Devil's number in some spreadsheet output.
# Is vectorized.
if (exists("na.to.zero.and.zero.to.small.devil")) rm(na.to.zero.and.zero.to.small.devil);
na.to.zero.and.zero.to.small.devil=function(x) {
   w0 = (x==0)
   w.na = is.na(x)
   x[w.na] = 0
   x[w0] = 666e-8
   x
}

#### create a permuted data set, by shuffling labels on the mice between
#### snps and traits
#### uses supplied parameter block.col.cols.to.permute.labels.on, and
#### uses pm$block.to.permute2, and
####      pm$block.to.permute3 if either is non-NULL
#### to enable the user to generate additional independences between
#### different blocks of data at once, if so desired.
#### The various blocks must be disjoint.
####
####
if (exists("permute.snps.traits")) rm(permute.snps.traits);
permute.snps.traits=function(pm, datCombined, block.of.cols.to.permute.labels.on, save.perm.to.file=FALSE) {
   block = block.of.cols.to.permute.labels.on
   nr=nrow(datCombined)
   if (nr <2) stop("bad call with 1 row of data to permute.snps.traits()")

   # implement additional blocks to do additional simultaneous permutation too
   # block2.to.permute
   # block3.to.permute

   # sanity checks-- blocks should not overlap, or else within block structure
   # will not be preserved.
   #
   if (length(intersect(block.of.cols.to.permute.labels.on,pm$block.to.permute2)) != 0) {
      stop("block.to.permute and block.to.permute2 cannot overlap!")
   }
   if (length(intersect(block.of.cols.to.permute.labels.on,pm$block.to.permute3)) != 0) {
      stop("block.to.permute and block.to.permute3 cannot overlap!")
   }
   if (length(intersect(pm$block.to.permute2,pm$block.to.permute3)) != 0) {
      stop("block.to.permute2 and block.to.permute3 cannot overlap!")
   }

   # sanity check range
   if (any(!(block.of.cols.to.permute.labels.on %in% (1:ncol(datCombined))))) {
      stop("Bad specification of block.of.cols.to.permute.labels.on: out of range")
   }
   if (any(!(pm$block.to.permute2 %in% (1:ncol(datCombined))))) {
      stop("Bad specification of pm$block.to.permute2: out of range")
   }
   if (any(!(pm$block.to.permute3 %in% (1:ncol(datCombined))))) {
      stop("Bad specification of pm$block.to.permute3: out of range")
   }


   perm = sample(nr)
   if (save.perm.to.file) {
      write.csv(file=paste(sep="",pm$neo.log.file,".permute.rows.",pm$perm.count-1,".txt"),perm,row.names=FALSE)
   }

  if (!is.null(pm$block.to.perm2)) {
     perm2 = sample(nr)
     datCombined[,pm$block.to.perm2]=datCombined[perm2,pm$block.to.perm2]
  }

  if (!is.null(pm$block.to.perm3)) {
     perm3 = sample(nr)
     datCombined[,pm$block.to.perm3]=datCombined[perm3,pm$block.to.perm3]
  }

   permutedSNPs=datCombined[perm,block]
   datCombined[,block]=permutedSNPs
   datCombined
}


# apply above function to multiple columns
if (exists("zero.dev.on.cols")) rm(zero.dev.on.cols);
zero.dev.on.cols=function(x,which.columns) {
   for (i in which.columns) {
      x[,i] = na.to.zero.and.zero.to.small.devil(x[,i])
   }
   x
}


# multiplier to convert natural log to log10
log.to.log10 = log10(exp(1));

# make it easy to reload this file after any changes.
if (exists("reload")) rm(reload);
if (exists("reneo")) rm(reneo);
reload=reneo=function() { source("neo.txt"); }


# A frequent user error is to input columns that are
# completely redundant with each other. This multicollinearity
# means it is impossible to distinguish between 
# two different redundant or even highly correlated variables.

# check.multi.collinearity(): checks for redundancy and
# stops with a complaint to fix it first before continuing.
# We try to do so in O(n log n) time for n columns to check,
# rather than in O(n^2) time, by comparing correlations
# with the first column and only investigating those columns with
# very high correlation with that first column, and those
# that share a very similar correlation with the first column.
#
#
if (exists("check.multi.collinearity")) { rm(check.multi.collinearity) }
check.multi.collinearity=function(x,pm=neo.get.param()) {

  r=list() # results
  r$has.redundant = FALSE # default
  r$redundant.cols = c() # default
  if (ncol(x) < 2) return(r);

  corx = cor(x[,1],x[,-1],use="p")

  # check #1: any highly correlated with first column?
  w = which(abs(corx) >= pm$datCombined.column.redundancy.reject.cor.level) 
  if (any(w)) {
      r$has.redundant = TRUE
      r$redundant.cols = c(1,w+1)
      return(r)
  }

  # if only 2 columns then we are done
  if (ncol(x) < 3) return(r);

  # IVAR: we have 3 or more columns, so length(corx) >= 2
  if (length(corx) < 2) stop("logic error in check.multi.collinearity")

  # check #2: investigate more closely any columns who
  # share the same correlation or close correlation with
  # that first column. Designed to avoid the all-by-all
  # time consuming O(n^2) correlation matrix computation.

  # remember that our indices start with the 2nd column of x being the 1st correlation in w  

  sx = sort(corx, decreasing=T, index.return=T)
  scorx = sx$x
  scorx.index.ordered = sx$ix

  # This is a heuristic fast approximation to detect multicolinearity, not a catch all.
  #
  # Analysis:
  #
  # We ignore pairs of columns if their sign of correlation is different.
  #
  # Assuming that the sign of the correlation is the same,
  # by the triangle inequality we know that |cor(a,b) - cor(a,c)| <= |1-cor(b,c)|.
  # 
  # If we assume that cor(b,c) is positive, then
  #     cor(b,c) <= 1 - |cor(a,b) - cor(a,c)|
  #
  # Since we want cor(b,c) < pm$datCombined.column.redundancy.reject.cor.level,
  #   we only compute cor(b,c) if we first find that 1 - |cor(a,b) - cor(a,c)| >= pm$datCombined.column.redundancy.reject.cor.level.
  #   or in other words, 1 - pm$datCombined.column.redundancy.reject.cor.level >= |cor(a,b) - cor(a,c)|.
  #
  # 
  flag.if.less.than = 1 - abs(pm$datCombined.column.redundancy.reject.cor.level)

  for (i in 1:(length(scorx)-1)) {
     j = i+1
     ab = scorx[i]
     ac = scorx[j]
     if (ab*ac <= 0) next; # skip if different signs or if one cor is zero.

     # with luck, few columns will pass this check and need explicit checking...
     # But also the sort may not have brought the closest pairs adjecent to each other, which
     # is why this is a hueristic function.
     # 
     if (abs(ab-ac) <= flag.if.less.than) {
        # compute actual correlation to check it
        bc.cor = cor(x[,i+1],x[,j+1],use="p") # add one to i and j b/c we never compute the correlation of the first column with itself.
        if (abs(bc.cor) >= pm$datCombined.column.redundancy.reject.cor.level) {
           r$has.redundant = TRUE
           r$redundant.cols = c(i+1,j+1)

           msg = paste("datCombined input problem detected. Columns numbered (",
                  paste(co.lin$redundant.cols,sep=",",collapse=","),") were redundant (correlated at |",signif(bc.cor,3),"| >= pm$datCombined.column.redundancy.reject.cor.level(",signif(pm$datCombined.column.redundancy.reject.cor.level,3),")) in datCombined. These were the columns labelled: ",
                  paste(collapse=",",sep=",",colnames(datCombined)[co.lin$redundant.cols]))

           #msg2 = paste("Here is their correlation matrix: ",
           #        stringify(cor(datCombined[,co.lin$redundant.cols],use="p")))
           r$msg = paste(msg,collapse="\n",sep="\n")

           return(r)
        }
     }
  } # end for i

  r
}


#
# Since Excel can be installed in many places, try hard to
# find it if not specified correctly by pm$excel.path.
#
# Then we can automatically open the spreadsheet at the
# end of the run.
#
if (exists("locate.excel")) { rm(locate.excel) }
locate.excel=function(pm=neo.get.param()) {

   if(!is.null(pm$excel.path)) {
      if (file.exists(pm$excel.path)) {
      return(pm$excel.path)
   }}
   # INVAR: pm$excel.path did not work, now try to guess...

   #helper function for pasting compactly
   pp=function(...) paste(sep="",...)

   # try to figure the system drive
   dr=Sys.getenv("SystemDrive")
   if (dr=="") { dr="C:" } # default

   excel.path.try = c(pp(dr,"\\Program Files\\Microsoft Office\\Office\\EXCEL.EXE"),pp(dr,"\\Program Files\\Microsoft Office\\Office10\\EXCEL.EXE"),pp(dr,"\\Program Files\\Microsoft Office\\Office11\\EXCEL.EXE"),pp(dr,"\\Program Files\\Microsoft Office\\Office12\\EXCEL.EXE"),pp(dr,"\\Program Files\\Microsoft Office\\Office13\\EXCEL.EXE"),pp(dr,"\\Program Files\\Microsoft Office\\Office14\\EXCEL.EXE"))

   for (i in 1:length(excel.path.try)) {
      if (file.exists(excel.path.try[i])) {
         return(excel.path.try[i])
      }
   }

   # default
   excel.path.try[1]
} # end locate.excel

#
# open.excel.file(): only works on Windows.
#
if (exists("open.excel.file")) { rm(open.excel.file) }
open.excel.file=function(fname,pm) {
   wd=getwd()
   path=gsub("/","\\\\",wd)
   file.path=paste(sep=""," \"",path,"\\",fname,"\"")
   file.path=path.os.convert(file.path); # customize to OS; for windows change / to \\

   excel.path=locate.excel(pm)

   dome=paste(sep="","\"",excel.path,"\" ",file.path)

   if (.Platform$OS.type == "windows") {
      # debug
      print(paste("Trying:",dome))

      system(command=dome,wait=FALSE,invisible=FALSE,show.output.on.console = FALSE)
    }
}

# helper function to convert "" (factor) entries to NA
if (exists("make.dataframe.numeric")) { rm(make.dataframe.numeric) }
make.dataframe.numeric=function(df,do.cols=1:ncol(df)) {
  newdf = df
  for (j in do.cols) {
     newdf[,j] = as.numeric(as.character(df[,j]))
  }
  newdf
}

if (exists("new.open.excel.file")) { rm(new.open.excel.file) }
new.open.excel.file=function(fname,pm) {
   wd=getwd()
   path=gsub("/","\\\\",wd)
   file.path=paste(sep=""," \"",path,"\\",fname,"\"")
   file.path=path.os.convert(file.path); # customize to OS

   # try hard to get a good excel path
   excel.path = pm$excel.path

   if(is.null(excel.path)) {
      excel.path=path.os.convert("C:\\Program Files\\Microsoft Office\\Office\\EXCEL.EXE");
   }

   len.xl.path = nchar(excel.path)
   if (len.xl.path < 9) excel.path = "C:\\Program Files\\Microsoft Office\\Office11\\EXCEL.EXE"

   excel.path.no.quotes = excel.path
   if (substring(excel.path,1,1) == "\"") {
      excel.path.no.quotes = substring(excel.path,2,nchar(excel.path)-1) 
   }

   if (file.exists(excel.path)) {

   }

   dome=paste(sep="",excel.path,file.path)

   if (.Platform $OS.type == "windows") {
      system(dome,wait=FALSE)
    }
}

if (exists("write.html.log")) { rm(write.html.log) }
write.html.log=function(pm,fname="default.log.html", title.1="My Title.1",title.2="(more detail)", title.3="title3", my.text="My test content.", html.log.dir="html.logfiles",parent.dir=getwd()) {
   dir.create(html.log.dir,showWarnings=FALSE,recursive=TRUE)
   wd=getwd()
   path=gsub("/","\\\\",wd)
   file.path=paste(sep="",path,"\\",html.log.dir,"\\",fname)
   file.path=path.os.convert(file.path); # customize to OS
   
  zz <- file(file.path, "wt")
     
  a1= "<html>\n<head>\n<title>";
  a2= "</title>\n</head>\n<body><PRE>\n\n";
  a3= "\n\n</PRE></body>\n</html>\n";
  
  title.string=paste(title.2,"...",title.1,"  ",title.3);

  writeChar(a1, con=zz,eos=NULL)
  writeChar(title.string, con=zz,eos=NULL)
  writeChar(a2, con=zz,eos=NULL)
  writeChar(my.text, con=zz,eos=NULL)
  writeChar(a3, con=zz,eos=NULL)
  close(zz);

  # Use basename() to trim the html directory name so that it is relative to the output directory--so
  # we can copy/email this directory and still have the hyperlink work!
  relocatable.file.path = paste(sep="",basename(html.log.dir),"\\",basename(file.path));
   
  my.link=paste(sep="","=HYPERLINK(\"",relocatable.file.path,"\",\"",title.1,"\")")
  path.os.convert(my.link); # customize to OS
}

# take the entries in a vector and convert the NAs to blanks: makes a character vector
if(exists("convert.na.to.blank") ) rm(convert.na.to.blank);
convert.na.to.blank=function(pm,v) {
   if (pm$turn.NA.scores.to.zero.and.zero.to.666) {v=na.to.zero.and.zero.to.small.devil(v) }
   v[which(is.na(as.character(v)))]=""
   v
}

# linab(): map from a pair indices (a,b) drawn from pm$A and pm$B (restricted choices
# of network edges to be analyzed rather than analyzing the complete graph by default)
# into a single LINEAR row number for a result set data frame.
# See function enumab() below for the reverse mapping.
#
if(exists("linab") ) rm(linab);
linab = function(a,b,pm) {
   lenA = length(pm$A)
   lenB = length(pm$B)
   if (lenA < 1 | lenB < 1) stop("bad call to linab(): lenA or lenB < 1.")

   addone = 0; # increment to 1 if a comes from pm$B and b comes from pm$A.
   wa = match(a,pm$A)
   wb = match(b,pm$B)

   if (is.na(wa)) {
      if (!is.na(wb)) stop("bad call to linab(): ")
      wa = match(b,pm$A)
      wb = match(a,pm$B)
      addone = 1
   }

  if (is.na(wa) | is.na(wb) ) {
     stop(paste(sep="","Bad call to linab(): wa[",wa,"] or wb[",wb,"] NOT FOUND in pm$A[",pm$A,"] or pm$B[",pm$B,"]"))
  }

  index = 1 + ((wb + (wa-1)*lenB)-1)*2 + addone
  index
}

# enumab(): opposite mapping of linab. Here we
# map from a single row number i to a pair of indices
# (a,b) where "a" is in pm$A and "b" is in pm$B.
#
if(exists("enumab") ) rm(enumab);
enumab = function(k,pm) {

   lenA = length(pm$A)
   lenB = length(pm$B)
   if (lenA < 1 | lenB < 1) stop("bad call to enumab(): lenA or lenB < 1.")
   if (k > 2* lenA * lenB) stop("bad call to enumab(): k > 2 * lenA * lenB.")

   i = (k+1) %/% 2 # irrespective of order, which pm$A,pm$B do we want?

   i.by.lb= ((i-1) %/% lenB) +1
   i.mod.lb = ((i-1) %% lenB) + 1

   if (k %% 2 == 0) {
      # even k
      ret = c(pm$B[i.mod.lb], pm$A[i.by.lb])
   } else { # odd
      ret = c(pm$A[i.by.lb], pm$B[i.mod.lb])
   }
   ret
}

# enumab.cn(): Same function as above, but return column
# names instead of numbers to make matching more robust.
#
if(exists("enumab.cn") ) rm(enumab.cn);
enumab.cn = function(k,pm) {

   lenA = length(pm$A.cn)
   lenB = length(pm$B.cn)
   if (lenA < 1 | lenB < 1) stop("bad call to enumab.cn(): length(pm$A.cn) or length(pm$B.cn) < 1.")

   i = (k+1) %/% 2 # irrespective of order, which pm$A.cn,pm$B.cn do we want?

   i.by.lb= ((i-1) %/% lenB) +1
   i.mod.lb = ((i-1) %% lenB) + 1

   if (k %% 2 == 0) {
      # even k
      ret = c(pm$B.cn[i.mod.lb], pm$A.cn[i.by.lb])
   } else { # odd
      ret = c(pm$A.cn[i.by.lb], pm$B.cn[i.mod.lb])
   }
   ret
}


# 
# Get 2 vectors of edges to scan from pm$A and pm$B; store in pm;
#    so pm=get.enumab.set(pm) is the call.
# 
# Use: for (k9 in 1:length(pm$enumerated.A)) { # use pm$enumerate.A[k9], pm$enumerated.B[k9] }
# 
if(exists("get.enumab.set") ) rm(get.enumab.set);
get.enumab.set=function(pm) {
   lenA = length(pm$A)
   lenB = length(pm$B)
   if (lenA < 1 | lenB < 1) stop("bad call to get.enumab.set(): lenA or lenB < 1")
  
   count = lenA*lenB

   enumerated.A = rep(NA,count)
   enumerated.B = rep(NA,count)

   for (k in 1:count) {
      e = enumab(1+(k-1)*2,pm)
      enumerated.A[k] = e[1]
      enumerated.B[k] = e[2]      
   }

   pm$enumerated.A = enumerated.A
   pm$enumerated.B = enumerated.B

   pm
}


# same as above, but use actual column names instead
# of indices.
# 
# Get 2 vectors of edges to scan from pm$A and pm$B; store in pm;
#    so pm=get.enumab.set(pm) is the call.
# 
# Use: for (k9 in 1:length(pm$enumerated.A)) { # use pm$enumerate.A[k9], pm$enumerated.B[k9] }
# 
if(exists("get.enumab.cn.set") ) rm(get.enumab.cn.set);
get.enumab.cn.set=function(pm) {
   lenA = length(pm$A.cn)
   lenB = length(pm$B.cn)
   if (lenA < 1 | lenB < 1) stop("bad call to get.enumab.cn.set(): lenA or lenB < 1")
  
   count = lenA*lenB

   enumerated.A.cn = rep(NA,count)
   enumerated.B.cn = rep(NA,count)

   for (k in 1:count) {
      e = enumab.cn(1+(k-1)*2,pm)
      enumerated.A.cn[k] = e[1]
      enumerated.B.cn[k] = e[2]      
   }

   pm$enumerated.A.cn = enumerated.A.cn
   pm$enumerated.B.cn = enumerated.B.cn

   pm
}



# convert from covariance matrix to correlation matrix
#  reduced from ggm::correlations()
if(exists("covx2cor") ) rm(covx2cor);
covx2cor=function(covx) {
   Dg <- 1/sqrt(diag(covx))
   r <- covx * outer(Dg, Dg)
   r
 }

# convert from covariance matrix to partial correlation matrix
# use ggm library function
if(exists("covx2pcor") ) rm(covx2pcor);
covx2pcor=function(covx) {
      parcor(covx)
   }

if(exists("try.failed") ) rm(try.failed);
try.failed=function(x) { inherits(x, "try-error") }


#
# since the sem() calls sometimes fail, log them, with the data, to reproduce later if we wish.
# Also log warning model did not converge calls.
#
if(exists("autolog.sem") ) rm(autolog.sem);
autolog.sem=function(pm,ram,covx,N,...) {
  last.warning=NULL;

  # IMPORTANT CONVERSION to CORRELATION MATRIX from COVARIANCE MATRIX
  # We were getting singular Hessians with some covariance matrices fed to sem(),
  # so we will try just feeding the correlation matrix and see if that
  # improves matters.
  if (pm$sem.fit.correlation.instead.of.covariance) { covx=covx2cor(covx); }
 
# PL: wrap the   r=try(sem(ram,covx,N,...)) in a suppressWarnings() and silent;
 
  suppressWarnings({ r=try(sem(ram,covx,N,...), silent = TRUE) });
  if(inherits(r, "try-error") || !is.null(last.warning)) {
     if(inherits(r, "try-error")) {
           r$had.error=geterrmessage();
           if (!exists(".Global.SEM.fail.count")) {
              .Global.SEM.fail.count <<- 1;
            } else {
              .Global.SEM.fail.count <<- .Global.SEM.fail.count + 1;
            }
           logfile = paste(sep="",pm$neo.log.file,".badsem.",.Global.SEM.fail.count,".rdat");
           #print(paste("sem() call failed...saving call arguments to file: ",logfile));
     } else {
           r$had.warning=last.warning;
           if (!exists(".Global.SEM.warn.count")) {
              .Global.SEM.warn.count <<- 1;
            } else {
              .Global.SEM.warn.count <<- .Global.SEM.warn.count + 1;
            }
           logfile = paste(sep="",pm$neo.log.file,".warnsem.",.Global.SEM.warn.count,".rdat");
           #print(paste("sem() call generated warning...saving call arguments to file: ",logfile));
     }
     badsem=list();
     badsem$pm=pm;
     badsem$ram=ram;
     badsem$covx=covx;
     
     # get human readable correlation and partial correlation matrices.
     badsem$cor=covx2cor(covx);
     badsem$pcor=covx2pcor(covx);
     
     badsem$N=N;
     badsem$full.call=match.call(expand.dots = FALSE); # document any ... passed arguments.
     badsem$reproduce="q=sem(ram=badsem$ram,S=badsem$covx,N=badsem$N)";
     if (!pm$no.log) save(badsem,file=logfile);
   } # end if try-error
  r
}


if(exists("try.sem") ) rm(try.sem);
try.sem=function(pm,ram,covx,N,ana.grad=TRUE) {
  a=try(autolog.sem(pm,ram,covx,N,analytic.gradient=ana.grad));
  if (try.failed(a)) {
     print(paste("try.sem: sem() call failed with analytic.gradient=",ana.grad));
     print("trying the other way (reversing analytical.gradient flag). Model(ram):");
     print(ram);
     a=try(autolog.sem(pm,ram,covx,N,analytic.gradient=!ana.grad));
     if (try.failed(a)) {
        print("sem() call failed both ways...Setting edge orienting scores to NA for ram:");
        cat(ram);
        a=NA;
        class(a)="FailedSEM";
     } else { cat("\n"); print(paste("try.sem() suceeded with analytical.gradient=",!ana.grad)); }
  }
  a
}

#
# Set up the mechanism by which we pass back NA for our sem() failed calls, from try.sem();
#
if(exists("summary.FailedSEM") ) rm(summary.FailedSEM);
summary.FailedSEM=function(a) {
  z=list()
  z$chisq=NA
  z
}


# fit an arbitrary sem model, given a M[from,to] matrix
# where which.x.cols specifies the order of the columns
# and rows of the matrix M and picks out the columns
# of the data frame x.
#
if(exists("fit.sem.from.matrix") ) rm(fit.sem.from.matrix);
fit.sem.from.matrix=function(which.x.cols,M,x,pm=neo.get.param(),use.ranks=FALSE) {

   cn=colnames(M)
   if (length(cn) > 0) {
      if (!all(cn==colnames(x)[which.x.cols])) {
         stop("fit.sem.from.matrix(): colnames(x)[which.x.cols] must match the column labels on the matrix M")
      }
   } else {
       stop("fit.sem.from.matrix(): must have column labels on supplied matrix M to report SEM model");
   }

   no.obs=nrow(x);    # allow no.obs to override the nrow(x), if specified.
   if (is.null(pm$no.obs.Z)) {pm$no.obs.Z = no.obs} # so we get correct zeo2() scores

   subx=data.frame(x[,which.x.cols])

   if (use.ranks) {
    for (g in 1:ncol(subx)) {
         subx[,g]=rank(subx[,g])
    }}

   covx=cov(subx);

   sem.M = make.ram(M)

   fit.M = try.sem(pm,sem.M$the.ram,covx,N=no.obs)
   fit.M.summary = summary(fit.M)

   z=list()
   z$sem = fit.M
   z$sem.summary = fit.M.summary
   z
}



#
# compare.local.sems: return the SEM fitting indices for 
#  a causal, reactive, and confounded model
#
# is safe if M.col is a vector too...
#
# either supply dataframe x directly....or give both covx and no.obs
#
if(exists("compare.local.sems") ) rm(compare.local.sems);
compare.local.sems=function(M.col,A.col,B.col,x=NULL,covx=NULL,no.obs=NULL,pm) {

   z=list() # return value

   if(is.null(x) && is.null(covx)) { 
      stop("compare.local.sems() called without supply either data matrix x or covariance matrix cov. Aborting.")
   }
   if (is.null(x) && is.null(no.obs)) {
      stop("compare.local.sems() called with covariance matrix but without giving no.obs. Aborting.")
   }
 
   # necessary for AhiddenB code below to always work...
   if(length(A.col) !=1 || length(B.col) != 1) {
       stop("compare.local.sems called with illegal input: A or B had more or less than 1 variable specified.")
   }


   if (!is.null(x)) {
     if (is.null(no.obs)) {
         no.obs=nrow(x);    # allow no.obs to override the nrow(x), if specified.
     }
     cn=colnames(x);
     subx=data.frame(x[,c(M.col,A.col,B.col)])
     covx=cov(subx);

   } else {
     cn=colnames(covx);
   }

   if (is.null(pm$no.obs.Z)) {pm$no.obs.Z = no.obs} # so we get correct zeo2() scores

   dn=cn[c(M.col,A.col,B.col)];
   nr=length(dn);
   M.empty=matrix(rep(0,nr*nr),nrow=nr,dimnames=list(dn,dn));

   # adjust indices to point into new matrix
   M.col=match(cn[M.col],dn)
   A.col  =match(cn[A.col],dn)
   B.col  =match(cn[B.col],dn)

   # model 1
   M.AtoB=M.empty;
   M.AtoB[M.col,A.col]=1;
   M.AtoB[A.col,B.col]=1;   
#   attr(M.AtoB,"model")=1;

   # model 2
   M.BtoA=M.empty;
   M.BtoA[M.col,B.col]=1;
   M.BtoA[B.col,A.col]=1;
#   attr(M.BtoA,"model")=2;

   # model 3
   M.conf=M.empty;
   M.conf[M.col,A.col]=1;
   M.conf[M.col,B.col]=1;
#   attr(M.conf,"model")=3;

   # model 4
   M.AcollideB=M.empty;
   M.AcollideB[M.col,A.col]=1;
   M.AcollideB[B.col,A.col]=1;
#   attr(M.AcollideB,"model")=4;

   # model 5
   M.BcollideA=M.empty;
   M.BcollideA[M.col,B.col]=1;
   M.BcollideA[A.col,B.col]=1;
#   attr(M.BcollideA,"model")=5;

   # model 6
   M.AhiddenB=M.empty;
   M.AhiddenB[M.col,A.col]=1;
#   attr(M.AhiddenB,"model")=6;

   intra.M.list = generate.intra.ma.pairlist(M.col,c(),pm)
   
   sem.M.conf = make.ram(M.conf, intra.M.list)
   sem.M.AtoB = make.ram(M.AtoB, intra.M.list)
   sem.M.BtoA = make.ram(M.BtoA, intra.M.list)
   sem.M.BcollideA = make.ram(M.BcollideA, intra.M.list)
   sem.M.AcollideB = make.ram(M.AcollideB, intra.M.list)

   intra.M.list[[ length(intra.M.list)+1]] = c(A.col,B.col)
   sem.M.AhiddenB = make.ram(M.AhiddenB, intra.M.list)

   fit.M.conf = summary(try.sem(pm,sem.M.conf$the.ram,covx,N=no.obs))
   fit.M.AtoB = summary(try.sem(pm,sem.M.AtoB$the.ram,covx,N=no.obs))
   fit.M.BtoA = summary(try.sem(pm,sem.M.BtoA$the.ram,covx,N=no.obs))
   fit.M.BcollideA = summary(try.sem(pm,sem.M.BcollideA$the.ram,covx,N=no.obs))
   fit.M.AcollideB = summary(try.sem(pm,sem.M.AcollideB$the.ram,covx,N=no.obs))

   # fit.M.AhiddenB always takes like 37 iterations and returns the same
   # value as AcollideB anyway, i.e. the models are not identifiable.
   #   fit.M.AhiddenB = summary(try.sem(pm,sem.M.AhiddenB$the.ram,covx,N=no.obs))

   


   z$title=paste(sep="","M(",paste(dn[M.col],collapse=","),"),A(",paste(dn[A.col],collapse=","),"),B(",paste(dn[B.col],collapse=","),") local SEM model comparisons to M->A->B.")

   if (!is.null(x)) {
     # we want to know how much of A the markers explain, and how much of B
     form1=paste(dn[A.col],"~",paste(collapse=" + ",dn[M.col]))
     summy1=summary(lm(as.formula(form1),data=subx))
     z$A.predicted.by.markers.R.squared = summy1$r.squared
     z$A.predicted.by.markers.Adjusted.R.squared = summy1$adj.r.squared
     z$A.predicted.by.markers.formula = form1

     form2 = paste(dn[B.col],"~",paste(collapse=" + ",dn[M.col]))
     summy2=summary(lm(as.formula(form2),data=subx))
     z$B.predicted.by.markers.R.squared = summy2$r.squared
     z$B.predicted.by.markers.Adjusted.R.squared = summy2$adj.r.square
     z$B.predicted.by.markers.formula = form2

     form3 = paste(dn[B.col],"~",dn[A.col])
     summy3=summary(lm(as.formula(form3),data=subx))
     z$B.predicted.by.A.gives.R.squared = summy3$r.squared
     z$B.predicted.by.A.gives.Adjusted.R.squared = summy3$adj.r.squared
     z$B.predicted.by.A.formula = form3
   }

   z$M.AtoB=fit.M.AtoB
   z$M.BtoA=fit.M.BtoA
   z$M.conf=fit.M.conf
   z$M.AcollideB=fit.M.AcollideB
   z$M.BcollideA=fit.M.BcollideA

#   z$M.AhiddenB=fit.M.AhiddenB

   # for confounding checking...
   if (length(M.col) ==1) {
      # single marker version
      z$zeo.M.A.given.B=zeo2detail(M.col,A.col,B.col,covx,pm) # $BLV is the final score
      z$zeo.M.B.given.A=zeo2detail(M.col,B.col,A.col,covx,pm) # $BLV
   } else {
      # average the BLVs for each marker
      z$zeo.M.A.given.B=multimarker.zeo2detail(M.col,A.col,B.col,covx,pm) # $BLV is the final score
      z$zeo.M.B.given.A=multimarker.zeo2detail(M.col,B.col,A.col,covx,pm) # $BLV
   }

   # convert to log10
   z$mlogp.M.AtoB=-pchisq(fit.M.AtoB$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10
   z$mlogp.M.BtoA=-pchisq(fit.M.BtoA$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10
   z$mlogp.M.conf=-pchisq(fit.M.conf$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10
   
   z$mlogp.M.AcollideB=-pchisq(fit.M.AcollideB$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10
   z$mlogp.M.BcollideA=-pchisq(fit.M.BcollideA$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10

#   z$mlogp.M.AhiddenB=-pchisq(fit.M.AhiddenB$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10

   z$mlogp.in.model.order=c(z$mlogp.M.AtoB,z$mlogp.M.BtoA,z$mlogp.M.conf,z$mlogp.M.AcollideB,z$mlogp.M.BcollideA) #,z$mlogp.M.AhiddenB)
   so=sort(z$mlogp.in.model.order,decreasing=FALSE,index.return=TRUE)
   z$ranked.models=so$ix
   z$ranked.models.mlogp=so$x

   # so far AcollideB and AhiddenB look the same???
   # For now, then, we EXCLUDE mlogp.M.AhiddenB from the min() below.


   ##
   ## the new statistics ZPathAB, ZPathBA, z$LEO.NB.BtoA,
   ## BLV, PearsonCor, PearsonCorP
   ##


   z$PathAB = fit.M.AtoB$coeff[4,1]
   z$PathBA = fit.M.BtoA$coeff[3,1]

   z$SEPathAB = fit.M.AtoB$coeff[4,2]
   z$SEPathBA = fit.M.BtoA$coeff[3,2]

   z$ZPathAB = fit.M.AtoB$coeff[4,3]
   z$ZPathBA = fit.M.BtoA$coeff[3,3]

   z$PPathAB = fit.M.AtoB$coeff[4,4]
   z$PPathBA = fit.M.BtoA$coeff[3,4]

   # Prior to 4 Feb 2008, the right-hand-sides of these two were swapped.
   # We now switch them to create consistency with SH implementation of BLV.
   z$BLV.AtoB = z$zeo.M.B.given.A$BLV
   z$BLV.BtoA = z$zeo.M.A.given.B$BLV

   # for M->A->B
   z$alt.model.best.min.mlogp = min(z$mlogp.M.conf,z$mlogp.M.BtoA,z$mlogp.M.BcollideA,z$mlogp.M.AcollideB)

   # for M->A<-B
   z$next.best.min.rev.mlogp = min(z$mlogp.M.conf,z$mlogp.M.BtoA,z$mlogp.M.AtoB,z$mlogp.M.BcollideA)

   # for M->B->A
   z$next.best.min.M.BtoA.numerator = min(z$mlogp.M.conf, z$mlogp.M.AtoB, z$mlogp.M.AcollideB, z$mlogp.M.BcollideA)

   z$LEO.NB.AtoB = -(z$mlogp.M.AtoB - z$alt.model.best.min.mlogp)

   z$LEO.NB.AcollideB = -(z$mlogp.M.AcollideB - z$next.best.min.rev.mlogp)

   z$LEO.NB.BtoA = -(z$mlogp.M.BtoA - z$next.best.min.M.BtoA.numerator)
   
   z$main.model.A.to.B.mlogp=z$mlogp.M.AtoB
   z$p.conf.over.p.AtoB = pchisq(fit.M.conf$chisq,1,lower.tail=FALSE)/pchisq(fit.M.AtoB$chisq,1,lower.tail=FALSE)
   z$p.conf.over.p.BtoA = pchisq(fit.M.conf$chisq,1,lower.tail=FALSE)/pchisq(fit.M.BtoA$chisq,1,lower.tail=FALSE)

    # the -log10 space version of p.conf.over.p.AtoB becomes LEO.I.AtoB (I for independent (confounded) model comparison)
    z$LEO.I.AtoB = z$mlogp.M.conf - z$mlogp.M.AtoB
    z$LEO.I.BtoA = z$mlogp.M.conf - z$mlogp.M.BtoA

    z$LEO.O.AtoB = z$mlogp.M.AcollideB - z$mlogp.M.AtoB
    z$LEO.O.BtoA = z$mlogp.M.BcollideA - z$mlogp.M.BtoA

   
   # for symmetry of A->B vs B->A to work...we should compare to M.AcollideB
   z$eo.losem.lod = -(z$mlogp.M.AtoB - z$mlogp.M.AcollideB);

   z
}


if(exists("compare.local.sems14") ) rm(compare.local.sems14);
compare.local.sems14=function(M.col,A.col,B.col,x=NULL,covx=NULL,no.obs=NULL) {

   if(is.null(x) && is.null(covx)) { 
      stop("compare.local.sems() called without supply either data matrix x or covariance matrix cov. Aborting.")
   }
   if (is.null(x) && is.null(no.obs)) {
      stop("compare.local.sems() called with covariance matrix but without giving no.obs. Aborting.")
   }
 
   # necessary for AhiddenB code below to always work...
   if(length(A.col) !=1 || length(B.col) != 1) {
       stop("compare.local.sems called with illegal input: A or B had more or less than 1 variable specified.")
   }

   if (!is.null(x)) {
     if (is.null(no.obs)) { no.obs=nrow(x); } # allow no.obs to override nrow(x), if specified.
     cn=colnames(x);
     subx=x[,c(M.col,A.col,B.col)];
     covx=cov(subx);
   } else {
     cn=colnames(covx);
   }

   dn=cn[c(M.col,A.col,B.col)];
   nr=length(dn);
   M.empty=matrix(rep(0,nr*nr),nrow=nr,dimnames=list(dn,dn));

   # adjust indices to point into new matrix
   M.col=match(cn[M.col],dn)
   A.col  =match(cn[A.col],dn)
   B.col  =match(cn[B.col],dn)

   # model 1
   M.AtoB=M.empty;
   M.AtoB[M.col,A.col]=1;
   M.AtoB[A.col,B.col]=1;   
#   attr(M.AtoB,"model")=1;


   # model 4
   M.AcollideB=M.empty;
   M.AcollideB[M.col,A.col]=1;
   M.AcollideB[B.col,A.col]=1;
#   attr(M.AcollideB,"model")=4;


   sem.M.AtoB = make.ram(M.AtoB)
   sem.M.AcollideB = make.ram(M.AcollideB)

   fit.M.AtoB = summary(try.sem(pm,sem.M.AtoB$the.ram,covx,N=no.obs))
   fit.M.AcollideB = summary(try.sem(pm,sem.M.AcollideB$the.ram,covx,N=no.obs))

   z=list()

   z$title=paste(sep="","M(",paste(dn[M.col],collapse=","),"),A(",paste(dn[A.col],collapse=","),"),B(",paste(dn[B.col],collapse=","),") local SEM model comparisons to M->A->B.")
   z$M.AtoB=fit.M.AtoB
   z$M.AcollideB=fit.M.AcollideB

   z$mlogp.M.AtoB=-pchisq(fit.M.AtoB$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10
   z$mlogp.M.AcollideB=-pchisq(fit.M.AcollideB$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10

   # for symmetry of A->B vs B->A to work...we should compare to M.AcollideB
   z$eo.losem.lod = -(z$mlogp.M.AtoB - z$mlogp.M.AcollideB);

   z
}


#
# minimal.sem.compare: return the LOD score for the
#   A.to.B versus AcollideB models as the edge score.
# 
# This is a version of compare.local.sem that has been
# stripped down for speed.
#
# Function is safe if M.col is a vector too...
#
# Either supply x directly....or give both covx and no.obs
#
if(exists("minimal.sem.compare") ) rm(minimal.sem.compare);
minimal.sem.compare=function(M.col,A.col,B.col,x=NULL,covx=NULL,no.obs=NULL) {

   if(is.null(x) && is.null(covx)) { 
      stop("minimal.sem.compare() called without supply either data matrix x or covariance matrix cov. Aborting.")
   }
   if (is.null(x) && is.null(no.obs)) {
      stop("minimal.sem.compare() called with covariance matrix but without giving no.obs. Aborting.")
   }
 
   if (!is.null(x)) {
     if (is.null(no.obs)) {
         no.obs=nrow(x);    # allow no.obs to override the nrow(x), if specified.
     }
     cn=colnames(x);
     subx=x[,c(M.col,A.col,B.col)];
     covx=cov(subx);
   } else {
     cn=colnames(covx);
   }

   dn=cn[c(M.col,A.col,B.col)];
   nr=length(dn);
   M.empty=matrix(rep(0,nr*nr),nrow=nr,dimnames=list(dn,dn));

   # adjust indices to point into new matrix
   M.col=match(cn[M.col],dn)
   A.col  =match(cn[A.col],dn)
   B.col  =match(cn[B.col],dn)

   # model 1
   M.AtoB=M.empty;
   M.AtoB[M.col,A.col]=1;
   M.AtoB[A.col,B.col]=1;   
#   attr(M.AtoB,"model")=1;

   # model 4
   M.AcollideB=M.empty;
   M.AcollideB[M.col,A.col]=1;
   M.AcollideB[B.col,A.col]=1;
#   attr(M.AcollideB,"model")=4;

   sem.M.AtoB = make.ram(M.AtoB)
   sem.M.AcollideB = make.ram(M.AcollideB)

   fit.M.AtoB = summary(try.sem(pm,sem.M.AtoB$the.ram,covx,N=no.obs))
   fit.M.AcollideB = summary(try.sem(pm,sem.M.AcollideB$the.ram,covx,N=no.obs))

   z=list()

   z$title=paste(sep="","M(",paste(dn[M.col],collapse=","),"),A(",paste(dn[A.col],collapse=","),"),B(",paste(dn[B.col],collapse=","),") local SEM model comparisons to M->A->B.")

   z$M.AtoB=fit.M.AtoB
   z$M.AcollideB=fit.M.AcollideB

   z$mlogp.M.AtoB=-pchisq(fit.M.AtoB$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10

   z$mlogp.M.AcollideB=-pchisq(fit.M.AcollideB$chisq,1,lower.tail=FALSE,log.p=TRUE)*log.to.log10

   # for symmetry of A->B vs B->A to work...we should compare to M.AcollideB
   z$eo.losem.lod = -(z$mlogp.M.AtoB - z$mlogp.M.AcollideB);

   z
} # end minimal.sem.compare()


if(exists("df.to.latex.table") ) rm(df.to.latex.table);
df.to.latex.table=function(df,sig.fig=2,use.row.names=TRUE,caption="Caption here.") {
 nr=nrow(df);
 nc=ncol(df);
 rn=rownames(df);
 if (use.row.names) { rn=paste(rn,"& "); extracol=" & "; c2="c|"; } else { rn[]=""; extracol=""; c2=""; }

 cn=colnames(df);
 gap="c|";
 x1 = paste(sep="","\\begin{table}\n \\centering\n\\begin{tabular}{|",c2,paste(rep(gap,nc),collapse=""),"}\n\\hline\n");
 x2 = paste(sep="","     ",extracol,paste(colnames(df),collapse=" & ")," \\\\\n\\hline\n");
 x=vector("list", length(nr))
 for (i in 1:nr) {
   x[[i]]=paste(sep="",rn[i],paste(signif(df[i,],digits=sig.fig),collapse=" & ")," \\\\\n");
 }

 x3=paste(sep="","\\hline\n\\end{tabular}\n  \\caption{",caption,"}\n  \\label{Ta:df.to.latex.table}\n\\end{table}\n");
 y=paste(x1,x2,paste(x,collapse=""),x3,collapse="")
 y
}


# verstion that puts x (sd) after each entry
if(exists("df.to.latex.table.sd") ) rm(df.to.latex.table.sd);
df.to.latex.table.sd=function(df,sig.fig=2,use.row.names=TRUE,sd.table,caption="Caption here.") {
 df=signif(df,sig.fig);
 sd.table=signif(sd.table,sig.fig)
 new.tab=df;
 nr=nrow(df);
 nc=ncol(df);
 for (i in 1:nr) { for (j in 1:nc) {
   new.tab[i,j] = paste(sep="",as.character(df[i,j])," (",sd.table[i,j],")");
 }}

 rn=rownames(df);
 if (use.row.names) { rn=paste(rn,"& "); extracol=" & "; c2="c|"; } else { rn[]=""; extracol=""; c2=""; }

 cn=colnames(df);
 gap="c|";
 x1 = paste(sep="","\\begin{table}\n \\centering\n\\begin{tabular}{|",c2,paste(rep(gap,nc),collapse=""),"}\n\\hline\n");
 x2 = paste(sep="","     ",extracol,paste(colnames(df),collapse=" & ")," \\\\\n\\hline\n");
 x=vector("list", length(nr))
 for (i in 1:nr) {
   x[[i]]=paste(sep="",rn[i],paste(new.tab[i,],collapse=" & ")," \\\\\n");
 }

 x3=paste(sep="","\\hline\n\\end{tabular}\n  \\caption{",caption,"}\n  \\label{Ta:df.to.latex.table}\n\\end{table}\n");
 y=paste(x1,x2,paste(x,collapse=""),x3,collapse="")
 y
}


# version that puts x (sd) after each entry, PLUS a 2 stand dev symbol [*]
if(exists("df.to.latex.table.sd.plus.test") ) rm(df.to.latex.table.sd.plus.test);
df.to.latex.table.sd.plus.test=function(df,sig.fig=2,use.row.names=TRUE,sd.table,caption="Caption here.",sd.away.zero=2,start.table.number=NA) {
 if(is.na(start.table.number) && !exists(".Global.table.num.df.to.latex")) {
    .Global.table.num.df.to.latex <<- 1;
 } else  if(is.na(start.table.number) && exists(".Global.table.num.df.to.latex")) {
    .Global.table.num.df.to.latex <<- .Global.table.num.df.to.latex + 1;
 } else {
    .Global.table.num.df.to.latex <<- start.table.number;
 }


 pos.mark="[$\\heartsuit$]";
 neg.mark="[$\\spadesuit$]";

 df.accur = df;
 sd.table.accur = sd.table;

 df=signif(df,sig.fig);
 sd.table=signif(sd.table,sig.fig)
 new.tab=df;
 nr=nrow(df);
 nc=ncol(df);

 for (i in 1:nr) { for (j in 1:nc) {
   if (df.accur[i,j] - sd.away.zero*sd.table.accur[i,j]>0) { mark=pos.mark; } else { mark=neg.mark; }
   new.tab[i,j] = paste(sep="",as.character(df[i,j])," (",sd.table[i,j],") ",mark);
 }}

 rn=rownames(df);
 if (use.row.names) { rn=paste(rn,"& "); extracol=" & "; c2="r|"; } else { rn[]=""; extracol=""; c2=""; }

 cn=colnames(df);
 gap="r|";
# x1 = paste(sep="","\\begin{table}\n \\centering\n\\begin{tabular}{|",c2,paste(rep(gap,nc),collapse=""),"}\n\\hline\n");

 x1=paste(sep="","\\begin{table}\n\\centering\n\\begin{small}\n\\begin{tabular}{|c||",paste(rep(gap,nc),collapse=""),"}\n\\hline\n\\small effect size  & \\multicolumn{",nc,"}{c|}{ ($1/\\sigma^2_{\\varepsilon}$) The ratio of genetic to environmental variance.} \\\\\n\\cline{2-6}\n($\\gamma$) & "); ## & 0.1 & 0.2 & 0.5 & 1 & 2 \\\\\n\\hline");

 # x2 = paste(sep="","     ",extracol,paste(colnames(df),collapse=" & ")," \\\\\n\\hline\n");
 x2 = paste(sep="",paste(colnames(df),collapse=" & ")," \\\\\n\\hline\n");
 
 x=vector("list", length(nr))
 for (i in 1:nr) {
   x[[i]]=paste(sep="",rn[i],paste(new.tab[i,],collapse=" & ")," \\\\\n");
 }

 x3=paste(sep="","\\hline\n\\end{tabular}\n\\end{small}\n  \\caption{",caption," ",pos.mark," means $> ",sd.away.zero,"$ standard deviations above zero, else ",neg.mark,"}\n  \\label{Ta:df.to.latex.table",.Global.table.num.df.to.latex,"}\n\\end{table}\n");
 y=paste(x1,x2,paste(x,collapse=""),x3,collapse="")
 y
}



# version that puts x (sd) after each entry, PLUS a 2 stand dev symbol [*]
if(exists("df.to.latex.typeIerror.print") ) rm(df.to.latex.typeIerror.print);
df.to.latex.typeIerror.print=function(df,sig.fig=2,use.row.names=TRUE,sd.table,caption="Caption here.",sd.away.zero=2,start.table.number=NA) {
 if(is.na(start.table.number) && !exists(".Global.table.num.df.to.latex")) {
    .Global.table.num.df.to.latex <<- 1;
 } else  if(is.na(start.table.number) && exists(".Global.table.num.df.to.latex")) {
    .Global.table.num.df.to.latex <<- .Global.table.num.df.to.latex + 1;
 } else {
    .Global.table.num.df.to.latex <<- start.table.number;
 }


 pos.mark="[$\\heartsuit$]";
 neg.mark="[$\\spadesuit$]";

 df.accur = df;
 sd.table.accur = sd.table;

 df=signif(df,sig.fig);
 sd.table=signif(sd.table,sig.fig)
 new.tab=df;
 nr=nrow(df);
 nc=ncol(df);

 row.head1="M-A coef"
 row.head2="$(\\alpha)$"
 col.head="M-B effect coefficient $(\\beta)$"

 for (i in 1:nr) { for (j in 1:nc) {
   if (df.accur[i,j] + 2*sd.table.accur[i,j]  >2) { mark=neg.mark; } else { mark=pos.mark; }
   new.tab[i,j] = paste(sep="",as.character(df[i,j])," (",sd.table[i,j],") ",mark);
 }}

 rn=rownames(df);
 if (use.row.names) { rn=paste(rn,"& "); extracol=" & "; c2="r|"; } else { rn[]=""; extracol=""; c2=""; }

 cn=colnames(df);
 gap="r|";
# x1 = paste(sep="","\\begin{table}\n \\centering\n\\begin{tabular}{|",c2,paste(rep(gap,nc),collapse=""),"}\n\\hline\n");

 x1=paste(sep="","\\begin{table}\n\\centering\n\\begin{small}\n\\begin{tabular}{|c||",paste(rep(gap,nc),collapse=""),"}\n\\hline\n\\small ",row.head1,"  & \\multicolumn{",nc,"}{c|}{ ",col.head,"} \\\\\n\\cline{2-6}\n",row.head2," & "); ## & 0.1 & 0.2 & 0.5 & 1 & 2 \\\\\n\\hline");

 # x2 = paste(sep="","     ",extracol,paste(colnames(df),collapse=" & ")," \\\\\n\\hline\n");
 x2 = paste(sep="",paste(colnames(df),collapse=" & ")," \\\\\n\\hline\n");
 
 x=vector("list", length(nr))
 for (i in 1:nr) {
   x[[i]]=paste(sep="",rn[i],paste(new.tab[i,],collapse=" & ")," \\\\\n");
 }

 x3=paste(sep="","\\hline\n\\end{tabular}\n\\end{small}\n  \\caption{",caption," ",neg.mark," means score $+ 2\\mbox{sd} > ",sd.away.zero,"$, else ",pos.mark,"}\n  \\label{Ta:df.to.latex.table",.Global.table.num.df.to.latex,"}\n\\end{table}\n");
 y=paste(x1,x2,paste(x,collapse=""),x3,collapse="")
 y
}


# 
if(exists("power.study") ) rm(power.study);
power.study=function() {
  start=-.8
  step=.5
  fin =.8
  print(paste("starting power study for step",step,date()))
  alpha = seq(start,fin,step); asize=length(alpha);
  beta  = seq(start-.1,fin,step); bsize=length(beta);
  delta = seq(start,fin,step); dsize=length(delta);
#  no.obs = c(20,50,100,200,500,1000,2000); nsize=length(no.obs);
  N.obs=c(50,100); nsize=length(N.obs)

  a=array(data = NA, dim = c(nsize,asize,bsize,dsize), dimnames = list(as.character(N.obs),as.character(alpha),as.character(beta),as.character(delta)))

  z=list();
  z$lod=a;
  z$mlogp.M.conf=a;
  z$mlogp.M.AtoB=a;
  z$mlogp.M.BtoA=a;
  z$mlogp.M.BcollideA=a;
  z$mlogp.M.AcollideB=a;
#  z$mlogp.M.AhiddenB=a;
  z$alt.model.best.min.mlogp=a;

  for (n in nsize) {
     for (i in 1:asize) {
        for (j in 1:bsize) {
           for (k in 1:dsize) {
              my.cor=matrix(nrow=3,c(1,alpha[i],beta[j],alpha[i],1,delta[k],beta[j],delta[k],1),dimnames=list(c("M","A","B"),c("M","A","B")))

              b=compare.local.sems(pm=pm,1,2,3,x=NULL,covx=my.cor,no.obs=N.obs[n]);
            
              z$mlogp.M.conf[n,i,j,k]=b$mlogp.M.conf;
              z$mlogp.M.AtoB[n,i,j,k]=b$mlogp.M.AtoB;
              z$mlogp.M.BtoA[n,i,j,k]=b$mlogp.M.BtoA;
              z$mlogp.M.BcollideA[n,i,j,k]=b$mlogp.M.BcollideA;
              z$mlogp.M.AcollideB[n,i,j,k]=b$mlogp.M.AcollideB;
#              z$mlogp.M.AhiddenB[n,i,j,k]=b$mlogp.M.AhiddenB;
              z$alt.model.best.min.mlogp[n,i,j,k]=b$alt.model.best.min.mlogp;
  }}}}

  print(paste("Finishing power study for step",step,date()));
  z  
}
#q=power.study()

sim.model.num=function(model.num,env.var,no.samples,MAeff,ABeff,MBeff) {
   s2=sqrt(2)
   M=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1

   if (model.num < 1 || model.num > 6) stop("model.num out of range in sim.model.num()")

   if(model.num==1) {
      A=MAeff*M+rnorm(no.samples,sd=sqrt(env.var))
      B=ABeff*A+rnorm(no.samples,sd=sqrt(env.var))
   }

   if(model.num==2) {
      B=MBeff*M+rnorm(no.samples,sd=sqrt(env.var))
      A=ABeff*B+rnorm(no.samples,sd=sqrt(env.var))
   }

   if(model.num==3) {
      B=MBeff*M+rnorm(no.samples,sd=sqrt(env.var))
      A=MAeff*M+rnorm(no.samples,sd=sqrt(env.var))
   }

   if(model.num==4) {
      B=rnorm(no.samples,sd=sqrt(env.var))
      A=MAeff*M+ABeff*B+rnorm(no.samples,sd=sqrt(env.var))
   }

   if(model.num==5) {
      A=rnorm(no.samples,sd=sqrt(env.var))
      B=MBeff*M+ABeff*A+rnorm(no.samples,sd=sqrt(env.var))
   }

   if(model.num==6) {
      H=rnorm(no.samples,sd=sqrt(env.var))
      A=MAeff*M+H+rnorm(no.samples,sd=sqrt(env.var))
      B=H+rnorm(no.samples,sd=sqrt(env.var))
   }

   data.frame(M,A,B)
}

if(exists("power.study2") ) rm(power.study2);
power.study2=function() {

   s2=sqrt(2)

   # simulate
   no.samples=400
   set.seed(1)
   gen.var=1;
   env.var=1;

   d1=date()

   g.to.e = c(.1,.2,.5,1,2);
   effect.coef=c(2,1,.5,.25,.1);
   n.choice=c(50,100,200,500,1000,2000,5000);

   lec=length(effect.coef)
   lge=length(g.to.e)
   lnc=length(n.choice)

   res=c();

   lod.mean=array(data = 0, dim = c(lec,lge,lnc), dimnames = list(as.character(effect.coef),as.character(g.to.e),as.character(n.choice)));
   lod.sd=lod.mean
   eo.mean=lod.mean
   eo.sd =lod.mean   

   # different true models simulated...
   models.to.do=1
   for(i in 1:models.to.do) {

   for (n in 1:lnc) {
     no.samples=n.choice[n]

   for (g in 1:lge) {
    env.var = 1/g.to.e[g]

   for (e in 1:lec) {
    ABeff=effect.coef[e]
    MAeff=1; #effect.coef[e] #, otherwise get a false reading on geneic:env variance
    MBeff=1; #effect.coef[e]

   num.sim=1000
   model.num=1
#   disc.matrix=matrix(nrow=6,ncol=6,data=0);
   eo.m = vector(length=num.sim);
   eo.m[]=0;
   lod.m=eo.m

      for (k in 1:num.sim) {

         mod=sim.model.num(model.num=i,env.var,no.samples,MAeff,ABeff,MBeff);

         covx=cov(mod)

         b=compare.local.sems14(1,2,3,x=NULL,covx=covx,no.obs=no.samples);

         best.mod=b$ranked.models[1]

         #disc.matrix[i,best.mod] = disc.matrix[i,b$ranked.models[1]] +1

#         if (best.mod == i) { alt.mod = b$ranked.models[2] } else { alt.mod = b$ranked.models[1] }
#         lod= -(b$ranked.models.mlogp[i] - b$ranked.models.mlogp[alt.mod])
#         lod.m[k]=lod

         eo.m[k]=b$eo.losem.lod
     } # k

#        lod.mean[e,g,n]=mean(lod.m)
#        lod.sd[e,g,n]=sd(lod.m)
        eo.mean[e,g,n]=mean(eo.m)
        eo.sd[e,g,n]=sd(eo.m)

        cat(paste(e,g,n,date(),"\n"))
     } # e
     } # g
     } # n

   d2=date()
  z=list()
#  z$lod.mean=lod.mean
#  z$lod.sd=lod.sd
  z$eo.mean=eo.mean
  z$eo.sd=eo.sd

  z$d1=d1
  z$d2=d2

#  z$lod.m=lod.m
  z$eo.m=eo.m
  save(z,file=paste(sep="","power.study2.i.is.",i,".reps.",num.sim,".rdat"))    
  res=c(z,res);

  } # end i over 6 models

  res
}


# do the ZEO stuff to compare
if(exists("power.study3") ) rm(power.study3);
power.study3=function() {

   s2=sqrt(2)

   # simulate
   no.samples=400
   set.seed(1)
   gen.var=1;
   env.var=1;

   d1=date()

   g.to.e = c(.1,.2,.5,1,2);
   effect.coef=c(2,1,.5,.25,.1);
   n.choice=c(50,100,200,500,1000,2000,5000);

   lec=length(effect.coef);
   lge=length(g.to.e);
   lnc=length(n.choice);

   res=c();

   lod.mean=array(data = 0, dim = c(lec,lge,lnc), dimnames = list(as.character(effect.coef),as.character(g.to.e),as.character(n.choice)));
   lod.sd=lod.mean
   eo.mean=lod.mean
   eo.sd =lod.mean   

   # different true models simulated...
   models.to.do=1
   for(i in 1:models.to.do) {

   for (n in 1:lnc) {
     no.samples=n.choice[n]

   for (g in 1:lge) {
    env.var = 1/g.to.e[g]

   for (e in 1:lec) {
    ABeff=effect.coef[e]
    MAeff=1; #effect.coef[e] #, otherwise get a false reading on geneic:env variance
    MBeff=1; #effect.coef[e]

   num.sim=10000
   model.num=1
#   disc.matrix=matrix(nrow=6,ncol=6,data=0);
   eo.m = vector(length=num.sim);
   eo.m[]=0;
      no.obs.Z=no.samples
      for (k in 1:num.sim) {

         mod=sim.model.num(model.num=i,env.var,no.samples,MAeff,ABeff,MBeff);

         covx=cov(mod)

            # twofer.eo and n.twofer.eo hold our results...
            z=zeo2(1,3,2,covx,pm); # z$eo

            #cor1=pcor(c(1,3,c()),covx)
            #Zm1B= sqrt(no.obs.Z-pm$fisher.dof.cor)*   fisher1(abs(cor1))

            #cor2=pcor(c(1,3,2),covx)
            #Zm1BgivenA= sqrt(no.obs.Z-pm$fisher.dof.pcor1)*   fisher1(abs(cor2))

            #s=1;
            #eo2 = Zm1B - s*Zm1BgivenA;
            eo2 = z$eo;
         eo.m[k]=eo2;
     } # k

        eo.mean[e,g,n]=mean(eo.m)
        eo.sd[e,g,n]=sd(eo.m)

        cat(paste(e,g,n,date(),"\n"))
     } # e
     } # g
     } # n

   d2=date()
  z=list()

  z$eo.mean=eo.mean
  z$eo.sd=eo.sd

  z$d1=d1
  z$d2=d2

  z$eo.m=eo.m
  save(z,file=paste(sep="","power.study3zeo.i.is.",i,".reps.",num.sim,".rdat"))    
  res=c(z,res);

  } # end i over 6 models

  res
}

# SimCR(): function used in false.positive.study1 below, to simulate
# causal and reactive PC variables around the trait T
#
# A-> y -> B where y is the trait.
# Which method is better at both detecting the causal arrow and the reactive arrow?
#
if(exists("SimCR") ) rm(SimCR);
SimCR=function(no.samples=100,use.traitsnp=T, nAsnp=3,nBsnp=3, nTsnp=3, env.var=.5,randomize.snps=FALSE){
   # for each SNP we assume Mendel's laws
   set.seed(1)
   s2=sqrt(2);

   SNPA=SNPB=SNPT=c();

   for (i in 1:nAsnp) {
       SNPA=cbind(SNPA,s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T )); colnames(SNPA)[i]=paste(sep="","SNP.A",i);}

   for (i in 1:nBsnp) {
       SNPB=cbind(SNPB,s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T )); colnames(SNPB)[i]=paste(sep="","SNP.B",i);}

   for (i in 1:nTsnp) {
       SNPT=cbind(SNPT,s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T )); colnames(SNPT)[i]=paste(sep="","SNP.T",i);}


             rSNPT=SNPT;
   if (use.traitsnp==FALSE) { rSNPT[,]=0;}

   A=rnorm(no.samples,sd=env.var)
   for (i in 1:nAsnp) { A = A + SNPA[,i] }
   A=scale(A);

   Trait=A+rnorm(no.samples,sd=env.var)
   for (i in 1:nTsnp) { Trait=Trait+rSNPT[,i]; }
   Trait=scale(Trait);

   B=Trait+rnorm(no.samples,sd=env.var)
   for (i in 1:nBsnp) { B=B+SNPB[,i] }
   B=scale(B); 

   #DANGEROUS CODE - scramble the relationship between SNPs and everything
   if (randomize.snps) {
      SNPA=sample(SNPA)
      SNPB=sample(SNPB)
   }

   #Trait=scale(Trait)
   Model1=data.frame(SNPA,SNPB,A,B,Trait)
   attr(Model1,"title")="SNPA -> A -> Trait -> B <- SNPB"

   Model1
}



if(exists("false.positive.study1") ) rm(false.positive.study1);
false.positive.study1=function() {

   s2=sqrt(2)

   # simulate
   no.samples=400
   set.seed(1)
   gen.var=1;
   env.var=.5;

   d1=date()

   num.sim=3
   model.num=1
   #   disc.matrix=matrix(nrow=6,ncol=6,data=0);
   eo.m = vector(length=num.sim);
   eo.m[]=0;
   lod.m=eo.m

   for (k in 1:num.sim) {

     mod=SimCR(no.samples=100,use.traitsnp=T, nAsnp=3,nBsnp=3, nTsnp=3, env.var=.5,randomize.snps=TRUE);

     q=neo(mod,skip.LEO=TRUE);
     
#     eo.m[k]=q$
     } # k

   d2=date()

  z=list();
  z$d1 = d1
  z$d2 = d2
  z$eo.m = eo.m

  save(z,file=paste(sep="","false.positive.study1.z.rdat"))    
}



if(exists("test.variance.logic") ) rm(test.variance.logic);
test.variance.logic=function() {
# todo for steve
#
# 3 excel files, one with snps
# one with details
# one with just LEO.I easy stuff.
#
# simulate 10 snps into A, 10 into B
# the noise SNPs -- simulate them as reactive to A
# so that they are correlated but not causal. (very clever).
#
# generate roc curve plots...
  #  for 0 <= omega/rho.ab <= 1
  # 
  # power at leo threshold should be very close to 100%, i.e. would like true positive rate to be close to 1.
  #
  # q25 = quantile(SNP, prob=.25)
  # q75 = quantile(SNP, prob=.75)
  # SNPprop = ifelse( SNPquant <= q25,0,ifelse(SNPquant <= q75,1,2))
  # 
  # --

 #
 # try with 10 true snps, 0 false ---->  1 true, 10 false.
 # using the false snps as simulated as reactive to A/B so they are correlated but not causal.
 # Thus we simulate LD.
#
# todo with jake: the srebp response pathway.
#
#

}

if(exists("generate.simulation.params.old") ) rm(generate.simulation.params.old);
generate.simulation.params.old=function(n.sample.points = 2) {

#  min.hb =  min.ha = .2 
#  max.hb =  max.ha = .8

  min.hb = min.ha = .5
  max.hb = max.ha = .5
  num.h.points= 1
  
  range.ha = max.ha - min.ha
  range.hb = max.hb - min.hb
  
  min.gamma2=0
  res=data.frame()

  ha.choices = seq(min.ha,max.ha,range.ha/num.h.points);
  for (ha in ha.choices) {
  #ha=ha.choices[1]
  
     hb.choices = seq(min.hb,max.hb,range.hb/num.h.points);
     for (hb in hb.choices) {
     #hb=hb.choices[1]
  
        min.errB = .05
        min.errA = .05
        max.gamma2 = round(min(1-hb-min.errB, 1-ha- min.errA, .6),3)

        gamma2.choices = seq(min.gamma2,max.gamma2,(max.gamma2-min.gamma2)/n.sample.points);
        for (gamma2 in gamma2.choices) {
        # gamma2=gamma2.choices[1]
  
           gamma=sqrt(gamma2)
           omega.min = 0;
           omega.max= round(-gamma2 + sqrt(gamma^4 +1 - gamma^2 - ha - min.errA),3);
           if (omega.max > .6) omega.max=.6
           omega.range = omega.max-omega.min;

           omega.choices = seq(omega.min, omega.max, omega.range/n.sample.points);
           for (omega in omega.choices) {
#omega = omega.choices[1]
  
              eA=1-2*omega*gamma^2 - gamma^2 - omega^2 - ha;
              eB=1-gamma2 - hb;
              rho.ab=gamma2+omega
              
              q=data.frame(ha,hb,gamma2,omega,eA,eB,rho.ab)
              res=rbind(res,q)
           }
        }
     }
   }
   round(res,3)
}



if(exists("generate.simulation.params") ) rm(generate.simulation.params);
generate.simulation.params=function() {

  min.hb =  min.ha = .2 
  max.hb =  max.ha = .6

#  min.hb = min.ha = .5
#  max.hb = max.ha = .5
  num.h.points= 3
  
  range.ha = max.ha - min.ha
  range.hb = max.hb - min.hb
  n.sample.points = 3
  omega.n.sample.points = 6

  min.gamma2=0
  res=data.frame()

  ha.choices = seq(min.ha,max.ha,range.ha/num.h.points);
  for (ha in ha.choices) {
  #ha=ha.choices[1]
  
     hb.choices = seq(min.hb,max.hb,range.hb/num.h.points);
     for (hb in hb.choices) {
     #hb=hb.choices[1]
  
        min.errB = .05
        min.errA = .05

        omega.top = .7
        
        omega.max = round(min(sqrt(1-ha-min.errA),omega.top),3)
        omega.min = 0
        
        omega.range = omega.max-omega.min;
        
        omega.choices = seq(omega.min, omega.max, omega.range/omega.n.sample.points)
        for (omega in omega.choices) {
          percent.max.omega = round(omega/omega.max,3)

          max.gamma2 = round(min((1-ha-omega^2-min.errA)/(1+2*omega),1-hb-min.errB,.6),3)
#          max.gamma2 = round(min(1-hb-min.errB, 1-ha- min.errA, .6),3)

        gamma2.choices = round(seq(min.gamma2,max.gamma2,(max.gamma2-min.gamma2)/n.sample.points),3);
        for (gamma2 in gamma2.choices) {
        # gamma2=gamma2.choices[1]
  
           gamma=sqrt(gamma2)
  
              eA=1-2*omega*gamma^2 - gamma^2 - omega^2 - ha;
              if (eA < 0) { stop("algebra error.") }
              eB=1-gamma2 - hb;
              if (eB < 0) { stop("algebra error.") }
              rho.ab=gamma2+omega
              if (rho.ab>1) { stop("algebra error.") }
              
              q=data.frame(ha,hb,gamma2,omega,eA,eB,rho.ab,percent.max.omega)
              res=rbind(res,q)
           }
        }
     }
   }
   round(res,3)
}

# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1)
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }

# ===================================================
# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)

err.bp<-function(daten,error,two.side=F){
 if(!is.numeric(daten)) {
      stop("All arguments must be numeric")}
 if(is.vector(daten)){ 
    xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
 }else{
    if (is.matrix(daten)){
      xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),dim=c(1,length(daten))))+0:(length(daten)-1)+.5
    }else{
      stop("First argument must either be a vector or a matrix") }
 }

 MW<-0.25*(max(xval)/length(xval))
 ERR1<-daten+error
 ERR2<-daten-error

 for(i in 1:length(daten)){
    segments(xval[i],daten[i],xval[i],ERR1[i])
    segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])

    if(two.side){
      segments(xval[i],daten[i],xval[i],ERR2[i])
      segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
    } 
 } 
}



if(exists("extract.chisq.right.log10prob") ) rm(extract.chisq.right.log10prob);
extract.chisq.right.log10prob=function(sem.obj) {
  if (!inherits(sem.obj,"summary.sem")) stop("bad object passed to extract.chisq.right.log10prob");

  -pchisq(sem.obj$chisq,sem.obj$df,lower.tail=FALSE,log.p=TRUE)*log.to.log10
}


if(exists("compare.single.vs.multiple.marker.simulations") ) rm(compare.single.vs.multiple.marker.simulations);
compare.single.vs.multiple.marker.simulations=function(param=NULL, no.samples=200, num.reps=50,pm=neo.get.param()) {

  #  param=generate.simulation.params.old(n.sample.points = 2)
  if (is.null(param)) { 
      param=generate.simulation.params();  
  } else {
      if (is.null(colnames(param))) { 
         colnames(param)=c("ha","hb","gamma2","omega","eA","eB","rho.ab","percent.max.omega")
         param=as.data.frame(param)
      }
  }

  N=no.samples
  set.seed(1)
  num.param=nrow(param)

  total.rows = num.reps*num.param

  param.keep=matrix(NA,nrow=total.rows,ncol=ncol(param))
  colnames(param.keep)=c("ha","hb","gamma2","omega","eA","eB","rho.ab","percent.max.omega")

 BLV.AtoB=rep(NA,total.rows)
 BLV.BtoA=rep(NA,total.rows)

 cor.ma.b=rep(NA,total.rows)
 pcor.ma.b.given.a=rep(NA,total.rows)
 cor.mb.a=rep(NA,total.rows)
 pcor.mb.a.given.b=rep(NA,total.rows)

 path.abs.z.f1ab.AtoB=rep(NA,total.rows)
 path.abs.z.f1ab.BtoA=rep(NA,total.rows)
 path.abs.z.f1ab.AcollideB=rep(NA,total.rows)
 path.abs.z.f1ab.BcollideA=rep(NA,total.rows)

 path.abs.z.f1ba.AtoB=rep(NA,total.rows)
 path.abs.z.f1ba.BtoA=rep(NA,total.rows)
 path.abs.z.f1ba.AcollideB=rep(NA,total.rows)
 path.abs.z.f1ba.BcollideA=rep(NA,total.rows)

 path.abs.z.f2.M.M1M2.AtoB=rep(NA,total.rows)
 path.abs.z.f2.M.M1M2.BtoA=rep(NA,total.rows)

 rmsea.f1ab.AtoB=rep(NA,total.rows)
 rmsea.f1ab.BtoA=rep(NA,total.rows)
 rmsea.f1ab.conf=rep(NA,total.rows)
 rmsea.f1ab.AcollideB=rep(NA,total.rows)
 rmsea.f1ab.BcollideA=rep(NA,total.rows)

 rmsea.f1ba.AtoB=rep(NA,total.rows)
 rmsea.f1ba.BtoA=rep(NA,total.rows)
 rmsea.f1ba.conf=rep(NA,total.rows)
 rmsea.f1ba.AcollideB=rep(NA,total.rows)
 rmsea.f1ba.BcollideA=rep(NA,total.rows)

 rmsea.f2.M.M1M2.AtoB=rep(NA,total.rows)
 rmsea.f2.M.M1M2.BtoA=rep(NA,total.rows)
 rmsea.f2.M.M1M2unresolv=rep(NA,total.rows)
 rmsea.f2.M.M1M2hidden.con=rep(NA,total.rows)

 rmsea90up.f1ab.AtoB=rep(NA,total.rows)
 rmsea90up.f1ab.BtoA=rep(NA,total.rows)
 rmsea90up.f1ab.conf=rep(NA,total.rows)
 rmsea90up.f1ab.AcollideB=rep(NA,total.rows)
 rmsea90up.f1ab.BcollideA=rep(NA,total.rows)

 rmsea90up.f1ba.AtoB=rep(NA,total.rows)
 rmsea90up.f1ba.BtoA=rep(NA,total.rows)
 rmsea90up.f1ba.conf=rep(NA,total.rows)
 rmsea90up.f1ba.AcollideB=rep(NA,total.rows)
 rmsea90up.f1ba.BcollideA=rep(NA,total.rows)

 rmsea90up.f2.M.M1M2.AtoB=rep(NA,total.rows)
 rmsea90up.f2.M.M1M2.BtoA=rep(NA,total.rows)
 rmsea90up.f2.M.M1M2unresolv=rep(NA,total.rows)
 rmsea90up.f2.M.M1M2hidden.con=rep(NA,total.rows)
  
 chisq.f1ab.AtoB=rep(NA,total.rows)
 chisq.f1ab.BtoA=rep(NA,total.rows)
 chisq.f1ab.conf=rep(NA,total.rows)
 chisq.f1ab.AcollideB=rep(NA,total.rows)
 chisq.f1ab.BcollideA=rep(NA,total.rows)

 chisq.f1ba.AtoB=rep(NA,total.rows)
 chisq.f1ba.BtoA=rep(NA,total.rows)
 chisq.f1ba.conf=rep(NA,total.rows)
 chisq.f1ba.AcollideB=rep(NA,total.rows)
 chisq.f1ba.BcollideA=rep(NA,total.rows)

 chisq.f2.M.M1M2.AtoB=rep(NA,total.rows)
 chisq.f2.M.M1M2.BtoA=rep(NA,total.rows)
 chisq.f2.M.M1M2unresolv=rep(NA,total.rows)
 chisq.f2.M.M1M2hidden.con=rep(NA,total.rows)
  
 chiprob.f1ab.AtoB=rep(NA,total.rows)
 chiprob.f1ab.BtoA=rep(NA,total.rows)
 chiprob.f1ab.conf=rep(NA,total.rows)
 chiprob.f1ab.AcollideB=rep(NA,total.rows)
 chiprob.f1ab.BcollideA=rep(NA,total.rows)

 chiprob.f1ba.AtoB=rep(NA,total.rows)
 chiprob.f1ba.BtoA=rep(NA,total.rows)
 chiprob.f1ba.conf=rep(NA,total.rows)
 chiprob.f1ba.AcollideB=rep(NA,total.rows)
 chiprob.f1ba.BcollideA=rep(NA,total.rows)

 chiprob.f2.M.M1M2.AtoB=rep(NA,total.rows)
 chiprob.f2.M.M1M2.BtoA=rep(NA,total.rows)
 chiprob.f2.M.M1M2unresolv=rep(NA,total.rows)
 chiprob.f2.M.M1M2hidden.con=rep(NA,total.rows)
  
  res.list=vector("list", total.rows)
  
 for (i in 1:num.param) {
   for (j in 1:num.reps) {
      k=(num.reps*(i-1))+j # where to store
#  for (i in 1:1) {
#    for (j in 1:2) {
   # q=data.frame(ha,hb,gamma2,omega,eA,eB,rho.ab)

     z =list()

     z$ha =   ha     = param$ha[i]
     z$hb =   hb     = param$hb[i]
     z$gamma2 =   gamma2 = param$gamma2[i]
     z$omega =   omega  = param$omega[i]
     z$eA =   eA =     param$eA[i]
     z$eB =   eB =     param$eB[i]
     z$rho.ab =   rho.ab = param$rho.ab[i]

     param.string=paste(sep="","ha=",ha, " hb=",hb, " g2=",gamma2, " eA=",eA," eB=",eB, " rho.ab=",rho.ab, " w=",omega)
   
#              varB= gamma2 + eB + hb;
#              varA = TA + 2*omega*gamma*gamma + gamma*gamma + omega*omega*varB + eA;

              # if B is standardized to var(B)=1, then....

#  std.varA = TA + gamma^2 + omega^2 + 2*gamma*gamma*omega + eA; # set == 1
#  covAB = omega*TB+gamma*gamma*(1+omega)+omega*eB
#  covAC = gamma+ omega*gamma
#  covBC = gamma

#  cov.SNPB.A = sqrt(TB)*omega

#  rho.AB = gamma^2 + omega

  SNPA=rnorm(no.samples);
  SNPB=rnorm(no.samples);

  C=rnorm(no.samples);

  errA=rnorm(no.samples,sd=sqrt(eA));
  errB=rnorm(no.samples,sd=sqrt(eB));

  B=errB+sqrt(gamma2)*C+ sqrt(hb)*SNPB;

  A=errA+omega*B+sqrt(gamma2)*C+sqrt(ha)*SNPA;

   
  x=data.frame(SNPA,SNPB,A,B);
   
    z$x=x
        z$fit1mb = compare.local.sems(pm=pm,M.col=2,A.col=3,B.col=4,x=x,no.obs=no.samples)
        z$fit1ma = compare.local.sems(pm=pm,M.col=1,A.col=3,B.col=4,x=x,no.obs=no.samples)
        z$fit2m =local.sem.four.var.m1m2(pm=pm,MA.col=1,A.col=3,B.col=4,x=x,MB.col=2,no.obs=N,fit.models=c(1,2,3,4));
    z$param.string = param.string;

 param.keep[k,]=as.matrix(param[i,]) # track params too


 BLV.AtoB[k]=z$fit2m$zeo.MA.B.given.A$BLV
 BLV.BtoA[k]=z$fit2m$zeo.MB.A.given.B$BLV

 cor.ma.b[k]         =z$fit2m$zeo.MA.B.given.A$cor
 pcor.ma.b.given.a[k]=z$fit2m$zeo.MA.B.given.A$pcor

 cor.mb.a[k]         =z$fit2m$zeo.MB.A.given.B$cor
 pcor.mb.a.given.b[k]=z$fit2m$zeo.MB.A.given.B$pcor


 



 path.abs.z.f1ab.AtoB[k]=abs(z$fit1ma$M.AtoB$coef["fr.A.B",3])
 path.abs.z.f1ab.BtoA[k]=abs(z$fit1ma$M.BtoA$coef["fr.B.A",3])
 path.abs.z.f1ab.AcollideB[k]=abs(z$fit1ma$M.AcollideB$coef["fr.B.A",3])
 path.abs.z.f1ab.BcollideA[k]=abs(z$fit1ma$M.BcollideA$coef["fr.A.B",3])

 path.abs.z.f1ba.AtoB[k]=abs(z$fit1mb$M.AtoB$coef["fr.A.B",3])
 path.abs.z.f1ba.BtoA[k]=abs(z$fit1mb$M.BtoA$coef["fr.B.A",3])
 path.abs.z.f1ba.AcollideB[k]=abs(z$fit1mb$M.BcollideA$coef["fr.A.B",3])
 path.abs.z.f1ba.BcollideA[k]=abs(z$fit1mb$M.AcollideB$coef["fr.B.A",3])

 path.abs.z.f2.M.M1M2.AtoB[k]=abs(z$fit2m$M.M1M2.AtoB$coef["fr.A.B",3])
 path.abs.z.f2.M.M1M2.BtoA[k]=abs(z$fit2m$M.M1M2.BtoA$coef["fr.B.A",3])

 rmsea.f1ab.AtoB[k]=z$fit1ma$M.AtoB$RMSEA[1]
 rmsea.f1ab.BtoA[k]=z$fit1ma$M.BtoA$RMSEA[1]
 rmsea.f1ab.conf[k]=z$fit1ma$M.conf$RMSEA[1]
 rmsea.f1ab.AcollideB[k]=z$fit1ma$M.AcollideB$RMSEA[1]
 rmsea.f1ab.BcollideA[k]=z$fit1ma$M.BcollideA$RMSEA[1]

 rmsea.f1ba.AtoB[k]=z$fit1mb$M.AtoB$RMSEA[1]
 rmsea.f1ba.BtoA[k]=z$fit1mb$M.BtoA$RMSEA[1]
 rmsea.f1ba.conf[k]=z$fit1mb$M.conf$RMSEA[1]
 rmsea.f1ba.AcollideB[k]=z$fit1mb$M.AcollideB$RMSEA[1]
 rmsea.f1ba.BcollideA[k]=z$fit1mb$M.BcollideA$RMSEA[1]

 rmsea.f2.M.M1M2.AtoB[k]=z$fit2m$M.M1M2.AtoB$RMSEA[1]
 rmsea.f2.M.M1M2.BtoA[k]=z$fit2m$M.M1M2.BtoA$RMSEA[1]
 rmsea.f2.M.M1M2unresolv[k]=z$fit2m$M.M1M2unresolv$RMSEA[1]
 rmsea.f2.M.M1M2hidden.con[k]=z$fit2m$M.M1M2hidden.con$RMSEA[1]

 rmsea90up.f1ab.AtoB[k]=z$fit1ma$M.AtoB$RMSEA[3]
 rmsea90up.f1ab.BtoA[k]=z$fit1ma$M.BtoA$RMSEA[3]
 rmsea90up.f1ab.conf[k]=z$fit1ma$M.conf$RMSEA[3]
 rmsea90up.f1ab.AcollideB[k]=z$fit1ma$M.AcollideB$RMSEA[3]
 rmsea90up.f1ab.BcollideA[k]=z$fit1ma$M.BcollideA$RMSEA[3]

 rmsea90up.f1ba.AtoB[k]=z$fit1mb$M.AtoB$RMSEA[3]
 rmsea90up.f1ba.BtoA[k]=z$fit1mb$M.BtoA$RMSEA[3]
 rmsea90up.f1ba.conf[k]=z$fit1mb$M.conf$RMSEA[3]
 rmsea90up.f1ba.AcollideB[k]=z$fit1mb$M.AcollideB$RMSEA[3]
 rmsea90up.f1ba.BcollideA[k]=z$fit1mb$M.BcollideA$RMSEA[3]

 rmsea90up.f2.M.M1M2.AtoB[k]=z$fit2m$M.M1M2.AtoB$RMSEA[3]
 rmsea90up.f2.M.M1M2.BtoA[k]=z$fit2m$M.M1M2.BtoA$RMSEA[3]
 rmsea90up.f2.M.M1M2unresolv[k]=z$fit2m$M.M1M2unresolv$RMSEA[3]
 rmsea90up.f2.M.M1M2hidden.con[k]=z$fit2m$M.M1M2hidden.con$RMSEA[3]
  
 chisq.f1ab.AtoB[k]=z$fit1ma$M.AtoB$chisq
 chisq.f1ab.BtoA[k]=z$fit1ma$M.BtoA$chisq
 chisq.f1ab.conf[k]=z$fit1ma$M.conf$chisq
 chisq.f1ab.AcollideB[k]=z$fit1ma$M.AcollideB$chisq
 chisq.f1ab.BcollideA[k]=z$fit1ma$M.BcollideA$chisq

 chisq.f1ba.AtoB[k]=z$fit1mb$M.AtoB$chisq
 chisq.f1ba.BtoA[k]=z$fit1mb$M.BtoA$chisq
 chisq.f1ba.conf[k]=z$fit1mb$M.conf$chisq
 chisq.f1ba.AcollideB[k]=z$fit1mb$M.AcollideB$chisq
 chisq.f1ba.BcollideA[k]=z$fit1mb$M.BcollideA$chisq

 chisq.f2.M.M1M2.AtoB[k]=z$fit2m$M.M1M2.AtoB$chisq
 chisq.f2.M.M1M2.BtoA[k]=z$fit2m$M.M1M2.BtoA$chisq
 chisq.f2.M.M1M2unresolv[k]=z$fit2m$M.M1M2unresolv$chisq
 chisq.f2.M.M1M2hidden.con[k]=z$fit2m$M.M1M2hidden.con$chisq

 chiprob.f1ab.AtoB[k]=extract.chisq.right.log10prob(z$fit1ma$M.AtoB)
 chiprob.f1ab.BtoA[k]=extract.chisq.right.log10prob(z$fit1ma$M.BtoA)
 chiprob.f1ab.conf[k]=extract.chisq.right.log10prob(z$fit1ma$M.conf)
 chiprob.f1ab.AcollideB[k]=extract.chisq.right.log10prob(z$fit1ma$M.AcollideB)
 chiprob.f1ab.BcollideA[k]=extract.chisq.right.log10prob(z$fit1ma$M.BcollideA)

 chiprob.f1ba.AtoB[k]=extract.chisq.right.log10prob(z$fit1mb$M.AtoB)
 chiprob.f1ba.BtoA[k]=extract.chisq.right.log10prob(z$fit1mb$M.BtoA)
 chiprob.f1ba.conf[k]=extract.chisq.right.log10prob(z$fit1mb$M.conf)
 chiprob.f1ba.AcollideB[k]=extract.chisq.right.log10prob(z$fit1mb$M.AcollideB)
 chiprob.f1ba.BcollideA[k]=extract.chisq.right.log10prob(z$fit1mb$M.BcollideA)

 chiprob.f2.M.M1M2.AtoB[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2.AtoB)
 chiprob.f2.M.M1M2.BtoA[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2.BtoA)
 chiprob.f2.M.M1M2unresolv[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2unresolv)
 chiprob.f2.M.M1M2hidden.con[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2hidden.con)
     
    res.list[[(i-1)*num.reps + j]]=z

   } # end j
#    save(res.list,file="my.res.list.rdat")

 } # end i over params

 scores=data.frame(param.keep,rmsea.f1ab.AtoB, rmsea.f1ab.BtoA, 
 rmsea.f1ab.conf, rmsea.f1ab.AcollideB, rmsea.f1ab.BcollideA,
 rmsea.f1ba.AtoB, rmsea.f1ba.BtoA, rmsea.f1ba.conf, rmsea.f1ba.AcollideB, 
 rmsea.f1ba.BcollideA, rmsea.f2.M.M1M2.AtoB, rmsea.f2.M.M1M2.BtoA, 
 rmsea.f2.M.M1M2unresolv, rmsea.f2.M.M1M2hidden.con, rmsea90up.f1ab.AtoB, 
 rmsea90up.f1ab.BtoA, rmsea90up.f1ab.conf, rmsea90up.f1ab.AcollideB, rmsea90up.f1ab.BcollideA,
 rmsea90up.f1ba.AtoB, rmsea90up.f1ba.BtoA, rmsea90up.f1ba.conf, rmsea90up.f1ba.AcollideB, rmsea90up.f1ba.BcollideA,
 rmsea90up.f2.M.M1M2.AtoB, rmsea90up.f2.M.M1M2.BtoA, rmsea90up.f2.M.M1M2unresolv, rmsea90up.f2.M.M1M2hidden.con,  
 chisq.f1ab.AtoB, chisq.f1ab.BtoA, chisq.f1ab.conf, chisq.f1ab.AcollideB, chisq.f1ab.BcollideA,
 chisq.f1ba.AtoB, chisq.f1ba.BtoA, chisq.f1ba.conf, chisq.f1ba.AcollideB, chisq.f1ba.BcollideA,
 chisq.f2.M.M1M2.AtoB, chisq.f2.M.M1M2.BtoA, chisq.f2.M.M1M2unresolv, chisq.f2.M.M1M2hidden.con,  
 chiprob.f1ab.AtoB, chiprob.f1ab.BtoA, chiprob.f1ab.conf, chiprob.f1ab.AcollideB, chiprob.f1ab.BcollideA,
 chiprob.f1ba.AtoB, chiprob.f1ba.BtoA, chiprob.f1ba.conf, chiprob.f1ba.AcollideB, chiprob.f1ba.BcollideA,
 chiprob.f2.M.M1M2.AtoB, chiprob.f2.M.M1M2.BtoA, chiprob.f2.M.M1M2unresolv, chiprob.f2.M.M1M2hidden.con,
 path.abs.z.f1ab.AtoB, path.abs.z.f1ab.BtoA,path.abs.z.f1ab.AcollideB,
 path.abs.z.f1ab.BcollideA,path.abs.z.f1ba.AtoB,path.abs.z.f1ba.BtoA,
 path.abs.z.f1ba.AcollideB,path.abs.z.f1ba.BcollideA,path.abs.z.f2.M.M1M2.AtoB,
 path.abs.z.f2.M.M1M2.BtoA, 
 BLV.AtoB, BLV.BtoA, 
 cor.ma.b, pcor.ma.b.given.a, cor.mb.a, pcor.mb.a.given.b)

 save(res.list,scores,file="my.new.res.list.scores.rdat")

 #### RETURN EARLY
 return(scores);

 # analysis before completion...
 # scores.new = scores ... b/c we loaded it separately.

# examples:
#     boxplot(split(log10(d$leo.i.correct[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",at = 1:length(unique(d$ntrack[all])) - 0.2,boxwex = 0.25,col = "yellow")#,ylim=c(-20,20))
#   boxplot(split(log10(d$leo.i[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",at = 1:length(unique(d$ntrack[all])) + 0.2,add=TRUE,boxwex = 0.25,col = "orange",show.names=FALSE)#,ylim=c(-20,20))

 file.out="power.study.rmsea.f2.AtoB.over5600.ps"
 postscript(file=file.out,horizontal=FALSE);
 boxplot(split(scores.new$rmsea.f2.M.M1M2.BtoA,scores.new$percent.max.omega),notch=TRUE,xlab="percent omega.max",ylab="RMSEA",col="yellow",at=1:length(unique(scores.new$percent.max.omega))-.2,boxwex=.25,ylim=c(-.15,.7))
 boxplot(split(scores.new$rmsea.f2.M.M1M2.AtoB,scores.new$percent.max.omega),notch=TRUE,col="violetred",at=1:length(unique(scores.new$percent.max.omega))+.2,boxwex=.25,add=TRUE,show.names=FALSE,main="Monte Carlo study: B causal for A" )
  legend(5.1,-.05, c("correct causal B to A", "incorrect A to B model"), fill = c("yellow", "violetred"))
  dev.off()

# things that vary
  # hb = .2, .333, .467, .6
  #
  # at hb = .2
  #    omega = 0, .117, .233, .350, .583, .700

  # start by studying the null model: no causal, no confounding
  # under different levels of genetic influence of HB

  no.gamma.no.omega=c(1:50,1401:1450,2801:2850,4201:4250)
  s=scores.new[no.gamma.no.omega,]

  file.out="power.study.rmsea.no.gamma.no.omega.hb.varies.ps"
  postscript(file=file.out,horizontal=FALSE);
  boxplot(split(s$rmsea.f2.M.M1M2.BtoA,s$hb),notch=TRUE,xlab="percent omega.max",ylab="RMSEA",col="yellow",at=1:length(unique(s$hb))-.2,boxwex=.25,log="y")#,ylim=c(-.15,.7))
  boxplot(split(s$rmsea.f2.M.M1M2.AtoB,s$percent.max.omega),notch=TRUE,col="violetred",at=1:length(unique(s$percent.max.omega))+.2,boxwex=.25,add=TRUE,show.names=FALSE,main="Monte Carlo study: B causal for A" )
  legend(5.1,-.05, c("correct causal B to A", "incorrect A to B model"), fill = c("yellow", "violetred"))
  dev.off()


  # make barplots with error bars
  # ex:
  # freq.bar(s,threshold=.1,th.col="rmsea.f2.M.M1M2.BtoA",group.col="percent.max.omega")

  
  
 rownames(scores) = as.character(1:nrow(scores))
 scores  
} # compare.single.vs.multiple.marker.simulations



if(exists("create.scores.from.res.list")) rm (create.scores.from.res.list);
create.scores.from.res.list=function(res.list, num.reps=50, no.samples=200) {

# we used to re-read the param from the generate.simulation.params function
# but now we read them exactly from the saved res.list to be accurate,
# since the res.list may be from an older generate.simulation.params call.
#
#  param=generate.simulation.params();
#   num.param=nrow(param)


  N=no.samples

 # count the good rows in res.list
 good.rows=0
 for (k in 1:length(res.list)) { if (!is.null(res.list[[k]])) good.rows=good.rows+1 }
 num.param = good.rows / num.reps

   param=matrix(NA,nrow=num.param,ncol=8)

  j=0;
        omega.top = .7
        min.errA=.05
  for (k in seq(1,good.rows,num.reps)) {
    j=j+1
    param[j,1]=res.list[[k]]$ha
    param[j,2]=res.list[[k]]$hb
    param[j,3]=res.list[[k]]$gamma2
    param[j,4]=res.list[[k]]$omega
    param[j,5]=res.list[[k]]$eA
    param[j,6]=res.list[[k]]$eB
    param[j,7]=res.list[[k]]$rho.ab

    # compute percent.max.omega
    omega.max = round(min(sqrt(1-res.list[[k]]$ha-min.errA),omega.top),3)

    param[j,8]= round(res.list[[k]]$omega / omega.max,3)
  }

 total.rows=num.reps*num.param

 param.keep=matrix(NA,nrow=total.rows,ncol=ncol(param))

 path.abs.z.f1ab.AtoB[k]=rep(NA,total.rows)
 path.abs.z.f1ab.BtoA[k]=rep(NA,total.rows)
 path.abs.z.f1ab.AcollideB[k]=rep(NA,total.rows)
 path.abs.z.f1ab.BcollideA[k]=rep(NA,total.rows)

 path.abs.z.f1ba.AtoB[k]=rep(NA,total.rows)
 path.abs.z.f1ba.BtoA[k]=rep(NA,total.rows)
 path.abs.z.f1ba.AcollideB[k]=rep(NA,total.rows)
 path.abs.z.f1ba.BcollideA[k]=rep(NA,total.rows)

 path.abs.z.f2.M.M1M2.AtoB[k]=rep(NA,total.rows)
 path.abs.z.f2.M.M1M2.BtoA[k]=rep(NA,total.rows)

 rmsea.f1ab.AtoB=rep(NA,total.rows)
 rmsea.f1ab.BtoA=rep(NA,total.rows)
 rmsea.f1ab.conf=rep(NA,total.rows)
 rmsea.f1ab.AcollideB=rep(NA,total.rows)
 rmsea.f1ab.BcollideA=rep(NA,total.rows)

 rmsea.f1ba.AtoB=rep(NA,total.rows)
 rmsea.f1ba.BtoA=rep(NA,total.rows)
 rmsea.f1ba.conf=rep(NA,total.rows)
 rmsea.f1ba.AcollideB=rep(NA,total.rows)
 rmsea.f1ba.BcollideA=rep(NA,total.rows)

 rmsea.f2.M.M1M2.AtoB=rep(NA,total.rows)
 rmsea.f2.M.M1M2.BtoA=rep(NA,total.rows)
 rmsea.f2.M.M1M2unresolv=rep(NA,total.rows)
 rmsea.f2.M.M1M2hidden.con=rep(NA,total.rows)

 rmsea90up.f1ab.AtoB=rep(NA,total.rows)
 rmsea90up.f1ab.BtoA=rep(NA,total.rows)
 rmsea90up.f1ab.conf=rep(NA,total.rows)
 rmsea90up.f1ab.AcollideB=rep(NA,total.rows)
 rmsea90up.f1ab.BcollideA=rep(NA,total.rows)

 rmsea90up.f1ba.AtoB=rep(NA,total.rows)
 rmsea90up.f1ba.BtoA=rep(NA,total.rows)
 rmsea90up.f1ba.conf=rep(NA,total.rows)
 rmsea90up.f1ba.AcollideB=rep(NA,total.rows)
 rmsea90up.f1ba.BcollideA=rep(NA,total.rows)

 rmsea90up.f2.M.M1M2.AtoB=rep(NA,total.rows)
 rmsea90up.f2.M.M1M2.BtoA=rep(NA,total.rows)
 rmsea90up.f2.M.M1M2unresolv=rep(NA,total.rows)
 rmsea90up.f2.M.M1M2hidden.con=rep(NA,total.rows)
  
 chisq.f1ab.AtoB=rep(NA,total.rows)
 chisq.f1ab.BtoA=rep(NA,total.rows)
 chisq.f1ab.conf=rep(NA,total.rows)
 chisq.f1ab.AcollideB=rep(NA,total.rows)
 chisq.f1ab.BcollideA=rep(NA,total.rows)

 chisq.f1ba.AtoB=rep(NA,total.rows)
 chisq.f1ba.BtoA=rep(NA,total.rows)
 chisq.f1ba.conf=rep(NA,total.rows)
 chisq.f1ba.AcollideB=rep(NA,total.rows)
 chisq.f1ba.BcollideA=rep(NA,total.rows)

 chisq.f2.M.M1M2.AtoB=rep(NA,total.rows)
 chisq.f2.M.M1M2.BtoA=rep(NA,total.rows)
 chisq.f2.M.M1M2unresolv=rep(NA,total.rows)
 chisq.f2.M.M1M2hidden.con=rep(NA,total.rows)
  
 chiprob.f1ab.AtoB=rep(NA,total.rows)
 chiprob.f1ab.BtoA=rep(NA,total.rows)
 chiprob.f1ab.conf=rep(NA,total.rows)
 chiprob.f1ab.AcollideB=rep(NA,total.rows)
 chiprob.f1ab.BcollideA=rep(NA,total.rows)

 chiprob.f1ba.AtoB=rep(NA,total.rows)
 chiprob.f1ba.BtoA=rep(NA,total.rows)
 chiprob.f1ba.conf=rep(NA,total.rows)
 chiprob.f1ba.AcollideB=rep(NA,total.rows)
 chiprob.f1ba.BcollideA=rep(NA,total.rows)

 chiprob.f2.M.M1M2.AtoB=rep(NA,total.rows)
 chiprob.f2.M.M1M2.BtoA=rep(NA,total.rows)
 chiprob.f2.M.M1M2unresolv=rep(NA,total.rows)
 chiprob.f2.M.M1M2hidden.con=rep(NA,total.rows)

 done=FALSE;  
 for (i in 1:num.param) {
   if (done) break;
   for (j in 1:num.reps) {

 z=    res.list[[(i-1)*num.reps + j]]
 if (is.null(z)) { done=TRUE; break; }

 path.abs.z.f1ab.AtoB[k]=abs(z$fit1ma$M.AtoB$coef["fr.A.B",3])
 path.abs.z.f1ab.BtoA[k]=abs(z$fit1ma$M.BtoA$coef["fr.B.A",3])
 path.abs.z.f1ab.AcollideB[k]=abs(z$fit1ma$M.AcollideB$coef["fr.B.A",3])
 path.abs.z.f1ab.BcollideA[k]=abs(z$fit1ma$M.BcollideA$coef["fr.A.B",3])

 path.abs.z.f1ba.AtoB[k]=abs(z$fit1mb$M.AtoB$coef["fr.A.B",3])
 path.abs.z.f1ba.BtoA[k]=abs(z$fit1mb$M.BtoA$coef["fr.B.A",3])
 path.abs.z.f1ba.AcollideB[k]=abs(z$fit1mb$M.BcollideA$coef["fr.A.B",3])
 path.abs.z.f1ba.BcollideA[k]=abs(z$fit1mb$M.AcollideB$coef["fr.B.A",3])

 path.abs.z.f2.M.M1M2.AtoB[k]=abs(z$fit2m$M.M1M2.AtoB$coef["fr.A.B",3])
 path.abs.z.f2.M.M1M2.BtoA[k]=abs(z$fit2m$M.M1M2.BtoA$coef["fr.B.A",3])


 param.keep[k,]=as.matrix(param[i,]) # track params too

 rmsea.f1ab.AtoB[k]=z$fit1ma$M.AtoB$RMSEA[1]
 rmsea.f1ab.BtoA[k]=z$fit1ma$M.BtoA$RMSEA[1]
 rmsea.f1ab.conf[k]=z$fit1ma$M.conf$RMSEA[1]
 rmsea.f1ab.AcollideB[k]=z$fit1ma$M.AcollideB$RMSEA[1]
 rmsea.f1ab.BcollideA[k]=z$fit1ma$M.BcollideA$RMSEA[1]

 rmsea.f1ba.AtoB[k]=z$fit1mb$M.AtoB$RMSEA[1]
 rmsea.f1ba.BtoA[k]=z$fit1mb$M.BtoA$RMSEA[1]
 rmsea.f1ba.conf[k]=z$fit1mb$M.conf$RMSEA[1]
 rmsea.f1ba.AcollideB[k]=z$fit1mb$M.AcollideB$RMSEA[1]
 rmsea.f1ba.BcollideA[k]=z$fit1mb$M.BcollideA$RMSEA[1]

 rmsea.f2.M.M1M2.AtoB[k]=z$fit2m$M.M1M2.AtoB$RMSEA[1]
 rmsea.f2.M.M1M2.BtoA[k]=z$fit2m$M.M1M2.BtoA$RMSEA[1]
 rmsea.f2.M.M1M2unresolv[k]=z$fit2m$M.M1M2unresolv$RMSEA[1]
 rmsea.f2.M.M1M2hidden.con[k]=z$fit2m$M.M1M2hidden.con$RMSEA[1]

 rmsea90up.f1ab.AtoB[k]=z$fit1ma$M.AtoB$RMSEA[3]
 rmsea90up.f1ab.BtoA[k]=z$fit1ma$M.BtoA$RMSEA[3]
 rmsea90up.f1ab.conf[k]=z$fit1ma$M.conf$RMSEA[3]
 rmsea90up.f1ab.AcollideB[k]=z$fit1ma$M.AcollideB$RMSEA[3]
 rmsea90up.f1ab.BcollideA[k]=z$fit1ma$M.BcollideA$RMSEA[3]

 rmsea90up.f1ba.AtoB[k]=z$fit1mb$M.AtoB$RMSEA[3]
 rmsea90up.f1ba.BtoA[k]=z$fit1mb$M.BtoA$RMSEA[3]
 rmsea90up.f1ba.conf[k]=z$fit1mb$M.conf$RMSEA[3]
 rmsea90up.f1ba.AcollideB[k]=z$fit1mb$M.AcollideB$RMSEA[3]
 rmsea90up.f1ba.BcollideA[k]=z$fit1mb$M.BcollideA$RMSEA[3]

 rmsea90up.f2.M.M1M2.AtoB[k]=z$fit2m$M.M1M2.AtoB$RMSEA[3]
 rmsea90up.f2.M.M1M2.BtoA[k]=z$fit2m$M.M1M2.BtoA$RMSEA[3]
 rmsea90up.f2.M.M1M2unresolv[k]=z$fit2m$M.M1M2unresolv$RMSEA[3]
 rmsea90up.f2.M.M1M2hidden.con[k]=z$fit2m$M.M1M2hidden.con$RMSEA[3]
  
 chisq.f1ab.AtoB[k]=z$fit1ma$M.AtoB$chisq
 chisq.f1ab.BtoA[k]=z$fit1ma$M.BtoA$chisq
 chisq.f1ab.conf[k]=z$fit1ma$M.conf$chisq
 chisq.f1ab.AcollideB[k]=z$fit1ma$M.AcollideB$chisq
 chisq.f1ab.BcollideA[k]=z$fit1ma$M.BcollideA$chisq

 chisq.f1ba.AtoB[k]=z$fit1mb$M.AtoB$chisq
 chisq.f1ba.BtoA[k]=z$fit1mb$M.BtoA$chisq
 chisq.f1ba.conf[k]=z$fit1mb$M.conf$chisq
 chisq.f1ba.AcollideB[k]=z$fit1mb$M.AcollideB$chisq
 chisq.f1ba.BcollideA[k]=z$fit1mb$M.BcollideA$chisq

 chisq.f2.M.M1M2.AtoB[k]=z$fit2m$M.M1M2.AtoB$chisq
 chisq.f2.M.M1M2.BtoA[k]=z$fit2m$M.M1M2.BtoA$chisq
 chisq.f2.M.M1M2unresolv[k]=z$fit2m$M.M1M2unresolv$chisq
 chisq.f2.M.M1M2hidden.con[k]=z$fit2m$M.M1M2hidden.con$chisq

 chiprob.f1ab.AtoB[k]=extract.chisq.right.log10prob(z$fit1ma$M.AtoB)
 chiprob.f1ab.BtoA[k]=extract.chisq.right.log10prob(z$fit1ma$M.BtoA)
 chiprob.f1ab.conf[k]=extract.chisq.right.log10prob(z$fit1ma$M.conf)
 chiprob.f1ab.AcollideB[k]=extract.chisq.right.log10prob(z$fit1ma$M.AcollideB)
 chiprob.f1ab.BcollideA[k]=extract.chisq.right.log10prob(z$fit1ma$M.BcollideA)

 chiprob.f1ba.AtoB[k]=extract.chisq.right.log10prob(z$fit1mb$M.AtoB)
 chiprob.f1ba.BtoA[k]=extract.chisq.right.log10prob(z$fit1mb$M.BtoA)
 chiprob.f1ba.conf[k]=extract.chisq.right.log10prob(z$fit1mb$M.conf)
 chiprob.f1ba.AcollideB[k]=extract.chisq.right.log10prob(z$fit1mb$M.AcollideB)
 chiprob.f1ba.BcollideA[k]=extract.chisq.right.log10prob(z$fit1mb$M.BcollideA)

 chiprob.f2.M.M1M2.AtoB[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2.AtoB)
 chiprob.f2.M.M1M2.BtoA[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2.BtoA)
 chiprob.f2.M.M1M2unresolv[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2unresolv)
 chiprob.f2.M.M1M2hidden.con[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2hidden.con)
     

   } # end j

 } # end i over params

     
 scores=data.frame(param.keep,rmsea.f1ab.AtoB, rmsea.f1ab.BtoA, rmsea.f1ab.conf, rmsea.f1ab.AcollideB, rmsea.f1ab.BcollideA,
 rmsea.f1ba.AtoB, rmsea.f1ba.BtoA, rmsea.f1ba.conf, rmsea.f1ba.AcollideB, rmsea.f1ba.BcollideA,
 rmsea.f2.M.M1M2.AtoB, rmsea.f2.M.M1M2.BtoA, rmsea.f2.M.M1M2unresolv, rmsea.f2.M.M1M2hidden.con,
 rmsea90up.f1ab.AtoB, rmsea90up.f1ab.BtoA, rmsea90up.f1ab.conf, rmsea90up.f1ab.AcollideB, rmsea90up.f1ab.BcollideA,
 rmsea90up.f1ba.AtoB, rmsea90up.f1ba.BtoA, rmsea90up.f1ba.conf, rmsea90up.f1ba.AcollideB, rmsea90up.f1ba.BcollideA,
 rmsea90up.f2.M.M1M2.AtoB, rmsea90up.f2.M.M1M2.BtoA, rmsea90up.f2.M.M1M2unresolv, rmsea90up.f2.M.M1M2hidden.con,  
 chisq.f1ab.AtoB, chisq.f1ab.BtoA, chisq.f1ab.conf, chisq.f1ab.AcollideB, chisq.f1ab.BcollideA,
 chisq.f1ba.AtoB, chisq.f1ba.BtoA, chisq.f1ba.conf, chisq.f1ba.AcollideB, chisq.f1ba.BcollideA,
 chisq.f2.M.M1M2.AtoB, chisq.f2.M.M1M2.BtoA, chisq.f2.M.M1M2unresolv, chisq.f2.M.M1M2hidden.con,  
 chiprob.f1ab.AtoB, chiprob.f1ab.BtoA, chiprob.f1ab.conf, chiprob.f1ab.AcollideB, chiprob.f1ab.BcollideA,
 chiprob.f1ba.AtoB, chiprob.f1ba.BtoA, chiprob.f1ba.conf, chiprob.f1ba.AcollideB, chiprob.f1ba.BcollideA,
 chiprob.f2.M.M1M2.AtoB, chiprob.f2.M.M1M2.BtoA, chiprob.f2.M.M1M2unresolv, chiprob.f2.M.M1M2hidden.con)

# save(scores,file="scores.rdat")

 scores
}

if(exists("check.variance.algebra") ) rm(check.variance.algebra);
check.variance.algebra=function(no.samples=100000) {

  TA=.3;
  TB=.5;
  gamma=.05;

  
  
  errA.min = .05
  omega.max= -gamma^2 + sqrt(gamma^4 +1 - gamma^2 - TA - errA.min)

  omega=.5;

  eA=1-2*omega*gamma^2 - gamma^2 - omega^2 - TA;
  eB=1-gamma^2 - TB;

  varB= gamma^2 + eB + TB;
  varA = TA + 2*omega*gamma*gamma + gamma*gamma + omega*omega*varB + eA;
  # if B is standardized to var(B)=1, then....

  std.varA = TA + gamma^2 + omega^2 + 2*gamma*gamma*omega + eA; # set == 1
  covAB = omega*TB+gamma*gamma*(1+omega)+omega*eB
  covAC = gamma+ omega*gamma
  covBC = gamma

  cov.SNPB.A = sqrt(TB)*omega

  rho.AB = gamma^2 + omega

  SNPA=rnorm(no.samples);
  SNPB=rnorm(no.samples);

  C=rnorm(no.samples);

  errA=rnorm(no.samples,sd=sqrt(eA));
  errB=rnorm(no.samples,sd=sqrt(eB));

  B=errB+gamma*C+ sqrt(TB)*SNPB;

  A=errA+omega*B+gamma*C+sqrt(TA)*SNPA;

  df=data.frame(A,B,C,SNPA,SNPB)
  cx=cor(df)
  

   # start the simulation
   ha =.5 # direct broad-sense heritability of A: the percentage of variance of A that is due to the genetic markers directly influence A, and not via B.
   hb =.5 # the broad-sense heritability of B: the percentage of variance of B that is due to the genetic loci influencing B.


   npoints = 5 # number of points at which to evaluate gamma^2
   min.err.var.b=0.05
   min.gamma2=.05
   max.gamma2=1-hb-min.err.var.b
   gamma2 = seq(min.gamma2,max.gamma2,(max.gamma2-min.gamma2)/npoints);
   rho = rev(gamma2);
   omega = matrix(nrow=length(rho),ncol=length(gamma2),dimnames=list(as.character(signif(rho,3)),as.character(signif(gamma2,3))))
   err.var.a=omega; err.var.a[]=NA; # get the matrix the right size

   for (i.g2 in 1:length(gamma2)) {
     for (i.r in 1:length(rho)) {
        w=omega[i.r,i.g2] = rho[i.r] - gamma2[i.g2]
        err.var.a[i.r,i.g2]=1-ha-(w + sqrt(gamma2[i.g2]))^2
     }
   }

   gamma = sqrt(gamma2)
   err.var.b = 1-hb-gamma2;

   # simulate
   ha=hb=.5
   set.seed(1);
#   Nchoices=20
#   for(N in Nchoices) {

   
#   for (g2 in gamma2) {
g2=gamma2[1]
     g = sqrt(g2)
     errb = 1-hb-g2
#     for (cor.ab in rho) {
cor.ab = rho[1]
        w= cor.ab - g2;
        erra = 1-ha-w^2-g^2-2*(g*g*w)

        param.string=paste(sep="","ha=",ha, " hb=",hb, " g2=",g2, " erra=",erra," errb",errb, " cor.ab=",cor.ab, " w=",w)
  
        # re-calibrate Chi-squared study:
        leo.i=c()
        leo.o=c()
        leo.r=c()
  ab=c()
        ba=c()
        unr=c() 
        hid=c()
        ntrack=c()

        chi2.ba=c()
        chi2.ab=c()
        chi2.unr=c()
        chi2.hid=c()

        numreps=200
        for (rangen in c(20,50,80,100,200,500,800,1000,2000)) {
           N=rangen
        for(reps in 1:numreps) {
        SNPA=rnorm(N);
        SNPB=rnorm(N);
        varSNPa=var(SNPA);
        varSNPb=var(SNPB);
        C=rnorm(N);
        B= (sqrt(hb)/sqrt(varSNPb)) * SNPB + rnorm(N,sd=sqrt(errb)) + g*C;
        A= (sqrt(ha)/sqrt(varSNPa)) * SNPA + rnorm(N,sd=sqrt(erra)) + g*C + w*B;

        x=data.frame(SNPA,SNPB,A,B);

        z=local.sem.four.var.m1m2(pm=pm,MA.col=1,A.col=3,B.col=4,x=x,MB.col=2,no.obs=N,fit.models=c(1,2,3,4));
        leo.i = c(leo.i,10^(z$mlogp.M.M1M2hidden.con - z$mlogp.M.M1M2.AtoB))
        leo.r  =  c(leo.r,10^(z$mlogp.M.M1M2unresolv - z$mlogp.M.M1M2.AtoB))
        leo.o = c(leo.o,10^(z$mlogp.M.M1M2.BtoA - z$mlogp.M.M1M2.AtoB))

        ab=c(ab,z$mlogp.M.M1M2.AtoB)
        ba=c(ba,z$mlogp.M.M1M2.BtoA)
        unr=c(unr,z$mlogp.M.M1M2unresolv)
        hid=c(hid,z$mlogp.M.M1M2hidden.con)
        ntrack=c(ntrack,N)
        
        chi2.ba=c(chi2.ba,z$M.M1M2.BtoA$chisq)
        chi2.ab=c(chi2.ab,z$M.M1M2.AtoB$chisq)
        chi2.unr=c(chi2.unr,z$M.M1M2unresolv$chisq)
        chi2.hid=c(chi2.hid,z$M.M1M2hidden.con$chisq)
        
        } # end reps
      } # end rangen

     save(leo.i,leo.o,leo.r,ab,ba,unr,hid,ntrack,chi2.ba,chi2.ab,chi2.unr,chi2.hid,file="plus.chi2.power.problem.through.with.hid.2K.study.rdat")

   # plot the raw scores...
   windows()
   rawf=-log10(d$chi2.ba/(d$ntrack-1))
   boxplot(split(rawf,d$ntrack),notch=T,xlab="N.observations",ylab="-log10(minimum Fit function)",main=paste(sep="","Fitting the correct model MB->B->A<-MA\n",param.string))

   # get the medians, fit a polynomial to them...
   split.raw=split(rawf,d$ntrack)
   N.str=names(split.raw)
   meds=vector(length=length(split.raw))
   for (i in 1:length(split.raw)) {
     meds[i]=median(split.raw[[i]])
   }

#polyroot(base)                      Find Zeros of a Real or Complex Polynomial
#locpoly(KernSmooth)                 Estimate Functions Using Local Polynomials
#fit.trend(sgeostat)                 Fit polynomial trend functions
#bs(splines)                         B-Spline Basis for Polynomial Splines
#polySpline(splines)                 Piecewise Polynomial Spline Representation
#loess(stats)                        Local Polynomial Regression Fitting
#poly(stats)                         Compute Orthogonal Polynomials
   
  
   windows()
   boxplot(split(d$chi2.ba,d$ntrack),notch=TRUE,xlab="N.observations",ylab="chisq statistics")
   title("Chi-square curve when Correct model is fit\nModel MB->B->A<-MA, Simulation params: h.a = h.b = .5\ncor(A,B)=.45, gamma^2=0.05 erra=0.25, errb=0.45")   
  
  windows()
     d=data.frame(ab,ba,unr,ntrack,leo.r,leo.o,chi2.ba,chi2.ab,chi2.unr)
     boxplot(split(d$ba,d$ntrack),notch=TRUE,xlab="N.observations",ylab="-log10(Prob(X^2>chisq))")
     title("Diagnosing the chi^2 fitting method Power Problem:\nthe fit of the correct model actually declines with N.\n50 trials, MA->A<-B<-B")
     title(paste(sep="","In ",numreps," trials, mlogp for MA->A<-B<-MB, as N varies"))

     windows()
     boxplot(split(d$ab,d$ntrack),notch=TRUE,xlab="N.observations",ylab="-log10(Prob(X^2>chisq))")
     title("Fit of the wrong model (MA->A->B<-MB) does worsen faster\n50 trials")

     windows()
     boxplot(split(d$unr,d$ntrack),notch=TRUE,xlab="N.observations",ylab="-log10(Prob(X^2>chisq))")
     title("Fit of the unresolvable model (MA->A->B<-MB) does worsen faster\n50 trials")

     windows()
     boxplot(split(d$hid,d$ntrack),notch=TRUE,xlab="N.observations",ylab="-log10(Prob(X^2>chisq))")
     title("Fit of the hidden confounder model (MA->A->B<-MB) does worsen faster\n50 trials")

#     } # end cor.ab
#  } # end g2
#} # end N in Nchoices

  
}  # end check.variance.algebra



if(exists("omega.gamma.model") ) rm(omega.gamma.model);
omega.gamma.model=function(sim.start=NULL) {

   if (is.null(sim.start)) {
      sim.num <<- 6 # increment each time
   } else  { sim.num <<- sim.start; }
   
   # start the simulation
   ha.start =.5 # direct broad-sense heritability of A: the percentage of variance of A that is due to the genetic markers directly influence A, and not via B.
   hb.start =.5 # the broad-sense heritability of B: the percentage of variance of B that is due to the genetic loci influencing B.
   for (minus.ha in c(0,.2,.4)) {
   for (minus.hb in c(0,.2,.4)) {

   ha = ha.start - minus.ha
   hb = hb.start - minus.hb 
     
   npoints = 3 # number of points at which to evaluate gamma^2
   min.err.var.b=0.05
   min.gamma2=.05
   max.gamma2=1-hb-min.err.var.b
   gamma2 = seq(min.gamma2,max.gamma2,(max.gamma2-min.gamma2)/npoints);
   rho = rev(gamma2);
   omega = matrix(nrow=length(rho),ncol=length(gamma2),dimnames=list(as.character(signif(rho,3)),as.character(signif(gamma2,3))))
   err.var.a=omega; err.var.a[]=NA; # get the matrix the right size

   for (i.g2 in 1:length(gamma2)) {
     for (i.r in 1:length(rho)) {
        w=omega[i.r,i.g2] = rho[i.r] - gamma2[i.g2]
        err.var.a[i.r,i.g2]=1-ha-(w + sqrt(gamma2[i.g2]))^2
     }
   }

   gamma = sqrt(gamma2)
   err.var.b = 1-hb-gamma2;

   interc=c()
   coef=c()
   sim.param=c()

   # simulate
   set.seed(1);
#   Nchoices=20
#   for(N in Nchoices) {

   
   for (g2 in gamma2) {
#g2=gamma2[1]
     g = sqrt(g2)
     errb = 1-hb-g2
     for (cor.ab in rho) {
#cor.ab = rho[1]
        w= cor.ab - g2;
        erra = 1-ha-w^2-g^2-2*(g*g*w)

        param.string=paste(sep="","ha=",ha, " hb=",hb, " g2=",g2, " erra=",erra," errb=",errb, " cor.ab=",cor.ab, " w=",w)
        param=list()
        param$ha=ha
        param$hb=hb
        param$g2=g2
        param$erra=erra
        param$errb=errb
        param$cor.ab=cor.ab
        param$w=w
  
        # re-calibrate Chi-squared study:
        leo.i=c()
        leo.o=c()
        leo.r=c()
  ab=c()
        ba=c()
        unr=c() 
        hid=c()
        ntrack=c()

        chi2.ba=c()
        chi2.ab=c()
        chi2.unr=c()
        chi2.hid=c()

        numreps=50
        for (rangen in c(20,50,80,100,200,500,800,1000,2000)) {
           N=rangen
        for(reps in 1:numreps) {
        SNPA=rnorm(N);
        SNPB=rnorm(N);
        varSNPa=var(SNPA);
        varSNPb=var(SNPB);
        C=rnorm(N);
        B= (sqrt(hb)/sqrt(varSNPb)) * SNPB + rnorm(N,sd=sqrt(errb)) + g*C;
        A= (sqrt(ha)/sqrt(varSNPa)) * SNPA + rnorm(N,sd=sqrt(erra)) + g*C + w*B;

        x=data.frame(SNPA,SNPB,A,B);

        z=local.sem.four.var.m1m2(pm=pm,MA.col=1,A.col=3,B.col=4,x=x,MB.col=2,no.obs=N,fit.models=c(1,2,3,4));
        leo.i = c(leo.i,10^(z$mlogp.M.M1M2hidden.con - z$mlogp.M.M1M2.BtoA))
        leo.r  =  c(leo.r,10^(z$mlogp.M.M1M2unresolv - z$mlogp.M.M1M2.BtoA))
        leo.o = c(leo.o,10^(z$mlogp.M.M1M2.AtoB - z$mlogp.M.M1M2.BtoA))

        ab=c(ab,z$mlogp.M.M1M2.AtoB)
        ba=c(ba,z$mlogp.M.M1M2.BtoA)
        unr=c(unr,z$mlogp.M.M1M2unresolv)
        hid=c(hid,z$mlogp.M.M1M2hidden.con)
        ntrack=c(ntrack,N)
        
        chi2.ba=c(chi2.ba,z$M.M1M2.BtoA$chisq)
        chi2.ab=c(chi2.ab,z$M.M1M2.AtoB$chisq)
        chi2.unr=c(chi2.unr,z$M.M1M2unresolv$chisq)
        chi2.hid=c(chi2.hid,z$M.M1M2hidden.con$chisq)
        
        } # end reps
      } # end rangen

   sim.num = sim.num +1

   n1=sort(unique(ntrack))
   
   just.f = split(chi2.ba/(ntrack-1),ntrack)
   just.chi2ba = split(chi2.ba,ntrack)

   f.means=vector(length=length(just.f))
   j.chi2ba=vector(length=length(just.f))

   for (i in 1:length(just.f)) {
     f.means[i]=mean(just.f[[i]])
     j.chi2ba[i]=mean(just.chi2ba[[i]])
   }

   my.f=j.chi2ba[4]/f.means
   postscript(file=paste(sep="","fit.multiplier.sim.",sim.num,".ps"))
   plot(sqrt(n1),my.f,xlab="sqrt(N)",ylab="Fit multiplier value to obtain corrected chi-squared")
   m1=lm(my.f~sqrt(n1))
   abline(m1$coefficients)
   title(paste(m1$coefficients,collapse=" "))
   dev.off()

   # store parameters in file too
   cat(param.string,file=paste(sep="","key.sim.num.",sim.num,".txt"))

   chi.corrector = (10/sqrt(ntrack-1))

   d=data.frame(ab,ba,unr,hid,ntrack,leo.r,leo.o,chi2.ba,chi2.ab,chi2.unr,chi2.hid,leo.i,chi.corrector)

   # recompute LEO.I and LEO.O, using new chisquared...

   chi2.ba.correct = d$chi2.ba * d$chi.corrector
   ba.correct=-pchisq(chi2.ba.correct,3,lower.tail=FALSE,log.p=TRUE)*log.to.log10

   chi2.ab.correct = d$chi2.ab * d$chi.corrector
   ab.correct=-pchisq(chi2.ab.correct,3,lower.tail=FALSE,log.p=TRUE)*log.to.log10

   chi2.hid.correct = d$chi2.hid * d$chi.corrector
   hid.correct=-pchisq(chi2.hid.correct,3,lower.tail=FALSE,log.p=TRUE)*log.to.log10

   leo.i=10^(d$hid - d$ba)
   leo.i.correct=10^(hid.correct - ba.correct)

   leo.o=10^(d$ab - d$ba)
   leo.o.correct=10^(ab.correct - ba.correct)

   d=data.frame(ab,ba,unr,hid,ntrack,leo.r,leo.o,chi2.ba,chi2.ab,chi2.unr,chi2.hid,leo.i,chi.corrector,leo.i.correct,leo.o.correct)
   attr(d,"params")=param.string
   save(d,file=paste(sep="","scan.power.simulate.with.correction.",sim.num,".rdat"))        
      } #cor.ab
   } # g2
 } # minus.hb
 } # minus.ha

 return(0);
   
   # end main simulation loops

   for (sim in 7:26) {
      file.str=paste(sep="","sims/scan.power.simulate.with.correction.",sim,".rdat")
      file.out=paste(sep="","sims/scan.power.simulate.with.correction.",sim,".")
      load(file.str)

      all=1:length(d$chi2.ba)
      param.string = attr(d,"params")

   # LEO.I corrected versus uncorrected.
#   postscript(file=paste(sep="",file.out,"LEO.I.correct.compare.ps"),horizontal=FALSE);
   pdf(file=paste(sep="",file.out,"LEO.I.correct.compare.pdf"));
   boxplot(split(log10(d$leo.i.correct[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",at = 1:length(unique(d$ntrack[all])) - 0.2,boxwex = 0.25,col = "yellow")#,ylim=c(-20,20))
   boxplot(split(log10(d$leo.i[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",at = 1:length(unique(d$ntrack[all])) + 0.2,add=TRUE,boxwex = 0.25,col = "orange",show.names=FALSE)#,ylim=c(-20,20))
  legend(0.9, 17, c("Corrected", "Uncorrected"), fill = c("yellow", "orange"))
  abline(a=1,b=0,lty=2)
  abline(a=0,b=0,lty=3)  
  title(paste(sep="","True postive detection above the dashed line at log10(LEO.I)=1\n",param.string)) 
  dev.off();

    }

   all=1:450
   all=1:length(d$chi2.ba)

   postscript(file=paste(sep="","advan.corrected.sim.num.",sim.num,".ps"))
   boxplot(split(d$chi2.ba[all] * d$chi.corrector[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="10 sqrt(N) corrected chisq statistic")
   dev.off()

   postscript(file=paste(sep="","comparison.corrected.vs.uncorrected.",sim.num,".ps"))
   boxplot(split(d$chi2.ba[all] * d$chi.corrector[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="Chisq statistic",at = 1:length(n1) - 0.2,boxwex = 0.25,col = "yellow")
   boxplot(split(d$chi2.ba[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="Chisq statistic",at = 1:length(n1) + 0.2,add=TRUE,boxwex = 0.25,col = "orange",show.names=FALSE)

  legend(2.4, 30, c("Corrected", "Uncorrected"), fill = c("yellow", "orange"))
  dev.off()

  windows()
   boxplot(split(log10(d$leo.i.correct[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",col = "yellow")
  abline(a=1,b=0,lty=2)
  abline(a=0,b=0,lty=3)  
  title(paste(sep="","True postive detection above the dashed line at log10(LEO.I)=1\n",param.string))

    windows()
   boxplot(split(log10(d$leo.o.correct[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.O)",col = "yellow")
  abline(a=log10(.01),b=0,lty=2)
#  abline(a=0,b=0,lty=3)  
#  title(paste(sep="","True postive detection above the dashed line at log10(LEO.O)=1\n",param.string))
  title(paste(sep="","True postive detection above the dashed line at LEO.0 > 2\n",param.string))

   # LEO.I corrected versus uncorrected.
   windows()
   boxplot(split(log10(d$leo.i.correct[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",at = 1:length(unique(d$ntrack[all])) - 0.2,boxwex = 0.25,col = "yellow")#,ylim=c(-20,20))
   boxplot(split(log10(d$leo.i[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",at = 1:length(unique(d$ntrack[all])) + 0.2,add=TRUE,boxwex = 0.25,col = "orange",show.names=FALSE)#,ylim=c(-20,20))
  legend(1.4, 17, c("Corrected", "Uncorrected"), fill = c("yellow", "orange"))
  abline(a=1,b=0,lty=2)
  abline(a=0,b=0,lty=3)  
  title(paste(sep="","True postive detection above the dashed line at log10(LEO.I)=1\n",param.string)) 

   windows()
   boxplot(split(-log10(d$leo.i[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="-log10(LEO.I)",col = "orange")
   boxplot(split(d$leo.o[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="LEO.O")

   windows()
   boxplot(split(hid[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="hid")

   windows()
   boxplot(split(d$ba[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="ba")

   windows()
   boxplot(split(d$ab[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="ab")

   windows()
   boxplot(split(d$chi2.ab[all] * d$chi.corrector[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="corrected chi2 ab")

   windows()
   boxplot(split(d$chi2.hid[all] * d$chi.corrector[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="corrected chi2 hid")

   windows()
   boxplot(split(d$chi2.unr[all] * d$chi.corrector[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="corrected chi2 unr")

   windows()
   boxplot(split(d$chi2.ab[all] * d$chi.corrector[all],d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="corrected chi2 ab")



   sim.param=c(sim.param,param.string)

#   a=lm(inv.meds ~ log10(n1))
   a=lm(1/f.means ~ log10(n1))

   postscript(paste(sep="","fit.function.correction.",sim.num,".ps"))
   plot(log10(n1),inv.means,xlab="log10(N)",ylab="1/(mean F*)",main=paste("Fit function correction\n",param.string))
   lines(log10(n1),a$fitted.values)
   dev.off()

   postscript(paste(sep="","no.title.fit.function.correction.",sim.num,".ps"))
   plot(log10(n1),inv.means,xlab="log10(N)",ylab="1/(mean F*)")
   lines(log10(n1),a$fitted.values)
   dev.off()

   interc=c(interc,a$coefficients[1])
   coef=c(coef,a$coefficients[2])

   # practice correction...
   

#polyroot(base)                      Find Zeros of a Real or Complex Polynomial
#locpoly(KernSmooth)                 Estimate Functions Using Local Polynomials
#fit.trend(sgeostat)                 Fit polynomial trend functions
#bs(splines)                         B-Spline Basis for Polynomial Splines
#polySpline(splines)                 Piecewise Polynomial Spline Representation
#loess(stats)                        Local Polynomial Regression Fitting
#poly(stats)                         Compute Orthogonal Polynomials

   
      windows()
      boxplot(split(d$chi2.ba,d$ntrack),notch=TRUE,xlab="N.observations",ylab="chisq statistics")
      title("Chi-square curve when Correct model is fit\nModel MB->B->A<-MA, Simulation params: h.a = h.b = .5\ncor(A,B)=.45, gamma^2=0.05 erra=0.25, errb=0.45")   
  
     windows()
     boxplot(split(d$ba,d$ntrack),notch=TRUE,xlab="N.observations",ylab="-log10(Prob(X^2>chisq))")
     title("Diagnosing the chi^2 fitting method Power Problem:\nthe fit of the correct model actually declines with N.\n50 trials, MA->A<-B<-B")
     title(paste(sep="","In ",numreps," trials, mlogp for MA->A<-B<-MB, as N varies"))

     windows()
     boxplot(split(d$ab,d$ntrack),notch=TRUE,xlab="N.observations",ylab="-log10(Prob(X^2>chisq))")
     title("Fit of the wrong model (MA->A->B<-MB) does worsen faster\n50 trials")

     windows()
     boxplot(split(d$unr,d$ntrack),notch=TRUE,xlab="N.observations",ylab="-log10(Prob(X^2>chisq))")
     title("Fit of the unresolvable model (MA->A->B<-MB) does worsen faster\n50 trials")

     windows()
     boxplot(split(d$hid,d$ntrack),notch=TRUE,xlab="N.observations",ylab="-log10(Prob(X^2>chisq))")
     title("Fit of the hidden confounder model (MA->A->B<-MB) does worsen faster\n50 trials")

#     } # end cor.ab
#  } # end g2
#} # end N in Nchoices

# barplots:

 library(gplots)
 # Example with confidence intervals and grid
 hh <- t(VADeaths)[, 5:1]
 mybarcol <- "gray20"
 ci.l <- hh * 0.85
 ci.u <- hh * 1.15
 mp <- barplot2(hh, beside = TRUE,
        col = c("lightblue", "mistyrose",
                "lightcyan", "lavender"),
        legend = colnames(VADeaths), ylim = c(0, 100),
        main = "Death Rates in Virginia", font.main = 4,
        sub = "Faked 95 percent error bars", col.sub = mybarcol,
        cex.names = 1.5, plot.ci = TRUE, ci.l = ci.l, ci.u = ci.u,
        plot.grid = TRUE)
 mtext(side = 1, at = colMeans(mp), line = 2,
      text = paste("Mean", formatC(colMeans(hh))), col = "red")
 box()

  
}  # end omega.gamma.model


library(gplots)

#
# make frequency barplot of th.col, height of the bar
#  is the frequency with which th.col is below threshold
#  group.col specifies the column which indicates the different groups.
#
# 
# ex:
# freq.bar(ext.scores,threshold=.1,th.col="rmsea.f2.M.M1M2.BtoA",group.col="percent.max.omega")
#
if(exists("freq.bar") ) rm(freq.bar);
freq.bar=function(dataIN, threshold, th.col, group.col, plot.freq.below.th=TRUE,...) {

   if (is.character(th.col)) {  th.col=which(colnames(dataIN)==th.col); }
   if (is.character(group.col)) {  group.col=which(colnames(dataIN)==group.col); }

   if (length(th.col) == 0) { stop("freq.bar error: th.col not found") }
   if (length(group.col) == 0) { stop("freq.bar error: group.col not found")}

   nm = colnames(dataIN)[c(group.col, th.col)]
   
   # Hint: choose the threshold so that the POWER for omega=0 should is smaller than 0.01.

   # assign a 1 if the rmsea is below threshold: make an indicator whose mean is the frequency of cases below the Threshold.
   # if threshold is Inf, then everything gets a 1 so we just are generating group means
   LEO.Icut=ifelse(dataIN[,th.col]< threshold,1,0)
   my.y.lab=paste("Frequency",nm[2],"<",threshold)

   if (!plot.freq.below.th) {
     LEO.Icut=ifelse(dataIN[,th.col]>threshold,1,0)
     my.y.lab=paste("Frequency",nm[2],">",threshold)
   }

   # group data by percent.max.omega, then get the frequency of rmsea being good in each group
   # ; the mean or expected value of an indicator is just the probability or frequency of that
   # indicator.

   # tapply automatically sorts the percent.max.omega groups into ascending order.
   POWER=as.vector(tapply(LEO.Icut,dataIN[,group.col],mean));
   group.names=round(sort(unique(dataIN[,group.col])),2)
   names(POWER)=as.character(group.names)
   
   SE.POWER= as.vector(tapply(LEO.Icut,dataIN[,group.col],stderr1))

   ci.l = POWER-1.96*SE.POWER
   ci.u = POWER+1.96*SE.POWER

   windows()
   par(mfrow=c(1,1))

   #barplot(POWER,names.arg=names(table(dataIN$percent.max.omega) ), ylab="POWER",main=paste("No. Samples=",50, "RMSEA Threshold", threshold))
   #err.bp(as.vector(POWER), as.vector(1.96*SE.POWER), two.side=T)

   mp <- barplot2(POWER,
             xpd=FALSE, # bars stay inside region
             beside = TRUE,
        col = "yellow",
#        col = c("yellow", "violet"),
#                "lightcyan", "orange"),
#        legend = 
#                  ylim = c(0, 1),
#        main = "Death Rates in Virginia",
        font.main = 4,
        sub = nm[1],
        ylab= my.y.lab,
#        col.sub = "yellow",
        cex.names = 1.25,
        plot.ci = TRUE,
               ci.l = ci.l,
               ci.u = ci.u,
        plot.grid = TRUE,...)

  invisible(mp)
}


if(exists("false.pos.pub") ) rm(false.pos.pub);
false.pos.pub=function() {

   s2=sqrt(2)

   # simulate
   no.samples=100
   set.seed(1)
   gen.var=1;
   env.var=1;

   d1=date()

   marker.pattern=list(c(1),c(1,1),c(1,1,1))
   g.to.e = c(10,1,.1)
   effect.coef=c(2,1,.5,.25,.1);
   n.choice=c(50,100,200);

   lmp=length(marker.pattern);
   lec=length(effect.coef);
   lge=length(g.to.e);
   lnc=length(n.choice);

   res=c();

   lod.mean=array(data = 0, dim = c(lec,lge,lnc), dimnames = list(as.character(effect.coef),as.character(g.to.e),as.character(n.choice)));
   lod.sd=lod.mean
   eo.mean=lod.mean
   eo.sd =lod.mean   

   # different true models simulated...
   models.to.do=3
   for(i in models.to.do) {

   for (n in 1:lnc) {
     no.samples=n.choice[n]

   for (g in 1:lge) {
    env.var = 1/g.to.e[g]

   for (e in 1:lec) {
    ABeff=effect.coef[e]
    MAeff=1; #effect.coef[e] #, otherwise get a false reading on geneic:env variance
    MBeff=1; #effect.coef[e]

   num.sim=10
   model.num=1
#   disc.matrix=matrix(nrow=6,ncol=6,data=0);
   eo.m = vector(length=num.sim);
   eo.m[]=0;
      no.obs.Z=no.samples
      for (k in 1:num.sim) {

         mod=sim.model.num(model.num=i,env.var,no.samples,MAeff,ABeff,MBeff);

         covx=cov(mod)

            # twofer.eo and n.twofer.eo hold our results...
            z=zeo2(1,3,2,covx,pm); # z$eo

            #cor1=pcor(c(1,3,c()),covx)
            #Zm1B= sqrt(no.obs.Z-pm$fisher.dof.cor)*   fisher1(abs(cor1))

            #cor2=pcor(c(1,3,2),covx)
            #Zm1BgivenA= sqrt(no.obs.Z-pm$fisher.dof.pcor1)*   fisher1(abs(cor2))

            #s=1;
            #eo2 = Zm1B - s*Zm1BgivenA;
            eo2 = z$eo;
         eo.m[k]=eo2;
     } # k

        eo.mean[e,g,n]=mean(eo.m)
        eo.sd[e,g,n]=sd(eo.m)

        cat(paste(e,g,n,date(),"\n"))
     } # e
     } # g
     } # n

   d2=date()
  z=list()

  z$eo.mean=eo.mean
  z$eo.sd=eo.sd

  z$d1=d1
  z$d2=d2

  z$eo.m=eo.m
  save(z,file=paste(sep="","power.study3zeo.i.is.",i,".reps.",num.sim,".rdat"))    
  res=c(z,res);

  } # end i over 6 models

  res
}

#
# compare SNP aggregation strategies with a fixed A and B
#
if(exists("aggregate.snps.fixed.A.B") ) rm(aggregate.snps.fixed.A.B);
aggregate.snps.fixed.A.B=function(datCombined,MA.snpcols,MB.snpcols,A.col, B.col, pm=neo.get.param())
{
  r=list() # results

  # condense matrix down
  x = datCombined[,c(MA.snpcols,MB.snpcols,A.col,B.col)]

  len.ma = length(MA.snpcols)
  len.mb = length(MB.snpcols)
  len.a  = length(A.col)
  len.b  = length(B.col)
  
  # adjust column numbers to reflect assignments in reduced x matrix
  MA.snpcols = 1:len.ma
  MB.snpcols = len.ma+ (1:len.mb)
  A.col = max(MB.snpcols) + (1:len.a)
  B.col = max(A.col) + (1:len.b)

  snpcols = c(MA.snpcols, MB.snpcols)

  covx=cov(x,use="pairwise.complete.obs")
  cx=abs(cov2cor(covx))

  pm$quiet=TRUE

  no.obs = nrow(x)
  cn=colnames(x)
  af=cx; af[]=0

  # set up the graph to analyze / aggregate SNPs on.
  af[MA.snpcols,A.col]=1
  af[MB.snpcols,B.col]=1
  af[A.col,B.col]=1  

  non.snp.cols=c(A.col,B.col) # very important this order, A.col first

  twofer.eo=af; twofer.eo[]=0
  twofer.log=matrix(data="",nrow=nrow(af),ncol=ncol(af))
  n.twofer.eo=matrix(NA,nrow=nrow(af),ncol=ncol(af))

  r$zeo.mylist.maxmax = walk.two.steps.max.max(pm,snpcols,non.snp.cols,af,no.obs,cx,covx,twofer.eo,n.twofer.eo,cn,twofer.log);
  r$zeo.mylist.mlreg = walk.two.steps.weighted.mean.mlreg(x,pm,snpcols,non.snp.cols,af,no.obs,cx,covx,twofer.eo,n.twofer.eo,cn,twofer.log);
  r$zeo.mylist.allvote = walk.two.steps.one.permutation.one.vote(x,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,n.twofer.eo,cn,twofer.log);

   r$leo.mylist.maxmax = losem.walk.two.steps.max.max(x,pm,snpcols,non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log);
   r$leo.mylist.mlreg = losem.walk.two.steps.weighted.mean.mlreg(x,pm,snpcols,non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log);
   r$leo.mylist.allvote = losem.walk.two.steps.one.permutation.one.vote(x,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log);



   invisible(r)
}


# permute.sma():
#    Do single.marker.analysis repeatedly on permutations of the data
#
# example call
#   permute.sma(datCombined=datCombined,block.to.permute=str.me$snpcols[1:5],num.perm=3,snpcols=str.me$snpcols[1:5],genecols=which.Insig1, traitcols=c(which.Fdft1,which.Dhcr7))
#
#    permute.sma(datCombined=dc,block.to.permute=1,num.perm=10,snpcols=1,genecols=2,traitcols=c(3))
#
# Since it calls permute.snps.traits(), then pm$block2.to.permute and pm$block3.to.permute
#  can be also used to specify addition blocks of columns to permute simultaneously.
#
# SunGridEngine usage: If the pm$ parameters are set correctly, then
# permutation jobs will be qsub-ed to SGE. Set pm$run.perms.on.sge = TRUE
# to accomplish this.
#
# Also, note that if num.perm overrides pm$number.BLOCK.permutations.
#
# Note that
# you should be running in a unique directory for this run, as in
# order to share data, the pre-permutation data is written to
# a well known (commonly named) file. If multiple runs share a
# directory, then this could get overwritten if the jobs overlap.
#
# Note that pm$sge.needs.neo.install.path, and
#           pm$sge.needs.r.program.install.path, must be set correctly for NEO/SGE to work.
#
# pm$just.do.permutation.no.unpermuted is honored, if pm$run.perms.on.sge == TRUE.
#
if(exists("permute.sma") ) rm(permute.sma);
permute.sma=function(datCombined,snpcols,genecols,traitcols,leo.i.th=1, leo.o.th=5, leo.nb.th=.3,pm=neo.get.param(),use.ranks=FALSE,build.multi.marker.to.gene.model=FALSE, impute.na=TRUE, num.perm=2, sma.perm.file="sma.perm.raw.csv", block.to.permute) {


   one.date=paste(sep="",gsub(":",".",gsub(" ","_",paste(sep="",date()))))
   print(paste(sep="","Starting permute.sma at ",one.date ))

   # use supplied parameter instead of pm version...but set pm$number.BLOCK.permutations so re-used code below works.
   if (pm$number.BLOCK.permutations!=0) {
      pm.print(pm,paste("over-riding pm$number.BLOCK.permutations with supplied num.perm = ",num.perm))
   }
   pm$number.BLOCK.permutations = num.perm

   if (any(!(block.to.permute %in% (1:ncol(datCombined))))) {
      stop("Bad specification of block.to.permute: out of range")
   }

   # auto-farm out permutations to SGE if requested
   if (pm$run.perms.on.sge & pm$number.BLOCK.permutations > 0) {
      # save the current environment
      # after changing pm$run.perms.on.sge to avoid infinite loop
      pm$run.perms.on.sge = FALSE

      # and give sge jobs their own directory, inside the main log directory
      pm$sge.run.dir = create.logpath.enclosed.sge.directory(pm)
      pm.print(pm,paste("Created directory for Sun Grid Engine job output: ",pm$sge.run.dir))      

      logfile.to.restore = pm$neo.log.file
      pm$neo.log.file = "sge_neo."

      # temporarily change directory into sge run directory
      cwd=getwd()
      setwd(pm$sge.run.dir)
      
      save(list=ls(), file="permute.sma.master.sge.image.rdat")
      pm$neo.log.file = logfile.to.restore

      
      # function neo.sge.perm loads the saved file "permute.sge.master.sge.image.rdat" and runs one permutation (could be more).
      #n.perm.job=ceiling(pm$number.BLOCK.permutations/10)

      n.perm.job=pm$number.BLOCK.permutations

      job.list=rep("sma.neo.sge.perm",n.perm.job)

      pm.print(pm,paste("Submitting ",n.perm.job," sma permutation jobs to Sun Grid Engine."))
      submit.to.sge(job.list=job.list,
                    neo.install.path=pm$sge.needs.neo.install.path,
                    r.program.install.path=pm$sge.needs.r.program.install.path)

      num.perm = pm$number.BLOCK.permutations = 0

      setwd(cwd) # restore directory, now that jobs are submitted.
      
      # terminate function call if we are halting early, only permutations requested.
      if (pm$just.do.permutation.no.unpermuted) return(invisible(c()))
    }

   sge.task = get.sge.task.or.zero()
   task.str=""
   if (sge.task > 0) task.str=paste(sep="",".",sge.task,".")

   # add timestamp (to allow re-runs without crashes) and .csv to the requested file
   sma.perm.file = paste(sep="",gsub(":",".",gsub(" ","_",paste(sep="",sma.perm.file,".",one.date))),".csv")

   
   # pre permutation, to compute p-values for observed!
   sma.observed= single.marker.analysis(datCombined=datCombined,snpcols=snpcols,genecols=genecols,traitcols=traitcols,leo.i.th=leo.i.th, leo.o.th=leo.o.th, leo.nb.th=leo.nb.th,pm=pm,use.ranks=use.ranks,build.multi.marker.to.gene.model=build.multi.marker.to.gene.model, impute.na=impute.na)

   sma.file = paste(sep="",gsub(":",".",gsub(" ","_",paste(sep="","sma.observed.unpermuted.",one.date))),".csv")
   write.csv(sma.observed,file=sma.file)

   print(paste("Done with unpermuted SMA analysis (see file",sma.file,")."))
   if (num.perm > 0) {
      print(paste(" Proceeding to do ",num.perm," permutations."))
   } else {
      return(invisible(c()))
   }

   obs.stat = sma.observed[,c(1,1,3,5,6,7,8,9,11,12:ncol(sma.observed))] 

      colnames(obs.stat)[1]="model"
      colnames(obs.stat)[2]="M"
      colnames(obs.stat)[3]="A"
      colnames(obs.stat)[4]="B"

      obs.stat[,2] = as.character(obs.stat[,2])
      obs.stat[,3] = as.character(obs.stat[,3])
      obs.stat[,4] = as.character(obs.stat[,4])

      # create unique row identifier
      obs.stat[,1]=paste(sep="->",obs.stat[,2],obs.stat[,3],obs.stat[,4])

   # Got the observed to compute p-values with...

   # Now start the permutations.

   # execute SMA
   for (i in 1:num.perm) {
      pm.print(pm,paste("single.marker.analysis permutations, doing : ",i," out of ",num.perm," at ",date()))

      datCombined=permute.snps.traits(pm=pm, datCombined, block.to.permute, save.perm.to.file=FALSE)

      sma.pre = single.marker.analysis(datCombined=datCombined,snpcols=snpcols,genecols=genecols,traitcols=traitcols,leo.i.th=leo.i.th, leo.o.th=leo.o.th, leo.nb.th=leo.nb.th,pm=pm,use.ranks=use.ranks,build.multi.marker.to.gene.model=build.multi.marker.to.gene.model, impute.na=impute.na)

           #old call: single.marker.analysis(datCombined=datCombined,snpcols=snpcols,pm=pm,...)

      # write out just statistics of interest...actually keep most everything now, adding :ncol
      sma = sma.pre[,c(1,1,3,5,6,7,8,9,11,12:ncol(sma.pre))] 
      colnames(sma)[1]="model"
      colnames(sma)[2]="M"
      colnames(sma)[3]="A"
      colnames(sma)[4]="B"

      sma[,2] = as.character(sma[,2])
      sma[,3] = as.character(sma[,3])
      sma[,4] = as.character(sma[,4])

      # create unique row identifier
      sma[,1]=paste(sep="->",sma[,2],sma[,3],sma[,4])

      if (i==1) {
         write.table(sma, file=sma.perm.file, col.names = NA, sep=",",dec=".",qmethod="double")
      } else {
         write.table(sma, file=sma.perm.file, append=TRUE, col.names = FALSE, sep=",",dec=".",qmethod="double")
      }
   } # end i loop

   # now do a variation on/modify create.permutation.report() to generate final permutation statistics.

   # these three variables accessed by create.permutation.report
   pm$perm.file = sma.perm.file
   pm$neo.log.file = gsub(":",".",gsub(" ","_",paste(sep="","sma.perms.",one.date)))
   pm$open.excel.at.end = TRUE

   #
   #
   # Adapt: create.permutation.report(pm,datCombined,snpcols) 
   #
   #

   # first 4 columns are row label stuff, skip to the 5th column for statistics.
   skipto=5

   perms = read.csv(file=pm$perm.file) # let rownames become their own column, since they will be repeated
                  # and non-unique rownames might be a problem.

   perms = perms[,-1] # get rid of 1,2,3,...,1,2,3,...,1,2,3,... and let edge be first column
   perms[,1]=as.character(perms[,1]) # avoid factor messes

   # and summarize
   aur=unique(perms[,1]) # first column has edge names
   num.perm.summary.rows = length(aur)

   perm.mean = perms[1:num.perm.summary.rows,]
   perm.mean[,skipto:ncol(perms)] = NA
   perm.mean[,1] = aur

   report = create.permutation.report.sma(pm,colnames(datCombined),snpcols,perms,skipto=5,num.stats=5,skipable=c(1:4),p.file=paste(sep="",pm$neo.log.file,".collated.sma.permutation.results.csv"), aur, perm.mean, obs.stat.use = obs.stat)

   print(paste(sep="","Done with permute.sma at ",date() ))

   invisible(report)
}


if (exists("create.permutation.report.sma")) rm(create.permutation.report.sma)
create.permutation.report.sma=function(pm,cn.datCombined,snpcols,perms,skipto=5,num.stats=5,skipable=c(1:4),p.file=paste(sep="",pm$neo.log.file,".collated.sma.permutation.results.csv"), aur, perm.mean, obs.stat.use) {

   # copy and reuse perm.mean as a labelled starting point matrix...

   # summary of good (not na) data
   perm.percent.good = perm.mean
   perm.count.good = perm.mean

   # sd
   perm.sd = perm.mean

   # percentage of time permutated statistic exceeded observed, empirically observed.
   perm.pvalue = perm.mean

   # theoretical p-value, based on assuming the null distribution of the statistics is Gaussian.
   perm.theory.pvalue = perm.mean

   # shapiro.test() results for Shapiro-Wilk test of normality
   perm.shapiro = perm.mean

   # max, skew, kurtosis permutation statistics
   perm.max = perm.mean
   perm.min = perm.mean
   perm.median = perm.mean
   perm.skewness = perm.mean
   perm.kurtosis = perm.mean


   # This may or may not be right to do, but simplifies the reporting for positive statistics.
   # --> Turn NA in permutation data into zero.
   ## perms[is.na(perms)]=0

   for (i in 1:length(aur)) {
      cur.row = aur[i]
      w = which(cur.row == perms[,1]) # first column has edge names
      if (length(w)>1) {
         perm.mean[i,skipto:ncol(perms)] = mean(perms[w,skipto:ncol(perms)],na.rm=TRUE)
         perm.max[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(max,x,na.rm=T))
         perm.min[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(min,x,na.rm=T))

         perm.median[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(median,x,na.rm=T))

         perm.skewness[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(skewness,x,na.rm=T))
         perm.kurtosis[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(kurtosis,x,na.rm=T))

         perm.sd[i,skipto:ncol(perms)] = fail.proof.FUN(sd,perms[w,skipto:ncol(perms)],na.rm=TRUE)
         perm.count.good[i,skipto:ncol(perms)]=apply(!is.na.or.nan(perms[w,skipto:ncol(perms)]),2,sum)
         perm.percent.good[i,skipto:ncol(perms)]=apply(!is.na.or.nan(perms[w,skipto:ncol(perms)]),2,sum)/pm$number.BLOCK.permutations
         perm.shapiro[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,fail.proof.shapiro.test)

         #perm.pvalue computed depends on the observed statistic in that column:
         for (j in skipto:ncol(perms)) {
             perm.pvalue[i,j] = sum(perms[w,j]>= obs.stat.use[i,j],na.rm=TRUE)/pm$number.BLOCK.permutations
             perm.theory.pvalue[i,j] = pnorm(obs.stat.use[i,j],mean=perm.mean[i,j],sd=perm.sd[i,j],lower.tail=FALSE)
         }
      } else {
         # handle special case of just one permutation gracefully
         perm.mean[i,skipto:ncol(perms)] = perms[w,skipto:ncol(perms)]
         perm.sd[i,skipto:ncol(perms)] = Inf #

         perm.max[i,skipto:ncol(perms)] = (perms[w,skipto:ncol(perms)])
         perm.skewness[i,skipto:ncol(perms)] = NA
         perm.kurtosis[i,skipto:ncol(perms)] = NA

         perm.min[i,skipto:ncol(perms)] = (perms[w,skipto:ncol(perms)])
         perm.median[i,skipto:ncol(perms)] = (perms[w,skipto:ncol(perms)])

         perm.percent.good[i,skipto:ncol(perms)]=apply(!is.na.or.nan(perms[w,skipto:ncol(perms)]),2,sum)/pm$number.BLOCK.permutations
         perm.count.good[i,skipto:ncol(perms)]=apply(!is.na.or.nan(perms[w,skipto:ncol(perms)]),2,sum)

         perm.pvalue[i,skipto:ncol(perms)] = NA # can't compute p-value given just a single permutation.
         perm.theory.pvalue[i,skipto:ncol(perms)] = NA # can't compute p-value given just a single permutation.
         perm.shapiro[i,skipto:ncol(perms)] = NA
      }
   }

   #skipable=c(1:4)

   # interleaf the results to provide a summary including p-value
   colnames(perm.pvalue)[-skipable] = paste(sep="","pval.empirical.",colnames(perm.pvalue)[-skipable])
   colnames(perm.theory.pvalue)[-skipable] = paste(sep="","pval.theory.",colnames(perm.theory.pvalue)[-skipable])
   colnames(perm.sd)[-skipable] = paste(sep="","sd.",colnames(perm.sd)[-skipable])
   colnames(perm.mean)[-skipable] = paste(sep="","mean.",colnames(perm.mean)[-skipable])
   colnames(perm.median)[-skipable] = paste(sep="","median.",colnames(perm.median)[-skipable])
   colnames(perm.min)[-skipable] = paste(sep="","min.",colnames(perm.min)[-skipable])
   colnames(perm.max)[-skipable] = paste(sep="","max.",colnames(perm.max)[-skipable])
   colnames(perm.skewness)[-skipable] = paste(sep="","skewness.",colnames(perm.skewness)[-skipable])
   colnames(perm.kurtosis)[-skipable] = paste(sep="","kurtosis.",colnames(perm.kurtosis)[-skipable])
   colnames(perm.percent.good)[-skipable] = paste(sep="","percent.good.",colnames(perm.percent.good)[-skipable])
   colnames(perm.count.good)[-skipable] = paste(sep="","count.good.",colnames(perm.count.good)[-skipable])
   colnames(perm.shapiro)[-skipable] = paste(sep="","shapiro.pval.",colnames(perm.shapiro)[-skipable])

   #num.stats = 5
   space = data.frame(matrix(data="",nrow=nrow(perm.mean),ncol=num.stats))

   interleaf = data.frame(cbind(obs.stat.use,  perm.pvalue[,-skipable],  perm.theory.pvalue[,-skipable],  perm.mean[,-skipable],  perm.sd[,-skipable],  perm.percent.good[,-skipable],  perm.count.good[,-skipable],  perm.shapiro[,-skipable],  perm.max[,-skipable],  perm.skewness[,-skipable],  perm.kurtosis[,-skipable],  perm.min[,-skipable], perm.median[,-skipable], space))

   cn.il = colnames(interleaf)

   i = 0
   for (k in skipto:ncol(perms)) {
      i=i+1
      j = (k-skipto)*14+skipto
      interleaf[,j] = obs.stat.use[,k]
      interleaf[,j+1] = perm.pvalue[,k]
      interleaf[,j+2] = perm.theory.pvalue[,k]
      interleaf[,j+3] = perm.mean[,k]
      interleaf[,j+4] = perm.sd[,k]
      interleaf[,j+5] = perm.max[,k]
      interleaf[,j+6] = perm.min[,k]
      interleaf[,j+7] = perm.median[,k]
      interleaf[,j+8] = perm.skewness[,k]
      interleaf[,j+9] = perm.kurtosis[,k]
      interleaf[,j+10] = perm.percent.good[,k]
      interleaf[,j+11] = perm.count.good[,k]
      interleaf[,j+12] = perm.shapiro[,k]
      interleaf[,j+13] = " "

      cn.il[j]=colnames(obs.stat.use)[k]
      cn.il[j+1]=colnames(perm.pvalue)[k]
      cn.il[j+2]=colnames(perm.theory.pvalue)[k]
      cn.il[j+3]=colnames(perm.mean)[k]
      cn.il[j+4]=colnames(perm.sd)[k]
      cn.il[j+5]=colnames(perm.max)[k]
      cn.il[j+6]=colnames(perm.min)[k]
      cn.il[j+7]=colnames(perm.median)[k]
      cn.il[j+8]=colnames(perm.skewness)[k]
      cn.il[j+9]=colnames(perm.kurtosis)[k]
      cn.il[j+10]=colnames(perm.percent.good)[k]
      cn.il[j+11]=colnames(perm.count.good)[k]
      cn.il[j+12]=colnames(perm.shapiro)[k]
      cn.il[j+13]=paste(rep("_",i+3),collapse="")
   }

   colnames(interleaf) = cn.il
   #p.file=paste(sep="",pm$neo.log.file,".collated.sma.permutation.results.csv")
   write.csv(interleaf,file=p.file)
   if (pm$open.excel.at.end) open.excel.file(p.file,pm)

   invisible(interleaf)

} # end create.permutation.report2


# create.GROUP.permutation.report(): a variation to compute FDR over a whole group of SNPs tested
# in one permutation.
#
# column one of perms should be shared by rows to be grouped.
# The aur variable specifies all unique row groupings (and
# so length(aur) is the number of rows in the report.
#
# Currently only a single significance level, sig.th.to.eval is
# evaluated in a single run.
#
# sig.th.to.eval supplies the nominal significance levels to
# count how many tests in a given permutation exceeded this level.
# We can then average the count over all permutations to get
# and FDR for the set of tests run, given the significance level.

if (exists("create.GROUP.permutation.report")) rm(create.GROUP.permutation.report)
create.GROUP.permutation.report=function(pm,cn.datCombined,snpcols,perms,skipto=5,num.stats=5,skipable=c(1:4),p.file=paste(sep="",pm$neo.log.file,".GROUP.permutation.results.csv"), aur, perm.mean, obs.stat.use,sig.th.to.eval = 1) {

   # copy and reuse perm.mean as a labelled starting point matrix...

   # summary of good (not na) data
   perm.percent.good = perm.mean
   perm.count.good = perm.mean

   # number of scores exceeding the supplied threshold, sig.th.to.eval
   perm.count.exceed.th =perm.mean

   # sd
   perm.sd = perm.mean

   # percentage of time permutated statistic exceeded observed, empirically observed.
   perm.pvalue = perm.mean

   # theoretical p-value, based on assuming the null distribution of the statistics is Gaussian.
   perm.theory.pvalue = perm.mean

   # shapiro.test() results for Shapiro-Wilk test of normality
   perm.shapiro = perm.mean

   # max, skew, kurtosis permutation statistics
   perm.max = perm.mean
   perm.min = perm.mean
   perm.median = perm.mean
   perm.skewness = perm.mean
   perm.kurtosis = perm.mean

   for (i in 1:length(aur)) {
      cur.row = aur[i]
      w = which(cur.row == perms[,1]) # first column has edge names
      if (length(w)>1) {
         perm.mean[i,skipto:ncol(perms)] = mean(perms[w,skipto:ncol(perms)],na.rm=TRUE)
         perm.max[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(max,x,na.rm=T))
         perm.min[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(min,x,na.rm=T))

         perm.median[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(median,x,na.rm=T))

         perm.skewness[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(skewness,x,na.rm=T))
         perm.kurtosis[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,function(x) fail.proof.FUN(kurtosis,x,na.rm=T))

         perm.sd[i,skipto:ncol(perms)] = fail.proof.FUN(sd,perms[w,skipto:ncol(perms)],na.rm=TRUE)
         perm.count.good[i,skipto:ncol(perms)]=apply(!is.na.or.nan(perms[w,skipto:ncol(perms)]),2,sum)
         perm.percent.good[i,skipto:ncol(perms)]=apply(!is.na.or.nan(perms[w,skipto:ncol(perms)]),2,sum)/pm$number.BLOCK.permutations
         perm.shapiro[i,skipto:ncol(perms)] = apply(perms[w,skipto:ncol(perms)],2,fail.proof.shapiro.test)

         #perm.pvalue computed depends on the observed statistic in that column:
         for (j in skipto:ncol(perms)) {
             perm.pvalue[i,j] = sum(perms[w,j]>= sig.th.to.eval,na.rm=TRUE)/pm$number.BLOCK.permutations
             perm.theory.pvalue[i,j] = pnorm(sig.th.to.eval,mean=perm.mean[i,j],sd=perm.sd[i,j],lower.tail=FALSE)
             perm.count.exceed.th[i,j] = sum(perms[w,j]>= sig.th.to.eval,na.rm=TRUE)
         }
      } else {
         # handle special case of just one permutation gracefully
         perm.mean[i,skipto:ncol(perms)] = perms[w,skipto:ncol(perms)]
         perm.sd[i,skipto:ncol(perms)] = Inf #

         perm.max[i,skipto:ncol(perms)] = (perms[w,skipto:ncol(perms)])
         perm.skewness[i,skipto:ncol(perms)] = NA
         perm.kurtosis[i,skipto:ncol(perms)] = NA

         perm.min[i,skipto:ncol(perms)] = (perms[w,skipto:ncol(perms)])
         perm.median[i,skipto:ncol(perms)] = (perms[w,skipto:ncol(perms)])

         perm.percent.good[i,skipto:ncol(perms)]=apply(!is.na.or.nan(perms[w,skipto:ncol(perms)]),2,sum)/pm$number.BLOCK.permutations
         perm.count.good[i,skipto:ncol(perms)]=apply(!is.na.or.nan(perms[w,skipto:ncol(perms)]),2,sum)
         perm.count.exceed.th[i,skipto:ncol(perms)]=apply(perms[w,skipto:ncol(perms)]>sig.th.to.eval,2,sum)
         
         perm.pvalue[i,skipto:ncol(perms)] = NA # can't compute p-value given just a single permutation.
         perm.theory.pvalue[i,skipto:ncol(perms)] = NA # can't compute p-value given just a single permutation.
         perm.shapiro[i,skipto:ncol(perms)] = NA
      }
   }

   #skipable=c(1:4)

   # interleaf the results to provide a summary including p-value
   colnames(perm.count.exceed.th)[-skipable] = paste(sep="","count.exceed_",sig.th.to.eval,"_",colnames(perm.count.exceed.th)[-skipable])
   colnames(perm.pvalue)[-skipable] = paste(sep="","pval.empirical.",colnames(perm.pvalue)[-skipable])
   colnames(perm.theory.pvalue)[-skipable] = paste(sep="","pval.theory.",colnames(perm.theory.pvalue)[-skipable])
   colnames(perm.sd)[-skipable] = paste(sep="","sd.",colnames(perm.sd)[-skipable])
   colnames(perm.mean)[-skipable] = paste(sep="","mean.",colnames(perm.mean)[-skipable])
   colnames(perm.median)[-skipable] = paste(sep="","median.",colnames(perm.median)[-skipable])
   colnames(perm.min)[-skipable] = paste(sep="","min.",colnames(perm.min)[-skipable])
   colnames(perm.max)[-skipable] = paste(sep="","max.",colnames(perm.max)[-skipable])
   colnames(perm.skewness)[-skipable] = paste(sep="","skewness.",colnames(perm.skewness)[-skipable])
   colnames(perm.kurtosis)[-skipable] = paste(sep="","kurtosis.",colnames(perm.kurtosis)[-skipable])
   colnames(perm.percent.good)[-skipable] = paste(sep="","percent.good.",colnames(perm.percent.good)[-skipable])
   colnames(perm.count.good)[-skipable] = paste(sep="","count.good.",colnames(perm.count.good)[-skipable])
   colnames(perm.shapiro)[-skipable] = paste(sep="","shapiro.pval.",colnames(perm.shapiro)[-skipable])

   #num.stats = 5
   space = data.frame(matrix(data="",nrow=nrow(perm.mean),ncol=num.stats))

   interleaf = data.frame(cbind(perm.count.exceed.th,  perm.pvalue[,-skipable],  perm.theory.pvalue[,-skipable],  perm.mean[,-skipable],  perm.sd[,-skipable],  perm.percent.good[,-skipable],  perm.count.good[,-skipable],  perm.shapiro[,-skipable],  perm.max[,-skipable],  perm.skewness[,-skipable],  perm.kurtosis[,-skipable],  perm.min[,-skipable], perm.median[,-skipable], space))

   cn.il = colnames(interleaf)

   i = 0
   for (k in skipto:ncol(perms)) {
      i=i+1
      j = (k-skipto)*14+skipto
      interleaf[,j] =   perm.count.exceed.th[,k]
      interleaf[,j+1] = perm.pvalue[,k]
      interleaf[,j+2] = perm.theory.pvalue[,k]
      interleaf[,j+3] = perm.mean[,k]
      interleaf[,j+4] = perm.sd[,k]
      interleaf[,j+5] = perm.max[,k]
      interleaf[,j+6] = perm.min[,k]
      interleaf[,j+7] = perm.median[,k]
      interleaf[,j+8] = perm.skewness[,k]
      interleaf[,j+9] = perm.kurtosis[,k]
      interleaf[,j+10] = perm.percent.good[,k]
      interleaf[,j+11] = perm.count.good[,k]
      interleaf[,j+12] = perm.shapiro[,k]
      interleaf[,j+13] = " "

      cn.il[j]=colnames(perm.count.exceed.th)[k]
      cn.il[j+1]=colnames(perm.pvalue)[k]
      cn.il[j+2]=colnames(perm.theory.pvalue)[k]
      cn.il[j+3]=colnames(perm.mean)[k]
      cn.il[j+4]=colnames(perm.sd)[k]
      cn.il[j+5]=colnames(perm.max)[k]
      cn.il[j+6]=colnames(perm.min)[k]
      cn.il[j+7]=colnames(perm.median)[k]
      cn.il[j+8]=colnames(perm.skewness)[k]
      cn.il[j+9]=colnames(perm.kurtosis)[k]
      cn.il[j+10]=colnames(perm.percent.good)[k]
      cn.il[j+11]=colnames(perm.count.good)[k]
      cn.il[j+12]=colnames(perm.shapiro)[k]
      cn.il[j+13]=paste(rep("_",i+3),collapse="")
   }

   colnames(interleaf) = cn.il
   #p.file=paste(sep="",pm$neo.log.file,".collated.sma.permutation.results.csv")

# for now we actually just want perm.count.exceed.th, rather than interleaf--b/c all those
# other statistics aren't really relevant to the grouping...at the moment.
   
#   write.csv(interleaf,file=p.file)

   write.csv(perm.count.exceed.th,file=p.file)
   if (pm$open.excel.at.end) open.excel.file(p.file,pm)

   invisible(interleaf)

} # end create.GROUP.permutation.report


#returns a dataframe fitting the 3 variable models SNP -> gene -> trait
#for each of the supplied snpcols, genecols, and traitcols
#
if(exists("single.marker.analysis") ) rm(single.marker.analysis);
single.marker.analysis=function(datCombined,snpcols,genecols,traitcols,leo.i.th=1, leo.o.th=5, leo.nb.th=.3,pm=neo.get.param(),use.ranks=FALSE,build.multi.marker.to.gene.model=FALSE, impute.na=TRUE)
{
 # should we convert the gene and trait cols to ranks first?
 if (use.ranks) {
    for (g in genecols) {
         datCombined[,g]=rank(datCombined[,g])
    }
    for (g in traitcols) {
         datCombined[,g]=rank(datCombined[,g])
    }
 }

 cn=colnames(datCombined)

 tsize=length(traitcols)*length(snpcols)*length(genecols)

 M.col=1
 A.col=2
 B.col=3

 mlogp.M.AtoB=rep(NA,tsize)
 mlogp.M.BtoA=rep(NA,tsize)
 mlogp.M.conf=rep(NA,tsize)
 mlogp.M.AcollideB=rep(NA,tsize)
 mlogp.M.BcollideA=rep(NA,tsize)

 RMSEA.AtoB=rep(NA,tsize)
 RMSEA.BtoA=rep(NA,tsize)
 RMSEA.conf=rep(NA,tsize)
 RMSEA.AcollideB=rep(NA,tsize)
 RMSEA.BcollideA=rep(NA,tsize)


 #mlogp.M.AhiddenB=rep(NA,tsize)

 leo.i = rep(NA,tsize)
 leo.o = rep(NA,tsize)
 leo.ab.over.ba = rep(NA,tsize)
 good.o=rep(NA,tsize)
 good.i=rep(NA,tsize)
 good.nb=rep(NA,tsize)

   BAc.vs.NextBest = rep(NA,tsize)

   PathAB = rep(NA,tsize)
   PathBA = rep(NA,tsize)
   SEPathAB = rep(NA,tsize)
   SEPathBA = rep(NA,tsize)
   ZPathAB = rep(NA,tsize)
   ZPathBA = rep(NA,tsize)
   PPathAB = rep(NA,tsize)
   PPathBA = rep(NA,tsize)
   BLV.AtoB = rep(NA,tsize)
   BLV.BtoA = rep(NA,tsize)
   leo.nb.AtoB = rep(NA,tsize)
   leo.nb.BtoA = rep(NA,tsize)
   AcollideB.vs.NextBest = rep(NA,tsize)


 rn= rep(NA,tsize)
 rn1= rep(NA,tsize)
 rn2= rep(NA,tsize)
 rn3= rep(NA,tsize)
 rn4= rep(NA,tsize)
 rn5= rep(NA,tsize)

 ii=0; # store at this row
 for (i in 1:length(traitcols) ) {
  t=traitcols[i];
  for (j in 1:length(snpcols)) {
    s = snpcols[j];
     for (k in 1:length(genecols)) {
       g = genecols[k];

  ii=ii+1
  xuse=datCombined[,c(s,g,t)]
  xuse.imp = xuse
  if (impute.na) { xuse.imp = impute(xuse,what="median") }
  z=compare.local.sems(pm=pm,M.col,A.col,B.col,xuse.imp)

  mlogp.M.AtoB[ii] = signif( z$mlogp.M.AtoB,3)
  mlogp.M.BtoA[ii] = signif( z$mlogp.M.BtoA,3)
  mlogp.M.conf[ii] = signif( z$mlogp.M.conf,3)
  mlogp.M.AcollideB[ii] = signif( z$mlogp.M.AcollideB,3)
  mlogp.M.BcollideA[ii] = signif( z$mlogp.M.BcollideA,3)
#  mlogp.M.AhiddenB[ii] = signif( z$mlogp.M.AhiddenB,3)
  rn[ii]=paste(sep="","M:",cn[s]," to A:",cn[g]," to B:",cn[t])

  rn1[ii]=paste(sep="","M:",cn[s])
  rn2[ii]="to"
  rn3[ii]=paste(sep="","A:",cn[g])
  rn4[ii]="to"
  rn5[ii]=paste(sep="","B:",cn[t])
  
  s1.i=z$mlogp.M.conf - z$mlogp.M.AtoB;
  s2.o=z$mlogp.M.AcollideB - z$mlogp.M.AtoB;
  leo.i[ii] = signif((s1.i),3)
  leo.o[ii] = signif((s2.o),3)
  leo.ab.over.ba[ii] = signif( z$mlogp.M.BtoA - z$mlogp.M.AtoB ,3)


   ##
   ## the new statistics ZPathAB, ZPathBA, z$LEO.NB.BtoA,
   ## BLV, PearsonCor, PearsonCorP
   ##

   PathAB[ii] = signif(z$PathAB,3)
   PathBA[ii] = signif(z$PathBA,3)

   SEPathAB[ii] = signif(z$SEPathAB,3)
   SEPathBA[ii] = signif(z$SEPathBA,3)

   ZPathAB[ii] = signif(z$ZPathAB,3)
   ZPathBA[ii] = signif(z$ZPathBA,3)

   PPathAB[ii] = signif(z$PPathAB,3)
   PPathBA[ii] = signif(z$PPathBA,3)

   BLV.AtoB[ii] = signif(z$BLV.AtoB,3)
   BLV.BtoA[ii] = signif(z$BLV.BtoA,3)

   RMSEA.AtoB[ii]=signif(z$M.AtoB$RMSEA[1],3)
   RMSEA.BtoA[ii]=signif(z$M.BtoA$RMSEA[1],3)
   RMSEA.conf[ii]=signif(z$M.conf$RMSEA[1],3)
   RMSEA.AcollideB[ii]=signif(z$M.AcollideB$RMSEA[1],3)
   RMSEA.BcollideA[ii]=signif(z$M.BcollideA$RMSEA[1],3)

    s3.nb = z$LEO.NB.AtoB
    leo.nb.AtoB[ii] = signif(s3.nb,3)
    s3.nb.BtoA = z$LEO.NB.BtoA
    leo.nb.BtoA[ii] = signif(s3.nb.BtoA,3)
    AcollideB.vs.NextBest[ii] = signif(z$LEO.NB.AcollideB,3)

    good.o[ii] = good.i[ii] = good.nb[ii] = " "
    if (s1.i >= (leo.i.th)) { good.i[ii]="*" }
    if (s2.o >= (leo.o.th)) { good.o[ii]="*" }
    if (s3.nb >= (leo.nb.th)) { good.nb[ii]="*" }

  }}}

  leo.nb.AtoB = signif(leo.nb.AtoB,3)
  leo.nb.BtoA = signif(leo.nb.BtoA,3)

  df=data.frame(rn1,rn2,rn3,rn4,rn5, leo.nb.AtoB, leo.nb.BtoA, leo.i, leo.o, good.o, leo.ab.over.ba, mlogp.M.AtoB,mlogp.M.BtoA,mlogp.M.conf,mlogp.M.AcollideB,mlogp.M.BcollideA, PathAB, SEPathAB, ZPathAB, PPathAB, PathBA, SEPathBA, ZPathBA, PPathBA, BLV.AtoB, BLV.BtoA, AcollideB.vs.NextBest, RMSEA.AtoB, RMSEA.BtoA, RMSEA.conf, RMSEA.AcollideB, RMSEA.BcollideA) 

  colnames(df)[1]="model"
  colnames(df)[2:5]=""
##  rownames(df) = rn
  df
}




#
# Generate a confusion matrix in place of type I error analysis:
#  superior in most senses...do models 1,2,3,4,5
#
# p-value and specificity analysis
# (type 1 error analysis): simulate under null hyp of model 3 (confounding)
# do the ZEO stuff to compare
if(exists("typeIconfounding.study") ) rm(typeIconfounding.study);
typeIconfounding.study=function() {

   s2=sqrt(2)

   # simulate
   no.samples=400
   set.seed(1)
   gen.var=1;
   env.var=1;

   d1=date()

   g.to.e = c(.1,.2,.5,1,2);
   effect.coef=c(2,1,.5,.25,.1);
   n.choice=c(50,100,200,500,1000,2000,5000);

   ma.choice=c(.1,.25,.5,.75,1)
   mb.choice=c(.1,.25,.5,.75,1)

   mal=length(ma.choice);
   mbl=length(mb.choice);

   lec=length(effect.coef);
   lge=length(g.to.e);
   lnc=length(n.choice);

   res=c();

   eo.mean=array(data = 0, dim = c(mal,mbl,lnc), dimnames = list(as.character(ma.choice),as.character(mb.choice),as.character(n.choice)));
   eo.sd =eo.mean

   ABeff=1;

   # different true models simulated...
   models.to.do=c(3)
   for(i in models.to.do) {

   for (n in 1:lnc) {
     no.samples=n.choice[n]

   for (a in 1:mal) {
      MAeff=ma.choice[a];

   for (b in 1:mbl) {
      MBeff=mb.choice[b];

   num.sim=10
   model.num=1
#   disc.matrix=matrix(nrow=6,ncol=6,data=0);
   eo.m = vector(length=num.sim);
   eo.m[]=0;

      no.obs.Z=no.samples
      for (k in 1:num.sim) {

         mod=sim.model.num(model.num=i,env.var,no.samples,MAeff,ABeff,MBeff);

         covx=cov(mod)

            # twofer.eo and n.twofer.eo hold our results...
            z=zeo2(1,3,2,covx,pm); # z$eo

            #cor1=pcor(c(1,3,c()),covx)
            #Zm1B= sqrt(no.obs.Z-pm$fisher.dof.cor)*   fisher1(abs(cor1))

            #cor2=pcor(c(1,3,2),covx)
            #Zm1BgivenA= sqrt(no.obs.Z-pm$fisher.dof.pcor1)*   fisher1(abs(cor2))

            #s=1;
            #eo2 = Zm1B - s*Zm1BgivenA;
            eo2 = z$eo;
         eo.m[k]=eo2;
     } # k

        eo.mean[a,b,n]=mean(eo.m)
        eo.sd[a,b,n]=sd(eo.m)

        cat(paste(a,b,n,date(),"\n"))
     } # e
     } # g
     } # n

  d2=date()
  z=list()

  z$eo.mean=eo.mean
  z$eo.sd=eo.sd

  z$d1=d1
  z$d2=d2

  z$eo.m=eo.m
  save(z,file=paste(sep="","typeIconfounding.study.zeo.i.is.",i,".reps.",num.sim,".rdat"))    
  res=c(z,res);

  } # end i over 6 models

  res
}

#
# implement option to add covariances within markers on the m1m2 models
# (see   pm$add.MA.MB.covar.in.local.sem.four.var.m1m2)
#
if(exists("generate.intra.ma.pairlist") ) rm(generate.intra.ma.pairlist);
generate.intra.ma.pairlist=function(MA.col,MB.col,pm) {
  list.of.hidden.confounded.pair.indices = NULL
  
  if (pm$add.MA.MB.covar.in.local.sem.four.var.m1m2 & (length(MA.col) > 1 | length(MB.col) > 1)) {
    list.of.hidden.confounded.pair.indices = list()
    list.pos = 1
    if (length(MA.col) > 1) {
      for (i in 1:(length(MA.col)-1)) {
        for (j in (i+1):(length(MA.col))) {
          list.of.hidden.confounded.pair.indices[[ list.pos ]] = c(MA.col[i],MA.col[j])
          list.pos = list.pos +1
        }
      }
    }
    
    if (length(MB.col) > 1) {
      for (i in 1:(length(MB.col)-1)) {
        for (j in (i+1):(length(MB.col))) {
          list.of.hidden.confounded.pair.indices[[ list.pos ]] = c(MB.col[i],MB.col[j])
          list.pos = list.pos +1
        }
      }
    }
    
  } # end if (pm$add.MA.MB...
  
  list.of.hidden.confounded.pair.indices
} # end generate.intra.ma.pairlist



# also safe if MA.col or MB.col is a vector
#
# fit.models can be reduced to save time if you only want to fit fewer models
# (confounded model takes a handful of iterations, generally)
#
if(exists("local.sem.four.var.m1m2") ) rm(local.sem.four.var.m1m2);
local.sem.four.var.m1m2=function(MA.col,A.col,B.col,x,MB.col,no.obs=NULL,fit.models=c(1,2,3,4),ten.correct.chisq=FALSE,pm) {
   pm.cat(pm,"\n");
   pm.print(pm,paste("Starting local.sem.four.var.m1m2()...",date()));

   # necessary for xr1 code below to always work...
   if(length(A.col) !=1 || length(B.col) != 1) {
       stop("local.sem.four.var.m1m2 called with illegal input: A or B had more or less than 1 variable specified.")
   }

   if (is.null(no.obs)) { no.obs=nrow(x); }

   if (is.null(pm$no.obs.Z)) {pm$no.obs.Z = no.obs} # so we get correct zeo2() scores

   # correct the chisquare values that get off as N varies from 100
   if (ten.correct.chisq && no.obs > 100) {no.obs = 10 * sqrt(no.obs)  }

   cn=colnames(x);

   if (is.character(MA.col)) { MA.col = which(cn==MA.col);  }
   if (is.character(MB.col)) { MB.col = which(cn==MB.col);  }

   dn=cn[c(MA.col,MB.col,A.col,B.col)];

   subx=data.frame(x[,c(MA.col,MB.col,A.col,B.col)])
   covx=cov(subx);

   # adjust indices to point into new matrix
   MA.col =match(cn[MA.col],dn)
   MB.col =match(cn[MB.col],dn)
   A.col  =match(cn[A.col],dn)
   B.col  =match(cn[B.col],dn)

     z=list()

     z$title=paste(sep="","MA(",paste(dn[MA.col],collapse=","),"),MB(",paste(dn[MB.col],collapse=","),"),A(",paste(dn[A.col],collapse=","),"),B(",paste(dn[B.col],collapse=","),") local SEM model fit of MA->A->B<-MB.")

     # we want to know how much of A the markers explain, and how much of B
     form1=paste(dn[A.col],"~",paste(collapse=" + ",dn[MA.col]))
     summy1=summary(lm(as.formula(form1),data=subx))
     z$A.predicted.by.MA.markers.R.squared = summy1$r.squared
     z$A.predicted.by.MA.markers.Adjusted.R.squared = summy1$adj.r.squared
     z$A.predicted.by.MA.markers.formula = form1

     form2=paste(dn[B.col],"~",paste(collapse=" + ",dn[MB.col]))
     summy2=summary(lm(as.formula(form2),data=subx))
     z$B.predicted.by.MB.markers.R.squared = summy2$r.squared
     z$B.predicted.by.MB.markers.Adjusted.R.squared = summy2$adj.r.squared
     z$B.predicted.by.MB.markers.formula = form2

     form3=paste(dn[B.col],"~",dn[A.col])
     summy3=summary(lm(as.formula(form3),data=subx))
     z$B.predicted.by.A.R.squared = summy3$r.squared
     z$B.predicted.by.A.Adjusted.R.squared = summy3$adj.r.squared
     z$B.predicted.by.A.formula = form3


     # 4 variable fit
     nr=length(dn);
     M.empty=matrix(rep(0,nr*nr),nrow=nr,dimnames=list(dn,dn));

     intra.M.list = generate.intra.ma.pairlist(MA.col,MB.col,pm)


     ana.grad=T; # otherwise we sometimes get: Error in nlm(if (analytic.gradient) objective.2 else objective.1, start,  : probable coding error in analytic gradient

     if (1 %in% fit.models) {

     M.M1M2AtoB=M.empty;
     M.M1M2AtoB[MA.col,A.col]=1;
     M.M1M2AtoB[MB.col,B.col]=1;
     M.M1M2AtoB[A.col,B.col]=1;
     
       sem.M.M1M2.AtoB = make.ram(M.M1M2AtoB, list.of.hidden.confounded.pair.indices=intra.M.list)
       full.sem.output.M1M2.AtoB = try.sem(pm,sem.M.M1M2.AtoB$the.ram,covx,N=no.obs,ana.grad)
       z$M.M1M2.AtoB = summary(full.sem.output.M1M2.AtoB)
       z$mlogp.M.M1M2.AtoB=-pchisq(z$M.M1M2.AtoB$chisq,df=z$M.M1M2.AtoB$df,lower.tail=FALSE,log.p=TRUE)*log.to.log10
     }
   
     if (2 %in% fit.models) {
        M.M1M2BtoA=M.empty;
        M.M1M2BtoA[MA.col,A.col]=1;
        M.M1M2BtoA[MB.col,B.col]=1;
        M.M1M2BtoA[B.col,A.col]=1;

       sem.M.M1M2.BtoA = make.ram(M.M1M2BtoA, list.of.hidden.confounded.pair.indices=intra.M.list)
       full.sem.output.M1M2.BtoA = try.sem(pm,sem.M.M1M2.BtoA$the.ram,covx,N=no.obs,ana.grad)
       z$M.M1M2.BtoA = summary(full.sem.output.M1M2.BtoA)
       z$mlogp.M.M1M2.BtoA=-pchisq(z$M.M1M2.BtoA$chisq,df=z$M.M1M2.BtoA$df,lower.tail=FALSE,log.p=TRUE)*log.to.log10
     }

     if (3 %in% fit.models) {
        M.M1M2unresolv=M.empty;
        M.M1M2unresolv[MA.col,A.col]=1;
        M.M1M2unresolv[MB.col,B.col]=1;
        M.M1M2unresolv[MA.col,B.col]=1;
        M.M1M2unresolv[MB.col,A.col]=1;
       
        sem.M.M1M2unresolv = make.ram(M.M1M2unresolv, list.of.hidden.confounded.pair.indices=intra.M.list)
        full.sem.output.M1M2.unresolv = try.sem(pm,sem.M.M1M2unresolv$the.ram,covx,N=no.obs,ana.grad)
        z$M.M1M2unresolv = summary(full.sem.output.M1M2.unresolv)
        z$mlogp.M.M1M2unresolv=-pchisq(z$M.M1M2unresolv$chisq,df=z$M.M1M2unresolv$df,lower.tail=FALSE,log.p=TRUE)*log.to.log10
      }

     if (4 %in% fit.models) {
        M.M1M2hidden.con=M.empty;
        M.M1M2hidden.con[MA.col,A.col]=1;
        M.M1M2hidden.con[MB.col,B.col]=1;

        sem.M.M1M2hidden.con = make.ram(M.M1M2hidden.con, list.of.hidden.confounded.pair.indices=intra.M.list)

        # manually add in the covariance between A.col and B.col that is due to hidden var
        xr1a=paste(dn[A.col],"<->",dn[B.col]); # the extra row in the making
        xr1b=paste(sep="","hidden.confound.of.",dn[A.col],".",dn[B.col])
        xr1c="NA";
        xr1=cbind(xr1a,xr1b,xr1c);
        dimnames(xr1) = NULL

        df = sem.M.M1M2hidden.con$the.ram;
        df2=cbind(df[,1],df[,2],df[,3])
     
        sem.M.M1M2hidden.con$the.ram = rbind(df2,xr1);
        class(sem.M.M1M2hidden.con$the.ram) <- "mod"

        # end manual addition of hidden covar between A and B

        full.sem.output.M1M2.hidden.con = try.sem(pm,sem.M.M1M2hidden.con$the.ram,covx,N=no.obs,ana.grad)
        z$M.M1M2hidden.con = summary(full.sem.output.M1M2.hidden.con)
        z$mlogp.M.M1M2hidden.con=-pchisq(z$M.M1M2hidden.con$chisq,df=z$M.M1M2hidden.con$df,lower.tail=FALSE,log.p=TRUE)*log.to.log10
      }

     # additionally combine the unresolvable and the confounded: comboNull. The
     # problem with this model is that it is so saturated that it can fit everything well. And few df to test.
     if (5 %in% fit.models) {
        M.M1M2comboNull=M.empty;
        M.M1M2comboNull[MA.col,A.col]=1;
        M.M1M2comboNull[MB.col,B.col]=1;
        M.M1M2comboNull[MA.col,B.col]=1;
        M.M1M2comboNull[MB.col,A.col]=1;
       
        sem.M.M1M2comboNull = make.ram(M.M1M2comboNull, list.of.hidden.confounded.pair.indices=intra.M.list)

        # manually add in the covariance between A.col and B.col that is due to hidden var
        xr1a=paste(dn[A.col],"<->",dn[B.col]); # the extra row in the making
        xr1b=paste(sep="","hidden.confound.of.",dn[A.col],".",dn[B.col])
        xr1c="NA";
        xr1=cbind(xr1a,xr1b,xr1c);
        dimnames(xr1) = NULL

        df = sem.M.M1M2comboNull$the.ram;
        df2=cbind(df[,1],df[,2],df[,3])
     
        sem.M.M1M2comboNull$the.ram = rbind(df2,xr1);
        class(sem.M.M1M2comboNull$the.ram) <- "mod"

        # end manual addition of hidden covar between A and B

        full.sem.output.M1M2.comboNull = try.sem(pm,sem.M.M1M2comboNull$the.ram,covx,N=no.obs,ana.grad)
        z$M.M1M2comboNull = summary(full.sem.output.M1M2.comboNull)
        z$mlogp.M.M1M2comboNull=-pchisq(z$M.M1M2comboNull$chisq,df=z$M.M1M2comboNull$df,lower.tail=FALSE,log.p=TRUE)*log.to.log10
      }



     # zeo scores
     z$zeo.MA.B.given.A=multimarker.zeo2detail(MA.col,B.col,A.col,covx,pm) # $BLV is the final score
     z$zeo.MB.A.given.B=multimarker.zeo2detail(MB.col,A.col,B.col,covx,pm) # $BLV


  if(length(fit.models)>=4) {
    # report all the summary stats too...
    z$alt.model.best.min.mlogp = min(z$mlogp.M.M1M2.BtoA,z$mlogp.M.M1M2unresolv, z$mlogp.M.M1M2hidden.con);
    z$main.model.M1M2.A.to.B.mlogp=z$mlogp.M.M1M2.AtoB;

#    z$LEO.NB.AtoB = -(z$mlogp.M.M1M2.AtoB - z$alt.model.best.min.mlogp);

    z$p.conf.over.p.AtoB = pchisq(z$M.M1M2hidden.con$chisq,1,lower.tail=FALSE)/pchisq(z$M.M1M2.AtoB$chisq,1,lower.tail=FALSE);
    z$p.conf.over.p.BtoA = pchisq(z$M.M1M2hidden.con$chisq,1,lower.tail=FALSE)/pchisq(z$M.M1M2.BtoA$chisq,1,lower.tail=FALSE);

    # the -log10 space version of p.conf.over.p.AtoB becomes LEO.I.AtoB (I for independent (confounded) model comparison)
    z$LEO.I.AtoB = z$mlogp.M.M1M2hidden.con - z$mlogp.M.M1M2.AtoB
    z$LEO.I.BtoA = z$mlogp.M.M1M2hidden.con - z$mlogp.M.M1M2.BtoA

    # new stats

    # locate the A->B path coefficient row in the output...

    ab.row = match(paste(sep="",dn[A.col]," -> ",dn[B.col]),sem.M.M1M2.AtoB$the.ram[,1])
    ba.row = match(paste(sep="",dn[B.col]," -> ",dn[A.col]),sem.M.M1M2.BtoA$the.ram[,1])

   z$PathAB = z$M.M1M2.AtoB$coeff[ab.row,1]
   z$SEPathAB = z$M.M1M2.AtoB$coeff[ab.row,2]
   z$ZPathAB = z$M.M1M2.AtoB$coeff[ab.row,3]
   z$PPathAB = z$M.M1M2.AtoB$coeff[ab.row,4]
   z$BLV.AtoB = z$zeo.M.A.given.B$BLV

   z$PathBA = z$M.M1M2.BtoA$coeff[ba.row,1]
   z$SEPathBA = z$M.M1M2.BtoA$coeff[ba.row,2]
   z$ZPathBA = z$M.M1M2.BtoA$coeff[ba.row,3]
   z$PPathBA = z$M.M1M2.BtoA$coeff[ba.row,4]
   z$BLV.BtoA = z$zeo.M.B.given.A$BLV

   z$LEO.NB.AtoB = z$AB.vs.NextBest = -(z$mlogp.M.M1M2.AtoB - min(z$mlogp.M.M1M2.BtoA,z$mlogp.M.M1M2unresolv, z$mlogp.M.M1M2hidden.con)) # aka LEO.NB.AtoB
   z$LEO.NB.BtoA = z$BA.vs.NextBest = -(z$mlogp.M.M1M2.BtoA - min(z$mlogp.M.M1M2.AtoB, z$mlogp.M.M1M2unresolv, z$mlogp.M.M1M2hidden.con)) # aka LEO.NB.BtoA
   z$BLV.AtoB = z$zeo.MA.B.given.A$mean.BLV
   z$BLV.BtoA = z$zeo.MB.A.given.B$mean.BLV
    
    # instead of next best, to get a symmetric score...we go AtoB against BtoA, in the M1M2 two marker models.
    z$eo.losem.lod = -(z$mlogp.M.M1M2.AtoB - z$mlogp.M.M1M2.BtoA);

    # give the modern names...LEO.O
    z$LEO.O.AtoB = -(z$mlogp.M.M1M2.AtoB - z$mlogp.M.M1M2.BtoA);
    z$LEO.O.BtoA = -(z$mlogp.M.M1M2.BtoA - z$mlogp.M.M1M2.AtoB);

    z$RMSEA.AtoB = z$M.M1M2.AtoB$RMSEA
    z$RMSEA.BtoA = z$M.M1M2.BtoA$RMSEA

  }
  z
}

# stringify()
#
# Print an object into a string, and return the string.
# Useful for capturing the results of nice formating that
# the print and summary functions do for the screen.
#
# Used alot for logging in our Local SEM model fitting functions.
#
if(exists("stringify") ) rm(stringify);
stringify=function(x,...) {
   paste(capture.output(print(x)),collapse="\n")
}

# version that cat x instead
if(exists("stringify.cat") ) rm(stringify.cat);
stringify.cat=function(x,...) {
   paste(capture.output(cat(x)),collapse="\n")
}

# check if we should ignore this edge
if(exists("edge.is.on.ignore.list")) rm(edge.is.on.ignore.list);
edge.is.on.ignore.list=function(i,j,pm) {
   if (is.null(pm$ignorable.edge.list) || length(pm$ignorable.edge.list) == 0) { return(FALSE) }
   for (k in 1:length(pm$ignorable.edge.list)) {
       a=pm$ignorable.edge.list[[k]]
       if ((a[1]==i && a[2]==j) || a[1]==j && a[2]==i) return(TRUE);
   }
   FALSE
}

# check if we should use this edge
if(exists("is.A.B.edge")) rm(is.A.B.edge);
is.A.B.edge=function(i,j,pm) {

   if (is.null(pm$A) | is.null(pm$B)) return(TRUE); # filter not active

   if (i %in% pm$A.cn & j %in% pm$B.cn) return(TRUE)
   if (j %in% pm$A.cn & i %in% pm$B.cn) return(TRUE)

   FALSE
}



if(exists("marker.to.string.of.marker.and.prob") ) rm(marker.to.string.of.marker.and.prob);
marker.to.string.of.marker.and.prob=function(Mset, cn, m1ab.log.list,wh) {
   m.cn = cn[Mset[wh]]
   model.probs = rep(NA,length(m.cn))
   model.probs.sortable = model.probs
   k=1
   for (w in wh) {
     model.probs.sortable[k]=m1ab.log.list[[w]]$mlogp.M.AtoB
     model.probs[k]=paste("(ModelProb:",signif(10^(-(m1ab.log.list[[w]]$mlogp.M.AtoB)),3),")",sep="")
     k=k+1
   }
   ix=sort(model.probs.sortable,index.return=TRUE)$ix
   paste(sep="",m.cn[ix],model.probs[ix],collapse=",")
}



########################################
# LOCAL SEM VERSION
# Max-max two-step strategy is different:
# we compute only one total ZEO. From two
# halfs, each of which is maximized independently.
########################################

if(exists("losem.walk.two.steps.max.max") ) rm(losem.walk.two.steps.max.max);
losem.walk.two.steps.max.max=function(x,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log,zeo.proxy=NA) {

   if (pm$do.leo.max) { pm.print(pm,paste("Starting losem.walk.two.steps.max.max()...",date())); }

   # initializer for LEO.I, BLV, ZPathAB, AB.vs.NextBest
   na.init=twofer.eo; na.init[]=NA;

   # track the relative probability of the confounded model
   LEO.I = na.init

   # track the newly added statistics, initialized to NA
   BLV = na.init
   ZPathAB = na.init
   AB.vs.NextBest = na.init
   
   # best single marker a->b RMSEA
   # best single marker b->a RMSEA
   BestAB.single.marker.RMSEA = na.init
   BestAB.single.marker.RMSEA.the.marker = na.init

   Simple.Max.vs.Max = na.init

  # setup and return a list with NULL entries if we aren't doing this function
   if (!pm$do.leo.max) {
      ret=list()
      ret$twofer.eo = twofer.eo # aka LEO.O
      ret$twofer.log = twofer.log
      ret$LEO.I = LEO.I
      ret$BLV = BLV
      ret$ZPathAB = ZPathAB
      ret$AB.vs.NextBest = AB.vs.NextBest #aka LEO.NB
      ret$BestAB.single.marker.RMSEA = BestAB.single.marker.RMSEA
      ret$BestAB.single.marker.RMSEA.the.marker = BestAB.single.marker.RMSEA.the.marker
      ret$Simple.Max.vs.Max = Simple.Max.vs.Max
      return(ret)
   }

   
  for (i in 1:(length(non.snp.cols)-1)) {
   for (j in (i+1):length(non.snp.cols)) {
      A= non.snp.cols[i];
      B= non.snp.cols[j];

      # use the A,B edge filters if present
      if (!is.null(pm$A) & !is.A.B.edge(cn[A],cn[B],pm)) next;

      if(edge.is.on.ignore.list(A,B,pm)) { twofer.log[A,B]="Ignored edge b/c on ignore list";  next; }
      if(!is.na(pm$zeo.proxy.for.leo.th) && !is.na(zeo.proxy[A,B])) {
         if (abs(zeo.proxy[A,B]) < pm$zeo.proxy.for.leo.th) { twofer.log[A,B]=paste("Edge eliminated by zeo.proxy score of",signif(zeo.proxy[A,B],3),"with elimination th:",pm$zeo.proxy.for.leo.th);  next; }
      }

      pm.cat(pm,paste(sep="",cn[A],"->",cn[B],"...")); # give a sense of progress
      
      if (af[A,B]==0 && af[B,A]==0) next; 
      # INVAR: we have an edge from A to B possible.
      M.A = setdiff(which(af[,A]==1),non.snp.cols); # the Markers into A
      M.B = setdiff(which(af[,B]==1),non.snp.cols); # the Markers into B


      # make sure any overlapping belong only to the set with the strongest correlation...
      MAB=union(M.A,M.B);
      M.A=M.B=c();
      for (m in MAB) {
         if (cx[m,A] >= cx[m,B]) { M.A=c(M.A,m); } else { M.B = c(M.B,m); }
      }

      Mset=c(); # grandparent candidates
      parent=NA;
      child=NA;
      singleAB=1; # 1 for M->A->B; -1 for M->B->A
      if (length(M.A) == 0 || length(M.B) == 0)
      {
#         if (!pm$work.hard.at.EO.contingencies.if.not.both.M.B.and.M.A.markers) {
#              twofer.log[A,B]=paste("Max.max contingency: markers only into one of A and B and not stronger into one. Since pm$work.hard.at.EO.contingencies.if.not.both.M.B.and.M.A.markers is FALSE, we don't report a single marker Max score that would be on a different scale from the two-marker scores.")
#
#              next;
#         }

          if(length(M.A) > 0) {  
              # do a single marker score, one two step, or half a ZEO: M.A -> A -> B
              Mset=M.A;
              parent=A;
              child=B;
              singleAB=1;

              # now compute A-B edge orientation using just the one marker

              # space for results:
              m1ab = cx[Mset,parent]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[Mset];}
              m1ab.log.list=vector("list", length(m1ab)) # create empty list of given length
              m1ab.rmsea = m1ab
              max.k = NA

              for (k in 1:length(Mset)) {
                  m1 = Mset[k];
                  z=compare.local.sems(pm=pm,m1,parent,child,x);
                  m1ab[k] = z$eo.losem.lod;
                  m1ab.rmsea[k] = z$M.AtoB$RMSEA[1]
                  m1ab.log.list[[k]]=z;
              }
              m1ab.names.presort=names(m1ab);
              m1ab=sort(m1ab,decreasing=T);

              top.m1 = which(m1ab.names.presort==names(m1ab)[1]); # figure out which log entry to keep.

              twofer.eo[A,B]=(m1ab[1]);
              twofer.eo[B,A]=-twofer.eo[A,B];

        BestAB.single.marker.RMSEA[A,B] = min(m1ab.rmsea)
        BestAB.single.marker.RMSEA.the.marker[A,B] = marker.to.string.of.marker.and.prob(Mset,cn,m1ab.log.list,wh=which(m1ab.rmsea == min(m1ab.rmsea)))  # get all the markers with minimal RMSEA, if multiple.

              BestAB.single.marker.RMSEA[B,A] = NA # leave blank b/c no markers stronger into B than A

              twofer.log[A,B]=paste("Max.max contingency: Only M.A markers, no M.B markers (all markers stronger into A than B)\n",stringify(m1ab.log.list[[top.m1]]))

              # recompute z for the top pick, to get all the right summary statistics
              z=compare.local.sems(pm=pm,Mset[top.m1],parent,child,x);

              # new stats 
              LEO.I[A,B]=z$LEO.I.AtoB;
              LEO.I[B,A]=z$LEO.I.BtoA;
              
              AB.vs.NextBest[A,B] = z$LEO.NB.AtoB
              AB.vs.NextBest[B,A] = z$LEO.NB.AcollideB
              
              BLV[A,B] = z$BLV.AtoB
              BLV[B,A] = z$BLV.BtoA
              
              ZPathAB[A,B] = z$ZPathAB
              ZPathAB[B,A] = z$ZPathBA
              
          } else if (length(M.B) > 0) {
              # interesting situation here. We have no markers into A.
              # We have marker(s) into B. We wish to return a score
              # for the A->B edge into twofer.eo[A,B] that reflects
              # M->B<-A situation. Hence we added a score to compare.local.sems()
              # to reflect this situation; we want to use here 
              #
              # compare.local.sems()$LEO.NB.AcollideB
              #
              #     rather than
              #
              # compare.local.sems()$eo.losem.lod
              #

              # do a single marker score
              Mset=M.B
              parent=B;
              child=A;
              singleAB=-1;

              # now compute A-B edge orientation using just the one marker

              # space for results:
              m1ab = cx[Mset,parent]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[Mset];}
              m1ab.log.list=vector("list", length(m1ab)) # create empty list of given length
              m1ab.rmsea = m1ab

              for (k in 1:length(Mset)) {
                  m1 = Mset[k];
                  z=compare.local.sems(pm=pm,m1,parent,child,x);
                  m1ab[k] = z$LEO.NB.AcollideB;
                  m1ab.rmsea[k] = z$M.AtoB$RMSEA[1]
                  m1ab.log.list[[k]]=z;
              }
              m1ab.names.presort=names(m1ab);
              m1ab = sort(m1ab,decreasing=T) # index.return=T will loose the labels we need
              top.m1 = which(m1ab.names.presort==names(m1ab)[1]); #figure out which m1ab.log.list entry to keep

              #top.ba = which.max(m1ba)

              twofer.eo[A,B]=m1ab[1];
              twofer.eo[B,A]=-twofer.eo[A,B];

              twofer.log[A,B]=paste("Max-max contingency: Only M.B markers, no M.A markers (all markers stronger into B than A).\n",stringify(m1ab.log.list[[top.m1]]))

              BestAB.single.marker.RMSEA[A,B] = NA # no M.A markers, so we leave blank

        BestAB.single.marker.RMSEA[B,A] = min(m1ab.rmsea)

        BestAB.single.marker.RMSEA.the.marker[B,A] = marker.to.string.of.marker.and.prob(Mset,cn,m1ab.log.list,wh=which(m1ab.rmsea == min(m1ab.rmsea)))  # get all the markers with minimal RMSEA, if multiple.

              # recompute z for the top pick, to get all the right summary statistics
              z=compare.local.sems(pm=pm,Mset[top.m1],parent,child,x);

              # new stats 
              LEO.I[A,B]=z$LEO.I.BtoA;
              LEO.I[B,A]=z$LEO.I.AtoB;
              
              AB.vs.NextBest[A,B] = z$LEO.NB.BtoA
              AB.vs.NextBest[B,A] = z$LEO.NB.AtoB
              
              BLV[A,B] = z$BLV.BtoA
              BLV[B,A] = z$BLV.AtoB
              
              ZPathAB[A,B] = z$ZPathBA
              ZPathAB[B,A] = z$ZPathAB
              
          } else next;

          next; # b/c rest of code is for markers into both A and B
      }

      if (length(M.A) == 0 || length(M.B) == 0) next; # require markers into both A and B.

      # INVAR: M.A has markers with stronger (>=) connection to A than B
      #    and M.B has markers with stronger connection to B than A

      # compute the Z.2steps
      # reserve space, and get the names() right
      m1ab = cx[M.A,A]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[M.A];}
      m2ba = cx[M.B,B]; m2ba[]=NA; if (length(m2ba)==1) { names(m2ba)=rownames(cx)[M.B];}

      m1ab.rmsea = m1ab
      m2ba.rmsea = m2ba

      m1ab.AB.prob = m1ab
      m2ba.AB.prob = m2ba

      m1ab.log.list=vector("list", length(m1ab)) # create empty list of given length
      m2ba.log.list=vector("list", length(m2ba)) # create empty list of given length      

      for (k in 1:length(M.A)) {
          m1 = M.A[k];

          z=compare.local.sems(pm=pm,m1,A,B,x);
          m1ab[k] = z$eo.losem.lod;
          m1ab.rmsea[k] = z$M.AtoB$RMSEA[1]
          m1ab.AB.prob[k] = z$mlogp.M.AtoB
          m1ab.log.list[[k]]=z;
      }

      for (k in 1:length(M.B)) {
          m2 = M.B[k];

          z=compare.local.sems(pm=pm,m2,B,A,x);
          m2ba[k] = z$eo.losem.lod;
          m2ba.rmsea[k] = z$M.AtoB$RMSEA[1]
          m2ba.AB.prob[k] = z$mlogp.M.AtoB
          m2ba.log.list[[k]]=z;
      }

      m1ab.unsorted = m1ab
      m2ba.unsorted = m2ba

      m1ab=sort(m1ab,decreasing=T)
      m2ba=sort(m2ba,decreasing=T)

      BestAB.single.marker.RMSEA[A,B] = min(m1ab.rmsea)

      BestAB.single.marker.RMSEA.the.marker[A,B] = marker.to.string.of.marker.and.prob(M.A,cn,m1ab.log.list,wh=which(m1ab.rmsea == min(m1ab.rmsea)))  # get all the markers with minimal RMSEA, if multiple.

      BestAB.single.marker.RMSEA[B,A] = min(m2ba.rmsea)
      BestAB.single.marker.RMSEA.the.marker[B,A] = marker.to.string.of.marker.and.prob(M.B,cn,m2ba.log.list,wh=which(m2ba.rmsea == min(m2ba.rmsea)))  # get all the markers with minimal RMSEA, if multiple.

      m1ab.names=names(m1ab);
      m2ba.names=names(m2ba);

      found.indep.markers=F;
      top.m1 = which(cn==names(m1ab)[1]); 
      top.m2 = which(cn==names(m2ba)[1]);

      top.m1.log.index = match(cn[top.m1], names(m1ab.unsorted))
      top.m2.log.index = match(cn[top.m2], names(m2ba.unsorted))

      # also return the Naive non-independent log10(Pr(M.A -> A -> B) / Pr(M.B -> B -> A)) for best M.A and best M.B
      min.prob.m1ab = which.min(m1ab.AB.prob)
      min.prob.m2ba = which.min(m2ba.AB.prob)
      Simple.Max.vs.Max[A,B] = - m1ab.log.list[[min.prob.m1ab]]$mlogp.M.AtoB + m2ba.log.list[[min.prob.m2ba]]$mlogp.M.AtoB
      Simple.Max.vs.Max[B,A] = - Simple.Max.vs.Max[A,B]
      
      if (cx[top.m1,top.m2] < pm$cor.ind.th) {
         #yee-haw, our top ranked markers are indeed independent.

         z=local.sem.four.var.m1m2(pm=pm,MA.col=top.m1,A.col=A,B.col=B,x,MB.col=top.m2);
         twofer.eo[A,B]=z$eo.losem.lod;
         twofer.eo[B,A]=-twofer.eo[A,B];
         twofer.log[A,B]=paste("Max-max had both M.A and M.B markers. Three models follow.\n\nTwo marker model:\n",stringify(z),"\n\nMAX One marker model forward direction:\n",stringify(m1ab.log.list[[top.m1.log.index]]),"\n\nMAX One marker model reverse direction:\n",stringify(m2ba.log.list[[top.m2.log.index]]));

         # new stats 
         LEO.I[A,B]=z$LEO.I.AtoB;
         LEO.I[B,A]=z$LEO.I.BtoA;

         AB.vs.NextBest[A,B] = z$LEO.NB.AtoB
         AB.vs.NextBest[B,A] = z$LEO.NB.BtoA

         BLV[A,B] = z$BLV.AtoB
         BLV[B,A] = z$BLV.BtoA
         
         ZPathAB[A,B] = z$ZPathAB
         ZPathAB[B,A] = z$ZPathBA

      } else {
         # top ranking markers are not independent...so scan down both lists...

         # =================================
         # pass A: M.A has priority, scan down M.B
         # =================================
         Aprior=NA
         Aprior.m1=NA;
         Aprior.m2=NA;
         Aprior.m1ab = Aprior.m2ba = NA;
         done=F;
         for (k1 in 1:length(m1ab)) {
             if (done) break;
             m1=m1ab.names[k1];
             for (k2 in 1:length(m2ba)) {
                 m2=m2ba.names[k2];
                 if (cx[m1,m2] < pm$cor.ind.th) {
                     Aprior=m1ab[k1]-m2ba[k2];
                     Aprior.m1=m1; Aprior.m2=m2;
                     Aprior.k1=k1; Aprior.k2=k2;
                     Aprior.m1ab = m1ab[k1]; Aprior.m2ba = m2ba[k2];
                     done=T;
                     break;
                 }
             }
         }
         # ===============================
         # pass B: M.B has priority, scan down M.A
         # ===============================
         Bprior=NA
         Bprior.m1=NA;
         Bprior.m2=NA;
         Bprior.m1ab = Bprior.m2ba = NA;
         done=F;
         for (k2 in 1:length(m2ba)) {
             if (done) break;
             m2=m2ba.names[k2];
             for (k1 in 1:length(m1ab)) {
                 m1=m1ab.names[k1];
                 if (cx[m1,m2] < pm$cor.ind.th) {
                     Bprior=m1ab[k1]-m2ba[k2];
                     Bprior.m1=m1; Bprior.m2=m2;
                     Bprior.k1=k1; Bprior.k2=k2;
                     Bprior.m1ab = m1ab[k1]; Bprior.m2ba = m2ba[k2];
                     done=T;
                     break;
                 }
             }
         }
         # ===============================
         # handle the 3 possible cases, with full logging so we can understand the results.
         # Either 1) Not so good: we never found independent M1 and M2, resort to single marker analysis.
         #        2) Good: found M1,M2 when M.A had priority and M.B had priority and they agree 
         #           on signs (or there was only one pair of markers). Report success at orienting the edge.
         #        3) Bad: there were two pairs of markers, differing when M.A had priority and when
         #           M.B had priority, and they differ in edge orienting signs. Report evidence is inconclusive.
         # ===============================

         if (is.na(Aprior) || is.na(Bprior)) {
             # no pair of markers in M.A and M.B was independent enough to meet our criteria
             twofer.log[A,B]=paste(sep="","no indep markers found at pm$cor.ind.th=",signif(pm$cor.ind.th,3));
         } else if (Aprior*Bprior > 0) {
             # we agreement in terms of sign-yes! very strong evidence
             # show the stronger of the two pieces of evidence: i.e. prefer the markers m1 and m2 they give strongest edge orienting support.
             if (abs(Aprior) > abs(Bprior)) { # use A

                 z=local.sem.four.var.m1m2(pm=pm,MA.col=Aprior.m1,A.col=A,B.col=B,x,MB.col=Aprior.m2);
                 twofer.eo[A,B]=z$eo.losem.lod;
                 twofer.eo[B,A]=-twofer.eo[A,B];
                 twofer.log[A,B]=paste("Max-max search for independent markers succeeded with agreement between the chosen M.A and M.B marker sign, using Aprior marker as the lead.\n\n",stringify(z))

                 # new stats 
                 LEO.I[A,B]=z$LEO.I.AtoB;
                 LEO.I[B,A]=z$LEO.I.BtoA;

                 AB.vs.NextBest[A,B] = z$LEO.NB.AtoB
                 AB.vs.NextBest[B,A] = z$LEO.NB.BtoA

                 BLV[A,B] = z$BLV.AtoB
                 BLV[B,A] = z$BLV.BtoA
         
                 ZPathAB[A,B] = z$ZPathAB
                 ZPathAB[B,A] = z$ZPathBA
                 
             } else { # use B

                 z=local.sem.four.var.m1m2(pm=pm,MA.col=Bprior.m1,A.col=A,B.col=B,x,MB.col=Bprior.m2);
                 twofer.eo[A,B]=z$eo.losem.lod;
                 twofer.eo[B,A]=-twofer.eo[A,B];
                 twofer.log[A,B]=paste("Max-max search for independent markers succeeded with agreement between the chosen M.A and M.B marker sign, using Bprior marker as the lead.\n\n",stringify(z))

                 # new stats 
                 LEO.I[A,B]=z$LEO.I.AtoB;
                 LEO.I[B,A]=z$LEO.I.BtoA;

                 AB.vs.NextBest[A,B] = z$LEO.NB.AtoB
                 AB.vs.NextBest[B,A] = z$LEO.NB.BtoA

                 BLV[A,B] = z$BLV.AtoB
                 BLV[B,A] = z$BLV.BtoA
         
                 ZPathAB[A,B] = z$ZPathAB
                 ZPathAB[B,A] = z$ZPathBA
                 
             }
         } else {
             # disagreement in terms of signs.
             twofer.eo[A,B]=twofer.eo[B,A]=NA;

       Aprior.m1.num=which(cn==Aprior.m1);
       Aprior.m2.num=which(cn==Aprior.m2);
       Bprior.m1.num=which(cn==Bprior.m1);
       Bprior.m2.num=which(cn==Bprior.m2);

             # A: display the component Z scores and correlations.
             A.sz1=paste(sep=""," [[Zm1B(",display.Zm1B(Aprior.m1.num,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(Aprior.m1.num,B,A,covx,pm),")]] \n");
             A.sz2=paste(sep=""," [[Zm2A(",display.Zm2A(Aprior.m2.num,A,covx,pm),")-Zm2A|B(",display.Zm2AgivenB(Aprior.m2.num,A,B,covx,pm),")]] \n");

             # B: display the component Z scores and correlations.
             B.sz1=paste(sep=""," [[Zm1B(",display.Zm1B(Bprior.m1.num,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(Bprior.m1.num,B,A,covx,pm),")]] \n");
             B.sz2=paste(sep=""," [[Zm2A(",display.Zm2A(Bprior.m2.num,A,covx,pm),")-Zm2A|B(",display.Zm2AgivenB(Bprior.m2.num,A,B,covx,pm),")]] \n");

             part1= paste(sep="","\nXXX Disagreement in terms of signs of M.A priority vs. M.B priority. \nXXX M1(",Aprior.m1,")->A(",cn[A],")->B(",cn[B],") had Z.2ab=",signif(Aprior.m1ab,3),A.sz1,"\nXXX M2(",Aprior.m2,")->B(",cn[B],")->A(",cn[A],") had Z.2ba=",signif(Aprior.m2ba,3),A.sz2,"\nXXX     WITH M.A ranked ",Aprior.k1," (out of ",length(m1ab),"), AND M.B ranked ",Aprior.k2," (out of ",length(m2ba),").\nXXX ZEO[A,B]=",signif(twofer.eo[A,B],3),"\nXXX Evidence from M1(",Aprior.m1,") and M2(",Aprior.m2,") with cor: ", signif(cx[Aprior.m1.num,Aprior.m2.num],3)," < pm$cor.ind.th(",signif(pm$cor.ind.th,3),")\nXXX");

             part2=paste(sep="","\nXXX M1(",Bprior.m1,")->A(",cn[A],")->B(",cn[B],") had Z.2ab=",signif(Bprior.m1ab,3),B.sz1,"\nXXX M2(",Bprior.m2,")->B(",cn[B],")->A(",cn[A],") had Z.2ba=",signif(Bprior.m2ba,3),B.sz2,"\nXXX     WITH M.A ranked ",Bprior.k1," (out of ",length(m1ab),"), AND M.B ranked ",Bprior.k2," (out of ",length(m2ba),").\nXXX ZEO[A,B]=",signif(twofer.eo[A,B],3),"\nXXX Evidence from M1(",Bprior.m1,") and M2(",Bprior.m2,") with cor: ", signif(cx[Bprior.m1.num,Bprior.m2.num],3)," < pm$cor.ind.th(",signif(pm$cor.ind.th,3),")\n");

             twofer.log[A,B]=paste(sep="","\n**[",cn[A],":",cn[B],"]**",part1,part2);
           }
        } # end top ranking markers are not independent.
   } # end i or A
} # end j or B

  pm.cat(pm,"\n");

  ret=list()
  ret$twofer.eo = twofer.eo # aka LEO.O
  ret$twofer.log = twofer.log
  ret$LEO.I = LEO.I
  ret$BLV = BLV
  ret$ZPathAB = ZPathAB
  ret$AB.vs.NextBest = AB.vs.NextBest #aka LEO.NB

  ret$BestAB.single.marker.RMSEA = BestAB.single.marker.RMSEA
  ret$BestAB.single.marker.RMSEA.the.marker = BestAB.single.marker.RMSEA.the.marker
  ret$Simple.Max.vs.Max = Simple.Max.vs.Max

  ret
}
# end losem.walk.two.steps.max.max
########################################






########################################
# LOCAL SEM VERSION
# Multiple Linear Regression two-step strategy builds a multiple regression model for A (and
#  separately B), using the markers in the model more strongly associated with A
#  than with B, by forward stepwise regression. Then we compare how well this model
#  does on B and on A|B.  Likewise for B. Combine these two R^2.
#  
########################################

if(exists("losem.walk.two.steps.weighted.mean.mlreg") ) rm(losem.walk.two.steps.weighted.mean.mlreg);
losem.walk.two.steps.weighted.mean.mlreg=function(x,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log,zeo.proxy=NA) {

   pm.print(pm,paste("Starting losem.walk.two.steps.weighted.mean.mlreg()...",date()));

   # allocate log if not already
   if (length(twofer.log)==0) { twofer.log=vector("list", length(non.snp.cols)*2) } # create empty list of given length

   # initializer for LEO.I, BLV, ZPathAB, AB.vs.NextBest
   na.init=twofer.eo; na.init[]=NA;

   # track the relative probability of the confounded model
   LEO.I = na.init

   # track the final SNPs chosen, without having to open the log file
   Final.SNPs.LEO.NB.OCA = na.init
   
   # Return model probs/p-values in spreadsheet, so we can sort without looking at the log
   Minus.log.Model.P.value.AtoB = na.init

   # track the newly added statistics, initialized to NA
   BLV = na.init
   ZPathAB = na.init
   AB.vs.NextBest = na.init
   RMSEA.2m = na.init
   
   for (i in 1:(length(non.snp.cols)-1)) {
      for (j in (i+1):length(non.snp.cols)) {
      A= non.snp.cols[i];
      B= non.snp.cols[j];

      # use the A,B edge filters if present
      if (!is.null(pm$A) & !is.A.B.edge(cn[A],cn[B],pm)) next;

      if(edge.is.on.ignore.list(A,B,pm)) { twofer.log[A,B]="Ignored edge b/c on ignore list";  next; }
      if(!is.na(pm$zeo.proxy.for.leo.th) && !is.na(zeo.proxy[A,B])) {
         if (abs(zeo.proxy[A,B]) < pm$zeo.proxy.for.leo.th) { twofer.log[A,B]=paste("Edge eliminated by zeo.proxy score of",signif(zeo.proxy[A,B],3),"with elimination th:",pm$zeo.proxy.for.leo.th);  next; }
      }

         pm.cat(pm,paste(sep="",cn[A],"->",cn[B],"...")); # give a sense of progress

      if (af[A,B]==0 && af[B,A]==0) next; 
      # INVAR: we have an edge from A to B possible.

      # Make sure any overlapping belong only to the set with the strongest correlation...M.A or M.B.
      # And make sure any SNP passes our minimum dependence threshold.

      # big Marker Assignment Consistency function moved below, so can be reused by CPA function too...
      mac=obtain.marker.assignment.consistency(pm, snpcols, cx, cn, A, B)

      M.A = mac$M.A
      M.B = mac$M.B

     if (length(M.A) == 0 && length(M.B) == 0) next;

      # handle the one side only cases...
      Mset=c(); # grandparent candidates
      parent=NA;
      child=NA;
      singleAB=1; # 1 for M->A->B; -1 for M->B->A
      if (length(M.A) == 0 || length(M.B) == 0)
      {
          if(length(M.A) > 0) {  
              # do a single marker score, one two step, or half a ZEO: M.A -> A -> B
              Mset=M.A;
              parent=A;
              child=B;
              singleAB=1;
              step.1 = forward.step.multimarker.twostep(x,cn[A],cn[M.A],cn[M.A][1]);

#             #A.coef = step.1$coefficients[-1];
              A.coef = step.1$percent.of.explained.variance.explained;
              A.who.cn  = names(step.1$coefficients)[-1]; # names
              A.who.index  = match(A.who.cn,cn); # get column index within Mset.cn

              log.summary=paste("Forward step-wise regression model built (when no M.B markers) to predict:",cn[parent]," includes these terms:\n",stringify(A.who.cn),"with coefficients",stringify(A.coef));

              z=compare.local.sems(pm=pm,A.who.index,parent,child,x);
              twofer.eo[A,B]=z$eo.losem.lod;
              twofer.eo[B,A]=-twofer.eo[A,B];
              twofer.log[A,B]=paste(stringify(z),"\n",log.summary);

              LEO.I[A,B]=z$LEO.I.AtoB;
              LEO.I[B,A]=z$LEO.I.BtoA;

              AB.vs.NextBest[A,B] = z$LEO.NB.AtoB
              AB.vs.NextBest[B,A] = z$LEO.NB.BtoA

        BLV[A,B] = z$BLV.AtoB
        BLV[B,A] = z$BLV.BtoA

              ZPathAB[A,B] = z$ZPathAB
              ZPathAB[B,A] = z$ZPathBA

              # track in spreadsheet which SNPs we ended up using
              fsnp.cpa=paste(sep="","c(\"",paste(sep="",A.who.cn,collapse="\",\""),"\")")
              fsnp.oca="c()"      
              Final.SNPs.LEO.NB.OCA[A,B]=paste(sep="","fsnp.cpa=",fsnp.cpa,";fsnp.oca=",fsnp.oca);
              Final.SNPs.LEO.NB.OCA[B,A]=paste(sep="","fsnp.cpa=",fsnp.oca,";fsnp.oca=",fsnp.cpa);

              Minus.log.Model.P.value.AtoB[A,B]=z$mlogp.M.AtoB
              Minus.log.Model.P.value.AtoB[B,A]=z$mlogp.M.BtoA
              
          } else if (length(M.B) > 0) {
              # do a single marker score, one two step, or half a ZEO: M.B -> B -> A
              Mset=M.B
              parent=B;
              child=A;
              singleAB=-1;
              step.1 = forward.step.multimarker.twostep(x,cn[B],cn[M.B],cn[M.B][1]);

#              B.coef = step.1$coefficients[-1];
              B.coef = step.1$percent.of.explained.variance.explained;
              B.who.cn  = names(step.1$coefficients)[-1]; # names
              B.who.index  = match(B.who.cn,cn); # get column index within Mset.cn

              log.summary=paste("Forward step-wise regression model built (when no M.A markers) to predict:",cn[parent]," includes these terms:\n",stringify(B.who.cn),"with coefficients",stringify(B.coef));

              z=compare.local.sems(pm=pm,B.who.index,parent,child,x);
              twofer.eo[A,B]=z$LEO.NB.AcollideB; # use reversed result here
              twofer.eo[B,A]=-twofer.eo[A,B];
              twofer.log[A,B]=paste(stringify(z),"\n",log.summary)

              LEO.I[B,A]=z$LEO.I.AtoB
              LEO.I[A,B]=z$LEO.I.BtoA

              AB.vs.NextBest[B,A] = z$LEO.NB.AtoB
              AB.vs.NextBest[A,B] = z$LEO.NB.AcollideB

        BLV[B,A] = z$BLV.AtoB
        BLV[A,B] = z$BLV.BtoA

              ZPathAB[B,A] = z$ZPathAB
              ZPathAB[A,B] = z$ZPathBA

              # track in spreadsheet which SNPs we ended up using
              fsnp.cpa=paste(sep="","c(\"",paste(sep="",B.who.cn,collapse="\",\""),"\")")
              fsnp.oca="c()"      
              Final.SNPs.LEO.NB.OCA[A,B]=paste(sep="","fsnp.cpa=",fsnp.cpa,";fsnp.oca=",fsnp.oca);
              Final.SNPs.LEO.NB.OCA[B,A]=paste(sep="","fsnp.cpa=",fsnp.oca,";fsnp.oca=",fsnp.cpa);

              Minus.log.Model.P.value.AtoB[B,A]=z$mlogp.M.AtoB
              Minus.log.Model.P.value.AtoB[A,B]=z$mlogp.M.BtoA
              
          } else next;

          next; # b/c rest of code is for markers into both A and B
      }


      if (length(M.A) == 0 || length(M.B) == 0) next; # require markers into both A and B.

      # now...both markers...
      M.A.cn = cn[M.A]
      M.B.cn = cn[M.B]

      # do the forward stepwise SNP selection for traits A and B
      step.A = forward.step.multimarker.twostep(x,cn[A],M.A.cn,M.A.cn[1],view.step=F);
      step.B = forward.step.multimarker.twostep(x,cn[B],M.B.cn,M.B.cn[1],view.step=F);

      #A.coef = step.A$coefficients[-1];
      A.coef = step.A$percent.of.explained.variance.explained;
      A.who.cn  = names(step.A$coefficients)[-1]; # names
      A.who.index  = match(A.who.cn, M.A.cn); # get column index within M.A
      A.who.index.cn = M.A[A.who.index]; # get column index within x
      A.coef.normalized = A.coef / sum(A.coef)

      #B.coef = step.B$coefficients[-1];
      B.coef = step.B$percent.of.explained.variance.explained;
      B.who.cn  = names(step.B$coefficients)[-1]; # names
      B.who.index  = match(B.who.cn, M.B.cn); # get column index within M.A
      B.who.index.cn = M.B[B.who.index]; # get column index within x
      B.coef.normalized = B.coef / sum(B.coef)

              log.summaryA=paste("Forward step-wise regression model built (when both M.A and M.B markers) to predict:",cn[A]," includes these terms:",paste(A.who.cn,collapse=","),"with coefficients",paste(signif(A.coef,3),collapse=","));
              log.summaryB=paste("Forward step-wise regression model built (when both M.A and M.B markers) to predict:",cn[B]," includes these terms:",paste(B.who.cn,collapse=","),"with coefficients",paste(signif(B.coef,3),collapse=","));
              log.summary=paste(sep="\n",log.summaryA, log.summaryB,"\n",stringify(summary(step.A)),stringify(summary(step.B)),collapse="\n");

      # reassing M.A and M.B, so we can re-use the following code, inserting the weights....
  old.M.A = M.A
      old.M.B = M.B
      M.A = A.who.index.cn
      M.B = B.who.index.cn
      M.A.cn = cn[M.A]
      M.B.cn = cn[M.B]

      # INVAR: M.A has markers with stronger (>=) connection to A than B
      #    and M.B has markers with stronger connection to B than A

      # For now, let the SEM model fitting find the balance between markers into A
      # and for markers into B.

      z=local.sem.four.var.m1m2(pm=pm,MA.col=M.A,A.col=A,B.col=B,x,MB.col=M.B);
      twofer.eo[A,B]=z$eo.losem.lod;
      twofer.eo[B,A]=-twofer.eo[A,B];

      twofer.log[A,B]=paste(stringify(z),"\n",log.summary,collapse="\n");

      # and save the SNPs in column, so we can draw the graph

      # track in spreadsheet which SNPs we ended up using
      fsnp.cpa=paste(sep="","c(\"",paste(sep="",M.A.cn,collapse="\",\""),"\")")
      fsnp.oca=paste(sep="","c(\"",paste(sep="",M.B.cn,collapse="\",\""),"\")")      
      Final.SNPs.LEO.NB.OCA[A,B]=paste(sep="","fsnp.cpa=",fsnp.cpa,";fsnp.oca=",fsnp.oca);
      Final.SNPs.LEO.NB.OCA[B,A]=paste(sep="","fsnp.cpa=",fsnp.oca,";fsnp.oca=",fsnp.cpa);
      
      LEO.I[A,B]=z$LEO.I.AtoB;
      LEO.I[B,A]=z$LEO.I.BtoA;

      # new stats
      AB.vs.NextBest[A,B] = z$LEO.NB.AtoB
      AB.vs.NextBest[B,A] = z$LEO.NB.BtoA

      BLV[A,B] = z$BLV.AtoB
      BLV[B,A] = z$BLV.BtoA

      ZPathAB[A,B] = z$ZPathAB
      ZPathAB[B,A] = z$ZPathBA

      RMSEA.2m[A,B] = z$RMSEA.AtoB[1]
      RMSEA.2m[B,A] = z$RMSEA.BtoA[1]

      Minus.log.Model.P.value.AtoB[A,B]=z$mlogp.M.M1M2.AtoB
      Minus.log.Model.P.value.AtoB[B,A]=z$mlogp.M.M1M2.BtoA


   } # end i or A
} # end j or B

  pm.cat(pm,"\n");

  ret=list()
  ret$twofer.eo = twofer.eo # aka LEO.O
  ret$twofer.log = twofer.log
  ret$LEO.I = LEO.I # aka 1/LEO.I
  ret$BLV = BLV
  ret$ZPathAB = ZPathAB
  ret$AB.vs.NextBest = AB.vs.NextBest #aka LEO.NB
  ret$RMSEA.2m = RMSEA.2m
  ret$Final.SNPs.LEO.NB.OCA = Final.SNPs.LEO.NB.OCA
  ret$Minus.log.Model.P.value.AtoB = Minus.log.Model.P.value.AtoB   

  ret
}
# end losem.walk.two.steps.weighted.mean.mlreg
########################################


#
# helper function for obtain.marker.assignment.consistency()
#
if(exists("snps.on.same.chrom")) rm(snps.on.same.chrom);
snps.on.same.chrom=function(MA.snp.cn, other.snp.cn) {

      chrm1=sapply(strsplit(other.snp.cn,".chr",fixed=TRUE),function(x) x[2])
      chrm=sapply(strsplit(chrm1,".",fixed=TRUE),function(x) x[1])

      MA.chrm1=sapply(strsplit(MA.snp.cn,".chr",fixed=TRUE),function(x) x[2])
      MA.chrm=sapply(strsplit(MA.chrm1,".",fixed=TRUE),function(x) x[1])

      flag.same.chrom = rep(FALSE,length(other.snp.cn))

      # exclude remaining snps if they fall on the same chromosome as the forced M.A set
      for (i in 1:length(MA.snp.cn)) {
         flag.same.chrom[MA.chrm[i]==chrm] = TRUE
      }

   flag.same.chrom
}




# obtain.marker.assignment.consistency():
#
# snp selection function used in both 
#     losem.walk.two.steps.weighted.mean.mlreg
# and multione.losem.walk.two.steps.weighted.mean.mlreg
#
# INVAR: we have an edge from A to B possible.
#
# Make sure any overlapping belong only to the set with the strongest correlation...M.A or M.B.
# And make sure any SNP passes our minimum dependence threshold
#
# Update, 4 March 2008: we also now implement forcing of SNP assignment
#                       here, referring to pm$forced.MA.cn and pm$forced.MB.cn
#                       for the forced assignment of SNPs to A and B, if specified (not null).
#

if(exists("obtain.marker.assignment.consistency")) rm(obtain.marker.assignment.consistency);
obtain.marker.assignment.consistency=function(pm,snpcols,cx,cn, A, B) {

   z=list() # return this

      MAB= snpcols;

      M.A=M.B=c();
      top.A=top.B=NA;
      top.A.cx=top.B.cx=-2;

      # reduce this set as we assign snps
      MAB.remaining = MAB

      # try to get the chromsomes if possible
      got.chrom = FALSE # default
      chrm1=sapply(strsplit(cn[MAB],".chr",fixed=TRUE),function(x) x[2])
      chrm=sapply(strsplit(chrm1,".",fixed=TRUE),function(x) x[1])
      if (all(!is.na(chrm)) & all(chrm != "")) { got.chrom = TRUE; }


   #
   #  Handle forcing of markers in here too, so that even with forced M.A markers,
   #   the M.B markers can be less correlated, and vice-versa.
   # 
   # if we have either M.A and M.B forced...

   if (!is.null(pm$forced.MA.cn) | !is.null(pm$forced.MB.cn)) {

      # just MA forced....
      if (!is.null(pm$forced.MA.cn) & is.null(pm$forced.MB.cn)) {
          M.A = match(pm$forced.MA.cn, cn)
          if (any(is.na(M.A))) {
              print("got a problem: no M.A matches found for forced MA cols...dropping into debugger.");
              browser()
          }
          MAB.remaining = setdiff(MAB, M.A)

          if (got.chrom & pm$chromsome.wide.strong.marker.consistency) {
             MAB.remaining = MAB.remaining[!snps.on.same.chrom(MA.snp.cn=pm$forced.MA.cn, other.snp.cn=cn[MAB.remaining])]
         }

         for (m in MAB.remaining) {
             if (cx[m,B] >= cx[m,A]) { # still enforcing consistency...
                 if (cx[m,B] > pm$cor.dep.th) {                        
                       M.B=c(m,M.B)
                 }
             }
         }
      } 

      # just MB forced
      if (!is.null(pm$forced.MB.cn) & is.null(pm$forced.MA.cn)) {
         M.B = match(pm$forced.MB.cn, cn)
         MAB.remaining = setdiff(MAB, M.B)

          if (got.chrom & pm$chromsome.wide.strong.marker.consistency) {
             MAB.remaining = MAB.remaining[!snps.on.same.chrom(MA.snp.cn=pm$forced.MB.cn, other.snp.cn=cn[MAB.remaining])]
         } 

         for (m in MAB.remaining) {
             if (cx[m,A] >= cx[m,B]) { # still enforcing consistency...
                 if (cx[m,A] > pm$cor.dep.th) { # only care if its a significant snp correlation...
                    M.A=c(m,M.A)
                 }
             }
          }

      }

     # both MA and MB forced
     if (!is.null(pm$forced.MA.cn) & !is.null(pm$forced.MB.cn)) {
         M.A = match(pm$forced.MA.cn, cn)
         M.B = match(pm$forced.MB.cn, cn)
     }

      # now sort M.A and M.B so they are ordered according to absolute marginal correlation
      sort.a.cx = sort(cx[M.A,A],decreasing=TRUE,index.return=TRUE)
      sort.b.cx = sort(cx[M.B,B],decreasing=TRUE,index.return=TRUE)
      M.A = M.A[sort.a.cx$ix]
      M.B = M.B[sort.b.cx$ix]

      z$M.A=M.A
      z$M.B=M.B
      return(z)
   } # end if either marker set is forced.



###      #
###      # Marker assingment consistency principle implementation follows.
###      #


      if (!got.chrom | (got.chrom & !pm$chromsome.wide.strong.marker.consistency)) {
        for (m in MAB) {
         if (cx[m,A] >= cx[m,B]) { 
              if (cx[m,A] > pm$cor.dep.th) { 
                  M.A=c(m,M.A)
              }
         } else {
              if (cx[m,B] > pm$cor.dep.th) {
                  M.B = c(m,M.B)
              }
         }
        }

      } else { # use chrom info

         # INVAR: pm$chromsome.wide.strong.marker.consistency == TRUE & got.chrom == TRUE

         # break ties on the same chromsome first
         chrm.table = table(chrm)  # group by same chromosome
         competing.chrom = names(which(chrm.table>1))
         if (length(competing.chrom)>0) {
             # have some ties to break
             for (i in 1:length(competing.chrom)) {
                  chrom.with.tie = competing.chrom[i]

                  snps.with.tie = MAB[chrom.with.tie==chrm]

                  # does the chrom go to A or B? which has the stronger cors?
                  top.cor.A = max(cx[snps.with.tie,A])
                  top.cor.B = max(cx[snps.with.tie,B])

                  w.top.cor.A = which.max(cx[snps.with.tie,A])
                  w.top.cor.B = which.max(cx[snps.with.tie,B])
                  snp.top.cor.A = snps.with.tie[w.top.cor.A]
                  snp.top.cor.B = snps.with.tie[w.top.cor.B]

                  # take out the snps on this chrom, leaving the remaining ones for next
                  MAB.remaining = setdiff(MAB.remaining, snps.with.tie)

                  if (top.cor.A > top.cor.B) {
                      m=snp.top.cor.A
                      if (cx[m,A] > pm$cor.dep.th) { 
                         M.A=c(m,M.A)
                      }
                  } else {
                      m = snp.top.cor.B
                      if (cx[m,B] > pm$cor.dep.th) {
                         M.B = c(m,M.B)
                      }
                  }


             } # end i loop over all repeated chromosomes...

         } # end if length(competing.chrom) > 0
         # all remaining snps in MAB.remaining are from distinct chromosomes

         for (m in MAB.remaining) {
            if (cx[m,A] >= cx[m,B]) { 
                 if (cx[m,A] > pm$cor.dep.th) { 
                     # simultaneously semi-sort so strongest correlation is first...
                     if (cx[m,A] > top.A.cx) { top.A = m; top.A.cx=cx[m,A]; M.A=c(m,M.A); } else { M.A=c(M.A,m); } 
                 }
            } else {
                 if (cx[m,B] > pm$cor.dep.th) {
                     # simultaneously semi-sort so strongest correlated marker is frist
                     if (cx[m,B] > top.B.cx) { top.B = m; top.B.cx=cx[m,B]; M.B = c(m,M.B);} else { M.B = c(M.B,m); }
                 }
            }
        } # end for m in MAB.remaining
         
      } # end else use chrom info


      # now sort M.A and M.B so they are ordered according to absolute marginal correlation
      sort.a.cx = sort(cx[M.A,A],decreasing=TRUE,index.return=TRUE)
      sort.b.cx = sort(cx[M.B,B],decreasing=TRUE,index.return=TRUE)
      M.A = M.A[sort.a.cx$ix]
      M.B = M.B[sort.b.cx$ix]

      if (!is.null(pm$max.snps.in.OCA.CPA.sets) & pm$max.snps.in.OCA.CPA.sets > 0) {
         a.endpt = min(length(M.A),pm$max.snps.in.OCA.CPA.sets)
         b.endpt = min(length(M.B),pm$max.snps.in.OCA.CPA.sets)

         if (length(M.A) > 0) { M.A = M.A[1:a.endpt]} else {M.A=c()}
         if (length(M.B) > 0) { M.B = M.B[1:b.endpt]} else {M.B=c()}
      }

   z$M.A=M.A
   z$M.B=M.B
   z
}


# Testing forced MA assignment on the MALE data, using the female best markers.
if (exists("unit.test.forced.marker.assignment")) { rm(unit.test.forced.marker.assignment) }
unit.test.forced.marker.assignment=function() { 

  a=load("/home/jaten/dev/peculiar/bxh/liver.1146snps.23388mrna.21clinical.bxh.apoe.null.rdat")
  liver.bxh.male.female=liver.bxh.male.female.cm
  x=liver.bxh.male.female.cm

  female=liver.bxh.male.female$sex==2
  male=liver.bxh.male.female$sex==1
#  datCombined=liver.bxh.male.female[female,-str.me$sex]
  datCombined=liver.bxh.male.female[male,-str.me$sex]

  pm=neo.get.param()  

  which.Insig1=pmatch("Insig1",colnames(liver.bxh.male.female))
  which.Fdft1=pmatch("Fdft1.MMT00082285",colnames(liver.bxh.male.female))
  which.Dhcr7=pmatch("Dhcr7",colnames(liver.bxh.male.female))
  which.B4301=pmatch("B430110G05Rik",colnames(liver.bxh.male.female))
  
  pm$run.title="FORCED.male.chr16.17.8.Insig1"

  pm$forced.MA.colnum=match(c("rs3670293.chr8.bp100063765.Morgans0.341901","rs3722983.chr17.bp7478326.Morgans0.019370","rs3705921.chr16.bp47999105.Morgans0.236810"),colnames(datCombined))

   one.date=date() 

  my.new.dir=gsub(":",".",gsub(" ","_",paste(sep="","seed.insig1.",one.date)))
  dir.create(my.new.dir)
  setwd(my.new.dir)
  
  seed.root.leaf.neo(pm, datCombined, snpcols=str.me$snpcols, known.upstream=which.Insig1, known.downstream=which.Fdft1, genecols=str.me$genecols, require.sma.model.prob.over=.5, require.sma.leo.nb.over=.75, num.seed.snps=2, require.cor.with.known.up.over=.5) 

  
}




#
# multione.losem.walk.two.steps.weighted.mean.mlreg
#
if(exists("multione.losem.walk.two.steps.weighted.mean.mlreg") ) rm(multione.losem.walk.two.steps.weighted.mean.mlreg);
multione.losem.walk.two.steps.weighted.mean.mlreg=function(x,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log,zeo.proxy=NA) {

   pm.print(pm,paste("Starting multione.losem.walk.two.steps.mean.mlreg()...",date()));

   # allocate log if not already
   if (length(twofer.log)==0) { twofer.log=vector("list", length(non.snp.cols)*2) } # create empty list of given length

   # initializer for LEO.I, BLV, ZPathAB, AB.vs.NextBest
   na.init=twofer.eo; na.init[]=NA;

   # track the relative probability of the confounded model
   LEO.I = na.init

   # track the newly added statistics, initialized to NA
   BLV = na.init
   ZPathAB = na.init
   AB.vs.NextBest = na.init


   for (i in 1:(length(non.snp.cols)-1)) {
      for (j in (i+1):length(non.snp.cols)) {
      A= non.snp.cols[i];
      B= non.snp.cols[j];

      # use the A,B edge filters if present
      if (!is.null(pm$A) & !is.A.B.edge(cn[A],cn[B],pm)) next;

      if(edge.is.on.ignore.list(A,B,pm)) { twofer.log[A,B]="Ignored edge b/c on ignore list";  next; }
      if(!is.na(pm$zeo.proxy.for.leo.th) && !is.na(zeo.proxy[A,B])) {
         if (abs(zeo.proxy[A,B]) < pm$zeo.proxy.for.leo.th) { twofer.log[A,B]=paste("Edge eliminated by zeo.proxy score of",signif(zeo.proxy[A,B],3),"with elimination th:",pm$zeo.proxy.for.leo.th);  next; }
      }

         pm.cat(pm,paste(sep="",cn[A],"->",cn[B],"...")); # give a sense of progress

      if (af[A,B]==0 && af[B,A]==0) next; 
      # INVAR: we have an edge from A to B possible.

      # Make sure any overlapping belong only to the set with the strongest correlation...M.A or M.B.
      # And make sure any SNP passes our minimum dependence threshold.

      mac=obtain.marker.assignment.consistency(pm, snpcols, cx, cn, A, B)
      M.A = mac$M.A
      M.B = mac$M.B

     if (length(M.A) == 0 && length(M.B) == 0) next;

           ########################################
           if (length(M.A) > 0) {
       
              #### multi-marker into A only

              # do a single marker score, one two step, or half a ZEO: M.A -> A -> B
              Mset=M.A;
              parent=A;
              child=B;
              singleAB=1;
              step.1 = forward.step.multimarker.twostep(x,cn[A],cn[M.A],cn[M.A][1]);

              A.coef = step.1$percent.of.explained.variance.explained;
              A.who.cn  = names(step.1$coefficients)[-1]; # names
              A.who.index  = match(A.who.cn,cn); # get column index within Mset.cn

              log.summary=paste("Forward step-wise regression model built to predict:",cn[parent]," includes these terms:\n",stringify(A.who.cn),"with coefficients",stringify(A.coef));

              z=compare.local.sems(pm=pm,A.who.index,parent,child,x);
              twofer.eo[A,B]=z$LEO.O.AtoB
              twofer.log[A,B]=paste(stringify(z),"\n",log.summary);

              LEO.I[A,B]=z$LEO.I.AtoB;
              AB.vs.NextBest[A,B] = z$LEO.NB.AtoB
        BLV[A,B] = z$BLV.AtoB
              ZPathAB[A,B] = z$ZPathAB
         }
      
         ###################################################

         #### multi-marker into B only

         if (length(M.B) > 0) {
              # do a single marker score, one two step, or half a ZEO: M.B -> B -> A
              Mset=M.B
              parent=B;
              child=A;
              singleAB=-1;
              step.1 = forward.step.multimarker.twostep(x,cn[B],cn[M.B],cn[M.B][1]);

              B.coef = step.1$percent.of.explained.variance.explained;
              B.who.cn  = names(step.1$coefficients)[-1]; # names
              B.who.index  = match(B.who.cn,cn); # get column index within Mset.cn

              log.summary=paste("Forward step-wise regression model built to predict:",cn[parent]," includes these terms:\n",stringify(B.who.cn),"with coefficients",stringify(B.coef));

              z=compare.local.sems(pm=pm,B.who.index,parent,child,x);
              twofer.eo[B,A]=z$LEO.O.AtoB
              twofer.log[B,A]=paste(stringify(z),"\n",log.summary)

              LEO.I[B,A]=z$LEO.I.AtoB
              AB.vs.NextBest[B,A] = z$LEO.NB.AtoB
        BLV[B,A] = z$BLV.AtoB
              ZPathAB[B,A] = z$ZPathAB
            }
      
      } # end i or A
   } # end j or B

  pm.cat(pm,"\n");

  ret=list()
  ret$twofer.eo = twofer.eo # aka LEO.O
  ret$twofer.log = twofer.log
  ret$LEO.I = LEO.I # aka 1/LEO.I
  ret$BLV = BLV
  ret$ZPathAB = ZPathAB
  ret$AB.vs.NextBest = AB.vs.NextBest #aka LEO.NB

  ret
}
# end multione.losem.walk.two.steps.weighted.mean.mlreg
########################################

########################################
#
# Steve requested this  score that fits all vs. all m1m2 models and then averages...
#
########################################

if(exists("m1m2.average.losem.walk.two.steps") ) rm(m1m2.average.losem.walk.two.steps);
m1m2.average.losem.walk.two.steps=function(x,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log,zeo.proxy=NA) {

   pm.print(pm,paste("Starting  m1m2.average.losem.walk.two.steps.mean.mlreg()...",date()));

   # allocate log if not already
   if (length(twofer.log)==0) { stop("Must pass valid twofer.log into function.") } # create empty list of given length

   # initialize
   twofer.eo[]=NA

   for (i in 1:(length(non.snp.cols)-1)) {
      for (j in (i+1):length(non.snp.cols)) {
      A= non.snp.cols[i];
      B= non.snp.cols[j];

      # use the A,B edge filters if present
      if (!is.null(pm$A) & !is.A.B.edge(cn[A],cn[B],pm)) next;

      if(edge.is.on.ignore.list(A,B,pm)) { twofer.log[A,B]="Ignored edge b/c on ignore list";  next; }
      if(!is.na(pm$zeo.proxy.for.leo.th) && !is.na(zeo.proxy[A,B])) {
         if (abs(zeo.proxy[A,B]) < pm$zeo.proxy.for.leo.th) { twofer.log[A,B]=paste("Edge eliminated by zeo.proxy score of",signif(zeo.proxy[A,B],3),"with elimination th:",pm$zeo.proxy.for.leo.th);  next; }
      }

         pm.cat(pm,paste(sep="",cn[A],"->",cn[B],"... at ",date(), "...\n")); # give a sense of progress

      if (af[A,B]==0 && af[B,A]==0) next; 
      # INVAR: we have an edge from A to B possible.

      # Make sure any overlapping belong only to the set with the strongest correlation...M.A or M.B.
      # And make sure any SNP passes our minimum dependence threshold.
      MAB= snpcols;
      M.A=M.B=c();
      top.A=top.B=NA;
      top.A.cx=top.B.cx=-2;

      # keep the correlations too, so we can sort by strongest absolute correlation, and order M.A and M.B
      M.A.abs.cor=c()
      M.B.abs.cor=c()

      for (m in MAB) {
         if (cx[m,A] >= cx[m,B]) { 
              if (cx[m,A] > pm$cor.dep.th) { 
                  
                  # simultaneously semi-sort so strongest correlation is first...
                  if (cx[m,A] > top.A.cx) { 
                      top.A = m; top.A.cx=cx[m,A]; 
                      M.A=c(m,M.A); 
                      M.A.abs.cor = c(abs(cx[m,A]),M.A.abs.cor)
                  } else { 
                      M.A=c(M.A,m); 
                      M.A.abs.cor = c(M.A.abs.cor,abs(cx[m,A]))
                  } 
              }
         } else {
              if (cx[m,B] > pm$cor.dep.th) {
                  # simultaneously semi-sort so strongest correlated marker is frist
                  if (cx[m,B] > top.B.cx) { 
                       top.B = m; 
                       top.B.cx=cx[m,B]; 
                       M.B = c(m,M.B);
                       M.B.abs.cor = c(abs(cx[m,B]),M.B.abs.cor)
                   } else { 
                       M.B = c(M.B,m)
                       M.B.abs.cor = c(M.B.abs.cor,abs(cx[m,B]))
                   }
              }
         }
      }

      if (length(M.A) == 0 || length(M.B) == 0) { # require markers into both A and B.
          twofer.log[A,B] = "Didn't find both M.A and M.B sets, so stopping."
          twofer.log[B,A] = twofer.log[A,B]
          next
      }

      # sort by M.A.abs.cor so that M.A is ordered from strongest correlation to weakest
      ma.sort.order = sort(M.A.abs.cor,decreasing=T,index.return=T)$ix
      M.A.abs.cor = M.A.abs.cor[ma.sort.order]
      M.A = M.A[ma.sort.order]
      names(M.A.abs.cor) = cn[M.A]

      # same for M.B.abs.cor
      mb.sort.order = sort(M.B.abs.cor,decreasing=T,index.return=T)$ix
      M.B.abs.cor = M.B.abs.cor[mb.sort.order]
      M.B = M.B[mb.sort.order]
      names(M.B.abs.cor) = cn[M.B]

      # for logging
      abs.correlations = list()
      abs.correlations$M.A.abs.cor = M.A.abs.cor
      abs.correlations$M.B.abs.cor = M.B.abs.cor
      abs.cor.txt = stringify(abs.correlations)

      # now...both markers...
      M.A.cn = cn[M.A]
      M.B.cn = cn[M.B]

      # INVAR: M.A has markers with stronger (>=) connection to A than B
      #    and M.B has markers with stronger connection to B than A

      # For now, let the SEM model fitting find the balance between markers into A
      # and for markers into B.

      average.orthomarker = 0
      average.orthomarker.BA = 0
      a.o.count = 0
      total.count = length(M.A) * length(M.B)
      log.vec = rep(NA,total.count)

      m1m2.score.AB.labels = rep("",total.count)
      m1m2.score.BA.labels = rep("",total.count)
      m1m2.scores.AB = rep(NA,total.count)
      m1m2.scores.BA = rep(NA,total.count)

      for (one.marker.a in M.A) {
         for (one.marker.b in M.B) {
            z=local.sem.four.var.m1m2(pm=pm,MA.col=one.marker.a,A.col=A,B.col=B,x,MB.col=one.marker.b);
            a.o.count = a.o.count + 1
            average.orthomarker = average.orthomarker + z$LEO.NB.AtoB
            average.orthomarker.BA = average.orthomarker.BA + z$LEO.NB.BtoA

            log.vec[a.o.count] = stringify(z)

            m1m2.scores.AB[a.o.count] = signif(z$LEO.NB.AtoB,3)
            m1m2.scores.BA[a.o.count] = signif(z$LEO.NB.BtoA,3)
            m1m2.score.AB.labels[a.o.count] = paste(sep="","A(",cn[A],").to.B(",cn[B],").with.MA(",cn[one.marker.a],").and.MB(",cn[one.marker.b],")")
            m1m2.score.BA.labels[a.o.count] = paste(sep="","A(",cn[B],").to.B(",cn[A],").with.MA(",cn[one.marker.b],").and.MB(",cn[one.marker.a],")")
         }
      }

      #names(m1m2.scores.AB) = m1m2.score.AB.labels
      #names(m1m2.scores.BA) = m1m2.score.BA.labels

      m1m2.scores.AB.stringified = paste(sep="","m1m2.scores.AB = c(",paste(collapse=",",m1m2.scores.AB),")")
      m1m2.scores.BA.stringified = paste(sep="","m1m2.scores.BA = c(",paste(collapse=",",m1m2.scores.BA),")")

      m1m2.labels.AB  = paste(sep="","m1m2.labels.AB = c(\"",paste(collapse="\",\"",m1m2.score.AB.labels),"\")")
      m1m2.labels.BA  = paste(sep="","m1m2.labels.BA = c(\"",paste(collapse="\",\"",m1m2.score.BA.labels),"\")")

      details.AB = paste(m1m2.labels.AB, m1m2.scores.AB.stringified,"\n", sep="\n",collapse="\n")
      details.BA = paste(m1m2.labels.BA, m1m2.scores.BA.stringified,"\n", sep="\n",collapse="\n")

      details = paste(details.AB, details.BA, abs.cor.txt, collapse="\n\n",sep="\n\n")

      average.orthomarker = average.orthomarker / a.o.count
      average.orthomarker.BA = average.orthomarker.BA / a.o.count

      twofer.eo[A,B]= average.orthomarker
      twofer.eo[B,A]= average.orthomarker.BA

      preamble.txt = paste(paste("M.A markers were: ",cn[M.A],collapse="\n"),"\n",paste("M.B markers were: ",cn[M.B],collapse="\n"),"\n",details,sep="\n",collapse="\n")

      title.om.avg = paste("\n==== NEXT PAIR of MARKERS ===== In the full log for the M1.M2.Average over ",total.count," orthomarker model(s) from: ",cn[A]," to ",cn[B])

      A.to.B.intro.text = paste("M1M2.Average LOG for A(",cn[A],") -> B(",cn[B],") with M1M2.AVG score of ",signif(twofer.eo[A,B],3),"\n")
      B.to.A.intro.text = paste("M1M2.Average LOG for B(",cn[B],") -> A(",cn[A],") with M1M2.AVG score of ",signif(twofer.eo[B,A],3),"\n")

      twofer.log[A,B]=paste(A.to.B.intro.text,"\n",preamble.txt,paste(title.om.avg,"\n\n",log.vec,"\n",collapse="\n"),sep="\n",collapse="\n")
      twofer.log[B,A]=paste(B.to.A.intro.text,"\n",preamble.txt,paste(title.om.avg,"\n\n",log.vec,"\n",collapse="\n"),sep="\n",collapse="\n")
 
   } # end i or A
} # end j or B


  pm.cat(pm,"\n");

  ret=list()
  ret$M1M2.AVG = twofer.eo
  ret$M1M2.AVG.log = twofer.log

  ret
}
# end m1m2.average.losem.walk.two.steps
########################################





# zlist.mean: helper function for losem.walk.two.steps.one.permutation.one.vote
# EFFECT: given a list of named lists, and a fieldname good for all the contained lists,
# compute the mean of the field over all the contained lists.
#
if(exists("zlist.mean") ) rm(zlist.mean)
zlist.mean=function(zlist, fieldname) {
   res=0;
   for (i in 1:length(zlist)) {
     a=zlist[[i]];
     res=res+ a[[which(names(a)==fieldname)]];
   }
   res/length(zlist)
 }

########################################
# "All vote" a.k.a. "Democracy" a.k.a. "One permutation, One vote."
# LOCAL SEM VERSION: do the model fitting using compare.local.sems() for
#                    each of the SNPs, then average at the end. 
#  One permutation-one vote two-step strategy 
#  lets all marker combinations vote, with equal weighting
#  
#
# Possible Todo later: Implement model where all the SNPs go into one huge SEM model
#   and fit that??? ie. 50 SNPs -> A -> B <- 50 SNPs  VS.  50 SNPs -> A <- B <- 50 SNPs
# This procedure is not that. Instead, it computes the scores for each 3 variable
# model (using a different SNP each time) seperately, then averages at the end.
#
########################################
if(exists("losem.walk.two.steps.one.permutation.one.vote") ) rm(losem.walk.two.steps.one.permutation.one.vote);
losem.walk.two.steps.one.permutation.one.vote=function(x,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log,zeo.proxy=NA) {

pm.print(pm,paste("Starting losem.walk.two.steps.one.permutation.one.vote()...",date()));

# allocate log if not already
if (length(twofer.log)==0) { twofer.log=vector("list", length(non.snp.cols)*2) } # create empty list of given length

   # initializer for LEO.I, BLV, ZPathAB, AB.vs.NextBest
   na.init=twofer.eo; na.init[]=NA;

   # track the relative probability of the confounded model
   LEO.I = na.init

   # track the newly added statistics, initialized to NA
   BLV = na.init
   ZPathAB = na.init
   AB.vs.NextBest = na.init

   # return NULL if not doing this function
   if (!pm$do.leo.all) {
       ret=list()
       ret$twofer.eo = twofer.eo # aka LEO.O
       ret$twofer.log = twofer.log
       ret$LEO.I = LEO.I
       ret$BLV = BLV
       ret$ZPathAB = ZPathAB
       ret$AB.vs.NextBest = AB.vs.NextBest #aka LEO.NB
       return(ret)
   }

for (i in 1:(length(non.snp.cols)-1)) {
   for (j in (i+1):length(non.snp.cols)) {
      A= non.snp.cols[i];
      B= non.snp.cols[j];

      # use the A,B edge filters if present
      if (!is.null(pm$A) & !is.A.B.edge(cn[A],cn[B],pm)) next;

      if(edge.is.on.ignore.list(A,B,pm)) { twofer.log[A,B]="Ignored edge b/c on ignore list";  next; }
      if(!is.na(pm$zeo.proxy.for.leo.th) && !is.na(zeo.proxy[A,B])) {
         if (abs(zeo.proxy[A,B]) < pm$zeo.proxy.for.leo.th) { twofer.log[A,B]=paste("Edge eliminated by zeo.proxy score of",signif(zeo.proxy[A,B],3),"with elimination th:",pm$zeo.proxy.for.leo.th);  next; }
      }

      pm.cat(pm,paste(sep="",cn[A],"->",cn[B],"...")); # give a sense of progress

      if (af[A,B]==0 && af[B,A]==0) next; 
      # INVAR: we have an edge from A to B possible.
      M.A = setdiff(which(af[,A]==1),non.snp.cols); # the Markers into A
      M.B = setdiff(which(af[,B]==1),non.snp.cols); # the Markers into B

      # make sure any overlapping belong only to the set with the strongest correlation...
      MAB=union(M.A,M.B);
      M.A=M.B=c();
      for (m in MAB) {
         if (cx[m,A] >= cx[m,B]) { M.A=c(M.A,m); } else { M.B = c(M.B,m); }
      }

      Mset=c(); # grandparent candidates
      parent=NA;
      child=NA;
      singleAB=1; # 1 for M->A->B; -1 for M->B->A
      if (length(M.A) == 0 || length(M.B) == 0)
      {
          if(length(M.A) > 0) {  
              # do a single marker score, one two step, or half a ZEO: M.A -> A -> B
              Mset=M.A;
              parent=A;
              child=B;

              # now compute A-B edge orientation using just the one marker
              # space for results:
              m1ab = cx[Mset,parent]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[Mset];}
              zlist = list() # store the z's, so we can compute the rest of the new summary stats too
              m1ab.str=m1ab; m1ab.str[]=NA;

              for (k in 1:length(Mset)) {
                  m1 = Mset[k];
                  z=zlist[[k]]=compare.local.sems(pm=pm,m1,parent,child,x);
                  m1ab[k] = z$eo.losem.lod;
                  m1ab.str[[k]]=stringify(z);
              }

              mean.m1ab = mean(m1ab);
  
              twofer.eo[A,B]=mean.m1ab;
              twofer.eo[B,A]=-twofer.eo[A,B];
              twofer.log[A,B]=stringify.cat(paste(sep="\n",m1ab.str,collapse="\n(Mean component)\n"));

              # new stats 
              LEO.I[A,B]=zlist.mean(zlist,"LEO.I.AtoB")
              LEO.I[B,A]=zlist.mean(zlist,"LEO.I.BtoA")
              
              AB.vs.NextBest[A,B] = zlist.mean(zlist,"LEO.NB.AtoB")
              AB.vs.NextBest[B,A] = zlist.mean(zlist,"LEO.NB.AcollideB")
              
              BLV[A,B] = zlist.mean(zlist,"BLV.AtoB")
              BLV[B,A] = zlist.mean(zlist,"BLV.BtoA")
              
              ZPathAB[A,B] = zlist.mean(zlist,"ZPathAB")
              ZPathAB[B,A] = zlist.mean(zlist,"ZPathBA")
              
              } else if (length(M.B) > 0) {

                  # interesting situation here. We have no markers into A.
                  # We have marker(s) into B. We wish to return a score
                  # for the A->B edge into twofer.eo[A,B] that reflects
                  # M->B<-A situation. Hence we added a score to compare.local.sems()
                  # to reflect this situation; we want to use here 
                  #
                  # compare.local.sems()$LEO.NB.AcollideB  (which is $LEO.A'collideB' vs next best: when A' is our B; and B' is our A)
                  #
                  #     rather than
                  #
                  # compare.local.sems()$eo.losem.lod
                  #

                  Mset=M.B
                  parent=B;
                  child=A;

                  # now compute A-B edge orientation using just the one marker
                  # space for results:
                  m1ab = cx[Mset,parent]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[Mset];}
                  m1ab.str=m1ab; m1ab.str[]=NA;
                  zlist = list() # store the z's, so we can compute the rest of the new summary stats too

                  for (k in 1:length(Mset)) {
                      m1 = Mset[k];
                      z=zlist[[k]]=compare.local.sems(pm=pm,m1,parent,child,x);
                      m1ab[k] = z$LEO.NB.AcollideB;
                      m1ab.str[[k]]=stringify.cat(paste("\n\n|||*** Using LEO.NB.AcollideB (REVERSED) ***|||\n",stringify(z)));
                  }

                  mean.m1ab = mean(m1ab);
  
                  twofer.eo[A,B]=mean.m1ab;
                  twofer.eo[B,A]=-twofer.eo[A,B];
                  twofer.log[A,B]=stringify.cat(paste(sep="\n",stringify.cat(m1ab.str),collapse="\n(Mean REV component)\n"));

                  # new stats...punt
                  LEO.I[A,B]= zlist.mean(zlist,"mlogp.M.conf") - zlist.mean(zlist,"mlogp.M.AcollideB")
                  LEO.I[B,A]= zlist.mean(zlist,"LEO.I.AtoB")
              
                  AB.vs.NextBest[A,B] = twofer.eo[A,B]
                  AB.vs.NextBest[B,A] = -twofer.eo[A,B]
              
                  BLV[A,B] = NA #zlist.mean(zlist,"BLV.AtoB")
                  BLV[B,A] = NA #zlist.mean(zlist,"BLV.BtoA")
                  
                  ZPathAB[A,B] = NA #zlist.mean(zlist,"ZPathAB")
                  ZPathAB[B,A] = NA #zlist.mean(zlist,"ZPathBA")

              } else next;

          next; # b/c rest of code is for markers into both A and B
      }

      if (length(M.A) == 0 || length(M.B) == 0) next; # require markers into both A and B.

      # INVAR: M.A has markers with stronger (>=) connection to A than B
      #    and M.B has markers with stronger connection to B than A

      # compute the Z.2steps
      # reserve space, and get the names() right
      m1ab = cx[M.A,A]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[M.A];}
      m2ba = cx[M.B,B]; m2ba[]=NA; if (length(m2ba)==1) { names(m2ba)=rownames(cx)[M.B];}

      # store the string output too
      m1ab.str=m1ab;
      m2ba.str=m2ba;

      zlist.ma = list()
      zlist.mb = list()
      
      for (k in 1:length(M.A)) {
           m1 = M.A[k];
           zlist.ma[[k]]=z=compare.local.sems(pm=pm,m1,A,B,x);
           m1ab[k] = z$eo.losem.lod;
           m1ab.str[k]=stringify(z);
           m1ab.str[k]=paste(sep="","All vote component log, direction Forward ___ M.A(",cn[m1],") -> A(",cn[A],") -> B(",cn[B],") score used here after: ___\n",stringify(z));
      }

      for (k in 1:length(M.B)) {
           m2 = M.B[k];
           zlist.mb[[k]]=z=compare.local.sems(pm=pm,m2,B,A,x);
           m2ba[k] = z$eo.losem.lod;
           m2ba.str[k]=paste(sep="","All vote component log, direction Backwards ___ M.B(",cn[m2],") -> B(",cn[B],") -> A(",cn[A],") score used here after: ___\n",stringify(z));
      }

      # now do the weighting
      z.mean.a.b = mean(m1ab);
      z.mean.b.a = mean(m2ba);
      z.total = z.mean.a.b - z.mean.b.a;

        for.log1=paste("(FWD: mean of ",stringify(m1ab),")=",signif(z.mean.a.b,3))
        for.log2=paste("(BACK: mean of ",stringify(m2ba),")=",signif(z.mean.b.a,3))
        for.log=paste(sep="\n minus \n",for.log1,for.log2)
        log.entry.summary=paste(for.log,"=TOTAL:",signif(z.total,3),"\n");

         twofer.eo[A,B]=z.total;
         twofer.eo[B,A]=-twofer.eo[A,B];
         twofer.log[A,B]=stringify.cat(paste("All vote log for: A(",cn[A],") -> B(",cn[B],")\n",paste(collapse="(NEXT FWD LOG)\n",m1ab.str),paste(collapse="(NEXT BACKWARDS LOG)\n",m2ba.str),log.entry.summary));


      # new stats 
      LEO.I[A,B]=zlist.mean(zlist.ma,"LEO.I.AtoB") - zlist.mean(zlist.mb,"LEO.I.AtoB")
      LEO.I[B,A]= - LEO.I[A,B]

      AB.vs.NextBest[A,B] = zlist.mean(zlist.ma,"LEO.NB.AtoB") -  zlist.mean(zlist.mb,"LEO.NB.AtoB")
      AB.vs.NextBest[B,A] = -AB.vs.NextBest[A,B]

      BLV[A,B] = zlist.mean(zlist.ma,"BLV.AtoB") -  zlist.mean(zlist.mb,"BLV.AtoB")
      BLV[B,A] = zlist.mean(zlist.ma,"BLV.BtoA") -  zlist.mean(zlist.mb,"BLV.BtoA")

      ZPathAB[A,B] = zlist.mean(zlist.ma,"ZPathAB") -  zlist.mean(zlist.mb,"ZPathAB")
      ZPathAB[B,A] = zlist.mean(zlist.ma,"ZPathBA") -  zlist.mean(zlist.mb,"ZPathBA")

   } # end i or A
} # end j or B

pm.cat(pm,"\n");
#old: list(twofer.eo,twofer.log)

# old  ret=list(twofer.eo,twofer.log,LEO.I)

  ret=list()
  ret$twofer.eo = twofer.eo # aka LEO.O
  ret$twofer.log = twofer.log
  ret$LEO.I = LEO.I
  ret$BLV = BLV
  ret$ZPathAB = ZPathAB
  ret$AB.vs.NextBest = AB.vs.NextBest #aka LEO.NB

  ret


}
# end losem.walk.two.steps.one.permutation.one.vote
########################################



#
#  stepAIC can't efficiently handle itself on the big data, take 2x number of steps
#   say 10 instead of 5 here, then let stepAIC find it's own cutoff using AIC
#
if(exists("manual.stepwise.regression") ) rm(manual.stepwise.regression);
manual.stepwise.regression=function(yname,snpcolnames,x,steps=10) {
   if (steps > ncol(x)) { steps = ncol(x)-1; }

   # do this up front now....
   # fill in missing data, or else we crash; SNPs are already imputed by impute.nn, so this is just for y
   #x = as.data.frame(impute(x,what="median"))

   cn=colnames(x);
   ynum = which(cn == yname);
   snpcols.num = match(snpcolnames,cn)
   ynew = x[,ynum];
   found=c();
   found.num=c();

   for (i in 1:steps) {
      cor.y.snp = cor(ynew,x[,snpcols.num],method="pearson",use="pairwise.complete.obs");
      wmax = which.max(abs(cor.y.snp))
      found=c(found,snpcolnames[wmax]);
      found.num=c(found.num,snpcols.num[wmax])
      
      # get the residuals after wmax is put into the model
      df=data.frame(ynew,x[,found.num]); # put impute here
      ynew = lm(ynew~.,data=df)$residuals;
    }
  cn[unique(found.num)]
}
# test call: msr=manual.stepwise.regression(yname="PC2",snpcolnames=colnames(x)[c(1,2,3,4,5,6)],x=x,steps=10)

#
# may be wise to limit the number of steps to 5-7 or so (to only pick 5-7 variables--b/c otherwise
#  the sem/nlm optimization routines can get overwhelmed.
#
# Use forward step-wise regression to estimate a model for the impact of multiple markers on A
if(exists("forward.step.multimarker.twostep") ) rm(forward.step.multimarker.twostep);
forward.step.multimarker.twostep=function(df,ycol,snpcols,single.start.col,nsteps=5,view.step=F) {
   new.snpcols = manual.stepwise.regression(ycol,snpcols,df,steps=nsteps)

   # take care of wierd bug situation when single.start.col is not chosen by
   # the greedy manual.stepwise.regression function....but single.start.col
   # must also be in the f.upper(in addition to always being the only term in f.lower)
   # or else stepAIC crashes.
   if (!(single.start.col %in% new.snpcols)) { new.snpcols = c(single.start.col, new.snpcols) }

   f.upper <- as.formula(paste(ycol, "~", paste(new.snpcols, collapse = " + ")))
   f.lower <- as.formula(paste(ycol, "~", single.start.col))

   f0 <- lm(f.lower, data.frame(df))
   f.step <- stepAIC(f0, list(upper = f.upper, lower = f.lower), direction="forward",trace = F, steps=nsteps)
   if (view.step) {
      print("===========================================================");
      print(paste(sep="","SNP forward stepwise selection for non-snp: ",ycol))
      print("===========================================================");
      print(paste(sep="","f.upper was: ",stringify(f.upper)));
      print(paste(sep="","f.lower was: ",stringify(f.lower)));
      print(f.step$anova)
   }

   # To weight each Z-score according to the percentage of variance explained
   # we need to find the residual sum of squares for variable in the model: anova.
   #  ....therefore...
   # compute percentage of variance for each chosen variable
   anova.step.A=anova(f.step);
   sum.sq.step.A=anova.step.A$"Sum Sq"
   len.sum.sq.step.A=length(sum.sq.step.A);
   total.sum.A.without.residuals = sum(sum.sq.step.A[1:(len.sum.sq.step.A-1)]);
   A.percent.explained.var = sum.sq.step.A[1:(len.sum.sq.step.A-1)]/total.sum.A.without.residuals

   names(A.percent.explained.var)=rownames(anova.step.A)[1:(len.sum.sq.step.A-1)]

   f.step$percent.of.explained.variance.explained=A.percent.explained.var

   f.step
}


# pick SNPs by forward step-wise regression...
if(exists("forward.step.pick.snp") ) rm(forward.step.pick.snp);
forward.step.pick.snp=function(df,ycol,snpcols,single.start.col,nsteps=5,view.step=F) {

   if (nsteps < 1) stop("must call forward.step.pick.snp with at least 1 snp requested!")

   new.snpcols = manual.stepwise.regression(ycol,snpcols,df,steps=nsteps*3)

   # take care of wierd bug situation when single.start.col is not chosen by
   # the greedy manual.stepwise.regression function....but single.start.col
   # must also be in the f.upper(in addition to always being the only term in f.lower)
   # or else stepAIC crashes.
   if (!(single.start.col %in% new.snpcols)) { new.snpcols = c(single.start.col, new.snpcols) }

   f.upper <- as.formula(paste(ycol, "~", paste(new.snpcols, collapse = " + ")))
   f.lower <- as.formula(paste(ycol, "~", single.start.col))

   f0 <- lm(f.lower, df)
   f.step <- stepAIC(f0, list(upper = f.upper, lower = f.lower), direction="forward",trace = F, steps=nsteps)

   if (view.step) {
      print("===========================================================");
      print(paste(sep="","SNP forward stepwise selection for non-snp: ",ycol))
      print("===========================================================");
      print(f.step$anova)
   }
   chosen.snps=names(f.step$coefficients)[-1]
}


#
# zeo2(): Centralize ZEO computation: so we can alter/play with it, if so desired.
#
# N.B.: the four display.Zm1B() like functions have the original ZEO hard coded in them.
#
if(exists("zeo2") ) rm(zeo2);
zeo2=function(gp,gk,pa,covx,pm) {
   z=list();
   z$Zm1B= sqrt(pm$no.obs.Z-pm$fisher.dof.cor)*   fisher1(abs(pcor(c(gp,gk,c()),covx)))
   z$Zm1BgivenA= sqrt(pm$no.obs.Z-pm$fisher.dof.pcor1)*   fisher1(abs(pcor(c(gp,gk,pa),covx)))
   z$eo = z$Zm1B -z$Zm1BgivenA
   z
}
# example: zeo2(gp,gk,pa,covx,pm)

# more details version
if(exists("zeo2detail") ) rm(zeo2detail);
zeo2detail=function(gp,gk,pa,covx,pm) {
   z=list();
   z$cor = pcor(c(gp,gk,c()),covx)
   z$pcor = pcor(c(gp,gk,pa),covx)
   z$Zm1B= sqrt(pm$no.obs.Z-pm$fisher.dof.cor)*   fisher1(abs(pcor(c(gp,gk,c()),covx)))
   z$Zm1BgivenA= sqrt(pm$no.obs.Z-pm$fisher.dof.pcor1)*   fisher1(abs(pcor(c(gp,gk,pa),covx)))
   z$BLV = z$Zm1B -z$Zm1BgivenA
   
   z
}

# Handles multiple markers as Grandparents and averages version
if(exists("multimarker.zeo2detail") ) rm(multimarker.zeo2detail);
multimarker.zeo2detail=function(gp,gk,pa,covx,pm) {
   z=list();

   cn = colnames(covx)
   z$multimarker.zeo2detail.grandparents = cn[gp]
   z$multimarker.zeo2detail.parent = cn[pa]
   z$multimarker.zeo2detail.kid = cn[gk]

   z$v.cor=rep(NA,length(gp))
   z$v.pcor=rep(NA,length(gp))
   z$v.Zm1B=rep(NA,length(gp))
   z$v.Zm1BgivenA=rep(NA,length(gp))
   z$v.BLV=rep(NA,length(gp))

   k=1
   for (grandpa in gp) {
     z$v.cor[k] = pcor(c(grandpa,gk,c()),covx)
     z$v.pcor[k] = pcor(c(grandpa,gk,pa),covx)
     z$v.Zm1B[k] = sqrt(pm$no.obs.Z-pm$fisher.dof.cor)*   fisher1(abs(pcor(c(grandpa,gk,c()),covx)))
     z$v.Zm1BgivenA[k] = sqrt(pm$no.obs.Z-pm$fisher.dof.pcor1)*   fisher1(abs(pcor(c(grandpa,gk,pa),covx)))
     z$v.BLV[k] = z$v.Zm1B[k] -z$v.Zm1BgivenA[k]

     k=k+1
   }

   z$mean.cor  = mean(z$v.cor)
   z$mean.pcor = mean(z$v.pcor)
   z$mean.Zm1B = mean(z$v.Zm1B)
   z$mean.Zm1BgivenA = mean(z$v.Zm1BgivenA)
   z$mean.BLV = mean(z$v.BLV)

   z$BLV = z$mean.BLV # for back compatibility with single marker utilizing code; use the mean.

   z
}



# convenience Reporting versions of the four Z-score components, to ease the logging.
# ; return the correlation and the Z-score of the abs value of the correlation (or partial correlation)

if(exists("display.Zm1B") ) rm(display.Zm1B);
display.Zm1B=function(m1,B,covx,pm) { rho=pcor(c(m1,B,c()),covx); Z=sqrt(pm$no.obs.Z-pm$fisher.dof.cor)*   fisher1(abs(rho)); paste(sep="","(rho=",signif(rho,2),";Zabs=",signif(Z,2),")") }

if(exists("display.Zm1BgivenA") ) rm(display.Zm1BgivenA);
display.Zm1BgivenA=function(m1,B,A,covx,pm) { rho=pcor(c(m1,B,A),covx); Z=sqrt(pm$no.obs.Z-pm$fisher.dof.pcor1)*   fisher1(abs(rho)); paste(sep="","(rho=",signif(rho,2),";Zabs=",signif(Z,2),")") }

if(exists("display.Zm2A") ) rm(display.Zm2A);
display.Zm2A=function(m2,A,covx,pm) { rho=pcor(c(m2,A,c()),covx); Z=sqrt(pm$no.obs.Z-pm$fisher.dof.cor)*   fisher1(abs(rho)); paste(sep="","(rho=",signif(rho,2),";Zabs=",signif(Z,2),")") }

if(exists("display.Zm2AgivenB") ) rm(display.Zm2AgivenB);
display.Zm2AgivenB=function(m2,A,B,covx,pm) { rho=pcor(c(m2,A,B),covx); Z=sqrt(pm$no.obs.Z-pm$fisher.dof.pcor1)*   fisher1(abs(rho)); paste(sep="","(rho=",signif(rho,2),";Zabs=",signif(Z,2),")") }

# sample calls:
# display.Zm1B(m1,B,covx,pm)
# display.Zm1BgivenA(m1,B,A,covx,pm)
# display.Zm2A(m2,A,covx,pm)
# display.Zm2AgivenB(m2,A,B,covx,pm)



########################################
# Max-max two-step strategy is different:
# we compute only one total ZEO. From two
# halfs, each of which is maximized independently.
########################################

if(exists("walk.two.steps.max.max") ) rm(walk.two.steps.max.max);
walk.two.steps.max.max=function(pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,n.twofer.eo,cn,twofer.log) {

pm.print(pm,paste("Starting walk.two.steps.max.max()...",date()));

for (i in 1:(length(non.snp.cols)-1)) {
   for (j in (i+1):length(non.snp.cols)) {
      A= non.snp.cols[i];
      B= non.snp.cols[j];

      if (af[A,B]==0 && af[B,A]==0) next; 
      # INVAR: we have an edge from A to B possible.
      M.A = setdiff(which(af[,A]==1),non.snp.cols); # the Markers into A
      M.B = setdiff(which(af[,B]==1),non.snp.cols); # the Markers into B

      # make sure any overlapping belong only to the set with the strongest correlation...
      MAB=union(M.A,M.B);
      M.A=M.B=c();
      for (m in MAB) {
         if (cx[m,A] >= cx[m,B]) { M.A=c(M.A,m); } else { M.B = c(M.B,m); }
      }

      Mset=c(); # grandparent candidates
      parent=NA;
      child=NA;
      singleAB=1; # 1 for M->A->B; -1 for M->B->A
      if (length(M.A) == 0 || length(M.B) == 0)
      {
          if(length(M.A) > 0) {  
              # do a single marker score, one two step, or half a ZEO: M.A -> A -> B
              Mset=M.A;
              parent=A;
              child=B;
              singleAB=1;              
          } else if (length(M.B) > 0) {
              # do a single marker score, one two step, or half a ZEO: M.B -> B -> A
              Mset=M.B
              parent=B;
              child=A;
              singleAB=-1;
          } else next;

          # now compute A-B edge orientation using just the one marker

          # space for results:
          m1ab = cx[Mset,parent]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[Mset];}

          for (k in 1:length(Mset)) {
              m1 = Mset[k];
              z2=zeo2(m1,child,parent,covx,pm); # z2$eo
              m1ab[k] = z2$eo;
          }
          m1ab=sort(m1ab,decreasing=T);
          m1ab.names=names(m1ab);
          top.m1 = which(cn==names(m1ab)[1]); 

          twofer.eo[A,B]=singleAB*(m1ab[1]);
          twofer.eo[B,A]=-twofer.eo[A,B];
          n.twofer.eo[A,B] = n.twofer.eo[B,A] = 2; # vestigial legacy from previous methodology

         # all the details of the components of the computation are here:
         sz1=paste(sep=""," [[Zm1B(",display.Zm1B(top.m1,child,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(top.m1,child,parent,covx,pm),")]] \n");

         twofer.log[A,B]=paste(sep="","\nNo markers into B-only A; for: M1(",names(m1ab)[1],")->A(",cn[parent],")->B(",cn[child],") had Z.2ab=",signif(m1ab[1],3),sz1,"\n")

          next; # b/c rest of code is for markers into both A and B
      }

      if (length(M.A) == 0 || length(M.B) == 0) next; # require markers into both A and B.

      # INVAR: M.A has markers with stronger (>=) connection to A than B
      #    and M.B has markers with stronger connection to B than A

      # compute the Z.2steps
      # reserve space, and get the names() right
      m1ab = cx[M.A,A]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[M.A];}
      m2ba = cx[M.B,B]; m2ba[]=NA; if (length(m2ba)==1) { names(m2ba)=rownames(cx)[M.B];}

      for (k in 1:length(M.A)) {
          m1 = M.A[k];
          z2=zeo2(m1,B,A,covx,pm);
          m1ab[k] = z2$eo; 
      }

      for (k in 1:length(M.B)) {
          m2 = M.B[k];
          z2=zeo2(m2,A,B,covx,pm);
          m2ba[k] = z2$eo; 
      }

      m1ab=sort(m1ab,decreasing=T);
      m2ba=sort(m2ba,decreasing=T);
      m1ab.names=names(m1ab);
      m2ba.names=names(m2ba);

      found.indep.markers=F;
      top.m1 = which(cn==names(m1ab)[1]); 
      top.m2 = which(cn==names(m2ba)[1]);

      if (cx[top.m1,top.m2] < pm$cor.ind.th) {
         #yee-haw, our top ranked markers are indeed independent.
         twofer.eo[A,B]=(m1ab[1]-m2ba[1])/2;
         twofer.eo[B,A]=-twofer.eo[A,B];
         n.twofer.eo[A,B] = n.twofer.eo[B,A] = 4; # vestigial legacy from previous methodology

         # all the details of the components of the computation are here:
         sz1=paste(sep=""," [[Zm1B(",display.Zm1B(top.m1,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(top.m1,B,A,covx,pm),")]] \n");
         sz2=paste(sep=""," [[Zm2A(",display.Zm2A(top.m2,A,covx,pm),")-Zm2A|B(",display.Zm2AgivenB(top.m2,A,B,covx,pm),")]] \n");

         twofer.log[A,B]=paste(sep="","\nM1(",names(m1ab)[1],")->A(",cn[A],")->B(",cn[B],") had Z.2ab=",signif(m1ab[1],3),sz1," M2(",names(m2ba)[1],")->B(",cn[B],")->A(",cn[A],") had Z.2ba=",signif(m2ba[1],3),sz2,"\n*** WITH TOP ranked M1 (out of ",length(M.A),"), WITH TOP ranked M2 (out of ",length(M.B),"). \n*** ZEO[A,B]=", signif(2*twofer.eo[A,B],3),"/2 = ",signif(twofer.eo[A,B],3),"\n*** Evidence from M1(",names(m1ab)[1],") and M2(",names(m2ba)[1],") which are correlated at: ", signif(cx[names(m1ab)[1],names(m2ba)[1]],3)," < pm$cor.ind.th(",signif(pm$cor.ind.th,3),")\n");

      } else {
         # top ranking markers are not independent...so scan down both lists...

         # =================================
         # pass A: M.A has priority, scan down M.B
         # =================================
         Aprior=NA
         Aprior.m1=NA;
         Aprior.m2=NA;
         Aprior.m1ab = Aprior.m2ba = NA;
         done=F;
         for (k1 in 1:length(m1ab)) {
             if (done) break;
             m1=m1ab.names[k1];
             for (k2 in 1:length(m2ba)) {
                 m2=m2ba.names[k2];
                 if (cx[m1,m2] < pm$cor.ind.th) {
                     Aprior=m1ab[k1]-m2ba[k2];
                     Aprior.m1=m1; Aprior.m2=m2;
                     Aprior.k1=k1; Aprior.k2=k2;
                     Aprior.m1ab = m1ab[k1]; Aprior.m2ba = m2ba[k2];
                     done=T;
                     break;
                 }
             }
         }
         # ===============================
         # pass B: M.B has priority, scan down M.A
         # ===============================
         Bprior=NA
         Bprior.m1=NA;
         Bprior.m2=NA;
         Bprior.m1ab = Bprior.m2ba = NA;
         done=F;
         for (k2 in 1:length(m2ba)) {
             if (done) break;
             m2=m2ba.names[k2];
             for (k1 in 1:length(m1ab)) {
                 m1=m1ab.names[k1];
                 if (cx[m1,m2] < pm$cor.ind.th) {
                     Bprior=m1ab[k1]-m2ba[k2];
                     Bprior.m1=m1; Bprior.m2=m2;
                     Bprior.k1=k1; Bprior.k2=k2;
                     Bprior.m1ab = m1ab[k1]; Bprior.m2ba = m2ba[k2];
                     done=T;
                     break;
                 }
             }
         }
         # ===============================
         # handle the 3 possible cases, with full logging so we can understand the results.
         # Either 1) Not so good: we never found independent M1 and M2, resort to single marker analysis.
         #        2) Good: found M1,M2 when M.A had priority and M.B had priority and they agree 
         #           on signs (or there was only one pair of markers). Report success at orienting the edge.
         #        3) Bad: there were two pairs of markers, differing when M.A had priority and when
         #           M.B had priority, and they differ in edge orienting signs. Report evidence is inconclusive.
         # ===============================

         if (is.na(Aprior) || is.na(Bprior)) {
             # no pair of markers in M.A and M.B was independent enough to meet our criteria
             twofer.log[A,B]=paste(sep="","no indep markers found at pm$cor.ind.th=",signif(pm$cor.ind.th,3));
         } else if (Aprior*Bprior > 0) {
             # we agreement in terms of sign-yes! very strong evidence
             # show the stronger of the two pieces of evidence: i.e. prefer the markers m1 and m2 they give strongest edge orienting support.
             if (abs(Aprior) > abs(Bprior)) { # use A
                 twofer.eo[A,B]=Aprior/2;
                 twofer.eo[B,A]=-twofer.eo[A,B];
                 n.twofer.eo[A,B] = n.twofer.eo[B,A] = 4; # vestigial legacy from previous methodology

                 Aprior.m1.num=which(cn==Aprior.m1);
           Aprior.m2.num=which(cn==Aprior.m2);

                 # display the component Z scores and correlations.
                 sz1=paste(sep=""," [[Zm1B(",display.Zm1B(which(cn==Aprior.m1),B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(which(cn==Aprior.m1),B,A,covx,pm),")]] \n");
                 sz2=paste(sep=""," [[Zm2A(",display.Zm2A(which(cn==Aprior.m2),A,covx,pm),")-Zm2A|B(",display.Zm2AgivenB(which(cn==Aprior.m2),A,B,covx,pm),")]] \n");

                 twofer.log[A,B]= paste(sep="","\n**[",cn[A],":",cn[B],"]**","\nM1(",Aprior.m1,")->A(",cn[A],")->B(",cn[B],") had Z.2ab=",signif(Aprior.m1ab,3),sz1,"\n", "M2(",Aprior.m2,")->B(",cn[B],")->A(",cn[A],") had Z.2ba=",signif(Aprior.m2ba,3),sz2,"\n     WITH M.A ranked ",Aprior.k1," (out of ",length(m1ab),"), AND M.B ranked ",Aprior.k2," (out of ",length(m2ba),"). \n*** ZEO[A,B]=", signif(2*twofer.eo[A,B],3), " /2 = ",signif(twofer.eo[A,B],3),"\n*** Evidence from M1(",Aprior.m1,") and M2(",Aprior.m2,") which are correlated at: ", signif(cx[Aprior.m1.num,Aprior.m2.num],3)," < pm$cor.ind.th(",signif(pm$cor.ind.th,3),")\n");
                 
             } else { # use B
                 twofer.eo[A,B]=Bprior/2;
                 twofer.eo[B,A]=-twofer.eo[A,B];
                 n.twofer.eo[A,B] = n.twofer.eo[B,A] = 4; # vestigial legacy from previous methodology

                 # display the component Z scores and correlations.
                 Bprior.m1.num=which(cn==Bprior.m1);
                 Bprior.m2.num=which(cn==Bprior.m2)
                 sz1=paste(sep=""," [[Zm1B(",display.Zm1B(Bprior.m1.num,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(Bprior.m1.num, B,A,covx,pm),")]] \n");
                 sz2=paste(sep=""," [[Zm2A(",display.Zm2A(Bprior.m2.num,A,covx,pm),")-Zm2A|B(",display.Zm2AgivenB(Bprior.m2.num, A,B,covx,pm),")]] \n");

                 twofer.log[A,B]= paste(sep="","\n**[",cn[A],":",cn[B],"]**","M1(",Bprior.m1,")->A(",cn[A],")->B(",cn[B],") had Z.2ab=",signif(Bprior.m1ab,3),sz1,"\n", "M2(",Bprior.m2,")->B(",cn[B],")->A(",cn[A],") had Z.2ba=",signif(Bprior.m2ba,3),sz2,"\n     WITH M.A ranked ",Bprior.k1," (out of ",length(m1ab),"), AND M.B ranked ",Bprior.k2," (out of ",length(m2ba),").\n*** ZEO[A,B]=", signif(2*twofer.eo[A,B],3)," /2 = ",signif(twofer.eo[A,B],3),"\n*** Evidence from M1(",Bprior.m1,") and M2(",Bprior.m2,") which are correlated at: ", signif(cx[Bprior.m1,Bprior.m2],3)," < pm$cor.ind.th(",signif(pm$cor.ind.th,3),")\n");

             }
         } else {
             # disagreement in terms of signs.
             twofer.eo[A,B]=twofer.eo[B,A]=NA;
             n.twofer.eo[A,B]=n.twofer.eo[B,A]=NA;

                 Aprior.m1.num=which(cn==Aprior.m1);
           Aprior.m2.num=which(cn==Aprior.m2);
           Bprior.m1.num=which(cn==Bprior.m1);
           Bprior.m2.num=which(cn==Bprior.m2);

             # A: display the component Z scores and correlations.
             A.sz1=paste(sep=""," [[Zm1B(",display.Zm1B(Aprior.m1.num,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(Aprior.m1.num,B,A,covx,pm),")]] \n");
             A.sz2=paste(sep=""," [[Zm2A(",display.Zm2A(Aprior.m2.num,A,covx,pm),")-Zm2A|B(",display.Zm2AgivenB(Aprior.m2.num,A,B,covx,pm),")]] \n");

             # B: display the component Z scores and correlations.
             B.sz1=paste(sep=""," [[Zm1B(",display.Zm1B(Bprior.m1.num,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(Bprior.m1.num,B,A,covx,pm),")]] \n");
             B.sz2=paste(sep=""," [[Zm2A(",display.Zm2A(Bprior.m2.num,A,covx,pm),")-Zm2A|B(",display.Zm2AgivenB(Bprior.m2.num,A,B,covx,pm),")]] \n");

             part1= paste(sep="","\nXXX Disagreement in terms of signs of M.A priority vs. M.B priority. \nXXX M1(",Aprior.m1,")->A(",cn[A],")->B(",cn[B],") had Z.2ab=",signif(Aprior.m1ab,3),A.sz1,"\nXXX M2(",Aprior.m2,")->B(",cn[B],")->A(",cn[A],") had Z.2ba=",signif(Aprior.m2ba,3),A.sz2,"\nXXX     WITH M.A ranked ",Aprior.k1," (out of ",length(m1ab),"), AND M.B ranked ",Aprior.k2," (out of ",length(m2ba),").\nXXX ZEO[A,B]=", signif(2*twofer.eo[A,B],3), " /2 = ",signif(twofer.eo[A,B],3),"\nXXX Evidence from M1(",Aprior.m1,") and M2(",Aprior.m2,") with cor: ", signif(cx[Aprior.m1.num,Aprior.m2.num],3)," < pm$cor.ind.th(",signif(pm$cor.ind.th,3),")\nXXX");

             part2=paste(sep="","\nXXX M1(",Bprior.m1,")->A(",cn[A],")->B(",cn[B],") had Z.2ab=",signif(Bprior.m1ab,3),B.sz1,"\nXXX M2(",Bprior.m2,")->B(",cn[B],")->A(",cn[A],") had Z.2ba=",signif(Bprior.m2ba,3),B.sz2,"\nXXX     WITH M.A ranked ",Bprior.k1," (out of ",length(m1ab),"), AND M.B ranked ",Bprior.k2," (out of ",length(m2ba),").\nXXX ZEO[A,B]=", signif(2*twofer.eo[A,B],3), " /2 = ",signif(twofer.eo[A,B],3),"\nXXX Evidence from M1(",Bprior.m1,") and M2(",Bprior.m2,") with cor: ", signif(cx[Bprior.m1,Bprior.m2],3)," < pm$cor.ind.th(",signif(pm$cor.ind.th,3),")\n");

             twofer.log[A,B]=paste(sep="","\n**[",cn[A],":",cn[B],"]**",part1,part2);
           }
        } # end top ranking markers are not independent.
   } # end i or A
} # end j or B

list(twofer.eo,n.twofer.eo,twofer.log)
}
# end walk.two.steps.max.max
########################################






########################################
# Multiple Linear Regression two-step strategy builds a multiple regression model for A (and
#  separately B), using the markers in the model more strongly associated with A
#  than with B, by forward stepwise regression. Then we compare how well this model
#  does on B and on A|B.  Likewise for B. Combine these two R^2.
#  
########################################

if(exists("walk.two.steps.weighted.mean.mlreg") ) rm(walk.two.steps.weighted.mean.mlreg);
walk.two.steps.weighted.mean.mlreg=function(x,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,n.twofer.eo,cn,twofer.log) {

pm.print(pm,paste("Starting walk.two.steps.weighted.mean.mlreg()...",date()));

for (i in 1:(length(non.snp.cols)-1)) {
   for (j in (i+1):length(non.snp.cols)) {
      A= non.snp.cols[i];
      B= non.snp.cols[j];

      if (af[A,B]==0 && af[B,A]==0) next; 
      # INVAR: we have an edge from A to B possible.

      # Make sure any overlapping belong only to the set with the strongest correlation...M.A or M.B.
      # And make sure any SNP passes our minimum dependence threshold.
      MAB= snpcols;
      M.A=M.B=c();
      top.A=top.B=NA;
      top.A.cx=top.B.cx=-2;
      for (m in MAB) {
         if (cx[m,A] >= cx[m,B]) { 
              if (cx[m,A] > pm$cor.dep.th) { 
                  # simultaneously semi-sort so strongest correlation is first...
                  if (cx[m,A] > top.A.cx) { top.A = m; top.A.cx=cx[m,A]; M.A=c(m,M.A); } else { M.A=c(M.A,m); } 
              }
         } else {
              if (cx[m,B] > pm$cor.dep.th) {
                  # simultaneously semi-sort so strongest correlated marker is frist
                  if (cx[m,B] > top.B.cx) { top.B = m; top.B.cx=cx[m,B]; M.B = c(m,M.B);} else { M.B = c(M.B,m); }
              }
         }
      }

     if (length(M.A) == 0 && length(M.B) == 0) next;

      # handle the one side only cases...
      Mset=c(); # grandparent candidates
      parent=NA;
      child=NA;
      singleAB=1; # 1 for M->A->B; -1 for M->B->A
      if (length(M.A) == 0 || length(M.B) == 0)
      {
          if(length(M.A) > 0) {  
              # do a single marker score, one two step, or half a ZEO: M.A -> A -> B
              Mset=M.A;
              parent=A;
              child=B;
              singleAB=1;
              step.1 = forward.step.multimarker.twostep(x,cn[A],cn[M.A],cn[M.A][1]);
          } else if (length(M.B) > 0) {
              # do a single marker score, one two step, or half a ZEO: M.B -> B -> A
              Mset=M.B
              parent=B;
              child=A;
              singleAB=-1;
              step.1 = forward.step.multimarker.twostep(x,cn[B],cn[M.B],cn[M.B][1]);
          } else next;

          # now compute A-B edge orientation using just the one marker
          Mset.cn=cn[Mset];

#         #A.coef = step.1$coefficients[-1];
          A.coef = step.1$percent.of.explained.variance.explained;
          A.who.cn  = names(step.1$coefficients)[-1]; # names
          A.who.index  = match(A.who.cn, Mset.cn); # get column index within Mset.cn
          A.who.index.cn = Mset[A.who.index]; # get column index within x
          A.coef.normalized = A.coef / sum(A.coef)

          # reassigning M.A, so we can re-use the following code, inserting the weights....
          old.M.A = M.A
          M.A = A.who.index.cn
          M.A.cn = cn[M.A]
          Mset = A.who.index.cn; # the new, reduced size Mset, from the forward.step selection process
          Mset.cn = cn[A.who.index.cn]

          # space for results:
          m1ab = cx[Mset,parent]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[Mset];}
          m1ab.str=m1ab;

          for (k in 1:length(Mset)) {
              m1 = Mset[k];
              z2=zeo2(m1,child,parent,covx,pm); # z2$eo
              m1ab[k] = z2$eo;

              m1ab.str[k]=paste(sep="","beta_",k," (",signif(A.coef.normalized[k],2),")*[[m1=(",cn[m1],");Zm1B(",display.Zm1B(m1,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(m1,B,A,covx,pm),")]] \n");
          }

          # now do the weighting
          z.mean.a.b = sum(m1ab * A.coef.normalized);
          z.total = z.mean.a.b;

          twofer.eo[A,B]=(singleAB*z.total)/sqrt(2);
          twofer.eo[B,A]=-twofer.eo[A,B];
          n.twofer.eo[A,B] = n.twofer.eo[B,A] = 2; # vestigial legacy from previous methodology

          sz1=paste(m1ab.str,collapse="");
          M.A.cn.together=paste(M.A.cn,collapse=",");
 
          twofer.log[A,B]=paste(sep="","\n**[",cn[A],":",cn[B],"]**","\n No markers into B-only A; for: M1(",M.A.cn.together,")->A(",cn[parent],")->B(",cn[child],") had weighted mean score:",signif(sqrt(2)*twofer.eo[A,B],3),"| components: Z.2ab=",signif(z.mean.a.b,3),"\n",sz1, "\n\n*** twostep.eo[A,B]=", signif(sqrt(2)*twofer.eo[A,B],3)," /sqrt(2) = ",signif(twofer.eo[A,B],3),"\n");

          next; # b/c rest of code is for markers into both A and B
      }


      if (length(M.A) == 0 || length(M.B) == 0) next; # require markers into both A and B.

      # now...both markers...
      M.A.cn = cn[M.A]
      M.B.cn = cn[M.B]

      # do the forward stepwise SNP selection for traits A and B
      step.A = forward.step.multimarker.twostep(x,cn[A],M.A.cn,M.A.cn[1]);
      step.B = forward.step.multimarker.twostep(x,cn[B],M.B.cn,M.B.cn[1]);

      # To weight each Z-score according to the percentage of variance explained
      A.percent.explained.var = step.A$percent.of.explained.variance.explained;

      A.who.cn=names(A.percent.explained.var); #old: rownames(anova.step.A)[1:(len.sum.sq.step.A-1)]
      A.who.index  = match(A.who.cn, M.A.cn); # get column index within M.A
      A.who.index.cn = M.A[A.who.index]; # get column index within x

      B.percent.explained.var = step.B$percent.of.explained.variance.explained;

      B.who.cn=names(B.percent.explained.var); # old: rownames(anova.step.B)[1:(len.sum.sq.step.B-1)]
      B.who.index  = match(B.who.cn, M.B.cn); # get column index within M.A
      B.who.index.cn = M.B[B.who.index]; # get column index within x

      # reassing M.A and M.B, so we can re-use the following code, inserting the weights....
  old.M.A = M.A
      old.M.B = M.B
      M.A = A.who.index.cn
      M.B = B.who.index.cn
      M.A.cn = cn[M.A]
      M.B.cn = cn[M.B]
      

      # INVAR: M.A has markers with stronger (>=) connection to A than B
      #    and M.B has markers with stronger connection to B than A

      # compute the Z.2steps
      # reserve space, and get the names() right
      m1ab = cx[M.A,A]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[M.A];}
      m2ba = cx[M.B,B]; m2ba[]=NA; if (length(m2ba)==1) { names(m2ba)=rownames(cx)[M.B];}

      # store the string output too
      m1ab.str=m1ab;
      m2ba.str=m2ba;

      for (k in 1:length(M.A)) {
          m1 = M.A[k];
          z2=zeo2(m1,B,A,covx,pm); # z2$eo
          m1ab[k] = z2$eo;

          m1ab.str[k]=paste(sep="","percent.explained.var_",k," (",signif(A.percent.explained.var[k],2),")*[[m1=(",cn[m1],");Zm1B(",display.Zm1B(m1,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(m1,B,A,covx,pm),")]] \n");
          
      }

      for (k in 1:length(M.B)) {
          m2 = M.B[k];
          z2=zeo2(m2,A,B,covx,pm);
          m2ba[k] = z2$eo;

          m2ba.str[k]=paste(sep=""," percent.explained.var_",k,"(",signif(B.percent.explained.var[k],2),")*[[m2=(",cn[m2],");Zm2A(",display.Zm2A(m2,A,covx,pm),")-Zm2A|B(",display.Zm2AgivenB(m2,A,B,covx,pm),")]] \n");
      }

      # now do the weighting
      z.mean.a.b = sum(m1ab * A.percent.explained.var);
      z.mean.b.a = sum(m2ba * B.percent.explained.var);
      z.total = z.mean.a.b - z.mean.b.a;

         twofer.eo[A,B]=z.total/2;
         twofer.eo[B,A]=-twofer.eo[A,B];
         n.twofer.eo[A,B] = n.twofer.eo[B,A] = 4; # vestigial legacy from previous methodology

         # all the details of the components of the computation are here:
         sz1=paste(m1ab.str,collapse="");
         sz2=paste(m2ba.str,collapse="");
         details = paste(sz1,"\nMINUS\n",sz2);

         M.A.cn.together=paste(M.A.cn,collapse=",");
         M.B.cn.together=paste(M.B.cn,collapse=",");

         twofer.log[A,B]=paste(sep="","\n**[",cn[A],":",cn[B],"]**","\nM1(",M.A.cn.together,")->A(",cn[A],")->B(",cn[B],") had weighted mean score: ",signif(2*twofer.eo[A,B],3),"|Components: Z.2ab=",signif(z.mean.a.b,3),"\n",sz1,"\n\n*MINUS*\n M2(",M.B.cn.together,")->B(",cn[B],")->A(",cn[A],") had weighted mean Z.2ba=",signif(z.mean.b.a,3),"\n",sz2,"\n\n*** ZEO[A,B]=", signif(2*twofer.eo[A,B],3)," /2 = ",signif(twofer.eo[A,B]/2,3),"\n");
 
   } # end i or A
} # end j or B

list(twofer.eo,n.twofer.eo,twofer.log)
}
# end walk.two.steps.weighted.mean.mlreg
########################################



########################################
#  One permutation-one vote two-step strategy 
#  lets all marker combinations vote, with equal weighting
#  
########################################
if(exists("walk.two.steps.one.permutation.one.vote") ) rm(walk.two.steps.one.permutation.one.vote);
walk.two.steps.one.permutation.one.vote=function(x,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,n.twofer.eo,cn,twofer.log) {

pm.print(pm,paste("Starting walk.two.steps.one.permutation.one.vote()...",date()));

for (i in 1:(length(non.snp.cols)-1)) {
   for (j in (i+1):length(non.snp.cols)) {
      A= non.snp.cols[i];
      B= non.snp.cols[j];
      if (af[A,B]==0 && af[B,A]==0) next; 
      # INVAR: we have an edge from A to B possible.
      M.A = setdiff(which(af[,A]==1),non.snp.cols); # the Markers into A
      M.B = setdiff(which(af[,B]==1),non.snp.cols); # the Markers into B

      # make sure any overlapping belong only to the set with the strongest correlation...
      MAB=union(M.A,M.B);
      M.A=M.B=c();
      for (m in MAB) {
         if (cx[m,A] >= cx[m,B]) { M.A=c(M.A,m); } else { M.B = c(M.B,m); }
      }

      Mset=c(); # grandparent candidates
      parent=NA;
      child=NA;
      singleAB=1; # 1 for M->A->B; -1 for M->B->A
      if (length(M.A) == 0 || length(M.B) == 0)
      {
          if(length(M.A) > 0) {  
              # do a single marker score, one two step, or half a ZEO: M.A -> A -> B
              Mset=M.A;
              parent=A;
              child=B;
              singleAB=1;              
          } else if (length(M.B) > 0) {
              # do a single marker score, one two step, or half a ZEO: M.B -> B -> A
              Mset=M.B
              parent=B;
              child=A;
              singleAB=-1;
          } else next;

          # now compute A-B edge orientation using just the one marker

          # space for results:
          m1ab = cx[Mset,parent]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[Mset];}
          m1ab.str=m1ab; m1ab.str[]=NA;

          for (k in 1:length(Mset)) {
              m1 = Mset[k];
              z2=zeo2(m1,child,parent,covx,pm); # z2$eo
              m1ab[k] = z2$eo;

              m1ab.str[k]=paste(sep="","[[m1=(",cn[m1],");Zm1B(",display.Zm1B(m1,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(m1,B,A,covx,pm),")]] \n");
          }

          mean.m1ab = mean(m1ab);

          twofer.eo[A,B]=singleAB*(mean.m1ab)/sqrt(2);
          twofer.eo[B,A]=-twofer.eo[A,B];
          n.twofer.eo[A,B] = n.twofer.eo[B,A] = 2; # vestigial legacy from previous methodology

         # all the details of the components of the computation are here:

          sz1=paste(m1ab.str,collapse="");
         twofer.log[A,B]=paste(sep="","\n**[",cn[A],":",cn[B],"]**","\nNo markers into B-only A; for: M1(",names(m1ab)[1],")->A(",cn[parent],")->B(",cn[child],") had Z.2ab=",signif(sqrt(2)*twofer.eo[A,B],3)," /sqrt(2) = ",signif(twofer.eo[A,B],3),"...sum of:\n",sz1,"\n");

          next; # b/c rest of code is for markers into both A and B
      }


      if (length(M.A) == 0 || length(M.B) == 0) next; # require markers into both A and B.

      # INVAR: M.A has markers with stronger (>=) connection to A than B
      #    and M.B has markers with stronger connection to B than A

      # compute the Z.2steps
      # reserve space, and get the names() right
      m1ab = cx[M.A,A]; m1ab[]=NA; if (length(m1ab)==1) { names(m1ab)=rownames(cx)[M.A];}
      m2ba = cx[M.B,B]; m2ba[]=NA; if (length(m2ba)==1) { names(m2ba)=rownames(cx)[M.B];}

      # store the string output too
      m1ab.str=m1ab;
      m2ba.str=m2ba;

      for (k in 1:length(M.A)) {
          m1 = M.A[k];
          z2=zeo2(m1,B,A,covx,pm); # z2$eo
          m1ab[k] = z2$eo;

          m1ab.str[k]=paste(sep="","[[m1=(",cn[m1],");Zm1B(",display.Zm1B(m1,B,covx,pm),")-Zm1B|A(",display.Zm1BgivenA(m1,B,A,covx,pm),")]=  ",signif(z2$eo,2),"] \n");
      }

      for (k in 1:length(M.B)) {
          m2 = M.B[k];

          z2=zeo2(m2,A,B,covx,pm); # z2$eo
          m2ba[k] = z2$eo;

          m2ba.str[k]=paste(sep="","[[m2=(",cn[m2],");Zm2A(",display.Zm2A(m2,A,covx,pm),")-Zm2A|B(",display.Zm2AgivenB(m2,A,B,covx,pm),")]=  ",signif(z2$eo,2),"] \n");
      }

      # now do the weighting
      z.mean.a.b = mean(m1ab);
      z.mean.b.a = mean(m2ba);
      z.total = z.mean.a.b - z.mean.b.a;

         twofer.eo[A,B]=z.total/2;
         twofer.eo[B,A]=-twofer.eo[A,B];
         n.twofer.eo[A,B] = n.twofer.eo[B,A] = 4; # vestigial legacy from previous methodology

         # all the details of the components of the computation are here:
         sz1=paste(m1ab.str,collapse="");
         sz2=paste(m2ba.str,collapse="");
         details = paste(sz1,"\nMINUS\n",sz2);

         M.A.cn.together=paste(cn[M.A],collapse=",");
         M.B.cn.together=paste(cn[M.B],collapse=",");

         twofer.log[A,B]=paste(sep="","\n**[",cn[A],":",cn[B],"]**","\nOverall mean before division=",signif(z.total,3)," from Z.mean.a->b(",signif(z.mean.a.b,3),") MINUS Z.mean.ba(",signif(z.mean.b.a,3),") = ",signif(z.total,3),"\nM1(",M.A.cn.together,")->A(",cn[A],")->B(",cn[B],") had weighted mean score: ",signif(2*twofer.eo[A,B],3),"|Components: Z.2ab=",signif(z.mean.a.b,3),"\n",sz1,"\n\n*MINUS*\n M2(",M.B.cn.together,")->B(",cn[B],")->A(",cn[A],") had weighted mean Z.2ba=",signif(z.mean.b.a,3),"\n",sz2,"\n\n*** ZEO[A,B]=", signif(2*twofer.eo[A,B],3)," /2 = ",signif(twofer.eo[A,B],3),"\n");

   } # end i or A
} # end j or B

list(twofer.eo,n.twofer.eo,twofer.log)
}
# end walk.two.steps.one.permutation.one.vote
########################################




if(exists("viewnodes") ) rm(viewnodes);
viewnodes=function(eigendag,from.nodes=c(),to.nodes=c()) {
   eigendag$stats[grep(paste(sep="",from.nodes,".* -> .*",to.nodes,".*"),as.character(rownames(eigendag$stats))),] 
}


# make.ram: turn a dag matrix into a ram specification for MaxLik model fitting.
if(exists("make.ram") ) rm(make.ram);
make.ram=function (m,list.of.hidden.confounded.pair.indices=NULL)
{
    cn=colnames(m);
    nc=length(cn);
    edgelist=list();
    my.par=list();
    my.start=list();

    hidden.count=0;
    if (!is.null(list.of.hidden.confounded.pair.indices)) {
       hidden.count=length(list.of.hidden.confounded.pair.indices);
    }

    # track which variances we have put in the model already
    var.already.in=vector(length=nc);
    var.already.in[]=F;

    used.cols=c(); # track the nodes that are used, so we can omit others from the covar matrix and avoid warnings
    for (j in 1:nc) { for (i in 1:nc) {
        if (i==j && any(m[,i]!=0, m[i,]!=0) && !var.already.in[i] ) {
                       edgelist=c(edgelist,paste(cn[i],"<->",cn[i])); 
                       my.par=c(my.par,paste(sep="","var.",cn[i])); 
                       my.start=c(my.start,NA);
                       var.already.in[i]=TRUE;
                  } else { 
                       if(m[i,j]==1) { edgelist=c(edgelist,paste(cn[i],"->",cn[j])); 
                                       my.par=c(my.par,paste(sep="","fr.",cn[i],".",cn[j])); 
                                       my.start=c(my.start,NA);
                                       used.cols=c(i,j,used.cols);
                                       }}
    }}

     # manually add in the covariance between A.col and B.col that is due to hidden var
     if (!is.null(list.of.hidden.confounded.pair.indices)) {
        pair.count=1;
        for (pair in list.of.hidden.confounded.pair.indices) {
             # old way, which was same as M.AhiddenB (model 4 and 6 were not identifiable)
              # ...is the best way...!
             edgelist=c(edgelist,paste(cn[pair[1]],"<->",cn[pair[2]]));
             my.par=c(my.par,paste(sep="","covar.",cn[pair[1]],".",cn[pair[2]]));
             my.start=c(my.start,NA);
             used.cols=c(pair[1],pair[2],used.cols);

#              # new way,  works, but same as above
#              latent = paste(sep="","Lat",pair.count);
#       
#              edgelist=c(edgelist,paste(latent,"<->",latent));
#              #my.par=c(my.par,paste(sep="","Var.",latent));
#              my.par=c(my.par,NA);
#              my.start=c(my.start,1);
#
#              edgelist=c(edgelist,paste(sep="",latent," -> ",cn[pair[1]]));
#              edgelist=c(edgelist,paste(sep="",latent," -> ",cn[pair[2]]));
#              my.par=c(my.par,paste(sep="",latent,".",cn[pair[1]]));
#              my.par=c(my.par,paste(sep="",latent,".",cn[pair[2]]));
#              #my.par=c(my.par,NA,NA)
#              my.start=c(my.start,1,1);
#              used.cols=c(pair[1],pair[2],used.cols);

            if (!var.already.in[pair[1]]) {
                       edgelist=c(edgelist,paste(cn[pair[1]],"<->",cn[pair[1]])); 
                       my.par=c(my.par,paste(sep="","var.",cn[pair[1]])); 
                       my.start=c(my.start,NA);
                       var.already.in[pair[1]]=TRUE;
            }

            if (!var.already.in[pair[2]]) {
                       edgelist=c(edgelist,paste(cn[pair[2]],"<->",cn[pair[2]])); 
                       my.par=c(my.par,paste(sep="","var.",cn[pair[2]])); 
                       my.start=c(my.start,NA);
                       var.already.in[pair[2]]=TRUE;
            }
             pair.count=pair.count+1;
        }
     }
     # end manual addition of hidden covar between A and B

    suppressWarnings({the.ram <- cbind(paste(edgelist), paste(my.par), as.numeric(paste(my.start)))});
    # important: NA in my.par indicates a fixed variable
    if (any(is.na(my.par))){ the.ram[which(is.na(my.par)),2]=NA; } # make real NA, not string "NA"
    class(the.ram) <- "mod"
    my.list=list();
    my.list$the.ram=the.ram
    my.list$used.cols=unique(used.cols)
    my.list
}

# pcor(): partial correlation, copied from ggm so we don't  necessarily need the whole ggm library

#if(exists("pcor") ) rm(pcor);
pcor=function (u, S) 
{
    k <- solve(S[u, u])
    -k[1, 2]/sqrt(k[1, 1] * k[2, 2])
}

# fisher's transform to take the correlation coefficient into a normal distribution
if(exists("fisher1") ) rm(fisher1);
fisher1=function(r) {
   0.5*log((1+r)/(1-r))
}

# DFScycle(): use Depth-first-search to
#        check if a cycle is created by adding edge from i->j # i.e. starting from j can we return then to i?
#         This function neede in the original FTC algorithm by Cohen et al.
# PRE: diagonal of dagmat should be 0; no self-loops.
if(exists("DFScycle") ) rm(DFScycle);
DFScycle=function(dagmat, i,j) {
    if (dagmat[j,i]) { return(TRUE); }
    for(k in which(dagmat[j,] !=0)) {
        if (DFScycle(dagmat,i,k)) { return(TRUE); }
    }
    return(FALSE);
}

# check22(): helper function to check 2 quantitative traits (A and B) and 2 SNPs for
#                     A->B or B->A or neither. 
# covx is the covariance matrix of the data matrix x
# cx is abs value of correlation matrix of data matrix x
# A, B, m1 and m2 are column numbers, referring to x. PRE: m1 and m2 specify SNP columns
# pm is a list of parameters
# return value $edge: -1 => A->B, 0 => no edge, 1=> B->A
# return value $zeo = sum of the z-scores for the A->B edge correlation propositions; the positive value => A->B

if(exists("check22") ) rm(check22);
check22=function(covx, cx, A, B, m1, m2, pm, no.obs,zeo.log) {
     ab=list();
     ab$edge=0;
     ab$zeo = 0;
     ab$zeo.log = zeo.log
     ab$why.not

     # avoid double counting: only do A->B if A<B
     if (A > B) { ab$why.not="we only check A,B when the numbering is A<B"; return(ab); }

    # are entry conditions satisfied?
    if (cx[m1,m2] > pm$cor.ind.th) { ab$why.not = "cx[m1,m2] > pm$cor.ind.th"; return(ab);}
    if (cx[m1,A] < pm$cor.dep.th) {  ab$why.not = "cx[m1,A] < pm$cor.dep.th"; return(ab); }
    if (cx[m2,B] < pm$cor.dep.th) {  ab$why.not = "cx[m2,B] < pm$cor.dep.th"; return(ab); }
    if (cx[A,B]  < pm$cor.dep.th) {  ab$why.not = "cx[A,B]  < pm$cor.dep.th"; return(ab); }

    # do we discard b/c looks like Model 4: confounding?
    if (cx[m2,A] < pm$cor.ind.th && cx[m1,B] < pm$cor.ind.th) { ab$why.not=" cx[m2,A] < pm$cor.ind.th && cx[m1,B] < pm$cor.ind.th"; return(ab); }
    
    # do we discard b/c looks like epistasis/both markers interact with both traits?
    if (cx[m2,A] > pm$cor.dep.th && cx[m1,B] > pm$cor.dep.th) { ab$why.not = "cx[m2,A] > pm$cor.dep.th && cx[m1,B] > pm$cor.dep.th "; return(ab)};
    
    # compute ZEO: Z-score Edge Orienting statistic: + meaning A-> B 
    z=zeo2(m1,B,A,covx,pm); # z$eo
    w=zeo2(m2,A,B,covx,pm); # w$eo
    ab$zeo = z$eo - w$eo;

    cn=colnames(cx);
    print(digits=3,paste("Setting ZEO=",signif(ab$zeo,3),"for m1,m2,a,b:",cn[m1],cn[m2],cn[A],cn[B]," of --+ sum of these ",signif(z$Zm1B,3),signif(w$Zm1B,3),signif(z$Zm1BgivenA,3),signif(w$Zm1BgivenA,3)));

    str.Zm1B=paste(sep="",cn[m1],".",cn[B],"(",signif(z$Zm1B,2),")");
    str.Zm1BgivenA= paste(sep="",cn[m1],".",cn[B],"|",cn[A],"(",signif(z$Zm1BgivenA,2),")");
    str.Zm2A=paste(sep="",cn[m2],".",cn[A],"(",signif(w$Zm1B,2),")");
    str.Zm2AgivenB= paste(sep="",cn[m2],".",cn[A],"|",cn[B],"(",signif(w$Zm1BgivenA,2),")");
    str.zeo = paste(sep="","  [", str.Zm1B,"-",str.Zm2A,"-",str.Zm1BgivenA, "+", str.Zm2AgivenB,"=",signif(ab$zeo,2),"]  ");

    # log it
    ab$zeo.log[A,B]=paste(ab$zeo.log[A,B],"+",str.zeo,"\n");
    ab$zeo.log[B,A]=paste(ab$zeo.log[B,A],"-",str.zeo,"\n");

    # if both of eqns 8.6 and 8.8 favor one model then conclude that model
    if (abs(pcor(c(m1,B,A),covx)) < pm$cor.ind.th && abs(pcor(c(m2,A,B),covx)) > pm$cor.dep.th) {ab$edge=-1; return(ab); }
#A->B
    if (abs(pcor(c(m1,B,A),covx)) > pm$cor.dep.th && abs(pcor(c(m2,A,B),covx)) < pm$cor.ind.th) {ab$edge=1; return(ab); } 
#B->A
    print(paste("check22 concluding neither b/c abs(pcor(c(m1,B,A),covx))=", signif(abs(pcor(c(m1,B,A),covx)),3),"and abs(pcor(c(m2,A,B),covx))=", signif(abs(pcor(c(m2,A,B),covx)),3)))
    return(ab); # default
}

# checkLCD(): helper function to check 2 quantitative traits (A and B) and 1 SNPs for
#                     SNP->A->B or SNP->B->A
# covx is the covariance matrix of the data matrix x
# cx is abs value of correlation matrix of data matrix x
# A, B, m1 are column numbers, referring to x. PRE: m1 specifies a SNP column
# pm is a list of parameters
# return value: 1 => A->B, 0 => no edge or B->A

if(exists("checkLCD") ) rm(checkLCD);
checkLCD=function(covx, cx, A, B, m1, pm) {
    r=list(); r$edge = FALSE;

    # avoid double counting: only do A->B if A<B
    if (A > B) return(r);

    # are entry conditions satisfied?
    if (cx[m1,A] >= pm$cor.dep.th) {
        if (cx[A,B]  >= pm$cor.dep.th) {

        # LCD is a simple conditional independence
            apc= abs(pcor(c(m1,B,A),covx));
            if (apc <= pm$cor.ind.th) { r$edge=TRUE; r$drop = cx[m1,B]-apc;}  #A->B
        }
    }
    r; # default 
}











# checkV()
#
#  check for a V-structure between pa1 and pa2 and their child, sink.
#
if(exists("checkV") ) rm(checkV);
checkV=function(covx, cx, sink, pa1, pa2, pm) {
    r=list(); r$v = FALSE; r$pa.cor= cx[pa1,pa2]; r$sink.pa1 = cx[sink,pa1]; r$sink.pa2 = cx[sink,pa2];
    apc= abs(pcor(c(pa1,pa2,sink),covx));
    r$cond.dep.pa=apc;
    if (cx[pa1,pa2] <= pm$cor.ind.th) { # parents must be uncorrelated
        if (cx[sink,pa1] >= pm$cor.dep.th) {
            if (cx[sink,pa2]  >= pm$cor.dep.th) {
                if (apc >= pm$induced.dep.th) { r$v=TRUE; } # print(paste("checkV says YES with sink",cn[sink],"and parents",cn[pa1],cn[pa2]));}
            }
        }
    }
    r; # default 
}

# checkRVV()
#
#  pa1 can be either marker or cts trait (like gene expression)
# 
#  gpa2.1 and gpa2.2 are are potential grandparents (they are parents of pa2)
#
#  return 1 if pa2 -> sink is verified. else 0.
if(exists("checkRVV") ) rm(checkRVV);
checkRVV=function(covx, cx, sink, pa1, pa2, gpa2.1, gpa2.2, pm, cn) {
      #print(paste("checkRVV checking ",cn[pa2],"->",cn[sink],"using",cn[gpa2.1],cn[gpa2.2]));
      r=list(); r$rvv=FALSE;
  v1=checkV(covx,cx, sink, pa1, pa2, pm);
      if (!v1$v) { return(r); } #print("aborting checkRVV:failed first v-structure"); return(r); }
      
  v2=checkV(covx, cx, pa2, gpa2.1, gpa2.2, pm);
      if (!v2$v) { return(r); } # print("aborting checkRVV:failed second v-structure"); return(r); }

      apc1 = abs(pcor(c(gpa2.1,sink,pa2),covx));
      apc2 = abs(pcor(c(gpa2.2,sink,pa2),covx));

      gp.gc1 = cx[gpa2.1,sink];
      gp.gc2 = cx[gpa2.2,sink];

  #print(paste("RVV diag: cor between",cn[gpa2.1],"and",cn[sink]," given ",cn[pa2],":",apc1));
  #print(paste("RVV diag: cor between",cn[gpa2.2],"and",cn[sink]," given ",cn[pa2],":",apc2));

  #print(paste("RVV diag: grandparent-child cor between",cn[gpa2.1],"and",cn[sink],":",gp.gc1));
  #print(paste("RVV diag: grandparent-child cor between",cn[gpa2.2],"and",cn[sink],":",gp.gc2));

      if (apc1 <= pm$cor.ind.th && gp.gc1 > pm$cor.dep.th) {
          r$rvv=TRUE;
          r$sink = sink;
          r$pa = pa2;
          r$gpa = gpa2.1;
          r$gpa.sink.cor = gp.gc1;
          r$cond.dep.gpa.sink = apc1;
          #print(paste("RVV accepting on ",cn[gpa2.1]));
          return(r);
      }

      if (apc2 <= pm$cor.ind.th && gp.gc2 > pm$cor.dep.th) {
          r$rvv=TRUE;
          r$sink = sink;
          r$pa = pa2;
          r$gpa = gpa2.2;
          r$gpa.sink.cor = gp.gc2;
          r$cond.dep.gpa.sink = apc2;
          #print(paste("RVV accepting on ",cn[gpa2.2]));
          return(r);
      }
   #print("RVV rejecting.");
   r;
}

# make drawing the returned RVV structure a seperate function for ease of repetition
if(exists("graph") ) rm(neo.graph);
neo.graph=function(z, adjustme=F, nowind=F, noparms=F, no.edge.labels=F){
if (!nowind) { cross.platform.windows(); }

res=drawGraph(z$M, coor=z$coor,adjust=adjustme,beta=0);

if (!noparms) {
 mtext(paste("cor.dep.th =",z$pm$cor.dep.th,
            "cor.ind.th =",z$pm$cor.ind.th,
          "  pcor.th =",z$pm$pcor.th, 
          " SNP.below =",z$pm$snp.below),line=2);
          
 mtext(paste("omega.th =",z$pm$omega.th,
            "  k.pcor.th =",z$pm$k.pcor.th,
            "  ind.dep.th =", z$pm$induced.dep.th,"  edge.th =", z$pm$graph.arrow.edge.th),line=0);
}

#sem=summary(z$sem.fit)

#mtext(paste("RMSEA=",sem$RMSEA),line=0);

if (!is.null(z$my.title)) { title(sub=z$my.title); }

# scores are no longer symmetric, so redo this logic:
if (!no.edge.labels) {
  cn=colnames(z$M);
  # res has final.coordinates
  for (ii in 1:(length(z$non.snp.cols)-1)) { 
      for (jj in (ii+1):length(z$non.snp.cols)) {
         i = z$non.snp.cols[ii]
         j = z$non.snp.cols[jj]

      # in the lower triangle, if the entry is na or negative, then check the upper triangle
      if (is.nan(z$final.eo[i,j]) | is.na(z$final.eo[i,j]) | (z$final.eo[i,j] <= 0)) {
            tmp=j; j=i; i=tmp; # swap i and j
      }

      if (!is.nan(z$final.eo[i,j]) && !is.na(z$final.eo[i,j])) { 
            x=(res[i,1]+res[j,1])/2;
            y=(res[i,2]+res[j,2])/2;
            i.tmp=i; j.tmp=j;
            score= z$final.eo[i,j];
            if (score < 0) {i.tmp=j; j.tmp=i; score = z$final.eo[j,i]; } # try the opposite one if neg
            if (score <= z$pm$graph.arrow.edge.th/2) next; # skip if still negative or 0 or small

            # supernode stuff, only if available
            sup=""
            if ("super" %in% names(z)) {
                 if (z$super[i]==z$super[j]) { sup="[S]"; }
            }

            arr="->"; if (abs(score) < z$pm$graph.arrow.edge.th) { arr="-]";}
            my.legend= paste(sep="",sup,cn[i.tmp],arr,cn[j.tmp],"=",signif(score,2));
            legend(x,y,legend=my.legend,bty="n",xjust=0.5,yjust=.25);            
#          } else if (z$final.eo[i,j] < (-pm$graph.arrow.edge.th)) { 
#            
#          }
      }
  }}
}

return(res);
}


# mapstar: function to allow us to print nice [*] vs. [ ] outputs
mapstar=function(n) { if (is.na(n) || (!is.finite(n)) || n==0) { return(' '); } else { return('*'); }}

# helper function to enumerate all pairs of a list of variables into a matrix, so we can
#  check possible parents and possible SNP grandparents
#
# L = a list of 2 or more numbers
if(exists("enumpairs") ) rm(enumpairs);
enumpairs=function(L) {
    if(length(L) < 2) stop("enumpairs called with list of length < 2");
    r=c()
    lenl = length(L);
    if (lenl==2) return(rbind(c(),L));

    for(i in 1:(lenl-1)) {
        for (j in (i+1):lenl) {
            r=rbind(r,c(L[i],L[j]));
        }}
      return(r);
    }

# filter SNPS by hclustering
if(exists("filterSNP") ) rm(filterSNP);
filterSNP=function(snps, cor.below.cluster.def=.3) {
    #h1=hclust(dist(t(snps)),method="average");
    tmp.snp.cor = cor(snps, method="p",use="pairwise.complete.obs")
    h1= hclust(as.dist( 1-abs(t(tmp.snp.cor))), method="average" );
    #plot(h1);
    a=cutree(h1,h=cor.below.cluster.def);
    cna=names(a);
    my.clusters=unique(a);
    my.snp.choices=my.clusters; my.snp.choices[]=0;
    for (i in my.clusters) {
        cluster.i = cna[a %in% i];
        # choose SNP closest to mean SNP response for this cluster.
        data.i = snps[,cluster.i];
        mean.i = mean(as.data.frame(t(data.i)));
        cor.i=cor(mean.i,data.i, method="p",use="pairwise.complete.obs");
        my.snp.choices[i]=which(cluster.i[which.max(cor.i)]==cna);
    }
    my.snp.choices
}

# filter first by hclustering SNPs, then by picking the
# one in each cluster with the strongest non-SNP correlation.
if(exists("filterSNP.strongPC") ) rm(filterSNP.strongPC);
filterSNP.strongPC=function(snps, non.snp.cols, x, cor.below.cluster.def=.3, force.cutree.count=NA) {
    #h1=hclust(dist(t(snps)),method="average");

    # handle NA and infinite values....
    ddd=as.dist( 1-abs(t(cor(snps, method="p",use="pairwise.complete.obs"))));
    ddd[which(is.na(ddd))] = 1;
    ddd[which(is.infinite(ddd))] = 1;

    h1= hclust(ddd, method="average" );
    #plot(h1);
    if (is.na(force.cutree.count)) {
         a=cutree(h1,h=cor.below.cluster.def);
    } else {
         a=cutree(h1,k=min(ncol(snps),force.cutree.count)); # force 20 groups for the 3000 SNP dataset
    }

    cna=names(a);
    non.cna=colnames(x)[non.snp.cols];
    my.clusters=unique(a);
    my.snp.choices=my.clusters; my.snp.choices[]=0;
    my.snp.abscor=my.snp.choices; # store abs cor level for PC chosen.
    my.snp.pc =c(); # store pc chosen.

    num.non = length(non.snp.cols)
    for (i in my.clusters) {
        cluster.i = cna[a %in% i];
        cluster.i.size = length(cluster.i);

        # choose SNP with maximum response in the non-snps
        data.i = snps[,cluster.i];
        data.i.x = cbind(data.i, x[,non.snp.cols]);
        cor.data.i.x = abs(cor(data.i.x, method="p",use="pairwise.complete.obs"));

        # zero out the missing and infinite bad boys
        cor.data.i.x[which(is.na(cor.data.i.x))] = 0;
        cor.data.i.x[which(is.infinite(cor.data.i.x))] = 0;

        ncolhere = ncol(cor.data.i.x);

        # zero the SNP:SNP and pc:pc correlations, so they don't compete
        cor.data.i.x[1:cluster.i.size,1:cluster.i.size]=0;
        non.snp.zeroable = (cluster.i.size+1):(cluster.i.size+num.non)
        cor.data.i.x[non.snp.zeroable, non.snp.zeroable]=0;

        # zero one part of the symmetric matrix, so we get a SNP chosen
        cor.data.i.x[1:cluster.i.size, non.snp.zeroable]=0;

        wm = which.max(cor.data.i.x);
        best.row = ((wm-1)%% ncolhere)+1;
        best.col = 1 + ((wm-1)  %/% ncolhere);

        my.snp.choices[i]=cluster.i[best.col];
        my.snp.abscor[i]=signif(cor.data.i.x[wm],3);
        my.snp.pc = c(my.snp.pc, non.cna[best.row - cluster.i.size]);
        #print(my.snp.pc)
    }
    snp.name = cna[my.snp.choices];
    data.frame(my.snp.choices,my.snp.abscor,my.snp.pc)
}

# filterSNP.maxpriority.PC()
#
# Restrict chosen SNPs to having minor allele frequency above a specified number.
#
# Replace first and second steps below with the following: for each non-snp trait, do
#  foward stepwise regression to select SNPs from the big pool of SNPs, so that we can
#  and will catch the multiple SNPs that may influence each non-SNP, especially those
#  that may offer orthogonal influence.
#
# First, allow each non-snp priority: pick the SNP most strongly associated with each PC/eigengene (in general, a non-SNP).
# Second, condition each non-snp on the chosen snp and search again for next most significant snp. Repeat to get 3-4 most strong SNPs.
# Thirdly, cluster the SNPs, and for each of say 10 clusters, pick the snp from the cluster that is most strongly
# associated with over all the candiate non-SNPs. This should get us a sense of any other strong
# clusters of SNPs that can act as nice somewhat orthogonal controls.

#
#

if(exists("filterSNP.maxpriority.PC") ) rm(filterSNP.maxpriority.PC);
filterSNP.maxpriority.PC=function(snpcols, non.snp.cols, x, cor.below.cluster.def=.3, force.cutree.count=NA, orthog.search.depth=5,view.step=F,top.N.snps.per.trait=12,pm) {
    #h1=hclust(dist(t(snps)),method="average");

    cn=colnames(x)

    # single snp needs special handling b/c we don't get matrices like we expect
    if (length(snpcols) ==1) {
       my.snp.choice=cn[snpcols];
       signed.c = cor(x[snpcols],x[non.snp.cols],method="pearson",use="pairwise.complete.obs") 
       ac= abs(signed.c)
       wmax= which.max(ac);
       my.snp.abscor=signif(ac[wmax],3);
       my.snp.cor=signif(signed.c[wmax],3);
       my.snp.pc=colnames(ac)[wmax];
       my.pass="greedy";
       non.redundant=TRUE
       redundant.note = " "
       return(data.frame(my.snp.choice,my.snp.abscor,my.snp.cor,my.snp.pc,my.pass,non.redundant,redundant.note)) 
    }

    
    x=data.frame(x)
    snps = data.frame(x[,snpcols])
    cn.snp = colnames(snps);

    # sanity
    s=colSds(as.matrix(snps))
    if(any(s==0)) {
       stop(paste("Fatal error: Some SNPs have no variation: ",cn.snp[s==0]))
    }
    # ===========================================================
    # find any SNPs with minor allele frequency below our cutoff
    # ===========================================================
    too.small=c();
    if (!is.null(pm$minimum.minor.allele.freq)) {
       for (i in 1:ncol(snps)) {
           ta= table(snps[,i]); # tables ignores NA
           w.min= which.min(ta)
           ith.min.allele=names(w.min);
           allele.freq = ta[w.min]/sum(ta)
           if (length(ta) < 2) { too.small=c(too.small,i); pm.print(pm,paste(sep="","Eliminating non-varying SNP ",cn.snp[i]));} # ignore non-varying SNPs
           if (allele.freq < pm$minimum.minor.allele.freq) {
                pm.print(pm,paste(sep="","dropping SNP '",cn.snp[i],"' because minor allele ",ith.min.allele," had frequency ",allele.freq," below our cutoff of ",pm$minimum.minor.allele.freq));
                too.small=c(too.small,i);
           }
       }
    }
 
    # reassign snp and cn.snp if we are shrinking our SNP selection
    if (length(too.small) > 0)
    {
        snps=snps[,-too.small];
        cn.snp=colnames(snps);
        snpcols = setdiff(snpcols,too.small)
    }
    
    # =================================================================
    # first pick out the (top.N.snps.per.trait) strongest correlated SNP for each non.snp trait
    # =================================================================

    pm.print(pm,paste("Starting greedy snp selection: ",date()))

    # can't pick more SNPs than we have.
    if (top.N.snps.per.trait > length(snpcols)) { top.N.snps.per.trait=length(snpcols); }

    signed.cor.snp.non = cor(snps,x[,non.snp.cols],method="pearson",use="pairwise.complete.obs")
    cor.snp.non = abs(signed.cor.snp.non)

    rownames(cor.snp.non) = cn.snp;
    cn.non= colnames(x)[non.snp.cols];
    colnames(cor.snp.non) = cn.non;

    rownames(signed.cor.snp.non) = cn.snp
    colnames(signed.cor.snp.non) = cn.non

    # zero out the missing and infinite bad boys
    cor.snp.non[which(is.na(cor.snp.non))] = 0;
    cor.snp.non[which(is.infinite(cor.snp.non))] = 0;

    signed.cor.snp.non[which(is.na(cor.snp.non))] = 0;
    signed.cor.snp.non[which(is.infinite(cor.snp.non))] = 0;

    max.snp = cor.snp.non[1:top.N.snps.per.trait,]; max.snp[]=NA; rownames(max.snp)=NULL;
    # fix the above line when only 1 row...
    if (1 == top.N.snps.per.trait) { 
        max.snp = matrix(nrow=top.N.snps.per.trait,cor.snp.non[1:top.N.snps.per.trait,]);
        max.snp[]=NA; 
  rownames(max.snp)=NULL;
        colnames(max.snp)=colnames(cor.snp.non)
    }

    max.snp.cor = max.snp; # empty with NA but labelled matrix
    max.snp.signed.cor = max.snp; # empty but labelled matrix
    max.names=max.snp; # empty but labelled matrix
    pc.namer=max.snp; # make linearization easy later.


    for (i in 1:ncol(cor.snp.non)) {

         # use the A,B edge filters if present to avoid picking SNPs for irrelevant traits.
         if (!is.null(pm$A) & !is.null(pm$B) & !(non.snp.cols[i] %in% pm$A | non.snp.cols[i] %in% pm$B)) {
            pm.print(pm,paste("Skipping greedy SNP selection for trait",cn.non[i],"since it wasn't on the pm$A or pm$B lists."))
            next;
         }

         # Are we at a forced SNP assignment?
         if (!is.null(pm$forced.MA.colnum) & !is.null(pm$A) & non.snp.cols[i] %in% pm$A) {
            nforced = length(pm$forced.MA.colnum)
            
            # and eliminate NA values--if, for instance, a forced snp had no variation or too small minor allele freq.
            tmp = match(pm$forced.MA.cn, cn.snp)
            if (any(is.na(tmp))) {
                  tmp = tmp[!is.na(tmp)];
                  pm.print(pm,paste("Warning: non-variation or below minor allele freq threshold causing us to eliminate a pm$forced.MA.colnum SNP assignment for ",cn.non[i]));
            }
            if (length(tmp)) {
               max.snp[1:nforced,i]=tmp
               max.names[1:nforced,i]=cn.snp[tmp]
            }

           pm.print(pm,paste("=================== CONCLUSION OF (Forced) SNP SELECTION DETAILS FOR ",cn.non[i]," === at ",date()))
           pm.cat(pm,paste(sep="","snps.for.",cn.non[i],"=c(\"",paste(collapse="\",\"",cn.snp[max.snp[,i]]),"\")\n"))
           pm.print(pm,      "================================================================================")            
            
            next;
         }
         # if in pm$A, it can't be in pm$B. Check if B has forced SNP assignment too.
         if (!is.null(pm$forced.MB.colnum) & !is.null(pm$B) & non.snp.cols[i] %in% pm$B) {
            nforced = length(pm$forced.MB.colnum)

            tmp = match(pm$forced.MB.cn, cn.snp)
            if (any(is.na(tmp))) {
                  tmp = tmp[!is.na(tmp)];
                  pm.print(pm,paste("Warning: non-variation or below minor allele freq threshold causing us to eliminate a pm$forced.MB.colnum SNP assignment for ",cn.non[i]));
            }
            if (length(tmp)) {
               max.snp[1:nforced,i]=tmp
               max.names[1:nforced,i]=cn.snp[tmp]
            }
            
           pm.print(pm,paste("=================== CONCLUSION OF (Forced) SNP SELECTION DETAILS FOR ",cn.non[i]," === at ",date()))
           pm.cat(pm,paste(sep="","snps.for.",cn.non[i],"=c(\"",paste(collapse="\",\"",cn.snp[max.snp[,i]]),"\")\n"))
           pm.print(pm,      "================================================================================")            
            
           next;
         }
         
         # reset eliminated.snp each time through so each variable gets to pick its own best greedy SNPs
         eliminated.snp=cor.snp.non; eliminated.snp[]=FALSE # track if SNP already eliminated

         scsn=sort(cor.snp.non[,i],index.return=TRUE,decreasing=TRUE);

         # now we need to weed out the redundant snps
         
         if (pm$skip.greedy.snps.if.redundant) {

           max.snp[1,i] = scsn$ix[1]
           if (pm$show.snp.selection.log) {
              pm.print(pm,paste("========================================================== SNP SELECTION DETAILS FOR ",cn.non[i]," ==="))
              pm.print(pm,paste("filterSNP.maxpriority.PC(): For trait",cn.non[i],
                             ": Choosing STRONGEST greedy SNP",cn.snp[scsn$ix[1]],
                             "with snp-trait cor =",signif(signed.cor.snp.non[scsn$ix[1],i],3)))
           }
           next.non.redundant = 2
           if (top.N.snps.per.trait > 1) {

           for(ki in 2:nrow(cor.snp.non)) {
               current.snp.index = scsn$ix[ki]
               if (pm$show.snp.selection.log) {
                  pm.print(pm,paste("considering: ",cn.snp[current.snp.index])) #debug
               }
               flag.ki.redun=FALSE
               closest.cor = -2
               closest.cor.index = NA

               for(kj in 1:(ki-1)) {
                 kj.idx = scsn$ix[kj]

                 # exlude checking already eliminated SNPS...
                 if (eliminated.snp[kj.idx,i]) next;

                 cor.cur.ki.with.kj = cor(x[,current.snp.index],x[,kj.idx],method="pearson",use="pairwise.complete.obs")

                 if (abs(cor.cur.ki.with.kj) >= pm$snp.redundancy.cor.level) {
                       # don't want this ki
                       if (pm$show.snp.selection.log) {
                           pm.print(pm,paste("filterSNP.maxpriority.PC(): Skipping greedy SNP",cn.snp[current.snp.index],
                                  "because it has abs cor =",signif(cor.cur.ki.with.kj,3),">= ",pm$snp.redundancy.cor.level,
                                  "(pm$snp.redundancy.cor.level), with previously chosen SNP",cn.snp[kj.idx]))
                       }         
                       eliminated.snp[current.snp.index,i]=TRUE
                       flag.ki.redun=TRUE
                       break;
                 }

                 # track how close we got to the redundancy threshold
                 if (abs(cor.cur.ki.with.kj) > closest.cor) {
                     closest.cor = abs(cor.cur.ki.with.kj)
                     closest.cor.index = kj.idx
                 }
                 
                } # end loop over the kj: now, did we get through them all without hitting redundancy?

               if (!flag.ki.redun) {
                   # Yes. So: we want this ki
                   max.snp[next.non.redundant,i]=current.snp.index
                   if (pm$show.snp.selection.log) {
                         pm.print(pm,paste("filterSNP.maxpriority.PC(): For trait",cn.non[i],
                                     ": Choosing (number ",next.non.redundant,") non-redundant greedy SNP",cn.snp[current.snp.index],
                                     "with snp-trait cor =",signif(signed.cor.snp.non[current.snp.index,i],3),
                                     "; The closest SNP was",signif(closest.cor,3),
                                     "correlated, that was SNP ",cn.snp[kj.idx]))
                   }
                   next.non.redundant = next.non.redundant+1
                   if (next.non.redundant > top.N.snps.per.trait) break; # we are done
                 }
             } # end for ki

           pm.print(pm,paste("=================== CONCLUSION OF (Greedy) SNP SELECTION DETAILS FOR ",cn.non[i]," === at ",date()))
           pm.cat(pm,paste(sep="","snps.for.",cn.non[i],"=c(\"",paste(collapse="\",\"",cn.snp[max.snp[,i]]),"\")\n"))
           pm.print(pm,      "================================================================================")

         } # end, if (top.N.snps.per.trait > 1)
                                        #
         } else {
            max.snp[,i] = scsn$ix[1:top.N.snps.per.trait];
         }

         max.snp.cor[,i]=cor.snp.non[max.snp[,i],i];
         max.snp.signed.cor[,i]=signed.cor.snp.non[max.snp[,i],i];
         max.names[,i]=cn.snp[max.snp[,i]];
         pc.namer[,i]=cn.non[i];
       }


        
    # max.snp has the numeric index of maximally correlated SNPs for each non.snp
    # and max.names names them.
    # max.names = cn.snp[max.snp]
    #

    nextra = max(orthog.search.depth - 1,0); # number of forward step-wise selected SNPs total
    
    pm.print(pm,paste("Done with greedy snp selection: ",date()))

    if (nextra > 0) {
       pm.print(pm,paste("Starting forward snp selection: ",date()))

       # ===============================================================
       # forward stepwise regression to pick up to (nextra) additional SNPs...
       # ===============================================================
       all.step.chosen.pc=c()
       all.step.chosen=c();
       all.step.chosen.cor=c();
       all.step.chosen.signed.cor=c();

       for (i in 1:ncol(cor.snp.non)) {

         # use the A,B edge filters if present to avoid picking SNPs for irrelevant traits.
         if (!is.null(pm$A) & !is.null(pm$B) & !(non.snp.cols[i] %in% pm$A | non.snp.cols[i] %in% pm$B)) {
            pm.print(pm,paste("Skipping Forward SNP selection for trait",cn.non[i],"since it wasn't on the pm$A or pm$B lists."))
            next;
         }

         pm.print(pm,paste(sep="","picking SNPs for trait ",i," out of ",ncol(cor.snp.non)," at ",date()))
         chosen= forward.step.pick.snp(df=x,ycol=cn.non[i],snpcols=cn.snp,single.start.col=max.names[1,i],nsteps=nextra, view.step);
         chosen.cor = cor.snp.non[chosen,i];
         chosen.signed.cor = signed.cor.snp.non[chosen,i];
         my.non = rep(cn.non[i],length(chosen));

         # and accumulate master lists
         all.step.chosen.pc=c(all.step.chosen.pc,my.non);
         all.step.chosen=c(all.step.chosen,chosen);
         all.step.chosen.cor = c(all.step.chosen.cor,chosen.cor);
         all.step.chosen.signed.cor = c(all.step.chosen.signed.cor,chosen.signed.cor);
       }

       pm.print(pm,paste("Done with forward snp selection: ",date()))
     } # end if nextra > 0
    
    # ==========================================================
    # now be sure to include max.snp as well as cluster "centers"...
    # ==========================================================

    # back to our regularly scheduled clustering....
    if (length(snpcols) > 2 && pm$do.snp.pick.hcluster) { # need minimum of 3 snps to have stuff to cluster

     # handle NA and infinite values....
     ddd=as.dist( 1-abs(t(cor(snps, method="p",use="pairwise.complete.obs"))));
     ddd[which(is.na(ddd))] = 1;
     ddd[which(is.infinite(ddd))] = 1;

       h1= hclust(ddd, method="average")
       #plot(h1);
       if (is.na(force.cutree.count)) {
            a=cutree(h1,h=cor.below.cluster.def);
            if (!is.unix()) {
                    cross.platform.windows();
                    plot(h1,main=paste("SNP selection: clusters below",cor.below.cluster.def,"line.")); 
                    abline(a=cor.below.cluster.def,b=0);
            }
       } else {
            a=cutree(h1,k=min(ncol(snps),force.cutree.count)); # force a specific number of groups of SNPs
            if (!is.unix()) {
               cross.platform.windows();
               plot(h1,main=paste("SNP selection: cluster forced count: ",force.cutree.count))
            }
       }

    cna=names(a);
    non.cna=colnames(x)[non.snp.cols];
    my.clusters=unique(a);
    my.snp.choices=my.clusters; my.snp.choices[]=0;
    my.snp.abscor=my.snp.choices; # store abs cor level for PC chosen.
    my.snp.signedcor=my.snp.choices; # store signed cor level for PC chosen.
    my.snp.pc =my.snp.choices; # store pc chosen.

    if (pm$do.snp.pick.hcluster) {
       pm.print(pm,paste(sep="","filterSNP.maxpriority.PC: choose ",length(my.clusters)," SNP clusters."));
    }

    num.non = length(non.snp.cols);
    x.non= x[,non.snp.cols];
    cn.x.non=colnames(x.non);
    k=1
    for (i in my.clusters) {
        cluster.i = cna[a %in% i];
        cluster.i.size = length(cluster.i);

        # choose SNP with maximum response in the non-snps
        data.i = snps[,cluster.i];
        signed.cor.data.i.x = cor(data.i, x.non, method="p",use="pairwise.complete.obs")
        cor.data.i.x = abs(signed.cor.data.i.x)

        # zero out the missing and infinite bad boys
        cor.data.i.x[which(is.na(cor.data.i.x))] = 0;
        cor.data.i.x[which(is.infinite(cor.data.i.x))] = 0;

        # now we want to find the single SNP with highest absolute correlation
        nrowhere = nrow(cor.data.i.x);

        wm = which.max(cor.data.i.x); # in column major (Fortran) matrix ordering
        best.row = ((wm-1)%% nrowhere)+1;
        best.col = 1 + ((wm-1)  %/% nrowhere);

        my.snp.choices[k]=cluster.i[best.row];
        my.snp.abscor[k]=cor.data.i.x[wm];
        my.snp.signedcor[k]=signed.cor.data.i.x[wm];

        my.snp.pc[k] = cn.x.non[best.col]
        #print(my.snp.pc)
        k=k+1
      } # end for i in my.clusters
    } else { # end bracket from if (length(snpcols) > 2 && pm$do.snp.pick.hcluster)
         a=snpcols; names(a)=colnames(x)[snpcols]; # just one SNP
    }

    # add the maximally correlated SNP for each trait, plus the cluster snps with maximal correlation...
    #   with max SNPs to come first.

    # accumulate in well named (or at least differently named) vectors
     pass        = c()
     snp.choice = c() # final version; versus with an s on the end, my.snp.choices<- was for the hclustering phase.
     snp.abscor  = c()
     snp.cor     = c()
     snp.pc      = c()

    # put fstep SNPs first, so that they are retained during redundancy checking
    # in preferance to the greedily chosen.
    # This shouldn't make much of a difference if our redundancy threshold is
    # low enough, but might if we aren't stringent and only reject redundants at .95 and above.
    #
    if(pm$do.snp.pick.fstep & nextra > 0) {
     pass   = c(pass,rep("fstep",length(all.step.chosen)))
     snp.choice = c(snp.choice,all.step.chosen)
     snp.abscor  = c(snp.abscor,signif(c(all.step.chosen.cor),3))
     snp.cor     = c(snp.cor,signif(c(all.step.chosen.signed.cor),3))
     snp.pc      = c(snp.pc,all.step.chosen.pc)
    }

    if(pm$do.snp.pick.greedy) {
     pass   = c(pass,rep("greedy",length(max.names)))
     snp.choice = c(snp.choice,max.names[1:length(max.names)])
     snp.abscor  = c(snp.abscor,signif(c(max.snp.cor[1:length(max.snp.cor)]),3))
     snp.cor     = c(snp.cor,signif(c(max.snp.signed.cor[1:length(max.snp.cor)]),3))
     snp.pc      = c(snp.pc,pc.namer[1:length(pc.namer)])
    } 

    if(pm$do.snp.pick.hcluster) {
     pass   = c(pass,rep("hclustSNPs",length(my.snp.choices))) 
     snp.choice = c(snp.choice, my.snp.choices)
     snp.abscor  = c(snp.abscor,signif(c(my.snp.abscor),3))
     snp.cor     = c(snp.cor,signif(c(my.snp.signedcor),3))
     snp.pc      = c(snp.pc,my.snp.pc)
    } 

    # now we can rename for backwards compatbility
    my.pass        = pass
    my.snp.choice  = snp.choice
    my.snp.abscor  = snp.abscor
    my.snp.cor     = snp.cor
    my.snp.pc      = snp.pc

    nrow.my.snp = length(my.pass)
    redundant=rep(FALSE,nrow.my.snp)
    redundant.note=rep(" ",nrow.my.snp)
    redundant.prior.snp = rep(NA,nrow.my.snp) # the one we keep instead

    # ==========================================================================
    # We may have some SNPs that are still so highly correlated with each other
    # that we need to discard one or more to avoid multi-collinearity.
    # Toward that end, discard any later SNPs that have too strong a correlation (>= pm$snp.redundancy.cor.level)
    # with earlier SNPs in our selection process;
    # ==========================================================================

    # To get a unique list of chosen SNPs, don't do this:
    # uniq.snp = unique(my.snp.choice); 
    # Instead, we prefer to use duplicated() so we can keep the ordering of the SNPs stable.
    # We don't know that unique() is stable (Maybe it is, maybe it's not. It's not documented in the R documentation, so
    # we shouldn't count on it.)
    has.dup = duplicated(my.snp.choice)
    redundant[has.dup]=TRUE # repeats found by both greedy and forward (or hcluster) not needed.
    uniq.snp = my.snp.choice[!has.dup]

    uniq.snp = uniq.snp[!is.na(uniq.snp)] # remove NA
    uniq.snp = uniq.snp[uniq.snp!="0"] # remove "0"

    usnp = x[uniq.snp];
    ucn=colnames(usnp);
    cor.usnp = abs(cor(usnp,method="pearson",use="pairwise.complete.obs"));
    if (ncol(cor.usnp) > 1) {
    for (i in 1:(nrow(cor.usnp)-1)) {
        # and don't bother checking for redundancy with an already redundant SNP...
        if (redundant[i]) next;
        for (j in (i+1):(ncol(cor.usnp))) {
            if (cor.usnp[i,j] >= pm$snp.redundancy.cor.level) {  
               # j is redundant with i
               primary.row = which(my.snp.choice==ucn[i]);
               red.result.rows=which(my.snp.choice==ucn[j]);
               for (k in red.result.rows) {
                   redundant.note[k] = paste(sep="","Retracted: ",ucn[j]," was redundant w/ ",ucn[i]," at cor: ",signif(cor.usnp[i,j],2));
                   redundant[k]=TRUE
                   redundant.prior.snp[k] = ucn[i]
               }
               if (!pm$skip.snp.retraction) { my.snp.choice[red.result.rows]=NA; }
            }
    }}
  } # end if ncol(cor.usnp) > 1

 non.redundant = !redundant

 df = data.frame(my.snp.choice,my.snp.abscor, my.snp.cor,my.snp.pc,my.pass,non.redundant,redundant.note,redundant.prior.snp)
 # remove any NA that snuck in and return that frame
 df[!is.na(my.snp.choice),]
}




# find the supernodes in the non.snps (here called traits, but includes PCs and Traits.)
#
# PRE: ncol(traits) >= 2
#
if(exists("supernode.traits") ) rm(supernode.traits);
supernode.traits=function(traits, cor.above.supernode.def=.7,pm) {
    if (ncol(traits) < 2) stop("supernode.traits() called with fewer than 2 column traits. Can't do that.")
    bel=1-cor.above.supernode.def;
    h1= hclust(as.dist( 1-abs(t(cor(traits, method="p",use="p")))), method="average" );
    if (ncol(traits) >2 & !pm$quiet) {
         cross.platform.windows();
         if (!is.unix()) {
             plot(h1,main=paste("Supernodes below",bel,"line."));
             abline(a=bel,b=0);
         }
         a=cutree(h1,h=bel);
       } else { # handle the case of 2 traits
             a=1:ncol(traits);
             attr(a,"names")=colnames(traits)
             if (cor(traits[,1],traits[,2],method="p",use="p") >= cor.above.supernode.def) {
               a=c(1,1); # put the two columns in the same supernode
               attr(a,"names")=colnames(traits)
             }
       }
    a
}


if(exists("standalone.pick.snp") ) rm(standalone.pick.snp);
standalone.pick.snp=function(x,snpcols=c(), traitcols=c(),pm) {

   nc=ncol(x);
   cn=c(colnames(x)," "); # add extra blank column for printing null LCD outputs

   # Guess at Trait if not given.
   # We want traits that 
   if (length(traitcols)==0) {
       traitcols = c(which(tolower(cn)=="trait"),which(tolower(cn)=="y"));
       if (length(traitcols)==0) {
          pm.print(pm,"Hmmm: traitcols not given and no data columns start with 'Trait' or 'y'"); }
       pm.print(pm,"traitcols not specificed manually, looking for 'trait' or 'y' column-names. Using:");
       pm.print(pm,cn[traitcols]);
   }

   traits = cn[traitcols];


   # guess at snpcols if not given
   if (length(snpcols)==0) { snpcols=which(tolower(substr(cn,1,3)) == "snp");
               if (length(snpcols)==0) { 
                 stop("snpcols not given and no data columns start with 'SNP': Must have SNPs to run current codebase without elaborate VCdep computation.");} else {
                 pm.print(pm,"snpcols not specified manually, looking for 'SNP' prefixed column names"); 
                 if (length(snpcols) < 100) {
                    pm.print(pm,paste("   ...we'll suppose that these ",length(snpcols)," are our markers/SNPS:"));
                    pm.print(pm,cn[snpcols]); 
                  } else {
                    pm.print(pm,"...more than 100 snpcols....skipping printout.")
                 }
}}

   numsnps=length(snpcols);

   non.snp.cols = setdiff(1:nc,snpcols);
   num.non.snps = length(non.snp.cols);

  if (is.null(pm$top.N.snps.per.trait)) pm$top.N.snps.per.trait=12;

  my.snp = filterSNP.maxpriority.PC(snpcols, non.snp.cols, x, cor.below.cluster.def=pm$snp.below, force.cutree.count=pm$force.cutree.count,orthog.search.depth=pm$orthog.search.depth,view.step=F,top.N.snps.per.trait=pm$top.N.snps.per.trait,pm=pm);
  save(my.snp,file="standalone.pick.snp.most_recent.rdat");
  my.snp
}

# test data for impute.nn
#x1=matrix(c(1,2,3,1,2,3,1,2,3),nrow=9);
#x2=matrix(c(1,2,NA,1,NA,3,1,2,3),nrow=9);
#x3=matrix(c(1,2,NA,4,NA,6,7,8,9),nrow=9);
#x4=matrix(c(1,2,3,1,2,3,1,2,4),nrow=9);
#x5=matrix(c(1,2,NA,1,NA,3,1,2,5),nrow=9);
#x6=matrix(c(1,2,NA,4,NA,6,7,8,6),nrow=9);
#x=data.frame(x1,x2,x3,x4,x5,x6)
#snpcols=1:2

#IMPUTE MISSING DATA USING NEAREST.NEIGHBORS -- only impute snpcols
if(exists("impute.nn") ) rm(impute.nn);
impute.nn=function(x,snpcols,pm) {
    nc=ncol(x);
    cn=colnames(x);
    rn=rownames(x);

    if (length(snpcols) < 1) stop("snpcols cannot be empty in impute.nn().")

    # block out the non-snpcols so we don't accidentally choose them
    # as the nearest-neighbors for computation if x has lousy nearest-neighbors
    x2=x[,snpcols]

    # screen out factors in snpcols, because we want alleles numeric and ordered for F2 crosses
    # and factors as default read in may not be ordered
    for (i in snpcols) {
       if (is.factor(x[,i])) stop("The specified snpcols must be numeric rather than factors so that correlations are computed correctly in impute.nn().")
    }

    # rpart imputation was toooooo sloooooowwww. Do Nearest-Neighbor instead.
    corx2=cor(x2, method="p",use="pairwise.complete.obs");

    # mock up the full correlation matrix including non-snps (even though we didn't compute it)
    # just so we have *sane* indexing again
    corx = matrix(nrow=nc, ncol=nc, data=-100) # put non-SNP columns last, with correlation of -100
    colnames(corx) = rownames(corx) = cn
    corx[snpcols,snpcols] = corx2

    impute.log=gsub(":",".",gsub(" ","_",paste(sep="","nn.imputation.na.log.",attr(x,"title"),date(),".txt")))
    sink(file=impute.log,type="output");
    pm.print(pm,paste(sep="","'",impute.log,"' is this file:"));
    for (j in 1:nc) {
          if ((j %in% snpcols) && any(is.na(x[,j]))) {
              pm.print(pm,paste("impute.nn(): imputing column: ",cn[j],date()));
        my.neighbors = corx[j,];
              my.neighbors[j]=-2; # put self last
              ix = sort(my.neighbors,decreasing=T,index.return=T)$ix;
              for(i in which(is.na(x[,j]))) {
                  #INVAR: x[i,j] is NA. So first check x[i,closest neighbor of column j]
                  ki=0;
                  for (k in ix) {
                      ki=ki+1;
                      if (!is.na(x[i,k])) {
                          if (my.neighbors[k] < pm$impute.nn.minimum.correlation.th) {
                               pm.print(pm,paste(sep="","SKIPPING IMPUTATION by NN because correlation ",signif(my.neighbors[k],3)," < ",
                                        pm$impute.nn.minimum.correlation.th," pm$impute.nn.minimum.correlation.th wasn't ",
                                        "reached for the nearest neighbor. For x[",
                                        i,"(",rn[i],"),",j,"(",cn[j],")]=NA from ",
                                        ki,"-th NN (cor=",signif(my.neighbors[k],3),"): x[",i,"(",rn[i],"),",k,"(",cn[k],")]=",x[i,k]));
                               break;                    
                          }
                          x[i,j] = x[i,k]; 
                          pm.print(pm,paste(sep="","NN imputing x[",i,"(",rn[i],"),",j,"(",cn[j],")]=NA from ",
                               ki,"-th NN (cor=",signif(my.neighbors[k],3),"): x[",i,"(",rn[i],"),",k,"(",cn[k],")]=",x[i,k]));
                          break;
                      }
                  }
              }
          }
    } # end j in nc
    sink(); # manually turn off log before showing file; otherwise redundant.
    file.show(impute.log);
    data.frame(x)
}
#xn=impute.nn(x)

#IMPUTE MISSING DATA USING RPART
if(exists("impute.rpart") ) rm(impute.rpart);
impute.rpart=function(x,snpcols) {
    x=data.frame(x); # label unlabelled x columns by converting to data frame
    cn=colnames(x);

    if (length(snpcols)==0) {
       snpcols=which(tolower(substr(cn,1,3)) == "snp");
    }

    # test: x=matrix(c(1,2,NA,4,NA,6,7,8,9),nrow=3);
    nc=ncol(x);
    cn=colnames(x);

#    # DON'T DO THIS: it massively slows down rpart(), and isn't needed
#    print(paste("impute.rpart(): FACTOR-izing columns ",date()));
#    # make sure SNP columns are treated as factors
#    for (i in 1:nc) {
#        print(paste("Factorizing column",i));
#        if (i %in% snpcols) {
#     if (i==1) {
#                x=data.frame(factor(x[,i]),x[,-i]);
#            } else if (i < nc) {
#                x=data.frame(x[,1:(i-1)],as.factor(x[,i]),x[,(i+1):nc]);
#            } else {
#                x=data.frame(x[,1:(i-1)],as.factor(x[,i]));
#            }
#        }
#    }
#    colnames(x)=cn; # restore in case messed up by factor-izing
#    print(paste("impute.rpart(): DONE FACTOR-izing columns! ",date()));

    xnew=x; # before NA fill, and don't use filled values that were NA to predict other NA values (too much momentum in one direction)
    # fill in missing SNPs
      for (i in 1:length(cn)) 
          if (any(is.na(x[,i]))) {
              pm.print(pm,paste("impute.rpart(): imputing column: ",cn[i],date()));
              if (i %in% snpcols) {
                 # predict this SNP column
                 f.upper <- as.formula(paste(cn[i], "~", paste(cn[-i], collapse = " + ")))
                 rp=rpart(f.upper,data=x,method="class",control=rpart.control(cp=.0001,minsplit=2,minbucket=1,maxsurrogate=10))
           xnew[,i]=predict(rp,newdata=x,type="class"); #set newdata=x to get predictions on the NA columns too
              } else {
                 # predict this qualitative column
                 f.upper <- as.formula(paste(cn[i], "~", paste(cn[-i], collapse = " + ")))
                 rp=rpart(f.upper,data=x,method="anova",control=rpart.control(cp=.0001,minsplit=2,minbucket=1,maxsurrogate=10))
           xnew[,i]=predict(rp,newdata=x,type="vector"); #set newdata=x to get predictions on the NA columns too
              }
      }
      
   colnames(xnew)=cn; # restore in case messed up by assignments from predict()
   xnew
}
# test calls for impute.rpart
# x=impute.rpart(x,snpcols;)
# if (any(is.na(x))) { x=impute.rpart(x,snpcols); }


# Network Edge Orienting. (NEO)
# 
# We break down neo() into five component functions. 
#
# Each of these is a separate function:
#
# 1. SNP primary, neo.snp.primary()
#
# 2. SNP election, neo.snp.election()
#
# 3. SNP based edge orienting, neo.snp.edge.orient()
#
# 4. Collider-based edge orienting, neo.collider.orient()
#
# 5. Model Fitting Indices and Statistics, neo.overall.fit()

#
# CLASS: "NEO" objects will hold intermediate results along the way
#



# 1. SNP primary, neo.snp.primary()
#
# Here we do the initial SNP picking, but not the final SNP election.
#
# If specified, no.obs over-rides nrow(x)
#
if(exists("neo.snp.primary") ) rm(neo.snp.primary);
neo.snp.primary=function(datCombined, snpcols=neo.guess.snps(datCombined,pm), traitcols=neo.guess.traits(datCombined,pm), pm=neo.get.param(),
                         skip.snp.select=F,no.obs=NULL,no.snp.impute=F, neo.log.file="neo.log",...) {

   r=list(); # our returned and circulated result list, holds parameters to be passed between the neo. functions

   r$neo.log.file = neo.log.file

   nc=ncol(datCombined);
   cn=c(colnames(datCombined)," "); # add extra blank column for printing null LCD outputs

   orig.traitcols.cn=cn[traitcols];
   traits = cn[traitcols];
   orig.snpcols.cn=cn[snpcols];
   orig.cn = cn

   datCombined.snpfilled=datCombined; # only use imputed SNPs for SNP selection, nothing else...

   # impute missing data by predicting using the other columns...in the non-ignored columns.
   if (!no.snp.impute && !skip.snp.select && any(is.na(datCombined[,snpcols]))) {
       pm.print(pm,paste("imputing missing values in SNPs...",date()));
       #datCombined=impute.rpart(datCombined,snpcols); 
       datCombined.snpfilled=impute.nn(datCombined,snpcols,pm); 
       pm.print(pm,paste("done with imputation of NA for SNPs...",date()));
      }


   numsnps=length(snpcols);

   non.snp.cols = setdiff(1:nc,snpcols);
   num.non.snps = length(non.snp.cols);

   ### ======================================
   ### use SNP clustering function
   ### ======================================

   snps=datCombined[,snpcols]

   # handle the case where snpcols==1 and we need to re-convert to matrix
   if (numsnps == 1) {
     snps = matrix(snps,nrow=nrow(datCombined))
   }

   if (numsnps > 0) {
     if (skip.snp.select) {
        pm.print(pm,"Skipping SNP selection and using any (all) provided SNPs.")
        my.snp=NULL;
     } else {

     pm.print(pm,paste("choosing SNPs...",date()));

     if (is.null(pm$top.N.snps.per.trait)) pm$top.N.snps.per.trait=12; # reasonable default

     my.snp = filterSNP.maxpriority.PC(snpcols, non.snp.cols, datCombined.snpfilled, cor.below.cluster.def=pm$snp.below, force.cutree.count=pm$force.cutree.count,view.step=F,orthog.search.depth=pm$orthog.search.depth,top.N.snps.per.trait=pm$top.N.snps.per.trait,pm=pm);

     if (pm$do.snp.pick.hcluster && !pm$do.snp.pick.greedy && !pm$do.snp.pick.fstep) {
        pm.print(pm,"************ SPECIAL: just doing hclustering SNPS ONLY **********");
        my.snp=my.snp[(my.snp$my.pass == "hclustSNPs"),];
        pm.print(pm,paste(sep="","Using hclustSNPs, count=",nrow(my.snp)));
     }

     pm.print(pm,paste("Done choosing SNPs...",date()));

     if (!pm$no.log) save(my.snp,file=paste(sep="",neo.log.file,".my.snp.choices.rdat"));
     my.snp.choices = sort(unique(my.snp$my.snp.choice));

     my.snp.choices=my.snp.choices[!is.na(my.snp.choices)]; # eliminate NA

     my.snp.real = match(my.snp.choices, cn); # re-index by the original columns

     pm.print(pm,my.snp)

     snpcols=my.snp.real; # replace old defn of snpcols with cluster center SNPs
     numsnps=length(snpcols);

    } # else not skip.snp.select

     # As long as numsnps > 0, 
     # rearrange so SNPs are all at the beginning, just to make plotting easier later...
     new.datCombined = datCombined[,c(snpcols,non.snp.cols)]; 
     datCombined=new.datCombined; rm(new.datCombined);
   } # if numsnps > 0


   # ================================================================================
   # re-generate the column info that changed...when we re-arranged the columns of datCombined
   #=================================================================================

   # and pack r
   r$nc = nc = ncol(datCombined);
   if (is.null(no.obs)) { r$no.obs=nrow(datCombined); } else { r$no.obs = no.obs}
   r$cn = cn = c(colnames(datCombined)," ");
   snpcols = match(orig.snpcols.cn,cn);
   snpcols = snpcols[!is.na(snpcols)];
   r$snpcols = snpcols;
   r$non.snp.cols = non.snp.cols = setdiff(1:nc,snpcols);
   r$num.non.snps = length(non.snp.cols);

   # update the pm$ignorable.edge.list if present
   if (!is.null(pm$ignorable.edge.list)) {
      for (i in 1:length(pm$ignorable.edge.list)) {
         a=pm$ignorable.edge.list[[i]]
         if (length(a)==2) {
            a[1] = match(orig.cn[a[1]],cn)
            a[2] = match(orig.cn[a[2]],cn)
            pm$ignorable.edge.list[[i]] = a
         }
      }
   }

   traitcols = match(orig.traitcols.cn,cn);
   r$traitcols = traitcols[!is.na(traitcols)];

   r$snps = snps;
   r$datCombined = datCombined;
   r$pm = pm;
   r$my.snp = my.snp
   r$numsnps = numsnps

   class(r)="NEO.after.snp.primary";
   r
} # end neo.snp.primary


# 2. SNP supernoding, neo.create.super.nodes()
if(exists("neo.create.super.nodes") ) rm(neo.create.super.nodes);
neo.create.super.nodes=function(r,...) {

   if (!inherits(r,"NEO.after.snp.primary")) {
      stop("neo.create.super.nodes() stopping because r (1st argument) must inherit from class NEO.after.snp.primary");
   }

    # unpack r
    nc=r$nc
    cn=r$cn
    no.obs=r$no.obs
    snpcols=r$snpcols

    non.snp.cols=r$non.snp.cols
    traitcols=r$traitcols
    numsnps = r$numsnps
    num.non.snps = r$num.non.snps
    snps=r$snps

    pm=r$pm
    my.snp=r$my.snp
    datCombined=r$datCombined

    # ========================================
    # the the correlation matrix
    # ========================================

    # check for completely redundant columns

    xcor=cor(datCombined, use="pairwise.complete.obs",method="p");
    ##xcor=cor(datCombined,method="spearman",use="pairwise"); # outlier robust version
    cx=abs(xcor);

    # =======================================
    # find any SuperNodes among the traits, where we suspend the conditional independence
    #  filters to avoid mutual cancellation of reasonable edges.
    # =======================================
    pctraits= datCombined[,non.snp.cols];
    pre.super = supernode.traits(pctraits, cor.above.supernode.def=pm$cor.above.supernode.def,pm=pm) 

    if (numsnps > 0) {  
        super = c(1:length(snpcols), pre.super+length(snpcols)); # put each snpcol in it's own supernode, too.
    } else {
        super=pre.super;
    }

    pm.print(pm,"SUPERNODE SUMMARY:");
    for (sup in unique(super)) {
        wsup = which(sup == super);
        pm.print(pm,paste("Supernode ", sup,":", paste(cn[wsup],collapse=", ")));
    }

    # ==============================================================
    # create an un.filter to later restore the connections between SuperNode members
    # ==============================================================
    un.super=cx; un.super[]=0;

    super.assigns=names(pre.super);
    my.supers=unique(pre.super);
    for (i in my.supers) {
        super.i = super.assigns[pre.super %in% i];
        super.i.wch = which(cn %in% super.i);
        un.super[super.i.wch, super.i.wch]=1;
    }
    # and zero the diagonal
    for (i in 1:ncol(un.super)) {un.super[i,i]=0;}

    # ===========================================================
    # Exclude the trait column(s) from any supernode.
    # ===========================================================
    un.super[,traitcols]=0;
    un.super[traitcols,]=0;
    super[traitcols]=max(pre.super)+length(snpcols)+traitcols; # some high number so traits will each have their own supernode.

    # =====================================================
    # Done with supernodes.
    # =====================================================
    pm.print(pm,paste("Done creating supernodes...",date()));

    r$super=super
    r$un.super=un.super
    r$cx = cx

    class(r) = "NEO.after.create.super.nodes"
    r

} # end neo.create.super.nodes


# 2-5 Filtering
if(exists("neo.edge.filters") ) rm(neo.edge.filters);
neo.edge.filters=function(r,...) {

   if (!inherits(r,"NEO.after.create.super.nodes")) {
      stop("neo.edge.filters() stopping because r (1st argument) must inherit from class NEO.after.create.super.nodes");
   }

    # unpack r
    nc=r$nc
    cn=r$cn
    no.obs=r$no.obs
    snpcols=r$snpcols

    non.snp.cols=r$non.snp.cols
    traitcols=r$traitcols
    numsnps = r$numsnps

    num.non.snps = r$num.non.snps
    snps=r$snps

    pm=r$pm
    my.snp=r$my.snp
    datCombined=r$datCombined

    super=r$super
    un.super=r$un.super
    cx = r$cx

   # =====================================================
   # start creating and using the filters
   # =====================================================

## redone above: pm$A and pm$A.cn will be consistent, but
##  we are tracking the global batch processing in pm$enumerated.A.cn
##  (and pm$enumerated.B.cn) so that we don't want to do that again here.
##
   ABfilter = (cx > 100) # start with an all false, labelled matrix
   ABfilter[] = TRUE # filter out nothing # by default

##    ABfilter[] = FALSE
##    if (!is.null(pm$A) & !is.null(pm$B)) {
##       # rematch after shrinking of matrix
##       pm$A = match(pm$A.cn,colnames(cx))
##       pm$B = match(pm$B.cn,colnames(cx))

##       pm.print(pm,"Using pm$A and pm$B filters.");

##       # next line depends upon prior call to the effect of: pm=get.enumab.cn.set(pm)
##       # Now we enumerated specifically the entries to consider.
##       if (is.null(pm$enumerated.A.cn)) {
##          stop("pm$enumerated.A.cn not set, so pm=get.enumab.cn.set(pm) not called prior to here. Halting!")
##       }
##       for (i in 1:length(pm$enumerated.A.cn)) {
##         ABfilter[pm$enumerated.A.cn[i],pm$enumerated.B.cn[i]]=TRUE
##         ABfilter[pm$enumerated.B.cn[i],pm$enumerated.A.cn[i]]=TRUE
##       }
##    } else {
##       ABfilter[] = TRUE # filter out nothing
##    }

   # useful default filter value: ones
   ones = matrix(ncol=nc,nrow=nc);
   ones[,]=1;

   # a logical filter lcx
   lcx = (cx>pm$cor.dep.th); # keep large correlation variables

   # snps have no children, so zero those columns.
   # Later, use exog to toss out variables.
   exog = cx; exog[]=1;
   if (numsnps > 0) { exog[,snpcols] = 0; }# ; zeroing column j means that var j has no parents. 
   for (i in 1:nc) { exog[i,i]=0; }

   # create a snp.selection.filter that restricts the SNP parents of a variable
   # to those chosen by the SNP selection process, as indicated by the r$my.snp data frame
   chosen.snp.filter=cx; chosen.snp.filter[]=1
   if (!is.null(my.snp)) {
     for (v in non.snp.cols) {
        v.snp = unique(my.snp$my.snp.choice[which(my.snp$my.snp.pc==cn[v])])
        w.row = match(v.snp,cn)
        non.w.row = setdiff(snpcols,w.row) # of all the snp cols, pick out the non-selected ones for this variable v
        chosen.snp.filter[non.w.row,v]=0   # and zero those rows for this column v.
     } 
   }

   # fill in any NA values before calling pcor.shrink, because it will crash
   # use the median, from package e1071
   if(any(is.na(datCombined))) {
      datCombined = as.data.frame(impute(datCombined,what="median"))
   }

   px=cx
   px[]=0 # off by default
   if (pm$skip.pcor.filters.currently) {
      px=abs(pcor.shrink(datCombined,verbose=FALSE)) # use Strimmer's well-behaved function based on pseudo-inverse
   }
   lpx = (px >= pm$pcor.th) # keep variables whose partials are large
   
   #compute Cohen's omega coefficient filter
   pxpc = abs(px)/abs(cx);
   omega = abs(ones-pxpc);
   lomega = (omega <= pm$omega.th); # keep small omega variables=>small indirect influences
   
   ### Jason's omega or jomega will be the ratio of the minimum k pcor to the total correlation
   jomega = omega; jomega[]=0;
   
   # computing R.squared on 50K column data frames doesn't work
   small.data=FALSE
   if (small.data) {
     # establish the queue priority based on largest r.squared. 
     R=vector(length=nc)
     for (i in 1:nc) {
              R[i]=summary(lm(as.formula(paste(cn[i],"~.")),data=datCombined))$r.squared; 
             #holds the r.squared for each variables regressed on all others
     }   
   
     q=sort(R,index.return=T,decreasing=T) # q$ix  lists the variables in order of priority 
     # for being considered as child nodes....using Cohen's heuristic that it's better to have 
     # children that have large R-squared when regressed on the rest of the variables.
   
     rev=sort(q$ix,index.return=T); # give us a map to where the strongest R-squared vars are
     # by indexing:  rev$ix[i] gives the rank of i in terms of its R-squared when regressed on all others.
   
     Rrank = substr(paste(sep="","[",rev$ix,"] ",R),1,10);
     R.sq= matrix(nrow=nc,byrow=TRUE,rep(Rrank,nc))
   } else {
     R.sq = NULL
   }

   
   # filter on single variable conditional independence checks: elminate the direct
   # link between i and j if there is a variable k which renders cor(i,j|k) almost zero, or below k.pcor.th
   
   covx=cov(datCombined)
   pcor.filtered = covx; pcor.filtered[,]=1;
   
   #Vcdep needed for non-SNP RVV checks, but not otherwise
   #Vcdep=array(dim=c(nc,nc,nc)); Vcdep[]=0; # also check for induced 
                      # conditional dependence, while we are enumerating triples,
   # and computing partial correlations, re-cycle the work in pcijk by noting at the same time if we could 
   # have an  unshielded collider (v-structure) present by setting cdep[i,j].
   
   min.k.abs.pcor=cx;
   min.k.abs.pcor[]=1;
   min.k.abs.pcor.ignore.supernode = min.k.abs.pcor;
   min.k=cx; min.k[]=NA;
   min.k.nosup=min.k;
   jomega.nosup = jomega;
   
   if (pm$k.pcor.th <= 0) {
      pm.print(pm,paste("skipping kmin/Vcdep computation loops...b/c pm$k.pcor.th <=0 ...",date()));
   } else {
      pm.print(pm,paste("starting kmin/Vcdep computation loops...",date()));
      for(i in 1:(nc-1)) {
        pm.cat(pm,paste(sep="",cn[i],"..."));
        for (j in (i+1):nc) {
          for (k in 1:nc) {
            if (k !=i && k != j) {
   
                  # do a version of min.k that ignores supernodes too.
                  pcijk2 = abs(pcor(c(i,j,k),covx));
                  if(pcijk2 < min.k.abs.pcor.ignore.supernode[i,j]) {
                      min.k.abs.pcor.ignore.supernode[i,j] = min.k.abs.pcor.ignore.supernode[j,i] = pcijk2; 
                      min.k.nosup[i,j]=k; min.k.nosup[j,i]=k;
                      jomega.nosup[i,j] = jomega.nosup[j,i] = pcijk2 / cx[i,j]; 
                 }
   
                  # don't test traits from the same SuperNode against each other
                  if (super[k] == super[i] || super[k]==super[j]) next;
   
                  pcijk = abs(pcor(c(i,j,k),covx));
                  if(pcijk < min.k.abs.pcor[i,j]) {
                      min.k.abs.pcor[i,j] = min.k.abs.pcor[j,i] = pcijk; 
                      min.k[i,j]=k; min.k[j,i]=k;
                      jomega[i,j] = jomega[j,i] = pcijk / cx[i,j]; 
                 }
      # take out VCdep for now
      #  we ARE using the min.k.abs.pcor, AND the pcor.filtered[,] so save the loopage.
             #print(paste(i,j,pcijk, pcijk > pm$induced.dep.th, 
      #                     #cx[i,j] <= pm$cor.ind.th, cx[i,j]));
                  if (pcijk < pm$k.pcor.th) {
                 pcor.filtered[i,j] = pcor.filtered[j,i]=0;
                  }
   #               if (pcijk >= pm$induced.dep.th && cx[i,j] <= pm$cor.ind.th) {
   #                   Vcdep[i,j,k]= rev$ix[k]; Vcdep[j,i,k]= rev$ix[k];
   #                   # possible v-structure parents of k,
   #                   # we do a symmetric pair of assignments 
   #                   # becasuse i-j order is irrelevant; 
   #                   # and the rev$ix[k] value gives the rank 
   #                   # in terms of preferred child status
   #               }
                }
             }
        }
    }
    pm.cat(pm,"\n");
    pm.print(pm,paste("Done with kmin/Vcdep computation loops...",date()));
    # INVAR: Vcdep holds all possible v-structure (unshielded collider) triples
   }
   
   ljomega = (jomega <= 1- pm$omega.th); # keep large jomega variables=>small indirect influences


   # now we want to keep all edges that are into and out of trait columns if they pass the
   # lcx filter--and then ignore the lpx and lomega and pcor.filtered filters for traits
   restore.trait.edges=lcx;
   restore.trait.edges[]=0;
   restore.trait.edges[traitcols,]=lcx[traitcols,];
   restore.trait.edges[,traitcols]=lcx[,traitcols];
   
   restore.trait.edges = restore.trait.edges & exog; # don't restore edges into SNPs
   
   # combine all filters; convert from logical to 0-1 matrix; zero the diagonal
   if (!pm$skip.pcor.filters.currently) {
      after.filter=(((lcx & lpx & lomega & pcor.filtered & exog & ABfilter & chosen.snp.filter) | restore.trait.edges) | un.super )+0; 
   } else {
      # much simpler scheme, no correlation filtering...
      after.filter=(((exog & ABfilter & chosen.snp.filter) | restore.trait.edges) | un.super )+0; 
   }
 
   af=after.filter;
   for (i in 1:nc) { af[i,i]=0; }

   # now apply the pm$only.score.edges.into  requirements: only bother to score
   # edges into this set, if given, i.e. 

   # override everything except SNP selection if we have pm$only.score.edges.into  columns
   if (!is.null(pm$only.score.edges.into)) {
      # rows=from, cols=to; for each for, the only column(s) that should have any 1's are
      # given by the pm$only.score.edges.into set.
      set.only.edges.into.matrix=lcx;
      set.only.edges.into.matrix[]=0;
      w = match(pm$only.score.edges.into,cn)
      w = w[!is.na(w)]
      af[non.snp.cols,w]=1
      af[,-w]=0
   }

   r$af = af;
   r$exog = exog;  
   r$lomega = lomega;
   r$omega = omega;
   r$lcx = lcx;
   r$lpx = lpx;
   r$px = px;
   r$R.sq = R.sq;
   r$pcor.filtered = pcor.filtered;
   r$min.k.abs.pcor = min.k.abs.pcor;
   r$ljomega = ljomega;
   r$jomega = jomega;
   r$min.k = min.k;
   r$min.k.abs.pcor.ignore.supernode = min.k.abs.pcor.ignore.supernode;
   r$min.k.nosup = min.k.nosup;
   r$covx = covx;
   r$datCombined = datCombined; # keep imputed data version

   class(r)="NEO.after.neo.edge.filters";
   r   

} # neo.edge.filters


# 3. SNP based edge orienting, neo.snp.edge.orient()
if(exists("neo.edge.orient") ) rm(neo.edge.orient);
neo.edge.orient=function(r,fit.full.model=FALSE, repos.nodes=FALSE, skip.LEO=FALSE, ignorable.edge.list, ...) {

   if (!inherits(r,"NEO.after.neo.edge.filters")) {
      stop("neo.edge.orient() stopping because r (1st argument) must inherit from class NEO.after.neo.edge.filters");
   }

    # unpack r
    datCombined=r$datCombined
    nc=r$nc
    cn=r$cn
    no.obs=r$no.obs
    snpcols=r$snpcols

    non.snp.cols=r$non.snp.cols
    numsnps = r$numsnps

    num.non.snps = r$num.non.snps;
    pm=r$pm;

    super=r$super;
    un.super=r$un.super;
    cx = r$cx;
    af=r$af;

    exog = r$exog;
    lomega = r$lomega;
    omega = r$omega;
    lcx = r$lcx;
    lpx = r$lpx;
    px = r$px;
    R.sq = r$R.sq;
    pcor.filtered = r$pcor.filtered;
    min.k.abs.pcor = r$min.k.abs.pcor;
    ljomega = r$ljomega;
    jomega = r$jomega;
    min.k = r$min.k;
    min.k.abs.pcor.ignore.supernode = r$min.k.abs.pcor.ignore.supernode;
    min.k.nosup = r$min.k.nosup;
    covx = r$covx;
    neo.log.file = r$neo.log.file

    # done unpacking r

   # no, let user decide...  options(warn=1); # print warnings as the occur.

   # choose and log which orienting we are using
   eo.choice=c("forward.step.selection","all.snp.permutations.equal.votes","max.vs.max")
   eo.type=eo.choice[pmatch(pm$eo.type,eo.choice)];

   if (is.na(eo.type)) { eo.type= "all.snp.permutations.equal.votes"; }
   pm.print(pm,paste("eo.type = ",eo.type))

   pm.print(pm,paste("All starting neo.edge.orient()...",date()));

   op = options(); # to restore digits at the end

   # grab title from x before any conversion
   z=list(); # the results at the end; what is returned.
   z$call = match.call(); # save call for later reference

   if (!is.null(pm$run.title)) { #(!is.null(attr(datCombined,"title"))) {
      z$my.title=pm$run.title # attr(datCombined,"title");
      title.string=paste("================================================= Model title: ",z$my.title," ===");
      pm.print(pm,paste(rep("=",nchar(title.string)),collapse=""));
      pm.print(pm,title.string);
      pm.print(pm,paste(rep("=",nchar(title.string)),collapse=""));
   } # name the model, so neo.graph can too.



M=af*0; # our empty model

zeo=M; zeo[]=0; # track the sum of all the z-scores for the component edge tests in CCO/CEO
nzeo=M; nzeo[]=0; # count total number of z-scores we add

lcd.edge =M; lcd.edge[]=nc+1; # note where LCD finds edges
 lcd.drop =M; lcd.drop[]=0;# and note value of conditional indepence drop: cor(m1,B) - cor(m1,B | A)
lcd.possibly.erase = M; lcd.possibly.erase[]=FALSE; # note where the assymetry should be if it's not learned any
 # other way.

# implement SNP->one->two checks....twofer checks, and store stats here
twofer.eo=M; twofer.eo[]=0; # track the sum of all the z-scores for the component edge tests in CCO/CEO
n.twofer.eo=M; n.twofer.eo[]=0; # count total number of z-scores we add
twofer.log=matrix(as.character(n.twofer.eo), nrow=nrow(M),dimnames=dimnames(M));
twofer.log[]=""

zeo.log=twofer.log; # matrix of character strings, for the ZEO logging

opp.arrow = M; opp.arrow[]=0; # set to 1 if we remove edge due to opposite arrow being found

edge.log = c(); # track the edges we find

# initialize M only after initializing the other matrices above that need to be 0 (zero-ed!!!)
M[snpcols,]=af[snpcols,]; # that has the straight up QTL to begin with.

# ======================================================
# Begin edge orienting - RVV check first...
# ======================================================
pm.print(pm,paste("Begin edge orienting...",date()));

### search for unshielded colliders upstream from our 
# sink (first variable/column)
# alteratively, start recursion with rev$ix[1], the best child candidate

#  Without SNP information, follow FTC in terms of child priority but add V-structure check.
#list the edges into our sink: this will be the column of our sink var

pm.print(pm,paste(sep="","Length of snpcols is: ",length(snpcols)));

if (length(snpcols) < 1) {
   next.edge = -1; # track priority of edge assignment, and FTC => edge < 0, RVV => edge >0
   for (sink  in q$ix) { # later make this a recursive call starting 
                         #with Specified sink var?
      sink.pac = which(af[,sink] != 0); #PArent Candidates for sink
      if (length(sink.pac) < 2) next;
      # parsel out the PACs, two at a time
      pa1pa2=enumpairs(sink.pac)
      for (i in 1:nrow(pa1pa2)) {
           a = pa1pa2[i,1];
           b = pa1pa2[i,2];
           if (Vcdep[a, b, sink] !=0) {
                if(!DFScycle(M,a,sink) && !DFScycle(M,b,sink)) {
                M[a,sink]=next.edge;  M[b,sink]=next.edge; next.edge=next.edge-1;
                # and prevent the other ways, since we have directionality
                M[sink,a]=NaN; M[sink,b]=NaN;
               }
           }
      }
   }
   # fill in the rest from af with undirected links that we didn't resolve:
   for (i in 1:nc) {
       for (j in 1:nc) {
           if (af[i,j] !=0 && M[i,j] != -1 && !is.nan(M[i,j])) { M[i,j]=R[j]; }
   }}
} # end if length(snpcols) < 2

# ======================================================
# Begin SNP-based edge orienting: to avoid double counting,
#  only check the direction A->B where A is smaller index than B
#  (arbitrary but avoids double counting). Save the negative of
#   the edge score of A->B as the edge score for B->A.
# ======================================================

# WITH SNP information: compare and contrast RVV to FTC: 
# Do the RVV and use the SNPs in check22() to establish directionality

if (length(snpcols) >= 1) {

########################################
# begin new Two-step computation: where we maximimize Z_A->B and
# seperately maximum Z_B->A and then make sure the markers are
# independent.
########################################


########################################
# begin new Two-fer SNP walking - for EO building blocks...
########################################

# eo.type is one of c("max.vs.max","forward.step.selection","all.snp.permutations.equal.votes")


mylist.maxmax = walk.two.steps.max.max(pm,snpcols,non.snp.cols,af,no.obs,cx,covx,twofer.eo,n.twofer.eo,cn,twofer.log);
mylist.mlreg = walk.two.steps.weighted.mean.mlreg(datCombined,pm,snpcols,non.snp.cols,af,no.obs,cx,covx,twofer.eo,n.twofer.eo,cn,twofer.log);
mylist.allvote = walk.two.steps.one.permutation.one.vote(datCombined,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,n.twofer.eo,cn,twofer.log);

r$M1M2.AVG = NA # default
r$M1M2.AVG.log = NA

if (pm$do.m1m2.average) {
   mylist.avg.over.all.m1m2 = m1m2.average.losem.walk.two.steps(datCombined,pm,snpcols, non.snp.cols, af, no.obs, cx, covx, twofer.eo, cn, twofer.log)
   r$M1M2.AVG     = mylist.avg.over.all.m1m2$M1M2.AVG
   r$M1M2.AVG.log = mylist.avg.over.all.m1m2$M1M2.AVG.log

  # and save the individual scores
  r$m1m2.scores.AB = mylist.avg.over.all.m1m2$m1m2.scores.AB
  r$m1m2.scores.BA = mylist.avg.over.all.m1m2$m1m2.scores.BA
}

# save stuff so we might recover later if fail after this
r$zeo.mylist.maxmax = mylist.maxmax
r$zeo.mylist.mlreg = mylist.mlreg
r$zeo.mylist.allvote = mylist.allvote

   logfile = paste(sep="",pm$neo.log.file,".r.after.zeo.walking.rdat");
   if (!pm$no.log) { 
         save(r,file=logfile);
         pm.print(pm,paste("Wrote ZEO scores to r$zeo.mylist.* in file: ",logfile));
   }

if (skip.LEO) {
   # make up NA lists
   print("skip.LEO is TRUE. Cirumventing losem.* calls()");
   empty.matrix=matrix(NA,nrow=nrow(af),ncol=ncol(af));
   losem.mylist.maxmax = list(empty.matrix,empty.matrix);
   losem.mylist.mlreg = list(empty.matrix,empty.matrix);
   losem.mylist.allvote = list(empty.matrix,empty.matrix);
} else {
   losem.mylist.maxmax = losem.walk.two.steps.max.max(datCombined,pm,snpcols,non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log,zeo.proxy=mylist.maxmax[[1]]);
   r$leo.mylist.maxmax = losem.mylist.maxmax;
   if (!pm$no.log) save(r,file=paste(sep="",pm$neo.log.file,".r.after.leo.max.rdat"))

   losem.mylist.mlreg = losem.walk.two.steps.weighted.mean.mlreg(datCombined,pm,snpcols,non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log,zeo.proxy=mylist.mlreg[[1]]);
   r$leo.mylist.mlreg = losem.mylist.mlreg;
   if (!pm$no.log) save(r,file=paste(sep="",pm$neo.log.file,".r.after.leo.for.rdat"))


   r$leo.mylist.multione.mlreg = multione.losem.walk.two.steps.weighted.mean.mlreg(x=datCombined,pm=pm,snpcols=snpcols,non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log,zeo.proxy=mylist.mlreg[[1]]);
   if (!pm$no.log) save(r,file=paste(sep="",pm$neo.log.file,".r.after.leo.mulitone.rdat"))


   losem.mylist.allvote = losem.walk.two.steps.one.permutation.one.vote(datCombined,pm,snpcols, non.snp.cols,af,no.obs,cx,covx,twofer.eo,cn,twofer.log,zeo.proxy=mylist.allvote[[1]]);
   r$leo.mylist.allvote = losem.mylist.allvote;
   if (!pm$no.log) save(r,file=paste(sep="",pm$neo.log.file,".r.after.leo.for.max.all.rdat"))

   logfile = paste(sep="",pm$neo.log.file,".r.after.leo.walking.rdat");
   if (!pm$no.log) save(r,file=logfile);
   pm.print(pm,paste("Wrote LEO scores to r$leo.mylist.* in file: ",logfile));
}

} # end if length(snpcols) >1
   
 # pack r
 r$skip.LEO=skip.LEO;
 r$fit.full.model=fit.full.model;
 r$repos.nodes=repos.nodes;
 r$M = M;
 r$eo.type=eo.type;
 r$opp.arrow = opp.arrow
 r$z=z

r$zeo=zeo
r$nzeo=nzeo
r$lcd.edge=lcd.edge
r$lcd.drop=lcd.drop 
r$lcd.possibly.erase=lcd.possibly.erase
r$twofer.eo=twofer.eo
r$n.twofer.eo=n.twofer.eo
r$twofer.log=twofer.log
r$zeo.log=zeo.log
r$edge.log=edge.log




 class(r) = c("NEO","NEO.after.neo.edge.orient");
 r

} # end neo.edge.orient


if(exists("summary.NEO") ) rm(summary.NEO);
summary.NEO=function(r) {

  # unpack r

   if (!inherits(r,"NEO.after.neo.edge.orient")) {
      stop("summary.NEO() stopping because r (1st argument) must inherit from class NEO.after.neo.edge.orient");
   }

    # unpack r
    opp.arrow = r$opp.arrow
    z = r$z
    skip.LEO = r$skip.LEO
    fit.full.model = r$fit.full.model
    repos.nodes = r$repos.nodes
    M = r$M
    eo.type = r$eo.type

    zeo=r$zeo
    nzeo=r$nzeo
    lcd.edge=r$lcd.edge
    lcd.drop =r$lcd.drop 
    lcd.possibly.erase=r$lcd.possibly.erase
    twofer.eo=r$twofer.eo
    n.twofer.eo=r$n.twofer.eo
    twofer.log=r$twofer.log
    zeo.log=r$zeo.log
    edge.log=r$edge.log
   
    datCombined=r$datCombined
    nc=r$nc
    cn=r$cn
    no.obs=r$no.obs
    snpcols=r$snpcols

    non.snp.cols=r$non.snp.cols
    numsnps = r$numsnps

    num.non.snps = r$num.non.snps;
    pm=r$pm;

    super=r$super;
    un.super=r$un.super;
    cx = r$cx;
    af=r$af;

    exog = r$exog;
    lomega = r$lomega;
    omega = r$omega;
    lcx = r$lcx;
    lpx = r$lpx;
    px = r$px;
    R.sq = r$R.sq;
    pcor.filtered = r$pcor.filtered;
    min.k.abs.pcor = r$min.k.abs.pcor;
    ljomega = r$ljomega;
    jomega = r$jomega;
    min.k = r$min.k;
    min.k.abs.pcor.ignore.supernode = r$min.k.abs.pcor.ignore.supernode;
    min.k.nosup = r$min.k.nosup;
    covx = r$covx;
#    neo.log.file = r$neo.log.file

    mylist.maxmax = r$zeo.mylist.maxmax
    mylist.mlreg = r$zeo.mylist.mlreg
    mylist.allvote = r$zeo.mylist.allvote

    losem.mylist.maxmax = r$leo.mylist.maxmax;
    losem.mylist.mlreg = r$leo.mylist.mlreg;
    losem.mylist.allvote = r$leo.mylist.allvote;

    # done unpacking r

   options(warn=1); # print warnings as the occur.

   if (length(snpcols) > 0) {
     
   # choose and log which orienting we are using
   eo.choice=c("forward.step.selection","all.snp.permutations.equal.votes","max.vs.max")
   eo.type=eo.choice[pmatch(pm$eo.type,eo.choice)];

 
  
# check for the consistency of the Edge Orienting techniques, and report any differences in sign on an edge.
# if so, issue a: minority report.

if (pm$do.minority.report.max.for.all) {
for (i in 1:(nrow(M)-1)) { for (j in (i+1):nrow(M)) {
          eo.max.max = mylist.maxmax[[1]][i,j]
          eo.for.step = mylist.mlreg[[1]][i,j]
          eo.allvote = mylist.allvote[[1]][i,j]

          if (is.na(eo.max.max) && is.na(eo.for.step) && is.na(eo.allvote)) next;

          if (is.na(eo.max.max)) eo.max.max=0;
          if (is.na(eo.for.step)) eo.for.step=0;
          if (is.na(eo.allvote)) eo.allvote=0;

          if (eo.max.max > 0 && eo.for.step > 0 && eo.allvote > 0) next;
          if (eo.max.max < 0 && eo.for.step < 0 && eo.allvote < 0) next;
          if (eo.max.max == 0 && eo.for.step == 0 && eo.allvote == 0) next;
          msg=paste(">>>> See ZEO MINORITY REPORT on edge",cn[i],"->",cn[j]," b/c scores were (max,for,all):",signif(eo.max.max,3),signif(eo.for.step,3),signif(eo.allvote,3),"  ;in model:",z$my.title);
          print(" ");
          print(msg);
          print(" ");
          warning(msg);
}}

## repeat this Minority Reporting for the LOSEM statistics

if (!skip.LEO) {
for (i in 1:(nrow(M)-1)) { for (j in (i+1):nrow(M)) {
          eo.max.max = losem.mylist.maxmax$twofer.eo[i,j]
          eo.for.step = losem.mylist.mlreg$twofer.eo[i,j]
          eo.allvote = losem.mylist.allvote$twofer.eo[i,j]

          if (is.na(eo.max.max) && is.na(eo.for.step) && is.na(eo.allvote)) next;

          if (is.na(eo.max.max)) eo.max.max=0;
          if (is.na(eo.for.step)) eo.for.step=0;
          if (is.na(eo.allvote)) eo.allvote=0;

          if (eo.max.max > 0 && eo.for.step > 0 && eo.allvote > 0) next;
          if (eo.max.max < 0 && eo.for.step < 0 && eo.allvote < 0) next;
          if (eo.max.max == 0 && eo.for.step == 0 && eo.allvote == 0) next;

          msg=paste(">>>> See LOSEM MINORITY REPORT on edge",cn[i],"->",cn[j]," b/c scores were (max,for,all):",signif(eo.max.max,3),signif(eo.for.step,3),signif(eo.allvote,3),"  ;in model:",z$my.title);
          pm.print(pm," ");
          pm.print(pm,msg);
          pm.print(pm," ");
          if (!pm$quiet) warning(msg);
}}
}
} # end if (pm$do.minority.report.max.for.all)

## Report any difference of opinion between the LOSEM and the ZEO

if (!skip.LEO & pm$warn.if.zeo.and.leo.disagree) {
for (i in 1:(nrow(M)-1)) { for (j in (i+1):nrow(M)) {
          zeo.eo.max.max = mylist.maxmax[[1]][i,j]
          zeo.eo.for.step = mylist.mlreg[[1]][i,j]
          zeo.eo.allvote = mylist.allvote[[1]][i,j]

          losem.eo.max.max = losem.mylist.maxmax$twofer.eo[i,j]
          losem.eo.for.step = losem.mylist.mlreg$twofer.eo[i,j]
          losem.eo.allvote = losem.mylist.allvote$twofer.eo[i,j]

          if (is.na(zeo.eo.max.max) && is.na(zeo.eo.for.step) && is.na(zeo.eo.allvote)) next;

          if (is.na(zeo.eo.max.max)) zeo.eo.max.max=0;
          if (is.na(zeo.eo.for.step)) zeo.eo.for.step=0;
          if (is.na(zeo.eo.allvote)) zeo.eo.allvote=0;

          if (is.na(losem.eo.max.max) && is.na(losem.eo.for.step) && is.na(losem.eo.allvote)) next;

          if (is.na(losem.eo.max.max)) losem.eo.max.max=0;
          if (is.na(losem.eo.for.step)) losem.eo.for.step=0;
          if (is.na(losem.eo.allvote)) losem.eo.allvote=0;

          if (zeo.eo.max.max > 0 && zeo.eo.for.step > 0 && zeo.eo.allvote > 0 && 
               losem.eo.max.max > 0 && losem.eo.for.step > 0 && losem.eo.allvote > 0) next;
          if (zeo.eo.max.max < 0 && zeo.eo.for.step < 0 && zeo.eo.allvote < 0 && 
               losem.eo.max.max < 0 && losem.eo.for.step < 0 && losem.eo.allvote < 0) next;
          if (zeo.eo.max.max == 0 && zeo.eo.for.step == 0 && zeo.eo.allvote == 0 && 
               losem.eo.max.max == 0 && losem.eo.for.step == 0 && losem.eo.allvote == 0) next;

          msg=paste(">>>> LEO and ZEO disagree. On edge",cn[i],"->",cn[j]," b/c scores were ZEO(max,for,all):",signif(zeo.eo.max.max,3),signif(zeo.eo.for.step,3),signif(zeo.eo.allvote,3),"  ;in model:",z$my.title,"  LOSEM(max,for,all):",signif(losem.eo.max.max,3),signif(losem.eo.for.step,3),signif(losem.eo.allvote,3));
          pm.print(pm," ");
          pm.print(pm,msg);
          pm.print(pm," ");
          if (!pm$quiet) warning(msg);
}}
}
 

mylist=mylist.maxmax; # default

losem.mylist=losem.mylist.maxmax; #default

if (eo.type == "max.vs.max") {
     mylist= mylist.maxmax;
     losem.mylist= losem.mylist.maxmax;
} else if (eo.type == "forward.step.selection") {
     mylist = mylist.mlreg;
     losem.mylist = losem.mylist.mlreg;
} else if (eo.type == "all.snp.permutations.equal.votes") {
     mylist = mylist.allvote;
     losem.mylist = losem.mylist.allvote;
} else warning(paste("unrecognized eo.type argument: ",eo.type,"\nmust be one \nof c(\"max.vs.max\",\"forward.step.selection\",\"all.snp.permutations.equal.votes\")\nDefaulting to max.vs.max."))

pm.print(pm,"--------------------------------------------------");
pm.print(pm,paste("EDGE Orienting using criteria: ",eo.type));
pm.print(pm,"--------------------------------------------------");

twofer.eo = mylist[[1]]
n.twofer.eo = mylist[[2]]
twofer.log = mylist[[3]]


#if (!skip.LEO) losem.twofer.eo = losem.mylist$twofer.eo # LEO.O
 if (!skip.LEO) {
   # this is the scores we use for the graph.
   losem.twofer.eo = losem.mylist$AB.vs.NextBest # LEO.NB, while FOR,MAX,ALL depends on pm$eo.type
   # and fill in NA values with 0's ... so we can neo.graph() okay.
   na.mat = is.na(losem.twofer.eo)
   losem.twofer.eo[na.mat]=0
    }

#sink.cand = setdiff(q$ix,snpcols); #SNPs can't be sinks && still use the child priority suggested by q$ix 
# Modification by Peter Langfelder as suggested by JEA: q is not defined, so treat q$ix as empty
sink.cand = NULL;
for (sink  in sink.cand) { 
   sink.pac = which(af[,sink] !=0); # PArent Candidates for sink

   snp.pac = snpcols[snpcols %in% sink.pac] # list parents of sink that are SNPs

   non.snp.pac = setdiff(sink.pac,snp.pac);

   # do the LCD checks
   for (A in non.snp.pac) {    if (!pm$no.log) save(z,file=z$post.save.fname);
         # setdiff so we exclude sink as grandparent of itself.
         gpa = setdiff(which(af[,A] !=0),sink);

             for (j in gpa) {
                  if (j %in% snpcols) {
                             ab=checkLCD(covx, cx, A, sink, j, pm);
                             if (ab$edge) {
                                 msg = paste("found LCD directional edge from", cn[A], "to", cn[sink],"supported by SNP:",cn[j]);
                                 pm.print(pm,msg);
                                 edge.log=rbind(edge.log,msg);
                                 M[A,sink] = 1;
                                 opp.arrow[sink,A]=1;
                                 lcd.edge[A,sink] = j;
                                 lcd.drop[A,sink] = ab$drop; 
                             }
                      }
                 } #end for j
   } # end for A in non.snp.pac / end LCD checks
   
   
   if (length(sink.pac) < 2) next;
   
   # do the RVV checks on the current sink
   # parsel out the PACs, two at a time
   pa1pa2=enumpairs(sink.pac)
   for (i in 1:nrow(pa1pa2)) {
          pa1 = pa1pa2[i,1];
          pa2 = pa1pa2[i,2];
          cv=checkV(covx, cx, sink, pa1, pa2, pm);
          if (cv$v == FALSE) {
              #pm.print(pm,paste("diag: checkV answered NO on sink",cn[sink],"and parents",cn[pa1],cn[pa2],"details:",cv));
              next; }

          gpa1 = setdiff(which(af[,pa1] !=0),sink);
          if (length(gpa1) >= 2) {
              gg=enumpairs(gpa1);
              for (j in 1:nrow(gg)) {
                 gpa1.1 = gg[j,1];
                 gpa1.2 = gg[j,2];

                 # avoid double counting: only do A->B if A<B
                 if (pa1 > sink) next;

                 cvv = checkRVV(covx, cx, sink, pa2, pa1, gpa1.1, gpa1.2, pm,cn);
                 if (cvv$rvv) {
                     M[pa1,sink]=1;
                     opp.arrow[sink,pa1]=1;
                     msg = paste("found RVV directional edge from", cn[pa1], "to", cn[sink],"supported by Grandparents of",cn[pa1],":",cn[gpa1.1],cn[gpa1.2],"and parent:",cn[pa2]);
                     pm.print(pm,msg);
                     edge.log=rbind(edge.log,msg);
                 }
              }
          }

          gpa2 = setdiff(which(af[,pa2] !=0),sink);
          if (length(gpa2) < 2) next;

          gg=enumpairs(gpa2);
          for (j in 1:nrow(gg)) {
               gpa2.1 = gg[j,1];
               gpa2.2 = gg[j,2];

               # avoid double counting: only do A->B if A<B
               if (pa2 > sink) next;

               cvv = checkRVV(covx, cx, sink, pa1, pa2, gpa2.1, gpa2.2, pm,cn);
               if (cvv$rvv) {
                   M[pa2,sink]=1;
                   opp.arrow[sink,pa2]=1;
                   msg = paste("found RVV directional edge from", cn[pa2], "to", cn[sink],"supported by Grandparents of",cn[pa2],":",cn[gpa2.1],cn[gpa2.2],"and parent:",cn[pa1]);
                   pm.print(pm,msg);
                   edge.log=rbind(edge.log,msg);
               }
          }
    } # end i in 1:nrow(pa1pa2) / RVV checks


   # do the 2 SNP check, at this point this should certainly be redundant!
   if (length(snp.pac) < 1) next;

   # parsel out the PACs, two at a time
   pa1pa2=enumpairs(sink.pac)
   for (i in 1:nrow(pa1pa2)) {
         B=0;
         if (pa1pa2[i,1] %in% snpcols  && !(pa1pa2[i,2] %in% snpcols)) {
                 m1 = pa1pa2[i,1];
                 B = pa1pa2[i,2];
         }

         if (pa1pa2[i,2] %in% snpcols && !(pa1pa2[i,1] %in% snpcols)) {
                 m1 = pa1pa2[i,2];
                 B = pa1pa2[i,1];
         }

         # avoid double counting: only do A->B if A<B
         if (sink >= B) next;

         if (B != 0) {
                 # setdiff so we exclude sink as grandparent of itself.
                 gpa = setdiff(which(af[,B] !=0),c(sink,m1)); 
                 for (j in gpa) {
                      if (j %in% snpcols) {
                             ab=check22(covx, cx, sink, B, m1, j, pm, no.obs,zeo.log);
                             zeo.log=ab$zeo.log; # update log after the call.
                             if (ab$zeo > 0) { 
                                 zeo[sink,B] = zeo[sink,B] + ab$zeo; nzeo[sink,B]=nzeo[sink,B]+4;
                                 zeo[B,sink] = zeo[B,sink] - ab$zeo; nzeo[B,sink]=nzeo[B,sink]+4;

                                 msg = paste("found CCO directional edge from", cn[sink], "to", cn[B],"supported by SNPs:",cn[j],"and",cn[m1]);
                                 pm.print(pm,msg);
                                 edge.log=rbind(edge.log,msg);
                                 M[sink,B]=1; 
                                 #af[B,sink]=0; # comment out: leave af as is! add stuff to M instead.
                                 opp.arrow[B,sink]=1;} 
                             else if (ab$zeo < 0) { 
                                 # CCO == Controlled Competitive Orienting of edges.
                                 zeo[B,sink] = zeo[B,sink] - ab$zeo; nzeo[B,sink]=nzeo[B,sink]+4; # zeo is already opposite, so flip the sign
                                 zeo[sink,B] = zeo[sink,B] + ab$zeo; nzeo[sink,B]=nzeo[sink,B]+4;

                                 msg = paste("found CCO directional edge from", cn[B], "to", cn[sink],"supported by SNPs:",cn[j],"and",cn[m1]);
             pm.print(pm,msg);
                                 edge.log=rbind(edge.log,msg);
                                 M[B,sink]=1; 
                                 #af[sink,B]=0; # comment out: leave af is is! ADD stuff to M instead.
                                 opp.arrow[sink,B]=1; } else {
                                       pm.print(pm,paste("ab$zeo==0: no ZEO edge found among (A,B,m1,m2):",cn[sink],cn[B],cn[m1],cn[j],"\n")); pm.print(pm,ab$why.not)}
                      }
                 } #end for j
         } #end if B != 0
   } # end for i
  } # end for sink

#### Force the primary M result to be due to two-step scores...

#for those edges that have twostep.edge orienting done on them, use the positive direction...
final.eo = twofer.eo 
if (!skip.LEO && !is.null(pm$use.leo) && pm$use.leo==TRUE) { 
   final.eo = losem.twofer.eo; 
   pm.print(pm,"--------------------------------------------------");
   pm.print(pm,paste("EDGE Orienting using Local SEM (LEO score) models. "));
   pm.print(pm,"--------------------------------------------------");
 } else {
   pm.print(pm,"--------------------------------------------------");
   pm.print(pm,paste("EDGE Orienting using ZEO. "));
   pm.print(pm,"--------------------------------------------------");
}

#pm.print(pm,"M before twostep orienting:");
#pm.print(pm,M)

# modified to work for non-symmetric final.eo scores...
for (i in 1:(nrow(M)-1)) { for (j in (i+1):nrow(M)) {
      if (!is.nan(final.eo[i,j]) && !is.na(final.eo[i,j])) { 
          if (final.eo[i,j] > pm$graph.arrow.edge.th & final.eo[j,i] < pm$graph.arrow.edge.th) { 
             M[i,j]=1; M[j,i]=0; opp.arrow[j,i]=1; 
          } else if (final.eo[j,i] > pm$graph.arrow.edge.th & final.eo[i,j] < pm$graph.arrow.edge.th) { 
           M[j,i]=1; M[i,j]=0; opp.arrow[i,j]=1;
          } else if (final.eo[i,j] != 0) {
             # edge, but indeterminate....
             M[i,j]=M[j,i]=1;
          }
      }
}}
#pm.print(pm,"M after twostep orienting:");
#pm.print(pm,M)

     # end if length(snpcols) > 0
 } else {
    pm.print(pm,"No SNPs: current configuration does no final edge orienting.");
 }
   
pm.print(pm,paste("Done with edge orienting...",date()));

pm.print(pm,paste("Writing logs and computing final statistics...",date()));

if (pm$do.final.colliders) {
# ========================================
# Compute and save just the two-step score for all found colliders
# ========================================
collider.name=c();
marg.cor=c();
cond.cor=c();
Z.marg=c();
Z.cond=c();
Z.cond.minus.Z.marg =c(); # Z.parents.correlation - Z.parents.given.child.correlation
child=c()

for (j in 1:(ncol(M))) { 
# for each column j, evaluate each pair of rows that has a 1 (collider at j)
      for (i in 1:(nrow(M)-1)) {
      for (k in (i+1):nrow(M)) {
          if (M[i,j] && M[k,j]) { 
              child=c(child,j);
              collider.name=c(collider.name,paste(cn[i],"->",cn[j],"<-",cn[k]));

              rho= cx[i,k];
              marg.cor=c(marg.cor, rho);

              z2= NA # init to try-failed defaults.
              rho.given.j = NA; # init to try-failed defaults. avoid crashing pcor with singular values.
              if (rho < 1 ) {
                  pcor.tmp = try(pcor(c(i,k,j),covx));
                  if (!inherits(pcor.tmp, "try-error")) {
                     rho.given.j= pcor.tmp;
                     z2=zeo2(i,k,j,covx,pm);
                  }
               }
              cond.cor=c(cond.cor, rho.given.j);
              
              Zpa=sqrt(pm$no.obs.Z-pm$fisher.dof.cor)*   fisher1(abs(rho));
              Z.marg=c(Z.marg,Zpa);

              Zpgc=sqrt(pm$no.obs.Z-pm$fisher.dof.pcor1)*   fisher1(abs(rho.given.j));
              Z.cond=c(Z.cond,Zpgc);

              if (!is.na(z2[1])) {
                   Z.cond.minus.Z.marg=c(Z.cond.minus.Z.marg,-z2$eo);
                 } else {
                   Z.cond.minus.Z.marg=c(Z.cond.minus.Z.marg,NA);
                 }
          }
      }
}}

   
if (length(marg.cor) > 0) {
   marg.cor=signif(marg.cor,3);
   cond.cor=signif(cond.cor,3);
   Z.marg=signif(Z.marg,3);
   Z.cond=signif(Z.cond,3);
   Z.cond.minus.Z.marg=signif(Z.cond.minus.Z.marg,3);
} else {
   warning("no edges oriented by ZEO");
   pm.print(pm,"No edges oriented by ZEO");
}

collider.zz = data.frame(child,collider.name, marg.cor, cond.cor, Z.marg, Z.cond, Z.cond.minus.Z.marg);
} else { collider.zz = NA; pm.print(pm,"Skipping final colliders since pm$do.final.colliders was FALSE.");  } # end do.final.colliders   


# compute our direction only and full model
dir.only = M;

# OLD: M = ((af - opp.arrow) | un.super)+0;
# NEW: orient using two-step modifications to M above...
M = (M | un.super)+0;

# =================================================
# Post processing: present the data nicely.
# =================================================

# create an index to linearize the matrices and turn them into the
# dataframe which summarizes each edge.
linearizer=matrix(1:(nc*nc),nrow=nc);
qlin=linearizer[as.logical(exog)]; # use exog to discard self-edges and edges into SNPs
 # details for the spr.nd column
super2=list();

# create a dataframe with summary for each possible edge
edgelist=list();
edgelist[[nc*nc]]=NA; # pre-allocate the list
super2=edgelist # pre-allocate
edgelist.plus.num=edgelist; # version of edgelist with (e.g. [3,6]) numeric indices, for ease finding log entries.

k=1
for (j in 1:nc) { for (i in 1:nc) { 
    edgelist[[k]]=paste(cn[i],"->",cn[j]) 
    edgelist.plus.num[[k]]=paste(cn[i],"[",i,"] ->",cn[j]," [",j,"]")
    if (super[i]==super[j]) { super2[[k]]="[*]"; } else { super2[[k]]="[ ]"} 
    k=k+1
}}

#edgelist[qlin];

# we allowed LCD to generate arcs possibly in both directions 
# if LCD found an arc in only one direction, then only it will be in M.
# There is no further need to clean up 


stringify.qlin=function(cx,lcx,qlin) {
     st=lapply(lcx[qlin]+0,mapstar);
     substr(paste(sep="","[",st,"] ",cx[qlin]),1,8) } # 8 = number of characters to keep in a column.

edge= paste(edgelist[qlin]);
edge.plus.num = paste(edgelist.plus.num[qlin]);
keep=M[qlin];
cohen.omega= stringify.qlin((1-omega),lomega,qlin);
#abs.cor = cx[qlin];
abs.cor = stringify.qlin(cx,lcx,qlin);
abs.pcor = stringify.qlin(px,lpx,qlin);
# not computed anymore: Rsq = R.sq[qlin];
k.min.abs.pcor = stringify.qlin(min.k.abs.pcor, pcor.filtered,qlin);


jomeg = stringify.qlin(jomega, ljomega , qlin);
k.min = cn[min.k[qlin]];

k.min.nosup.pcor = stringify.qlin(min.k.abs.pcor.ignore.supernode,(min.k.abs.pcor.ignore.supernode > pm$k.pcor.th),qlin);
k.min.nosup = cn[min.k.nosup[qlin]];

Rm.op = paste(sep="","[",lapply(1-opp.arrow[qlin],mapstar),"]");
a.flt = af[qlin];
Spr=as.character(super2)[qlin];
ZEO = (zeo / sqrt(nzeo))[qlin];
n.z = nzeo[qlin]

step2.eo = final.eo[qlin]
step2.eo = signif(step2.eo,3);
n.2 = n.twofer.eo[qlin]

# rank the most significant edges:
step2.eo.rank=sort(step2.eo,decreasing=TRUE,index.return=TRUE)$ix


lcd=substr(paste(sep="","[",cn[lcd.edge[qlin]],"] ",lcd.drop[qlin]),1,12);

#options(digits = 4);

zeo.max=mylist.maxmax[[1]][qlin]
zeo.for=mylist.mlreg[[1]][qlin]
zeo.all=mylist.allvote[[1]][qlin]

 # keep just 3 sig figures
sig.fig.keep=3

# make it easy to pair up the best single marker RMSEAs for an edge, by transposing
BestAB.BA = t(losem.mylist.maxmax$BestAB.single.marker.RMSEA)
BestAB.BA.marker = t(losem.mylist.maxmax$BestAB.single.marker.RMSEA.the.marker)

SIMPLE.MAX.MAX=convert.na.to.blank(pm,signif(losem.mylist.maxmax$Simple.Max.vs.Max[qlin],sig.fig.keep))
   
# 
BEST.1M.RMSEA.AB  = signif(losem.mylist.maxmax$BestAB.single.marker.RMSEA[qlin],sig.fig.keep)
BEST.1M.RMSEA.BA  = signif(BestAB.BA[qlin],sig.fig.keep)
BEST.1M.MARKER.AB = losem.mylist.maxmax$BestAB.single.marker.RMSEA.the.marker[qlin]
BEST.1M.MARKER.BA = BestAB.BA.marker[qlin]


LEO.O.MAX=losem.mylist.maxmax$twofer.eo[qlin]
LEO.O.FOR=losem.mylist.mlreg$twofer.eo[qlin]
LEO.O.ALL=losem.mylist.allvote$twofer.eo[qlin]

LEO.I.FOR = signif(losem.mylist.mlreg$LEO.I[qlin],sig.fig.keep)
LEO.I.MAX = signif(losem.mylist.maxmax$LEO.I[qlin],sig.fig.keep)
LEO.I.ALL = signif(losem.mylist.allvote$LEO.I[qlin],sig.fig.keep)

 RMSEA.2m=losem.mylist.mlreg$RMSEA.2m[qlin]

 # LEO.MULTIONE.NB is now called LEO.NB.CPA too.
 LEO.MULTIONE.NB = convert.na.to.blank(pm,signif(r$leo.mylist.multione.mlreg$AB.vs.NextBest[qlin],sig.fig.keep))

if (pm$do.m1m2.average) {
   M1M2.AVG =  convert.na.to.blank(pm,signif(r$M1M2.AVG[qlin],sig.fig.keep))
   # already there: M1M2.AVG.log = r$M1M2.AVG.log
} else {
   M1M2.AVG = rep(" ",length(qlin))
   r$M1M2.AVG.log = rep("",length(qlin))
}



 # to sort by LEO.I, descending, use just before making the final data frames.
# sort.by=LEO.I.FOR
#sort.by =-RMSEA.2m

sort.by=losem.mylist.mlreg$AB.vs.NextBest[qlin] # sort by LEO.NB.OCA, but not called LEO.NB.OCA yet b/c not in final sorted order.


 # default order is no order...
 priority.sort=1:length(sort.by)  

 if (any(!is.na(sort.by))) {
    # Now compute priority.sort to order the rows of the spreadsheet csv/excel file.
    # We do special handling to get the NA to sort to the end AND get an ix sort index back....
    #
    w.is.na=which(is.na(sort.by))
    min.tmp = min(sort.by,na.rm=TRUE)
    num.na=length(w.is.na)
    if (num.na) {
      sort.by[w.is.na]=(min.tmp-1):(min.tmp-num.na)
    }
    priority.sort=sort(sort.by,decreasing=TRUE,index.return=TRUE)$ix
 } else {
   pm.print(pm,"Hmmm: all ordering (sort.by) scores were NA! Spreadsheet will be unsorted.");
 }

 # after determining the priority sort, then get rid of NAs
 LEO.I.FOR = convert.na.to.blank(pm,LEO.I.FOR)
 LEO.I.MAX = convert.na.to.blank(pm,LEO.I.MAX)
 LEO.I.ALL = convert.na.to.blank(pm,LEO.I.ALL)

 RMSEA.2m = convert.na.to.blank(pm,signif(RMSEA.2m,sig.fig.keep))

 BLV = convert.na.to.blank(pm,signif(losem.mylist.mlreg$BLV[qlin],sig.fig.keep))
 ZPathAB = convert.na.to.blank(pm,signif(losem.mylist.mlreg$ZPathAB[qlin],sig.fig.keep))

 LEO.NB.FOR = convert.na.to.blank(pm,signif(losem.mylist.mlreg$AB.vs.NextBest[qlin],sig.fig.keep))
 LEO.NB.MAX = convert.na.to.blank(pm,signif(losem.mylist.maxmax$AB.vs.NextBest[qlin],sig.fig.keep))
 LEO.NB.ALL = convert.na.to.blank(pm,signif(losem.mylist.allvote$AB.vs.NextBest[qlin],sig.fig.keep))

 Final.SNPs.LEO.NB.OCA = convert.na.to.blank(pm,losem.mylist.mlreg$Final.SNPs.LEO.NB.OCA[qlin])
 Minus.log.Model.P.value.AtoB = convert.na.to.blank(pm,signif(losem.mylist.mlreg$Minus.log.Model.P.value.AtoB[qlin],sig.fig.keep))
 Model.P.value.AtoB = convert.na.to.blank(pm,signif(10^-(losem.mylist.mlreg$Minus.log.Model.P.value.AtoB[qlin]),sig.fig.keep))
 P.weighted.LEO.NB.OCA = convert.na.to.blank(pm,signif(losem.mylist.mlreg$AB.vs.NextBest[qlin]*(10^-(losem.mylist.mlreg$Minus.log.Model.P.value.AtoB[qlin])),sig.fig.keep))

 PearsonCor = signif(cor(datCombined)[qlin],2)

 # c(BLV,ZPathAB,AB.vs.NextBest,PearsonCor)

 zeo.max=convert.na.to.blank(pm,signif(zeo.max,sig.fig.keep))
 zeo.for=convert.na.to.blank(pm,signif(zeo.for,sig.fig.keep))
 zeo.all=convert.na.to.blank(pm,signif(zeo.all,sig.fig.keep))

 LEO.O.MAX=convert.na.to.blank(pm,signif(LEO.O.MAX,sig.fig.keep))
 LEO.O.FOR=convert.na.to.blank(pm,signif(LEO.O.FOR,sig.fig.keep))
 LEO.O.ALL=convert.na.to.blank(pm,signif(LEO.O.ALL,sig.fig.keep))

# hyperlinks to log entries...
max.zeo.log.link=mylist.maxmax[[3]]
all.zeo.log.link=mylist.allvote[[3]]
for.zeo.log.link=mylist.mlreg[[3]]

if (!skip.LEO) {
max.leo.log.link=losem.mylist.maxmax$twofer.log
all.leo.log.link=losem.mylist.allvote$twofer.log
for.leo.log.link=losem.mylist.mlreg$twofer.log
multione.leo.log.link=r$leo.mylist.multione.mlreg$twofer.log
m1m2.avg.log.link=r$M1M2.AVG.log
}

#
# put the actual score as title1 == what shows up in the spreadsheet column
#   ...but with a * after it to indicate its got a log...eventually replace
#   the non-log entries too?
#
if (!pm$no.log) {
for (i in 1:nc) { for (j in 1:nc) {

###   Already taken care of by the log.link == "" :   
###   if (!is.null(pm$A) & !is.A.B.edge(cn[i],cn[j],pm)) next;

   if (max.zeo.log.link[i,j] != "") { max.zeo.log.link[i,j]=write.html.log(pm,fname=paste(sep="",cn[i],".to.",cn[j],".html"),title.1=paste(sep="",signif(mylist.maxmax[[1]][i,j],sig.fig.keep)," *"),title.2=paste(sep="",cn[i]," -> ",cn[j]), title.3="MAX.ZEO", my.text=max.zeo.log.link[i,j],html.log.dir=paste(sep="",pm$neo.log.file,".zeo.max"),parent.dir=); }
   if (all.zeo.log.link[i,j] != "") { all.zeo.log.link[i,j]=write.html.log(pm,fname=paste(sep="",cn[i],".to.",cn[j],".html"),title.1=paste(sep="",signif(mylist.allvote[[1]][i,j],sig.fig.keep)," *"),title.2=paste(sep="",cn[i]," -> ",cn[j]), title.3="ALL.ZEO", my.text=all.zeo.log.link[i,j],html.log.dir=paste(sep="",pm$neo.log.file,".zeo.all")); }
   if (for.zeo.log.link[i,j] != "") { for.zeo.log.link[i,j]=write.html.log(pm,fname=paste(sep="",cn[i],".to.",cn[j],".html"),title.1=paste(sep="",signif(mylist.mlreg[[1]][i,j],sig.fig.keep)," *"),title.2=paste(sep="",cn[i]," -> ",cn[j]), title.3="FOR.ZEO", my.text=for.zeo.log.link[i,j],html.log.dir=paste(sep="",pm$neo.log.file,".zeo.for")); }

  if (!skip.LEO) {
   if (max.leo.log.link[i,j] != "") { max.leo.log.link[i,j]=write.html.log(pm,fname=paste(sep="",cn[i],".to.",cn[j],".html"),title.1=paste(sep="",signif(losem.mylist.maxmax$twofer.eo[i,j],sig.fig.keep)," *"),title.2=paste(sep="",cn[i]," -> ",cn[j]), title.3="MAX.LEO",my.text=max.leo.log.link[i,j],html.log.dir=paste(sep="",pm$neo.log.file,".leo.max")); }
   if (all.leo.log.link[i,j] != "") { all.leo.log.link[i,j]=write.html.log(pm,fname=paste(sep="",cn[i],".to.",cn[j],".html"),title.1=paste(sep="",signif(losem.mylist.allvote$twofer.eo[i,j],sig.fig.keep)," *"),title.2=paste(sep="",cn[i]," -> ",cn[j]), title.3="ALL.LEO",my.text=all.leo.log.link[i,j],html.log.dir=paste(sep="",pm$neo.log.file,".leo.all")); }

### old, FOR logging:   if (for.leo.log.link[i,j] != "") { for.leo.log.link[i,j]=write.html.log(pm,fname=paste(sep="",cn[i],".to.",cn[j],".html"),title.1=paste(sep="",signif(losem.mylist.mlreg$twofer.eo[i,j],sig.fig.keep)," *"),title.2=paste(sep="",cn[i]," -> ",cn[j]), title.3="FOR.LEO",my.text=for.leo.log.link[i,j],html.log.dir=paste(sep="",pm$neo.log.file,".leo.for")); }
### new, OCA specific logging:
   if (for.leo.log.link[i,j] != "") { for.leo.log.link[i,j]=write.html.log(pm,fname=paste(sep="",cn[i],".to.",cn[j],".html"),title.1=paste(sep="",signif(losem.mylist.mlreg$AB.vs.NextBest[i,j],sig.fig.keep)," *"),title.2=paste(sep="",cn[i]," -> ",cn[j]), title.3="LEO.NB.OCA",my.text=for.leo.log.link[i,j],html.log.dir=paste(sep="",pm$neo.log.file,".leo.nb.oca")); }

   if (r$leo.mylist.multione.mlreg$twofer.log[i,j] != "") { multione.leo.log.link[i,j]=write.html.log(pm,fname=paste(sep="",cn[i],".to.",cn[j],".html"),title.1=paste(sep="",signif(r$leo.mylist.multione.mlreg$AB.vs.NextBest[i,j],sig.fig.keep)," *"),title.2=paste(sep="",cn[i]," -> ",cn[j]), title.3="MULTIONE.LEO",my.text=multione.leo.log.link[i,j],html.log.dir=paste(sep="",pm$neo.log.file,".leo.cpa")); }
   if (pm$do.m1m2.average) {   
      if (r$M1M2.AVG.log[i,j] != "") { m1m2.avg.log.link[i,j]=write.html.log(pm,fname=paste(sep="",cn[i],".to.",cn[j],".html"),title.1=paste(sep="",signif(r$M1M2.AVG[i,j],sig.fig.keep)," *"),title.2=paste(sep="",cn[i]," -> ",cn[j]), title.3="M1M2.AVG",my.text=m1m2.avg.log.link[i,j],html.log.dir=paste(sep="",pm$neo.log.file,".m1m2.avg")); 
      } 
   }
 }

 }}
} # end if (!pm$no.log)
   

# correct the log entry string for the opposite edge.
if (exists("flip.sign.on.log.string")) rm(flip.sign.on.log.string);
flip.sign.on.log.string=function(s) {
   len = nchar(s,type="c");
   start = regexpr("\",\"",s);
   if (start == -1) return("");
   start=start+3;
   if (substr(s,start,start)=="-") {
      return(paste(sep="",substr(s,1,start-1),substr(s,start+1,len-2)," [op]\")"));
   } else {
      return(paste(sep="",substr(s,1,start-1),"-",substr(s,start,len-2)," [op]\")"));
   }
}

# now fill in logs for the opposite edges....
for (i in 1:nc) { for (j in 1:nc) {
   if (max.zeo.log.link[i,j] == "" && max.zeo.log.link[j,i] != "") { max.zeo.log.link[i,j]=flip.sign.on.log.string(max.zeo.log.link[j,i]); }
   if (all.zeo.log.link[i,j] == "" && all.zeo.log.link[j,i] != "") { all.zeo.log.link[i,j]=flip.sign.on.log.string(all.zeo.log.link[j,i]); }
   if (for.zeo.log.link[i,j] == "" && for.zeo.log.link[j,i] != "") { for.zeo.log.link[i,j]=flip.sign.on.log.string(for.zeo.log.link[j,i]); }

  if (!skip.LEO) {
    if (max.leo.log.link[i,j] == "" && max.leo.log.link[j,i] != "") { max.leo.log.link[i,j]=flip.sign.on.log.string(max.leo.log.link[j,i]); }
    if (all.leo.log.link[i,j] == "" && all.leo.log.link[j,i] != "") { all.leo.log.link[i,j]=flip.sign.on.log.string(all.leo.log.link[j,i]); }
    if (for.leo.log.link[i,j] == "" && for.leo.log.link[j,i] != "") { for.leo.log.link[i,j]=flip.sign.on.log.string(for.leo.log.link[j,i]); }
   }
 }}


# reduce to just the linearized....
max.zeo.log=max.zeo.log.link[qlin];
all.zeo.log=all.zeo.log.link[qlin];
for.zeo.log=for.zeo.log.link[qlin];

if (!skip.LEO) {
   max.leo.log=max.leo.log.link[qlin];
   all.leo.log=all.leo.log.link[qlin];
   for.leo.log=for.leo.log.link[qlin];
   multione.leo.log=multione.leo.log.link[qlin];
   LINK.to.LOG.M1M2.AVG = m1m2.avg.log.link[qlin]

} else {
   max.leo.log=rep(NA,length(qlin));
   all.leo.log=rep(NA,length(qlin));
   for.leo.log=rep(NA,length(qlin));
   multione.leo.log=rep(NA,length(qlin));
   LINK.to.LOG.M1M2.AVG = rep(NA,length(qlin));
}


I1=rep("|zeo>",length(qlin));
I2=rep("<LEO|",length(qlin));
#I3=rep("|>",length(qlin));
I4=rep("<ZEO|",length(qlin));
II=rep("     ",length(qlin)); # blank spacer


# just do one spreadsheet now, just df.for.excel, and make df the same.

abs.pcor.best = k.min.abs.pcor 
best.conditioning.node = k.min
abs.pcor.best.outside.supernode =  k.min.nosup.pcor
best.conditioning.node.outside.supernode = k.min.nosup
relative.change.cor.pcor.best = jomeg
relative.change.cor.pcor = cohen.omega
Together.in.supernode = Spr

LINK.to.LOG.LEO.O.FOR = for.leo.log
LINK.to.LOG.LEO.O.MAX = max.leo.log
LINK.to.LOG.LEO.O.ALL = all.leo.log
LINK.to.LOG.MULTIONE  = multione.leo.log

LINK.to.LOG.OCA = for.leo.log

BLV.or.BilayerZscore = BLV
ZtestPathCoefficientAB = ZPathAB
EdgeSurvives.First.filters = a.flt
EdgeSurvives.Final.filter = keep

# rename to reflect the improved names in Paper.
LEO.NB.CPA = LEO.MULTIONE.NB
LINK.to.LOG.CPA = LINK.to.LOG.MULTIONE
LEO.NB.OCA = LEO.NB.FOR

# took out: , Minus.log.Model.P.value.AtoB # b/c kind of redundant with Model.P.value.AtoB

df.for.excel= data.frame(LEO.NB.OCA, LINK.to.LOG.OCA, Final.SNPs.LEO.NB.OCA, Model.P.value.AtoB, P.weighted.LEO.NB.OCA, II, LEO.NB.CPA, LINK.to.LOG.CPA, II, SIMPLE.MAX.MAX, LEO.NB.MAX, LEO.NB.ALL, II, M1M2.AVG, LINK.to.LOG.M1M2.AVG, II, LEO.I.FOR, LEO.I.MAX, LEO.I.ALL, II, LEO.O.FOR, LEO.O.MAX, LEO.O.ALL, II, PearsonCor, II, LINK.to.LOG.LEO.O.FOR, LINK.to.LOG.LEO.O.MAX, LINK.to.LOG.LEO.O.ALL, II, zeo.for,zeo.max,zeo.all, II, BLV.or.BilayerZscore, ZtestPathCoefficientAB, II, EdgeSurvives.First.filters, EdgeSurvives.Final.filter, II, all.zeo.log,max.zeo.log,for.zeo.log, II, ZEO, abs.cor, abs.pcor , abs.pcor.best, best.conditioning.node, abs.pcor.best.outside.supernode, best.conditioning.node.outside.supernode,  relative.change.cor.pcor , relative.change.cor.pcor.best, Together.in.supernode,RMSEA.2m,BEST.1M.RMSEA.AB,BEST.1M.RMSEA.BA, II, BEST.1M.MARKER.AB,BEST.1M.MARKER.BA)
rownames(df.for.excel)=edge;

df=df.for.excel

# clean up the colname blanks
w.blank = which(substr(colnames(df.for.excel),1,2)=="II")
num.underscores=5+(1:length(w.blank))
for (i in 1:length(w.blank)) { 
    colnames(df.for.excel)[w.blank[i]] = paste(rep("_",num.underscores[i]),collapse=""); 
}


df.for.excel=df.for.excel[priority.sort,]

# might not always work, but we want the date and run title
# to be on the csv file to for sane viewing/keeping track of things.
excel.file=get.excel.file(pm) # path.os.convert(paste(sep="",dirname(pm$neo.log.file),"\\",basename(dirname(pm$neo.log.file)),".csv"))
z$excel.file=excel.file

if (!pm$doing.permutation) { # don't overwrite our original excel results with the permutations!!!
   if (pm$first.batch.pass) {
    # col.names=NA to get the space over row.names in the spreadsheet.

    if (!pm$no.log) write.table(df.for.excel, file=excel.file, col.names = NA, sep=",",dec=".",qmethod="double" )
   } else {
    if (!pm$no.log) write.table(df.for.excel, file=excel.file,append=TRUE, col.names = FALSE, sep=",",dec=".",qmethod="double")
   }
}

if (pm$doing.permutation) {
    # note if doing permutations, we must save the spreadsheet...no other record of the data is kept.

   # report the TO and FROM columns for ease of sorting later, if someone wants to excel it.
   splt=strsplit(edge,split=" -> ",fixed=TRUE)
   FROM.A=sapply(splt,function(x) x[1])
   TO.B  =sapply(splt,function(x) x[2])

   # construct different dataframe, becuase the links to logs were mistakenly containing full actually logs! (Ugh!)
    df.for.perm= data.frame(FROM.A,TO.B,LEO.NB.OCA, LEO.NB.CPA, SIMPLE.MAX.MAX, LEO.NB.MAX, LEO.NB.ALL, M1M2.AVG, LEO.I.FOR, LEO.I.MAX, LEO.I.ALL, LEO.O.FOR, LEO.O.MAX, LEO.O.ALL, PearsonCor, zeo.for,zeo.max,zeo.all, BLV.or.BilayerZscore)
    rownames(df.for.perm)=edge;

    # convert empty factors to NA
    df.for.perm=make.dataframe.numeric(df.for.perm,do.cols=3:ncol(df.for.perm))

    # only write headers on the first pass
    if (pm$perm.count == 2) {
           write.table(df.for.perm, file=pm$perm.file, col.names = NA, sep=",",dec=".",qmethod="double" )
    } else {
           # INVAR: pm$perm.count > 2
           write.table(df.for.perm, file=pm$perm.file,append=TRUE, col.names = FALSE, sep=",",dec=".",qmethod="double")
    }
}

   #
   #  Leave openning of excel file to upper level neo() wrapper now.
   #
   ### # open excel only after permutations done too.
   ### 
   ### if (pm$open.excel.at.end) {
   ###    if (pm$open.excel.now) { #  | pm$number.BLOCK.permutations==0) # comment out b/c was openning after each batch 
   ###        open.excel.file(excel.file,pm)
   ###    }
   ### }


# arrange snps in a line across the middle,
# and arrange non-snps in a half-circle underneath.

# PRE: the SNPs are the first 1:numsnps columns of datCombined
  
# start mostly inside the unit circle centered at zero
snp.y=rep(.7,numsnps)
# bump the evens up to try to get the labels visible again when we have lots of SNPS
odd.snps=(0:ceiling((numsnps/2)-1))*2+1;
snp.y[odd.snps] = snp.y[odd.snps]-.08;
snp.x=-1+2*(0:(numsnps-1))/(numsnps-1);
if (numsnps==1) {
   snp.x=0;
} 

radius=nc/2;
ang.inc = pi/(num.non.snps-1);
ang = (0:(num.non.snps-1))* ang.inc + pi;
non.snp.place=complex(argument=ang,modulus=1);
non.snp.x = Re(non.snp.place);
non.snp.y = Im(non.snp.place);

place=matrix((10+(c(snp.x,non.snp.x, snp.y, non.snp.y)+1)*40),nrow=nc,ncol=2)

# now, in case our PRE-condition that the SNPs were not first
# in terms of order, re-do placement to correct for that...
new.place=place;
new.place[snpcols,]=place[1:length(snp.x),]
new.place[non.snp.cols,]=place[(length(snp.x)+1):nrow(place),]
place=new.place


z$M = M;
z$super=super
z$dir.only = dir.only;
z$af = af; z$qlin = qlin;
#z$Vcdep = Vcdep;
z$cx = cx;
z$pcor = px;
z$covx = covx;
z$numsnps = numsnps
z$lcd.edge = lcd.edge
z$lcd.drop = lcd.drop
z$coor = place;
z$pm = pm;
z$non.snp.cols=non.snp.cols;
z$snpcols = snpcols;
z$edge.log = edge.log;
z$zeo = zeo;
z$nzeo = nzeo;
z$ZEO=zeo/sqrt(nzeo);
z$stats=df;
z$datCombined = datCombined; #retain imputed data

z$twostep.log=twofer.log;
z$final.eo=final.eo;
z$zeo.log=zeo.log;

z$eo.all.log=mylist.allvote[[3]]
z$eo.max.log=mylist.maxmax[[3]]
z$eo.for.log=mylist.mlreg[[3]]

z$losem.eo.all.log=losem.mylist.allvote$twofer.log
z$losem.eo.max.log=losem.mylist.maxmax$twofer.log
z$losem.eo.for.log=losem.mylist.mlreg$twofer.log

z$my.snp = r$my.snp

z$min.k.abs.pcor = min.k.abs.pcor
#z$min.k = min.k
z$min.k = matrix(cn[min.k],nrow=nrow(min.k),dimnames=dimnames(min.k.abs.pcor))


z$min.k.abs.pcor.ignore.supernode = min.k.abs.pcor.ignore.supernode;
z$min.k.nosup = matrix(cn[min.k.nosup],nrow=nrow(min.k.nosup),dimnames=dimnames(min.k.abs.pcor.ignore.supernode))
z$collider.zz = collider.zz

#cross.platform.windows();
#drawGraph(af, coor=place,adjust=T,beta=5);

#final.coor = neo.graph(z,adjust=repos.nodes);
#z$final.coor = final.coor;

topMax = min(10,length(qlin));
z$top10=z$stats[step2.eo.rank[1:topMax],]
z$bot10=z$stats[step2.eo.rank[(length(step2.eo.rank)):(length(step2.eo.rank)-topMax)],]

# since the composite can take a long time or crash...save before that.
z$pre.save.fname=paste(sep="",pm$neo.log.file,".z.PRE.COMPOSITE.FIT.rdat");
if (!pm$no.log) save(z,file=z$pre.save.fname);

# compute sem model fitting indices
if (exists("fit.full.model") && fit.full.model==TRUE) {
   pm.print(pm,paste("starting composite model sem() run...",date()));
   my.sem.struct = make.ram(M)
   sem.used.cols = sort(unique(my.sem.struct$used.cols));
   covx.sem=covx[sem.used.cols,sem.used.cols];
   sem.fit = try.sem(pm,my.sem.struct$the.ram,covx.sem,N=no.obs)
   z$sem.fit = sem.fit;
}
# now we open the excel file instead of printing the full thing
# pm.print(pm,z$stats)

#options(op) # reset (all) initial options 

z$post.save.fname=paste(sep="",pm$neo.log.file,".z.rdat");
if (!pm$no.log) save(z,file=z$post.save.fname);
if (!pm$no.log) pm.print(pm,paste(sep="","See file '",excel.file,"' and that directory for records of the run."));
pm.print(pm,paste("All done...",date()));


return(invisible(z));
} # end summary.NEO


####
# extract and re-use placement function for drawing graphs at the end.
####

if(exists("do.graph.placement") ) rm(do.graph.placement);
do.graph.placement=function(snpcols,numsnps, nc, num.non.snps,non.snp.cols) {

   # arrange snps in a line across the middle,
   # and arrange non-snps in a half-circle underneath.

   # PRE: the SNPs are the first 1:numsnps columns of datCombined
  
   # start mostly inside the unit circle centered at zero

  snp.y=rep(.7,numsnps)

  # bump the evens up to try to get the labels visible again when we have lots of SNPS
  odd.snps=(0:ceiling((numsnps/2)-1))*2+1;
  snp.y[odd.snps] = snp.y[odd.snps]-.08;
  snp.x=-1+2*(0:(numsnps-1))/(numsnps-1);

  if (numsnps==1) {
     snp.x=0;
     } 

     radius=nc/2;
     ang.inc = pi/(num.non.snps-1);
     ang = (0:(num.non.snps-1))* ang.inc + pi;
     non.snp.place=complex(argument=ang,modulus=1);
     non.snp.x = Re(non.snp.place);
     non.snp.y = Im(non.snp.place);

     place=matrix((10+(c(snp.x,non.snp.x, snp.y, non.snp.y)+1)*40),nrow=nc,ncol=2)

     # now, in case our PRE-condition that the SNPs were not first
     # in terms of order, re-do placement to correct for that...
     new.place=place;
     new.place[snpcols,]=place[1:length(snp.x),]
     new.place[non.snp.cols,]=place[(length(snp.x)+1):nrow(place),]

   place=new.place
}

# unit/sanity test
if(exists("neo.test") ) rm(neo.test);
neo.test=function(no.samples=500,...){
  set.seed(1); # make it reproducible
  s2=sqrt(2);
  SNP=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);

  SNP2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);
  SNP3=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);
  SNP4=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);

  env.var=0.5;
  A=SNP+rnorm(no.samples,sd=sqrt(env.var));
  B=A+rnorm(no.samples,sd=sqrt(env.var));

  datCombined=data.frame(SNP,A,B,SNP2,SNP3,SNP4);
  #attr(datCombined,"title")="unit test SNP.to.A.to.B"
  pm=neo.get.param()
  pm$run.title="unit test SNP.to.A.to.B"
  pm$no.log=FALSE
  pm$quiet=FALSE
  pm$A=2
  pm$B=3
  pm$top.N.snps.per.trait=1

  z=neo(datCombined,pm=pm,...);
} # end neo.test

if(exists("neo.test.perm") ) rm(neo.test.perm);
neo.test.perm=function(no.samples=500,...){
  set.seed(1); # make it reproducible
  s2=sqrt(2);
  SNP=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);

  SNP2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);
  SNP3=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);
  SNP4=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);

  env.var=0.5;
  A=SNP+rnorm(no.samples,sd=sqrt(env.var));
  B=A+rnorm(no.samples,sd=sqrt(env.var));

  datCombined=data.frame(SNP,A,B,SNP2,SNP3,SNP4);
  #attr(datCombined,"title")="unit test SNP.to.A.to.B"
  pm=neo.get.param()
  pm$run.title="unit test SNP.to.A.to.B"
  pm$no.log=FALSE
  pm$quiet=FALSE
  pm$A=2
  pm$B=3
  pm$top.N.snps.per.trait=2

  pm$number.BLOCK.permutations = 2
  pm$run.perms.on.sge = TRUE

  z=neo(datCombined,pm=pm,...);
} # end neo.test.perm


if(exists("neo.test.m1m2.avg") ) rm(neo.test.m1m2.avg);
neo.test.m1m2.avg=function(no.samples=500,...){
  set.seed(1); # make it reproducible
  s2=sqrt(2);
  SNP=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);

  SNP2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);
  SNP3=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);
  SNP4=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T);

  env.var=0.5;
  A=scale(SNP+SNP2+rnorm(no.samples,sd=sqrt(env.var)))
  B=scale(A+SNP3+SNP4+ rnorm(no.samples,sd=sqrt(env.var)))

  datCombined=data.frame(SNP,A,B,SNP2,SNP3,SNP4);

  pm=neo.get.param()
  pm$run.title="m1m2.avg test SNP.to.A.to.B"
  pm$no.log=FALSE
  pm$quiet=FALSE
  pm$A=2
  pm$B=3
  pm$top.N.snps.per.trait=3
  pm$do.m1m2.average = TRUE

  z=neo(datCombined,pm=pm,...);
} # end neo.test.m1m2.avg


if(exists("neo.test10") ) rm(neo.test10);
neo.test10=function(no.samples=500,...){
  pm=neo.get.param()
  pm$run.title="ChallengeData10run"
  datCombined=read.csv("ChallengeData10.csv")
  z=neo(datCombined,pm=pm);
} # end neo.test10


if(exists("neo.test2") ) rm(neo.test2);
neo.test2=function(no.samples=500,gene.to.env.var=10,generate=TRUE,pm=neo.get.param(),...){

   if (generate) {
     set.seed(1);
     numSNPgen=60;

   datSNP=data.frame(matrix(NA, nrow=no.samples, numSNPgen))

   names(datSNP)=paste("SNP",1:numSNPgen,"a", sep="");

   s2=sqrt(2)
   e2g=sqrt(1/gene.to.env.var);
   for (i in 1:numSNPgen) {
      SNPa=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T)
      datSNP[,i]=SNPa;
    } # end of for loop

    PCblue=scale(datSNP$SNP20a+ datSNP$SNP21a+ datSNP$SNP22a +datSNP$SNP23a + datSNP$SNP24a - datSNP$SNP25a - datSNP$SNP26a - datSNP$SNP27a - datSNP$SNP28a - datSNP$SNP29a)+e2g*rnorm(no.samples);

    PCturquoise=scale(PCblue+ datSNP$SNP10a + datSNP$SNP11a+ datSNP$SNP12a+ datSNP$SNP13a + datSNP$SNP14a - datSNP$SNP15a - datSNP$SNP16a - datSNP$SNP17a - datSNP$SNP18a - datSNP$SNP19a)+e2g*rnorm(no.samples);

    PCyellow=scale(PCblue+ datSNP$SNP50a+ datSNP$SNP51a+ datSNP$SNP52a+ datSNP$SNP53a + datSNP$SNP54a + datSNP$SNP55a - datSNP$SNP56a - datSNP$SNP57a - datSNP$SNP58a - datSNP$SNP59a)+e2g*rnorm(no.samples);

    PCbrown= scale(datSNP$SNP30a+ datSNP$SNP31a+ datSNP$SNP32a+ datSNP$SNP33a + datSNP$SNP34a - datSNP$SNP35a - datSNP$SNP36a - datSNP$SNP37a - datSNP$SNP38a - datSNP$SNP39a)+e2g*rnorm(no.samples);

    Trait=scale(PCbrown + datSNP$SNP1a+ datSNP$SNP2a+ datSNP$SNP3a + datSNP$SNP4a + datSNP$SNP5a - datSNP$SNP6a - datSNP$SNP7a - datSNP$SNP8a - datSNP$SNP9a)+e2g*rnorm(no.samples)+e2g*rnorm(no.samples);

    PCgreen=scale(Trait  + datSNP$SNP40a + datSNP$SNP41a+ datSNP$SNP42a+ datSNP$SNP43a + datSNP$SNP44a - datSNP$SNP45a - datSNP$SNP46a - datSNP$SNP47a - datSNP$SNP48a - datSNP$SNP49a)+e2g*rnorm(no.samples);

#    datCombined=data.frame(Trait, PCturquoise,PCblue,PCbrown, PCgreen, PCyellow,datSNP)

     SNP20a=datSNP$SNP20a;
     SNP21a=datSNP$SNP21a;
#     datCombined=data.frame(Trait, PCturquoise,PCblue,PCbrown, PCgreen, PCyellow,SNP20a,SNP21a);
     datCombined=data.frame(cbind(Trait, PCturquoise,PCblue,PCbrown, PCgreen, PCyellow,datSNP))

     write.csv(datCombined,"neo.test2.csv",row.names=F);
     } else {
        datCombined=read.csv("neo.test2.csv");
     }

     pm$run.title="unit.test2"
     z=neo(datCombined,pm=pm,...);

  z
} # end neo.test2


# changing neo.test2 will probably break this little test of marker assignment forcing...but it's easily fixed.
if(exists("neo.test.force") ) rm(neo.test.force);
neo.test.force=function(){
  pm=neo.get.param()
  pm$A = 1
  pm$forced.MA.colnum = c(7,8) # make this point to a SNP!
  neo.test2(pm=pm)
}

if(exists("neo.test2ab") ) rm(neo.test2ab);
neo.test2ab=function(no.samples=500,gene.to.env.var=10,generate=TRUE,...){

   if (generate) {
     set.seed(1);
     numSNPgen=60;

   datSNP=data.frame(matrix(NA, nrow=no.samples, numSNPgen))

   names(datSNP)=paste("SNP",1:numSNPgen,"a", sep="");

   s2=sqrt(2)
   e2g=sqrt(1/gene.to.env.var);
   for (i in 1:numSNPgen) {
      SNPa=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T)
      datSNP[,i]=SNPa;
    } # end of for loop

    PCblue=scale(datSNP$SNP20a+ datSNP$SNP21a+ datSNP$SNP22a +datSNP$SNP23a + datSNP$SNP24a - datSNP$SNP25a - datSNP$SNP26a - datSNP$SNP27a - datSNP$SNP28a - datSNP$SNP29a)+e2g*rnorm(no.samples);

    PCturquoise=scale(PCblue+ datSNP$SNP10a + datSNP$SNP11a+ datSNP$SNP12a+ datSNP$SNP13a + datSNP$SNP14a - datSNP$SNP15a - datSNP$SNP16a - datSNP$SNP17a - datSNP$SNP18a - datSNP$SNP19a)+e2g*rnorm(no.samples);

    PCyellow=scale(PCblue+ datSNP$SNP50a+ datSNP$SNP51a+ datSNP$SNP52a+ datSNP$SNP53a + datSNP$SNP54a + datSNP$SNP55a - datSNP$SNP56a - datSNP$SNP57a - datSNP$SNP58a - datSNP$SNP59a)+e2g*rnorm(no.samples);

    PCbrown= scale(datSNP$SNP30a+ datSNP$SNP31a+ datSNP$SNP32a+ datSNP$SNP33a + datSNP$SNP34a - datSNP$SNP35a - datSNP$SNP36a - datSNP$SNP37a - datSNP$SNP38a - datSNP$SNP39a)+e2g*rnorm(no.samples);

    Trait=scale(PCbrown + datSNP$SNP1a+ datSNP$SNP2a+ datSNP$SNP3a + datSNP$SNP4a + datSNP$SNP5a - datSNP$SNP6a - datSNP$SNP7a - datSNP$SNP8a - datSNP$SNP9a)+e2g*rnorm(no.samples)+e2g*rnorm(no.samples);

    PCgreen=scale(Trait  + datSNP$SNP40a + datSNP$SNP41a+ datSNP$SNP42a+ datSNP$SNP43a + datSNP$SNP44a - datSNP$SNP45a - datSNP$SNP46a - datSNP$SNP47a - datSNP$SNP48a - datSNP$SNP49a)+e2g*rnorm(no.samples);

#    datCombined=data.frame(Trait, PCturquoise,PCblue,PCbrown, PCgreen, PCyellow,datSNP)

     SNP20a=datSNP$SNP20a;
     SNP21a=datSNP$SNP21a;
     datCombined=data.frame(Trait, PCturquoise,PCblue,PCbrown, PCgreen, PCyellow,SNP20a,SNP21a);

     write.csv(datCombined,"neo.test2ab.csv",row.names=F);
     } else {
        datCombined=read.csv("neo.test2ab.csv");
     }
     pm=neo.get.param()
     cn=colnames(datCombined)
     pm$A = match("PCbrown",cn)
     pm$B = match("Trait",cn)
     z=neo(datCombined,pm=pm,...);

} # end neo.test2ab


#
# Simulate the full gamma-omega model with weighted discrete simulated SNPs into
# A with alpha weights and into B with beta weights.
#
#
#
if(exists("neo.test3") ) rm(neo.test3);
neo.test3=function(no.samples=200,gene.to.env.var=10,generate=TRUE,...){

   set.seed(1);
   N=no.samples

   Ktop = 10
   nsim=10
  # wrap with a sim loop
  for (sim in 1:nsim) {

   confounded=FALSE # set to TRUE to simulate under null hypothesis
    
   do.a.snp=5
   do.b.snp=5
   noise.a.snp = 10
   noise.b.snp = 10

   # simulate SNPs nearby a QTL peak, for each do.a.snp and for each do.b.snp add this many:
   nearby.a.snps = 6
   nearby.b.snps = 6

   numAsnp=do.a.snp+noise.a.snp 
   numBsnp=do.b.snp+noise.b.snp 

   # the weights for the SNPs coming into A
   alpha.A.snp = c(rep(sqrt(1/do.a.snp),do.a.snp),rep(0,noise.a.snp) )

   # the weights for the SNPs coming into B
   beta.B.snp  = c(rep(sqrt(1/do.b.snp),do.b.snp), rep(0,noise.b.snp) )

   datSNP=data.frame(matrix(NA, nrow=no.samples, numAsnp+numBsnp))

   names(datSNP)[1:numAsnp] = paste("SNP.A.",1:numAsnp, sep="");
   Bspan = (1:numBsnp) + numAsnp
   names(datSNP)[Bspan] = paste("SNP.B.",1:numBsnp, sep="");
   SNP.A=SNP.B=0

   for (i in 1:numAsnp) {
      SNPa=scale(sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T))
      SNP.A=SNP.A+ (alpha.A.snp[i] * SNPa)
      datSNP[,i]=SNPa;
    } # end of for loop

   k=1
   for (i in Bspan) {
      SNPa=scale(sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T))
      SNP.B=SNP.B+ (beta.B.snp[k] * SNPa)
      datSNP[,i]=SNPa;
      k=k+1
    } # end of for loop

   # now generate the nearby SNPs too
   nearB = do.b.snp*nearby.b.snps 
   nearA = do.a.snp*nearby.a.snps
   nearTot = nearA+nearB
   nearbySNP=data.frame(matrix(NA, nrow=no.samples, nearTot))
   names(nearbySNP)[1:nearA] = paste("SNP.near.A.",1:nearA, sep="");
   nearBspan = (1:nearB) + nearA
   names(nearbySNP)[nearBspan] = paste("SNP.near.B.",1:nearB, sep="");

   percentage.off.range = c(.75,.95) # approximate correlation values with the neighboring true signal SNP
   if (nearA > 0) {
    for (i in 1:nearA) {
      which.a = 1+floor((i-1)/nearby.a.snps)
      p.off=runif(1,min=percentage.off.range[1], max=percentage.off.range[2])
      num.to.alter = floor(N*(1-p.off))
      SNPa=datSNP[,which.a]
      
      new.val=sample(c(0,1,2), size=num.to.alter, prob=c(.25,.5,.25) ,replace=T)
      SNPa[sample(1:N,size=num.to.alter,replace=F)] = new.val

      nearbySNP[,i]=SNPa;
      colnames(nearbySNP)[i]=paste(sep="","SNP.neighbor.",i,".near.",colnames(datSNP)[which.a])
   }
   }

   if (nearB > 0 ) {
    for (i in nearBspan) {
      which.b = 1+floor((i-(nearA+1))/nearby.b.snps)
      p.off=runif(1,min=percentage.off.range[1], max=percentage.off.range[2])

      num.to.alter = floor(N*(1-p.off))
      SNPa=datSNP[,Bspan[which.b]]
      
      new.val=sample(c(0,1,2), size=num.to.alter, prob=c(.25,.5,.25) ,replace=T)
      SNPa[sample(1:N,size=num.to.alter,replace=F)] = new.val

      nearbySNP[,i]=SNPa;
      colnames(nearbySNP)[i]=paste(sep="","SNP.neighbor.",i,".near.",colnames(datSNP)[Bspan[which.b]])
    }
   }
   
     ha=hb=.5
     if (confounded) { g2=0.4 } else { g2=0}
     g = sqrt(g2)
     errb = 1-hb-g2
     cor.ab = .4
     w= cor.ab - g2;
     erra = 1-ha-w^2-g^2-2*(g*g*w)
     C = rnorm(N);

    B= sqrt(hb) * scale(SNP.B) + rnorm(N,sd=sqrt(errb)) + g*C;
    A= sqrt(ha) * scale(SNP.A) + rnorm(N,sd=sqrt(erra)) + g*C + w*B;

#     datCombined=data.frame(cbind(datSNP,A,B));
     datCombined=data.frame(cbind(datSNP,nearbySNP,A,B));

#     write.csv(datCombined,"neo.test3.csv",row.names=F);

pm=neo.get.param()
pm$no.log=FALSE
pm$quiet=TRUE
pm$top.N.snps.per.trait=1
pm$orthog.search.depth=0
pm$run.title="test3run"

#     z=neo(datCombined,pm=pm);

# ========================================
# evaluate forward first:
# ========================================
 forsnp.leo.nb.for=matrix(NA,ncol=Ktop,nrow=nsim)
 forsnp.leo.nb.max=matrix(NA,ncol=Ktop,nrow=nsim)
 forsnp.leo.nb.all=matrix(NA,ncol=Ktop,nrow=nsim)
 zlist=list()

 for (k in 1:Ktop) {
 pm$do.snp.pick.greedy = FALSE
 pm$do.snp.pick.fstep  = TRUE

 pm$top.N.snps.per.trait=1
 pm$orthog.search.depth=k-1

 #  pm$force.cutree.count = 2
 #  pm$orthog.search.depth = k
  z1 = neo(datCombined,pm=pm)
#  zlist[[k]]=z1
  forsnp.leo.nb.for[sim,k]=anac(z1$stats["B -> A","LEO.NB.FOR"])
  forsnp.leo.nb.max[sim,k]=anac(z1$stats["B -> A","LEO.NB.MAX"])
  forsnp.leo.nb.all[sim,k]=anac(z1$stats["B -> A","LEO.NB.ALL"])
 }

## # results: with no neighboring SNPs, 5 SNPs into B trait:
## #forsnp.leo.nb.for = c(0.647,1.780,3.780,2.950,3.960,3.960,3.960,3.960,3.960,3.960)
## #forsnp.leo.nb.max = c(0.647,1.570,1.570,1.270,1.270,1.270,1.270,1.270,1.270,1.270)
## #forsnp.leo.nb.all = c(3.33,3.10,3.09,2.22,2.25,2.25,2.25,2.25,2.25,2.25)

## # results: with 6 neighboring SNPs:
## #forsnp.leo.nb.for= c(0.59,1.36,1.51,1.49,2.59,2.59,2.59,2.59,2.59,2.59)
## #forsnp.leo.nb.max = c(0.590,0.316,0.328,0.328,0.328,0.328,0.328,0.328,0.328,0.328)
## #forsnp.leo.nb.all = c(4.04,2.20,1.68,1.15,1.13,1.13,1.13,1.13,1.13,1.13)


## postscript(file="forward.step.snp.impact.with.6.neighbors.ps",horizontal=FALSE) 
## plot(1:length(forsnp.leo.nb.for),forsnp.leo.nb.for,type="l",ylab="LEO.NB",xlab="Number of Forward-step SNPs used on a 5 SNP trait",ylim=c(0,5),cex.lab=1.5)
## points(1:length(forsnp.leo.nb.for),forsnp.leo.nb.for,pch=21,bg="red",cex=1.5)

## lines(1:length(forsnp.leo.nb.all),forsnp.leo.nb.all)
## points(1:length(forsnp.leo.nb.all),forsnp.leo.nb.all,pch=22,cex=1.5,bg="yellow") # square

## lines(1:length(forsnp.leo.nb.max),forsnp.leo.nb.max)
## points(1:length(forsnp.leo.nb.max),forsnp.leo.nb.max,pch=23,cex=1.5,bg="blue") # diamond

## legend(7,4.5, c("LEO.NB.FOR", "LEO.NB.ALL", "LEO.NB.MAX"), cex=1.2, pch=c(21,22,23),pt.bg=c("red","yellow","blue"))
## title("With six neighboring SNPs, impact on \nForward-step SNP selection",cex.main=1.5,sub="ha=hb=0.5; omega=0.4; gamma=0",cex.sub=.75)
## dev.off()
   

# ========================================
# evaluate greedy next
# ========================================

 greedysnp.leo.nb.for=matrix(NA,ncol=Ktop,nrow=nsim)
 greedysnp.leo.nb.max=matrix(NA,ncol=Ktop,nrow=nsim)
 greedysnp.leo.nb.all=matrix(NA,ncol=Ktop,nrow=nsim)
 zlist=list()

 for (k in 1:Ktop) {
 pm$do.snp.pick.greedy = TRUE
 pm$do.snp.pick.fstep  = FALSE

 pm$top.N.snps.per.trait=k
 pm$orthog.search.depth=0

 #  pm$force.cutree.count = 2
 #  pm$orthog.search.depth = k
  z1 = neo(datCombined,pm=pm)
#  zlist[[k]]=z1
  greedysnp.leo.nb.for[sim,k]=anac(z1$stats["B -> A","LEO.NB.FOR"])
  greedysnp.leo.nb.max[sim,k]=anac(z1$stats["B -> A","LEO.NB.MAX"])
  greedysnp.leo.nb.all[sim,k]=anac(z1$stats["B -> A","LEO.NB.ALL"])
 }

## # no neighbors:
## greedysnp.leo.nb.for =c(0.647,1.760,2.010,2.950,3.960,3.960,3.960,3.960,3.960,3.960)
## greedysnp.leo.nb.all=c(3.33,3.56,2.56,2.22,2.25,2.25,2.25,2.25,2.25,2.25)
## greedysnp.leo.nb.max=c(0.647,1.580,1.570,1.270,1.270,1.270,1.270,1.270,1.270,1.270)

## postscript(file="greedy.snp.impact.no.neighbors.ps",horizontal=FALSE) 
## plot(1:length(greedysnp.leo.nb.for),greedysnp.leo.nb.for,type="l",ylab="LEO.NB",xlab="Number of Greedy SNPs used on a 5 SNP trait",ylim=c(0,5),cex.lab=1.5)
## points(1:length(greedysnp.leo.nb.for),greedysnp.leo.nb.for,pch=21,bg="red",cex=1.5)

## lines(1:length(greedysnp.leo.nb.all),greedysnp.leo.nb.all)
## points(1:length(greedysnp.leo.nb.all),greedysnp.leo.nb.all,pch=22,cex=1.5,bg="yellow") # square

## lines(1:length(greedysnp.leo.nb.max),greedysnp.leo.nb.max)
## points(1:length(greedysnp.leo.nb.max),greedysnp.leo.nb.max,pch=23,cex=1.5,bg="blue") # diamond

## legend(7,5, c("LEO.NB.FOR", "LEO.NB.ALL", "LEO.NB.MAX"), cex=1.2, pch=c(21,22,23),pt.bg=c("red","yellow","blue"))
## title("No neighboring SNPs, impact of \nGreedy SNP selection",cex.main=1.5,sub="ha=hb=0.5; omega=0.4; gamma=0",cex.sub=.75)
## dev.off()


## # with 6 neighbors:
## greedysnp.leo.nb.for=c(0.590,0.834,1.920,2.220,2.220,3.010,2.970,2.970,1.830,1.830)
## greedysnp.leo.nb.max=c(0.5900,0.7450,0.7450,0.3160,-0.1280,0.3280,0.0354,0.0354,0.3280,0.3280)
## greedysnp.leo.nb.all=c(4.04,3.40,3.25,2.34,1.56,1.86,1.45,1.54,1.30,1.38)


## postscript(file="greedy.snp.impact.with.6.neighbors.ps",horizontal=FALSE) 
## plot(1:length(greedysnp.leo.nb.for),greedysnp.leo.nb.for,type="l",ylab="LEO.NB",xlab="Number of Greedy SNPs used on a 5 SNP trait",ylim=c(0,5),cex.lab=1.5)
## points(1:length(greedysnp.leo.nb.for),greedysnp.leo.nb.for,pch=21,bg="red",cex=1.5)

## lines(1:length(greedysnp.leo.nb.all),greedysnp.leo.nb.all)
## points(1:length(greedysnp.leo.nb.all),greedysnp.leo.nb.all,pch=22,cex=1.5,bg="yellow") # square

## lines(1:length(greedysnp.leo.nb.max),greedysnp.leo.nb.max)
## points(1:length(greedysnp.leo.nb.max),greedysnp.leo.nb.max,pch=23,cex=1.5,bg="blue") # diamond

## legend(7,5, c("LEO.NB.FOR", "LEO.NB.ALL", "LEO.NB.MAX"), cex=1.2, pch=c(21,22,23),pt.bg=c("red","yellow","blue"))
## title("With six neighboring SNPs, impact of \nGreedy SNP selection",cex.main=1.5,sub="ha=hb=0.5; omega=0.4; gamma=0",cex.sub=.75)
## dev.off()


# ========================================
# evaluate forward and greedy together
# ========================================
 for.gre.snp.leo.nb.for=matrix(NA,ncol=Ktop,nrow=nsim)
 for.gre.snp.leo.nb.max=matrix(NA,ncol=Ktop,nrow=nsim)
 for.gre.snp.leo.nb.all=matrix(NA,ncol=Ktop,nrow=nsim)
 zlist=list()

 for (k in 1:Ktop) {
 pm$do.snp.pick.greedy = TRUE
  pm$do.snp.pick.fstep  = TRUE

 pm$top.N.snps.per.trait=k
 pm$orthog.search.depth=k-1

  z1 = neo(datCombined,pm=pm)
#  zlist[[k]]=z1
  for.gre.snp.leo.nb.for[sim,k]=anac(z1$stats["B -> A","LEO.NB.FOR"])
  for.gre.snp.leo.nb.max[sim,k]=anac(z1$stats["B -> A","LEO.NB.MAX"])
  for.gre.snp.leo.nb.all[sim,k]=anac(z1$stats["B -> A","LEO.NB.ALL"])
 }

 print(paste(sep="","on ",sim," pass, out of ",nsim," at ",date()))
 } # end sim loop

   r=list()
r$forsnp.leo.nb.for=forsnp.leo.nb.for
r$forsnp.leo.nb.max=forsnp.leo.nb.max
r$forsnp.leo.nb.all   =forsnp.leo.nb.all   
r$for.gre.snp.leo.nb.for=for.gre.snp.leo.nb.for
r$for.gre.snp.leo.nb.max=for.gre.snp.leo.nb.max
r$for.gre.snp.leo.nb.all=for.gre.snp.leo.nb.all
r$greedysnp.leo.nb.for=greedysnp.leo.nb.for
r$greedysnp.leo.nb.max=greedysnp.leo.nb.max
r$greedysnp.leo.nb.all=greedysnp.leo.nb.all

   return(r)
   
## # no neighbors:
## for.gre.snp.leo.nb.for=c(0.59,1.36,1.51,2.57,2.57,2.57,2.57,2.57,2.57,2.57)
## for.gre.snp.leo.nb.max=c(0.590,0.316,0.328,0.328,0.328,0.328,0.328,0.328,0.328,0.328)
## for.gre.snp.leo.nb.all=c(4.04,2.20,1.68,1.22,1.22,1.22,1.22,1.22,1.22,1.22)

## postscript(file="forward.and.greedy.snp.impact.no.neighbors.ps",horizontal=FALSE) 
## plot(1:length(for.gre.snp.leo.nb.for),for.gre.snp.leo.nb.for,type="l",ylab="LEO.NB",xlab="Number of Forward and Greedy SNPs used on a 5 SNP trait",ylim=c(0,5),cex.lab=1.4)
## points(1:length(for.gre.snp.leo.nb.for),for.gre.snp.leo.nb.for,pch=21,bg="red",cex=1.5)

## lines(1:length(for.gre.snp.leo.nb.all),for.gre.snp.leo.nb.all)
## points(1:length(for.gre.snp.leo.nb.all),for.gre.snp.leo.nb.all,pch=22,cex=1.5,bg="yellow") # square

## lines(1:length(for.gre.snp.leo.nb.max),for.gre.snp.leo.nb.max)
## points(1:length(for.gre.snp.leo.nb.max),for.gre.snp.leo.nb.max,pch=23,cex=1.5,bg="blue") # diamond

## legend(7,5, c("LEO.NB.FOR", "LEO.NB.ALL", "LEO.NB.MAX"), cex=1.2, pch=c(21,22,23),pt.bg=c("red","yellow","blue"))
## title("No neighboring SNPs, impact of \nForward and Greedy SNP selection",cex.main=1.5,sub="ha=hb=0.5; omega=0.4; gamma=0",cex.sub=.75)
## dev.off()



## # with neighbors
## for.gre.snp.leo.nb.for=c(0.59,2.22,3.01,2.89,2.59,2.59,2.59,2.59,2.59,2.59)
## for.gre.snp.leo.nb.max=c(0.5900,0.3160,0.3280,0.3280,0.0354,0.0354,0.0354,0.3280,0.0812,0.0812)
## for.gre.snp.leo.nb.all=c(4.040,2.340,1.960,1.340,0.969,0.969,0.969,1.300,1.160,1.270)

## postscript(file="forward.and.greedy.snp.impact.with.6.neighbors.ps",horizontal=FALSE) 
## plot(1:length(for.gre.snp.leo.nb.for),for.gre.snp.leo.nb.for,type="l",ylab="LEO.NB",xlab="Number of Forward and Greedy SNPs used on a 5 SNP trait",ylim=c(0,5),cex.lab=1.4)
## points(1:length(for.gre.snp.leo.nb.for),for.gre.snp.leo.nb.for,pch=21,bg="red",cex=1.5)

## lines(1:length(for.gre.snp.leo.nb.all),for.gre.snp.leo.nb.all)
## points(1:length(for.gre.snp.leo.nb.all),for.gre.snp.leo.nb.all,pch=22,cex=1.5,bg="yellow") # square

## lines(1:length(for.gre.snp.leo.nb.max),for.gre.snp.leo.nb.max)
## points(1:length(for.gre.snp.leo.nb.max),for.gre.snp.leo.nb.max,pch=23,cex=1.5,bg="blue") # diamond

## legend(7,5, c("LEO.NB.FOR", "LEO.NB.ALL", "LEO.NB.MAX"), cex=1.2, pch=c(21,22,23),pt.bg=c("red","yellow","blue"))
## title("With six neighboring SNPs, impact of \nForward and Greedy SNP selection",cex.main=1.5,sub="ha=hb=0.5; omega=0.4; gamma=0",cex.sub=.75)
## dev.off()

# with 6 neighbors, and with 0% causal signal, i.e. total confounding:
##  forsnp.leo.nb.for    =c(-1.88,-3.44,-4.39,-7.47,-13,-13,-13,-13,-13,-13)
##  forsnp.leo.nb.max    =c(-1.88,-2.38,-1.04,-1.04,-1.04,-1.04,-1.04,-1.04,-1.04,-1.04)
##  forsnp.leo.nb.all    =c(-1.26,-0.601,-1.06,-1.25,-1.46,-1.35,-1.26,-1.53,-1.53,-1.53)
##  for.gre.snp.leo.nb.for =c(-1.88,-3.44,-4.39,-12.4,-13,-13,-13,-13,-13.9,-13.9)
##  for.gre.snp.leo.nb.max =c(-1.88,-2.38,-1.04,-1.04,-1.04,-1.04,-1.04,-1.04,-1.04,-1.04)
##  for.gre.snp.leo.nb.all =c(-1.26,-0.601,-1.28,-1.39,-1.5,-1.4,-1.19,-1.36,-1.25,-1.25)
##  greedysnp.leo.nb.for =c(-1.88,-1.49,-1.43,-3.04,-5.83,-5.76,-7.42,-12.2,-10.2,-10.2)
##  greedysnp.leo.nb.max =c(-1.88,-1.88,-1.28,-1.28,-1.28,-1.65,-1.65,-1.65,-1.04,-1.04)
##  greedysnp.leo.nb.all =c(-1.26,-0.815,-0.92,-1.05,-1.1,-1.75,-1.89,-1.79,-1.44,-1.56)



} # end neo.test3


#
# make a directory for SGE output (very long) within the logfile directory
#   as specified by dirname(pm$neo.log.file)
# 
#
if(exists("create.logpath.enclosed.sge.directory") ) rm(create.logpath.enclosed.sge.directory);
create.logpath.enclosed.sge.directory=function(pm) {

   starting.dir = dirname(pm$neo.log.file)
   new.sge.dir = path.os.convert(paste(sep="",getwd(),"/",starting.dir,"/sge.neo.job.output/"))

   dir.create(new.sge.dir,showWarnings=FALSE,recursive=TRUE)

   new.sge.dir
}


#
# neo.mkdir.logfile(): put all log files in a new directory with timestamp
#
if(exists("neo.mkdir.logfile") ) rm(neo.mkdir.logfile);
neo.mkdir.logfile=function(run.title="") {

   one.date=date() 
   neo.logs.dir=gsub(":",".",gsub(" ","_",paste(sep="","neo.logs.",one.date)))
   run.title=gsub(" ","_",run.title)
   neo.logs.dir = paste(sep="",neo.logs.dir,".",run.title)

   dir.create(neo.logs.dir,showWarnings=FALSE,recursive=TRUE)

# old way, paths too long:   neo.log.file=paste(sep="\\",neo.logs.dir, gsub(":",".",gsub(" ","_",paste(sep="","neo.logfile.",one.date,".",run.title))))
   neo.log.file=paste(sep="\\",neo.logs.dir, "neolog")
   path.os.convert(neo.log.file); # customize to OS
}

# convert a path to unix "/" or dos "\\" depending on the OS detected.
if(exists("path.os.convert") ) rm(path.os.convert);
path.os.convert=function(path) {
  if (.Platform $OS.type == "windows") {
      return(gsub("/","\\\\",path));
  }
  gsub("\\\\","/",path)
}

if(exists("neo.guess.snps") ) rm(neo.guess.snps);
neo.guess.snps=function(x,pm) {
   cn=colnames(x);
   # guess at snpcols if not given
   snpcols=which(tolower(substr(cn,1,3)) == "snp");

   if (length(snpcols)==0) { 
         stop("snpcols not given and no data columns start with 'SNP': Must have SNPs, at the moment.");
   } else {
            pm.print(pm,"snpcols not specified manually, looking for 'SNP' prefixed column names"); 
         if (length(snpcols) < 100) {
         pm.print(pm,paste("   ...we'll suppose that these ",length(snpcols)," are our markers/SNPS:"));
         pm.print(pm,cn[snpcols]);
         } else {
         pm.print(pm, "too many (> 100) snpcols, skipping their printout.")
         }
   }
   snpcols
}

# helper function to find excel path name
if(exists("get.excel.file") ) rm(get.excel.file);
get.excel.file=function(pm) {
   if (dirname(pm$neo.log.file)==".") {
     # not working inside a log directory??? try to return something sensible nonetheless...
     return(paste(sep="",pm$neo.log.file,".csv"))
   }
   excel.file=path.os.convert(paste(sep="",dirname(pm$neo.log.file),"\\",basename(dirname(pm$neo.log.file)),".csv"))
}

# helper function to silence printing based on pm$quiet
if(exists("pm.print") ) rm(pm.print);
pm.print=function(pm,x) {
   if (pm$quiet) return(NULL);
   print(x)
}

if(exists("pm.cat") ) rm(pm.cat);
pm.cat=function(pm,x) {
   if (pm$quiet) return(NULL);
   cat(x)
}

if(exists("pm.warning") ) rm(pm.warning);
pm.warning=function(pm,x) {
   if (pm$quiet) return(NULL);
   warning(x)
}


# Guess at Trait if not given.
if(exists("neo.guess.traits") ) rm(neo.guess.traits);
neo.guess.traits=function(datCombined,pm) {

       if (!pm$use.traitcols.guess) return(c()) # bail early if we don't want to guess.

       cn=colnames(datCombined);
 
       traitcols = c(which(tolower(substr(cn,1,1))=="t"),which(tolower(substr(cn,1,1))=="y"));
       if (length(traitcols)==0) {
           pm.print(pm,"Hmmm...traitcols not given and no data columns start with 't' or 'y' : not using any special trait treatment."); 
       }

       pm.print(pm,"traitcols not specificed manually, looking for 'trait' or 'y' column-names. Using:");
       pm.print(pm,cn[traitcols]);

       traitcols
}

# detect bad column names that through off the as.formula construction in forward.step.pick.snp
if(exists("starts.with.zero.thru.nine") ) rm(starts.with.zero.thru.nine);
starts.with.zero.thru.nine=function(a) {
  if (substr(a,1,1) == "0") return(TRUE);
  if (substr(a,1,1) == "1") return(TRUE);
  if (substr(a,1,1) == "2") return(TRUE);
  if (substr(a,1,1) == "3") return(TRUE);
  if (substr(a,1,1) == "4") return(TRUE);

  if (substr(a,1,1) == "5") return(TRUE);
  if (substr(a,1,1) == "6") return(TRUE);
  if (substr(a,1,1) == "7") return(TRUE);
  if (substr(a,1,1) == "8") return(TRUE);
  if (substr(a,1,1) == "9") return(TRUE);
  FALSE;
}

# helper function to convert from factors to numeric via character--affected by options(digits = 4); getOptions("digits") to check.
if(exists("anac") ) rm(anac);
anac=function(x) as.numeric(as.character(x))

if(exists("noise.snp.sim") ) rm(noise.snp.sim);
noise.snp.sim=function(param=NULL, no.samples=200, num.reps=50,pm=neo.get.param()) {


ha=hb=.6

rho.ab=.5
percent.max.omega=1; # .6667
omega=percent.max.omega*rho.ab
gamma2=rho.ab-omega

gamma=sqrt(gamma2)

eA=1-2*omega*gamma^2 - gamma^2 - omega^2 - ha

eB=1-gamma2 - hb;


  q=data.frame(ha,hb,gamma2,omega,eA,eB,rho.ab,percent.max.omega)

# just for development DEBUG, take out later, so function works...
param=q
no.samples=200
pm=neo.get.param()
pm$skip.snp.select = TRUE
pm$quiet = TRUE
pm$run.title="trashme"
pm$cor.ind.th = .3 # so max doesn't exit away?
pm$zeo.proxy.for.leo.th = NA


  num.reps=5
  num.param=1
  total.rows = num.reps*num.param

  BLV.AtoB=rep(NA,total.rows)
  BLV.BtoA=rep(NA,total.rows)

  leo.nb.for = rep(NA,total.rows)
  leo.i.for  = rep(NA,total.rows)
  leo.o.for  = rep(NA,total.rows)

  zeo.for = rep(NA,total.rows)
  zeo.max  = rep(NA,total.rows)
  zeo.all  = rep(NA,total.rows)

  leo.nb.max = rep(NA,total.rows)
  leo.i.max  = rep(NA,total.rows)
  leo.o.max  = rep(NA,total.rows)

  zeo.nb.max = rep(NA,total.rows)
  zeo.i.max  = rep(NA,total.rows)
  zeo.o.max  = rep(NA,total.rows)

  leo.nb.all = rep(NA,total.rows)
  leo.i.all  = rep(NA,total.rows)
  leo.o.all  = rep(NA,total.rows)

  zeo.nb.all = rep(NA,total.rows)
  zeo.i.all  = rep(NA,total.rows)
  zeo.o.all  = rep(NA,total.rows)

  # needed for anac() to not round too much....
  options(digits = 4); # affects the conversion of spreadsheet columns from factors to character to numeric, at the character to numeric step.

  
 for (i in 1:num.param) {
   for (j in 1:num.reps) {
      k=(num.reps*(i-1))+j # where to store

     z =list()

     z$ha =   ha     = param$ha[i]
     z$hb =   hb     = param$hb[i]
     z$gamma2 =   gamma2 = param$gamma2[i]
     z$omega =   omega  = param$omega[i]
     z$eA =   eA =     param$eA[i]
     z$eB =   eB =     param$eB[i]
     z$rho.ab =   rho.ab = param$rho.ab[i]

     param.string=paste(sep="","ha=",ha, " hb=",hb, " g2=",gamma2, " eA=",eA," eB=",eB, " rho.ab=",rho.ab, " w=",omega)
   
#              varB= gamma2 + eB + hb;
#              varA = TA + 2*omega*gamma*gamma + gamma*gamma + omega*omega*varB + eA;

              # if B is standardized to var(B)=1, then....

#  std.varA = TA + gamma^2 + omega^2 + 2*gamma*gamma*omega + eA; # set == 1
#  covAB = omega*TB+gamma*gamma*(1+omega)+omega*eB
#  covAC = gamma+ omega*gamma
#  covBC = gamma

#  cov.SNPB.A = sqrt(TB)*omega

#  rho.AB = gamma^2 + omega

  SNPA=rnorm(no.samples);
  SNPB=rnorm(no.samples);

   s2=sqrt(2)
   SNPA1=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPA2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPA3=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPA4=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPA5=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPA6=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPA7=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPA8=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPA9=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPA10=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1

  SNPA = scale(SNPA1+SNPA2+SNPA3+SNPA4+SNPA5+SNPA6+SNPA7+SNPA8+SNPA9+SNPA10)

#   MB1=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
#   MB2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
#   MB2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
#   SNPB = scale(MB1+MB2)


   SNPB1=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPB2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPB3=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPB4=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPB5=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPB6=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPB7=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPB8=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPB9=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
   SNPB10=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1

  SNPB = scale(SNPB1+SNPB2+SNPB3+SNPB4+SNPB5+SNPB6+SNPB7+SNPB8+SNPB9+SNPB10)


  C=rnorm(no.samples);

  errA=rnorm(no.samples,sd=sqrt(eA));
  errB=rnorm(no.samples,sd=sqrt(eB));

  B=errB+sqrt(gamma2)*C+ sqrt(hb)*SNPB;

  A=errA+omega*B+sqrt(gamma2)*C+sqrt(ha)*SNPA;


     x=data.frame(SNPA1,SNPA2,SNPA3,SNPA4,SNPA5,SNPA6,SNPA7,SNPA8,SNPA9,SNPA10,SNPB,A,B);
#  x=data.frame(SNPA,SNPB,A,B);

  z$neo =neo(x,pm=pm,traitcols=3)   
#    z$x=x
##        z$fit1mb = compare.local.sems(pm=pm,M.col=2,A.col=3,B.col=4,x=x,no.obs=no.samples)
#        z$fit1ma = compare.local.sems(pm=pm,M.col=1,A.col=3,B.col=4,x=x,no.obs=no.samples)
#        z$fit2m =local.sem.four.var.m1m2(pm=pm,MA.col=1,A.col=3,B.col=4,x=x,MB.col=2,no.obs=N,fit.models=c(1,2,3,4));
    z$param.string = param.string;



  b.a.row=which(rownames(z$neo$stats)=="B -> A")

  leo.nb.for[k] = anac(z$neo$stats$LEO.NB.FOR[b.a.row])
  leo.nb.max[k] = anac(z$neo$stats$LEO.NB.MAX[b.a.row])
  leo.nb.all[k] = anac(z$neo$stats$LEO.NB.ALL[b.a.row])


  zeo.for[k] = anac(z$neo$stats$zeo.for[b.a.row])
  zeo.max[k] = anac(z$neo$stats$zeo.max[b.a.row])
  zeo.all[k] = anac(z$neo$stats$zeo.all[b.a.row])
  
  leo.i.for[k] = anac(z$neo$stats$LEO.I.FOR[b.a.row])
  leo.i.max[k] = anac(z$neo$stats$LEO.I.MAX[b.a.row])
  leo.i.all[k] = anac(z$neo$stats$LEO.I.ALL[b.a.row])

  leo.o.for[k] = anac(z$neo$stats$LEO.O.FOR[b.a.row])
  leo.o.max[k] = anac(z$neo$stats$LEO.O.MAX[b.a.row])
  leo.o.all[k] = anac(z$neo$stats$LEO.O.ALL[b.a.row])



 #param.keep[k,]=as.matrix(param[i,]) # track params too


# BLV.AtoB[k]=z$fit2m$zeo.MA.B.given.A$BLV
# BLV.BtoA[k]=z$fit2m$zeo.MB.A.given.B$BLV

print(paste(date(),"    j=",j))

}} # end i,j


  ###### old stuff:

  #  param=generate.simulation.params.old(n.sample.points = 2)
  if (is.null(param)) { 
      param=generate.simulation.params();  
  } else {
      if (is.null(colnames(param))) { 
         colnames(param)=c("ha","hb","gamma2","omega","eA","eB","rho.ab","percent.max.omega")
         param=as.data.frame(param)
      }
  }

  N=no.samples
  set.seed(1)
  num.param=nrow(param)

  total.rows = num.reps*num.param

  param.keep=matrix(NA,nrow=total.rows,ncol=ncol(param))
  colnames(param.keep)=c("ha","hb","gamma2","omega","eA","eB","rho.ab","percent.max.omega")

 BLV.AtoB=rep(NA,total.rows)
 BLV.BtoA=rep(NA,total.rows)

 cor.ma.b=rep(NA,total.rows)
 pcor.ma.b.given.a=rep(NA,total.rows)
 cor.mb.a=rep(NA,total.rows)
 pcor.mb.a.given.b=rep(NA,total.rows)

 path.abs.z.f1ab.AtoB=rep(NA,total.rows)
 path.abs.z.f1ab.BtoA=rep(NA,total.rows)
 path.abs.z.f1ab.AcollideB=rep(NA,total.rows)
 path.abs.z.f1ab.BcollideA=rep(NA,total.rows)

 path.abs.z.f1ba.AtoB=rep(NA,total.rows)
 path.abs.z.f1ba.BtoA=rep(NA,total.rows)
 path.abs.z.f1ba.AcollideB=rep(NA,total.rows)
 path.abs.z.f1ba.BcollideA=rep(NA,total.rows)

 path.abs.z.f2.M.M1M2.AtoB=rep(NA,total.rows)
 path.abs.z.f2.M.M1M2.BtoA=rep(NA,total.rows)

 rmsea.f1ab.AtoB=rep(NA,total.rows)
 rmsea.f1ab.BtoA=rep(NA,total.rows)
 rmsea.f1ab.conf=rep(NA,total.rows)
 rmsea.f1ab.AcollideB=rep(NA,total.rows)
 rmsea.f1ab.BcollideA=rep(NA,total.rows)

 rmsea.f1ba.AtoB=rep(NA,total.rows)
 rmsea.f1ba.BtoA=rep(NA,total.rows)
 rmsea.f1ba.conf=rep(NA,total.rows)
 rmsea.f1ba.AcollideB=rep(NA,total.rows)
 rmsea.f1ba.BcollideA=rep(NA,total.rows)

 rmsea.f2.M.M1M2.AtoB=rep(NA,total.rows)
 rmsea.f2.M.M1M2.BtoA=rep(NA,total.rows)
 rmsea.f2.M.M1M2unresolv=rep(NA,total.rows)
 rmsea.f2.M.M1M2hidden.con=rep(NA,total.rows)

 rmsea90up.f1ab.AtoB=rep(NA,total.rows)
 rmsea90up.f1ab.BtoA=rep(NA,total.rows)
 rmsea90up.f1ab.conf=rep(NA,total.rows)
 rmsea90up.f1ab.AcollideB=rep(NA,total.rows)
 rmsea90up.f1ab.BcollideA=rep(NA,total.rows)

 rmsea90up.f1ba.AtoB=rep(NA,total.rows)
 rmsea90up.f1ba.BtoA=rep(NA,total.rows)
 rmsea90up.f1ba.conf=rep(NA,total.rows)
 rmsea90up.f1ba.AcollideB=rep(NA,total.rows)
 rmsea90up.f1ba.BcollideA=rep(NA,total.rows)

 rmsea90up.f2.M.M1M2.AtoB=rep(NA,total.rows)
 rmsea90up.f2.M.M1M2.BtoA=rep(NA,total.rows)
 rmsea90up.f2.M.M1M2unresolv=rep(NA,total.rows)
 rmsea90up.f2.M.M1M2hidden.con=rep(NA,total.rows)
  
 chisq.f1ab.AtoB=rep(NA,total.rows)
 chisq.f1ab.BtoA=rep(NA,total.rows)
 chisq.f1ab.conf=rep(NA,total.rows)
 chisq.f1ab.AcollideB=rep(NA,total.rows)
 chisq.f1ab.BcollideA=rep(NA,total.rows)

 chisq.f1ba.AtoB=rep(NA,total.rows)
 chisq.f1ba.BtoA=rep(NA,total.rows)
 chisq.f1ba.conf=rep(NA,total.rows)
 chisq.f1ba.AcollideB=rep(NA,total.rows)
 chisq.f1ba.BcollideA=rep(NA,total.rows)

 chisq.f2.M.M1M2.AtoB=rep(NA,total.rows)
 chisq.f2.M.M1M2.BtoA=rep(NA,total.rows)
 chisq.f2.M.M1M2unresolv=rep(NA,total.rows)
 chisq.f2.M.M1M2hidden.con=rep(NA,total.rows)
  
 chiprob.f1ab.AtoB=rep(NA,total.rows)
 chiprob.f1ab.BtoA=rep(NA,total.rows)
 chiprob.f1ab.conf=rep(NA,total.rows)
 chiprob.f1ab.AcollideB=rep(NA,total.rows)
 chiprob.f1ab.BcollideA=rep(NA,total.rows)

 chiprob.f1ba.AtoB=rep(NA,total.rows)
 chiprob.f1ba.BtoA=rep(NA,total.rows)
 chiprob.f1ba.conf=rep(NA,total.rows)
 chiprob.f1ba.AcollideB=rep(NA,total.rows)
 chiprob.f1ba.BcollideA=rep(NA,total.rows)

 chiprob.f2.M.M1M2.AtoB=rep(NA,total.rows)
 chiprob.f2.M.M1M2.BtoA=rep(NA,total.rows)
 chiprob.f2.M.M1M2unresolv=rep(NA,total.rows)
 chiprob.f2.M.M1M2hidden.con=rep(NA,total.rows)
  
  res.list=vector("list", total.rows)
  
 for (i in 1:num.param) {
   for (j in 1:num.reps) {
      k=(num.reps*(i-1))+j # where to store
#  for (i in 1:1) {
#    for (j in 1:2) {
   # q=data.frame(ha,hb,gamma2,omega,eA,eB,rho.ab)

     z =list()

     z$ha =   ha     = param$ha[i]
     z$hb =   hb     = param$hb[i]
     z$gamma2 =   gamma2 = param$gamma2[i]
     z$omega =   omega  = param$omega[i]
     z$eA =   eA =     param$eA[i]
     z$eB =   eB =     param$eB[i]
     z$rho.ab =   rho.ab = param$rho.ab[i]

     param.string=paste(sep="","ha=",ha, " hb=",hb, " g2=",gamma2, " eA=",eA," eB=",eB, " rho.ab=",rho.ab, " w=",omega)
   
#              varB= gamma2 + eB + hb;
#              varA = TA + 2*omega*gamma*gamma + gamma*gamma + omega*omega*varB + eA;

              # if B is standardized to var(B)=1, then....

#  std.varA = TA + gamma^2 + omega^2 + 2*gamma*gamma*omega + eA; # set == 1
#  covAB = omega*TB+gamma*gamma*(1+omega)+omega*eB
#  covAC = gamma+ omega*gamma
#  covBC = gamma

#  cov.SNPB.A = sqrt(TB)*omega

#  rho.AB = gamma^2 + omega

  SNPA=rnorm(no.samples);
  SNPB=rnorm(no.samples);

  C=rnorm(no.samples);

  errA=rnorm(no.samples,sd=sqrt(eA));
  errB=rnorm(no.samples,sd=sqrt(eB));

  B=errB+sqrt(gamma2)*C+ sqrt(hb)*SNPB;

  A=errA+omega*B+sqrt(gamma2)*C+sqrt(ha)*SNPA;

   
  x=data.frame(SNPA,SNPB,A,B);
   
    z$x=x
#        z$fit1mb = compare.local.sems(pm=pm,M.col=2,A.col=3,B.col=4,x=x,no.obs=no.samples)
#        z$fit1ma = compare.local.sems(pm=pm,M.col=1,A.col=3,B.col=4,x=x,no.obs=no.samples)
#        z$fit2m =local.sem.four.var.m1m2(pm=pm,MA.col=1,A.col=3,B.col=4,x=x,MB.col=2,no.obs=N,fit.models=c(1,2,3,4));
    z$param.string = param.string;

 param.keep[k,]=as.matrix(param[i,]) # track params too


 BLV.AtoB[k]=z$fit2m$zeo.MA.B.given.A$BLV
 BLV.BtoA[k]=z$fit2m$zeo.MB.A.given.B$BLV

 cor.ma.b[k]         =z$fit2m$zeo.MA.B.given.A$cor
 pcor.ma.b.given.a[k]=z$fit2m$zeo.MA.B.given.A$pcor

 cor.mb.a[k]         =z$fit2m$zeo.MB.A.given.B$cor
 pcor.mb.a.given.b[k]=z$fit2m$zeo.MB.A.given.B$pcor


 



 path.abs.z.f1ab.AtoB[k]=abs(z$fit1ma$M.AtoB$coef["fr.A.B",3])
 path.abs.z.f1ab.BtoA[k]=abs(z$fit1ma$M.BtoA$coef["fr.B.A",3])
 path.abs.z.f1ab.AcollideB[k]=abs(z$fit1ma$M.AcollideB$coef["fr.B.A",3])
 path.abs.z.f1ab.BcollideA[k]=abs(z$fit1ma$M.BcollideA$coef["fr.A.B",3])

 path.abs.z.f1ba.AtoB[k]=abs(z$fit1mb$M.AtoB$coef["fr.A.B",3])
 path.abs.z.f1ba.BtoA[k]=abs(z$fit1mb$M.BtoA$coef["fr.B.A",3])
 path.abs.z.f1ba.AcollideB[k]=abs(z$fit1mb$M.BcollideA$coef["fr.A.B",3])
 path.abs.z.f1ba.BcollideA[k]=abs(z$fit1mb$M.AcollideB$coef["fr.B.A",3])

 path.abs.z.f2.M.M1M2.AtoB[k]=abs(z$fit2m$M.M1M2.AtoB$coef["fr.A.B",3])
 path.abs.z.f2.M.M1M2.BtoA[k]=abs(z$fit2m$M.M1M2.BtoA$coef["fr.B.A",3])

 rmsea.f1ab.AtoB[k]=z$fit1ma$M.AtoB$RMSEA[1]
 rmsea.f1ab.BtoA[k]=z$fit1ma$M.BtoA$RMSEA[1]
 rmsea.f1ab.conf[k]=z$fit1ma$M.conf$RMSEA[1]
 rmsea.f1ab.AcollideB[k]=z$fit1ma$M.AcollideB$RMSEA[1]
 rmsea.f1ab.BcollideA[k]=z$fit1ma$M.BcollideA$RMSEA[1]

 rmsea.f1ba.AtoB[k]=z$fit1mb$M.AtoB$RMSEA[1]
 rmsea.f1ba.BtoA[k]=z$fit1mb$M.BtoA$RMSEA[1]
 rmsea.f1ba.conf[k]=z$fit1mb$M.conf$RMSEA[1]
 rmsea.f1ba.AcollideB[k]=z$fit1mb$M.AcollideB$RMSEA[1]
 rmsea.f1ba.BcollideA[k]=z$fit1mb$M.BcollideA$RMSEA[1]

 rmsea.f2.M.M1M2.AtoB[k]=z$fit2m$M.M1M2.AtoB$RMSEA[1]
 rmsea.f2.M.M1M2.BtoA[k]=z$fit2m$M.M1M2.BtoA$RMSEA[1]
 rmsea.f2.M.M1M2unresolv[k]=z$fit2m$M.M1M2unresolv$RMSEA[1]
 rmsea.f2.M.M1M2hidden.con[k]=z$fit2m$M.M1M2hidden.con$RMSEA[1]

 rmsea90up.f1ab.AtoB[k]=z$fit1ma$M.AtoB$RMSEA[3]
 rmsea90up.f1ab.BtoA[k]=z$fit1ma$M.BtoA$RMSEA[3]
 rmsea90up.f1ab.conf[k]=z$fit1ma$M.conf$RMSEA[3]
 rmsea90up.f1ab.AcollideB[k]=z$fit1ma$M.AcollideB$RMSEA[3]
 rmsea90up.f1ab.BcollideA[k]=z$fit1ma$M.BcollideA$RMSEA[3]

 rmsea90up.f1ba.AtoB[k]=z$fit1mb$M.AtoB$RMSEA[3]
 rmsea90up.f1ba.BtoA[k]=z$fit1mb$M.BtoA$RMSEA[3]
 rmsea90up.f1ba.conf[k]=z$fit1mb$M.conf$RMSEA[3]
 rmsea90up.f1ba.AcollideB[k]=z$fit1mb$M.AcollideB$RMSEA[3]
 rmsea90up.f1ba.BcollideA[k]=z$fit1mb$M.BcollideA$RMSEA[3]

 rmsea90up.f2.M.M1M2.AtoB[k]=z$fit2m$M.M1M2.AtoB$RMSEA[3]
 rmsea90up.f2.M.M1M2.BtoA[k]=z$fit2m$M.M1M2.BtoA$RMSEA[3]
 rmsea90up.f2.M.M1M2unresolv[k]=z$fit2m$M.M1M2unresolv$RMSEA[3]
 rmsea90up.f2.M.M1M2hidden.con[k]=z$fit2m$M.M1M2hidden.con$RMSEA[3]
  
 chisq.f1ab.AtoB[k]=z$fit1ma$M.AtoB$chisq
 chisq.f1ab.BtoA[k]=z$fit1ma$M.BtoA$chisq
 chisq.f1ab.conf[k]=z$fit1ma$M.conf$chisq
 chisq.f1ab.AcollideB[k]=z$fit1ma$M.AcollideB$chisq
 chisq.f1ab.BcollideA[k]=z$fit1ma$M.BcollideA$chisq

 chisq.f1ba.AtoB[k]=z$fit1mb$M.AtoB$chisq
 chisq.f1ba.BtoA[k]=z$fit1mb$M.BtoA$chisq
 chisq.f1ba.conf[k]=z$fit1mb$M.conf$chisq
 chisq.f1ba.AcollideB[k]=z$fit1mb$M.AcollideB$chisq
 chisq.f1ba.BcollideA[k]=z$fit1mb$M.BcollideA$chisq

 chisq.f2.M.M1M2.AtoB[k]=z$fit2m$M.M1M2.AtoB$chisq
 chisq.f2.M.M1M2.BtoA[k]=z$fit2m$M.M1M2.BtoA$chisq
 chisq.f2.M.M1M2unresolv[k]=z$fit2m$M.M1M2unresolv$chisq
 chisq.f2.M.M1M2hidden.con[k]=z$fit2m$M.M1M2hidden.con$chisq

 chiprob.f1ab.AtoB[k]=extract.chisq.right.log10prob(z$fit1ma$M.AtoB)
 chiprob.f1ab.BtoA[k]=extract.chisq.right.log10prob(z$fit1ma$M.BtoA)
 chiprob.f1ab.conf[k]=extract.chisq.right.log10prob(z$fit1ma$M.conf)
 chiprob.f1ab.AcollideB[k]=extract.chisq.right.log10prob(z$fit1ma$M.AcollideB)
 chiprob.f1ab.BcollideA[k]=extract.chisq.right.log10prob(z$fit1ma$M.BcollideA)

 chiprob.f1ba.AtoB[k]=extract.chisq.right.log10prob(z$fit1mb$M.AtoB)
 chiprob.f1ba.BtoA[k]=extract.chisq.right.log10prob(z$fit1mb$M.BtoA)
 chiprob.f1ba.conf[k]=extract.chisq.right.log10prob(z$fit1mb$M.conf)
 chiprob.f1ba.AcollideB[k]=extract.chisq.right.log10prob(z$fit1mb$M.AcollideB)
 chiprob.f1ba.BcollideA[k]=extract.chisq.right.log10prob(z$fit1mb$M.BcollideA)

 chiprob.f2.M.M1M2.AtoB[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2.AtoB)
 chiprob.f2.M.M1M2.BtoA[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2.BtoA)
 chiprob.f2.M.M1M2unresolv[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2unresolv)
 chiprob.f2.M.M1M2hidden.con[k]=extract.chisq.right.log10prob(z$fit2m$M.M1M2hidden.con)
     
    res.list[[(i-1)*num.reps + j]]=z

   } # end j
#    save(res.list,file="my.res.list.rdat")

 } # end i over params

 scores=data.frame(param.keep,rmsea.f1ab.AtoB, rmsea.f1ab.BtoA, 
 rmsea.f1ab.conf, rmsea.f1ab.AcollideB, rmsea.f1ab.BcollideA,
 rmsea.f1ba.AtoB, rmsea.f1ba.BtoA, rmsea.f1ba.conf, rmsea.f1ba.AcollideB, 
 rmsea.f1ba.BcollideA, rmsea.f2.M.M1M2.AtoB, rmsea.f2.M.M1M2.BtoA, 
 rmsea.f2.M.M1M2unresolv, rmsea.f2.M.M1M2hidden.con, rmsea90up.f1ab.AtoB, 
 rmsea90up.f1ab.BtoA, rmsea90up.f1ab.conf, rmsea90up.f1ab.AcollideB, rmsea90up.f1ab.BcollideA,
 rmsea90up.f1ba.AtoB, rmsea90up.f1ba.BtoA, rmsea90up.f1ba.conf, rmsea90up.f1ba.AcollideB, rmsea90up.f1ba.BcollideA,
 rmsea90up.f2.M.M1M2.AtoB, rmsea90up.f2.M.M1M2.BtoA, rmsea90up.f2.M.M1M2unresolv, rmsea90up.f2.M.M1M2hidden.con,  
 chisq.f1ab.AtoB, chisq.f1ab.BtoA, chisq.f1ab.conf, chisq.f1ab.AcollideB, chisq.f1ab.BcollideA,
 chisq.f1ba.AtoB, chisq.f1ba.BtoA, chisq.f1ba.conf, chisq.f1ba.AcollideB, chisq.f1ba.BcollideA,
 chisq.f2.M.M1M2.AtoB, chisq.f2.M.M1M2.BtoA, chisq.f2.M.M1M2unresolv, chisq.f2.M.M1M2hidden.con,  
 chiprob.f1ab.AtoB, chiprob.f1ab.BtoA, chiprob.f1ab.conf, chiprob.f1ab.AcollideB, chiprob.f1ab.BcollideA,
 chiprob.f1ba.AtoB, chiprob.f1ba.BtoA, chiprob.f1ba.conf, chiprob.f1ba.AcollideB, chiprob.f1ba.BcollideA,
 chiprob.f2.M.M1M2.AtoB, chiprob.f2.M.M1M2.BtoA, chiprob.f2.M.M1M2unresolv, chiprob.f2.M.M1M2hidden.con,
 path.abs.z.f1ab.AtoB, path.abs.z.f1ab.BtoA,path.abs.z.f1ab.AcollideB,
 path.abs.z.f1ab.BcollideA,path.abs.z.f1ba.AtoB,path.abs.z.f1ba.BtoA,
 path.abs.z.f1ba.AcollideB,path.abs.z.f1ba.BcollideA,path.abs.z.f2.M.M1M2.AtoB,
 path.abs.z.f2.M.M1M2.BtoA, 
 BLV.AtoB, BLV.BtoA, 
 cor.ma.b, pcor.ma.b.given.a, cor.mb.a, pcor.mb.a.given.b)

 save(res.list,scores,file="my.new.res.list.scores.rdat")

 #### RETURN EARLY
 return(scores);

 # analysis before completion...
 # scores.new = scores ... b/c we loaded it separately.

# examples:
#     boxplot(split(log10(d$leo.i.correct[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",at = 1:length(unique(d$ntrack[all])) - 0.2,boxwex = 0.25,col = "yellow")#,ylim=c(-20,20))
#   boxplot(split(log10(d$leo.i[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",at = 1:length(unique(d$ntrack[all])) + 0.2,add=TRUE,boxwex = 0.25,col = "orange",show.names=FALSE)#,ylim=c(-20,20))

 file.out="power.study.rmsea.f2.AtoB.over5600.ps"
 postscript(file=file.out,horizontal=FALSE);
 boxplot(split(scores.new$rmsea.f2.M.M1M2.BtoA,scores.new$percent.max.omega),notch=TRUE,xlab="percent omega.max",ylab="RMSEA",col="yellow",at=1:length(unique(scores.new$percent.max.omega))-.2,boxwex=.25,ylim=c(-.15,.7))
 boxplot(split(scores.new$rmsea.f2.M.M1M2.AtoB,scores.new$percent.max.omega),notch=TRUE,col="violetred",at=1:length(unique(scores.new$percent.max.omega))+.2,boxwex=.25,add=TRUE,show.names=FALSE,main="Monte Carlo study: B causal for A" )
  legend(5.1,-.05, c("correct causal B to A", "incorrect A to B model"), fill = c("yellow", "violetred"))
  dev.off()

# things that vary
  # hb = .2, .333, .467, .6
  #
  # at hb = .2
  #    omega = 0, .117, .233, .350, .583, .700

  # start by studying the null model: no causal, no confounding
  # under different levels of genetic influence of HB

  no.gamma.no.omega=c(1:50,1401:1450,2801:2850,4201:4250)
  s=scores.new[no.gamma.no.omega,]

  file.out="power.study.rmsea.no.gamma.no.omega.hb.varies.ps"
  postscript(file=file.out,horizontal=FALSE);
  boxplot(split(s$rmsea.f2.M.M1M2.BtoA,s$hb),notch=TRUE,xlab="percent omega.max",ylab="RMSEA",col="yellow",at=1:length(unique(s$hb))-.2,boxwex=.25,log="y")#,ylim=c(-.15,.7))
  boxplot(split(s$rmsea.f2.M.M1M2.AtoB,s$percent.max.omega),notch=TRUE,col="violetred",at=1:length(unique(s$percent.max.omega))+.2,boxwex=.25,add=TRUE,show.names=FALSE,main="Monte Carlo study: B causal for A" )
  legend(5.1,-.05, c("correct causal B to A", "incorrect A to B model"), fill = c("yellow", "violetred"))
  dev.off()


  # make barplots with error bars
  # ex:
  # freq.bar(s,threshold=.1,th.col="rmsea.f2.M.M1M2.BtoA",group.col="percent.max.omega")

  
  
 rownames(scores) = as.character(1:nrow(scores))
 scores  
} # noise.snp.sim




if(exists("compare.for.max.all") ) rm(compare.for.max.all);
compare.for.max.all=function(param=NULL, no.samples=200, num.reps=50,pm=neo.get.param()) {


ha=hb=.6

rho.ab=.5
percent.max.omega=1; # .6667
omega=percent.max.omega*rho.ab
gamma2=rho.ab-omega

gamma=sqrt(gamma2)

eA=1-2*omega*gamma^2 - gamma^2 - omega^2 - ha

eB=1-gamma2 - hb;


  q=data.frame(ha,hb,gamma2,omega,eA,eB,rho.ab,percent.max.omega)

# just for development DEBUG, take out later, so function works...
param=q
no.samples=200
pm=neo.get.param()
pm$skip.snp.select = TRUE
pm$quiet = TRUE
pm$run.title="trashme"
pm$cor.ind.th = .3 # so max doesn't exit away?
pm$zeo.proxy.for.leo.th = NA
pm$no.log = TRUE

  num.reps=1
  num.param=1
  total.rows = num.reps*num.param

  param.keep=matrix(NA,nrow=total.rows,ncol=ncol(param))
  colnames(param.keep)=c("ha","hb","gamma2","omega","eA","eB","rho.ab","percent.max.omega")


  BLV.AtoB=rep(NA,total.rows)
  BLV.BtoA=rep(NA,total.rows)

  leo.nb.for = rep(NA,total.rows)
  leo.i.for  = rep(NA,total.rows)
  leo.o.for  = rep(NA,total.rows)

  zeo.for = rep(NA,total.rows)
  zeo.max  = rep(NA,total.rows)
  zeo.all  = rep(NA,total.rows)

  leo.nb.max = rep(NA,total.rows)
  leo.i.max  = rep(NA,total.rows)
  leo.o.max  = rep(NA,total.rows)

  zeo.nb.max = rep(NA,total.rows)
  zeo.i.max  = rep(NA,total.rows)
  zeo.o.max  = rep(NA,total.rows)

  leo.nb.all = rep(NA,total.rows)
  leo.i.all  = rep(NA,total.rows)
  leo.o.all  = rep(NA,total.rows)

  zeo.nb.all = rep(NA,total.rows)
  zeo.i.all  = rep(NA,total.rows)
  zeo.o.all  = rep(NA,total.rows)

  # needed for anac() to not round too much....
  options(digits = 4); # affects the conversion of spreadsheet columns from factors to character to numeric, at the character to numeric step.

  
 for (i in 1:num.param) {
   for (j in 1:num.reps) {
      k=(num.reps*(i-1))+j # where to store

      z =list()

     z$ha =   ha     = param$ha[i]
     z$hb =   hb     = param$hb[i]
     z$gamma2 =   gamma2 = param$gamma2[i]
     z$omega =   omega  = param$omega[i]
     z$eA =   eA =     param$eA[i]
     z$eB =   eB =     param$eB[i]
     z$rho.ab =   rho.ab = param$rho.ab[i]

     param.string=paste(sep="","ha=",ha, " hb=",hb, " g2=",gamma2, " eA=",eA," eB=",eB, " rho.ab=",rho.ab, " w=",omega)
   
#              varB= gamma2 + eB + hb;
#              varA = TA + 2*omega*gamma*gamma + gamma*gamma + omega*omega*varB + eA;

              # if B is standardized to var(B)=1, then....

#  std.varA = TA + gamma^2 + omega^2 + 2*gamma*gamma*omega + eA; # set == 1
#  covAB = omega*TB+gamma*gamma*(1+omega)+omega*eB
#  covAC = gamma+ omega*gamma
#  covBC = gamma

#  cov.SNPB.A = sqrt(TB)*omega

#  rho.AB = gamma^2 + omega


    s2=sqrt(2)
    SNPA1=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
    SNPA2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPA3=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPA4=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPA5=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPA6=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPA7=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPA8=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPA9=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPA10=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1

##   SNPA = scale(SNPA1+SNPA2+SNPA3+SNPA4+SNPA5+SNPA6+SNPA7+SNPA8+SNPA9+SNPA10)
   SNPA = scale(SNPA1)

## #   MB1=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
## #   MB2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
## #   MB2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
## #   SNPB = scale(MB1+MB2)


    SNPB1=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
    SNPB2=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPB3=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPB4=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPB5=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPB6=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPB7=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPB8=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPB9=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1
##    SNPB10=s2*sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T ); # variance=1

##   SNPB = scale(SNPB1+SNPB2+SNPB3+SNPB4+SNPB5+SNPB6+SNPB7+SNPB8+SNPB9+SNPB10)

   SNPB = scale(.8*SNPB1+.2*SNPB2)


# GAUSSIAN SNPA and SNPB for the moment

#  SNPA=rnorm(no.samples);
#  SNPB=rnorm(no.samples);


  C=rnorm(no.samples);

  errA=rnorm(no.samples,sd=sqrt(eA));
  errB=rnorm(no.samples,sd=sqrt(eB));

  B=errB+sqrt(gamma2)*C+ sqrt(hb)*SNPB;

  A=errA+omega*B+sqrt(gamma2)*C+sqrt(ha)*SNPA;


#     x=data.frame(SNPA1,SNPA2,SNPA3,SNPA4,SNPA5,SNPA6,SNPA7,SNPA8,SNPA9,SNPA10,SNPB,A,B);
     x=data.frame(SNPA1,SNPB1,SNPB2,A,B);
#  x=data.frame(SNPA,SNPB,A,B);

  z$neo =neo(x,pm=pm,traitcols=4)
#    z$x=x
##        z$fit1mb = compare.local.sems(pm=pm,M.col=2,A.col=3,B.col=4,x=x,no.obs=no.samples)
#        z$fit1ma = compare.local.sems(pm=pm,M.col=1,A.col=3,B.col=4,x=x,no.obs=no.samples)
#        z$fit2m =local.sem.four.var.m1m2(pm=pm,MA.col=1:10,A.col=12,B.col=13,x=x,MB.col=11,no.obs=N,fit.models=c(1,2,3,4));
#    z$param.string = param.string;


  b.a.row=which(rownames(z$neo$stats)=="B -> A")
  a.b.row=which(rownames(z$neo$stats)=="A -> B")

  leo.nb.for[k] = anac(z$neo$stats$LEO.NB.FOR[b.a.row])
  leo.nb.max[k] = anac(z$neo$stats$LEO.NB.MAX[b.a.row])
  leo.nb.all[k] = anac(z$neo$stats$LEO.NB.ALL[b.a.row])


  zeo.for[k] = anac(z$neo$stats$zeo.for[b.a.row])
  zeo.max[k] = anac(z$neo$stats$zeo.max[b.a.row])
  zeo.all[k] = anac(z$neo$stats$zeo.all[b.a.row])
  
  leo.i.for[k] = anac(z$neo$stats$LEO.I.FOR[b.a.row])
  leo.i.max[k] = anac(z$neo$stats$LEO.I.MAX[b.a.row])
  leo.i.all[k] = anac(z$neo$stats$LEO.I.ALL[b.a.row])

  leo.o.for[k] = anac(z$neo$stats$LEO.O.FOR[b.a.row])
  leo.o.max[k] = anac(z$neo$stats$LEO.O.MAX[b.a.row])
  leo.o.all[k] = anac(z$neo$stats$LEO.O.ALL[b.a.row])

 param.keep[k,]=as.matrix(param[i,]) # track params too

 BLV.AtoB[k]=anac(z$neo$stats$BilayerZscore[a.b.row])
 BLV.BtoA[k]=anac(z$neo$stats$BilayerZscore[b.a.row])

  print(paste(date(),"    j=",j))

}} # end i,j

  r=list()
  r$param = param
  r$leo.nb.for = leo.nb.for
  r$leo.nb.max = leo.nb.max
  r$leo.nb.all = leo.nb.all
  r$zeo.for = zeo.for
  r$zeo.max = zeo.max
  r$zeo.all = zeo.all
  r$leo.i.for = leo.i.for
  r$leo.i.max = leo.i.max
  r$leo.i.all = leo.i.all
  r$leo.o.for = leo.o.for
  r$leo.o.max = leo.o.max
  r$leo.o.all = leo.o.all
  r$BLV.AtoB = BLV.AtoB
  r$BLV.BtoA = BLV.BtoA

  r$num.reps = num.reps
  r$num.param = num.param
  r$total.rows = total.rows
  r$param.keep = param.keep

  r

} # compare.for.max.all



## [R] the chi-square test for trend

##     * This message: [ Message body ] [ More options ]
##     * Related messages: [ Next message ] [ Previous message ]

## From: Eric Lecoutre <lecoutre>
## Date: Thu Apr 10 11:26:08 2003

## Hi,

## Please find here a function to compute Mantel-Haenszel Chi squared as well
## as Cochran-Armitage test for trend (suitable when one of the dimensions is 2).
## Both functions take as argument a contingency table.
## By the way, the MH Chi? is very easy to compute as it is defined as n*rho.

## Eric Lecoutre

## ------------------------------------------------------------------

tablepearson=function(x,scores.type="table")
{

# Statistic
        sR=scores(x,1,scores.type)
        sC=scores(x,2,scores.type)
        n=sum(data)
        Rbar=sum(apply(x,1,sum)*sR)/n
        Cbar=sum(apply(x,2,sum)*sC)/n
        ssr=sum(x*(sR-Rbar)^2)
        ssc=sum(t(x)* (sC-Cbar)^2)
        tmpij=outer(sR,sC,FUN=function(a,b) return((a-Rbar)*(b-Cbar)))
        ssrc= sum(x*tmpij)
        v=ssrc
        w=sqrt(ssr*ssc)
        r=v/w
# ASE
        bij=outer(sR,sC, FUN=function(a,b)return((a-Rbar)^2*ssc + (b-Cbar)^2*ssr))
        tmp1=1/w^2
        tmp2=x*(w*tmpij - (bij*v)/(2*w))^2
        tmp3=sum(tmp2)
        ASE=tmp1*sqrt(tmp3)
# Test
        var0= (sum(x*tmpij) - (ssrc^2/n))/ (ssr*ssc)
        tb=r/sqrt(var0)
        p.value=2*(1-pnorm(tb))
# Output
        out=list(estimate=r,ASE=ASE,statistic=tb,p.value=p.value,name="Pearson Correlation",bornes=c(-1,1))
        class(out)="ordtest"
        return(out)
}

tableChisqMH=function(x)
{
        n=sum(x)
        G2=n*(tablepearson(x)^2)
        dll=1
        p.value=1-pchisq(G2,dll)
        out=list(estimate=G2,dll=dll,p.value=p.value,dim=dim(x),name="Mantel-Haenszel Chi-square")
        return(out)

}

tabletrend=function(x,transpose=FALSE)
{
        if (any(dim(x)==2))
        {
        if (transpose==TRUE) {
        x=t(x)
        }
        
        if (dim(x)[2]!=2){stop("Cochran-Armitage test for trend must be used with a (R,2) table. Use transpose argument",call.=FALSE) }
        
        nidot=apply(x,1,sum)
        n=sum(nidot)

        Ri=scores(x,1,"table")
        Rbar=sum(nidot*Ri)/n
        
        s2=sum(nidot*(Ri-Rbar)^2)
        pdot1=sum(x[,1])/n
        T=sum(x[,1]*(Ri-Rbar))/sqrt(pdot1*(1-pdot1)*s2)
        p.value.uni=1-pnorm(abs(T))
        p.value.bi=2*p.value.uni
        out=list(estimate=T,dim=dim(x),p.value.uni=p.value.uni,p.value.bi=p.value.bi,name="Cochran-Armitage test for trend")
        return(out)
        
        }
        else {stop("Cochran-Armitage test for trend must be used with a (2,C) or a (R,2) table",call.=FALSE) }
} 


test.robust=function() {

  num.reps=100

  sge.task = get.sge.task()

  set.seed(554324)
  r=list()

  r$forward.null = list()
  r$greedy.null = list()
  r$greedy.and.forward.null = list()

  r$forward.50.causal = list()
  r$greedy.50.causal = list()
  r$greedy.and.forward.50.causal = list()

  r$forward.100.causal = list()
  r$greedy.100.causal = list()
  r$greedy.and.forward.100.causal = list()

  for (k in 1:num.reps) {

   restricted.heritability.a = .4
   restricted.heritability.b = .4
   observed.a.b.correlation  = .3
   topK = 20

#     dat=neo.monte.carlo.data(do.a.snp=1,do.b.snp=1, percent.causal.signal=0,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=4, nearby.b.snps=4, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))
#
#     b=ncol(dat)
#     a=b-1
#
#     datSNP =dat[,-c(a,b)]
#     datExpr=dat[, c(a,b)]


##      r$greedy.null[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
##      r$forward.null[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
##      r$greedy.and.forward.null[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)
##
##      dat=neo.monte.carlo.data(do.a.snp=4,do.b.snp=4, percent.causal.signal=0.5,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=0, nearby.b.snps=0, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

##      datSNP =dat[,-c(a,b)]
##      datExpr=dat[, c(a,b)]
##
##      r$greedy.50.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
##      r$forward.50.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
##      r$greedy.and.forward.50.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)

      dat=neo.monte.carlo.data(do.a.snp=1,do.b.snp=1, percent.causal.signal=1,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=4, nearby.b.snps=4, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

      b=ncol(dat)
      a=b-1

     datSNP =dat[,-c(a,b)]
     datExpr=dat[, c(a,b)]

#if (sge.task==1)     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
#if (sge.task==2)     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
#if (sge.task==3)     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)

    r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
    r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
    r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)

# just do topK snps...
#     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=topK, TopForward=1, no.log.quiet=TRUE)
#     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=1, TopForward=topK, no.log.quiet=TRUE)
#     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=topK, TopGreedy=topK, no.log.quiet=TRUE)

     print(paste(sep="","Done with k=",k," of out ",num.reps," at ",date()))

     sge.task = get.sge.task()

     the.file=gsub(":",".",gsub(" ","_",paste(sep="","r.test.robust.task",sge.task,".results.",k,".",date(),".rdat")))
     save(r,file=the.file)

 } # end for num.reps


  r
}



test.robust.1a.1b=function() {

  num.reps=100

  sge.task = get.sge.task()

  set.seed(sge.task * 101)
  r=list()

  r$forward.null = list()
  r$greedy.null = list()
  r$greedy.and.forward.null = list()

  r$forward.50.causal = list()
  r$greedy.50.causal = list()
  r$greedy.and.forward.50.causal = list()

  r$forward.100.causal = list()
  r$greedy.100.causal = list()
  r$greedy.and.forward.100.causal = list()

  for (k in 1:num.reps) {

   restricted.heritability.a = .4
   restricted.heritability.b = .4
   observed.a.b.correlation  = .3
   topK = 20

      dat=neo.monte.carlo.data(do.a.snp=1,do.b.snp=1, percent.causal.signal=1,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=0, nearby.b.snps=0, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

      b=ncol(dat)
      a=b-1

     datSNP =dat[,-c(a,b)]
     datExpr=dat[, c(a,b)]

if (sge.task==1)     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
if (sge.task==2)     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
if (sge.task==3)     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)

# just do topK snps...
#     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=topK, TopForward=1, no.log.quiet=TRUE)
#     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=1, TopForward=topK, no.log.quiet=TRUE)
#     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=topK, TopGreedy=topK, no.log.quiet=TRUE)

     print(paste(sep="","Done with k=",k," of out ",num.reps," at ",date()))

     sge.task = get.sge.task()

     the.file=gsub(":",".",gsub(" ","_",paste(sep="","r.test.robust.task",sge.task,".results.",k,".",date(),".rdat")))
     save(r,file=the.file)

 } # end for num.reps


  r
}




test.robust.0a.1b.5causal=function() {

  num.reps=100

  sge.task = get.sge.task()

  set.seed(sge.task * 101)
  r=list()

  r$forward.null = list()
  r$greedy.null = list()
  r$greedy.and.forward.null = list()

  r$forward.50.causal = list()
  r$greedy.50.causal = list()
  r$greedy.and.forward.50.causal = list()

  r$forward.100.causal = list()
  r$greedy.100.causal = list()
  r$greedy.and.forward.100.causal = list()

  for (k in 1:num.reps) {

   restricted.heritability.a = .4
   restricted.heritability.b = .4
   observed.a.b.correlation  = .3
   topK = 20

      dat=neo.monte.carlo.data(do.a.snp=0,do.b.snp=1, percent.causal.signal=.5,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=0, nearby.b.snps=0, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

      b=ncol(dat)
      a=b-1

     datSNP =dat[,-c(a,b)]
     datExpr=dat[, c(a,b)]

if (sge.task==1)     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
if (sge.task==2)     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
if (sge.task==3)     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)

     print(paste(sep="","Done with k=",k," of out ",num.reps," at ",date()))

     sge.task = get.sge.task()

     the.file=gsub(":",".",gsub(" ","_",paste(sep="","r.test.robust.task",sge.task,".results.",k,".",date(),".rdat")))
     save(r,file=the.file)

 } # end for num.reps


  r
}



test.robust.1a.1b.5causal=function() {

  num.reps=100

  sge.task = get.sge.task()

  set.seed(sge.task * 101)
  r=list()

  r$forward.null = list()
  r$greedy.null = list()
  r$greedy.and.forward.null = list()

  r$forward.50.causal = list()
  r$greedy.50.causal = list()
  r$greedy.and.forward.50.causal = list()

  r$forward.100.causal = list()
  r$greedy.100.causal = list()
  r$greedy.and.forward.100.causal = list()

  for (k in 1:num.reps) {

   restricted.heritability.a = .4
   restricted.heritability.b = .4
   observed.a.b.correlation  = .3
   topK = 20

      dat=neo.monte.carlo.data(do.a.snp=1,do.b.snp=1, percent.causal.signal=.5,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=0, nearby.b.snps=0, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

      b=ncol(dat)
      a=b-1

     datSNP =dat[,-c(a,b)]
     datExpr=dat[, c(a,b)]

if (sge.task==1)     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
if (sge.task==2)     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
if (sge.task==3)     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)


     print(paste(sep="","Done with k=",k," of out ",num.reps," at ",date()))

     sge.task = get.sge.task()

     the.file=gsub(":",".",gsub(" ","_",paste(sep="","r.test.robust.task",sge.task,".results.",k,".",date(),".rdat")))
     save(r,file=the.file)

 } # end for num.reps


  r
}



test.robust.1a.1b.0causal=function() {

  num.reps=100

  sge.task = get.sge.task()

  set.seed(sge.task * 101)
  r=list()

  r$forward.null = list()
  r$greedy.null = list()
  r$greedy.and.forward.null = list()

  r$forward.50.causal = list()
  r$greedy.50.causal = list()
  r$greedy.and.forward.50.causal = list()

  r$forward.100.causal = list()
  r$greedy.100.causal = list()
  r$greedy.and.forward.100.causal = list()

  for (k in 1:num.reps) {

   restricted.heritability.a = .4
   restricted.heritability.b = .4
   observed.a.b.correlation  = .3
   topK = 20

      dat=neo.monte.carlo.data(do.a.snp=1,do.b.snp=1, percent.causal.signal=0,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=0, nearby.b.snps=0, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

      b=ncol(dat)
      a=b-1

     datSNP =dat[,-c(a,b)]
     datExpr=dat[, c(a,b)]

if (sge.task==1)     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
if (sge.task==2)     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
if (sge.task==3)     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)


     print(paste(sep="","Done with k=",k," of out ",num.reps," at ",date()))

     sge.task = get.sge.task()

     the.file=gsub(":",".",gsub(" ","_",paste(sep="","r.test.robust.task",sge.task,".results.",k,".",date(),".rdat")))
     save(r,file=the.file)

 } # end for num.reps


  r
}



test.robust.0a.1b.0causal=function() {

  num.reps=100

  sge.task = get.sge.task()

  set.seed(sge.task * 101)
  r=list()

  r$forward.null = list()
  r$greedy.null = list()
  r$greedy.and.forward.null = list()

  r$forward.50.causal = list()
  r$greedy.50.causal = list()
  r$greedy.and.forward.50.causal = list()

  r$forward.100.causal = list()
  r$greedy.100.causal = list()
  r$greedy.and.forward.100.causal = list()

  for (k in 1:num.reps) {

   restricted.heritability.a = .4
   restricted.heritability.b = .4
   observed.a.b.correlation  = .3
   topK = 20

      dat=neo.monte.carlo.data(do.a.snp=0,do.b.snp=1, percent.causal.signal=0,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=0, nearby.b.snps=0, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

      b=ncol(dat)
      a=b-1

     datSNP =dat[,-c(a,b)]
     datExpr=dat[, c(a,b)]

if (sge.task==1)     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
if (sge.task==2)     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
if (sge.task==3)     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)


     print(paste(sep="","Done with k=",k," of out ",num.reps," at ",date()))

     sge.task = get.sge.task()

     the.file=gsub(":",".",gsub(" ","_",paste(sep="","r.test.robust.task",sge.task,".results.",k,".",date(),".rdat")))
     save(r,file=the.file)

 } # end for num.reps


  r
}


test.robust.1a.1b.1causal.4an.4bn=function() {

  num.reps=100

  sge.task = get.sge.task()

  set.seed(sge.task * 101)
  r=list()

  r$forward.null = list()
  r$greedy.null = list()
  r$greedy.and.forward.null = list()

  r$forward.50.causal = list()
  r$greedy.50.causal = list()
  r$greedy.and.forward.50.causal = list()

  r$forward.100.causal = list()
  r$greedy.100.causal = list()
  r$greedy.and.forward.100.causal = list()

  for (k in 1:num.reps) {

   restricted.heritability.a = .4
   restricted.heritability.b = .4
   observed.a.b.correlation  = .3
   topK = 20

      dat=neo.monte.carlo.data(do.a.snp=1,do.b.snp=1, percent.causal.signal=1,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=4, nearby.b.snps=4, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

      b=ncol(dat)
      a=b-1

     datSNP =dat[,-c(a,b)]
     datExpr=dat[, c(a,b)]

if (sge.task==1)     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
if (sge.task==2)     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
if (sge.task==3)     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)


     print(paste(sep="","Done with k=",k," of out ",num.reps," at ",date()))

     sge.task = get.sge.task()

     the.file=gsub(":",".",gsub(" ","_",paste(sep="","r.test.robust.task",sge.task,".results.",k,".",date(),".rdat")))
     save(r,file=the.file)

 } # end for num.reps


  r
}



test.robust.1a.4b.4causal.4an.4bn=function() {

  num.reps=100

  sge.task = get.sge.task()

  set.seed(sge.task * 101)
  r=list()

  r$forward.null = list()
  r$greedy.null = list()
  r$greedy.and.forward.null = list()

  r$forward.50.causal = list()
  r$greedy.50.causal = list()
  r$greedy.and.forward.50.causal = list()

  r$forward.100.causal = list()
  r$greedy.100.causal = list()
  r$greedy.and.forward.100.causal = list()

  for (k in 1:num.reps) {

   restricted.heritability.a = .4
   restricted.heritability.b = .4
   observed.a.b.correlation  = .3
   topK = 20

      dat=neo.monte.carlo.data(do.a.snp=1,do.b.snp=4, percent.causal.signal=1,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=4, nearby.b.snps=4, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

      b=ncol(dat)
      a=b-1

     datSNP =dat[,-c(a,b)]
     datExpr=dat[, c(a,b)]

if (sge.task==1)     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
if (sge.task==2)     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
if (sge.task==3)     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)


     print(paste(sep="","Done with k=",k," of out ",num.reps," at ",date()))

     sge.task = get.sge.task()

     the.file=gsub(":",".",gsub(" ","_",paste(sep="","r.test.robust.task",sge.task,".results.",k,".",date(),".rdat")))
     save(r,file=the.file)

 } # end for num.reps


  r
}




test.robust.0a.4b.4causal.4an.4bn=function() {

  num.reps=100

  sge.task = get.sge.task()

  set.seed(sge.task * 101)
  r=list()

  r$forward.null = list()
  r$greedy.null = list()
  r$greedy.and.forward.null = list()

  r$forward.50.causal = list()
  r$greedy.50.causal = list()
  r$greedy.and.forward.50.causal = list()

  r$forward.100.causal = list()
  r$greedy.100.causal = list()
  r$greedy.and.forward.100.causal = list()

  for (k in 1:num.reps) {

   restricted.heritability.a = .4
   restricted.heritability.b = .4
   observed.a.b.correlation  = .3
   topK = 20

      dat=neo.monte.carlo.data(do.a.snp=0, do.b.snp=4, percent.causal.signal=1,    noise.a.snp=50, noise.b.snp=50, nearby.a.snps=4, nearby.b.snps=4, restricted.heritability.a = restricted.heritability.a, restricted.heritability.b = restricted.heritability.b, observed.a.b.correlation  = observed.a.b.correlation, percentage.range.for.neighbor.cor = c(.40,.70))

      b=ncol(dat)
      a=b-1

     datSNP =dat[,-c(a,b)]
     datExpr=dat[, c(a,b)]

if (sge.task==1)     r$greedy.100.causal[[k]] = NEOrobustness( datSNP, datExpr, TopGreedy=1:topK, TopForward=rep(1,topK), no.log.quiet=TRUE)
if (sge.task==2)     r$forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopGreedy=rep(1,topK), TopForward=1:topK, no.log.quiet=TRUE)
if (sge.task==3)     r$greedy.and.forward.100.causal[[k]] = NEOrobustness(datSNP, datExpr, TopForward=1:topK, TopGreedy=1:topK, no.log.quiet=TRUE)


     print(paste(sep="","Done with k=",k," of out ",num.reps," at ",date()))

     sge.task = get.sge.task()

     the.file=gsub(":",".",gsub(" ","_",paste(sep="","r.test.robust.task",sge.task,".results.",k,".",date(),".rdat")))
     save(r,file=the.file)

 } # end for num.reps


  r
}




# helper function for robustness.plot
   global.get.ylim = function(r,col.name) {
      ymin=Inf; ymax=-Inf;
      for (i in 1:length(r) ) { 
        ymin=min(ymin, as.numeric(r[[i]][,col.name]),na.rm=T)
        ymax=max(ymax, as.numeric(r[[i]][,col.name]),na.rm=T)
      }
      c(ymin,ymax)
    }



# r=cholesterol.pathway.robust()
cholesterol.pathway.robust=function() {

   # BxH Apoe null mouse cross.
   bxh.female.liver.cholesterol.genes =read.csv("bxh.female.liver.cholesterol.genes.csv")

   x=bxh.female.liver.cholesterol.genes[,-1] # get rid of mice names in first column
   
   datSNP.female = x[,1:1278]
   datExpr.female = x[,1279:1283]


   r=list() # hold the results in a descriptively named list
   
   top.Num.SNPs = 20
   
                                        # Greedy SNP selection
   r$greedy.Insig1.Mvd = NEOrobustnessGreedy( datSNP.female, datExpr.female[,c("Mvd","Insig1")], TopGreedySNPs=1:top.Num.SNPs, no.log.quiet=T)

   r$greedy.Insig1.Sqle = NEOrobustnessGreedy( datSNP.female, datExpr.female[,c("Sqle","Insig1")], TopGreedySNPs=1:top.Num.SNPs, no.log.quiet=T)

   r$greedy.Insig1.Dhcr7 = NEOrobustnessGreedy( datSNP.female, datExpr.female[,c("Dhcr7","Insig1")], TopGreedySNPs=1:top.Num.SNPs, no.log.quiet=T)

   r$greedy.Insig1.Fdft1 = NEOrobustnessGreedy( datSNP.female, datExpr.female[,c("Fdft1","Insig1")], TopGreedySNPs=1:top.Num.SNPs, no.log.quiet=T)

                                        # Forward SNP selection
   r$fwd.Insig1.Mvd = NEOrobustnessForward( datSNP.female, datExpr.female[,c("Mvd","Insig1")], TopForwardSNPs=1:top.Num.SNPs, no.log.quiet=T)

   r$fwd.Insig1.Sqle = NEOrobustnessForward( datSNP.female, datExpr.female[,c("Sqle","Insig1")], TopForwardSNPs=1:top.Num.SNPs, no.log.quiet=T)

   r$fwd.Insig1.Dhcr7 = NEOrobustnessForward( datSNP.female, datExpr.female[,c("Dhcr7","Insig1")], TopForwardSNPs=1:top.Num.SNPs, no.log.quiet=T)

   r$fwd.Insig1.Fdft1 = NEOrobustnessForward( datSNP.female, datExpr.female[,c("Fdft1","Insig1")], TopForwardSNPs=1:top.Num.SNPs, no.log.quiet=T)

   save(r,file="r.robust.cholesterol.with.CORRECT.snp.selection.20.rdat")


# given a r output list, and a col.name (such as "LEO.NB.FOR", or "LEO.NB.ALL", or zeo.for, zeo.all),
# find the ymin

   
   par(mfrow=c(4,2))

   robustness.plot(r$greedy.Insig1.Mvd,main="Insig1 to Mvd: Greedy SNPs")
   robustness.plot(r$fwd.Insig1.Mvd,main="Insig1 to Mvd: Forward SNPs")
   
   robustness.plot(r$greedy.Insig1.Sqle,main="Insig1 to Sqle : Greedy SNPs ")
   robustness.plot(r$fwd.Insig1.Sqle,main="Insig1 to Sqle : Forward SNPs ")
   
   robustness.plot(r$greedy.Insig1.Dhcr7,main="Insig1 to Dhcr7: Greedy SNPs")
   robustness.plot(r$fwd.Insig1.Dhcr7,main="Insig1 to Dhcr7: Forward SNPs")
   
   robustness.plot(r$greedy.Insig1.Fdft1,main="Insig1 to Fdft1: Greedy SNPs")
   robustness.plot(r$fwd.Insig1.Fdft1,main="Insig1 to Fdft1: Forward SNPs")
   





   r
}


post.test.robust=function() {
   x=r$greedy
   mn=c(-1.6376,.8388,.4952) #apply(x,2,mean)
   x[,1]=x[,1]-mn[1]
   x[,2]=x[,2]-mn[2]
   x[,3]=x[,3]-mn[3]

   matplot(x[,1:3],type="l",lty=1:4,lwd=3,col=1:4,ylab="LEO scores",xlab="Number of SNPs")
   legend(3.3,.8, c("LEO.NB.FORWARD"," SIMPLE.MAX.MAX"," LEO.NB.ALL"), col=1:3, lty=1:3,lwd=3, pch = "*",ncol = 1, cex=1.3)
   title("Simulation Study")

  r
}


# wrapper functions, to titrate the number of SNPs used
#
# datExpr is exactly two columns
# datSNP is greedy snps in order in which they are to be added, first 
# row first, then 1st and 2nd rows, then 1-3 rows, ..., until we get to TopNumberSNPs
#
# The model evaluated is column.2 -> column.1 of datExpr.
#
if(exists("NEOrobustnessGreedy") ) rm(NEOrobustnessGreedy);
NEOrobustnessGreedy=function(datSNP, datExpr, TopGreedySNPs=1:5, no.log.quiet = TRUE,   pm=neo.get.param()) {

   # replace by correct function instead of messing around
   return(NEOrobustness(datSNP, datExpr, TopGreedy=TopGreedySNPs, TopForward=rep(1,length(TopGreedySNPs)), no.log.quiet = no.log.quiet,   pm=pm))
   ############################################### DONE ##############################

   if (nrow(datSNP) != nrow(datExpr)) stop("dim on datSNP and datExpr must agree")

   rnS = rownames(datSNP)
   rnE = rownames(datExpr)

   if ((length(rnS)>0 & length(rnE)>0) & any(rnS != rnE)) stop("rownames of datSNP and datExpr must agree. Else cbind will crash.")

   # sanity check the supplied data frame -- it should not have any NA or non-varying columns
   if(any(is.na(datSNP))) { warning(paste("impute NA in datSNP before calling. We are repairing ",sum(is.na(datSNP))," NAs using: library(e1071); datSNP=impute(datSNP,what=\"median\");")); datSNP=impute(datSNP,what="median")}
   if(any(is.na(datExpr))) { warning(paste("impute NA in datExpr before calling. We are repairing ",sum(is.na(datExpr))," NAs using: library(e1071); datExpr=impute(datExpr,what=\"median\");")); datExpr=impute(datExpr,what="median") }
   if(ncol(datExpr) != 2) stop("datExpr must have exactly two columns.")
   if(nrow(datSNP) != nrow(datExpr)) stop("datExpr and datSNP must have matching rows.")

   pm$no.log=no.log.quiet
   pm$quiet =no.log.quiet
   pm$run.title = run.title = "robustness.greedy"
   pm$skip.snp.retraction=FALSE

   nruns = length(TopGreedySNPs)
   r=matrix(NA,ncol=16,nrow=nruns)
   colnames(r)=c("Num.SNP","c2.to.c1","LEO.NB.FOR","SIMPLE.MAX.MAX","LEO.NB.ALL","zeo.for","zeo.max","zeo.all","c1.to.c2","LEO.NB.FOR","SIMPLE.MAX.MAX","LEO.NB.ALL","zeo.for","zeo.max","zeo.all","SNPs")

   rownames(r)= paste(TopGreedySNPs)
   r[,1]=TopGreedySNPs

   snpcols = 1:ncol(datSNP)
   ###pm$turn.NA.scores.to.zero.and.zero.to.666 = TRUE
   pm$do.snp.pick.hcluster=FALSE

   # ======================================
   # important adjustments for greedy/fstep
   # ======================================
   pm$do.snp.pick.greedy=TRUE
   pm$do.snp.pick.fstep=FALSE



   # may need a few more in case of redundant SNPs
   pm$top.N.snps.per.trait = min(floor(2*max(TopGreedySNPs)),ncol(datSNP)) # take 200% more than needed, if we can.
   pm$orthog.search.depth = 1

   pm$skip.snp.select=FALSE
   cn=colnames(datExpr)

   #necessary
   pm$no.obs.Z = nrow(datExpr); # If you want, you can over-ride no.obs for Z-score computations here/and below.

   # INVAR: we have pm, set true no.obs rows
   pm$no.obs = nrow(datExpr);
   traitcols=c()

   if (pm$no.log) {
      neo.log.file="no.neo.log"
   } else {
     if (pm$keep.existing.neo.log.file) {
        if (length(pm$neo.log.file) < 1) stop("pm$keep.existing.neo.log.file was TRUE but pm$neo.log.file not specified. Halting.")
        # leave pm$neo.log.file intact
        neo.log.file = pm$neo.log.file
     } else {
       # setup logging - put all log files in a new, timestamped directory
       neo.log.file = neo.mkdir.logfile(pm$run.title);
       pm$neo.log.file=neo.log.file
     }
     sink(file = paste(sep="",neo.log.file,".txt"), append=TRUE, type="output",split = !pm$quiet);
     on.exit(while (sink.number()) sink()); # turn off log in case of crash, and upon normal exit.
     pm.print(pm,paste("Logfile:",neo.log.file));
   }

   datEarly =cbind(datSNP,datExpr)

   # find all the SNPs we'll ever want at once.
   r1 = neo.snp.primary(datEarly, snpcols=snpcols, traitcols=traitcols, 
                        pm=pm, 
                        skip.snp.select=FALSE,
                        neo.log.file=neo.log.file);

   # seperate out the SNPs into their respective categories
   my.snp = r1$my.snp

   if(any(!my.snp$non.redundant)) {
      warning("Note: Found some redundant SNPs. As a result, the number of SNPs utlized/the count may be off.")
      for (i in which(!my.snp$non.redundant)) {
          print(my.snp$redundant.note[i])
      }
      print("Here is my.snp: now we NEVER eliminate redundant SNPs because we were eliminating the best SNPs for B when they were previously chosen for A as well.:")
      print(my.snp)

      # we were eliminating the best SNPs for B when they were previously chosen for A as well;
      # so now we skip making them non-redundant
      # i.e. we don't do this: my.snp = my.snp[my.snp$non.redundant,]

      #print("AND here is my.snp AFTER we eliminate redundant SNPs:")
      #print(my.snp)
   }

   the.greedy = which(my.snp$my.pass == "greedy")
   the.fstep  = which(my.snp$my.pass == "fstep")

   col1rows   = which(my.snp$my.snp.pc == cn[1])
   col2rows   = which(my.snp$my.snp.pc == cn[2])

   # sort is probably redudant, but we want to be
   # sure we are taking SNPs in order.
   col1greedy = sort(intersect(the.greedy,col1rows))
   col2greedy = sort(intersect(the.greedy,col2rows))

   col1fstep = sort(intersect(the.fstep,col1rows))
   col2fstep = sort(intersect(the.fstep,col2rows))

   # get the ordered names of SNPs to index datSNP by
   col1greedy.cn = as.character(my.snp$my.snp.choice[col1greedy])
   col2greedy.cn = as.character(my.snp$my.snp.choice[col2greedy])

   col1fstep.cn = as.character(my.snp$my.snp.choice[col1fstep])
   col2fstep.cn = as.character(my.snp$my.snp.choice[col2fstep])

   available.col1 = length(col1greedy)
   available.col2 = length(col2greedy)

   max.found = max(available.col1, available.col2)

   for (i in 1:length(TopGreedySNPs)) {

     sim=TopGreedySNPs[i]

     if (sim > max.found) {
           print("breaking early b/c no more SNPs available")
           return(r) # no more SNPs available
     }

     pm$top.N.snps.per.trait = sim

     # INVAR: sim <= max(available.col1,available.col2)
     use.col1 = min(sim,available.col1)
     use.col2 = min(sim,available.col2)

     pm$skip.snp.select=TRUE

     final.chosen.snp.cols = c()
     if (use.col1 > 0) { final.chosen.snp.cols = col1greedy.cn[1:use.col1] }
     if (use.col2 > 0) { final.chosen.snp.cols = c(final.chosen.snp.cols,col2greedy.cn[1:use.col2]) }

     # make sure we aren't requesting the same set twice
     final.chosen.snp.cols = unique(final.chosen.snp.cols)

     snpcols=1:length(final.chosen.snp.cols)

     print(paste("         Greedy tried to pick: ",sim," SNPs and at this point am using ",length(final.chosen.snp.cols)," SNPs"))     

     if (length(final.chosen.snp.cols) < 1) {
        warning(paste("NEOrobustness function: Zero SNPs found at sim=",sim,"....skipping to next SNP count level."))
        next
     }
     
     # next 2 lines needed to handle one SNP situation
     datSNP.small = data.frame(datSNP[,final.chosen.snp.cols])
     colnames(datSNP.small) = final.chosen.snp.cols

     datCombined =cbind(datSNP.small,datExpr)
     snp.names = colnames(datCombined)[1:(ncol(datCombined)-ncol(datExpr))]

     pm$run.title = paste(sep="",run.title,sim)
     pm$quiet = F
     pm$no.log =F
     
     z = neo(datCombined,pm=pm,snpcols=snpcols)

    res.colname = paste(sep="",cn[2]," -> ",cn[1])
    res.colname.opp = paste(sep="",cn[1]," -> ",cn[2])
    r[i,2]=res.colname
    r[i,3]=anac(z$stats[res.colname,"LEO.NB.FOR"])
    r[i,4]=anac(z$stats[res.colname,"SIMPLE.MAX.MAX"])
    r[i,5]=anac(z$stats[res.colname,"LEO.NB.ALL"])
    r[i,6]=anac(z$stats[res.colname,"zeo.for"])
    r[i,7]=anac(z$stats[res.colname,"zeo.max"])
    r[i,8]=anac(z$stats[res.colname,"zeo.all"])

    r[i,9]=res.colname.opp
    r[i,10]=anac(z$stats[res.colname.opp,"LEO.NB.FOR"])
    r[i,11]=anac(z$stats[res.colname.opp,"SIMPLE.MAX.MAX"])
    r[i,12]=anac(z$stats[res.colname.opp,"LEO.NB.ALL"])
    r[i,13]=anac(z$stats[res.colname.opp,"zeo.for"])
    r[i,14]=anac(z$stats[res.colname.opp,"zeo.max"])
    r[i,15]=anac(z$stats[res.colname.opp,"zeo.all"])
    r[i,16]=paste(collapse=",",snp.names)

    print(paste(sep="","Progress through the NEOrobustnessGreedy loop: done with ",i," out of ",length(TopGreedySNPs)," at ",date()))     
  }

  #  Output=data.frame(LEO.forward,  LEO.max, LEO.all,  ZEO.forward,  ZEO.max, ZEO.all )

  # Where the rows of the data frame output correspond to TopGreedySNPs that are assigned to each network node (greedy search for SNPs).

  r
}

# version that doesn't try to factor out snp selection to save time, 
# b/c we had worries about correctness when redundant SNPs enter the picture
#
# specify the number of greedy and forward snps manually to use for each run
#
# Set back.compatible.column.names to TRUE to re-use old scripts that expect
#     LEO.NB.OCA to be called LEO.NB.FOR
# and LEO.NB.CPA to be called LEO.MULTIONE.NB
# 
if(exists("NEOrobustness") ) rm(NEOrobustness);
NEOrobustness=function(datSNP, datExpr, TopGreedy=1:5, TopForward=1:5, no.log.quiet = TRUE,   pm=neo.get.param(), back.compatible.column.names=F) {

   if (!is.null(pm$A) | !is.null(pm$B)) {
      pm.print(pm,"Ignoring pm$A and pm$B specifications during NEOrobustness analysis--as is only reasonable.")
      pm$A = NULL
      pm$B = NULL
   }

   if (length(TopGreedy) != length(TopForward)) stop("length of TopGreedy and TopForward must agree (but can have all 1's for instance, if need be.")

   if (nrow(datSNP) != nrow(datExpr)) stop("dim on datSNP and datExpr must agree")

   rnS = rownames(datSNP)
   rnE = rownames(datExpr)

   if ((length(rnS)>0 & length(rnE)>0) & any(rnS != rnE)) stop("rownames of datSNP and datExpr must agree. Else cbind will crash.")

   # sanity check the supplied data frame -- it should not have any NA or non-varying columns
   if(any(is.na(datSNP))) { warning(paste("impute NA in datSNP before calling. We are repairing ",sum(is.na(datSNP))," NAs using: library(e1071); datSNP=impute(datSNP,what=\"median\");")); datSNP=impute(datSNP,what="median")}
   if(any(is.na(datExpr))) { warning(paste("impute NA in datExpr before calling. We are repairing ",sum(is.na(datExpr))," NAs using: library(e1071); datExpr=impute(datExpr,what=\"median\");")); datExpr=impute(datExpr,what="median") }
   if(ncol(datExpr) != 2) stop("datExpr must have exactly two columns.")
   if(nrow(datSNP) != nrow(datExpr)) stop("datExpr and datSNP must have matching rows.")

   pm$no.log=no.log.quiet
   pm$quiet =no.log.quiet
   pm$run.title = run.title = "robustness.greedy"
   pm$skip.snp.retraction=FALSE

   nruns = length(TopGreedy)
   r=matrix(NA,ncol=21,nrow=nruns)

   if (back.compatible.column.names) {
      colnames(r)=c("Num.Greedy.SNP","c2.to.c1","LEO.NB.FOR","SIMPLE.MAX.MAX","LEO.NB.MAX","LEO.NB.ALL","LEO.MULTIONE.NB","zeo.for","zeo.max","zeo.all","c1.to.c2","LEO.NB.FOR","SIMPLE.MAX.MAX","LEO.NB.MAX","LEO.NB.ALL","LEO.MULTIONE.NB","zeo.for","zeo.max","zeo.all","SNPs","Num.Forward.SNP")
   } else {
      colnames(r)=c("Num.Greedy.SNP","c2.to.c1","LEO.NB.OCA","SIMPLE.MAX.MAX","LEO.NB.MAX","LEO.NB.ALL","LEO.NB.CPA","zeo.for","zeo.max","zeo.all","c1.to.c2","LEO.NB.OCA","SIMPLE.MAX.MAX","LEO.NB.MAX","LEO.NB.ALL","LEO.NB.CPA","zeo.for","zeo.max","zeo.all","SNPs","Num.Forward.SNP")
   }

   rownames(r)= paste(sep=",",paste(TopGreedy),paste(TopForward))
   r[,1]=TopGreedy
   r[,21] = TopForward

   snpcols = 1:ncol(datSNP)
   ###pm$turn.NA.scores.to.zero.and.zero.to.666 = TRUE
   pm$do.snp.pick.hcluster=FALSE

   # ======================================
   # important adjustments for greedy/fstep
   # ======================================
   pm$do.snp.pick.greedy=TRUE
   pm$do.snp.pick.fstep=TRUE
 
   pm$skip.snp.select=FALSE
   cn=colnames(datExpr)

   #necessary
   pm$no.obs.Z = nrow(datExpr); # If you want, you can over-ride no.obs for Z-score computations here/and below.

   # INVAR: we have pm, set true no.obs rows
   pm$no.obs = nrow(datExpr);
   traitcols=c()
   datCombined =cbind(datSNP,datExpr)
   snpcols=1:ncol(datSNP)

   for (i in 1:length(TopGreedy)) {

     sim=TopGreedy[i]
     sim.f=TopForward[i]

     pm$top.N.snps.per.trait = sim
     pm$orthog.search.depth = sim.f
     pm$skip.snp.select=FALSE

     pm$run.title = paste(sep="",run.title,sim)

     trait.col1 = ncol(datSNP) + 1
     trait.col2 = ncol(datSNP) + 2
     traitcols = c(trait.col1,trait.col2)

     z = neo(datCombined,pm=pm,snpcols=snpcols,traitcols=traitcols)

     snp.names = "" # punt on reporting this for now

     res.colname = paste(sep="",cn[2]," -> ",cn[1])
     res.colname.opp = paste(sep="",cn[1]," -> ",cn[2])
    r[i,2]=res.colname
    r[i,3]=anac(z$stats[res.colname,"LEO.NB.OCA"])
    r[i,4]=anac(z$stats[res.colname,"SIMPLE.MAX.MAX"])
    r[i,5]=anac(z$stats[res.colname,"LEO.NB.MAX"])
    r[i,6]=anac(z$stats[res.colname,"LEO.NB.ALL"])
    r[i,7]=anac(z$stats[res.colname,"LEO.NB.CPA"])

    r[i,8]=anac(z$stats[res.colname,"zeo.for"])
    r[i,9]=anac(z$stats[res.colname,"zeo.max"])
    r[i,10]=anac(z$stats[res.colname,"zeo.all"])

    r[i,11]=res.colname.opp
    r[i,12]=anac(z$stats[res.colname.opp,"LEO.NB.OCA"])
    r[i,13]=anac(z$stats[res.colname.opp,"SIMPLE.MAX.MAX"])
    r[i,14]=anac(z$stats[res.colname.opp,"LEO.NB.MAX"])
    r[i,15]=anac(z$stats[res.colname.opp,"LEO.NB.ALL"])
    r[i,16]=anac(z$stats[res.colname.opp,"LEO.NB.CPA"])     
    r[i,17]=anac(z$stats[res.colname.opp,"zeo.for"])
    r[i,18]=anac(z$stats[res.colname.opp,"zeo.max"])
    r[i,19]=anac(z$stats[res.colname.opp,"zeo.all"])
    r[i,20]=paste(collapse=",",snp.names)

    print(paste(sep="","Progress through the NEOrobustness loop: done with ",i," out of ",length(TopGreedy)," at ",date()))     
  }

  #  Output=data.frame(LEO.forward,  LEO.max, LEO.all,  ZEO.forward,  ZEO.max, ZEO.all )

  # Where the rows of the data frame output correspond to TopGreedy that are assigned to each network node (greedy search for SNPs).

  r
}
 


 

NEOrobustnessForward=function(datSNP, datExpr, TopForwardSNPs=1:6, no.log.quiet=TRUE,pm=neo.get.param()){

   # replace by correct function instead of messing around
   return(NEOrobustness(datSNP, datExpr, TopForward=TopForwardSNPs, TopGreedy=rep(1,length(TopForwardSNPs)), no.log.quiet = no.log.quiet,   pm=pm))
   ############################################### DONE ##############################


   # sanity check the supplied data frame -- it should not have any NA or non-varying columns

   if(any(is.na(datSNP))) { warning(paste("impute NA in datSNP before calling. We are repairing ",sum(is.na(datSNP))," NAs using: library(e1071); datSNP=impute(datSNP,what=\"median\");")); datSNP=impute(datSNP,what="median")}
   if(any(is.na(datExpr))) { warning(paste("impute NA in datExpr before calling. We are repairing ",sum(is.na(datExpr))," NAs using: library(e1071); datExpr=impute(datExpr,what=\"median\");")); datExpr=impute(datExpr,what="median") }
   if(ncol(datExpr) != 2) stop("datExpr must have exactly two columns.")
   if(nrow(datSNP) != nrow(datExpr)) stop("datExpr and datSNP must have matching rows.")

   pm$no.log=no.log.quiet
   pm$quiet =no.log.quiet
   pm$run.title = run.title = "robustness.forward"
   pm$skip.snp.retraction=TRUE

   nruns = length(TopForwardSNPs)
   r=matrix(NA,ncol=16,nrow=nruns)
   colnames(r)=c("Num.SNP","c2.to.c1","LEO.NB.FOR","SIMPLE.MAX.MAX","LEO.NB.ALL","zeo.for","zeo.max","zeo.all","c1.to.c2","LEO.NB.FOR","SIMPLE.MAX.MAX","LEO.NB.ALL","zeo.for","zeo.max","zeo.all","SNPs")
   rownames(r)= paste(TopForwardSNPs)
   r[,1]=TopForwardSNPs

   snpcols = 1:ncol(datSNP)
###   pm$turn.NA.scores.to.zero.and.zero.to.666 = TRUE
   pm$do.snp.pick.hcluster=FALSE

   # ======================================
   # important adjustments for greedy/fstep
   # ======================================
   pm$do.snp.pick.greedy=FALSE
   pm$do.snp.pick.fstep=TRUE


   # may need a few more in case of redundant SNPs
   pm$top.N.snps.per.trait = 1
   pm$orthog.search.depth = min(floor(2*max(TopForwardSNPs)),ncol(datSNP)) # take 200% more than needed, if we can.

   pm$skip.snp.select=FALSE
   cn=colnames(datExpr)

   #necessary
   pm$no.obs.Z = nrow(datExpr); # If you want, you can over-ride no.obs for Z-score computations here/and below.

   # INVAR: we have pm, set true no.obs rows
   pm$no.obs = nrow(datExpr);
   traitcols=c()

   if (pm$no.log) {
      neo.log.file="no.neo.log"
   } else {
     # setup logging - put all log files in a new, timestamped directory
     neo.log.file = neo.mkdir.logfile(pm$run.title);
     pm$neo.log.file=neo.log.file
     sink(file = paste(sep="",neo.log.file,".txt"), append=TRUE, type="output",split = !pm$quiet);
     on.exit(while (sink.number()) sink()); # turn off log in case of crash, and upon normal exit.
     pm.print(pm,paste("Logfile:",neo.log.file));
   }

   datEarly =cbind(datSNP,datExpr)

   # find all the SNPs we'll ever want at once.
   r1 = neo.snp.primary(datEarly, snpcols=snpcols, traitcols=traitcols, 
                        pm=pm, 
                        skip.snp.select=FALSE,
                        neo.log.file=neo.log.file);

   # seperate out the SNPs into their respective categories
   my.snp = r1$my.snp

   if(any(!my.snp$non.redundant)) {
      warning("Note: Found some redundant SNPs. As a result, the number of SNPs utlized/the count may be off.")
      for (i in which(!my.snp$non.redundant)) {
          print(my.snp$redundant.note[i])
      }
      print("Here is my.snp: now we NEVER eliminate redundant SNPs because we were eliminating the best SNPs for B when they were previously chosen for A as well.:")
      print(my.snp)

      # we were eliminating the best SNPs for B when they were previously chosen for A as well;
      # so now we skip making them non-redundant
      # i.e. we don't do this: my.snp = my.snp[my.snp$non.redundant,]

      #print("AND here is my.snp AFTER we eliminate redundant SNPs:")
      #print(my.snp)
   }

   the.greedy = which(my.snp$my.pass == "greedy")
   the.fstep  = which(my.snp$my.pass == "fstep")

   col1rows   = which(my.snp$my.snp.pc == cn[1])
   col2rows   = which(my.snp$my.snp.pc == cn[2])

   # sort is probably redudant, but we want to be
   # sure we are taking SNPs in order.
   col1greedy = sort(intersect(the.greedy,col1rows))
   col2greedy = sort(intersect(the.greedy,col2rows))

   col1fstep = sort(intersect(the.fstep,col1rows))
   col2fstep = sort(intersect(the.fstep,col2rows))

   # get the ordered names of SNPs to index datSNP by
   col1greedy.cn = as.character(my.snp$my.snp.choice[col1greedy])
   col2greedy.cn = as.character(my.snp$my.snp.choice[col2greedy])

   col1fstep.cn = as.character(my.snp$my.snp.choice[col1fstep])
   col2fstep.cn = as.character(my.snp$my.snp.choice[col2fstep])

   available.col1 = length(col1fstep)
   available.col2 = length(col2fstep)

   max.found = max(available.col1, available.col2)

   print(paste(sep="","starting NEOrobustnessForward at ",date()))

   for (i in 1:length(TopForwardSNPs)) {

     sim=TopForwardSNPs[i]

     if (sim > max.found) {
           print("breaking early b/c no more SNPs available")
           return(r) # no more SNPs available
     }

     pm$orthog.search.depth = sim
#     pm$top.N.snps.per.trait = sim

     # INVAR: sim <= max(available.col1,available.col2)
     use.col1 = min(sim,available.col1)
     use.col2 = min(sim,available.col2)

     pm$skip.snp.select=TRUE

     # fstep version
     final.chosen.snp.cols = c()
     if (use.col1 > 0) { final.chosen.snp.cols = col1fstep.cn[1:use.col1] }
     if (use.col2 > 0) { final.chosen.snp.cols = c(final.chosen.snp.cols,col2fstep.cn[1:use.col2]) }

     # make sure we aren't requesting the same set twice
     final.chosen.snp.cols = unique(final.chosen.snp.cols)

     snpcols=1:length(final.chosen.snp.cols)

     print(paste("         Forward tried to pick: ",sim," SNPs and at this point am using ",length(final.chosen.snp.cols)," SNPs"))     

     if (length(final.chosen.snp.cols) < 1) {
        warning(paste("NEOrobustness function: Zero SNPs found at sim=",sim,"....skipping to next SNP count level."))
        next
     }
     
     # next 2 lines needed to handle one SNP situation
     datSNP.small = data.frame(datSNP[,final.chosen.snp.cols])
     colnames(datSNP.small) = final.chosen.snp.cols

     datCombined =cbind(datSNP.small,datExpr)
     snp.names = colnames(datCombined)[1:(ncol(datCombined)-ncol(datExpr))]

     pm$run.title = paste(sep="",run.title,sim)
     z = neo(datCombined,pm=pm,snpcols=snpcols)

    res.colname = paste(sep="",cn[2]," -> ",cn[1])
    res.colname.opp = paste(sep="",cn[1]," -> ",cn[2])
    r[i,2]=res.colname
    r[i,3]=anac(z$stats[res.colname,"LEO.NB.FOR"])
    r[i,4]=anac(z$stats[res.colname,"SIMPLE.MAX.MAX"])
    r[i,5]=anac(z$stats[res.colname,"LEO.NB.ALL"])
    r[i,6]=anac(z$stats[res.colname,"zeo.for"])
    r[i,7]=anac(z$stats[res.colname,"zeo.max"])
    r[i,8]=anac(z$stats[res.colname,"zeo.all"])

    r[i,9]=res.colname.opp
    r[i,10]=anac(z$stats[res.colname.opp,"LEO.NB.FOR"])
    r[i,11]=anac(z$stats[res.colname.opp,"SIMPLE.MAX.MAX"])
    r[i,12]=anac(z$stats[res.colname.opp,"LEO.NB.ALL"])
    r[i,13]=anac(z$stats[res.colname.opp,"zeo.for"])
    r[i,14]=anac(z$stats[res.colname.opp,"zeo.max"])
    r[i,15]=anac(z$stats[res.colname.opp,"zeo.all"])
    r[i,16]=paste(collapse=",",snp.names)

    print(paste(sep="","Progress through the NEOrobustnessForward loop: done with ",i," out of ",length(TopForwardSNPs)," at ",date()))
  }

  r
}



NEOrobustnessForwardAndGreedy=function(datSNP, datExpr, TopForwardAndGreedySNPs=1:6, no.log.quiet=TRUE,pm=neo.get.param()){

   # replace by correct function instead of messing around
   return(NEOrobustness(datSNP, datExpr, TopGreedy=TopForwardAndGreedySNPs, TopForward=TopForwardAndGreedySNPs, no.log.quiet = no.log.quiet,   pm=pm))
   ############################################### DONE ##############################


   # sanity check the supplied data frame -- it should not have any NA or non-varying columns

   if(any(is.na(datSNP))) { warning(paste("impute NA in datSNP before calling. We are repairing ",sum(is.na(datSNP))," NAs using: library(e1071); datSNP=impute(datSNP,what=\"median\");")); datSNP=impute(datSNP,what="median")}
   if(any(is.na(datExpr))) { warning(paste("impute NA in datExpr before calling. We are repairing ",sum(is.na(datExpr))," NAs using: library(e1071); datExpr=impute(datExpr,what=\"median\");")); datExpr=impute(datExpr,what="median") }

   if(ncol(datExpr) != 2) stop("datExpr must have exactly two columns.")
   if(nrow(datSNP) != nrow(datExpr)) stop("datExpr and datSNP must have matching rows.")

   pm$no.log=no.log.quiet
   pm$quiet =no.log.quiet
   pm$run.title = run.title = "robustness.greedy.and.forward"
   pm$skip.snp.retraction=TRUE

   nruns = length(TopForwardAndGreedySNPs)
   r=matrix(NA,ncol=16,nrow=nruns)
   colnames(r)=c("Num.SNP","c2.to.c1","LEO.NB.FOR","SIMPLE.MAX.MAX","LEO.NB.ALL","zeo.for","zeo.max","zeo.all","c1.to.c2","LEO.NB.FOR","SIMPLE.MAX.MAX","LEO.NB.ALL","zeo.for","zeo.max","zeo.all","SNPs")

   rownames(r)= paste(TopForwardAndGreedySNPs)
   r[,1]=TopForwardAndGreedySNPs

   snpcols = 1:ncol(datSNP)
###   pm$turn.NA.scores.to.zero.and.zero.to.666 = TRUE
   pm$do.snp.pick.hcluster=FALSE

   # ======================================
   # important adjustments for greedy/fstep
   # ======================================
   pm$do.snp.pick.greedy=TRUE
   pm$do.snp.pick.fstep=TRUE



   # may need a few more in case of redundant SNPs
   pm$top.N.snps.per.trait = min(floor(2*max(TopForwardAndGreedySNPs)),ncol(datSNP)) # take 200% more than needed, if we can.
   pm$orthog.search.depth  = min(floor(2*max(TopForwardAndGreedySNPs)),ncol(datSNP)) # take 200% more than needed, if we can.


   pm$skip.snp.select=FALSE
   cn=colnames(datExpr)

   #necessary
   pm$no.obs.Z = nrow(datExpr); # If you want, you can over-ride no.obs for Z-score computations here/and below.

   # INVAR: we have pm, set true no.obs rows
   pm$no.obs = nrow(datExpr);
   traitcols=c()

   if (pm$no.log) {
      neo.log.file="no.neo.log"
   } else {
     # setup logging - put all log files in a new, timestamped directory
     neo.log.file = neo.mkdir.logfile(pm$run.title);
     pm$neo.log.file=neo.log.file
     sink(file = paste(sep="",neo.log.file,".txt"), append=TRUE, type="output",split = !pm$quiet);
     on.exit(while (sink.number()) sink()); # turn off log in case of crash, and upon normal exit.
     pm.print(pm,paste("Logfile:",neo.log.file));
   }

   datEarly =cbind(datSNP,datExpr)

   # find all the SNPs we'll ever want at once.
   r1 = neo.snp.primary(datEarly, snpcols=snpcols, traitcols=traitcols, 
                        pm=pm, 
                        skip.snp.select=FALSE,
                        neo.log.file=neo.log.file);

   # seperate out the SNPs into their respective categories
   my.snp = r1$my.snp


   if(any(!my.snp$non.redundant)) {
      warning("Note: Found some redundant SNPs. As a result, the number of SNPs utlized/the count may be off.")
      for (i in which(!my.snp$non.redundant)) {
          print(my.snp$redundant.note[i])
      }
      print("Here is my.snp: now we NEVER eliminate redundant SNPs because we were eliminating the best SNPs for B when they were previously chosen for A as well.:")
      print(my.snp)

      # we were eliminating the best SNPs for B when they were previously chosen for A as well;
      # so now we skip making them non-redundant
      # i.e. we don't do this: my.snp = my.snp[my.snp$non.redundant,]

      #print("AND here is my.snp AFTER we eliminate redundant SNPs:")
      #print(my.snp)
   }

   the.greedy = which(my.snp$my.pass == "greedy")
   the.fstep  = which(my.snp$my.pass == "fstep")

   col1rows   = which(my.snp$my.snp.pc == cn[1])
   col2rows   = which(my.snp$my.snp.pc == cn[2])

   # eliminate NA before proceeding
   

   # sort is probably redudant, but we want to be
   # sure we are taking SNPs in order.
   col1greedy = sort(intersect(the.greedy,col1rows))
   col2greedy = sort(intersect(the.greedy,col2rows))

   col1fstep = sort(intersect(the.fstep,col1rows))
   col2fstep = sort(intersect(the.fstep,col2rows))

   # get the ordered names of SNPs to index datSNP by
   col1greedy.cn = as.character(my.snp$my.snp.choice[col1greedy])
   col2greedy.cn = as.character(my.snp$my.snp.choice[col2greedy])

   col1fstep.cn = as.character(my.snp$my.snp.choice[col1fstep])
   col2fstep.cn = as.character(my.snp$my.snp.choice[col2fstep])

   available.col1.f = length(col1fstep)
   available.col1.g = length(col1greedy)

   available.col2.f = length(col2fstep)
   available.col2.g = length(col2greedy)

   available.col1 = max(available.col1.f,available.col1.g)
   available.col2 = max(available.col2.f,available.col2.g)

   max.found = max(available.col1, available.col2) 

   for (i in 1:length(TopForwardAndGreedySNPs)) {

     sim=TopForwardAndGreedySNPs[i]

     if (sim > max.found) {
           print("breaking early b/c no more SNPs available")
           return(r) # no more SNPs available
     }

     fwd.depth = ceiling((sim+1)/2)
     greedy.depth = ceiling(sim/2)

     pm$top.N.snps.per.trait = greedy.depth
     pm$orthog.search.depth =  fwd.depth

     # INVAR: sim <= max(available.col1,available.col2)
     use.col1.g = min(sim,available.col1.g)
     use.col1.f = min(sim,available.col1.f)

     use.col2.g = min(sim,available.col2.g)
     use.col2.f = min(sim,available.col2.f)

     # fstep and greedy version
     pm$skip.snp.select=TRUE

     final.chosen.snp.cols = c()
     if (use.col1.g > 0) { final.chosen.snp.cols = col1greedy.cn[1:use.col1.g] }
     if (use.col2.g > 0) { final.chosen.snp.cols = c(final.chosen.snp.cols,col2greedy.cn[1:use.col2.g]) }
     if (use.col1.f > 0) { final.chosen.snp.cols = c(final.chosen.snp.cols,col1fstep.cn[1:use.col1.f]) }
     if (use.col2.f > 0) { final.chosen.snp.cols = c(final.chosen.snp.cols,col2fstep.cn[1:use.col2.f]) }

     # make sure we aren't requesting the same set twice
     final.chosen.snp.cols = unique(final.chosen.snp.cols)

     snpcols=1:length(final.chosen.snp.cols)

     print(paste("      sim=",sim,"   Greedy+Forward tried to pick: ",fwd.depth," fwd SNPs; and ",greedy.depth," greedy SNPs and at this point am using ",length(final.chosen.snp.cols)," SNPs"))     

     if (length(final.chosen.snp.cols) < 1) {
        warning(paste("NEOrobustness function: Zero SNPs found at sim=",sim,"....skipping to next SNP count level."))
        next
     }
     
     # next 2 lines needed to handle one SNP situation
     datSNP.small = data.frame(datSNP[,final.chosen.snp.cols])
     colnames(datSNP.small) = final.chosen.snp.cols

     datCombined =cbind(datSNP.small,datExpr)
     snp.names = colnames(datCombined)[1:(ncol(datCombined)-ncol(datExpr))]

     pm$run.title = paste(sep="",run.title,sim)
     z = neo(datCombined,pm=pm,snpcols=snpcols)


    res.colname = paste(sep="",cn[2]," -> ",cn[1])
    res.colname.opp = paste(sep="",cn[1]," -> ",cn[2])
    r[i,2]=res.colname
    r[i,3]=anac(z$stats[res.colname,"LEO.NB.FOR"])
    r[i,4]=anac(z$stats[res.colname,"SIMPLE.MAX.MAX"])
    r[i,5]=anac(z$stats[res.colname,"LEO.NB.ALL"])
    r[i,6]=anac(z$stats[res.colname,"zeo.for"])
    r[i,7]=anac(z$stats[res.colname,"zeo.max"])
    r[i,8]=anac(z$stats[res.colname,"zeo.all"])

    r[i,9]=res.colname.opp
    r[i,10]=anac(z$stats[res.colname.opp,"LEO.NB.FOR"])
    r[i,11]=anac(z$stats[res.colname.opp,"SIMPLE.MAX.MAX"])
    r[i,12]=anac(z$stats[res.colname.opp,"LEO.NB.ALL"])
    r[i,13]=anac(z$stats[res.colname.opp,"zeo.for"])
    r[i,14]=anac(z$stats[res.colname.opp,"zeo.max"])
    r[i,15]=anac(z$stats[res.colname.opp,"zeo.all"])
    r[i,16]=paste(collapse=",",snp.names)

    print(paste(sep="","Progress through the NEOrobustnessForwardAndGreedy loop: done with ",i," out of ",length(TopForwardAndGreedySNPs)," at ",date()))     
  }

  r
}


#
#
# neo.monte.carlo.data() : simulate test data under controlled, parameterized conditions
#
#  The model is M.A -> A <- B <- M.B  where M.A and M.B are sets of markers.
#
#    parameters: 
#
# Number of primary causal SNPs into A and B nodes
#    do.a.snp=5, so.b.snp=5,
#
# Number of back-ground non-causal (independent) SNPs
#    noise.a.snp=10, noise.b.snp=10,
#
# Number of SNPs nearby the QTL LOD peak causal SNPS, the number for each causal SNP.
# For example, if nearby.a.snps = 6 means that 6 neighbor SNPs are generated for each do.a.snp
#    nearby.a.snps=6, nearby.b.snps=6,
#
# the correlation of the neighbors (nearby SNPs) with the trait is roughly specified by this range
#    percentage.range.for.neighbor.cor = c(.75,.95),
#
# key parameters for the genetic versus environmental noise into A and B.
# Note that not all values of parameters are compatible.
#    restricted.heritability.a = .5,
#    restricted.heritability.b = .5,
#    observed.a.b.correlation  = .4,
#
# set percent.causal.signal to 0 to simulate the fully confounded model
#    percent.causal.signal     = 1,
#
#
# The (percentage) weights for the SNPs coming into A
# We must supply alpha's for noise SNPs too: 0 for true noise, as here.
#    alpha.A.snp = c(rep(sqrt(1/do.a.snp),do.a.snp),rep(0,noise.a.snp) ),
#
# The weights for the SNPs coming into B
#    beta.B.snp  = c(rep(sqrt(1/do.b.snp),do.b.snp), rep(0,noise.b.snp) )
#
#
if(exists("neo.monte.carlo.data") ) rm(neo.monte.carlo.data);
neo.monte.carlo.data=function(no.samples=200,
   do.a.snp=5,
   do.b.snp=5,
   noise.a.snp=0,
   noise.b.snp=0,
   nearby.a.snps=6,
   nearby.b.snps=6,
   percentage.range.for.neighbor.cor = c(.50,.80),
   restricted.heritability.a = .5,
   restricted.heritability.b = .5,
   observed.a.b.correlation  = .4,
   percent.causal.signal = 1,
   alpha.A.snp = c(rep(sqrt(1/do.a.snp),do.a.snp),rep(0,noise.a.snp) ), 
   beta.B.snp  = c(rep(sqrt(1/do.b.snp),do.b.snp), rep(0,noise.b.snp))
) {

     # check if we have feasible parameters
     ha=restricted.heritability.a
     hb=restricted.heritability.b
     cor.ab = observed.a.b.correlation
     w = cor.ab * percent.causal.signal
     g2= cor.ab - w

     if (g2 < 0) stop("gamma-squared is less than zero: non-feasable heritabilities, observed.a.b.correlation, and percent.causal.signal")

     g = sqrt(g2)
     errb = 1-hb-g2
     erra = 1-ha-w^2-g^2-2*(g*g*w)

     if (erra <=0) stop("erra should be non-negativity: non-feasable parameters specified.")
     if (errb <=0) stop("errb should be non-negativity: non-feasable parameters specified.")

     # now start simulation stuff
     #set.seed(1);
     N=no.samples

     nearby.yes = FALSE
     if (nearby.a.snps > 0 || nearby.b.snps > 0) { nearby.yes = TRUE  }

     numAsnp=do.a.snp+noise.a.snp 
     numBsnp=do.b.snp+noise.b.snp 

     datSNP=data.frame(matrix(NA, nrow=no.samples, numAsnp+numBsnp))

     if (noise.a.snp > 0) { 
          if (do.a.snp > 0) {
             names(datSNP)[1:numAsnp] = c(paste("SNP.A.",1:do.a.snp, sep=""),paste("SNP.A.NOISE.",(1:noise.a.snp)+do.a.snp, sep=""))
          } else {
             names(datSNP)[1:numAsnp] = c(paste("SNP.A.NOISE.",(1:noise.a.snp), sep=""))
          }
     } else {
          if (do.a.snp > 0) { names(datSNP)[1:do.a.snp] = paste("SNP.A.",1:do.a.snp, sep="") }
     }

     Bspan = (1:numBsnp) + numAsnp
     if (noise.b.snp > 0) { 
         if (do.b.snp > 0) {
            names(datSNP)[Bspan] = c(paste("SNP.B.",1:do.b.snp, sep=""),paste("SNP.B.NOISE.",(1:noise.b.snp)+do.b.snp, sep=""))
         } else {
            names(datSNP)[Bspan] = c(paste("SNP.B.NOISE.",(1:noise.b.snp), sep=""))
         }
     } else {
          if (do.b.snp > 0) { names(datSNP)[Bspan] = paste("SNP.B.",1:do.b.snp, sep="") }
     }

     SNP.A=SNP.B=0

    if (numAsnp > 0) {
     for (i in 1:numAsnp) {
       SNPa=sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T)
       SNP.A=SNP.A+ (alpha.A.snp[i] * c(scale(SNPa)))
       datSNP[,i]=SNPa
     } # end of for loop
    }

   k=1

   if (numBsnp > 0) {
    for (i in Bspan) {
      SNPa=sample(c(0,1,2), size=no.samples, prob=c(.25,.5,.25) ,replace=T)
      SNP.B=SNP.B+ (beta.B.snp[k] * c(scale(SNPa)))
      datSNP[,i]=SNPa
      k=k+1
    } # end of for loop
   }

   # now generate the nearby SNPs too

   nearB = do.b.snp*nearby.b.snps
   nearA = do.a.snp*nearby.a.snps
   nearTot = nearA+nearB
   nearbySNP=data.frame(matrix(NA, nrow=no.samples, nearTot))
   if (nearA > 0) { names(nearbySNP)[1:nearA] = paste("SNP.near.A.",1:nearA, sep="") }
   nearBspan = (1:nearB) + nearA
   if (nearB > 0) { names(nearbySNP)[nearBspan] = paste("SNP.near.B.",1:nearB, sep="") }

   if (nearA > 0) {
   for (i in 1:nearA) {
      which.a = 1+floor((i-1)/nearby.a.snps)
      p.off=runif(1,min=percentage.range.for.neighbor.cor[1], max=percentage.range.for.neighbor.cor[2])
      num.to.alter = floor(N*(1-p.off))
      SNPa=datSNP[,which.a]
      
      new.val=sample(c(0,1,2), size=num.to.alter, prob=c(.25,.5,.25) ,replace=T)
      SNPa[sample(1:N,size=num.to.alter,replace=F)] = new.val

      nearbySNP[,i]=(SNPa)
      colnames(nearbySNP)[i]=paste(sep="","SNP.neighbor.",i,".near.",colnames(datSNP)[which.a])
   }}

   if (nearB > 0) {
   for (i in nearBspan) {
      which.b = 1+floor((i-(nearA+1))/nearby.b.snps)
      p.off=runif(1,min=percentage.range.for.neighbor.cor[1], max=percentage.range.for.neighbor.cor[2])

      num.to.alter = floor(N*(1-p.off))
      SNPa=datSNP[,Bspan[which.b]]
      
      new.val=sample(c(0,1,2), size=num.to.alter, prob=c(.25,.5,.25) ,replace=T)
      SNPa[sample(1:N,size=num.to.alter,replace=F)] = new.val

      nearbySNP[,i]=(SNPa)
      colnames(nearbySNP)[i]=paste(sep="","SNP.neighbor.",i,".near.",colnames(datSNP)[Bspan[which.b]])
   }}

     C = rnorm(N);

     if (do.b.snp > 0) {
        B= c(scale(sqrt(hb) * scale(SNP.B) + rnorm(N,sd=sqrt(errb)) + g*C))
     } else {
        B= c(rnorm(N,sd=sqrt(errb)) + g*C)
     }

     if (do.a.snp == 0) {
         A= c(rnorm(N,sd=sqrt(erra)) + g*C + w*B)
     } else {
         A= c(scale(sqrt(ha) * scale(SNP.A) + rnorm(N,sd=sqrt(erra)) + g*C + w*B))
     }

     if (nearby.yes) {
        datCombined=data.frame(cbind(datSNP,nearbySNP,A,B));
     } else {
        datCombined=data.frame(cbind(datSNP,A,B));
     }    

     datCombined  
}

als.analysis=function() {


#  x=read.csv(file="c:/good/tova.als/pTables/turquoiseTable.maxrank.csv")
#    x=read.csv(file="c:/good/tova.als/pTables/turquoiseTable.avgrank.csv")
    x=read.csv(file="c:/good/tova.als/pTables/turquoiseTable.csv")

  snpnames=as.character(x[,1])

  pvalue.trait = x[,2]
  pvalue.module = x[,3]

  trait.order = sort(pvalue.trait,decreasing=F,index.return=T)$ix
  module.order = sort(pvalue.module,decreasing=F,index.return=T)$ix

  topNmarkers = 60

  snpnames.trait.rank = snpnames[trait.order[1:topNmarkers]]
  snpnames.module.rank= snpnames[module.order[1:topNmarkers]]

#  max.rank.index = x[,2]

  m=load("c:/good/tova.als/RobustnessAnalysis/turquoiseSNP1.Rdat")

  which.cols = pmatch(snpnames.trait.rank,colnames(datCombinedTurquoise1imp))
  which.cols.mod = pmatch(snpnames.module.rank,colnames(datCombinedTurquoise1imp))

  # some where NA, presumably data entry error or spreadsheet got trashed somewhere
  which.cols = which.cols[!is.na(which.cols)]
  which.cols.mod = which.cols.mod[!is.na(which.cols.mod)]

  datSNP = datCombinedTurquoise1imp[,c(which.cols,which.cols.mod)]
  datExpr = datCombinedTurquoise1imp[,1:2]

  datCombined=data.frame(cbind(datSNP,datExpr))
#  sma = single.marker.analysis(datCombined,snpcols=1:ncol(datSNP),genecols=ncol(datSNP)+2,traitcols=ncol(datSNP)+1)

#  write.csv(sma,file="c:/good/tova.als/sma.trait.priority.csv")
#  topKmarkers = 10

  ### my.preferred.markers = colnames(datSNP)[sort(sma$leo.nb,index.return=T,decreasing=T)$ix[1:topKmarkers]]
#  my.preferred.markers = colnames(datSNP)[1:topKmarkers]

  #  good.sma.snps = colnames(datSNP)[sma$leo.nb > .0]

  #  both.positive.and.top.ranked.leo.nb.markers = intersect(good.sma.snps, my.preferred.markers)


#  best20="rs764571.chr11.bp91192772"
#  n1=neo(data.frame(cbind(datSNP[,best20],datExpr)),snpcols=1)

  r=list()
  r$greedy = NEOrobustnessGreedy( datSNP, datExpr, TopGreedySNPs=1:20)
  r$forward = NEOrobustnessForward(datSNP, datExpr, TopForwardSNPs=1:20)
#  r$greedy.and.forward  = NEOrobustnessForwardAndGreedy(datSNP, datExpr, TopForwardAndGreedySNPs=1:20)
  r
}

null.dist=function(){

 a=load("r.test.robust.results.20.Mon_Jul_16_12.51.12_2007.rdat")

 gre.mean = r$greedy[[1]]
 for.mean = r$greedy[[1]]
 gf.mean  = r$greedy[[1]]

 gre.mean[]=0
 for.mean[]=0
 gf.mean[] =0

 gre.sd = gre.mean
 for.sd = gre.mean
 gf.sd  = gre.mean

  do.list=function(the.mean, the.sd, the.src) {
   for (j in 1:8) {
      for(co in 1:6) {
          tmp = rep(0,20)
          for (i in 1:20) {
              tmp[i]=the.src[[i]][j,co]
          }
          the.mean[j,co] = mean(tmp)
          the.sd[j,co]   = sd(tmp)
      }
   }
   a=list()
   a$the.mean = the.mean
   a$the.sd = the.sd
   a
}

 do.list(gre.mean, gre.sd, r$greedy)

}

# =================================================
# begin recall.barplots() - a long function not meant
# to be executed but to hold all the functions
# and work that went into comparing power and
# false positive rates on the single and orthomarker
# models.
# =================================================

## to read an outputted LabelledMatrix into R
## that is of size (147 x something), use:
## ## to read in data into R that is (147 x something):
## r=read.table(file="/home/jaten/dev/bionumerics/src/bngen/cart_test/qtl_test_raw_uarray_transposed_sex1.dat",header=T,row.names=1,nrows=147,skip=2,sep="\t",na.strings = "nan")
##
## calibration of file load times:
## on iggi: a 293x24967 matrix was 71MB matrix output text file was read into R in 22minutes using this command.
## once saved from R it loads/saves quickly and takes up 56MB on disk.

recall.barplots = function() {

#load("scores.new.rdat") # scores.new
#load("scores.extended.scores.rdat") # ext.scores
#load("scores.extended.goodorder.rdat") # ext.scores
#load("fullscore.rdat")
load("scores424.rdat") # scores424
load("scores.new.for.omega.zero.ha.ha.gamma2.vary.rdat") # scores.new, null.param

# add the splitter variable
# w = 1:nrow(scores424)
# u8=unique(scores424[w,1:8])
# num.reps = length(w)/nrow(u8); # automatically determine num.reps, rather than assume 50
# num.grp = nrow(u8)

# splitter = (0:(length(w)-1)) %/% num.reps;
# splitter = splitter+1 
# scores424$splitter = splitter
# save(scores424, file="scores424.rdat")
# above doesn't need repeating!

# the rows for splitter==1 are 1:50
# the rows for splitter==2 are 51:100
# ...
# the rows for splitter==k are (50*(k-1)+1):(50*k)

if (exists("splitter.rows")) rm(splitter.rows)
splitter.rows=function(k,groupsize=50) { (groupsize*(k-1)+1):(groupsize*k) }

load("scores424.rdat") # scores424

#load("scores.new.rdat")
#> dim(scores.new)
#[1] 5600   64
#load("scores.extended.goodorder.rdat")
#> dim(ext.scores)
#[1] 6550   64

#scores2=rbind(scores.new, ext.scores)
#rownames(scores2)=as.character(1:nrow(scores2))
#fullscore=scores2
#save(fullscore,file="fullscore.rdat")

rownames(ext.scores)=as.character(1:nrow(ext.scores))
unique(ext.scores[,1:8])
##         ha    hb gamma2 omega    eA    eB rho.ab percent.max.omega
## 1    0.467 0.467  0.139 0.463 0.051 0.394  0.602             0.667

## 51   0.467 0.467  0.000 0.579 0.198 0.533  0.579             0.833
## 101  0.467 0.467  0.023 0.579 0.148 0.510  0.602             0.833
## 151  0.467 0.467  0.046 0.579 0.099 0.487  0.625             0.833
## 201  0.467 0.467  0.069 0.579 0.049 0.464  0.648             0.833

# do a fixed signal strength ha, ha equal and fixed, and compare
# rmsea values at different omega values.

ha6hb6.gamma0=c(5301:5350,5501:5550,5701:5750,5901:5950,6101:6150,6301:6350,6501:6550)
freq.bar(scores424[ha6hb6.gamma0,],threshold=.1,th.col="rmsea.f2.M.M1M2.BtoA",group.col="omega",plot.freq.below.th=FALSE,ylim=c(0,.1))

# compute P-RMSEA too
scores424$P.rmsea.f1ab.AtoB = scores424$rmsea.f1ab.AtoB * scores424$ab.path.pv.f1ab.AtoB
scores424$P.rmsea.f1ab.BtoA = scores424$rmsea.f1ab.BtoA * scores424$ab.path.pv.f1ab.BtoA
scores424$P.rmsea.f1ab.conf = scores424$rmsea.f1ab.conf * scores424$ab.path.pv.f1ab.conf
scores424$P.rmsea.f1ab.AcollideB = scores424$rmsea.f1ab.AcollideB * scores424$ab.path.pv.f1ab.AcollideB
scores424$P.rmsea.f1ab.BcollideA = scores424$rmsea.f1ab.BcollideA * scores424$ab.path.pv.f1ab.BcollideA
scores424$P.rmsea.f1ba.AtoB = scores424$rmsea.f1ba.AtoB * scores424$ab.path.pv.f1ba.AtoB
scores424$P.rmsea.f1ba.BtoA = scores424$rmsea.f1ba.BtoA * scores424$ab.path.pv.f1ba.BtoA
scores424$P.rmsea.f1ba.conf = scores424$rmsea.f1ba.conf * scores424$ab.path.pv.f1ba.conf
scores424$P.rmsea.f1ba.AcollideB = scores424$rmsea.f1ba.AcollideB * scores424$ab.path.pv.f1ba.AcollideB
scores424$P.rmsea.f1ba.BcollideA = scores424$rmsea.f1ba.BcollideA * scores424$ab.path.pv.f1ba.BcollideA
scores424$P.rmsea.f2.M.M1M2.AtoB = scores424$rmsea.f2.M.M1M2.AtoB * scores424$ab.path.pv.f2.M.M1M2.AtoB
scores424$P.rmsea.f2.M.M1M2.BtoA = scores424$rmsea.f2.M.M1M2.BtoA * scores424$ab.path.pv.f2.M.M1M2.BtoA
scores424$P.rmsea.f2.M.M1M2unresolv = scores424$rmsea.f2.M.M1M2unresolv * scores424$ab.path.pv.f2.M.M1M2unresolv
scores424$P.rmsea.f2.M.M1M2hidden.con = scores424$rmsea.f2.M.M1M2hidden.con * scores424$ab.path.pv.f2.M.M1M2hidden.con


# summarize null model  -- at different heritability levels
rows.null = which(scores424$gamma2==0 & scores424$omega==0)
rows.null2= which(scores424$gamma2==0 & scores424$omega==0 & scores424$ha != .467)

a=apply(scores424[rows.null,],2,mean)
b=apply(scores424[rows.null,],2,stderr1)
d=rbind(a,a-1.96*b,a+1.96*b,b)
rownames(d)=c("mean","ci.l","ci.u","stderr")
names(a)=colnames(scores424)


  if(exists("m1.model.bar.plot") ) rm(m1.model.bar.plot);
  m1.model.bar.plot=function(my.index, my.y.lab,d,main,do.postscript=FALSE,psfile="bargraph.file.ps",x.h.add=.02,...) {
    x=d[1,c(my.index)]
    x.l=d[2,c(my.index)]
    x.h=d[3,c(my.index)]
    x.se=d[4,c(my.index)]
    names(x)=substr(names(x),7,50)
    
    names(x)=names(x.l)=names(x.h)=names(x.se)=c("(1) M=>A=>B","(2) M=>B=>A","(3) A<=M=>B","(4) M=>A<=B","(5) M=>B<=A")

    if (do.postscript) postscript(file=psfile,horizontal=FALSE);
    mp <- barplot2(x,
                   xpd=FALSE, # bars stay inside region
                   beside = TRUE,
                                        #        col = "yellow",
                   col = c("white","lightblue","tan","brown","black"),
#                   density=c(-1,10,100,-1,-1),
                                        #        legend = 
                   main = main,
                   font.main = 4,
                   sub = "One marker models. 450 simulations of N=200. Mean (SE) plotted.\n",
                   ylab= my.y.lab,
                                        #        col.sub = "yellow",
                   cex.names = .8,
                   plot.ci = TRUE,
                   ci.l = x.l, 
                   ci.u = x.h, 
                   plot.grid = TRUE,...)
    
    for (i in 1:length(x)) {
      text(x=mp[i],y=x.h[i]+4*(x.h-x)[i],labels=paste(sep="",signif(x[i],2)," (",signif(x.se[i],2),")"))
    }
    box()
    if (do.postscript) dev.off()
  } # end model.bar.plot

#
# m1.model.bar.plot.rmsea.plus
#
#
#
#
#
#
#
#
# define a version that compares RMSEA to RMSEA + BilayerVerify (earlier known as zeo2) + PathAbsZ, the p-value
# for the path coefficient for the models with A-B paths.
# 
# for one marker f1ab
cols.zeo=c(75:76)
cols.path.abs.z=c(65:68)

bpo=c("rmsea.f1ab.AtoB","chisq.f1ab.AtoB","path.abs.z.f1ab.AtoB","rmsea.f1ab.BtoA","chisq.f1ab.BtoA","path.abs.z.f1ab.BtoA","rmsea.f1ab.AcollideB","path.abs.z.f1ab.AcollideB","chisq.f1ab.AcollideB","rmsea.f1ab.BcollideA","path.abs.z.f1ab.BcollideA","chisq.f1ab.BcollideA")

#bpo.conf=c("rmsea.f1ab.conf","chisq.f1ab.conf","zeo.ma.b.given.a","zeo.mb.a.given.b","cor.ma.b","pcor.ma.b.given.a","cor.mb.a","pcor.mb.a.given.b")

# take out the chisq
bpo.conf=c("rmsea.f1ab.conf","zeo.ma.b.given.a","zeo.mb.a.given.b","cor.ma.b","pcor.ma.b.given.a","cor.mb.a","pcor.mb.a.given.b")


if(exists("m1.model.bar.plot.rmsea.plus") ) rm(m1.model.bar.plot.rmsea.plus);
m1.model.bar.plot.rmsea.plus=function(my.index, cols.path.abs.z, cols.zeo, my.y.lab,d,main,do.postscript=FALSE,psfile="bargraph.file.ps",x.h.add=.02,...) {
    x=d[1,my.index]
    x.l=d[2,my.index]
    x.h=d[3,my.index]
    x.se=d[4,my.index]

#    names(x)=names(x.l)=names(x.h)=names(x.se)=c("(1) M=>A=>B","(2) M=>B=>A","(3) A<=M=>B","(4) M=>A<=B","(5) M=>B<=A")

#   barplot(split(log10(d$leo.i.correct[all]),d$ntrack[all]),notch=TRUE,ylab="RMSEA",at = 1:length(x) - 0.2,boxwex = 0.2,col = "yellow")
    boxplot(split(log10(d$leo.i[all]),d$ntrack[all]),notch=TRUE,xlab="N.observations",ylab="log10(LEO.I)",at = 1:length(unique(d$ntrack[all])) + 0.2,add=TRUE,boxwex = 0.25,col = "orange",show.names=FALSE)#,ylim=c(-20,20))


    if (do.postscript) postscript(file=psfile,horizontal=FALSE);
    mp <- barplot2(x,
                   xpd=FALSE, # bars stay inside region
                   beside = TRUE, # don't stack columns, put them alongside each other
                                        #        col = "yellow",
                   col = c("white","lightblue","tan","brown","black"),
#                   density=c(-1,10,100,-1,-1),
                                        #        legend = 
                   main = main,
                   font.main = 4,
                   sub = "10 simulations of N=200. Mean (SE) plotted with 1.96*SE error bars.\n",
                   ylab= my.y.lab,
                                        #        col.sub = "yellow",
                   cex.names = .6,
                   plot.ci = TRUE,
                   ci.l = x.l,
                   ci.u = x.h,
                   plot.grid = TRUE)

    mp <- barplot2(x,add=TRUE,
                   xpd=FALSE, # bars stay inside region
                   beside = TRUE,
                                        #        col = "yellow",
                   col = c("white","lightblue","tan","brown","black"),
#                   density=c(-1,10,100,-1,-1),
                                        #        legend = 
                   main = main,
                   font.main = 4,
                   sub = "One marker models. 450 simulations of N=200. Mean (SE) plotted.\n",
                   ylab= my.y.lab,
                                        #        col.sub = "yellow",
                   cex.names = .8,
                   plot.ci = TRUE,
                   ci.l = x.l, 
                   ci.u = x.h, 
                   plot.grid = TRUE)
    
    for (i in 1:length(x)) {
      text(x=mp[i],y=x.h[i]+4*(x.h-x)[i],labels=paste(sep="",signif(x[i],2)," (",signif(x.se[i],2),")"))
    }
    box()
    if (do.postscript) dev.off()
  } # end model.bar.plot

 # for the null.param, scores.new under the null comparison
  m1.model.bar.plot(my.index=rmsea.f1ab, my.y.lab="RMSEA", d,main = "Monte Carlo study of the Null model: M=>A, while B independent",do.postscript=FALSE,psfile="mc.null.1ab.rmsea.ps",ylim = c(0, .8),density=c(-1,10,100,-1,-1));



rmsea.f1ab=c(9:13)
rmsea.f1ba=c(14:18)

  m1.model.bar.plot(my.index=rmsea.f1ab, my.y.lab="RMSEA", d,main = "Monte Carlo study of the Null model: M=>A, while B independent",do.postscript=TRUE,psfile="mc.null.1ab.rmsea.ps",ylim = c(0, .8),density=c(-1,10,100,-1,-1));
  m1.model.bar.plot(my.index=rmsea.f1ba, my.y.lab="RMSEA", d,main = "Monte Carlo study of the Null model: M=>B, while A independent",do.postscript=TRUE,psfile="mc.null.1ba.rmsea.ps",ylim = c(0, .8),density=c(-1,10,100,-1,-1));
  m1.model.bar.plot(my.index=51:55, my.y.lab="-log10 Probability(>chi-squared)", d,main = "Monte Carlo study of the Null model: M=>A, while B independent",do.postscript=T,psfile="mc.null.1ab.chiprob.ps",ylim = c(.1, 50),log="y");#,density=c(-1,10,100,-1,-1));
  m1.model.bar.plot(my.index=56:60, my.y.lab="-log10 Probability(>chi-squared)", d,main = "Monte Carlo study of the Null model: M=>B, while A independent",do.postscript=T,psfile="mc.null.1ba.chiprob.ps",ylim = c(.1, 50),log="y");#,density=c(-1,10,100,-1,-1));

  chisq.f1ab=c(37:41)
  chisq.f1ba=c(42:46)

  rmsea.f2=c(19:22)

  m1.model.bar.plot(my.index=chisq.f1ab, my.y.lab="Chi-squared statistic, 1 df", d,main = "Monte Carlo study of the Null model: M=>A, while B independent",do.postscript=TRUE,psfile="mc.null.1ab.chisq.ps",ylim=c(.1,200),x.h.add=.6,log="y");
    m1.model.bar.plot(my.index=chisq.f1ba, my.y.lab="Chi-squared statistic, 1 df", d,main = "Monte Carlo study of the Null model: M=>B, while A independent",do.postscript=TRUE,psfile="mc.null.1ba.chisq.ps",ylim=c(.1,200),log="y");

#### The four variable models, under the null.


  if(exists("m2.model.bar.plot") ) rm(m2.model.bar.plot);
  m2.model.bar.plot=function(my.index, my.y.lab,d,main,do.postscript=FALSE,psfile="bargraph.file.ps",...) {
    x=d[1,c(my.index)]
    x.l=d[2,c(my.index)]
    x.h=d[3,c(my.index)]
    x.se=d[4,c(my.index)]
    names(x)=substr(names(x),7,50)
    
    names(x)=names(x.l)=names(x.h)=names(x.se)=c("(1) MA=>A=>B<=MB","(2) MA=>A<=B<=MB","(3) MA=>{A,B};MB=>{A,B}","(4) MA=>A<=(H)=>B<=MB")

    if (do.postscript) postscript(file=psfile,horizontal=FALSE);
    mp <- barplot2(x,
                   xpd=FALSE, # bars stay inside region
                   beside = TRUE,
                                        #        col = "yellow",
                   col = c("white","lightblue","tan","brown","black"),
#                   density=c(-1,10,100,-1,-1),
                                        #        legend = 
                   main = main,
                   font.main = 4,
                   sub = "Two marker models. 450 simulations of N=200. Mean (SE) plotted.\n",
                   ylab= my.y.lab,
                                        #        col.sub = "yellow",
                   cex.names = .7,
                   plot.ci = TRUE,
                   ci.l = x.l, 
                   ci.u = x.h, 
                   plot.grid = TRUE,...)
    
    for (i in 1:length(x)) {
      text(x=mp[i],y=x.h[i]+(x.h-x)[i],labels=paste(sep="",signif(x[i],2)," (",signif(x.se[i],2),")"))
    }
    box()
    if (do.postscript) dev.off()
  } # end model.bar.plot

 
 m2.model.bar.plot(d=d,my.index=19:22, my.y.lab="RMSEA", main="Monte Carlo study of Null model: MA=>A; MB=>B; A and B indep.",ylim=c(0,.04),do.postscript=T,psfile="mc.null.2m.rmsea.ps")
 m2.model.bar.plot(d=d,my.index=47:50, my.y.lab="Chi-squared statistic, df={3,3,2,3}", main="Monte Carlo study of Null model: MA=>A; MB=>B; A and B indep.",ylim=c(0,5),do.postscript=T,psfile="mc.null.2m.chisq.ps")
 m2.model.bar.plot(d=d,my.index=61:64, my.y.lab="-log10 Probability(>chi-squared)", main="Monte Carlo study of Null model: MA=>A; MB=>B; A and B indep.",ylim=c(0,.6),do.postscript=T,psfile="mc.null.2m.chiprob.ps")

##############################
# begin power graphs -- with actual signal


if(exists("power.graph") ) rm(power.graph);
power.graph=function(ha,hb,scores424,do.postscript=FALSE,do.labels=TRUE,...) {
 
 # pick out an ha and hb level with w.
 if (do.postscript) { 
     psfile=paste(sep="","power.plot.ha",round(ha,2),".hb",round(hb,2),"..rhoab.vs.leo.i.ps");
     postscript(file=psfile,horizontal=FALSE) 
 } else {
    # windows()
 }

 w=which(scores424$ha==ha & scores424$hb==hb)

 u8=unique(scores424[w,1:8])
 num.reps = length(w)/nrow(u8); # automatically determine num.reps, rather than assume 50
 num.grp = nrow(u8)

 #plot(scores424$rho.ab[w], scores424$rmsea.f2.M.M1M2.BtoA[w])


 leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])

 # split leo.i scores by the different parameter variations
 #bplist=split(leo.i,scores424$rho.ab[w])
 bplist=split(leo.i,scores424$splitter[w])

 group.firsts= seq(1,length(w),num.reps)

 bps=rep(NA,length(bplist));
 rhos=rep(NA,length(bplist));
 for (i in 1:length(bplist)) {
    bps[i]=boxplot.stats(bplist[[i]])$stats[3]
    rhos[i]=scores424$rho.ab[w[group.firsts[i]]]
 }

 scale.factor=num.grp / max(rhos) 

 # manually split out the levels, so we can be sure the
 # split is correct (not combining close or duplicated rho levels), and
 # then label the split groups correctly with their rho.ab.level, rhos.
 man.split=split(leo.i,scores424$splitter[w])
 names(man.split)=as.character(rhos)

 if (do.labels) {
#  boxplot(leo.i ~ scores424$rho.ab[w],notch=F,ylab="LEO.I for MA,MB model",xlab=paste(sep="","Horizontal axis: observed rho(A,B) = gamma^2 + omega. hA=",round(ha,2),", hB=",round(hb,2),". Dashed lines at LEO.I=2"),col="yellow",at=rhos*scale.factor,main="Monte Carlo study: detecting causal influence under partial\nconfounding: 0%, 17%, 33%, 50%, 67%, 83%, 100% causal",sub="Group labels show (omega/rho(A,B)) or percent true causal influence. Boxes for 50 replicates at N=200.")#,at=1:length(unique(scores424$percent.max.omega))-.2,boxwex=.25,ylim=c(-.15,.7))

  boxplot(man.split,notch=F,ylab="LEO.I for MA,MB model",xlab=paste(sep="","Horizontal axis: observed rho(A,B) = gamma^2 + omega. hA=",round(ha,2),", hB=",round(hb,2),". Dashed lines at LEO.I=2"),col="yellow",at=rhos*scale.factor,main="Monte Carlo study: detecting causal influence under partial\nconfounding: 0%, 17%, 33%, 50%, 67%, 83%, 100% causal",sub="Group labels show (omega/rho(A,B)) or percent true causal influence. Boxes for 50 replicates at N=200.",...)#,at=1:length(unique(scores424$percent.max.omega))-.2,boxwex=.25)

 } else {
  boxplot(man.split,notch=F,col="yellow",at=rhos*scale.factor,ylab="LEO.I for MA, MB model",sub=paste(sep="","Observed correlation(A,B)\nhA=",ha,",  hB=",hb),...)#,at=1:length(unique(scores424$percent.max.omega))-.2,boxwex=.25,ylim=c(-20,20))
 }


 # hit the center of the boxplots like this:
 points(x=rhos*scale.factor,y=bps,pch=19,col="red")

 # connect the centers that have the same percent.max.omega

 # split by groups that have the same percent.max.omega
 percent.omega.groups = split(w,scores424$percent.max.omega[w])
 ngroup = length(percent.omega.groups)
 
 # get the labels for the groups
 pog.nm=as.character(round(as.numeric(names(percent.omega.groups)),2))

 # omit last group: only 1 point
 for (j in 1:ngroup) {
 
 w.grp=percent.omega.groups[[j]]

 leo.i.grp = (scores424$chiprob.f2.M.M1M2hidden.con[w.grp]  - scores424$chiprob.f2.M.M1M2.BtoA[w.grp])

 rhos.grp=sort(unique(scores424$rho.ab[w.grp]))

 bplist.grp=split(leo.i.grp,scores424$rho.ab[w.grp])
 bps.grp=rep(NA,length(bplist.grp));

 bn = names(bplist.grp)

 for (i in 1:length(bplist.grp)) {
    bps.grp[i]=boxplot.stats(bplist.grp[[i]])$stats[3]
 }

 # take the max point in the first group as the reference location to label above
   bps.grp.top=max(bplist.grp[[1]])
   text(x=rhos.grp[1]*scale.factor,y=bps.grp.top+2,labels=pog.nm[j])

 # if (j != ngroup) { # don't plot lines for last group of singleton point.
    lines(x=rhos.grp*scale.factor,y=bps.grp)
 # }
} # end j

abline(a=2,b=0,lty=2)
#abline(v=.3*scale.factor,lty=3)  

  if (do.postscript) { 
     dev.off(); 
     system(paste(sep="","ps2pdf ",psfile),wait=T);
     system(paste(sep="","chmod +x ",paste(sep="",substr(psfile,1,(nchar(psfile,type="chars")-3)),".pdf")),wait=F);
  }
 invisible(0)
} # end power.graph

for (ha in unique(scores424$ha)) {
  for (hb in unique(scores424$hb)) {
    power.graph(ha, hb, scores424=scores424,do.postscript=T)
  }}

k=0
for (ha in unique(scores424$ha)) {
  for (hb in unique(scores424$hb)) {
    if (k %% 4 == 0) { windows(); par(mfrow=c(2,2)); }
    power.graph(ha, hb, scores424=scores424,do.postscript=F,do.labels=FALSE,ylim=c(-20,20),sub="observed correlation(A,B)")
    k=k+1
  }
 par(mfrow=c(1,1)); title(paste(sep="","Monte Carlo study: the effect of direct genetic\ninfluence on power to detect causal influence. hA=",round(ha,2)));
}

# version to create output postscript files
k=0
for (ha in unique(scores424$ha)) {
  postscript(file=paste(sep="","fourplex.power.study.ha.",round(ha,2),".ps"),horizontal=FALSE); 
  par(mfrow=c(2,2));

  for (hb in unique(scores424$hb)) {
    power.graph(ha, hb, scores424=scores424,do.postscript=F,do.labels=FALSE,ylim=c(-20,20),sub="observed correlation(A,B)")
    k=k+1
  }
 par(mfrow=c(1,1)); 
 title(paste(sep="","Monte Carlo study: the effect of direct genetic\ninfluence on power to detect causal influence. hA=",round(ha,2)));
 dev.off()
}


# end power graphs
###############################

xtabs(chiprob.f1ab.conf ~ ha + hb,data=scores424[setdiff(rows.null,which(scores424$ha == .467)),])

dataIN=data.frame(ext.scores2$percent.max.omega, ext.scores2$rmsea.f2.M.M1M2.BtoA)
colnames(dataIN)=c("percent.max.omega","RMSEA.M1M2.BtoA")
rmsea.threshold=.05

# Hint: choose the threshold so that the POWER for omega=0 should is smaller than 0.01.

# assign a 1 if the rmsea is below threshold: make an indicator whose mean is the frequency of cases below the Threshold.
LEO.Icut=ifelse(dataIN[,2]<rmsea.threshold,1,0)

# group data by percent.max.omega, then get the frequency of rmsea being good in each group
# ; the mean or expected value of an indicator is just the probability or frequency of that
# indicator.

# tapply automatically sorts the percent.max.omega groups into ascending order.
POWER=as.vector(tapply(LEO.Icut,dataIN$percent.max.omega,mean));
group.names=round(sort(unique(dataIN$percent.max.omega)),2)
names(POWER)=as.character(group.names)

SE.POWER= as.vector(tapply(LEO.Icut,dataIN$percent.max.omega,stderr1))

ci.l = POWER-SE.POWER
ci.u = POWER+SE.POWER


par(mfrow=c(1,1))

#barplot(POWER,names.arg=names(table(dataIN$percent.max.omega) ), ylab="POWER",main=paste("No. Samples=",50, "RMSEA Threshold", rmsea.threshold))
#err.bp(as.vector(POWER), as.vector(1.96*SE.POWER), two.side=T)

# Produce the above barplots for different choices of the threshold

   windows()

#ci.l=rbind(ci.l,ci.l)
#ci.u=rbind(ci.u,ci.u)
#POWER2=rbind(POWER,POWER)

   mp <- barplot2(POWER,
             beside = TRUE,
#        col = "yellow",
        col = c("yellow", "violet"),
#                "lightcyan", "orange"),
#        legend = 
                  ylim = c(0, 1),
#        main = "Death Rates in Virginia",
        font.main = 4,
#        sub = "Faked 95 percent error bars",
#        col.sub = "yellow",
        cex.names = 1.25,
        plot.ci = TRUE,
               ci.l = ci.l,
               ci.u = ci.u,
        plot.grid = TRUE)



library(gplots)

 # Example with confidence intervals and grid
 hh <- t(VADeaths)[, 5:1] # each row is a group with a different color bar
                          # each column is a group that is spaced differently on the x axis


 mybarcol <- "gray20"
 ci.l <- hh * 0.85
 ci.u <- hh * 1.15
mp <- barplot2(hh, beside = TRUE,
        col = c("yellow", "violet",
                "lightcyan", "orange"),
        legend = colnames(VADeaths), ylim = c(0, 100),
        main = "Death Rates in Virginia", font.main = 4,
        sub = "Faked 95 percent error bars", col.sub = mybarcol,
        cex.names = 1.5, plot.ci = TRUE, ci.l = ci.l, ci.u = ci.u,
        plot.grid = TRUE)
 mtext(side = 1, at = colMeans(mp), line = 2,
      text = paste("Mean", formatC(colMeans(hh))), col = "red")
 box()

## > unique(scores424[,1:8])
##          ha    hb gamma2 omega    eA    eB rho.ab percent.max.omega
## 1     0.200 0.200  0.000 0.000 0.800 0.800  0.000             0.000
## 51    0.200 0.200  0.200 0.000 0.600 0.600  0.200             0.000
## 101   0.200 0.200  0.400 0.000 0.400 0.400  0.400             0.000
## 151   0.200 0.200  0.600 0.000 0.200 0.200  0.600             0.000

## 201   0.200 0.200  0.000 0.117 0.786 0.800  0.117             0.167
## 251   0.200 0.200  0.199 0.117 0.541 0.601  0.316             0.167
## 301   0.200 0.200  0.398 0.117 0.296 0.402  0.515             0.167
## 351   0.200 0.200  0.597 0.117 0.050 0.203  0.714             0.167

## 401   0.200 0.200  0.000 0.233 0.746 0.800  0.233             0.333
## 451   0.200 0.200  0.158 0.233 0.514 0.642  0.391             0.333
## 501   0.200 0.200  0.316 0.233 0.282 0.484  0.549             0.333
## 551   0.200 0.200  0.474 0.233 0.050 0.326  0.707             0.333

## 601   0.200 0.200  0.000 0.350 0.677 0.800  0.350             0.500
## 651   0.200 0.200  0.123 0.350 0.468 0.677  0.473             0.500
## 701   0.200 0.200  0.246 0.350 0.259 0.554  0.596             0.500
## 751   0.200 0.200  0.369 0.350 0.050 0.431  0.719             0.500

## 801   0.200 0.200  0.000 0.467 0.582 0.800  0.467             0.667
## 851   0.200 0.200  0.092 0.467 0.404 0.708  0.559             0.667
## 901   0.200 0.200  0.183 0.467 0.228 0.617  0.650             0.667
## 951   0.200 0.200  0.275 0.467 0.051 0.525  0.742             0.667

## 1001  0.200 0.200  0.000 0.583 0.460 0.800  0.583             0.833
## 1051  0.200 0.200  0.063 0.583 0.323 0.737  0.646             0.833
## 1101  0.200 0.200  0.126 0.583 0.187 0.674  0.709             0.833
## 1151  0.200 0.200  0.189 0.583 0.050 0.611  0.772             0.833

## 1201  0.200 0.200  0.000 0.700 0.310 0.800  0.700             1.000
## 1251  0.200 0.200  0.036 0.700 0.224 0.764  0.736             1.000
## 1301  0.200 0.200  0.072 0.700 0.137 0.728  0.772             1.000
## 1351  0.200 0.200  0.108 0.700 0.051 0.692  0.808             1.000

## 1401  0.200 0.333  0.000 0.000 0.800 0.667  0.000             0.000
## 1451  0.200 0.333  0.200 0.000 0.600 0.467  0.200             0.000
## 1501  0.200 0.333  0.400 0.000 0.400 0.267  0.400             0.000
## 1551  0.200 0.333  0.600 0.000 0.200 0.067  0.600             0.000

## 1601  0.200 0.333  0.000 0.117 0.786 0.667  0.117             0.167
## 1651  0.200 0.333  0.199 0.117 0.541 0.468  0.316             0.167
## 1701  0.200 0.333  0.398 0.117 0.296 0.269  0.515             0.167
## 1751  0.200 0.333  0.597 0.117 0.050 0.070  0.714             0.167

## 1801  0.200 0.333  0.000 0.233 0.746 0.667  0.233             0.333
## 1851  0.200 0.333  0.158 0.233 0.514 0.509  0.391             0.333
## 1901  0.200 0.333  0.316 0.233 0.282 0.351  0.549             0.333
## 1951  0.200 0.333  0.474 0.233 0.050 0.193  0.707             0.333

## 2001  0.200 0.333  0.000 0.350 0.677 0.667  0.350             0.500
## 2051  0.200 0.333  0.123 0.350 0.468 0.544  0.473             0.500
## 2101  0.200 0.333  0.246 0.350 0.259 0.421  0.596             0.500
## 2151  0.200 0.333  0.369 0.350 0.050 0.298  0.719             0.500

## 2201  0.200 0.333  0.000 0.467 0.582 0.667  0.467             0.667
## 2251  0.200 0.333  0.092 0.467 0.404 0.575  0.559             0.667
## 2301  0.200 0.333  0.183 0.467 0.228 0.484  0.650             0.667
## 2351  0.200 0.333  0.275 0.467 0.051 0.392  0.742             0.667

## 2401  0.200 0.333  0.000 0.583 0.460 0.667  0.583             0.833
## 2451  0.200 0.333  0.063 0.583 0.323 0.604  0.646             0.833
## 2501  0.200 0.333  0.126 0.583 0.187 0.541  0.709             0.833
## 2551  0.200 0.333  0.189 0.583 0.050 0.478  0.772             0.833

## 2601  0.200 0.333  0.000 0.700 0.310 0.667  0.700             1.000
## 2651  0.200 0.333  0.036 0.700 0.224 0.631  0.736             1.000
## 2701  0.200 0.333  0.072 0.700 0.137 0.595  0.772             1.000
## 2751  0.200 0.333  0.108 0.700 0.051 0.559  0.808             1.000

## 2801  0.200 0.467  0.000 0.000 0.800 0.533  0.000             0.000
## 2851  0.200 0.467  0.161 0.000 0.639 0.372  0.161             0.000
## 2901  0.200 0.467  0.322 0.000 0.478 0.211  0.322             0.000
## 2951  0.200 0.467  0.483 0.000 0.317 0.050  0.483             0.000

## 3001  0.200 0.467  0.000 0.117 0.786 0.533  0.117             0.167
## 3051  0.200 0.467  0.161 0.117 0.588 0.372  0.278             0.167
## 3101  0.200 0.467  0.322 0.117 0.389 0.211  0.439             0.167
## 3151  0.200 0.467  0.483 0.117 0.191 0.050  0.600             0.167

## 3201  0.200 0.467  0.000 0.233 0.746 0.533  0.233             0.333
## 3251  0.200 0.467  0.158 0.233 0.514 0.375  0.391             0.333
## 3301  0.200 0.467  0.316 0.233 0.282 0.217  0.549             0.333
## 3351  0.200 0.467  0.474 0.233 0.050 0.059  0.707             0.333

## 3401  0.200 0.467  0.000 0.350 0.677 0.533  0.350             0.500
## 3451  0.200 0.467  0.123 0.350 0.468 0.410  0.473             0.500
## 3501  0.200 0.467  0.246 0.350 0.259 0.287  0.596             0.500
## 3551  0.200 0.467  0.369 0.350 0.050 0.164  0.719             0.500

## 3601  0.200 0.467  0.000 0.467 0.582 0.533  0.467             0.667
## 3651  0.200 0.467  0.092 0.467 0.404 0.441  0.559             0.667
## 3701  0.200 0.467  0.183 0.467 0.228 0.350  0.650             0.667
## 3751  0.200 0.467  0.275 0.467 0.051 0.258  0.742             0.667

## 3801  0.200 0.467  0.000 0.583 0.460 0.533  0.583             0.833
## 3851  0.200 0.467  0.063 0.583 0.323 0.470  0.646             0.833
## 3901  0.200 0.467  0.126 0.583 0.187 0.407  0.709             0.833
## 3951  0.200 0.467  0.189 0.583 0.050 0.344  0.772             0.833

## 4001  0.200 0.467  0.000 0.700 0.310 0.533  0.700             1.000
## 4051  0.200 0.467  0.036 0.700 0.224 0.497  0.736             1.000
## 4101  0.200 0.467  0.072 0.700 0.137 0.461  0.772             1.000
## 4151  0.200 0.467  0.108 0.700 0.051 0.425  0.808             1.000

## 4201  0.200 0.600  0.000 0.000 0.800 0.400  0.000             0.000
## 4251  0.200 0.600  0.117 0.000 0.683 0.283  0.117             0.000
## 4301  0.200 0.600  0.233 0.000 0.567 0.167  0.233             0.000
## 4351  0.200 0.600  0.350 0.000 0.450 0.050  0.350             0.000

## 4401  0.200 0.600  0.000 0.117 0.786 0.400  0.117             0.167
## 4451  0.200 0.600  0.117 0.117 0.642 0.283  0.234             0.167
## 4501  0.200 0.600  0.233 0.117 0.499 0.167  0.350             0.167
## 4551  0.200 0.600  0.350 0.117 0.355 0.050  0.467             0.167

## 4601  0.200 0.600  0.000 0.233 0.746 0.400  0.233             0.333
## 4651  0.200 0.600  0.117 0.233 0.574 0.283  0.350             0.333
## 4701  0.200 0.600  0.233 0.233 0.404 0.167  0.466             0.333
## 4751  0.200 0.600  0.350 0.233 0.232 0.050  0.583             0.333

## 4801  0.200 0.600  0.000 0.350 0.677 0.400  0.350             0.500
## 4851  0.200 0.600  0.117 0.350 0.479 0.283  0.467             0.500
## 4901  0.200 0.600  0.233 0.350 0.281 0.167  0.583             0.500
## 4951  0.200 0.600  0.350 0.350 0.083 0.050  0.700             0.500

## 5001  0.200 0.600  0.000 0.467 0.582 0.400  0.467             0.667
## 5051  0.200 0.600  0.092 0.467 0.404 0.308  0.559             0.667
## 5101  0.200 0.600  0.183 0.467 0.228 0.217  0.650             0.667
## 5151  0.200 0.600  0.275 0.467 0.051 0.125  0.742             0.667

## 5201  0.200 0.600  0.000 0.583 0.460 0.400  0.583             0.833
## 5251  0.200 0.600  0.063 0.583 0.323 0.337  0.646             0.833
## 5301  0.200 0.600  0.126 0.583 0.187 0.274  0.709             0.833
## 5351  0.200 0.600  0.189 0.583 0.050 0.211  0.772             0.833

## 5401  0.200 0.600  0.000 0.700 0.310 0.400  0.700             1.000
## 5451  0.200 0.600  0.036 0.700 0.224 0.364  0.736             1.000
## 5501  0.200 0.600  0.072 0.700 0.137 0.328  0.772             1.000
## 5551  0.200 0.600  0.108 0.700 0.051 0.292  0.808             1.000

## 5601  0.333 0.200  0.000 0.000 0.667 0.800  0.000             0.000
## 5651  0.333 0.200  0.200 0.000 0.467 0.600  0.200             0.000
## 5701  0.333 0.200  0.400 0.000 0.267 0.400  0.400             0.000
## 5751  0.333 0.200  0.600 0.000 0.067 0.200  0.600             0.000

## 5801  0.333 0.200  0.000 0.117 0.653 0.800  0.117             0.167
## 5851  0.333 0.200  0.163 0.117 0.452 0.637  0.280             0.167
## 5901  0.333 0.200  0.326 0.117 0.251 0.474  0.443             0.167
## 5951  0.333 0.200  0.489 0.117 0.050 0.311  0.606             0.167

## 6001  0.333 0.200  0.000 0.233 0.612 0.800  0.233             0.333
## 6051  0.333 0.200  0.128 0.233 0.424 0.672  0.361             0.333
## 6101  0.333 0.200  0.255 0.233 0.238 0.545  0.488             0.333
## 6151  0.333 0.200  0.383 0.233 0.050 0.417  0.616             0.333

## 6201  0.333 0.200  0.000 0.350 0.544 0.800  0.350             0.500
## 6251  0.333 0.200  0.097 0.350 0.379 0.703  0.447             0.500
## 6301  0.333 0.200  0.194 0.350 0.214 0.606  0.544             0.500
## 6351  0.333 0.200  0.291 0.350 0.049 0.509  0.641             0.500

## 6401  0.333 0.200  0.000 0.467 0.449 0.800  0.467             0.667
## 6451  0.333 0.200  0.069 0.467 0.315 0.731  0.536             0.667
## 6501  0.333 0.200  0.137 0.467 0.184 0.663  0.604             0.667
## 6551  0.333 0.200  0.206 0.467 0.051 0.594  0.673             0.667

## 6601  0.333 0.200  0.000 0.583 0.326 0.800  0.583             0.833
## 6651  0.333 0.200  0.043 0.583 0.233 0.757  0.626             0.833
## 6701  0.333 0.200  0.085 0.583 0.142 0.715  0.668             0.833
## 6751  0.333 0.200  0.128 0.583 0.049 0.672  0.711             0.833

## 6801  0.333 0.200  0.000 0.700 0.177 0.800  0.700             1.000
## 6851  0.333 0.200  0.018 0.700 0.133 0.782  0.718             1.000
## 6901  0.333 0.200  0.035 0.700 0.093 0.765  0.735             1.000
## 6951  0.333 0.200  0.053 0.700 0.049 0.747  0.753             1.000

## 7001  0.333 0.333  0.000 0.000 0.667 0.667  0.000             0.000
## 7051  0.333 0.333  0.200 0.000 0.467 0.467  0.200             0.000
## 7101  0.333 0.333  0.400 0.000 0.267 0.267  0.400             0.000
## 7151  0.333 0.333  0.600 0.000 0.067 0.067  0.600             0.000

## 7201  0.333 0.333  0.000 0.117 0.653 0.667  0.117             0.167
## 7251  0.333 0.333  0.163 0.117 0.452 0.504  0.280             0.167
## 7301  0.333 0.333  0.326 0.117 0.251 0.341  0.443             0.167
## 7351  0.333 0.333  0.489 0.117 0.050 0.178  0.606             0.167

## 7401  0.333 0.333  0.000 0.233 0.612 0.667  0.233             0.333
## 7451  0.333 0.333  0.128 0.233 0.424 0.539  0.361             0.333
## 7501  0.333 0.333  0.255 0.233 0.238 0.412  0.488             0.333
## 7551  0.333 0.333  0.383 0.233 0.050 0.284  0.616             0.333

## 7601  0.333 0.333  0.000 0.350 0.544 0.667  0.350             0.500
## 7651  0.333 0.333  0.097 0.350 0.379 0.570  0.447             0.500
## 7701  0.333 0.333  0.194 0.350 0.214 0.473  0.544             0.500
## 7751  0.333 0.333  0.291 0.350 0.049 0.376  0.641             0.500

## 7801  0.333 0.333  0.000 0.467 0.449 0.667  0.467             0.667
## 7851  0.333 0.333  0.069 0.467 0.315 0.598  0.536             0.667
## 7901  0.333 0.333  0.137 0.467 0.184 0.530  0.604             0.667
## 7951  0.333 0.333  0.206 0.467 0.051 0.461  0.673             0.667

## 8001  0.333 0.333  0.000 0.583 0.326 0.667  0.583             0.833
## 8051  0.333 0.333  0.043 0.583 0.233 0.624  0.626             0.833
## 8101  0.333 0.333  0.085 0.583 0.142 0.582  0.668             0.833
## 8151  0.333 0.333  0.128 0.583 0.049 0.539  0.711             0.833

## 8201  0.333 0.333  0.000 0.700 0.177 0.667  0.700             1.000
## 8251  0.333 0.333  0.018 0.700 0.133 0.649  0.718             1.000
## 8301  0.333 0.333  0.035 0.700 0.093 0.632  0.735             1.000
## 8351  0.333 0.333  0.053 0.700 0.049 0.614  0.753             1.000

## 8401  0.333 0.467  0.000 0.000 0.667 0.533  0.000             0.000
## 8451  0.333 0.467  0.161 0.000 0.506 0.372  0.161             0.000
## 8501  0.333 0.467  0.322 0.000 0.345 0.211  0.322             0.000
## 8551  0.333 0.467  0.483 0.000 0.184 0.050  0.483             0.000

## 8601  0.333 0.467  0.000 0.117 0.653 0.533  0.117             0.167
## 8651  0.333 0.467  0.161 0.117 0.454 0.372  0.278             0.167
## 8701  0.333 0.467  0.322 0.117 0.256 0.211  0.439             0.167
## 8751  0.333 0.467  0.483 0.117 0.057 0.050  0.600             0.167

## 8801  0.333 0.467  0.000 0.233 0.612 0.533  0.233             0.333
## 8851  0.333 0.467  0.128 0.233 0.424 0.405  0.361             0.333
## 8901  0.333 0.467  0.255 0.233 0.238 0.278  0.488             0.333
## 8951  0.333 0.467  0.383 0.233 0.050 0.150  0.616             0.333

## 9001  0.333 0.467  0.000 0.350 0.544 0.533  0.350             0.500
## 9051  0.333 0.467  0.097 0.350 0.379 0.436  0.447             0.500
## 9101  0.333 0.467  0.194 0.350 0.214 0.339  0.544             0.500
## 9151  0.333 0.467  0.291 0.350 0.049 0.242  0.641             0.500

## 9201  0.333 0.467  0.000 0.467 0.449 0.533  0.467             0.667
## 9251  0.333 0.467  0.069 0.467 0.315 0.464  0.536             0.667
## 9301  0.333 0.467  0.137 0.467 0.184 0.396  0.604             0.667
## 9351  0.333 0.467  0.206 0.467 0.051 0.327  0.673             0.667

## 9401  0.333 0.467  0.000 0.583 0.326 0.533  0.583             0.833
## 9451  0.333 0.467  0.043 0.583 0.233 0.490  0.626             0.833
## 9501  0.333 0.467  0.085 0.583 0.142 0.448  0.668             0.833
## 9551  0.333 0.467  0.128 0.583 0.049 0.405  0.711             0.833

## 9601  0.333 0.467  0.000 0.700 0.177 0.533  0.700             1.000
## 9651  0.333 0.467  0.018 0.700 0.133 0.515  0.718             1.000
## 9701  0.333 0.467  0.035 0.700 0.093 0.498  0.735             1.000
## 9751  0.333 0.467  0.053 0.700 0.049 0.480  0.753             1.000

## 9801  0.333 0.600  0.000 0.000 0.667 0.400  0.000             0.000
## 9851  0.333 0.600  0.117 0.000 0.550 0.283  0.117             0.000
## 9901  0.333 0.600  0.233 0.000 0.434 0.167  0.233             0.000
## 9951  0.333 0.600  0.350 0.000 0.317 0.050  0.350             0.000

## 10001 0.333 0.600  0.000 0.117 0.653 0.400  0.117             0.167
## 10051 0.333 0.600  0.117 0.117 0.509 0.283  0.234             0.167
## 10101 0.333 0.600  0.233 0.117 0.366 0.167  0.350             0.167
## 10151 0.333 0.600  0.350 0.117 0.221 0.050  0.467             0.167

## 10201 0.333 0.600  0.000 0.233 0.612 0.400  0.233             0.333
## 10251 0.333 0.600  0.117 0.233 0.441 0.283  0.350             0.333
## 10301 0.333 0.600  0.233 0.233 0.270 0.167  0.466             0.333
## 10351 0.333 0.600  0.350 0.233 0.099 0.050  0.583             0.333

## 10401 0.333 0.600  0.000 0.350 0.544 0.400  0.350             0.500
## 10451 0.333 0.600  0.097 0.350 0.379 0.303  0.447             0.500
## 10501 0.333 0.600  0.194 0.350 0.214 0.206  0.544             0.500
## 10551 0.333 0.600  0.291 0.350 0.049 0.109  0.641             0.500

## 10601 0.333 0.600  0.000 0.467 0.449 0.400  0.467             0.667
## 10651 0.333 0.600  0.069 0.467 0.315 0.331  0.536             0.667
## 10701 0.333 0.600  0.137 0.467 0.184 0.263  0.604             0.667
## 10751 0.333 0.600  0.206 0.467 0.051 0.194  0.673             0.667

## 10801 0.333 0.600  0.000 0.583 0.326 0.400  0.583             0.833
## 10851 0.333 0.600  0.043 0.583 0.233 0.357  0.626             0.833
## 10901 0.333 0.600  0.085 0.583 0.142 0.315  0.668             0.833
## 10951 0.333 0.600  0.128 0.583 0.049 0.272  0.711             0.833

## 11001 0.333 0.600  0.000 0.700 0.177 0.400  0.700             1.000
## 11051 0.333 0.600  0.018 0.700 0.133 0.382  0.718             1.000
## 11101 0.333 0.600  0.035 0.700 0.093 0.365  0.735             1.000
## 11151 0.333 0.600  0.053 0.700 0.049 0.347  0.753             1.000

## 11201 0.467 0.200  0.000 0.000 0.533 0.800  0.000             0.000
## 11251 0.467 0.200  0.161 0.000 0.372 0.639  0.161             0.000
## 11301 0.467 0.200  0.322 0.000 0.211 0.478  0.322             0.000
## 11351 0.467 0.200  0.483 0.000 0.050 0.317  0.483             0.000

## 11401 0.467 0.200  0.000 0.116 0.520 0.800  0.116             0.167
## 11451 0.467 0.200  0.127 0.116 0.363 0.673  0.243             0.167
## 11501 0.467 0.200  0.255 0.116 0.206 0.545  0.371             0.167
## 11551 0.467 0.200  0.382 0.116 0.049 0.418  0.498             0.167

## 11601 0.467 0.200  0.000 0.232 0.480 0.800  0.232             0.333
## 11651 0.467 0.200  0.098 0.232 0.336 0.702  0.330             0.333
## 11701 0.467 0.200  0.196 0.232 0.193 0.604  0.428             0.333
## 11751 0.467 0.200  0.294 0.232 0.049 0.506  0.526             0.333

## 11801 0.467 0.200  0.000 0.347 0.413 0.800  0.347             0.500
## 11851 0.467 0.200  0.071 0.347 0.292 0.729  0.418             0.500
## 11901 0.467 0.200  0.143 0.347 0.170 0.657  0.490             0.500
## 11951 0.467 0.200  0.214 0.347 0.050 0.586  0.561             0.500

## 12001 0.467 0.200  0.000 0.463 0.319 0.800  0.463             0.667
## 12051 0.467 0.200  0.046 0.463 0.230 0.754  0.509             0.667
## 12101 0.467 0.200  0.093 0.463 0.139 0.707  0.556             0.667
## 12151 0.467 0.200  0.139 0.463 0.051 0.661  0.602             0.667

## 12201 0.467 0.200  0.000 0.579 0.198 0.800  0.579             0.833
## 12251 0.467 0.200  0.023 0.579 0.148 0.777  0.602             0.833
## 12301 0.467 0.200  0.046 0.579 0.099 0.754  0.625             0.833
## 12351 0.467 0.200  0.069 0.579 0.049 0.731  0.648             0.833

## 12401 0.467 0.200  0.000 0.695 0.050 0.800  0.695             1.000

## 12451 0.467 0.333  0.000 0.000 0.533 0.667  0.000             0.000
## 12501 0.467 0.333  0.161 0.000 0.372 0.506  0.161             0.000
## 12551 0.467 0.333  0.322 0.000 0.211 0.345  0.322             0.000
## 12601 0.467 0.333  0.483 0.000 0.050 0.184  0.483             0.000

## 12651 0.467 0.333  0.000 0.116 0.520 0.667  0.116             0.167
## 12701 0.467 0.333  0.127 0.116 0.363 0.540  0.243             0.167
## 12751 0.467 0.333  0.255 0.116 0.206 0.412  0.371             0.167
## 12801 0.467 0.333  0.382 0.116 0.049 0.285  0.498             0.167

## 12851 0.467 0.333  0.000 0.232 0.480 0.667  0.232             0.333
## 12901 0.467 0.333  0.098 0.232 0.336 0.569  0.330             0.333
## 12951 0.467 0.333  0.196 0.232 0.193 0.471  0.428             0.333
## 13001 0.467 0.333  0.294 0.232 0.049 0.373  0.526             0.333

## 13051 0.467 0.333  0.000 0.347 0.413 0.667  0.347             0.500
## 13101 0.467 0.333  0.071 0.347 0.292 0.596  0.418             0.500
## 13151 0.467 0.333  0.143 0.347 0.170 0.524  0.490             0.500
## 13201 0.467 0.333  0.214 0.347 0.050 0.453  0.561             0.500

## 13251 0.467 0.333  0.000 0.463 0.319 0.667  0.463             0.667
## 13301 0.467 0.333  0.046 0.463 0.230 0.621  0.509             0.667
## 13351 0.467 0.333  0.093 0.463 0.139 0.574  0.556             0.667
## 13401 0.467 0.333  0.139 0.463 0.051 0.528  0.602             0.667

## 13451 0.467 0.333  0.000 0.579 0.198 0.667  0.579             0.833
## 13501 0.467 0.333  0.023 0.579 0.148 0.644  0.602             0.833
## 13551 0.467 0.333  0.046 0.579 0.099 0.621  0.625             0.833
## 13601 0.467 0.333  0.069 0.579 0.049 0.598  0.648             0.833

## 13651 0.467 0.333  0.000 0.695 0.050 0.667  0.695             1.000

## 13701 0.467 0.467  0.000 0.000 0.533 0.533  0.000             0.000
## 13751 0.467 0.467  0.161 0.000 0.372 0.372  0.161             0.000
## 13801 0.467 0.467  0.322 0.000 0.211 0.211  0.322             0.000
## 13851 0.467 0.467  0.483 0.000 0.050 0.050  0.483             0.000

## 13901 0.467 0.467  0.000 0.116 0.520 0.533  0.116             0.167
## 13951 0.467 0.467  0.127 0.116 0.363 0.406  0.243             0.167
## 14001 0.467 0.467  0.255 0.116 0.206 0.278  0.371             0.167
## 14051 0.467 0.467  0.382 0.116 0.049 0.151  0.498             0.167

## 14101 0.467 0.467  0.000 0.232 0.480 0.533  0.232             0.333
## 14151 0.467 0.467  0.098 0.232 0.336 0.435  0.330             0.333
## 14201 0.467 0.467  0.196 0.232 0.193 0.337  0.428             0.333
## 14251 0.467 0.467  0.294 0.232 0.049 0.239  0.526             0.333

## 14301 0.467 0.467  0.000 0.347 0.413 0.533  0.347             0.500
## 14351 0.467 0.467  0.071 0.347 0.292 0.462  0.418             0.500
## 14401 0.467 0.467  0.143 0.347 0.170 0.390  0.490             0.500
## 14451 0.467 0.467  0.214 0.347 0.050 0.319  0.561             0.500

## 14501 0.467 0.467  0.000 0.463 0.319 0.533  0.463             0.667
## 14551 0.467 0.467  0.046 0.463 0.230 0.487  0.509             0.667
## 14601 0.467 0.467  0.093 0.463 0.139 0.440  0.556             0.667
## 14651 0.467 0.467  0.139 0.463 0.051 0.394  0.602             0.667

## 14701 0.467 0.467  0.000 0.579 0.198 0.533  0.579             0.833
## 14751 0.467 0.467  0.023 0.579 0.148 0.510  0.602             0.833
## 14801 0.467 0.467  0.046 0.579 0.099 0.487  0.625             0.833
## 14851 0.467 0.467  0.069 0.579 0.049 0.464  0.648             0.833

## 14901 0.467 0.467  0.000 0.695 0.050 0.533  0.695             1.000

## 14951 0.467 0.600  0.000 0.000 0.533 0.400  0.000             0.000
## 15001 0.467 0.600  0.117 0.000 0.416 0.283  0.117             0.000
## 15051 0.467 0.600  0.233 0.000 0.300 0.167  0.233             0.000
## 15101 0.467 0.600  0.350 0.000 0.183 0.050  0.350             0.000

## 15151 0.467 0.600  0.000 0.116 0.520 0.400  0.116             0.167
## 15201 0.467 0.600  0.117 0.116 0.376 0.283  0.233             0.167
## 15251 0.467 0.600  0.233 0.116 0.233 0.167  0.349             0.167
## 15301 0.467 0.600  0.350 0.116 0.089 0.050  0.466             0.167

## 15351 0.467 0.600  0.000 0.232 0.480 0.400  0.232             0.333
## 15401 0.467 0.600  0.098 0.232 0.336 0.302  0.330             0.333
## 15451 0.467 0.600  0.196 0.232 0.193 0.204  0.428             0.333
## 15501 0.467 0.600  0.294 0.232 0.049 0.106  0.526             0.333

## 15551 0.467 0.600  0.000 0.347 0.413 0.400  0.347             0.500
## 15601 0.467 0.600  0.071 0.347 0.292 0.329  0.418             0.500
## 15651 0.467 0.600  0.143 0.347 0.170 0.257  0.490             0.500
## 15701 0.467 0.600  0.214 0.347 0.050 0.186  0.561             0.500

## 15751 0.467 0.600  0.000 0.463 0.319 0.400  0.463             0.667
## 15801 0.467 0.600  0.046 0.463 0.230 0.354  0.509             0.667
## 15851 0.467 0.600  0.093 0.463 0.139 0.307  0.556             0.667
## 15901 0.467 0.600  0.139 0.463 0.051 0.261  0.602             0.667

## 15951 0.467 0.600  0.000 0.579 0.198 0.400  0.579             0.833
## 16001 0.467 0.600  0.023 0.579 0.148 0.377  0.602             0.833
## 16051 0.467 0.600  0.046 0.579 0.099 0.354  0.625             0.833
## 16101 0.467 0.600  0.069 0.579 0.049 0.331  0.648             0.833

## 16151 0.467 0.600  0.000 0.695 0.050 0.400  0.695             1.000

## 16201 0.600 0.200  0.000 0.000 0.400 0.800  0.000             0.000
## 16251 0.600 0.200  0.117 0.000 0.283 0.683  0.117             0.000
## 16301 0.600 0.200  0.233 0.000 0.167 0.567  0.233             0.000
## 16351 0.600 0.200  0.350 0.000 0.050 0.450  0.350             0.000

## 16401 0.600 0.200  0.000 0.099 0.390 0.800  0.099             0.167
## 16451 0.600 0.200  0.095 0.099 0.277 0.705  0.194             0.167
## 16501 0.600 0.200  0.189 0.099 0.164 0.611  0.288             0.167
## 16551 0.600 0.200  0.284 0.099 0.050 0.516  0.383             0.167

## 16601 0.600 0.200  0.000 0.197 0.361 0.800  0.197             0.333
## 16651 0.600 0.200  0.074 0.197 0.258 0.726  0.271             0.333
## 16701 0.600 0.200  0.149 0.197 0.153 0.651  0.346             0.333
## 16751 0.600 0.200  0.223 0.197 0.050 0.577  0.420             0.333

## 16801 0.600 0.200  0.000 0.296 0.312 0.800  0.296             0.500
## 16851 0.600 0.200  0.055 0.296 0.225 0.745  0.351             0.500
## 16901 0.600 0.200  0.110 0.296 0.137 0.690  0.406             0.500
## 16951 0.600 0.200  0.165 0.296 0.050 0.635  0.461             0.500

## 17001 0.600 0.200  0.000 0.395 0.244 0.800  0.395             0.667
## 17051 0.600 0.200  0.036 0.395 0.180 0.764  0.431             0.667
## 17101 0.600 0.200  0.073 0.395 0.114 0.727  0.468             0.667
## 17151 0.600 0.200  0.109 0.395 0.049 0.691  0.504             0.667

## 17201 0.600 0.200  0.000 0.493 0.157 0.800  0.493             0.833
## 17251 0.600 0.200  0.018 0.493 0.121 0.782  0.511             0.833
## 17301 0.600 0.200  0.036 0.493 0.085 0.764  0.529             0.833
## 17351 0.600 0.200  0.054 0.493 0.049 0.746  0.547             0.833

## 17401 0.600 0.200  0.000 0.592 0.050 0.800  0.592             1.000

## 17451 0.600 0.333  0.000 0.000 0.400 0.667  0.000             0.000
## 17501 0.600 0.333  0.117 0.000 0.283 0.550  0.117             0.000
## 17551 0.600 0.333  0.233 0.000 0.167 0.434  0.233             0.000
## 17601 0.600 0.333  0.350 0.000 0.050 0.317  0.350             0.000

## 17651 0.600 0.333  0.000 0.099 0.390 0.667  0.099             0.167
## 17701 0.600 0.333  0.095 0.099 0.277 0.572  0.194             0.167
## 17751 0.600 0.333  0.189 0.099 0.164 0.478  0.288             0.167
## 17801 0.600 0.333  0.284 0.099 0.050 0.383  0.383             0.167

## 17851 0.600 0.333  0.000 0.197 0.361 0.667  0.197             0.333
## 17901 0.600 0.333  0.074 0.197 0.258 0.593  0.271             0.333
## 17951 0.600 0.333  0.149 0.197 0.153 0.518  0.346             0.333
## 18001 0.600 0.333  0.223 0.197 0.050 0.444  0.420             0.333

## 18051 0.600 0.333  0.000 0.296 0.312 0.667  0.296             0.500
## 18101 0.600 0.333  0.055 0.296 0.225 0.612  0.351             0.500
## 18151 0.600 0.333  0.110 0.296 0.137 0.557  0.406             0.500
## 18201 0.600 0.333  0.165 0.296 0.050 0.502  0.461             0.500

## 18251 0.600 0.333  0.000 0.395 0.244 0.667  0.395             0.667
## 18301 0.600 0.333  0.036 0.395 0.180 0.631  0.431             0.667
## 18351 0.600 0.333  0.073 0.395 0.114 0.594  0.468             0.667
## 18401 0.600 0.333  0.109 0.395 0.049 0.558  0.504             0.667

## 18451 0.600 0.333  0.000 0.493 0.157 0.667  0.493             0.833
## 18501 0.600 0.333  0.018 0.493 0.121 0.649  0.511             0.833
## 18551 0.600 0.333  0.036 0.493 0.085 0.631  0.529             0.833
## 18601 0.600 0.333  0.054 0.493 0.049 0.613  0.547             0.833

## 18651 0.600 0.333  0.000 0.592 0.050 0.667  0.592             1.000

## 18701 0.600 0.467  0.000 0.000 0.400 0.533  0.000             0.000
## 18751 0.600 0.467  0.117 0.000 0.283 0.416  0.117             0.000
## 18801 0.600 0.467  0.233 0.000 0.167 0.300  0.233             0.000
## 18851 0.600 0.467  0.350 0.000 0.050 0.183  0.350             0.000

## 18901 0.600 0.467  0.000 0.099 0.390 0.533  0.099             0.167
## 18951 0.600 0.467  0.095 0.099 0.277 0.438  0.194             0.167
## 19001 0.600 0.467  0.189 0.099 0.164 0.344  0.288             0.167
## 19051 0.600 0.467  0.284 0.099 0.050 0.249  0.383             0.167

## 19101 0.600 0.467  0.000 0.197 0.361 0.533  0.197             0.333
## 19151 0.600 0.467  0.074 0.197 0.258 0.459  0.271             0.333
## 19201 0.600 0.467  0.149 0.197 0.153 0.384  0.346             0.333
## 19251 0.600 0.467  0.223 0.197 0.050 0.310  0.420             0.333

## 19301 0.600 0.467  0.000 0.296 0.312 0.533  0.296             0.500
## 19351 0.600 0.467  0.055 0.296 0.225 0.478  0.351             0.500
## 19401 0.600 0.467  0.110 0.296 0.137 0.423  0.406             0.500
## 19451 0.600 0.467  0.165 0.296 0.050 0.368  0.461             0.500

## 19501 0.600 0.467  0.000 0.395 0.244 0.533  0.395             0.667
## 19551 0.600 0.467  0.036 0.395 0.180 0.497  0.431             0.667
## 19601 0.600 0.467  0.073 0.395 0.114 0.460  0.468             0.667
## 19651 0.600 0.467  0.109 0.395 0.049 0.424  0.504             0.667

## 19701 0.600 0.467  0.000 0.493 0.157 0.533  0.493             0.833
## 19751 0.600 0.467  0.018 0.493 0.121 0.515  0.511             0.833
## 19801 0.600 0.467  0.036 0.493 0.085 0.497  0.529             0.833
## 19851 0.600 0.467  0.054 0.493 0.049 0.479  0.547             0.833

## 19901 0.600 0.467  0.000 0.592 0.050 0.533  0.592             1.000

##          ha    hb gamma2 omega    eA    eB rho.ab percent.max.omega
## 19951 0.600 0.600  0.000 0.000 0.400 0.400  0.000             0.000
## 20001 0.600 0.600  0.117 0.000 0.283 0.283  0.117             0.000
## 20051 0.600 0.600  0.233 0.000 0.167 0.167  0.233             0.000
## 20101 0.600 0.600  0.350 0.000 0.050 0.050  0.350             0.000

## 20151 0.600 0.600  0.000 0.099 0.390 0.400  0.099             0.167
## 20201 0.600 0.600  0.095 0.099 0.277 0.305  0.194             0.167
## 20251 0.600 0.600  0.189 0.099 0.164 0.211  0.288             0.167
## 20301 0.600 0.600  0.284 0.099 0.050 0.116  0.383             0.167

## 20351 0.600 0.600  0.000 0.197 0.361 0.400  0.197             0.333
## 20401 0.600 0.600  0.074 0.197 0.258 0.326  0.271             0.333
## 20451 0.600 0.600  0.149 0.197 0.153 0.251  0.346             0.333
## 20501 0.600 0.600  0.223 0.197 0.050 0.177  0.420             0.333

## 20551 0.600 0.600  0.000 0.296 0.312 0.400  0.296             0.500
## 20601 0.600 0.600  0.055 0.296 0.225 0.345  0.351             0.500
## 20651 0.600 0.600  0.110 0.296 0.137 0.290  0.406             0.500
## 20701 0.600 0.600  0.165 0.296 0.050 0.235  0.461             0.500

## 20751 0.600 0.600  0.000 0.395 0.244 0.400  0.395             0.667
## 20801 0.600 0.600  0.036 0.395 0.180 0.364  0.431             0.667
## 20851 0.600 0.600  0.073 0.395 0.114 0.327  0.468             0.667
## 20901 0.600 0.600  0.109 0.395 0.049 0.291  0.504             0.667

## 20951 0.600 0.600  0.000 0.493 0.157 0.400  0.493             0.833
## 21001 0.600 0.600  0.018 0.493 0.121 0.382  0.511             0.833
## 21051 0.600 0.600  0.036 0.493 0.085 0.364  0.529             0.833
## 21101 0.600 0.600  0.054 0.493 0.049 0.346  0.547             0.833

## 21151 0.600 0.600  0.000 0.592 0.050 0.400  0.592             1.000 


# (pseudo code)
#
# for N,h_a,cor(A,B),...fixed.
#
# The following contains R=200 repolicates for each choice of omega.

# Omega.scaled=omega/cor(A,B)

# Choose omega.scaled=0, 25,50 , 75, 100

# Omega= takes 4 or more distinct values including 0.

datIN=data.frame(omega.scaled, LEO.I)
LEO.threshold=2

# Hint: choose the threshold so that the POWER for omega=0 should is smaller than 0.01.

LEO.Icut=ifelse(LEO.I>LEO.threshold,1,0)

POWER=as.vector(tapply(LEO.Icut,omega.scaled,mean));

SE.POWER = as.vector(tapply(LEO.Icut,omega.scaled,stderr1))

par(mfrow=c(1,1))

barplot(POWER,names.arg=names(table(omega.scaled) ), ylab="POWER")

err.bp(as.vector(POWER), as.vector(1.96*SE.POWER), two.side=T, main=paste("No. Samples=", "LEO Threshold", LEO.threshold) )
# Produce the above barplots for different choices of the threshold


# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1)
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }

# ===================================================
# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)

err.bp<-function(daten,error,two.side=F){
 if(!is.numeric(daten)) {
      stop("All arguments must be numeric")}
 if(is.vector(daten)){ 
    xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
 }else{
    if (is.matrix(daten)){
      xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),dim=c(1,length(daten))))+0:(length(daten)-1)+.5
    }else{
      stop("First argument must either be a vector or a matrix") }
 }

 MW<-0.25*(max(xval)/length(xval))
 ERR1<-daten+error
 ERR2<-daten-error

 for(i in 1:length(daten)){
    segments(xval[i],daten[i],xval[i],ERR1[i])
    segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])

    if(two.side){
      segments(xval[i],daten[i],xval[i],ERR2[i])
      segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
    } 
 } 
}




###################
# percent causal omega on the x-axis plots...

bpo=c("rmsea.f1ab.AtoB","chisq.f1ab.AtoB","path.abs.z.f1ab.AtoB","rmsea.f1ab.BtoA","chisq.f1ab.BtoA","path.abs.z.f1ab.BtoA","rmsea.f1ab.AcollideB","path.abs.z.f1ab.AcollideB","chisq.f1ab.AcollideB","rmsea.f1ab.BcollideA","path.abs.z.f1ab.BcollideA","chisq.f1ab.BcollideA")

#bpo.conf=c("rmsea.f1ab.conf","chisq.f1ab.conf","zeo.ma.b.given.a","zeo.mb.a.given.b","cor.ma.b","pcor.ma.b.given.a","cor.mb.a","pcor.mb.a.given.b")

# take out the chisq
bpo.conf=c("rmsea.f1ab.conf","zeo.ma.b.given.a","zeo.mb.a.given.b","cor.ma.b","pcor.ma.b.given.a","cor.mb.a","pcor.mb.a.given.b")


 u8=unique(scores.new[w,1:8])
 num.reps = length(w)/nrow(u8); # automatically determine num.reps, rather than assume 50
 num.grp = nrow(u8)


w = 1:nrow(scores.new)
num.reps=10
splitter = (0:(length(w)-1)) %/% num.reps;
splitter = splitter+1 
scores.new$splitter = splitter



 bplist=split(scores.new$zeo.ma.b.given.a,scores.new$splitter[w])

 group.firsts= seq(1,length(w),num.reps)

 bps=rep(NA,length(bplist));
 gamma2=rep(NA,length(bplist));
 for (i in 1:length(bplist)) {
    bps[i]=boxplot.stats(bplist[[i]])$stats[3]
    gamma2[i]=scores.new$gamma2[w[group.firsts[i]]]
 }

  scale.factor=num.grp / max(gamma2)

 # combine the gamma2 that do coincide
 gamma2.splitter=rep(NA,nrow(scores.new))
 ug2=sort(unique(gamma2))
 gamma2.splitter=match(scores.new$gamma2,ug2)


 # a.man.split=split(scores.new$zeo.ma.b.given.a, scores.new$splitter[w]
 # names(a.man.split)=as.character(gamma2)
 man.split=split(scores.new$zeo.ma.b.given.a, gamma2.splitter)
 names(man.split)=as.character(ug2)

  postscript(file="bilav.under.one.marker.null.ps",horizontal=FALSE);  
#  boxplot(scores.new$zeo.ma.b.given.a ~ scores.new$gamma2,notch=TRUE,col="yellow",ylab="BILAV score = Z(|cor(MA,B)|)-Z(|cor(MA,B|A)|)",sub=paste(sep="","Gamma-squared: strength of A-B correlation when entirely due to confounding.\nSampling as ha and hb range from [0.2 to 0.6]."),main = "Bilayer verification: Negative BILAV score implies confounding.\nMonte Carlo study of Null: MA=>A, while B independent.")
#  boxplot(a.man.split,at=gamma2*scale.factor,notch=TRUE,col="yellow",ylab="BILAV score = Z(|cor(MA,B)|)-Z(|cor(MA,B|A)|)",sub=paste(sep="","Gamma-squared: strength of A-B correlation when entirely due to confounding.\nSampling as ha and hb range from [0.2 to 0.6]."),main = "Bilayer verification: Negative BILAV score implies confounding.\nMonte Carlo study of Null: MA=>A, while B independent.")
  boxplot(man.split,at=1+(ug2)*(length(ug2)-1)/max(ug2),notch=TRUE,col="yellow",ylab="BILAV score = Z(|cor(MA,B)|)-Z(|cor(MA,B|A)|)",sub=paste(sep="","Gamma-squared: strength of A-B correlation when entirely due to confounding.\nGrid sampling as ha and hb step through [0.2,0.33,0.47,0.6]."),main = "Bilayer verification: Negative BILAV score implies confounding.\nMonte Carlo study of Null: MA=>A, while B independent.")
  dev.off()


# sum over the correlations for each fixed percent.max.omega

if(exists("power.graph2") ) rm(power.graph2);
power.graph2=function(ha,hb,scores424,do.postscript=FALSE,do.labels=TRUE,...) {
 
 # pick out an ha and hb level with w.
 if (do.postscript) { 
     psfile=paste(sep="","power.plot.ha",round(ha,2),".hb",round(hb,2),"..rhoab.vs.leo.i.ps");
     postscript(file=psfile,horizontal=FALSE) 
 } else {
    # windows()
 }

 w=which(scores424$ha==ha & scores424$hb==hb)

 u8=unique(scores424[w,1:8])
 num.reps = length(w)/nrow(u8); # automatically determine num.reps, rather than assume 50
 num.grp = nrow(u8)

 #plot(scores424$rho.ab[w], scores424$rmsea.f2.M.M1M2.BtoA[w])


 leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])

 # split leo.i scores by the different parameter variations
 #bplist=split(leo.i,scores424$rho.ab[w])
 bplist=split(leo.i,scores424$splitter[w])

 group.firsts= seq(1,length(w),num.reps)

 bps=rep(NA,length(bplist));
 rhos=rep(NA,length(bplist));
 for (i in 1:length(bplist)) {
    bps[i]=boxplot.stats(bplist[[i]])$stats[3]
    rhos[i]=scores424$rho.ab[w[group.firsts[i]]]
 }

 scale.factor=num.grp / max(rhos) 

 # manually split out the levels, so we can be sure the
 # split is correct (not combining close or duplicated rho levels), and
 # then label the split groups correctly with their rho.ab.level, rhos.
 man.split=split(leo.i,scores424$splitter[w])
 names(man.split)=as.character(rhos)

 if (do.labels) {
  boxplot(leo.i ~ scores424$percent.max.omega[w],notch=F,ylab="LEO.I for MA,MB model",xlab=paste(sep="","percent max omega"),col="yellow")

#  boxplot(leo.i ~ scores424$rho.ab[w],notch=F,ylab="LEO.I for MA,MB model",xlab=paste(sep="","Horizontal axis: observed rho(A,B) = gamma^2 + omega. hA=",round(ha,2),", hB=",round(hb,2),". Dashed lines at LEO.I=2"),col="yellow",at=rhos*scale.factor,main="Monte Carlo study: detecting causal influence under partial\nconfounding: 0%, 17%, 33%, 50%, 67%, 83%, 100% causal",sub="Group labels show (omega/rho(A,B)) or percent true causal influence. Boxes for 50 replicates at N=200.")#,at=1:length(unique(scores424$percent.max.omega))-.2,boxwex=.25,ylim=c(-.15,.7))

  boxplot(man.split,notch=F,ylab="LEO.I for MA,MB model",xlab=paste(sep="","Horizontal axis: observed rho(A,B) = gamma^2 + omega. hA=",round(ha,2),", hB=",round(hb,2),". Dashed lines at LEO.I=2"),col="yellow",at=rhos*scale.factor,main="Monte Carlo study: detecting causal influence under partial\nconfounding: 0%, 17%, 33%, 50%, 67%, 83%, 100% causal",sub="Group labels show (omega/rho(A,B)) or percent true causal influence. Boxes for 50 replicates at N=200.")#,at=1:length(unique(scores424$percent.max.omega))-.2,boxwex=.25)

 } else {
  boxplot(man.split,notch=F,col="yellow",at=rhos*scale.factor,ylab="LEO.I for MA, MB model",sub=paste(sep="","Observed correlation(A,B)\nhA=",ha,",  hB=",hb),...)#,at=1:length(unique(scores424$percent.max.omega))-.2,boxwex=.25,ylim=c(-20,20))
 }


 # hit the center of the boxplots like this:
 points(x=rhos*scale.factor,y=bps,pch=19,col="red")

 # connect the centers that have the same percent.max.omega

 # split by groups that have the same percent.max.omega
 percent.omega.groups = split(w,scores424$percent.max.omega[w])
 ngroup = length(percent.omega.groups)
 
 # get the labels for the groups
 pog.nm=as.character(round(as.numeric(names(percent.omega.groups)),2))

 # omit last group: only 1 point
 for (j in 1:ngroup) {
 
 w.grp=percent.omega.groups[[j]]

 leo.i.grp = (scores424$chiprob.f2.M.M1M2hidden.con[w.grp]  - scores424$chiprob.f2.M.M1M2.BtoA[w.grp])

 rhos.grp=sort(unique(scores424$rho.ab[w.grp]))

 bplist.grp=split(leo.i.grp,scores424$rho.ab[w.grp])
 bps.grp=rep(NA,length(bplist.grp));

 bn = names(bplist.grp)

 for (i in 1:length(bplist.grp)) {
    bps.grp[i]=boxplot.stats(bplist.grp[[i]])$stats[3]
 }

 # take the max point in the first group as the reference location to label above
   bps.grp.top=max(bplist.grp[[1]])
   text(x=rhos.grp[1]*scale.factor,y=bps.grp.top+2,labels=pog.nm[j])

 # if (j != ngroup) { # don't plot lines for last group of singleton point.
    lines(x=rhos.grp*scale.factor,y=bps.grp)
 # }
} # end j

abline(a=2,b=0,lty=2)
#abline(v=.3*scale.factor,lty=3)  

  if (do.postscript) { 
     dev.off(); 
     system(paste(sep="","ps2pdf ",psfile),wait=T);
     system(paste(sep="","chmod +x ",paste(sep="",substr(psfile,1,(nchar(psfile,type="chars")-3)),".pdf")),wait=F);
  }
 invisible(0)
} # end power.graph2


for (ha in unique(scores424$ha)) {
  for (hb in unique(scores424$hb)) {
    power.graph2(ha, hb, scores424=scores424,do.postscript=T)
  }}

#####################################################
# LEO.I boxplots
#####################################################



# boxplot for percent max omega on the x axis
ha=hb=.6
     postscript(file="leo.i.vs.percent.true.signal.ps",horizontal=FALSE) 
 w=which(scores424$ha==ha & scores424$hb==hb)

 leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])

  boxplot(leo.i ~ scores424$percent.max.omega[w],notch=F,ylab="LEO.I for MA,MB model",col="yellow",boxwex=.25,ylim=c(-40,40))

 w2=which(scores424$ha==.2 & scores424$hb==.2)
 leo.i2 = (scores424$chiprob.f2.M.M1M2hidden.con[w2]  - scores424$chiprob.f2.M.M1M2.BtoA[w2])

 pmo=sort(unique(scores424$percent.max.omega[w2]))
 ngrp.pmo = length(pmo);
 boxplot(leo.i2 ~ scores424$percent.max.omega[w2],at=(1:ngrp.pmo)+.5,add=T,notch=F,ylab="LEO.I for MA,MB model",xlab=paste(sep="","Perent true signal (|omega|/|cor(A,B)|)"),col="violetred",main="Monte Carlo study: LEO.I detection of true signal as genetic signals\n(restricted heritabilities ha, hb) vary. Dashed line threshold: LEO.I=2.",boxwex=.25,show.names=FALSE)

 abline(a=2,b=0,lty=2)
  legend(5.1,-24.5, c("ha=0.6, hb=0.6", "ha=0.2, hb=0.2"), fill = c("yellow", "violetred"))
  dev.off()

# same thing, but with LEO.O
# same thing, but with AB.vs.NextBest

ha=hb=.6
 w=which(scores424$ha==ha & scores424$hb==hb)

 leo.o = (scores424$chiprob.f2.M.M1M2.AtoB[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])

# boxplot for percent max omega on the x axis
     postscript(file="leo.o.vs.percent.true.signal.ps",horizontal=FALSE) 
 w=which(scores424$ha==ha & scores424$hb==hb)
  boxplot(leo.o ~ scores424$percent.max.omega[w],notch=F,ylab="LEO.O for MA,MB model",col="yellow",boxwex=.25,ylim=c(-20,100))

 w2=which(scores424$ha==.2 & scores424$hb==.2)
 leo.o2 = (scores424$chiprob.f2.M.M1M2.AtoB[w2]  - scores424$chiprob.f2.M.M1M2.BtoA[w2])

 pmo=sort(unique(scores424$percent.max.omega[w2]))
 ngrp.pmo = length(pmo);
 boxplot(leo.o2 ~ scores424$percent.max.omega[w2],at=(1:ngrp.pmo)+.5,add=T,notch=F,ylab="LEO.O for MA,MB model",xlab=paste(sep="","Perent true signal (|omega|/|cor(A,B)|)"),col="violetred",main="Monte Carlo study: LEO.O detection of true signal as genetic signals\n(restricted heritabilities ha, hb) vary. Dashed line threshold: LEO.O=2.",boxwex=.25,show.names=FALSE)

 abline(h=2,lty=2)
  legend(1.4,67, c("ha=0.6, hb=0.6", "ha=0.2, hb=0.2"), fill = c("yellow", "violetred"))
  dev.off()



# same thing, but with AB.vs.NextBest

ha=hb=.6
 w=which(scores424$ha==ha & scores424$hb==hb)

next.best = pmin(scores424$chiprob.f2.M.M1M2.AtoB[w],
                scores424$chiprob.f2.M.M1M2unresolv[w],
                scores424$chiprob.f2.M.M1M2hidden.con[w])

 leo.nb = (next.best  - scores424$chiprob.f2.M.M1M2.BtoA[w])

# boxplot for percent max omega on the x axis
     postscript(file="leo.nb.vs.percent.true.signal.ps",horizontal=FALSE) 
 w=which(scores424$ha==ha & scores424$hb==hb)
  boxplot(leo.nb ~ scores424$percent.max.omega[w],notch=F,ylab="LEO.NB for MA,MB model",col="yellow",boxwex=.25,ylim=c(-20,40))

 w2=which(scores424$ha==.2 & scores424$hb==.2)
 leo.nb2 = (scores424$chiprob.f2.M.M1M2.AtoB[w2]  - scores424$chiprob.f2.M.M1M2.BtoA[w2])

 pmo=sort(unique(scores424$percent.max.omega[w2]))
 ngrp.pmo = length(pmo);
 boxplot(leo.nb2 ~ scores424$percent.max.omega[w2],at=(1:ngrp.pmo)+.5,add=T,notch=F,ylab="LEO.NB for MA,MB model",xlab=paste(sep="","Perent true signal (|omega|/|cor(A,B)|)"),col="violetred",main="Monte Carlo study: LEO.NB detection of true signal as genetic signals\n(restricted heritabilities ha, hb) vary. Dashed line threshold: LEO.NB=0.3",boxwex=.25,show.names=FALSE)

 abline(h=.3,lty=2)
  legend(1.2,29, c("ha=0.6, hb=0.6", "ha=0.2, hb=0.2"), fill = c("yellow", "violetred"))
  dev.off()







# do plots showing the advantage of BILAZ: start by giving the null distribution
scores.new$BILAZ=sqrt(197)*fisher1(abs(scores.new$cor.ma.b)) - sqrt(197)*fisher1(abs(scores.new$pcor.ma.b.given.a))

unique(scores.new[,1:8])

hahb6=scores.new[601:640,]
boxplot(hahb6$BILAZ ~ hahb6$gamma2 )

#k=0
#for (wg in seq(1,640,40)) {

# summarize the BILAZ studies....
for (ha in sort(unique(scores.new$ha))) {
#   windows(); 
 har = round(ha,2)
 postscript(file=paste(sep="","bilaz.under.null.ha",har,".ps"),horizontal=FALSE)
   par(mfrow=c(2,2))

   for (hb in sort(unique(scores.new$hb[scores.new$ha==ha]))) {
   w = which(scores.new$ha == ha & scores.new$hb == hb)

 hbr = round(hb,2)
 boxplot(scores.new$BILAZ[w] ~ scores.new$gamma2[w], ylab="BILAZ", xlab=paste(sep="","gamma^2 (confounded correlation);\n ha=",har,"; hb=",hbr),ylim=c(-8,2),col="yellow")
  abline(a=0,b=0,lty=2)
 }
 par(mfrow=c(1,1)); title(main="Monte Carlo study of BILAZ scores under Null model: no true A=>B signal",sub="10 simulations of N=200 at each level.")
 dev.off()

}

# Very start small, advnatage of using BILAZ score
ha=hb=.5

omega= c(0,  0,  0, .2, .2, .4)
gamma2=c(0, .2, .4,  0, .2,  0)

omega= c(.2, .4, .4)
gamma2=c(.4, .2, .4)


gamma=sqrt(gamma2)

eA=1-2*omega*gamma^2 - gamma^2 - omega^2 - ha

## > eA
## [1] 0.063

eB=1-gamma2 - hb;
## > eB
## [1] 0.25

rho.ab=gamma2+omega
## > rho.ab
## [1] 0.5

percent.max.omega=omega/rho.ab

q=data.frame(ha,hb,gamma2,omega,eA,eB,rho.ab,percent.max.omega)

## > q
##    ha  hb gamma2 omega    eA   eB rho.ab percent.max.omega
## 1 0.5 0.5   0.25  0.25 0.063 0.25    0.5               0.5
param=q

r=compare.single.vs.multiple.marker.simulations(param=param)

r$BILAZ.AtoB = r$zeo.ma.b.given.a
r$BILAZ.BtoA = r$zeo.mb.a.given.b

colnames(r)[1:8]=c("ha","hb","gamma2","omega","eA","eB","rho.ab","percent.max.omega")
save(r,file="r.small.bilaz.roc.curve.data.rdat")

param
#   ha  hb gamma2 omega   eA  eB rho.ab percent.max.omega
#1 0.5 0.5    0.4   0.0 0.10 0.1    0.4               0.0
#2 0.5 0.5    0.2   0.2 0.18 0.3    0.4               0.5
#3 0.5 0.5    0.0   0.4 0.34 0.5    0.4               1.0

w0g0=1:50
w0g2=51:100
w0g4=101:150
w2g0=151:200
w2g2=201:250
w4g0=251:300

w.alltrue = 101:150
w.allfalse = 1:50
w.halftrue = 51:100

# RMSEA power
windows()

postscript(file=paste(sep="","rmsea.3x3.gamma2.omega.isolated.ps"),horizontal=FALSE)
par(mfrow=c(3,3))

hist(r$rmsea.f2.M.M1M2.BtoA[w0g0],xlim=c(0,.6),main="",xlab=paste(sep="","omega=0  gamma^2=0"),col="yellow")
hist(r$rmsea.f2.M.M1M2.BtoA[w0g2],xlim=c(0,.6),main="",xlab=paste(sep="","omega=0  gamma^2=0.2"),col="yellow")
hist(r$rmsea.f2.M.M1M2.BtoA[w0g4],xlim=c(0,.6),main="",xlab=paste(sep="","omega=0  gamma^2=0.4"),col="orange")

hist(r$rmsea.f2.M.M1M2.BtoA[w2g0],xlim=c(0,.6),main="",xlab=paste(sep="","omega=0.2  gamma^2=0"),col="yellow")
hist(r$rmsea.f2.M.M1M2.BtoA[w2g2],xlim=c(0,.6),main="",xlab=paste(sep="","omega=0.2  gamma^2=0.2"),col="orange")
plot(c(),c())

hist(r$rmsea.f2.M.M1M2.BtoA[w4g0],xlim=c(0,.6),main="",xlab=paste(sep="","omega=0.4  gamma^2=0"),col="orange")


par(mfrow=c(1,1))
title("RMSEA in Monte Carlo study of Two marker MA=>A<=B<=MB model",sub="ha=hb=0.5; 50 replicates at N=200.")
dev.off()


# LEO.I
r$m2.leo.i = (r$chiprob.f2.M.M1M2hidden.con  - r$chiprob.f2.M.M1M2.BtoA)

windows()

postscript(file=paste(sep="","leo.i.3x3.gamma2.omega.isolated.ps"),horizontal=FALSE)
par(mfrow=c(3,3))


hist(r$m2.leo.i[w0g0],xlim=c(-40,20),main="",xlab=paste(sep="","omega=0  gamma^2=0"),col="yellow")
hist(r$m2.leo.i[w0g2],xlim=c(-40,20),main="",xlab=paste(sep="","omega=0  gamma^2=0.2"),col="yellow")
hist(r$m2.leo.i[w0g4],xlim=c(-40,20),main="",xlab=paste(sep="","omega=0  gamma^2=0.4"),col="orange")

hist(r$m2.leo.i[w2g0],xlim=c(-40,20),main="",xlab=paste(sep="","omega=0.2  gamma^2=0"),col="yellow")
hist(r$m2.leo.i[w2g2],xlim=c(-40,20),main="",xlab=paste(sep="","omega=0.2  gamma^2=0.2"),col="orange")
plot(c(),c())

hist(r$m2.leo.i[w4g0],xlim=c(-40,20),main="",xlab=paste(sep="","omega=0.4  gamma^2=0"),col="orange")

par(mfrow=c(1,1))
title("LEO.I in Monte Carlo study of Two marker MB=>B=>A<=MA model",sub="ha=hb=0.5; 50 replicates at N=200.")
dev.off()


windows()

postscript(file=paste(sep="","rmsea.power.full.half.no.omega.ps"),horizontal=FALSE)
par(mfrow=c(3,1))
#plot(c(),c())
hist(r$rmsea.f2.M.M1M2.BtoA[w.alltrue],xlim=c(0,.6),main="",sub="100% true signal",xlab="RMSEA",col="orange")
hist(r$rmsea.f2.M.M1M2.BtoA[w.halftrue],xlim=c(0,.6),main="",sub="50% true signal",xlab="RMSEA",breaks=9,col="orange")
hist(r$rmsea.f2.M.M1M2.BtoA[w.allfalse],xlim=c(0,.6),main="",sub="0% true signal",xlab="RMSEA",col="orange")

par(mfrow=c(1,1))
title("Monte Carlo power study of Two marker MB=>B=>A<=MA model\nha=hb=0.5; cor(A,B)=0.4; 50 replicates at N=200.")
dev.off()





# BILAZ
r$BILAZ.BtoA=sqrt(197)*fisher1(abs(r$cor.mb.a)) - sqrt(196)*fisher1(abs(r$pcor.mb.a.given.b))

windows()

postscript(file=paste(sep="","bilaz.3x3.gamma2.omega.isolated.ps"),horizontal=FALSE)
par(mfrow=c(3,3))


hist(r$BILAZ.BtoA[w0g0],xlim=c(-5,5),main="",xlab=paste(sep="","omega=0  gamma^2=0"),col="yellow")
hist(r$BILAZ.BtoA[w0g2],xlim=c(-5,5),main="",xlab=paste(sep="","omega=0  gamma^2=0.2"),col="yellow")
hist(r$BILAZ.BtoA[w0g4],xlim=c(-5,5),main="",xlab=paste(sep="","omega=0  gamma^2=0.4"),col="orange")

hist(r$BILAZ.BtoA[w2g0],xlim=c(-5,5),main="",xlab=paste(sep="","omega=0.2  gamma^2=0"),col="yellow")
hist(r$BILAZ.BtoA[w2g2],xlim=c(-5,5),main="",xlab=paste(sep="","omega=0.2  gamma^2=0.2"),col="orange")
plot(c(),c())

hist(r$BILAZ.BtoA[w4g0],xlim=c(-5,5),main="",xlab=paste(sep="","omega=0.4  gamma^2=0"),col="orange")

par(mfrow=c(1,1))
title("BILAZ in Monte Carlo study of Two marker MB=>B=>A<=MA model",sub="ha=hb=0.5; 50 replicates at N=200.")
dev.off()


# make barplots for LEO.I, LEO.O, LEO.NB that show the effect of the thresholding
# that depict the power of the method to detect causal true signal.

ha=hb=.6

 w=which(scores424$ha==ha & scores424$hb==hb)

 leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])

#dataIN=data.frame(ext.scores2$percent.max.omega, ext.scores2$rmsea.f2.M.M1M2.BtoA)
#colnames(dataIN)=c("percent.max.omega","RMSEA.M1M2.BtoA")
th=1

# Hint: choose the threshold so that the POWER for omega=0 should is smaller than 0.01.

# assign a 1 if the rmsea is below threshold: make an indicator whose mean is the frequency of cases below the Threshold.
LEO.Icut=ifelse(leo.i> th,1,0)

# group data by percent.max.omega, then get the frequency of rmsea being good in each group
# ; the mean or expected value of an indicator is just the probability or frequency of that
# indicator.

# tapply automatically sorts the percent.max.omega groups into ascending order.
POWER=as.vector(tapply(LEO.Icut,scores424$percent.max.omega[w],mean));
group.names=round(sort(unique(scores424$percent.max.omega[w])),2)
names(POWER)=as.character(group.names)

SE.POWER= as.vector(tapply(LEO.Icut,scores424$percent.max.omega[w],stderr1))

ci.l = POWER-1.96*SE.POWER
ci.u = POWER+1.96*SE.POWER


par(mfrow=c(1,1))

#barplot(POWER,names.arg=names(table(dataIN$percent.max.omega) ), ylab="POWER",main=paste("No. Samples=",50, "RMSEA Threshold", rmsea.threshold))
#err.bp(as.vector(POWER), as.vector(1.96*SE.POWER), two.side=T)

# Produce the above barplots for different choices of the threshold

   windows()
     postscript(file="leo.i.BOXPLOT.vs.percent.true.signal.ha.hb.0.6.ps",horizontal=FALSE)

#ci.l=rbind(ci.l,ci.l)
#ci.u=rbind(ci.u,ci.u)
#POWER2=rbind(POWER,POWER)

   mp <- barplot2(POWER,
             beside = TRUE,
        col = "yellow",
#        col = c("yellow", "violet"),
#                "lightcyan", "orange"),
#        legend = 
#                  ylim = c(0, 1),
        main = "Power of LEO.I in Two Marker MA=>A<=B<=MB model\nha=hb=0.6",
        font.main = 4,
                   ylab="Power",xlab=paste(sep="","Proportion causal at LEO.I significance threshold = ",signif(th,2)),
#        sub = "Faked 95 percent error bars",
#        col.sub = "yellow",
        cex.lab = 1.5,
        cex.names = 1.5,
        cex.axis = 1.5,
        plot.ci = TRUE,
               ci.l = ci.l,
               ci.u = ci.u,
        plot.grid = TRUE)

dev.off()


# repeat the same, but with different heritability

ha=hb=.2

 w=which(scores424$ha==ha & scores424$hb==hb)

 leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])

#dataIN=data.frame(ext.scores2$percent.max.omega, ext.scores2$rmsea.f2.M.M1M2.BtoA)
#colnames(dataIN)=c("percent.max.omega","RMSEA.M1M2.BtoA")
th=1

# Hint: choose the threshold so that the POWER for omega=0 should is smaller than 0.01.

# assign a 1 if the rmsea is below threshold: make an indicator whose mean is the frequency of cases below the Threshold.
LEO.Icut=ifelse(leo.i> th,1,0)

# group data by percent.max.omega, then get the frequency of rmsea being good in each group
# ; the mean or expected value of an indicator is just the probability or frequency of that
# indicator.

# tapply automatically sorts the percent.max.omega groups into ascending order.
POWER=as.vector(tapply(LEO.Icut,scores424$percent.max.omega[w],mean));
group.names=round(sort(unique(scores424$percent.max.omega[w])),2)
names(POWER)=as.character(group.names)

SE.POWER= as.vector(tapply(LEO.Icut,scores424$percent.max.omega[w],stderr1))

ci.l = POWER-1.96*SE.POWER
ci.u = POWER+1.96*SE.POWER


par(mfrow=c(1,1))

#barplot(POWER,names.arg=names(table(dataIN$percent.max.omega) ), ylab="POWER",main=paste("No. Samples=",50, "RMSEA Threshold", rmsea.threshold))
#err.bp(as.vector(POWER), as.vector(1.96*SE.POWER), two.side=T)

# Produce the above barplots for different choices of the threshold

   windows()
     postscript(file="leo.i.BOXPLOT.vs.percent.true.signal.ha.hb.0.2.ps",horizontal=FALSE)

#ci.l=rbind(ci.l,ci.l)
#ci.u=rbind(ci.u,ci.u)
#POWER2=rbind(POWER,POWER)

   mp <- barplot2(POWER,
             beside = TRUE,
        col = "violetred",
#        col = c("yellow", "violet"),
#                "lightcyan", "orange"),
#        legend = 
                  ylim = c(0, 1),
        main = "Power of LEO.I in Two Marker MA=>A<=B<=MB model\nha=hb=0.2",
        font.main = 4,
                   ylab="Power",xlab=paste(sep="","Proportion causal at LEO.I significance threshold = ",signif(th,2)),
#        sub = "Faked 95 percent error bars",
#        col.sub = "yellow",
        cex.lab = 1.5,
        cex.names = 1.5,
        cex.axis = 1.5,
        plot.ci = TRUE,
               ci.l = ci.l,
               ci.u = ci.u,
        plot.grid = TRUE)

dev.off()


##### funtion form

hahb=ha=hb=.6
 w=which(scores424$ha==ha & scores424$hb==hb)
score = (scores424$chiprob.f2.M.M1M2.AtoB[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])
score.name = "LEO.O"
th=1

power.plot(score.name="LEO.O",score=score,th=1,hahb=.6,scores424=scores424,w=w)



if(exists("power.plot") ) rm(power.plot);
power.plot=function(do.postscript=FALSE,hahb,score, score.name, th, my.color="yellow",scores424,w) {
 

#dataIN=data.frame(ext.scores2$percent.max.omega, ext.scores2$rmsea.f2.M.M1M2.BtoA)
#colnames(dataIN)=c("percent.max.omega","RMSEA.M1M2.BtoA")

# Hint: choose the threshold so that the POWER for omega=0 should is smaller than 0.01.

# assign a 1 if the rmsea is below threshold: make an indicator whose mean is the frequency of cases below the Threshold.
LEO.Icut=ifelse(score> th,1,0)

# group data by percent.max.omega, then get the frequency of rmsea being good in each group
# ; the mean or expected value of an indicator is just the probability or frequency of that
# indicator.

# tapply automatically sorts the percent.max.omega groups into ascending order.
POWER=as.vector(tapply(LEO.Icut,scores424$percent.max.omega[w],mean));
group.names=round(sort(unique(scores424$percent.max.omega[w])),2)
names(POWER)=as.character(group.names)

SE.POWER= as.vector(tapply(LEO.Icut,scores424$percent.max.omega[w],stderr1))

ci.l = POWER-1.96*SE.POWER
ci.u = POWER+1.96*SE.POWER


par(mfrow=c(1,1))

#barplot(POWER,names.arg=names(table(dataIN$percent.max.omega) ), ylab="POWER",main=paste("No. Samples=",50, "RMSEA Threshold", rmsea.threshold))
#err.bp(as.vector(POWER), as.vector(1.96*SE.POWER), two.side=T)

# Produce the above barplots for different choices of the threshold

  if (!do.postscript) windows()
  if (do.postscript)   postscript(file=paste(sep="",score.name,".BARPLOT.vs.percent.true.signal.ha.hb.",hahb,".ps"),horizontal=FALSE)

#ci.l=rbind(ci.l,ci.l)
#ci.u=rbind(ci.u,ci.u)
#POWER2=rbind(POWER,POWER)

   mp <- barplot2(POWER,
             beside = TRUE,
        col = my.color,
#        col = c("yellow", "violet"),
#                "lightcyan", "orange"),
#        legend = 
                  ylim = c(0, 1),
        main = paste(sep="","Power of ",score.name," in Two Marker MA=>A<=B<=MB model\nha=hb=",hahb),
        font.main = 4,
                   ylab="Power",xlab=paste(sep="","Proportion causal at ",score.name," threshold = ",signif(th,2)),
#        sub = "Faked 95 percent error bars",
#        col.sub = "yellow",
        cex.lab = 1.5,
        cex.names = 1.5,
        cex.axis = 1.5,
        plot.ci = TRUE,
               ci.l = ci.l,
               ci.u = ci.u,
        plot.grid = TRUE)

   for (k in 1:length(POWER)) {
      text(x=mp[k],y=(ci.u[k])+.02,labels=paste(sep="",as.character(signif(POWER[k],2))," (",as.character(signif(SE.POWER[k],2)),")"))
    }
  box()
  if (do.postscript) dev.off()

  invisible(data.frame(POWER,SE.POWER,ci.l,ci.u))
}



hahb=ha=hb=.6
 w=which(scores424$ha==ha & scores424$hb==hb)
score = (scores424$chiprob.f2.M.M1M2.AtoB[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])
score.name = "LEO.O"
th=5

power.plot(score.name="LEO.O",score=score,th=th,hahb=.6,scores424=scores424,w=w,do.postscript=T)

hahb=ha=hb=.2
 w=which(scores424$ha==ha & scores424$hb==hb)
score = (scores424$chiprob.f2.M.M1M2.AtoB[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])
score.name = "LEO.O"
th=5

z=power.plot(score.name="LEO.O",score=score,th=th,hahb=hahb,scores424=scores424,w=w,my.color="violetred",do.postscript=T)



hahb=ha=hb=.6
 w=which(scores424$ha==ha & scores424$hb==hb)


next.best = pmin(scores424$chiprob.f2.M.M1M2.AtoB[w],
                scores424$chiprob.f2.M.M1M2unresolv[w],
                scores424$chiprob.f2.M.M1M2hidden.con[w])

score= leo.nb = (next.best  - scores424$chiprob.f2.M.M1M2.BtoA[w])

score.name = "LEO.NB"
th=0.3

z=power.plot(score.name="LEO.NB",score=score,th=th,hahb=.6,scores424=scores424,w=w,do.postscript=TRUE)

####################

hahb=ha=hb=.2
w=which(scores424$ha==ha & scores424$hb==hb)
next.best = pmin(scores424$chiprob.f2.M.M1M2.AtoB[w],
                scores424$chiprob.f2.M.M1M2unresolv[w],
                scores424$chiprob.f2.M.M1M2hidden.con[w])

score= leo.nb = (next.best  - scores424$chiprob.f2.M.M1M2.BtoA[w])

score.name = "LEO.NB"
th=0.3

z=power.plot(score.name=score.name,score=score,th=th,hahb=hahb,scores424=scores424,w=w,my.color="violetred",do.postscript=TRUE)

################
score.name="LEO.I"
hahb=ha=hb=.6
w=which(scores424$ha==ha & scores424$hb==hb)
score=leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])
th=1
z=power.plot(score.name=score.name,score=score,th=th,hahb=hahb,scores424=scores424,w=w,do.postscript=T)

hahb=ha=hb=.2
w=which(scores424$ha==ha & scores424$hb==hb)
score = leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])
th=1
z=power.plot(score.name=score.name,score=score,th=th,hahb=hahb,scores424=scores424,w=w,my.color="violetred",do.postscript=T)

##################
#compare.single.vs.multiple.markers LEO.I, LEO.O, LEO.NB sensitivities

hahb=ha=hb=.6
w=which(scores424$ha==ha & scores424$hb==hb)


score=mm.leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w]  - scores424$chiprob.f2.M.M1M2.BtoA[w])

score=single.leo.i = scores424$chiprob.f1ba.conf[w] - scores424$chiprob.f1ba.BtoA[w]
single.leo.o = scores424$chiprob.f1ba.AtoB[w] - scores424$chiprob.f1ba.BtoA[w]
single.leo.nb = pmin(scores424$chiprob.f1ba.AtoB[w], scores424$chiprob.f1ba.conf[w], scores424$chiprob.f1ba.AcollideB[w], scores424$chiprob.f1ba.BcollideA[w]) - scores424$chiprob.f1ba.BtoA[w]

score.name="Single.Marker.LEO.I"
score=single.leo.i
th=8
z=power.plot(score.name=score.name,score=score,th=th,hahb=hahb,scores424=scores424,w=w,do.postscript=F)


hahb=ha=hb=.2
w=which(scores424$ha==ha & scores424$hb==hb)



score=single.leo.i = scores424$chiprob.f1ba.conf[w] - scores424$chiprob.f1ba.BtoA[w]
single.leo.o = scores424$chiprob.f1ba.AtoB[w] - scores424$chiprob.f1ba.BtoA[w]
single.leo.nb = pmin(scores424$chiprob.f1ba.AtoB[w], scores424$chiprob.f1ba.conf[w], scores424$chiprob.f1ba.AcollideB[w], scores424$chiprob.f1ba.BcollideA[w]) - scores424$chiprob.f1ba.BtoA[w]

th=24

z=power.plot(score.name=score.name,score=score,th=th,hahb=hahb,scores424=scores424,w=w,my.color="violetred",do.postscript=F)

#######################
# comparing single vs multiple markers
#######################
load("scores424.rdat") # scores424
library(gplots)

# this function computes the standard error
if (exists("stderr1")) rm(stderr1)
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }


######################

# function to compare side-by-side the single versus double marker in barplots

if(exists("side.by.side.power.plot") ) rm(side.by.side.power.plot);
side.by.side.power.plot=function(do.postscript=FALSE,hahb1, hahb2 ,score1, score2, score1.name, score2.name, th1, th2, my.color1="yellow",my.color2="violetred",scores424,w1,w2,main=paste(sep="","Power comparison: two marker versus one marker scores\nMonte Carlo study using MA=>A<=B<=MB model; ha=hb=",hahb),legend.cex=1,...) {
 
# assign a 1 if the rmsea is below threshold: make an indicator whose mean is the frequency of cases below the Threshold.
LEO.Icut1=ifelse(score1> th1,1,0)
LEO.Icut2=ifelse(score2> th2,1,0)

# group data by percent.max.omega, then get the frequency of rmsea being good in each group
# ; the mean or expected value of an indicator is just the probability or frequency of that
# indicator.

# tapply automatically sorts the percent.max.omega groups into ascending order.


# for 1

POWER1=as.vector(tapply(LEO.Icut1,scores424$percent.max.omega[w1],mean));
group.names1=round(sort(unique(scores424$percent.max.omega[w1])),2)
names(POWER1)=as.character(group.names1)

SE.POWER1 = as.vector(tapply(LEO.Icut1,scores424$percent.max.omega[w1],stderr1))

#ci.l1 = POWER1-1.96*SE.POWER1
#ci.u1 = POWER1+1.96*SE.POWER1

ci.l1 = POWER1-SE.POWER1
ci.u1 = POWER1+SE.POWER1

# for 2

POWER2=as.vector(tapply(LEO.Icut2,scores424$percent.max.omega[w2],mean));
group.names2=round(sort(unique(scores424$percent.max.omega[w2])),2)
names(POWER2)=as.character(group.names2)

SE.POWER2 = as.vector(tapply(LEO.Icut2,scores424$percent.max.omega[w2],stderr1))

ci.l2 = POWER2-SE.POWER2
ci.u2 = POWER2+SE.POWER2

par(mfrow=c(1,1))

  

  if (!do.postscript) windows()
  if (do.postscript)   postscript(file=paste(sep="",gsub(" ",".",score1.name),".vs.",gsub(" ",".",score2.name),".BARPLOT.vs.percent.true.signal.side.by.side.",hahb1,".vs.",hahb2,".ps"),horizontal=FALSE)

   mp1 <- barplot2(rbind(POWER1,POWER2),
             beside = TRUE,
        col = c(my.color1,my.color2),
                  ylim = c(0, 1),
        main = main,
        font.main = 4,
        ylab="statistical power",
        xlab=paste(sep="","Percent causal signal (|omega|/|cor(A,B)|)"),
        cex.lab = 1.8,
        cex.names = 1.5,
        cex.axis = 1.5,
        plot.ci = TRUE,
               ci.l = rbind(ci.l1,ci.l2),
               ci.u = rbind(ci.u1,ci.u2),
        plot.grid = TRUE,...)

   for (k in 1:length(POWER1)) {
      text(x=mp1[1,k],y=(ci.u1[k])+.02,cex=.7,labels=paste(sep="",as.character(signif(POWER1[k],2))," (",as.character(signif(SE.POWER1[k],2)),")"))
    }

##    mp2 <- barplot2(POWER2, add=TRUE,
##              beside = TRUE,
##         col = my.color2,
##                   ylim = c(0, 1),
## #        main = paste(sep="","Power of ",score.name," in Two Marker MA=>A<=B<=MB model\nha=hb=",hahb),
##         font.main = 4,
##                    ylab="Power",
##   xlab=paste(sep="","Proportion causal at ",score.name," threshold = ",signif(th2,2)),
##         cex.lab = 1.5,
##         cex.names = 1.5,
##         cex.axis = 1.5,
##         plot.ci = TRUE,
##                ci.l = ci.l2,
##                ci.u = ci.u2,
##         plot.grid = TRUE)

   for (k in 1:length(POWER2)) {
      text(x=mp1[2,k],y=(ci.u2[k])+.02+.02*(k==1),cex=.7,labels=paste(sep="",as.character(signif(POWER2[k],2))," (",as.character(signif(SE.POWER2[k],2)),")"))
    }



  box()

  legend(1.32,.946, c(paste(sep="",score1.name,", threshold ",th1),paste(sep="",score2.name,", threshold ",th2)), fill = c(my.color1, my.color2),cex=legend.cex)



  if (do.postscript) dev.off()

  invisible(data.frame(POWER1,SE.POWER1,ci.l1,ci.u1,POWER2,SE.POWER2,ci.l2,ci.u2))
}


##################


hahb1=ha=hb=.6
hahb2=ha=hb=.6
w1=which(scores424$ha==ha & scores424$hb==hb)
w2=which(scores424$ha==ha & scores424$hb==hb)

score1=mm.leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w1]  - scores424$chiprob.f2.M.M1M2.BtoA[w1])
score2=single.leo.i = scores424$chiprob.f1ba.conf[w2] - scores424$chiprob.f1ba.BtoA[w2]

score1.name="Two Marker LEO.I"
score2.name="One Marker LEO.I"

th1=.7
th2=8

zz=side.by.side.power.plot(do.postscript=T,hahb1, hahb2 ,score1, score2, score1.name, score2.name, th1, th2, my.color1="violetred",my.color2="palegreen",scores424,w1,w2,main=paste(sep="","Power comparison: two marker versus one marker scores\nMonte Carlo study using MA=>A<=B<=MB model; ha=hb=",hahb1))

####################


hahb1=ha=hb=.2
hahb2=ha=hb=.2
w1=which(scores424$ha==ha & scores424$hb==hb)
w2=which(scores424$ha==ha & scores424$hb==hb)

score1=mm.leo.i = (scores424$chiprob.f2.M.M1M2hidden.con[w1]  - scores424$chiprob.f2.M.M1M2.BtoA[w1])
score2=single.leo.i = scores424$chiprob.f1ba.conf[w2] - scores424$chiprob.f1ba.BtoA[w2]

score1.name="Two Marker LEO.I"
score2.name="One Marker LEO.I"

th1=.3
th2=24

zz=side.by.side.power.plot(do.postscript=T,hahb1, hahb2 ,score1, score2, score1.name, score2.name, th1, th2, my.color1="violetred",my.color2="palegreen",scores424,w1,w2,main=paste(sep="","Power comparison: two marker versus one marker scores\nMonte Carlo study using MA=>A<=B<=MB model; ha=hb=",hahb1))




##############
# more single vs multiple comparisons

hahb1=ha=hb=.6
hahb2=ha=hb=.6
w1=which(scores424$ha==ha & scores424$hb==hb)
w2=which(scores424$ha==ha & scores424$hb==hb)

score2 = single.leo.nb = pmin(scores424$chiprob.f1ba.AtoB[w2], scores424$chiprob.f1ba.conf[w2], scores424$chiprob.f1ba.AcollideB[w2], scores424$chiprob.f1ba.BcollideA[w2]) - scores424$chiprob.f1ba.BtoA[w2]

next.best = pmin(scores424$chiprob.f2.M.M1M2.AtoB[w1],
                scores424$chiprob.f2.M.M1M2unresolv[w1],
                scores424$chiprob.f2.M.M1M2hidden.con[w1])

score1 = mm.leo.nb = (next.best  - scores424$chiprob.f2.M.M1M2.BtoA[w1])

score1.name="Two Marker LEO.NB"
score2.name="One Marker LEO.NB"

th1=.3
th2=1

zz=side.by.side.power.plot(do.postscript=F,hahb1, hahb2 ,score1, score2, score1.name, score2.name, th1, th2, my.color1="violetred",my.color2="palegreen",scores424,w1,w2,main=paste(sep="","Power comparison: two marker versus one marker scores\nMonte Carlo study using MA=>A<=B<=MB model; ha=hb=",hahb1),cex.main=1.3)

#zz=side.by.side.power.plot(do.postscript=FALSE,hahb1, hahb2 ,score1, score2, score1.name, score2.name, th1, th2, my.color1="blue",my.color2="orange",scores424,w1,w2,main=paste(sep="","Power comparison: two marker versus one marker scores\nMonte Carlo study using MA=>A<=B<=MB model; ha=hb=",hahb1))
#  legend(1.32,.946, c(paste(score1.name,", threshold",th1),paste(score2.name,", threshold",th2)), fill = c("blue", "orange"))


# .2 version
hahb1=ha=hb=.2
hahb2=ha=hb=.2
w1=which(scores424$ha==ha & scores424$hb==hb)
w2=which(scores424$ha==ha & scores424$hb==hb)

score2 = single.leo.nb = pmin(scores424$chiprob.f1ba.AtoB[w2], scores424$chiprob.f1ba.conf[w2], scores424$chiprob.f1ba.AcollideB[w2], scores424$chiprob.f1ba.BcollideA[w2]) - scores424$chiprob.f1ba.BtoA[w2]

next.best = pmin(scores424$chiprob.f2.M.M1M2.AtoB[w1],
                scores424$chiprob.f2.M.M1M2unresolv[w1],
                scores424$chiprob.f2.M.M1M2hidden.con[w1])

score1 = mm.leo.nb = (next.best  - scores424$chiprob.f2.M.M1M2.BtoA[w1])

score1.name="Two Marker LEO.NB"
score2.name="One Marker LEO.NB"

th1=.3
th2=.7

zz=side.by.side.power.plot(do.postscript=F,hahb1, hahb2 ,score1, score2, score1.name, score2.name, th1, th2, my.color1="violetred",my.color2="palegreen",scores424,w1,w2,main=paste(sep="","Power comparison: two marker versus one marker scores\nMonte Carlo study using MA=>A<=B<=MB model; ha=hb=",hahb1),cex.main=1.3)





################

hahb1=ha=hb=.6
hahb2=ha=hb=.6
w1=which(scores424$ha==ha & scores424$hb==hb)
w2=which(scores424$ha==ha & scores424$hb==hb)

single.leo.o.using.mb = scores424$chiprob.f1ba.BcollideA[w2] - scores424$chiprob.f1ba.BtoA[w2] # correct, should be positive
single.leo.o.using.ma = scores424$chiprob.f1ab.AcollideB[w2] - scores424$chiprob.f1ab.AtoB[w2] # opposite, should be negative

score2 = single.leo.o.using.mb

score1 = mm.leo.o = (scores424$chiprob.f2.M.M1M2.AtoB[w1]  - scores424$chiprob.f2.M.M1M2.BtoA[w1])

score1.name="Two Marker LEO.O"
score2.name="One Marker LEO.O"

th1=6
th2=.7

zz=side.by.side.power.plot(do.postscript=T,hahb1, hahb2 ,score1, score2, score1.name, score2.name, th1, th2, my.color1="violetred",my.color2="palegreen",scores424,w1,w2,main=paste(sep="","Power comparison: two marker versus one marker scores\nMonte Carlo study using MA=>A<=B<=MB model; ha=hb=",hahb1),legend.cex=.8,cex.main=1.3)

##################################

hahb1=ha=hb=.2
hahb2=ha=hb=.2
w1=which(scores424$ha==ha & scores424$hb==hb)
w2=which(scores424$ha==ha & scores424$hb==hb)

single.leo.o.using.mb = scores424$chiprob.f1ba.BcollideA[w2] - scores424$chiprob.f1ba.BtoA[w2] # correct, should be positive
single.leo.o.using.ma = scores424$chiprob.f1ab.AcollideB[w2] - scores424$chiprob.f1ab.AtoB[w2] # opposite, should be negative

score2 = single.leo.o.using.mb

score1 = mm.leo.o = (scores424$chiprob.f2.M.M1M2.AtoB[w1]  - scores424$chiprob.f2.M.M1M2.BtoA[w1])

score1.name="Two Marker LEO.O"
score2.name="One Marker LEO.O"

th1=4
th2=.7

zz=side.by.side.power.plot(do.postscript=T,hahb1, hahb2 ,score1, score2, score1.name, score2.name, th1, th2, my.color1="violetred",my.color2="palegreen",scores424,w1,w2,main=paste(sep="","Power comparison: two marker versus one marker scores\nMonte Carlo study using MA=>A<=B<=MB model; ha=hb=",hahb1),legend.cex=.8,cex.main=1.3)

###########################

postscript(file="leo.o.predicted.spurious.from.zeo.ps",horizontal=FALSE)
 plot(leo.o.for,zeo.for,xlab="LEO.O",ylab="ZEO",main="Predicting LEO.O scores from ZEO scores: under spurious signal",cex.main=1.3,cex.lab=1.5)
 abline(lm(zeo.for~leo.o.for))
dev.off()


postscript(file="leo.nb.predicted.spurious.from.zeo.ps",horizontal=FALSE)
plot(leo.nb.for,zeo.for,xlab="LEO.NB",ylab="ZEO",main="Predicting LEO.NB scores from ZEO scores: under spurious signal",cex.main=1.3,cex.lab=1.5)
 abline(lm(zeo.for~leo.nb.for))
dev.off()


postscript(file="leo.i.predicted.spurious.from.zeo.ps",horizontal=FALSE)
 plot(leo.i.for,zeo.for,xlab="LEO.I",ylab="ZEO",main="Predicting LEO.I scores from ZEO scores: under spurious signal",cex.main=1.3,cex.lab=1.5)
 abline(lm(zeo.for~leo.i.for))
dev.off()




save(  leo.nb.for,  zeo.for,  leo.i.for, leo.o.for, file="leo.vs.zeo.data.under.null.0.4cor.rdat")


} # end recall.barplots

# =================================================
#
# end recall.barplots()
#
# =================================================




# helper function for make.rqtl.file
translate.012.AHB = function(x) {
  xn = x
  xn[] = "-"
  xn[x==0]="A"
  xn[x==1]="H"
  xn[x==2]="B"
  xn
}

# convert from datCombined to a csv file subtable for qtl::read.cross(); Karl Broman's R/qtl package.
make.rqtl.file=function(x,snpcols,traitcols,file.name="my.read.cross.file.csv", has.cm.not.bp=FALSE,has.bp.and.morgans=FALSE,convert.chr20.to.X=FALSE) {

  library(qtl)
  
  ## from library(qtl)
  ## ?read.cross
  
  ## The input file is a comma-delimited text file. A different field separator
  ## may be specified via the argument sep, which will be passed to the function
  ## read.table). For example, in Europe, it is common to use a comma in place of
  ## the decimal point in numbers and so a semi-colon in place of a comma as the
  ## field separator; such data may be read by using sep=";" and dec=",". 
  ##
  ## The first line should contain the phenotype names followed by the marker names.
  ## At least one phenotype must be included; for example, include a numerical index for each individual. 
  ##
  ## The second line should contain blanks in the phenotype columns, followed by
  ## chromosome identifiers for each marker in all other columns. If a chromosome has
  ## the identifier X or x, it is assumed to be the X chromosome; otherwise, it is assumed to be an autosome. 
  ##
  ## An optional third line should contain blanks in the phenotype columns, followed by marker positions, in cM. 
  ##
  ## Marker order is taken from the cM positions, if provided; otherwise, it is taken from the column order. 
  ##
  ## Subsequent lines should give the data, with one line for each individual,
  ## and with phenotypes followed by genotypes. If possible, phenotypes are
  ## made numeric; otherwise they are converted to factors. 
  ##
  ## The cross is determined to be a backcross if only the first two elements of the
  ## genotypes string are found; otherwise, it is assumed to be an intercross. 
  
  # make a character data frame

  ##   T264 D10M44  D1M3  D1M75 D1M215
  ##          1 1 1 1  <- chromosome
  ##          0 0.99675 24.84773  40.41361 <- centomorgans
  ## 118.317  B B B H
  ## 264  - B B B
  ## 194.917  - H H H
  ## 264  B B H H
  ## 145.417  H H H H

  nc=ncol(x)
  nr=nrow(x)
  cn=colnames(x)

  nt = length(traitcols)
  snp.start = nt + 1

  #split out the chromosome and position information
  snp.chr.pos = strsplit(cn[snpcols],split=".",fixed=TRUE)

  if (length(snp.chr.pos[[1]]) < 3) stop("snpcols column names must contain snpname.chr10.bp23423423 or snpname.10.23423423 to indicate chromosome and position. Or rs3691092.chr1.bp23027438.Morgans0.066997 to indicate both bp and cM/100 genetic map location.")
  
  cn.snpcols = sapply(snp.chr.pos,function(x) x[1])
  chr = sapply(snp.chr.pos,function(x) x[2])
  bp  = sapply(snp.chr.pos,function(x) x[3])
  if (has.bp.and.morgans & length(snp.chr.pos[[1]]) == 5) {
      Morgans = paste(sep="","0.",sapply(snp.chr.pos,function(x) x[5]))
  }
  
  chr2= strsplit(chr,split="chr",fixed=TRUE)
  bp2 = strsplit(bp,split="bp",fixed=TRUE)

  
  if (length(chr2[[1]]) > 1 ) { chr.good = sapply(chr2,function(x) x[2]) } else {chr.good = sapply(chr2,function(x) x[1])}
  if (length(bp2[[1]]) > 1 ) { bp.good = sapply(bp2,function(x) x[2]) } else {bp.good = sapply(bp2,function(x) x[1])}

  if (has.cm.not.bp) {
    cM = as.numeric(sapply(strsplit(cn[snpcols],split="cm",fixed=TRUE),function(x) x[2]))
  } else {
     cM = as.numeric(bp.good)/2e6 # mouse 1 cM = 2Mbp
   }

  if (has.bp.and.morgans) {
     cM = as.numeric(Morgans)*100
  }

  # cast them to numeric and throw out the markers without good info...

  chr.num = as.numeric(chr.good)
  cM.num = as.numeric(cM)

  good.snp = !(is.na(chr.num) & is.na(cM.num))

  # reduce snpcols
  snpcols = snpcols[good.snp]
  numsnp = length(snpcols)


  chr.good = chr.good[good.snp]
  cM = cM[good.snp]
  cn.snpcols = cn.snpcols[good.snp]

  mrow = nr+2
  m=matrix(nrow=nr+2,ncol=length(snpcols)+length(traitcols),data="")

  # convert chromosome 20 to X if requested (for mice!): must also include "sex" as a covariate
  if(convert.chr20.to.X) {
     if (is.na(match("sex",cn[traitcols]))) {
        warning("Didn't find 'sex' in the traitcols; skipping conversion of chr20 to chrX.")
     } else {
        w.20 = which(chr.good=="20")
        chr.good[w.20]="X"
     }
  }

  m[1,snp.start:(snp.start+numsnp-1)]=chr.good
  m[2,snp.start:(snp.start+numsnp-1)]=cM

  colnames(m) = c(cn[traitcols],cn.snpcols)

  m[3:mrow,snp.start:(snp.start+numsnp-1)] = sapply(data.frame(as.matrix(x[1:(mrow-2),snpcols])),translate.012.AHB)

  m[3:mrow, 1:nt] = as.matrix(x[1:(mrow-2),traitcols])

  write.table(m,file=file.name,sep=",",row.names=FALSE)

  print(paste("Wrote R/qtl CVS file to ",file.name," relative to getwd:",getwd()," at ",date()))
  print("Starting the read.cross() of the file")    
  rc=read.cross(format="csv", dir=".", file=file.name)

  print(paste("Done with read.cross() of the file",date()))
  invisible(rc)
}



#
# sample call:
#
# snp.snp.abs.cor.matrix = abs(cor(datCombined[,snpcols],use="p"))
#
# g=grafting(datCombined,snpcols,already.chosen.snpcols, genecol, traitcol,pm=pm,snp.snp.abs.cor.matrix=snp.snp.abs.cor.matrix)
#
# 
# 

# grafting:
#
# find snps from datCombined[,setdiff(snpcols,already.chosen.snpcols)]
#  that are independent of already.chosen.snpcols in that
#  they show correlation at a level of pm$cor.ind.th or less.
#
# Further check that these snps are independent of the genecol variable.
#
# If the snp passes these two tests, it is a good graft, and is returned.
# Otherwise return NA.
#
# 
if(exists("grafting")) rm(grafting);
grafting=function(datCombined,snpcols,already.chosen.snpcols,genecol, traitcol,pm=neo.get.param(),use.ranks=FALSE)
{
 if (length(traitcol) !=1) stop("only single traitcol allowed")
 if (length(genecol) !=1) stop("only single genecol allowed")
 if (length(already.chosen.snpcols) < 1) stop("Must specify at least 1 already.chosen.snpcols  to be grafting.")

 # should we convert the gene and trait cols to ranks first?
 if (use.ranks) {
    for (g in genecols) {
         datCombined[,g]=rank(datCombined[,g])
    }
    for (g in traitcol) {
         datCombined[,g]=rank(datCombined[,g])
    }
 }

 cn=colnames(datCombined)

 snp.cn = cn[snpcols]
 trait.cn = cn[traitcol]

 use.snps = setdiff(snpcols, already.chosen.snpcols)

 dc.snp = datCombined[,use.snps]
 dc.trait = datCombined[,traitcol]

 snp.trait.cor = cor(dc.snp, dc.trait, use="p")
 abs.trait.cor = abs(snp.trait.cor)

 # snps on rows, trait at the top of the column

 sort.snp =  sort(abs.trait.cor, decreasing=TRUE, index.return=TRUE)

  my.snp=list()

  for (i in sort.snp$ix) {

   my.snp$grafted.snp.index = use.snps[i]
   my.snp$grafted.snp = cn[my.snp$grafted.snp.index]
   my.snp$grafted.snp.trait.cor = snp.trait.cor[i]

   # my.snp has the next strongest SNP with this trait
   
   # is my.snp below pm$cor.dep.th correlated with all of the already.chosen.snpcols?

   #   w.row = match(my.snp$grafted.snp, rownames(snp.snp.abs.cor.matrix))
   check.me =  abs(cor(datCombined[,my.snp$grafted.snp.index], datCombined[,already.chosen.snpcols],use="p"))

   max.check.me = max(check.me)
   if (max.check.me < pm$cor.ind.th) {
      my.snp$closest.already.chosen.snp.to.grafted.snp = cn[already.chosen.snpcols[which.max(check.me)]]
      my.snp$closest.already.chosen.snp.to.grafted.snp.abs.cor = max.check.me

      # passed first check, now confirm that is also independent of genecol
      gene.new.snp.cor = abs(cor(datCombined[,my.snp$grafted.snp.index],datCombined[,genecol],use="p"))

      if (gene.new.snp.cor < pm$cor.ind.th) {
         my.snp$grafted.to.genecol.abs.cor = gene.new.snp.cor
         my.snp$genecol = cn[genecol]
         my.snp$traitcol = cn[traitcol]
         return(my.snp)
      }
   }

 } # end for i

 NA # default
} # end grafting

#
# do a bunch of grafting and fit the orthomarker models
#
if(exists("neo.graft")) rm(neo.graft)
neo.graft=function(datCombined,snpcols,already.chosen.snpcols,genecol, traitcol,pm=neo.get.param(),use.ranks=FALSE) {

  z=list()

  z$graft   = grafting(datCombined,snpcols,already.chosen.snpcols,genecol, traitcol,pm,use.ranks)

  MA.col = already.chosen.snpcols
  A.col  = genecol
  B.col  = traitcol
  MB.col = z$graft$grafted.snp.index

  z$orthomarker.fit = local.sem.four.var.m1m2(MA.col,A.col,B.col,x=datCombined,MB.col,no.obs=NULL,fit.models=c(1,2,3,4),ten.correct.chisq=FALSE,pm)

  z
}

################# Sun grid engine execution functions ##########################
################# Sun grid engine execution functions ##########################
################# Sun grid engine execution functions ##########################
################# Sun grid engine execution functions ##########################

# when running lots of jobs on the cluster, get a Sun Grid Engine task number
if(exists("get.sge.task")) rm(get.sge.task)
get.sge.task=function(quiet=TRUE) {
  task = Sys.getenv("SGE_TASK_ID")
  if (task=="" | is.na(task)) { ta = 1 } else {  ta = as.numeric(task) }
  if (is.na(ta)) { ta = 1 }
  if (!quiet) { print(paste("get.sge.task(): got SGE_TASK_ID task ",ta)) }
  ta
}

# same as above but return 0 if no SGE_TASK_ID
if(exists("get.sge.task.or.zero")) rm(get.sge.task.or.zero)
get.sge.task.or.zero=function(quiet=TRUE) {
  task = Sys.getenv("SGE_TASK_ID")
  if (task=="" | is.na(task) ) { ta = 0 } else {   ta = as.numeric(task) }
  if (is.na(ta)) { ta = 0}
  if (!quiet) { print(paste("get.sge.task(): got SGE_TASK_ID task ",ta)) }
  ta
}


if(exists("simple.test.sge.func")) rm(simple.test.sge.func)
simple.test.sge.func=function() {
   ta =  get.sge.task()
   print(paste("hello Sun Grid Engine!; task id:",ta))
}


#
# leave job.list empty to do a test run, otherwise give the function (no arguments, just the function
#     name) to execute for each task in the list
#
if(exists("submit.to.sge")) rm(submit.to.sge)
submit.to.sge=function(job.list=list(),neo.install.path="/home/jaten/dev/peculiar/neo.txt",r.program.install.path="/usr/local/bin/R") {

  num.tasks = length(job.list)
  if (num.tasks == 0) num.tasks = 2 # test situation
  
  if (num.tasks > 1) {
     cat(paste(sep="","qsub -cwd -t 1:",num.tasks," neo.sge.2.qsub.me.sh\n"), file = "neo.sge.1.source.me")
  } else {
     cat(paste(sep="","qsub -cwd  neo.sge.2.qsub.me.sh\n"), file = "neo.sge.1.source.me")
  }
  
  # biomath03.q defective, so we had to insist on submitting elsewhere...by naming all good queues as options
  #specify.allowable.queues="" # default
  specify.allowable.queues="#!/bin/sh\n#$ -q biomath02.q,biomath04.q,biomath05.q,cpu02.q,cpu03.q,cpu05.q,cpu06.q,cpu07.q,cpu08.q,cpu09.q,cpu10.q,cpu11.q,cpu12.q,cpu13.q,cpu14.q,cpu15.q\n"

   cat(paste(sep="",specify.allowable.queues,"\n",r.program.install.path," -q --vanilla < neo.sge.3.script.r\n"), file = "neo.sge.2.qsub.me.sh")

   cat(paste(sep="","source(\"",neo.install.path,"\"); task=get.sge.task(); if (is.na(task)) { task = 1 }; load(\"neo.sge.job.list.rdat\"); do.call(neo.sge.job.list[[task]],args=list());\n"), file = "neo.sge.3.script.r")

   neo.sge.job.list=list()

  if (length(job.list)==0) {
     # test situation
     neo.sge.job.list[[1]]="neo.test"
     neo.sge.job.list[[2]]="simple.test.sge.func"
  } else {
     # real jobs
     neo.sge.job.list = job.list
  }
  
   save(neo.sge.job.list,file="neo.sge.job.list.rdat")

   wd=getwd()

   the.cmd = paste(sep="","source ",wd,"/neo.sge.1.source.me")
   system(the.cmd)
}


################# END   Sun grid engine execution functions ##########################
################# END   Sun grid engine execution functions ##########################
################# END   Sun grid engine execution functions ##########################
################# END   Sun grid engine execution functions ##########################


# get.frames: used to extract useful data frames from the robustness studies files.
get.frames = function(rlist) {
  
  num.frames = length(rlist)
  blank.frame = matrix(nrow=num.frames, ncol=20)
  colnames(blank.frame) = paste(sep="","snp.count.",1:20)
  rownames(blank.frame) = paste(sep="","rep.",1:num.frames)  
  
  the.frames = list()
  
  the.frames$my.LEO.NB.FOR.snp = blank.frame
  the.frames$my.LEO.NB.MAX.snp = blank.frame
  the.frames$my.LEO.NB.ALL.snp = blank.frame
  the.frames$my.SIMPLE.MAX.MAX.snp = blank.frame
  the.frames$my.LEO.MULTIONE.NB.snp = blank.frame

  for (i in 1:num.frames) {
    
    df=data.frame(rlist[[i]])

    the.frames$my.LEO.NB.FOR.snp[i,]=anac(df$LEO.NB.FOR)
    the.frames$my.LEO.NB.MAX.snp[i,]=anac(df$LEO.NB.MAX)
    the.frames$my.LEO.NB.ALL.snp[i,]=anac(df$LEO.NB.ALL)
    the.frames$my.SIMPLE.MAX.MAX.snp[i,]=anac(df$SIMPLE.MAX.MAX)
    the.frames$my.LEO.MULTIONE.NB.snp[i,]=anac(df$LEO.MULTIONE.NB)   
  }
  
  the.frames
}

# work up the robustness study results
compare.robustness.studies=function() {

  my.27robust.files = c(
    "mc.robust.1a.1b/r.test.robust.task1.results.100.Sun_Aug_12_06.02.53_2007.rdat",
    "mc.robust.1a.1b/r.test.robust.task2.results.100.Sun_Aug_12_14.38.15_2007.rdat",
    "mc.robust.1a.1b/r.test.robust.task3.results.100.Sun_Aug_12_03.30.12_2007.rdat",
    
    "mc.robust.0a.1b.5causal/r.test.robust.task1.results.100.Tue_Aug_14_10.58.51_2007.rdat",
    "mc.robust.0a.1b.5causal/r.test.robust.task2.results.100.Sun_Aug_12_16.05.02_2007.rdat",
    "mc.robust.0a.1b.5causal/r.test.robust.task3.results.100.Sun_Aug_12_13.23.41_2007.rdat",
    
    "mc.robust.1a.1b.5causal/r.test.robust.task1.results.100.Sun_Aug_12_11.35.58_2007.rdat",
    "mc.robust.1a.1b.5causal/r.test.robust.task2.results.100.Sun_Aug_12_08.03.22_2007.rdat",
    "mc.robust.1a.1b.5causal/r.test.robust.task3.results.100.Sun_Aug_12_07.16.37_2007.rdat",
    
    "mc.robust.1a.1b.0causal/r.test.robust.task1.results.100.Thu_Aug_16_03.59.40_2007.rdat",
    "mc.robust.1a.1b.0causal/r.test.robust.task2.results.100.Sun_Aug_12_14.26.48_2007.rdat",
    "mc.robust.1a.1b.0causal/r.test.robust.task3.results.100.Sun_Aug_12_13.12.29_2007.rdat",
    
    "mc.robust.0a.1b.0causal/r.test.robust.task1.results.100.Mon_Aug_13_05.34.40_2007.rdat",
    "mc.robust.0a.1b.0causal/r.test.robust.task2.results.100.Sun_Aug_12_11.58.00_2007.rdat",
    "mc.robust.0a.1b.0causal/r.test.robust.task3.results.100.Mon_Aug_13_06.48.41_2007.rdat",
    
    "mc.robust.1a.1b.1causal.4an.4bn/r.test.robust.task1.results.100.Mon_Aug_13_18.29.13_2007.rdat",
    "mc.robust.1a.1b.1causal.4an.4bn/r.test.robust.task2.results.100.Mon_Aug_13_20.23.08_2007.rdat",
    "mc.robust.1a.1b.1causal.4an.4bn/r.test.robust.task3.results.100.Mon_Aug_13_19.32.40_2007.rdat",
    
    "mc.robust.1a.4b.4causal.4an.4bn/r.test.robust.task1.results.100.Tue_Aug_14_21.15.45_2007.rdat",
    "mc.robust.1a.4b.4causal.4an.4bn/r.test.robust.task2.results.100.Tue_Aug_14_12.20.56_2007.rdat",
    "mc.robust.1a.4b.4causal.4an.4bn/r.test.robust.task3.results.100.Fri_Aug_17_03.38.37_2007.rdat",
    
    "mc.robust.0a.4b.4causal.4an.4bn/r.test.robust.task1.results.100.Tue_Aug_14_14.09.52_2007.rdat",
    "mc.robust.0a.4b.4causal.4an.4bn/r.test.robust.task2.results.100.Tue_Aug_14_10.32.41_2007.rdat",
    "mc.robust.0a.4b.4causal.4an.4bn/r.test.robust.task3.results.100.Thu_Aug_16_22.19.47_2007.rdat",
    
    "mc.robust.0a.1b.no.neighbors.1causal/r.test.robust.task1.results.100.Sun_Aug_12_14.56.15_2007.rdat",
    "mc.robust.0a.1b.no.neighbors.1causal/r.test.robust.task2.results.100.Sun_Aug_12_10.06.59_2007.rdat",
    "mc.robust.0a.1b.no.neighbors.1causal/r.test.robust.task3.results.100.Sun_Aug_12_05.28.02_2007.rdat"
    )


### two functions used to summarize the robustness data files: compute.leo.stats() and summarize.robust.analysis.files()
  
#################### compute stats
    
    compute.leo.stats = function(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX) {
      
      LEO.NB.FOR.mean = mean(LEO.NB.FOR,na.rm=T)
      LEO.NB.ALL.mean = mean(LEO.NB.ALL,na.rm=T)
      LEO.NB.MAX.mean = mean(LEO.NB.MAX,na.rm=T)
      LEO.MULTIONE.NB.mean = mean(LEO.MULTIONE.NB,na.rm=T)
      SIMPLE.MAX.MAX.mean = mean(SIMPLE.MAX.MAX,na.rm=T)
      
      LEO.NB.FOR.se = stderr1(LEO.NB.FOR)
      LEO.NB.ALL.se = stderr1(LEO.NB.ALL)
      LEO.NB.MAX.se = stderr1(LEO.NB.MAX)
      LEO.MULTIONE.NB.se = stderr1(LEO.MULTIONE.NB)
      SIMPLE.MAX.MAX.se = stderr1(SIMPLE.MAX.MAX)  
      
      LEO.NB.FOR.sd = sd(LEO.NB.FOR,na.rm=T)
      LEO.NB.ALL.sd = sd(LEO.NB.ALL,na.rm=T)
      LEO.NB.MAX.sd = sd(LEO.NB.MAX,na.rm=T)
      LEO.MULTIONE.NB.sd = sd(LEO.MULTIONE.NB,na.rm=T)
      SIMPLE.MAX.MAX.sd = sd(SIMPLE.MAX.MAX,na.rm=T)
      
####################
      
      LEO.NB.FOR.mean = mean(LEO.NB.FOR,na.rm=T)
      LEO.NB.ALL.mean = mean(LEO.NB.ALL,na.rm=T)
      LEO.NB.MAX.mean = mean(LEO.NB.MAX,na.rm=T)
      LEO.MULTIONE.NB.mean = mean(LEO.MULTIONE.NB,na.rm=T)
      SIMPLE.MAX.MAX.mean = mean(SIMPLE.MAX.MAX,na.rm=T)
      
      LEO.NB.FOR.se = stderr1(LEO.NB.FOR)
      LEO.NB.ALL.se = stderr1(LEO.NB.ALL)
      LEO.NB.MAX.se = stderr1(LEO.NB.MAX)
      LEO.MULTIONE.NB.se = stderr1(LEO.MULTIONE.NB)
      SIMPLE.MAX.MAX.se = stderr1(SIMPLE.MAX.MAX)  
      
      LEO.NB.FOR.sd = sd(LEO.NB.FOR,na.rm=T)
      LEO.NB.ALL.sd = sd(LEO.NB.ALL,na.rm=T)
      LEO.NB.MAX.sd = sd(LEO.NB.MAX,na.rm=T)
      LEO.MULTIONE.NB.sd = sd(LEO.MULTIONE.NB,na.rm=T)
      SIMPLE.MAX.MAX.sd = sd(SIMPLE.MAX.MAX,na.rm=T)
      
      data.frame(LEO.NB.FOR.mean, LEO.NB.FOR.se, LEO.NB.FOR.sd, LEO.NB.ALL.mean, LEO.NB.ALL.se, LEO.NB.ALL.sd, LEO.NB.MAX.mean, LEO.NB.MAX.se, LEO.NB.MAX.sd, LEO.MULTIONE.NB.mean, LEO.MULTIONE.NB.se, LEO.MULTIONE.NB.sd, SIMPLE.MAX.MAX.mean, SIMPLE.MAX.MAX.se, SIMPLE.MAX.MAX.sd)
    } # end compute.leo.stats

  
  summarize.robust.analysis.files = function(my.27robust.files,f1,f2,f3,snp.number) {
  
    load(my.27robust.files[f1]) # r
    greedy.frames = get.frames(r$greedy.100.causal)

    load(my.27robust.files[f2])
    forward.frames = get.frames(r$forward.100.causal)
    
    load(my.27robust.files[f3])
    fwd.and.gr.frames = get.frames(r$greedy.and.forward.100.causal)
    
                                        # get the first SNP of each.
    LEO.NB.FOR = c(greedy.frames$my.LEO.NB.FOR.snp[,snp.number],forward.frames$my.LEO.NB.FOR.snp[,snp.number],fwd.and.gr.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(greedy.frames$my.LEO.NB.ALL.snp[,snp.number],forward.frames$my.LEO.NB.ALL.snp[,snp.number],fwd.and.gr.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(greedy.frames$my.LEO.NB.MAX.snp[,snp.number],forward.frames$my.LEO.NB.MAX.snp[,snp.number],fwd.and.gr.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(greedy.frames$my.LEO.MULTIONE.NB.snp[,snp.number],forward.frames$my.LEO.MULTIONE.NB.snp[,snp.number],fwd.and.gr.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(greedy.frames$my.SIMPLE.MAX.MAX.snp[,snp.number],forward.frames$my.SIMPLE.MAX.MAX.snp[,snp.number],fwd.and.gr.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    
    x = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)    
  x
  } # end summarize.robust.analysis.files


  
# Numbered file index
# > my.27robust.files 
#  [1] "mc.robust.1a.1b/r.test.robust.task1.results.100.Sun_Aug_12_06.02.53_2007.rdat"                     
#  [2] "mc.robust.1a.1b/r.test.robust.task2.results.100.Sun_Aug_12_14.38.15_2007.rdat"                     
#  [3] "mc.robust.1a.1b/r.test.robust.task3.results.100.Sun_Aug_12_03.30.12_2007.rdat"                     
# 
#  [4] "mc.robust.0a.1b.5causal/r.test.robust.task1.results.100.Tue_Aug_14_10.58.51_2007.rdat"             
#  [5] "mc.robust.0a.1b.5causal/r.test.robust.task2.results.100.Sun_Aug_12_16.05.02_2007.rdat"             
#  [6] "mc.robust.0a.1b.5causal/r.test.robust.task3.results.100.Sun_Aug_12_13.23.41_2007.rdat"             
# 
#  [7] "mc.robust.1a.1b.5causal/r.test.robust.task1.results.100.Sun_Aug_12_11.35.58_2007.rdat"             
#  [8] "mc.robust.1a.1b.5causal/r.test.robust.task2.results.100.Sun_Aug_12_08.03.22_2007.rdat"             
#  [9] "mc.robust.1a.1b.5causal/r.test.robust.task3.results.100.Sun_Aug_12_07.16.37_2007.rdat"             
# 
# [10] "mc.robust.1a.1b.0causal/r.test.robust.task1.results.100.Thu_Aug_16_03.59.40_2007.rdat"             
# [11] "mc.robust.1a.1b.0causal/r.test.robust.task2.results.100.Sun_Aug_12_14.26.48_2007.rdat"             
# [12] "mc.robust.1a.1b.0causal/r.test.robust.task3.results.100.Sun_Aug_12_13.12.29_2007.rdat"             
# 
# [13] "mc.robust.0a.1b.0causal/r.test.robust.task1.results.100.Mon_Aug_13_05.34.40_2007.rdat"             
# [14] "mc.robust.0a.1b.0causal/r.test.robust.task2.results.100.Sun_Aug_12_11.58.00_2007.rdat"             
# [15] "mc.robust.0a.1b.0causal/r.test.robust.task3.results.100.Mon_Aug_13_06.48.41_2007.rdat"             
# 
# [16] "mc.robust.1a.1b.1causal.4an.4bn/r.test.robust.task1.results.100.Mon_Aug_13_18.29.13_2007.rdat"     
# [17] "mc.robust.1a.1b.1causal.4an.4bn/r.test.robust.task2.results.100.Mon_Aug_13_20.23.08_2007.rdat"     
# [18] "mc.robust.1a.1b.1causal.4an.4bn/r.test.robust.task3.results.100.Mon_Aug_13_19.32.40_2007.rdat"     
# 
# [19] "mc.robust.1a.4b.4causal.4an.4bn/r.test.robust.task1.results.100.Tue_Aug_14_21.15.45_2007.rdat"     
# [20] "mc.robust.1a.4b.4causal.4an.4bn/r.test.robust.task2.results.100.Tue_Aug_14_12.20.56_2007.rdat"     
# [21] "mc.robust.1a.4b.4causal.4an.4bn/r.test.robust.task3.results.100.Fri_Aug_17_03.38.37_2007.rdat"     
# 
# [22] "mc.robust.0a.4b.4causal.4an.4bn/r.test.robust.task1.results.100.Tue_Aug_14_14.09.52_2007.rdat"     
# [23] "mc.robust.0a.4b.4causal.4an.4bn/r.test.robust.task2.results.100.Tue_Aug_14_10.32.41_2007.rdat"     
# [24] "mc.robust.0a.4b.4causal.4an.4bn/r.test.robust.task3.results.100.Thu_Aug_16_22.19.47_2007.rdat"     
# 
# [25] "mc.robust.0a.1b.no.neighbors.1causal/r.test.robust.task1.results.100.Sun_Aug_12_14.56.15_2007.rdat"
# [26] "mc.robust.0a.1b.no.neighbors.1causal/r.test.robust.task2.results.100.Sun_Aug_12_10.06.59_2007.rdat"
# [27] "mc.robust.0a.1b.no.neighbors.1causal/r.test.robust.task3.results.100.Sun_Aug_12_05.28.02_2007.rdat"

  #  1a.1b 100% causal, no neighbors, 
  x.1.2.3 = summarize.robust.analysis.files(my.27robust.files, 1,2,3,snp.number=1)
  x.4.5.6 = summarize.robust.analysis.files(my.27robust.files,4,5,6,snp.number=1)
  x.7.8.9 = summarize.robust.analysis.files(my.27robust.files,7,8,9,snp.number=1)
  x.10.11.12 = summarize.robust.analysis.files(my.27robust.files,10,11,12,snp.number=1)
  x.13.14.15 = summarize.robust.analysis.files(my.27robust.files,13,14,15,snp.number=1)
  x.16.17.18 = summarize.robust.analysis.files(my.27robust.files,16,17,18,snp.number=1)
  x.25.26.27 = summarize.robust.analysis.files(my.27robust.files, 25,26,27,snp.number=1)

  df= rbind(x.1.2.3, x.4.5.6, x.7.8.9, x.10.11.12, x.13.14.15, x.16.17.18, x.25.26.27)
  rownames(df) = c("1a.1b.100","0a.1b.50","1a.1b.50","1a.1b.0causal","0a.1b.0causal","1a.1b.100.4an.4bn","0a.1b.100")

  # get a better row order
  df2 = df[c(4,3,1,6,5,2,7),]
  write.csv(df2,file="single.snp.robustness.analysis.17aug2007.csv")
  

  # next we want to anlyze the two NULL models as the number of SNPs incorporated goes from 1 to 20

  f1=10
  f2=11
  f3=12 

  load(my.27robust.files[f1]) # r
  greedy.frames = get.frames(r$greedy.100.causal)
  
  load(my.27robust.files[f2])
  forward.frames = get.frames(r$forward.100.causal)
  
  load(my.27robust.files[f3])
  fwd.and.gr.frames = get.frames(r$greedy.and.forward.100.causal)
  
  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(greedy.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(greedy.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(greedy.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(greedy.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(greedy.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

  # plot each statistic under the null hypothesis

  for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
  for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)

  th = .3
  zscore = (th - for.mean) / for.sd
  for.fp = pnorm(zscore,lower.tail=FALSE)
  

  all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
  all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)

  barplot2(all.mean,plot.ci=T, ci.l=all.mean - all.sd/10, ci.u=all.mean+all.sd/10 ,ylim=c(-.7,.5),xlab="Greedy SNPs",ylab="LEO.NB.ALL, under Null",cex.axis=1.5,cex.names=1.5,cex.lab=1.5)
  box()


  # what choice of TH keeps the FP rate below 5% here...
  
  th.all = 3.5
  zscore.all = (th.all - all.mean) / all.sd
  all.fp = pnorm(zscore.all,lower.tail=FALSE)

  # how about the simple.max.max
 
  smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
  smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)

  barplot2(smaxmax.mean,plot.ci=T, ci.l=smaxmax.mean - smaxmax.sd/10, ci.u=smaxmax.mean+smaxmax.sd/10 ,ylim=c(-.7,.5),xlab="Greedy SNPs",ylab="MAX.MAX, under Null",cex.axis=1.5,cex.names=1.5,cex.lab=1.5)
  box()

############################
# start new analysis of the 4 snp scenario
# [19] "mc.robust.1a.4b.4causal.4an.4bn/r.test.robust.task1.results.100.Tue_Aug_14_21.15.45_2007.rdat"     
# [20] "mc.robust.1a.4b.4causal.4an.4bn/r.test.robust.task2.results.100.Tue_Aug_14_12.20.56_2007.rdat"     
# [21] "mc.robust.1a.4b.4causal.4an.4bn/r.test.robust.task3.results.100.Fri_Aug_17_03.38.37_2007.rdat"     
############################
  
f1=19
  f2=20
  f3=21 

  load(my.27robust.files[f1]) # r
  greedy.frames = get.frames(r$greedy.100.causal)
  
  load(my.27robust.files[f2])
  forward.frames = get.frames(r$forward.100.causal)
  
  load(my.27robust.files[f3])
  fwd.and.gr.frames = get.frames(r$greedy.and.forward.100.causal)
  
  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(greedy.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(greedy.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(greedy.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(greedy.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(greedy.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

gr.smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
gr.smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)
names(gr.smaxmax.mean)=names(gr.smaxmax.sd)=paste(1:20)

gr.all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
gr.all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)
names(gr.all.mean)=names(gr.all.sd)=paste(1:20)

gr.for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
gr.for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)
names(gr.for.mean)=names(gr.for.sd)=paste(1:20)

### forward

  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(forward.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(forward.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(forward.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(forward.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(forward.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

fw.smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
fw.smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)
names(fw.smaxmax.mean)=names(fw.smaxmax.sd)=paste(1:20)

fw.all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
fw.all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)
names(fw.all.mean)=names(fw.all.sd)=paste(1:20)

fw.for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
fw.for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)
names(fw.for.mean)=names(fw.for.sd)=paste(1:20)


### forward and greedy: fwd.and.gr.frames

  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(fwd.and.gr.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(fwd.and.gr.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(fwd.and.gr.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(fwd.and.gr.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(fwd.and.gr.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

fg.smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
fg.smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)
names(fg.smaxmax.mean)=names(fg.smaxmax.sd)=paste(1:20)

fg.all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
fg.all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)
names(fg.all.mean)=names(fg.all.sd)=paste(1:20)

fg.for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
fg.for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)
names(fg.for.mean)=names(fg.for.sd)=paste(1:20)

# plot LEO.NB.FOR

 for.mean.3 = rbind(gr.for.mean, fw.for.mean, fg.for.mean)
 for.3.ci.l = rbind(gr.for.mean - gr.for.sd/10, fw.for.mean - fw.for.sd/10, fg.for.mean - fg.for.sd/10)
 for.3.ci.u = rbind(gr.for.mean + gr.for.sd/10, fw.for.mean + fw.for.sd/10, fg.for.mean + fg.for.sd/10)
 rownames(for.3.ci.l) = rownames(for.mean.3)
 rownames(for.3.ci.u) = rownames(for.mean.3)

  barplot2(for.mean.3,plot.ci=T, ci.l=for.3.ci.l, ci.u=for.3.ci.u,xlab="Number of SNPs",ylab="LEO.NB.FOR scoring a 4 SNP trait",cex.axis=1.5,cex.names=1.5,cex.lab=1.5,beside=T,ylim=c(0,1.5),col=heat.colors(3))
  box()
  legend(36,1.4, c("Greedy SNP selection","Forward SNP selection", "Forward & Greedy SNP selection"), fill = heat.colors(3))
  abline(h=0.3,lty=2)


# plot LEO.NB.ALL

 all.mean.3 = rbind(gr.all.mean, fw.all.mean, fg.all.mean)
 all.3.ci.l = rbind(gr.all.mean - gr.all.sd/10, fw.all.mean - fw.all.sd/10, fg.all.mean - fg.all.sd/10)
 all.3.ci.u = rbind(gr.all.mean + gr.all.sd/10, fw.all.mean + fw.all.sd/10, fg.all.mean + fg.all.sd/10)
 rownames(all.3.ci.l) = rownames(all.mean.3)
 rownames(all.3.ci.u) = rownames(all.mean.3)

  barplot2(all.mean.3,plot.ci=T, ci.l=all.3.ci.l, ci.u=all.3.ci.u,xlab="Number of SNPs",ylab="LEO.NB.ALL scoring a 4 SNP trait",cex.axis=1.5,cex.names=1.5,cex.lab=1.5,beside=T,ylim=c(0,4.6),col=heat.colors(3))
  box()
  legend(33,3.6, c("Greedy SNP selection","Forward SNP selection", "Forward & Greedy SNP selection"), fill = heat.colors(3))
  abline(h=0.3,lty=2)

# simple MAX.MAX

 smaxmax.mean.3 = rbind(gr.smaxmax.mean, fw.smaxmax.mean, fg.smaxmax.mean)
 smaxmax.3.ci.l = rbind(gr.smaxmax.mean - gr.smaxmax.sd/10, fw.smaxmax.mean - fw.smaxmax.sd/10, fg.smaxmax.mean - fg.smaxmax.sd/10)
 smaxmax.3.ci.u = rbind(gr.smaxmax.mean + gr.smaxmax.sd/10, fw.smaxmax.mean + fw.smaxmax.sd/10, fg.smaxmax.mean + fg.smaxmax.sd/10)
 rownames(smaxmax.3.ci.l) = rownames(smaxmax.mean.3)
 rownames(smaxmax.3.ci.u) = rownames(smaxmax.mean.3)

  barplot2(smaxmax.mean.3,plot.ci=T, ci.l=smaxmax.3.ci.l, ci.u=smaxmax.3.ci.u,xlab="Number of SNPs",ylab="LEO.NB.FOR scoring a 4 SNP trait",cex.axis=1.5,cex.names=1.5,cex.lab=1.5,beside=T,ylim=c(0,5),col=heat.colors(3))
  box()
  legend(31,3, c("Greedy SNP selection","Forward SNP selection", "Forward & Greedy SNP selection"), fill = heat.colors(3))
  abline(h=4.21,lty=2)

#################
#################
####

############################
# start new analysis of the Null model 100% confounding, to show the SNP titrations...
# [10] "mc.robust.1a.1b.0causal/r.test.robust.task1.results.100.Thu_Aug_16_03.59.40_2007.rdat"             
# [11] "mc.robust.1a.1b.0causal/r.test.robust.task2.results.100.Sun_Aug_12_14.26.48_2007.rdat"             
# [12] "mc.robust.1a.1b.0causal/r.test.robust.task3.results.100.Sun_Aug_12_13.12.29_2007.rdat"             
############################
  
f1=10
  f2=11
  f3=12 

  load(my.27robust.files[f1]) # r
  greedy.frames = get.frames(r$greedy.100.causal)
  
  load(my.27robust.files[f2])
  forward.frames = get.frames(r$forward.100.causal)
  
  load(my.27robust.files[f3])
  fwd.and.gr.frames = get.frames(r$greedy.and.forward.100.causal)
  
  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(greedy.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(greedy.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(greedy.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(greedy.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(greedy.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

gr.smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
gr.smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)
names(gr.smaxmax.mean)=names(gr.smaxmax.sd)=paste(1:20)

gr.all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
gr.all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)
names(gr.all.mean)=names(gr.all.sd)=paste(1:20)

gr.for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
gr.for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)
names(gr.for.mean)=names(gr.for.sd)=paste(1:20)

### forward

  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(forward.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(forward.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(forward.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(forward.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(forward.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

fw.smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
fw.smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)
names(fw.smaxmax.mean)=names(fw.smaxmax.sd)=paste(1:20)

fw.all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
fw.all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)
names(fw.all.mean)=names(fw.all.sd)=paste(1:20)

fw.for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
fw.for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)
names(fw.for.mean)=names(fw.for.sd)=paste(1:20)


### forward and greedy: fwd.and.gr.frames

  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(fwd.and.gr.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(fwd.and.gr.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(fwd.and.gr.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(fwd.and.gr.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(fwd.and.gr.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

fg.smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
fg.smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)
names(fg.smaxmax.mean)=names(fg.smaxmax.sd)=paste(1:20)

fg.all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
fg.all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)
names(fg.all.mean)=names(fg.all.sd)=paste(1:20)

fg.for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
fg.for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)
names(fg.for.mean)=names(fg.for.sd)=paste(1:20)

# plot LEO.NB.FOR

 for.mean.3 = rbind(gr.for.mean, fw.for.mean, fg.for.mean)
 for.3.ci.l = rbind(gr.for.mean - gr.for.sd/10, fw.for.mean - fw.for.sd/10, fg.for.mean - fg.for.sd/10)
 for.3.ci.u = rbind(gr.for.mean + gr.for.sd/10, fw.for.mean + fw.for.sd/10, fg.for.mean + fg.for.sd/10)
 rownames(for.3.ci.l) = rownames(for.mean.3)
 rownames(for.3.ci.u) = rownames(for.mean.3)

  barplot2(for.mean.3,plot.ci=T, ci.l=for.3.ci.l, ci.u=for.3.ci.u,xlab="Number of SNPs",ylab="LEO.NB.FOR scoring a 100% confounded model",cex.axis=1.5,cex.names=1.5,cex.lab=1.5,beside=T,ylim=c(-5.5,1),col=heat.colors(3))
  box()
  legend(36,-3.6, c("Greedy SNP selection","Forward SNP selection", "Forward & Greedy SNP selection"), fill = heat.colors(3))
  abline(h=0.3,lty=2)


# plot LEO.NB.ALL

 all.mean.3 = rbind(gr.all.mean, fw.all.mean, fg.all.mean)
 all.3.ci.l = rbind(gr.all.mean - gr.all.sd/10, fw.all.mean - fw.all.sd/10, fg.all.mean - fg.all.sd/10)
 all.3.ci.u = rbind(gr.all.mean + gr.all.sd/10, fw.all.mean + fw.all.sd/10, fg.all.mean + fg.all.sd/10)
 rownames(all.3.ci.l) = rownames(all.mean.3)
 rownames(all.3.ci.u) = rownames(all.mean.3)

  barplot2(all.mean.3,plot.ci=T, ci.l=all.3.ci.l, ci.u=all.3.ci.u,xlab="Number of SNPs",ylab="LEO.NB.ALL scoring a 100% confounded model",cex.axis=1.5,cex.names=1.5,cex.lab=1.5,beside=T,ylim=c(-1,1),col=heat.colors(3))
  box()
  legend(25,-.5, c("Greedy SNP selection","Forward SNP selection", "Forward & Greedy SNP selection"), fill = heat.colors(3))
#  abline(h=0.3,lty=2)

# simple MAX.MAX

 smaxmax.mean.3 = rbind(gr.smaxmax.mean, fw.smaxmax.mean, fg.smaxmax.mean)
 smaxmax.3.ci.l = rbind(gr.smaxmax.mean - gr.smaxmax.sd/10, fw.smaxmax.mean - fw.smaxmax.sd/10, fg.smaxmax.mean - fg.smaxmax.sd/10)
 smaxmax.3.ci.u = rbind(gr.smaxmax.mean + gr.smaxmax.sd/10, fw.smaxmax.mean + fw.smaxmax.sd/10, fg.smaxmax.mean + fg.smaxmax.sd/10)
 rownames(smaxmax.3.ci.l) = rownames(smaxmax.mean.3)
 rownames(smaxmax.3.ci.u) = rownames(smaxmax.mean.3)

  barplot2(smaxmax.mean.3,plot.ci=T, ci.l=smaxmax.3.ci.l, ci.u=smaxmax.3.ci.u,xlab="Number of SNPs",ylab="MAX.MAX scoring a 100% confounded model",cex.axis=1.5,cex.names=1.5,cex.lab=1.5,beside=T,ylim=c(-1,1),col=heat.colors(3))
  box()
  legend(26,-.6, c("Greedy SNP selection","Forward SNP selection", "Forward & Greedy SNP selection"), fill = heat.colors(3))
  abline(h=4.21,lty=2)


##########################
##########################
##########################
#
# Repeat for one SNP trait
#
##########################
##########################
##########################



############################
# start new analysis of the Single marker model, to show the SNP titrations...
# [16] "mc.robust.1a.1b.1causal.4an.4bn/r.test.robust.task1.results.100.Mon_Aug_13_18.29.13_2007.rdat"     
# [17] "mc.robust.1a.1b.1causal.4an.4bn/r.test.robust.task2.results.100.Mon_Aug_13_20.23.08_2007.rdat"     
# [18] "mc.robust.1a.1b.1causal.4an.4bn/r.test.robust.task3.results.100.Mon_Aug_13_19.32.40_2007.rdat"     
############################
  
f1=16
  f2=17
  f3=18 

  load(my.27robust.files[f1]) # r
  greedy.frames = get.frames(r$greedy.100.causal)
  
  load(my.27robust.files[f2])
  forward.frames = get.frames(r$forward.100.causal)
  
  load(my.27robust.files[f3])
  fwd.and.gr.frames = get.frames(r$greedy.and.forward.100.causal)
  
  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(greedy.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(greedy.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(greedy.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(greedy.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(greedy.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

gr.smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
gr.smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)
names(gr.smaxmax.mean)=names(gr.smaxmax.sd)=paste(1:20)

gr.all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
gr.all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)
names(gr.all.mean)=names(gr.all.sd)=paste(1:20)

gr.for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
gr.for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)
names(gr.for.mean)=names(gr.for.sd)=paste(1:20)

### forward

  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(forward.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(forward.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(forward.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(forward.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(forward.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

fw.smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
fw.smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)
names(fw.smaxmax.mean)=names(fw.smaxmax.sd)=paste(1:20)

fw.all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
fw.all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)
names(fw.all.mean)=names(fw.all.sd)=paste(1:20)

fw.for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
fw.for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)
names(fw.for.mean)=names(fw.for.sd)=paste(1:20)


### forward and greedy: fwd.and.gr.frames

  x=list()
  for (snp.number in 1:20) {
    LEO.NB.FOR = c(fwd.and.gr.frames$my.LEO.NB.FOR.snp[,snp.number])
    LEO.NB.ALL = c(fwd.and.gr.frames$my.LEO.NB.ALL.snp[,snp.number])
    LEO.NB.MAX = c(fwd.and.gr.frames$my.LEO.NB.MAX.snp[,snp.number])
    LEO.MULTIONE.NB = c(fwd.and.gr.frames$my.LEO.MULTIONE.NB.snp[,snp.number])
    SIMPLE.MAX.MAX = c(fwd.and.gr.frames$my.SIMPLE.MAX.MAX.snp[,snp.number])
    
    x[[snp.number]] = compute.leo.stats(LEO.NB.FOR, LEO.NB.ALL, LEO.NB.MAX, LEO.MULTIONE.NB, SIMPLE.MAX.MAX)
  }

fg.smaxmax.mean = sapply(x, function(x) x$SIMPLE.MAX.MAX.mean)
fg.smaxmax.sd   = sapply(x, function(x) x$SIMPLE.MAX.MAX.sd)
names(fg.smaxmax.mean)=names(fg.smaxmax.sd)=paste(1:20)

fg.all.mean = sapply(x, function(x) x$LEO.NB.ALL.mean)
fg.all.sd   = sapply(x, function(x) x$LEO.NB.ALL.sd)
names(fg.all.mean)=names(fg.all.sd)=paste(1:20)

fg.for.mean = sapply(x, function(x) x$LEO.NB.FOR.mean)
fg.for.sd   = sapply(x, function(x) x$LEO.NB.FOR.sd)
names(fg.for.mean)=names(fg.for.sd)=paste(1:20)

# plot LEO.NB.FOR

 for.mean.3 = rbind(gr.for.mean, fw.for.mean, fg.for.mean)
 for.3.ci.l = rbind(gr.for.mean - gr.for.sd/10, fw.for.mean - fw.for.sd/10, fg.for.mean - fg.for.sd/10)
 for.3.ci.u = rbind(gr.for.mean + gr.for.sd/10, fw.for.mean + fw.for.sd/10, fg.for.mean + fg.for.sd/10)
 rownames(for.3.ci.l) = rownames(for.mean.3)
 rownames(for.3.ci.u) = rownames(for.mean.3)

  barplot2(for.mean.3,plot.ci=T, ci.l=for.3.ci.l, ci.u=for.3.ci.u,xlab="Number of SNPs",ylab="LEO.NB.FOR scoring a single SNP model, 4AN/4BN",cex.axis=1.5,cex.names=1.5,cex.lab=1.5,beside=T,col=heat.colors(3),ylim=c(0,3))
  box()
  legend(38,1.9, c("Greedy SNP selection","Forward SNP selection", "Forward & Greedy SNP selection"), fill = heat.colors(3))
  abline(h=0.3,lty=2)


# plot LEO.NB.ALL

 all.mean.3 = rbind(gr.all.mean, fw.all.mean, fg.all.mean)
 all.3.ci.l = rbind(gr.all.mean - gr.all.sd/10, fw.all.mean - fw.all.sd/10, fg.all.mean - fg.all.sd/10)
 all.3.ci.u = rbind(gr.all.mean + gr.all.sd/10, fw.all.mean + fw.all.sd/10, fg.all.mean + fg.all.sd/10)
 rownames(all.3.ci.l) = rownames(all.mean.3)
 rownames(all.3.ci.u) = rownames(all.mean.3)

  barplot2(all.mean.3,plot.ci=T, ci.l=all.3.ci.l, ci.u=all.3.ci.u,xlab="Number of SNPs",ylab="LEO.NB.ALL scoring a single SNP model, 4AN/4BN",cex.axis=1.5,cex.names=1.5,cex.lab=1.5,beside=T,col=heat.colors(3),ylim=c(0,5.5))
  box()
  legend(36,4.1, c("Greedy SNP selection","Forward SNP selection", "Forward & Greedy SNP selection"), fill = heat.colors(3))
  abline(h=4.53,lty=2)

# simple MAX.MAX

 smaxmax.mean.3 = rbind(gr.smaxmax.mean, fw.smaxmax.mean, fg.smaxmax.mean)
 smaxmax.3.ci.l = rbind(gr.smaxmax.mean - gr.smaxmax.sd/10, fw.smaxmax.mean - fw.smaxmax.sd/10, fg.smaxmax.mean - fg.smaxmax.sd/10)
 smaxmax.3.ci.u = rbind(gr.smaxmax.mean + gr.smaxmax.sd/10, fw.smaxmax.mean + fw.smaxmax.sd/10, fg.smaxmax.mean + fg.smaxmax.sd/10)
 rownames(smaxmax.3.ci.l) = rownames(smaxmax.mean.3)
 rownames(smaxmax.3.ci.u) = rownames(smaxmax.mean.3)

  barplot2(smaxmax.mean.3,plot.ci=T, ci.l=smaxmax.3.ci.l, ci.u=smaxmax.3.ci.u,xlab="Number of SNPs",ylab="MAX.MAX scoring a single SNP model, 4AN/4BN",cex.axis=1.5,cex.names=1.5,cex.lab=1.5,beside=T,col=heat.colors(3),ylim=c(0,5))
  box()
  legend(26,2.2, c("Greedy SNP selection","Forward SNP selection", "Forward & Greedy SNP selection"), fill = heat.colors(3))
  abline(h=4.21,lty=2)

} # end compare.robustness.studies


if(exists("sim.loop")) rm(sim.loop)
sim.loop=function() {

   set.seed(1)
   N=200

   SNP.A = rnorm(N)
   SNP.B = rnorm(N)
   SNP.C = rnorm(N)

   time.pts = 50
   A = B = C = matrix(NA,nrow=N,ncol=time.pts)   
   
   # initial values
   A[,1] = .2 + rnorm(N)
   B[,1] = .2 + rnorm(N)
   C[,1] = .2 + rnorm(N) 

   # path coefficients
   a.to.b = .3
   b.to.c = .3
   c.to.a = .3

   # find an equilibrium
   for (i in 2:time.pts) {
     A[,i] = c.to.a * C[,i-1] + SNP.A # + rnorm(N,sd=.1)
     B[,i] = a.to.b * A[,i-1] + SNP.B # + rnorm(N,sd=.1)
     C[,i] = b.to.c * B[,i-1] + SNP.C # + rnorm(N,sd=.1)

     print(paste(sep="","At i=",i,"   A[1]=",A[1,i],"   B[1]=",B[1,i],"   C[1]=",C[1,i]))
   }

   col.A = A[,time.pts]
   col.B = B[,time.pts]
   col.C = C[,time.pts]

   x=data.frame(SNP.A,SNP.B,SNP.C,col.A,col.B,col.C)

   pm=neo.get.param()
   pm$run.title="sim.loop"
   z=neo(x,snpcols=1:3)

   z
}


auto.snp.insig1 = function() {

#setwd("c:/good")
source("neo.txt")
a=load("liver.snps.23388genes.clinical.bxh.male.and.female.rdat")
# a contains two variables:
a
# [1] "str.me"                "liver.bxh.male.female"

#Let's look at the dimensions of the liver data:
dim(liver.bxh.male.female)
#[1]   265 24714

str(str.me)
#List of 4
# $ snpcosl      : int [1:1278] 1 2 3 4 5 6 7 8 9 10 ...
# $ genecols     : int [1:23388] 1279 1280 1281 1282 1283 1284 1285 #1286 1287 1288 ...
# $ sex          : num 24667
# $ clinical.cols: int [1:47] 24668 24669 24670 24671 24672 24673 24674 #24675 24676 24677 ...

# We identify the indices of several genes from the expression data that we are interested 
# in.
which.Insig1=pmatch("Insig1",colnames(liver.bxh.male.female))
which.Fdft1=pmatch("Fdft1.MMT00082285",colnames(liver.bxh.male.female))
which.Dhcr7=pmatch("Dhcr7",colnames(liver.bxh.male.female))

# datCombined contains the data we will analyze
female=liver.bxh.male.female$sex==2
datCombined=liver.bxh.male.female[female,]


downstream.of.insig1 = c("Acac.MMT00048501","Acas2.MMT00009506","Adipor2.MMT00046649", "Aqp8.MMT00004671","B430110G05Rik.MMT00004109","BC036563.MMT00012566","C130074G19Rik.MMT00017549","Cs.MMT00082255","Dbi.MMT00068995","Dhcr7.MMT00053523","Dia1.MMT00026453","Dlat.MMT00009011","Eaf2.MMT00060307","Elovl6.MMT00035545","Fads1.MMT00046994","Fads2.MMT00012675","Fasn.MMT00018029","Fdft1.MMT00082285","Fdps.MMT00063985","Gale.MMT00077659","Map17.MMT00058422","MMT00012966","MMT00060232","Mod1.MMT00046920","Mx1.MMT00026973","Pgd.MMT00038022","Qdpr.MMT00038017","Rdh11.MMT00006408","Sc4mol.MMT00059097","Scd1.MMT00022878","Slc23a1.MMT00022094","Slc25a1.MMT00054225","Slc41a3.MMT00029681","Stard4.MMT00024116","X0610030G03Rik.MMT00036912","X6030440G05Rik.MMT00031386","Zdhhc6.MMT00047560")

# So first a NEO multimarker analysis.
traits=cbind(datCombined[,which.Insig1],datCombined[,downstream.of.insig1])
colnames(traits)[1]= colnames(datCombined)[which.Insig1]
datC = data.frame(cbind(datCombined[,str.me$snpcols],traits))

pm=neo.get.param()
pm$A = 1279
pm$B=c(match("Dhcr7.MMT00053523",colnames(datC)),match("Fdft1.MMT00082285",colnames(datC)))

# impute medians of missing trait data...
pm$rough.and.ready.NA.imputation = TRUE
pm$quiet=FALSE
pm$run.title="insig1.to.dhcr7.fdft1.auto.snp.select"
x=neo(datC, pm=pm, snpcols=1:1278)
save(x,file="auto.snp.dhcr7.fdft1.rdat")

}



if(exists("error.returns.NA")) rm(error.returns.NA)
error.returns.NA=function(x) {
  NA
}

if(exists("fail.proof.shapiro.test")) rm(fail.proof.shapiro.test)
fail.proof.shapiro.test=function(x) {
  op=options()
  options(error=error.returns.NA) # turn off error debuggin during this function
  r=try(shapiro.test(x),silent=TRUE)
  options(op) #restore original options; debuggin error stuff, etc
  if(inherits(r, "try-error")) return(NA)
  return(r$p.value)
}


# try-wrapper for generic function FUNC:
# return NA if the function FUNC fails when applied to x
#
if(exists("fail.proof.FUN")) rm(fail.proof.FUN)
fail.proof.FUN=function(FUNC,...) {
  op=options()
  options(error=error.returns.NA) # turn off error debuggin during this function
  options(warn=-1) # turn off warnings
  r=try(FUNC(...),silent=TRUE)
  options(op) #restore original options; debuggin error stuff, etc
  if(inherits(r, "try-error")) return(NA)
  r
}


if(exists("create.permutation.report")) rm(create.permutation.report)
create.permutation.report=function(pm,cn.datCombined,snpcols) {

   # read back in the permutations...and report the distribution.
   perms = read.csv(file=pm$perm.file) # let rownames become their own column, since they will be repeated
                  # and non-unique rownames might be a problem.

   perms[,1]=as.character(perms[,1]) # avoid factor messes

   # and summarize
   aur=unique(perms[,1]) # first column has edge names
   num.perm.summary.rows = length(aur)

   perm.mean = perms[1:num.perm.summary.rows,]
   perm.mean[,1:3]="" # wipe annotation
   perm.mean[,4:ncol(perms)] = NA
   perm.mean[,1] = aur

   # and fix up from and to
   splt=strsplit(aur,split=" -> ",fixed=TRUE)
   perm.mean[,2]= FROM.A = sapply(splt,function(x) x[1])
   perm.mean[,3]= TO.B   = sapply(splt,function(x) x[2])

   # and eliminitate rows that have snpcols in FROM.A
   w.from.snp = which(FROM.A %in% cn.datCombined[snpcols])
   aur = aur[-w.from.snp]
   perm.mean = perm.mean[-w.from.snp,]

   # get the observed
   obs.stat = read.csv(file=get.excel.file(pm),row.names=1)
   mt = match(colnames(perm.mean)[-c(1:3)],colnames(obs.stat))
   mt.rows = match(perm.mean[,1], rownames(obs.stat))
   if (any(is.na(mt.rows))) { stop("bad correspondence between permutation data and observed LEO scores.")}

   # make the observed match up with our permutations, except it lacks the first 3 annotation columns
   obs.stat.use = data.frame(cbind(perm.mean[,1:3],obs.stat[aur,mt]))

   # see if we can replace this next code with a common call to create.permutation.report.sma

   report = create.permutation.report.sma(pm,cn.datCombined,snpcols,perms,skipto=4,num.stats=17,skipable=c(1:3),p.file=paste(sep="",pm$neo.log.file,".collated.permutation.results.csv"), aur, perm.mean, obs.stat.use)

  if (0) { # while trying out create.permutation.report.sma() instead

   # summary of good (not na) data
   perm.percent.good = perm.mean

   # sd
   perm.sd = perm.mean

   # percentage of time permutated statistic exceeded observed, empirically observed.
   perm.pvalue = perm.mean

   # theoretical p-value, based on assuming the null distribution of the statistics is Gaussian.
   perm.theory.pvalue = perm.mean

   # shapiro.test() results for Shapiro-Wilk test of normality
   perm.shapiro = perm.mean


   # This may or may not be right to do, but simplifies the reporting for positive statistics.
   # --> Turn NA in permutation data into zero.
   ## perms[is.na(perms)]=0

   for (i in 1:length(aur)) {
      cur.row = aur[i]
      w = which(cur.row == perms[,1]) # first column has edge names
      if (length(w)>1) {
         ##### FIXME: the line below appears wrong; did the author mean colMeans instead of mean?
         perm.mean[i,4:ncol(perms)] = mean(perms[w,4:ncol(perms)],na.rm=TRUE)
         perm.sd[i,4:ncol(perms)] = colSds(as.matrix(perms[w,4:ncol(perms)]),na.rm=TRUE)
         perm.percent.good[i,4:ncol(perms)]=apply(!is.na.or.nan(perms[w,4:ncol(perms)]),2,sum)/pm$number.BLOCK.permutations
         perm.shapiro[i,4:ncol(perms)] = apply(perms[w,4:ncol(perms)],2,fail.proof.shapiro.test)

         #perm.pvalue computed depends on the observed statistic in that column:
         for (j in 4:ncol(perms)) {
             perm.pvalue[i,j] = sum(perms[w,j]>= obs.stat.use[i,j],na.rm=TRUE)/pm$number.BLOCK.permutations
             perm.theory.pvalue[i,j] = pnorm(obs.stat.use[i,j],mean=perm.mean[i,j],sd=perm.sd[i,j],lower.tail=FALSE)
         }
      } else {
         # handle special case of just one permutation gracefully
         perm.mean[i,4:ncol(perms)] = perms[w,4:ncol(perms)]
         perm.sd[i,4:ncol(perms)] = Inf #
         perm.percent.good[i,4:ncol(perms)]=apply(!is.na.or.nan(perms[w,4:ncol(perms)]),2,sum)/pm$number.BLOCK.permutations
         perm.pvalue[i,4:ncol(perms)] = NA # can't compute p-value given just a single permutation.
         perm.theory.pvalue[i,4:ncol(perms)] = NA # can't compute p-value given just a single permutation.
         perm.shapiro[i,4:ncol(perms)] = NA
      }
   }

#   write.csv(perm.mean,file=paste(sep="",pm$neo.log.file,".perm.means.csv"))
#   write.csv(perm.sd,file=paste(sep="",pm$neo.log.file,".perm.sd.csv"))
#   write.csv(perm.pvalue,file=paste(sep="",pm$neo.log.file,".perm.pvalue.csv"))
#   write.csv(perm.theory.pvalue,file=paste(sep="",pm$neo.log.file,".perm.pvalue.csv"))
#   write.csv(perm.percent.good,file=paste(sep="",pm$neo.log.file,".perm.summary.percent.good.csv"))

   # interleaf the results to provide a summary including p-value
   colnames(perm.pvalue)[-c(1:3)] = paste(sep="","pval.empirical.",colnames(perm.pvalue)[-c(1:3)])
   colnames(perm.theory.pvalue)[-c(1:3)] = paste(sep="","pval.theory.",colnames(perm.theory.pvalue)[-c(1:3)])
   colnames(perm.sd)[-c(1:3)] = paste(sep="","sd.",colnames(perm.sd)[-c(1:3)])
   colnames(perm.mean)[-c(1:3)] = paste(sep="","mean.",colnames(perm.mean)[-c(1:3)])
   colnames(perm.percent.good)[-c(1:3)] = paste(sep="","percent.good.",colnames(perm.percent.good)[-c(1:3)])
   colnames(perm.shapiro)[-c(1:3)] = paste(sep="","shapiro.pval.",colnames(perm.shapiro)[-c(1:3)])

   space = data.frame(cbind(rep(" ",nrow(obs.stat.use))))

   interleaf = data.frame(cbind(obs.stat.use, space, perm.pvalue[,-c(1:3)], space, perm.theory.pvalue[,-c(1:3)], space, perm.mean[,-c(1:3)], space, perm.sd[,-c(1:3)], space, perm.percent.good[,-c(1:3)], space, perm.shapiro, space))

   cn.il = colnames(interleaf)

   i = 0
   for (k in 4:ncol(perms)) {
      i=i+1
      j = (k-4)*8+4
      interleaf[,j] = obs.stat.use[,k]
      interleaf[,j+1] = perm.pvalue[,k]
      interleaf[,j+2] = perm.theory.pvalue[,k]
      interleaf[,j+3] = perm.mean[,k]
      interleaf[,j+4] = perm.sd[,k]
      interleaf[,j+5] = perm.percent.good[,k]
      interleaf[,j+6] = perm.shapiro[,k]
      interleaf[,j+7] = " "

      cn.il[j]=colnames(obs.stat.use)[k]
      cn.il[j+1]=colnames(perm.pvalue)[k]
      cn.il[j+2]=colnames(perm.theory.pvalue)[k]
      cn.il[j+3]=colnames(perm.mean)[k]
      cn.il[j+4]=colnames(perm.sd)[k]
      cn.il[j+5]=colnames(perm.percent.good)[k]
      cn.il[j+6]=colnames(perm.shapiro)[k]
      cn.il[j+7]=paste(rep("_",i+3),collapse="")
   }

   colnames(interleaf) = cn.il
   #p.file=paste(sep="",pm$neo.log.file,".collated.permutation.results.csv")
   write.csv(interleaf,file=p.file)
   if (pm$open.excel.at.end) open.excel.file(p.file,pm)
  } # end if(0)

}


if(exists("get.computer.name")) rm(get.computer.name)
get.computer.name=function(pm,cn.datCombined,snpcols) {
  if (.Platform$OS.type == "windows") { 
    return(try(Sys.getenv("COMPUTERNAME"),silent=TRUE))
  } 
  a=try(system("uname -n",intern=TRUE),silent=TRUE)
}

if (!exists(".Global.NEO.version")) {
   cat(paste("NEO loaded at ",date()," ...on ",get.computer.name(),"\n\n"))
   cat("Hint: Ctrl-W toggles console updating (on Windows); to allow you to\n     monitor the progress of a NEO run.\n\n")
   .Global.NEO.version <<- neover() # side effect: prints version
} else {
   ver=neover(quiet=TRUE)
   if (ver != .Global.NEO.version) {
      .Global.NEO.version <<- ver
      cat(paste("NEO loaded at ",date()," ...on ",get.computer.name(),"\n"))
      cat(.Global.NEO.version)
      cat("\n")
   }
}

