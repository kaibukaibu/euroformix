#' @title calcCommonLR
#' @author Oyvind Bleka
#' @description Calculating the LR for common unknowns
#' @details LR calculation of whether two paired unknowns between two fitted models is identical (or not).
#' Based on formula by Slooten (2017): Identifying common donors in DNA mixtures, with applications to database searches
#' @param mlefit1 1st Fitted object using calcMLE/contLikMLE function.
#' @param mlefit2 2nd Fitted object using calcMLE/contLikMLE function.
#' @param DC1 Optional deconvolute object from mlefit1
#' @param DC2 Optional deconvolute object from mlefit2
#' @param returnPerMarker Whether returning LR per marker instead of across all.
#' @return A matrix of LR between all unknown pairs. NULL if some have no unknowns. A ists if returnPerMarker=TRUE
#' @export

calcCommonLR = function(mlefit1,mlefit2,DC1=NULL,DC2=NULL,returnPerMarker=FALSE){ 
  Qallele="99" #can be any other than those observed
  if(is.null(DC1)) DC1 = deconvolve(mlefit1,alpha = 1) #ensure to list all outcomes
  if(is.null(DC2)) DC2 = deconvolve(mlefit2,alpha = 1)  #ensure to list all outcomes

  getGenoMatHelper = function(X) {
    return(t(matrix(unlist(strsplit(X,"/")),nrow=2)))
  }          
  #genoSorter = function(X) sapply(strsplit(X,"/"),function(x) paste0(sort(x),collapse="/"))
  
  #Obtain list of marginal probs (already performed in MLE step)
  DCmarg1 = DC1$table3 #get marginal probs
  DCmarg2 = DC2$table3 #get marginal probs
  c1 <- mlefit1$prepareC #returned from prepareC
  c2 <- mlefit2$prepareC #returned from prepareC

  #Obtain necessary data from objects  
  nU1 = c1$nU #number of unknowns in model1
  nU2 = c2$nU #number of unknowns in model2
  if(nU1==0 || nU2==0) return(NULL) #NOT POSSIBLE TO CALCULATE
  locs1 = c1$markerNames
  locs2 = c2$markerNames
  commonLocs = intersect(locs1,locs2)
  fstVec = c1$fst #use fst from 1st model
  
  DCmarg1[,1] = gsub("C","",DCmarg1[,1]) #remove letter prefix letter C
  DCmarg2[,1] = gsub("C","",DCmarg2[,1]) #remove letter prefix letter C
  nK1 = c1$NOK #number of knowns (conditionals)
  nK2 = c2$NOK #number of knowns (conditionals)
  Urng1 = seq_len(nU1) #range of unknowns
  Urng2 = seq_len(nU2) #range of unknowns
  Uname1 = paste0("U",Urng1)
  Uname2 = paste0("U",Urng2)
  KUrng1 = setNames(nK1 + Urng1,Uname1) #range of unknowns adjusted for knowns
  KUrng2 = setNames(nK2 + Urng2,Uname2) #range of unknowns adjusted for knowns 
  
  #Traverse each marker
  commonLR = matrix(0,ncol=nU2,nrow=nU1,dimnames = list(Uname1,Uname2))  
  commonLRmarker = matrix(list(),ncol=nU2,nrow=nU1,dimnames = list(Uname1,Uname2))  
  for(loc in commonLocs) {
  #  loc=  commonLocs[16]
    
    #Get marker information for each model:
    genoMat1 = c1$genoList[[loc]]
    genoMat2 = c2$genoList[[loc]]
    mind1 = which(locs1==loc) #get marker index
    mind2 = which(locs2==loc) #get marker index
    mrng1 = c1$startIndMarker_nAlleles[mind1] + seq_len(c1$nAlleles[mind1])
    mrng2 = c2$startIndMarker_nAlleles[mind2] + seq_len(c2$nAlleles[mind2])
    alleles1 = c1$alleleNames[mrng1]
    alleles2 = c2$alleleNames[mrng2]
    freq1 = setNames(c1$freqs[mrng1],alleles1) #frequencies
    freq2 = setNames(c2$freqs[mrng2],alleles2) #frequencies 
    known1counts = setNames(c1$maTyped[mrng1],alleles1)
    known2counts = setNames(c2$maTyped[mrng2],alleles2)
    known1Alelles = rep(names(known1counts),known1counts)
    known2Alelles = rep(names(known2counts),known2counts)
    knownAlleles = c(c(known1Alelles),c(known2Alelles)) #obtain all counts
    
    #OBTAIN COMMON STUFF BETWEEN THE MODELS 
    allelesCommon = intersect(alleles1,alleles2)
    if(length(knownAlleles)>0) knownAlleles[!knownAlleles%in%allelesCommon] = Qallele
    
    #Prepare frequencies
    freqsCommon = c(freq1[allelesCommon],freq2[allelesCommon])
    freqsCommon = freqsCommon[!duplicated(names(freqsCommon))]
    freqsCommonNoQ = freqsCommon[names(freqsCommon)!=Qallele] #keep only those not Q
    if(length(freqsCommonNoQ)==0) {
      freqsCommon = setNames(1,Qallele) #assign
    } else {
      freqSum = sum(freqsCommonNoQ)
      if(freqSum>1) freqsCommonNoQ = freqsCommonNoQ/freqSum #normalize to not exceed 1
      Qfreq = setNames(1-sum(freqsCommonNoQ),Qallele)
      freqsCommon = c(freqsCommonNoQ,Qfreq) #always Add Q alleles as a possibility
    }
    order = order(names(freqsCommon),decreasing = FALSE) #sort wrt names
    freqsCommon = freqsCommon[order] #update order
    fstMarker = fstVec[mind1] #obtain fst for 1st model
    
    #Obtain outcome of genotypes (including frequencies) across both models
    GlistMarker = calcGjoint(freq = freqsCommon, nU = 1, fst = fstMarker, refK = knownAlleles)
    genosUnique = GlistMarker$G
    nGenos = nrow(genosUnique) 
    
    #traversing all pairwise combinations of unknowns (and store result per marker)
    genosUniqueCombined = paste0(genosUnique[,1],"/",genosUnique[,2]) #combined
    genosProb = setNames(GlistMarker$Gprob,genosUniqueCombined)
    
    #Prepare genotypes in DC table to be same as for the "common genotype" outcome
    DCtab1Marker = DCmarg1[DCmarg1[,2]==loc,,drop=FALSE]
    DCtab2Marker = DCmarg2[DCmarg2[,2]==loc,,drop=FALSE]
      
    #get updated genotype names due to having overlap
    getUpdatedGenos = function(genoVec) {
      genoMat = getGenoMatHelper(genoVec)
      swap = genoMat[,2]<genoMat[,1]
      if(any(swap)) genoMat[swap,] = genoMat[swap,2:1]
      for(g in 1:2) genoMat[!genoMat[,g]%in%genosUnique[,g],g] = Qallele
      genoVec = paste0(genoMat[,1],"/",genoMat[,2])
      return(genoVec)
    }
    updatedGenos1 = getUpdatedGenos(DCtab1Marker[,3]) 
    updatedGenos2 = getUpdatedGenos(DCtab2Marker[,3])
    
    for(u1 in Urng1) {
      for(u2 in Urng2) {
        #loc=commonLocs[1];
        #u1=2;u2=3
        indRows1 = DCtab1Marker[,1]==KUrng1[u1] 
        indRows2 = DCtab2Marker[,1]==KUrng2[u2] 
        
        #Obtain genotypes from DC (updated version)
        genosDC1 = updatedGenos1[indRows1] #obtain geno model 1
        genosDC2 = updatedGenos2[indRows2] #obtain geno model 2
        probs1 = as.numeric(DCtab1Marker[indRows1,4])
        probs2 = as.numeric(DCtab2Marker[indRows2,4])
        
        #Traverse all unique genotypes
        sumProd = 0 #get final sum
        for(genoIdx in 1:nGenos) {
        #  genoIdx=10
          priorProb = genosProb[genoIdx] #get genotype prior prob
          if(priorProb==0) next #skip if this is not possible (apriori)
          genoPropCombined = genosUniqueCombined[genoIdx] #proposed genotype
          
          #Find correct row for each DC
          DCrows1 = genosDC1==genoPropCombined
          DCrows2 = genosDC2==genoPropCombined
          prob1 = sum(probs1[DCrows1])
          prob2 = sum(probs2[DCrows2])
          term = prob1*prob2/priorProb
          sumProd = sumProd + term
          #print(genoPropCombined)
          #print(term)
        } #end all genotypes
        commonLR[u1,u2] = commonLR[u1,u2] + log10(sumProd)
        newListElem = setNames(log10(sumProd),loc)
        commonLRmarker[[u1,u2]] = c(commonLRmarker[[u1,u2]],newListElem) #add to list
      } #end u1
    }#end u2
  } #end for each marker
 # barplot(commonLRmarker[[1,2]])
  if(returnPerMarker) {
    return(commonLRmarker)
  } else {
    return(commonLR)
  }
} #end function