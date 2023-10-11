##metabolomics integration functions
##here we store some functions we've been developing for metabolomics integration


if(!require('KEGGgraph')){
  BiocManager::install("KEGGgraph")
  library(KEGGgraph)
}

##gets all kegg pathways that contain nodes in the nodelist
getMetNetworksForNodes<-function(nodelist,allpaths){
  allnodes<-lapply(allpaths,function(x){
     # print(x)
    tmp <- tempfile()
    retrieveKGML(x, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
    mapkG<-data.frame()
    try(mapkG <- parseKGML2DataFrame(tmp,expandGenes=TRUE,genesOnly=FALSE))
    if(nrow(mapkG)>0)
      mapkG$pathway=x
    if(length(intersect(union(mapkG$from,mapkG$to),nodelist))>0)
      return(mapkG)

  })
  fullnet<-do.call(rbind,allnodes)
  return(fullnet)

}


##this leverages the same logic that we have for the phospho data
##basically we take the network with gene names and compounds and adds them to the ppi network
addNetworkToProts<-function(metnet,exprProts=c(),
                            filterTypes=c('expression','ubiquitination','indirect effect','missing interaction',NA)){
  require(PCSF)
  data(STRING)
  if(length(exprProts)>0)
    STRING=STRING|>
    subset(from%in%exprProts)|>
    subset(to%in%exprProts)


  mval=mean(STRING$cost)
  newmet<-subset(metnet,!subtype%in%(filterTypes))
  res<-apply(newmet,1,function(x){
    fg=x[['newfrom']]
    tg=x[['newto']]
    ss=subset(STRING,from==fg)|>
      subset(to==tg)
    #print(ss)
    if(nrow(ss>0))
      return(list(from='',to='',cost=0)) ##basically the interaction is already in the network, do not add
    newx=list(from=fg,to=tg,cost=mval/2)
    return(newx)
  })

  adf = do.call(rbind,res)|>
    as.data.frame()|>
    subset(cost!=0)
  ##first get diffex proteins
  ppi <- construct_interactome(rbind(STRING,adf))
  return(ppi)
}



mapNetworkToHGNC<-function(fullNetwork,kegg_met){
  require(org.Hs.eg.db)

  ##now we map gene identifers and metabolit names
  fromids<-translateKEGGID2GeneID(fullNetwork$from)
  geneids<-sapply(mget(fromids, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
  toids<-translateKEGGID2GeneID(fullNetwork$to)
  tgeneids<-sapply(mget(toids,org.Hs.egSYMBOL,ifnotfound=NA),'[[',1)

  ##here is our complete mapping
  mapping<-data.frame(kegg=c(fullNetwork$from,fullNetwork$to),
                      gid=c(fromids,toids),
                      symbol=c(geneids,tgeneids))|>
    distinct()|>
    subset(!is.na(symbol))|>
    rbind(data.frame(kegg=kegg_met$compoundId,gid=kegg_met$KEGG,symbol=kegg_met$Metabolite))

  #the join introduces NAs but does a decent job
  newNetwork<-fullNetwork|>
    left_join(dplyr::rename(mapping,from='kegg',newfrom='symbol'))|>
    dplyr::select(-gid)|>
    left_join(dplyr::rename(mapping,to='kegg',newto='symbol'))|>
    dplyr::select(-gid)

  ##TODO make this metabolites but for now just use Kegg ids
  fixedNetwork<-apply(newNetwork,1,function(x){
    if(is.na(x[['newto']])){
      x[['newto']]=x[['to']]
    }
    if(is.na(x[['newfrom']])){
      x[['newfrom']]=x[['from']]
    }
    return(x)
  })|>
    t()|>
    as.data.frame()
  return(fixedNetwork)
}
