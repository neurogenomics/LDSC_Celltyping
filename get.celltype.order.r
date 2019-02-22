# get.celltype.order
# - Used when generating figures from LDSC results. Returns the order of celltypes (w.r.t. a particular gwas) based on either z-score ('z'), p-val ('p') or slope ('slope')
#' @param res The LDSC results dataframe, trimmed to only include data for a single GWAS and file
#' @param orderBY (optional) Only used if ct=NULL (showing all celltypes as facets). Determines whether to sort the facets by z-score or slope. Set "z" for z-score, "p" for p-value or "slope" for slope.
#' @param decreasing (optional) Should the sort order go from smallest to largest?
get.celltype.order <- function(res,orderBY,decreasing=FALSE){
    # Check arguments
    if(sum(orderBY %in% c("z","slope","p"))!=1){stop("ERROR: orderBY must be either 'z', 'p' or 'slope' if ct=NULL")}
    if(length(unique(res$gwas))>1){stop("ERROR: 'res' argument passed to get.celltype.order must only contain results from a single GWAS")}
    if(length(unique(res$file))>1){stop("ERROR: 'res' argument passed to get.celltype.order must only contain results from a single single cell dataset file")}
    if(length(unique(res$celltype))<2){stop("ERROR: 'res' argument passed to get.celltype.order must contain results from more than one celltype")}
    
    # Get order
    cts = unique(as.character(res$celltype))
    ct_slope = ct_zscore = ct_pval = rep(0,length(cts))
    names(ct_slope) = names(ct_zscore) = names(ct_pval) = cts
    for(ctI in cts){
        ct_res = res[res$celltype==ctI,]
        #ct_slope[ctI]  = lm(Coefficient_z.score~percentile,data=ct_res)$coefficients[2]
        ct_slope[ctI]  = lm(Coefficient_z.score~percentile,data=ct_res)$coefficients[2]
        ct_zscore[ctI] = ct_res[ct_res$percentile==1,"Coefficient_z.score"]
        ct_pval[ctI]   = ct_res[ct_res$percentile==1,"Enrichment_p"]
    }
    if(orderBY=="z"){ct_order=ct_zscore}
    if(orderBY=="slope"){ct_order=ct_slope}  
    if(orderBY=="p"){ct_order=ct_pval}  
    
    cto = names(sort(ct_order,decreasing=decreasing))
    ct_mapping = unique(res[,c("celltype","writtenName")])
    rownames(ct_mapping)=ct_mapping$celltype
    res$celltype = factor(as.character(res$celltype),levels=cto)
    res$writtenName = factor(as.character(res$writtenName),levels=ct_mapping[cto,]$writtenName)
    
    # If slope, then append to subres
    slope_data = data.frame(celltype=names(ct_slope),slope=ct_slope)
    res = merge(slope_data,res,by="celltype")
    
    return(res)
}