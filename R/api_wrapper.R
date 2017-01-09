db.list <- c("WikiPathways_2016", "KEGG_2016", "Biocarta_2016", "GO_Biological_Process_2015")


#' Perform functional enrichment on a pair of "up"- and "down"- regulated genes sets,
#' such as those generated in a differential expression analysis.
#'
#' This function wraps enrichGeneList to call Enrichr to perform functional enrichment tests on
#' the provided gene lists, for the databases specified in the databases argument.
#' Separate analyses are run on upregulated and downregulated genes, the returned results are rbind'ed.
#' Databases are specified as seen in the web interface, with underscores for spaces
#' (e.g. "WikiPathways_2016", "KEGG_2016", "GO_Biological_Process"). See \url{http://amp.pharm.mssm.edu/Enrichr/#stats}
#' for more databases.
#' @param up.genes a list of up-regulated gene symbols
#' @param dn.genes a list of down-regulated gene symbols
#' @param databases a character vector of Enrichr-fronted databases. Default: KEGG_2016
#' @param fdr.cutoff An FDR (adjusted p-value) threshold by which to limit the list of enriched pathways.
#' Default = 0.1
#' @keywords functional enrichment Enrichr
#' @export
enrichFullGeneList <- function(up.genes, dn.genes, databases = "KEGG_2016", fdr.cutoff = 0.1) {
    up.gene.res <- enrichGeneList(up.genes, databases, fdr.cutoff)
    dn.gene.res <- enrichGeneList(dn.genes, databases, fdr.cutoff)
    if (nrow(up.gene.res) > 0) {
      up.gene.res$direction <- "UP"
    }
    if (nrow(dn.gene.res) > 0) {
      dn.gene.res$direction <- "DN"
    }
    if (nrow(up.gene.res) == 0 & nrow(dn.gene.res) == 0) {
      return(data.frame(databases = databases, category = "Nothing significant", stringsAsFactors = FALSE))
    } else {
      return(rbind(up.gene.res, dn.gene.res))
    }
}


#' Perform functional enrichment on a set of genes.
#'
#' This function interacts with Enrichr's REST API in order to perform functional enrichment of a single
#' set of genes, for a set of specified databases which are already fronted by Enrichr.
#' Databases are specified as seen in the web interface, with underscores for spaces
#' (e.g. "WikiPathways_2016", "KEGG_2016", "GO_Biological_Process"). See \url{http://amp.pharm.mssm.edu/Enrichr/#stats}
#' for more databases.
#'
#' @param gene.list a character vector of gene symbols. Required
#' @param databases a character vector of Enrichr-fronted databases. Default: KEGG_2016
#' @param fdr.cutoff An FDR (adjusted p-value) threshold by which to limit the list of enriched pathways.
#' Default = 0.1
#' @keywords functional enrichment Enrichr
#' @export
#' @examples
#' \donotrun{
#' # Analysis of genes associated with asperger syndrome
#' genes <- c("DISC1", "ASPG4", "ASPG1", "ASPG2", "SLC6A4", "ASPG3", "FRAXE", "FRAXA", "FHIT", "NTM", "SLTM", "RASGRP4", "NOS2", "NOS1", "SHANK3", "DISC2", "TSNAX", "OXTR", "ARSD")
#' res <- enrichGeneList(gene.list = genes, databases = c("KEGG_2016", "WikiPathways_2016"), fdr.cutoff = 0.1)
#' DT::datatable(res)
#' }
##
enrichGeneList <- function(gene.list, databases = "KEGG_2016", fdr.cutoff = 0.1) {
    ######Step 1: Post gene list to EnrichR
    req.body <- list(list=paste(gene.list, collapse="\n"))
    post.req <- httr::POST("http://amp.pharm.mssm.edu/Enrichr/enrich", encode="multipart", body=I(req.body))

    #TODO: Real error handling..
    if (!grepl("success", httr::http_status(post.req)$category, ignore.case=T)) stop("Posting gene list to EnrichR failed")

    ######Step 2: Get results from posted gene list
    database.enrichments <- list()
    for (idx in 1:length(databases)) {
        database <- databases[idx]
        get.req <- httr::GET(paste("http://amp.pharm.mssm.edu/Enrichr/enrich?backgroundType=", database, sep=""))
        if (!grepl("success", httr::http_status(get.req)$category, ignore.case=T)) stop("Retrieving results from EnrichR failed")

        response.content <- mungeResponseContent(httr::content(get.req)[[database]])

        if (length(response.content) > 1) {
            database.res <- data.table::rbindlist(response.content)
            database.res[, 1] <- rep(database, nrow(database.res))
            database.enrichments[[idx]] <- database.res[, paste("V", c(1, 2, 3, 7, 6), sep=""), with=F]
        }
    }

    query.results <- as.data.frame(data.table::rbindlist(database.enrichments))
    colnames(query.results) <- c("database", "category", "pval", "qval", "genes")

    if (!is.null(fdr.cutoff)) {
        query.results <- query.results[query.results$qval < fdr.cutoff, ]
    }

    return(query.results)
}



#' Munge the Enrichr API response so it'll squeeze neatly (if untidily) into a dataframe.
#'
#' The response from the Enrichr API is a list of lists, where each nested list item represents an enriched
#' category. The 6th item of each category (i.e. response.content[[category.idx]][[6]]) corresponds to the
#' genes that overlapped with the gene set behind that category. This function bascically collapses that list of
#' genes into a single string.
#'
#' I'm sorry you ever had to look at this.
#'
#' @param response.content result of calling httr::content on the GET request to the Enrichr API, after submitting a list for enrichment.
#'
mungeResponseContent <- function(response.content) {
    munged.content <- response.content
    if (length(response.content) == 0) return(NA)

    for (idx in 1:length(response.content)) {
        munged.content[[idx]][[6]] <- paste(munged.content[[idx]][[6]], collapse=",")
    }

    return(munged.content)
}




