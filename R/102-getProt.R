#' @title Retrieve Protein Sequence in various Formats from Databases
#'
#' @description Retrieve Protein Sequence in various Formats from Databases(BMgetDrug)
#' 
#' @details This function retrieves protein sequence in various formats from three databases.
#' 
#' @param id A character vector, as the protein ID(s).
#' @param from The database, one of \code{'uniprot'}, \code{'kegg'}, \code{'pdb'}.
#' @param type The returned protein format, one of \code{fasta}, \code{pdb}, \code{aaseq}.
#' @param parallel An integer, the parallel parameter, indicates how many 
#'                 process the user would like to use for retrieving 
#'                 the data (using RCurl), default is \code{5}. 
#'                 For regular cases, we recommend a number less than \code{20}.
#' 
#' @return A length of \code{id} character list, each element 
#' containing the corresponding protein sequence(s) or file(s).
#' 
#' @keywords BMgetProt
#'
#' @aliases BMgetProt
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @seealso See \code{\link{BMgetDrug}} for retrieving drug molecules 
#' from five databases.
#' 
#' @export BMgetProt
#' 
#' @name getProt
#' 
#' @examples
#' \donttest{
#' # BMgetProt
#' id = c('P00750', 'P00751', 'P00752')
#' BMgetProt(id, from = 'uniprot', type = 'aaseq')}
#' 

BMgetProt = function (id,
                    from = c('uniprot', 'kegg', 'pdb'),
                    type = c('fasta', 'pdb', 'aaseq'), 
                    parallel = 5) {

    if (is.null(from)) stop('Must specify a data source')
    if (is.null(type)) stop('Must specify a data type')

    # Exclude 3 special case from total 9 possible combinations

    if (from == 'uniprot' & type == 'pdb') stop('UniProt only supports type = "fasta" or type = "aaseq"')
    if (from == 'kegg' & type == 'pdb') stop('KEGG only supports type = "fasta" or type = "aaseq"')
    if (from == 'pdb' & type == 'fasta') stop('RCSB PDB only supports type = "pdb" or type = "aaseq"')

    FromDict = c('uniprot' = 'UniProt', 'kegg' = 'KEGG', 'pdb' = 'RCSBPDB')  
    TypeDict = c('fasta' = 'FASTA', 'pdb' = 'PDB', 'aaseq' = 'Seq')

    NamePart1 = TypeDict[type]
    NamePart2 = FromDict[from]

    FunctionName = paste('get', NamePart1, 'From', NamePart2, sep = '')

    Prot = eval(parse(text = paste(FunctionName, '(', 
                                   gsub('\\"', '\'', capture.output(dput(id))), 
                                   ', ', parallel, ')', sep = '')))

    return(Prot)

}

#' @rdname getProt
#' 
#' @title Retrieve Protein Sequence (FASTA Format) from the UniProt Database
#'
#' @description Retrieve Protein Sequence (FASTA Format) from the UniProt Database(BMgetProt...UinProt)
#' 
#' @details This function retrieves protein sequences (FASTA format) 
#' from the UniProt database.
#' 
#' @keywords getProt BMgetProtFASTAUinProt UniProt
#'
#' @aliases BMgetProtFASTAUinProt
#' 
#' @export BMgetProtFASTAUinProt
#' 
#' @importFrom RCurl getURLAsynchronous
#' 
#' @references
#' UniProt. \url{https://www.uniprot.org/}
#' 
#' UniProt REST API Documentation. \url{https://www.uniprot.org/faq/28}
#' 
#' @examples
#' \donttest{
#' # BMgetProtFASTAUinProt
#' id = c('P00750', 'P00751', 'P00752')
#' BMgetProtFASTAUinProt(id)}
#' 

BMgetProtFASTAUinProt = function (id, parallel = 5) {
  
  # example id:  P00750
  # example url: http://www.uniprot.org/uniprot/P00750.fasta
  
  fastaURL = paste0('https://www.uniprot.org/uniprot/', id, '.fasta')
  
  fastaTxt = getURLAsynchronous(url = fastaURL, perform = parallel)
  
  return(fastaTxt)
  
}

#' @rdname getProt
#' 
#' @keywords getProt BMgetProtSeqUniProt UniProt
#'
#' @aliases BMgetProtSeqUniProt
#' 
#' @export BMgetProtSeqUniProt
#' 
#' @references
#' UniProt. \url{https://www.uniprot.org/}
#' 
#' UniProt REST API Documentation. \url{https://www.uniprot.org/faq/28}
#' 
#' @examples
#' \donttest{
#' # BMgetProtSeqUniProt
#' id = c('P00750', 'P00751', 'P00752')
#' BMgetProtSeqUniProt(id)}
#' 

BMgetProtSeqUniProt = function (id, parallel = 5) {
  
  # example id:  P00750
  # example url: https://www.uniprot.org/uniprot/P00750.fasta
  
  fastaTxt = BMgetProtFASTAUinProt(id, parallel)
  
  tmpfile = tempfile(pattern = paste0(id, '-'), fileext = 'fasta')
  for (i in 1:length(id)) write(fastaTxt[[i]], tmpfile[i])
  
  AASeq = lapply(tmpfile, readFASTA)
  
  unlink(tmpfile)
  
  return(AASeq)
  
}

#' @rdname getProt
#' 
#' @title Retrieve Protein Sequence (FASTA Format) from the KEGG Database
#'
#' @description Retrieve Protein Sequence (FASTA Format) from the KEGG Database(BMgetProt...KEGG)
#' 
#' @details This function retrieves protein sequences (FASTA format) 
#' from the KEGG database.
#' 
#' @keywords getProt BMgetProtFASTAKEGG KEGG
#'
#' @aliases BMgetProtFASTAKEGG
#' 
#' @export BMgetProtFASTAKEGG
#' 
#' @importFrom RCurl getURLAsynchronous
#' 
#' @examples
#' \donttest{
#' #  BMgetProtFASTAKEGG
#' id = c('hsa:10161', 'hsa:10162')
#' BMgetProtFASTAKEGG(id)}
#' 

BMgetProtFASTAKEGG = function (id, parallel = 5) {
  
  # example id : hsa:10161
  # example url: http://rest.kegg.jp/get/hsa:10161/aaseq
  
  fastaURL = paste0('http://rest.kegg.jp/get/', id, '/aaseq')
  
  fastaTxt = getURLAsynchronous(url = fastaURL, perform = parallel)
  
  return(fastaTxt)
  
}

#' @rdname getProt
#' 
#' @keywords getProt BMgetProtSeqKEGG KEGG
#' 
#' @aliases BMgetProtSeqKEGG
#' 
#' @export BMgetProtSeqKEGG
#' 
#' @examples
#' \donttest{
#' # BMgetProtSeqKEG
#' id = c('hsa:10161', 'hsa:10162')
#' BMgetProtSeqKEGG(id)}
#' 

BMgetProtSeqKEGG = function (id, parallel = 5) {
  
  # example id : hsa:10161
  # example url: http://rest.kegg.jp/get/hsa:10161/aaseq
  
  fastaTxt = BMgetProtFASTAKEGG(id, parallel)
  
  tmpfile = tempfile(pattern = paste0(id, '-'), fileext = 'fasta')
  for (i in 1:length(id)) write(fastaTxt[[i]], tmpfile[i])
  
  AASeq = lapply(tmpfile, readFASTA)
  
  unlink(tmpfile)
  
  return(AASeq)
  
}

#' @rdname getProt
#' 
#' @title Retrieve Protein Sequence (PDB Format) from RCSB PDB
#' 
#' @description Retrieve Protein Sequence (PDB Format) from RCSB PDB(BMgetProt...RCSBPDB)
#' 
#' @details This function retrieves protein sequences (PDB format) from RCSB PDB.
#' 
#' @keywords getProt BMgetProtPDBRCSBPDB PDB
#'
#' @aliases BMgetProtPDBRCSBPDB
#'
#' @export BMgetProtPDBRCSBPDB
#' 
#' @importFrom RCurl getURLAsynchronous
#' 
#' @examples
#' \donttest{
#' # BMgetProtPDBRCSBPDB
#' id = c('4HHB', '4FF9')
#' BMgetProtPDBRCSBPDB(id)}
#' 

BMgetProtPDBRCSBPDB = function (id, parallel = 5) {
  
  # example id : 4HHB
  # example url: http://www.rcsb.org/pdb/files/4HHB.pdb
  
  pdbURL = paste0('http://files.rcsb.org/view/', id, '.pdb')
  
  pdbTxt = getURLAsynchronous(url = pdbURL, perform = parallel)
  
  return(pdbTxt)
  
}

#' @rdname getProt
#' 
#' @keywords getProt BMgetProtSeqRCSBPDB PDB
#' 
#' @aliases BMgetProtSeqRCSBPDB
#' 
#' @export BMgetProtSeqRCSBPDB
#' 
#' @importFrom RCurl getURLAsynchronous
#' 
#' @examples
#' \donttest{
#' # BMgetProtSeqRCSBPDB
#' id = c('4HHB', '4FF9')
#' BMgetProtSeqRCSBPDB(id)}
#' 

BMgetProtSeqRCSBPDB = function (id, parallel = 5) {
  
  # example id : 4HHB
  # example url: http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=4HHB
  
  fastaURL = paste0('http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=', id)
  
  fastaTxt = getURLAsynchronous(url = fastaURL, perform = parallel)
  
  tmpfile = tempfile(pattern = paste0(id, '-'), fileext = 'fasta')
  for (i in 1:length(id)) write(fastaTxt[[i]], tmpfile[i])
  
  AASeq = lapply(tmpfile, readFASTA)
  
  unlink(tmpfile)
  
  return(AASeq)
  
}
