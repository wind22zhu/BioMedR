#' @title Retrieve Drug Molecules in MOL and SMILES Format from Databases
#'
#' @description Retrieve Drug Molecules in MOL and SMILES Format from Databases(BMgetDrug)
#' 
#' @details This function retrieves drug molecules in MOL and SMILES format from five databases.
#' 
#' @param id A character vector, as the drug ID(s).
#' @param from The database, one of \code{'pubchem'}, \code{'chembl'}, \code{'cas'}, 
#'             \code{'kegg'}, \code{'drugbank'}.
#' @param type The returned molecule format, \code{mol} or \code{smile}.
#' @param parallel An integer, the parallel parameter, indicates how many 
#'                 process the user would like to use for retrieving 
#'                 the data (using RCurl), default is \code{5}. 
#'                 For regular cases, we recommend a number less than \code{20}.
#' 
#' @return A length of \code{id} character vector, 
#' each element containing the corresponding drug molecule.
#' 
#' @keywords BMgetDrug
#'
#' @aliases BMgetDrug
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @seealso See \code{\link{BMgetProt}} for retrieving protein sequences 
#' from three databases.
#' 
#' @export BMgetDrug
#' 
#' @name getDrug
#' 
#' @examples
#' \donttest{
#' # BMgetDrug
#' id = c('DB00859', 'DB00860')
#' BMgetDrug(id, 'drugbank', 'smile')}
#' 

BMgetDrug = function (id, 
                    from = c('pubchem', 'chembl', 'cas', 'kegg', 'drugbank'), 
                    type = c('mol', 'smile'), 
                    parallel = 5) {

    if (is.null(from)) stop('Must specify a data source')
    if (is.null(type)) stop('Must specify a data type')

    # Exclude 1 special case from total 10 possible combinations

    if (from == 'cas' & type == 'smile') stop('CAS only supports type = "mol" (InChI)')

    FromDict = c('pubchem' = 'PubChem', 'chembl' = 'ChEMBL', 
                 'cas' = 'CAS', 'kegg' = 'KEGG', 'drugbank' = 'DrugBank')

    TypeDict = c('mol' = 'Mol', 'smile' = 'Smi')

    NamePart1 = TypeDict[type]
    NamePart2 = FromDict[from]

    FunctionName = paste('get', NamePart1, 'From', NamePart2, sep = '')

    Drug = eval(parse(text = paste(FunctionName, '(', 
                                   gsub('\\"', '\'', capture.output(dput(id))), 
                                   ', ', parallel, ')', sep = '')))

    return(Drug)

}


#' @rdname getDrug
#' 
#' @title Retrieve Drug Molecules in MOL and Smi Format from the PubChem Database
#'
#' @description Retrieve Drug Molecules in MOL and Smi Format from the PubChem 
#' Database(BMgetDrug...PubChem)
#' 
#' @details This function retrieves drug molecules in MOL format from the PubChem database.
#' 
#' @keywords getDrug BMgetDrugMolPubChem PubChem
#'
#' @aliases BMgetDrugMolPubChem
#' 
#' @export BMgetDrugMolPubChem
#' 
#' @importFrom RCurl getURLAsynchronous
#' 
#' @examples
#' \donttest{
#' # BMgetDrugMolPubChem
#' id = c('7847562', '7847563')  # Penicillamine
#' BMgetDrugMolPubChem(id)}
#' 

BMgetDrugMolPubChem = function (id, parallel = 5) {
  
  # example id : 7847562 (Penicillamine)
  # example url: http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=7847562&disopt=DisplaySDF
  
  SdfURL = paste0('http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=', 
                  id, '&disopt=DisplaySDF')
  
  SdfTxt = getURLAsynchronous(url = SdfURL, perform = parallel)
  
  return(SdfTxt)
  
}

#' @rdname getDrug
#' 
#' @keywords getDrug BMgetDrugSmiPubChem PubChem
#'
#' @aliases BMgetDrugSmiPubChem
#' 
#' @export BMgetDrugSmiPubChem
#' 
#' @importFrom rcdk load.molecules get.smiles
#' 
#' @examples
#' \donttest{
#' # BMgetDrugsmiPubChem
#' id = c('7847562', '7847563')  # Penicillamine
#' BMgetDrugSmiPubChem(id)}
#' 

BMgetDrugSmiPubChem = function (id, parallel = 5) {
  
  # example id : 7847562 (Penicillamine)
  # example url: http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=7847562&disopt=DisplaySDF
  
  SdfTxt = BMgetDrugMolPubChem(id, parallel)
  
  # rcdk::load.molecules() only loads files on the disk
  # so we have to do this
  tmpfile = tempfile(pattern = paste0(id, '-'), fileext = 'sdf')
  for (i in 1:length(id)) write(SdfTxt[[i]], tmpfile[i])
  Mol = load.molecules(tmpfile)
  Smi = sapply(Mol, get.smiles)
  unlink(tmpfile)
  
  return(Smi)
  
}

#' @rdname getDrug
#' 
#' @title Retrieve Drug Molecules in MOL and Smi Format from the ChEMBL Database
#'
#' @description Retrieve Drug Molecules in MOL and Smi Format from the ChEMBL 
#' Database(BMgetDrug...ChEMBL)
#' 
#' @details This function retrieves drug molecules in MOL format from the ChEMBL database.
#' 
#' @keywords getDrug BMgetDrugMolChEMBL ChEMBL
#'
#' @aliases BMgetDrugMolChEMBL
#' 
#' @export BMgetDrugMolChEMBL
#' 
#' @importFrom RCurl getURLAsynchronous
#' 
#' @examples
#' \donttest{
#' # BMgetDrugMolChEMBL
#' id = 'CHEMBL1430'  # Penicillamine
#' BMgetDrugMolChEMBL(id)}
#' 

BMgetDrugMolChEMBL = function (id, parallel = 5) {
  
  # example id : CHEMBL1430 (Penicillamine)
  # example url: https://www.ebi.ac.uk/chembldb/compound/inspect/CHEMBL1430
  # then we get: https://www.ebi.ac.uk/chembldb/download_helper/getmol/369179
  
  MolPageURL = paste0('https://www.ebi.ac.uk/chembldb/compound/inspect/', id)
  MolPageTxt = getURLAsynchronous(url = MolPageURL, perform = parallel)
  
  n = length(id)
  tmp1 = rep(NA, n)
  tmp2 = rep(NA, n)
  
  for (i in 1:n) {
    tmp1[i] = strsplit(MolPageTxt, 
                       "<a href='/chembldb/download_helper/getmol/")[[1]][2]
  }
  
  for (i in 1:n) {
    tmp2[i] = strsplit(tmp1[i], "'>Download MolFile")[[1]][1]
  }
  
  MolURL = paste0('https://www.ebi.ac.uk/chembldb/download_helper/getmol/', tmp2)
  MolTxt = getURLAsynchronous(url = MolURL, perform = parallel)
  
  return(MolTxt)
  
}

#' @rdname getDrug
#' 
#' @keywords getDrug BMgetDrugSmiChEMBL ChEMBL
#'
#' @aliases BMgetDrugSmiChEMBL
#' 
#' @export BMgetDrugSmiChEMBL
#' 
#' @importFrom RCurl getURLAsynchronous
#' @importFrom rjson fromJSON
#' 
#' @examples
#' \donttest{
#' # BMgetDrugSmiChEMBL
#' id = 'CHEMBL1430'  # Penicillamine
#' BMgetDrugSmiChEMBL(id)}
#' 

BMgetDrugSmiChEMBL = function (id, parallel = 5) {
  
  # example id : CHEMBL1430 (Penicillamine)
  # example url: https://www.ebi.ac.uk/chemblws/compounds/CHEMBL1430.json
  
  MolURL = paste0('https://www.ebi.ac.uk/chemblws/compounds/', id, '.json')
  
  MolTxt = getURLAsynchronous(url = MolURL, perform = parallel)
  
  SmiTxt = lapply(MolTxt, fromJSON)
  
  Smi = sapply(unlist(SmiTxt, recursive = FALSE), `[[`, 'smiles')
  
  names(Smi) = NULL
  
  return(Smi)
  
}

#' @rdname getDrug
#' @title Retrieve Drug Molecules in InChI Format from the CAS Database
#'
#' @description Retrieve Drug Molecules in InChI Format from the CAS Database(BMDrugMolCAS)
#' 
#' @details This function retrieves drug molecules in InChI format from the CAS database.
#' CAS database only provides InChI data, so here we return the molecule 
#' in InChI format, users could convert them to SMILES format using 
#' Open Babel (\url{http://openbabel.org/}) or other third-party tools.
#' 
#' @keywords getDrug BMDrugMolCAS CAS
#' 
#' @export BMDrugMolCAS
#' 
#' @importFrom RCurl getURLAsynchronous
#' 
#' @examples
#' \donttest{
#' # BMDrugMolCAS
#' id = '52-67-5'  # Penicillamine
#' BMDrugMolCAS(id)}
#' 

BMDrugMolCAS = function (id, parallel = 5) {
  
  # CAS only provides InChI, so here we return InChI 
  # users could convert to SMILES using Open Babel
  
  # example id : 52-67-5 (Penicillamine)
  # example url: http://www.chemnet.com/cas/supplier.cgi?terms=52-67-5&exact=dict
  
  InChIURL = paste0('http://www.chemnet.com/cas/supplier.cgi?terms=', id, 
                    '&exact=dict')
  
  InChITxt = getURLAsynchronous(url = InChIURL, perform = parallel)
  
  n = length(id)
  tmp1 = rep(NA, n)
  tmp2 = rep(NA, n)
  for (i in 1:n) tmp1[i] = strsplit(InChITxt[[i]], 'InChI=')[[1]][2]
  for (i in 1:n) tmp2[i] = strsplit(tmp1[i], '</td>')[[1]][1]
  
  InChI = paste0('InChI=', tmp2)
  
  return(InChI)
  
}

#' @rdname getDrug
#' 
#' @title Retrieve Drug Molecules in MOL and Smi Format from the KEGG Database
#'
#' @description Retrieve Drug Molecules in MOL and Smi Format from the KEGG 
#' Database(BMgetDrug...KEGG)
#' 
#' @details This function retrieves drug molecules in MOL format from the KEGG database.
#' 
#' @keywords getDrug BMgetDrugMolKEGG KEGG
#' 
#' @export BMgetDrugMolKEGG
#' 
#' @importFrom RCurl getURLAsynchronous
#' 
#' @examples
#' \donttest{
#' # BMgetDrugMolKEGG
#' id = 'D00496'  # Penicillamine
#' BMgetDrugMolKEGG(id)}
#' 

BMgetDrugMolKEGG = function (id, parallel = 5) {
  
  # example id     : D00496 (Penicillamine)
  # example url    : http://rest.kegg.jp/get/D00496/mol
  # KEGG API Intro : http://www.kegg.jp/kegg/rest/keggapi.html
  
  MolURL = paste0('http://rest.kegg.jp/get/', id, '/mol')
  
  MolTxt = getURLAsynchronous(url = MolURL, perform = parallel)
  
  return(MolTxt)
  
}

#' @rdname getDrug
#' 
#' @keywords BMgetDrug BMgetDrugSmiKEGG KEGG
#'
#' @aliases BMgetDrugSmiKEGG
#' 
#' @export BMgetDrugSmiKEGG
#' 
#' @importFrom rcdk load.molecules get.smiles
#' 
#' @examples
#' \donttest{
#' # BMgetDrugSmiKEGG
#' id = 'D00496'  # Penicillamine
#' BMgetDrugSmiKEGG(id)}
#' 

BMgetDrugSmiKEGG = function (id, parallel = 5) {
  
  # example id     : D00496 (Penicillamine)
  # example url    : http://rest.kegg.jp/get/D00496/mol
  # KEGG API Intro : http://www.kegg.jp/kegg/rest/keggapi.html
  
  MolTxt = BMgetDrugMolKEGG(id, parallel)
  
  # rcdk::load.molecules() only loads files on the disk
  # so we have to do this
  tmpfile = tempfile(pattern = paste0(id, '-'), fileext = 'mol')
  for (i in 1:length(id)) write(MolTxt[[i]], tmpfile[i])
  
  Mol = load.molecules(tmpfile)
  Smi = sapply(Mol, get.smiles)
  
  unlink(tmpfile)
  
  return(Smi)
  
}

#' @rdname getDrug
#' 
#' @title Retrieve Drug Molecules in MOL and Smi Format from the DrugBank Database
#'
#' @description Retrieve Drug Molecules in MOL and Smi Format from the DrugBank 
#' Database(BMgetDrug...DrugBank)
#' 
#' @details This function retrieves drug molecules in MOL format from the DrugBank database.
#' 
#' @keywords getDrug BMgetDrugMolDrugBank DrugBank
#'
#' @aliases BMgetDrugMolDrugBank
#' 
#' @export BMgetDrugMolDrugBank
#' 
#' @importFrom RCurl getURLAsynchronous
#' 
#' @examples
#' \donttest{
#' # BMgetDrugMolDrugBank
#' id = 'DB00859'  # Penicillamine
#' BMgetDrugMolDrugBank(id)}
#' 

BMgetDrugMolDrugBank = function (id, parallel = 5) {
  
  # example id : DB00859 (Penicillamine)
  # example url: http://www.drugbank.ca/drugs/DB00859.sdf
  
  SdfURL = paste0('http://www.drugbank.ca/drugs/', id, '.sdf')
  
  SdfTxt = getURLAsynchronous(url = SdfURL, perform = parallel)
  
  return(SdfTxt)
  
}

#' @rdname getDrug
#' 
#' @keywords getDrug BMgetDrugSmiDrugBank DrugBank
#'
#' @aliases BMgetDrugSmiDrugBank
#' 
#' @export BMgetDrugSmiDrugBank
#' 
#' @importFrom rcdk load.molecules get.smiles
#' 
#' @examples
#' \donttest{
#' # BMgetDrugSmiDrugBank
#' id = 'DB00859'  # Penicillamine
#' BMgetDrugSmiDrugBank(id)}
#' 

BMgetDrugSmiDrugBank = function (id, parallel = 5) {
  
  # example id : DB00859 (Penicillamine)
  # example url: http://www.drugbank.ca/drugs/DB00859.sdf
  
  SdfTxt = BMgetDrugMolDrugBank(id, parallel)
  
  # rcdk::load.molecules() only loads files on the disk
  # so we have to do this
  tmpfile = tempfile(pattern = paste0(id, '-'), fileext = 'sdf')
  for (i in 1:length(id)) write(SdfTxt[[i]], tmpfile[i])
  
  Mol = load.molecules(tmpfile)
  Smi = sapply(Mol, get.smiles)
  
  unlink(tmpfile)
  
  return(Smi)
  
}