library(readr)
library(stringr)

readRefTable=function(filename = "reftable.bin", header = "header.txt", N = 0)
{
##################################################################
    ## Headers parsing
    header = read_file(header)          # whole file in one string
    
### First part, scenarios
    ##  extraction of scenarios descriptions
    rescen = "\\bscenario\\s+(\\d+)\\s+.*\\n((?:(?!(?:scenario|\\n)).*\\n)+)"
    scenprematch = str_extract_all(header,rescen)[[1]]
    nscenh = length(scenprematch)       # Number of scenarios
    scenmatch = str_match_all(scenprematch,rescen)
    ## scendesc is the vector of scenarios descriptions
    scendesc <- vector(mode="character", length=nscenh)
    # loop on each scenario to store their descriptions
    for(i in 1:nscenh)
    {
        scendesc[as.integer(scenmatch[[i]][2])] = scenmatch[[i]][3]
    }

### Second part, hist params
    ## extraction of parameters names (and the corresponding law too,
    ## in order to check if it is a dirac base, aka constant, as we
    ## ignore it in the reftable)
    reparamtot = "historical parameters priors \\((\\d+)\\D"
    ## nparamtoth is the total number of hist params in header
    nparamtoth = as.integer(str_match(str_extract(header,reparamtot)[[1]],reparamtot)[2])
    reparamlist = paste0("\\bhistorical parameters priors.*\\n((?:\\w+\\W[^\\n]*\\n){",nparamtoth,"})")
    paramlistmatch = str_match(str_extract_all(header,reparamlist),reparamlist)[2]
    ## regexp to extract : param name, mini and maximum of the law
    reparam = "(\\w+)\\W+\\w\\W+\\w\\w\\[([^,\\]]+),([^,\\]]+)[,\\]][^\\n]*\\n"
    paramsh = str_match_all(str_extract_all(paramlistmatch,reparam)[[1]],reparam)
    paramsdesc = list()
    reali = 1;
    ## loop on all hist parameters to "filter" all constants
    ## parameters
    for(i in 1:nparamtoth)
    {
        mini = as.numeric(paramsh[[i]][3])
        maxi = as.numeric(paramsh[[i]][4])
        if(maxi != 0.0)
            if ((maxi-mini)/maxi > 0.000001)
                {
                    paramsdesc[paramsh[[i]][2]] = reali
                    reali = reali + 1
                }
    }
    realparamtot = reali - 1 # this is the real number of non constant
                                        # historical parameters

### Third part, extract parameters for each scenario
    # parambyscenh stores the indices of used parameters for each
    # scenario
    parambyscenh <- vector(mode="numeric",length=nscenh)
    # Loop for all scenarios
    for(i in 1:nscenh)
    {
        templist = list()
        # list of terms used in scenario descriptions
        listterms = strsplit(scendesc[i],"\\W")[[1]]
        m = 1
        for(j in 1:length(listterms)) 
        {
            if (!is.null(paramsdesc[listterms[j]][[1]]))
            {
                templist[m] = paramsdesc[listterms[j]][[1]]
                m = m + 1
            }
        }
        parambyscenh[i] = list((unique(unlist(templist))))
    }

    # Get all names
    restatsname = "\\n\\nscenario\\s+.*"
    allcolspre = tail(strsplit(str_extract(header,restatsname),"\\s+")[[1]],-2)

##################################################################
    ## Reftable feed
    # Stream from reftable file
    to.read = file(filename,"rb")
    # number of records
    realnrec = readBin(to.read, integer(), endian = "little")
    nrec = if (N > 0) { min(N,realnrec) } else { realnrec }
    # number of scenarios
    nscen = readBin(to.read, integer(), endian = "little")
    # number of records for each scenario
    nrecscen = readBin(to.read, integer(), n = nscen, endian = "little")
    # number of used parameters (non constant)
    nparam = readBin(to.read, integer(), n = nscen, endian = "little")
    # number of stats
    nstat = readBin(to.read, integer(), endian = "little")
    # Get all params
    paramsname = head(allcolspre,-nstat)
    # Number of mutation parameters
    nmutparams = length(paramsname) - realparamtot
    # Stats matrix
    stats = matrix(nrow = nrec, ncol = nstat)
    colnames(stats) <- tail(allcolspre,nstat)
    # Param matrix
    params = matrix(nrow = nrec, ncol = realparamtot + nmutparams)
    colnames(params) <- paramsname
    # Scenarios vector
    scenarios = vector(mode="numeric", length = nrec)
    # Parsing of records
    for(i in 1:nrec)
    {
        scen = readBin(to.read, integer(), endian = "little")
        scenarios[i] = scen
        lparams = readBin(to.read, numeric(), n = nparam[scen], size = 4, endian = "little")
        for(j in 1:length(parambyscenh[[scen]]))
        {
            params[i,parambyscenh[[scen]][j]] = lparams[j]
        }
        # Accounting for mutational parameters
        if (nmutparams > 0)
        {
            for (jm in 1:nmutparams)
            {
                params[i,realparamtot+jm] = lparams[nparam[scen]-nmutparams+jm]
            }
        }
        lstats = readBin(to.read, numeric(), n = nstat, size = 4, endian = "little")
        for(j in 1:nstat)
        {
            stats[i,j] = lstats[j]
        }
    }
    close(to.read)
    list(nrec = nrec, nscen = nscen, nrecscen = nrecscen, nparam = nparam, scenarios = as.factor(scenarios), params = params, stats = stats)
}
