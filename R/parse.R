## -*- coding: utf-8 -*-

stan.file = function(fname){
    sprintf('%s/stan_templates/%s', path.package('bhsdtr2'), fname)
}

## Parsing functions (model code from template)
parse.PAR = function(lines, par_types){
    parsed = NULL
    for(l in lines){
        if(rmatch('PAR', l)){
            for(par_type in par_types)
                parsed[length(parsed) + 1] = gsub('PAR', par_type, l)
        }else{
            parsed = c(parsed, l)
        }
    }
    parsed
}

parse.likelihood = function(lines, model){
    parsed = NULL
    for(l in lines){
        if(rmatch('//likelihood', l)){
            parsed = c(parsed, readLines(stan.file(sprintf('likelihood_%s.stan', model))))
        }else{
            parsed = c(parsed, l)
        }
    }
    parsed
}

parse.link = function(lines, par, link){
    parsed = NULL
    for(l in lines){
        if(rmatch(sprintf('//%s-link', par), l)){
            parsed = c(parsed, readLines(stan.file(sprintf('%s_link_%s.stan', par, link))))
        }else{
            parsed = c(parsed, l)
        }
    }
    parsed
}

parse.random = function(lines, random, par_types){
    parsed = NULL
    for(part in lines){## If this is part of the random effects' specification ...
        if(rmatch(sprintf('//random-(%s)', paste(par_types, collapse = '|')), part)){
            ## ... and there are random effects in the model ...
            if(length(random) > 0)
                for(par in names(random)){
                    for(g in 1:length(random[[par]])){
                        if(!is.null(random[[par]][[g]]) & rmatch(sprintf('//random-%s', par), part))
                            parsed[length(parsed) + 1] = gsub('%', g, part)
                    }
                }
        }else{
            ## If this is not a part of random effects structure then do not change it
            parsed = c(parsed, part)
        }
        ## if there are no random effects in the model and the
        ## template part defines part of the random effects structure then ignore this part
    }
    parsed
}
