## -*- coding: utf-8 -*-

parse.model.code = function(m){
    ## Bounds can be specificed only for fixed effects and random effects' standard deviations
    bounds.fe = bounds.sd = list()
    for(par in names(m$fixed)){
        ## gamma_prior_fixed_lb == 1, gamma_prior_fixed_ub == NA -> bounds.fe$gamma == 'l'
        bounds.fe[[par]] = ''
        for(type in c('l', 'u'))
            if(!is.na(m$sdata[[sprintf('%s_prior_fixed_%sb', par, type)]]))
                bounds.fe[[par]] = paste(bounds.fe[[par]], type, sep = '')
    }
    for(par in names(m$random)){
        bounds.sd[[par]] = rep('', length(m$random[[par]]))
        for(g in 1:length(m$random[[par]])){
            for(type in c('l', 'u'))
                if(!is.na(m$sdata[[sprintf('%s_prior_random_sd_%sb_%d', par, type, g)]]))
                    bounds.sd[[par]][g] = paste(bounds.sd[[par]][g], type, sep = '')
        }
    }
    parsed = ''
    lines = readLines(paste(path.package('bhsdtr2'), '/stan_templates/ordinal.stan', sep = ''))
    l = 1
    while(l <= length(lines)){
        if(grepl('//cb', lines[l])){
            cb = l
            for(k in (l + 1):length(lines))
                if(grepl('//ce', lines[k]))break
            ce = k
            header = lines[cb]
            chunk = paste(lines[(cb + 1):(ce - 1)], collapse = '\n')
            processed = process.chunk(header, chunk, m$model, m$links, m$fixed, m$random, m$sdata, bounds.fe, bounds.sd)
            if(processed != "")
                parsed = paste(parsed, processed, sep = '\n')
            l = k + 1
        }else{
            parsed = paste(parsed, lines[l], sep = '\n')
            l = l + 1
        }
    }
    parsed
}

## Bounds possible only for fixed effects and random effects'
## std. devs.
process.chunk = function(header, chunk, model, links, fixed, random, sdata, bounds.fe, bounds.sd){
    parsed = NULL
    test = regmatches(header, regexpr('\\{.*\\}', header))
    if(length(test) == 0){
        test = parse(text = 'TRUE')
    }else{
        test = parse(text = test)
    }
    process = T
    if(grepl('fpariter', header)){
        for(par in names(fixed)){
            if(eval(test))
                parsed = c(parsed, gsub('PAR', par, chunk))
        }
    }else if(grepl('rpariter', header)){
        for(par in names(random)){
            for(g in 1:length(random[[par]]))
                if(eval(test))
                    parsed = c(parsed, gsub('G', g, gsub('PAR', par, chunk)))
        }
    }else if(eval(test)){
        parsed = chunk
    }
    paste(parsed, collapse = '\n')
}
