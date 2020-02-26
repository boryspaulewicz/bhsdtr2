## -*- coding: utf-8 -*-

#' Creates the SDT stan model code.
#'
#' \code{make_stan_model} function creates the stan model code
#' defining the SDT model with additional regression structure and
#' random effects if these are specified using the \code{random}
#' argument. See \code{\link{make_stan_data}} for details on how to
#' specify a random effects structure. The model is described in detail
#' in the paper (TODO reference)
#'
#' @param random an optional list specifying random effects
#'     structure. See \code{\link{make_stan_data}} for details.
#' @param gamma_link can be either 'softmax' (described in the paper),
#'     'log_distance', or 'log_ratio' (See the Readme file in the
#'     github repository)
#' @param model can be either 'sdt' (the default), 'uvsdt', 'metad',
#'     'ordinal', or 'uvordinal'
#' @return a string containing the full model definition in the stan
#'     modelling language.
#' @examples
#' data(gabor)
#' model = make_stan_model(list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1)))
#' cat(model)
#' @export
make_stan_model = function(random = NULL, gamma_link = 'softmax', model = 'sdt', delta_link = 'log'){
    check_link(gamma_link, delta_link)
    check_model(model)
    if(model %in% c('ordinal', 'uvordinal')){
        par_types = c('eta', 'gamma')
    }else{
        par_types = c('delta', 'gamma')
    }
    if(model %in% c('uvsdt', 'uvordinal'))par_types = c(par_types, 'theta')
    parsed = parse_PAR(readLines(paste(path.package('bhsdtr'), '/stan_templates/ordinal_template.stan', sep = '')), par_types)
    parsed = parse_likelihood(parsed, model)
    parsed = parse_random(parsed, random, par_types)
    links = list(delta = delta_link, gamma = gamma_link)
    for(par in names(links))
        parsed = parse_link(parsed, par, links[[par]])
    paste(parsed, collapse = '\n')
}

parse_PAR = function(lines, par_types){
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

parse_likelihood = function(lines, model){
    parsed = NULL
    for(l in lines){
        if(rmatch('//likelihood', l)){
            parsed = c(parsed, readLines(sprintf('%s/stan_templates/likelihood_%s.stan',
                                                 paste(path.package('bhsdtr')), model)))
        }else{
            parsed = c(parsed, l)
        }
    }
    parsed
}

parse_random = function(lines, random, par_types){
    parsed = NULL
    for(part in lines){
        ## If this is part of the random effects' specification ...
        if(rmatch(sprintf('//random-(%s)', paste(c('common', par_types), collapse = '|')), part)){
            ## ... and there are random effects in the model ...
            if(length(random) > 0)
                for(l in 1:length(random)){
                    ## ... then add the line read from the template with % replaced by the grouping factor number ...
                    if(rmatch('//random-common', part))
                        parsed[length(parsed)+1] = gsub('%', l, part)
                    ## ... and do the same with parts of delta/gamma
                    ## random effects' specification if delta/gamma is
                    ## associated with random effects of the given
                    ## grouping factor ...
                    for(par in par_types)
                        if(!is.null(random[[l]][[par]]) & rmatch(sprintf('//random-%s', par), part))
                            parsed[length(parsed)+1] = gsub('%', l, part)
                }
        }else{
            parsed = c(parsed, part)
        }
    }
    parsed
}

parse_link = function(lines, par, link){
    parsed = NULL
    for(l in lines){
        if(rmatch(sprintf('//%s-link', par), l)){
            parsed = c(parsed, readLines(sprintf('%s/stan_templates/%s_link_%s.stan',
                                                 paste(path.package('bhsdtr')), par, link)))
        }else{
            parsed = c(parsed, l)
        }
    }
    parsed
}

rmatch = function (pattern, vector){
  res = TRUE
  for (i in 1:length(vector)) {
    if (length(grep(pattern, vector[i])) > 0) {
      res[i] = TRUE
    }
    else {
      res[i] = FALSE
    }
  }
  res
}
