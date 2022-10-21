extraDerivatives = function(){
  
  logit = function(x) log(x/(1-x))
  sigmoid = function(x) 1/(1+exp(-x))
  odds = function(x) x/(1-x)

  Deriv(quote(sigmoid(x)),cache.exp = F)
  Deriv(quote(logit(x)),cache.exp = F)
  Deriv(quote(odds(x)),cache.exp = F)
  
}
