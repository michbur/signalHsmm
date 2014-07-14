#viterbi

duration.viterbi <- function(probka, pipar, tpmpar, od, params){
  max.duration <- dim(params)[1]
  viterbi <- matrix(nrow=length(probka), ncol=4)
  psi <- matrix(nrow=length(probka), ncol=4)
  dura <- matrix(nrow=length(probka), ncol=4)
  for(j in 1:4){
    viterbi[1,j]=log(pipar[j])+log(od[j,probka[1]])
    psi[1,j] = 1
    dura[1,j] = 1
  }  
  for(i in 2:length(probka)){
    for(j in 1:4){
      max=-Inf
      max.i = 0
      max.dur = 0 
      for(k in 1:4){
        for(d in 1:min(max.duration, i)){
          # poprawka na wszystkie sygnaly
          #rozpisac na zmienne temp - krok po kroku
          if(i-d==0){
            if(j==1){
              previous = 1
              transition <- 1
            }
            else{
              previous = -Inf
              transition <- 1
            }
          }
          else{
            previous <- viterbi[i-d,k]
            transition <- log(tpmpar[k,j])
          }
          duration <- log(params[d,j])
          responses <- sum(log(od[j,probka[(i-d+1):i]]))
          #print(probka[(i-d+1):i])
          #print(od[j,probka[(i-d+1):i]])
          #print(previous)
          #print(transition)
          #print(duration)
          #print(responses)
          if(previous + transition + duration + responses>max){
            max = previous + transition + duration + responses
            max.i = k
            max.dur = d
          }
        }
        
      }
      viterbi[i,j]=max
      psi[i,j]=max.i
      dura[i,j]=max.dur
    }  
  }
  path = NULL
  path[length(probka)] = which.max(viterbi[length(probka),])
  i = length(probka)-1
  last = length(probka)
  while(i>1){
    if(last-i < dura[last, path[last]]){
      path[i] = path[last]
    }
    else{
      path[i] = psi[last, path[last]]
      last=i
    }
    i = i -1
  }
  path[1]=1

  return(list(path =path,
              viterbi = viterbi,
              psi = psi,
              duration = dura))
}
