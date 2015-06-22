<img src="https://github.com/michbur/signal.hsmm/blob/master/inst/logo.png" alt="signalHsmm" style="height: 200px;"/>

Predict Presence of Signal Peptides
-------------------------

signalHsmm predicts presence of signal peptides in eukaryotic proteins using hidden semi-Markov models. The implemented algorithm can be accessed as a web-based service www.smorfland.uni.wroc.pl/signalHsmm 

Local instance of signalHsmm
------------------------
signalHsmm can be also used locally as the R package. It can be installed from CRAN using:

```R
install.packages("signalHsmm")
```

You can install the latest development version of the code using the `devtools` R package.

```R
# Install devtools, if you haven't already.
install.packages("devtools")

library(devtools)
install_github("michbur/signalHsmm")
```

After installation GUI can be accessed locally:

```R
library(signalHsmm)
gui_signalHsmm()
```
All signalHsmm functionalities can be also invoked in batch mode, for example:

```R
run_signalHsmm(benchmark_dat[1:10])
```
