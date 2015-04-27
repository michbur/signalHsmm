![signalHsmm](https://github.com/michbur/signalHsmm/blob/master/inst/logo.png)

signal.hsmm
=========================
Predict Presence of Signal Peptides
-------------------------

signalHsmmpredicts presence of signal peptides in eukaryotic proteins using
hidden semi-Markov models. The implemented algorithm can be accessed as a web-based service http://michbur.shinyapps.io/signalhsmm/ 

Local instance of signal.hsmm
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
library(signal.hsmm)
gui_signal.hsmm()
```
All signalHsmm functionalities can be also invoked in batch mode, for example:

```R
run_signal.hsmm(benchmark_dat[1:10])
```
