

netTS: Introduction
------------

   The netTS package is meant for relational data that takes place through time, e.g., interactions between two individuals with a time stamp. The package takes a time-aggregated network approach and uses an adjustable moving window to measure how a social network changes through time.

   The moving window approach allows a user to define the size of a time window (e.g., windowsize = 1 month) and the amount to move the window (e.g., windowshift = 1 day). This moving window then subsets the relational data within a window, creates a network, and extracts a network measure. It then shifts over in time and repeats the process. By altering the size and shift of this moving window it is then possible to measure how networks change in time.
   
   The package also provides tools for choosing window sizes, comparing observed network measures to null models, and for simulating network data to aid in statistical model construction and testing.

  
![](inst/readme_figs/diag_movingWindow.png)
  

Install netTS, load some libraries
----------------------------------

``` r
devtools::install_github("tbonne/netTS")
library(netTS)
library(lubridate)
library(ggplot2)
library(igraph)
library(reshape2)
```


Tutorials for getting started
------------------------------------------

1. [Introduction](https://tbonne.github.io/netTS/inst/Tutorials/Intro_to_netTS.html)
2. [Choosing a window size](https://tbonne.github.io/netTS/inst/Tutorials/Intro_to_netTS_WindowSizeChoice.html)
3. [Controling for sampling effort](https://tbonne.github.io/netTS/inst/Tutorials/Control_for_sampling_effort.html)
4. [Controling for nodes entering and leaving the network.](https://tbonne.github.io/netTS/inst/Tutorials/Control_for_entering_and_leaving.Rmd.html)
5. [Network similarity](https://tbonne.github.io/netTS/inst/Tutorials/Network_Similarity.html)
6. [Using network permutations](https://tbonne.github.io/netTS/inst/Tutorials/Using_permutations.html)
7. [Using network simulations](https://tbonne.github.io/netTS/inst/Tutorials/Using_simulated_event_data.html)
