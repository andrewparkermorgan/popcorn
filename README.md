# popcorn

A lightweight parser for various file formats produced by software used for population genetic analyses, plus a handufl of plotting utilities.  The package makes use of tools from the `tidyverse` and generally tries to respect their idioms and vernacular.

Documentation is a work in progress.

Includes parsers for these formats, among others:

* `PLINK`
	* `*.fam` (sample/family metadata): `read_fam()`
	* `*.map`/`*.bim` (marker map): `read_map()`
* `ADMIXTURE`
	* `*.Q` + `*.Q_se` (estimated admixture proportions and standard errors): `read_Q_matrix()`
* `TreeMix`
	* `read_treemix()` for loading the population tree and log-liklihood
	* `read_f3stats()` for reading $f_3$ statistics estimated by `threepop`
	* `read_f4stats()` for reading $f_4$ statistics estimated by `fourpop`

## Credits
This package borrows internal functions from at least the following packages:

* `TreeMix` (https://bitbucket.org/nygcresearch/treemix/)
* `admixturegraph` (https://github.com/mailund/admixture_graph/)
* `pophelper` (https://github.com/royfrancis/pophelper/)

## Examples	
Suppose we ran `ADIXTURE` with $K = 3$ and bootstrap standard errors turned on like so:
```
$> admixture -B test.bed 3
```

We obtain these output files:
```
$> ls test*
test.3.Q   test.3.Q_bias   test.3.Q_se   test.3.P
```

The code below reads the result and makes the typical barplot of individual admixture proportions.
```
library(popcorn)
library(ggplot2)

pops <- read_fam("test.fam")
Q <- read_Q_matrix("test.3.Q")
so <- sort_by_cluster(Q)
QQ <- tidy(Q, pops)

plot_admixture(Q, label = TRUE, sort.order = so)
```