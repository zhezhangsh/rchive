# rchive
An comprehensive archive of genomic data structured and stored as R data objects




#### To install and load rchive library within your R console
It's highly recommended to store all Rchive data collection in a common folder <RCHIVE_HOME> and set a permanent environment variable as the path to this folder.

```
# Create a permanent environment variable, replace </path/to/rchive> with the folder of your own choice
# If the RCHIVE_HOME variable is not set, all Rchive data collection will be saved to the root directory by default
write('Sys.setenv(RCHIVE_HOME = "/path/to/rchive")', file = "~/.Rprofile", append = TRUE);

library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);
```

#### To load the functions without installing the rchive package
```
# Load R functions
library(devtools);
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");
```
