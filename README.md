# ndvi_ram

This code is for extracting phenology metrics from MODIS data using double-logistic curve fitting. First, MODIS data is downloaded/updated using the MODIStsp package. The selected options for the different time series downloaded are saved in the 4 .json files where dates, passwords and possibly paths have to be updated to download the latest files. The .json files can be opened through the Load Options button like this:

```r
library(MODIStsp)
MODIStsp()
```

The dates, the username, the password and the paths have to be updated in the .json file prior to loading the file (some functionalities do not work well, e.g. when searching for the path to the UdS NDVI server path). Some functions within the MODIStsp package could be used to directly download the data through a script, but some thing are buggy and do no seem to work well (dates, extent not respected).Thus the interactive UI has to be used for the download.

The rest is in the get_metrics.r script that will produce a .csv file with the phenology values.
