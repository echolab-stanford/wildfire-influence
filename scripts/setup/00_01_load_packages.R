library(lubridate)

library(sf) # Linking to GEOS 3.11.2, GDAL 3.6.4, PROJ 9.2.0
library(tigris)

library(magrittr)
library(tidyverse)
library(feather)

import::from(segmented, .except = "select")
library(fixest)
library(robslopes)
library(strucchange)

library(gtable)
library(cowplot)
library(geofacet)
library(ggforce)
library(ggridges)
library(MetBrewer)
library(usmap)
library(scales)

library(foreach)
library(parallel)
library(doParallel)
library(future.apply)
library(future.callr)
library(doRNG)
library(doFuture)
library(tictoc)
library(callr)
