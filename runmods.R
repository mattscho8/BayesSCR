#### Sequentially run the models
## 1. DA (x2: Stan and NIMBLE)
## 2. SCDL (x4: 2xStan, 2xNIMBLE)
## 3. ODL (x1: Stan)

#### Obtain timings and posterior summaries for each model
## Timing ignores compilation time

library(here)

source(here("DA","DA_Stan.R"))
source(here("DA","DA_NIMBLE.R"))
source(here("SCDL","SCDL_Stan_BPP.R"))
source(here("SCDL","SCDL_Stan_PPP.R"))
source(here("SCDL","SCDL_NIMBLE_BPP.R"))
source(here("SCDL","SCDL_NIMBLE_PPP.R"))
source(here("ODL","ODL_Stan_PPP.R"))
