# rstl
A Python port of R's stl function. Translated rather literally from the original R/Fortran source and vectorized with NumPy.

For more information see [the R manual](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/stl.html).

## Installation
```
pip install --user rstl
```

## Usage
Import the STL class and call its constructor.

```python
import numpy as np
from rstl import STL

ts = np.arange(144)
freq = 12

stl = STL(ts, freq, "periodic")

trend = stl.trend
```

## Documentation
```python
class STL(ts, freq, s_window, s_degree=0, t_window=None,
          t_degree=1, l_window=None, l_degree=None, s_jump=None,
          t_jump=None, l_jump=None, robust=False, inner=None, outer=None):
```

Note that a value of None means that the default R value will be used.
(These could not be specified in the signature because Python does not allow referencing other parameters).

Args:
* ts: The time series (numpy array).
* freq: The frequency of the time series.
* s_window: Either the character string "periodic" or the span (in lags) of the loess window for seasonal extraction, which should be odd and at least 7, according to Cleveland et al.

Optional args:
* s_degree, t_window, t_degree, l_window, l_degree, s_jump, t_jump, l_jump, robust, inner and outer. See [the R manual](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/stl.html) for their meanings and defaults.

Attributes:
* trend: The trend component of the time series (numpy array).
* seasonal: The seasonal component of the time series (numpy array).
* remainder: The remainder of the time series not explained by trend and seasonal components (numpy array).
* weights: The final robust weights (all one if fitting is not done robustly) (numpy array).

* s_window, s_degree, t_window, t_degree, l_window, l_degree, s_jump, t_jump, l_jump, inner and outer. Note that these may have been altered by the program.

## Performance
According to tests a runtime increase of factor ~3 should be expected.

## Copyright and License
Python port Copyright 2018 Eric Rauch.

Original source Copyright 2014 B.D. Ripley; Fortran code by Cleveland et al (1990) from ‘netlib’.

Licensed under the GPLv3.
