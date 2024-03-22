# EVE scores

This directory contains EVE scores we computed for the [experimental datasets](../dms_data).
The scores are stored as hdf5 files and can be read into dataframes using the `pandas` library.

## Example usage

```python
import pandas as pd

fn = "data/eve/pab1.h5"
df = pd.read_hdf(fn)

print(df.head())
```

This would print the dataframe:

|   | variant | eve        |
|---|---------|------------|
| 0 | G0A     | -12.307312 |
| 1 | G0C     | -12.496277 |
| 2 | G0D     | -12.437378 |
| 3 | G0E     | -12.308533 |
| 4 | G0N     | -12.576096 |

Note the variants are 0-indexed.