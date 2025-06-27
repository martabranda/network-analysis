# Network Analysis for Time Series Data

Implementation of network models for longitudinal data.

## Methods

- **Graphical VAR Models**: Temporal and contemporaneous network estimation
- **Time Series Detrending**: Removing day/beep effects before analysis
- **Network Preprocessing**: Data centering and stationarity checks
- **Bootstrap Stability**: Assessing edge and centrality reliability

## Requirements

```r
install.packages(c("qgraph", "graphicalVAR", "mlVAR", "bootnet", 
                  "dplyr", "lm.beta"))
