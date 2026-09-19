# Insulin-divisor sensitivity

AndroidAPS maps a nominal 75-minute rapid-insulin peak to divisor 55, a 55-minute ultra-rapid peak to divisor 65, and a 45-minute Lyumjev-style peak to divisor 75. The principal audit replay fixed divisor 75 because insulin type was not consistently available.

## Same-window endpoint prediction

| divisor | V1 MAE | V1 bias | V2 MAE | V2 bias |
|---:|---:|---:|---:|---:|
| 55 | 25.1 | -15.6 | 38.2 | +3.9 |
| 65 | 24.7 | -12.9 | 42.8 | +13.1 |
| 75 | 24.6 | -10.4 | 49.7 | +25.2 |

Static ISF MAE was 20.3 mg/dL on the same windows.

## Six-hour action-balance comparison

| divisor | V1 / proxy | V1 within 30% | V2 / proxy | V2 within 30% |
|---:|---:|---:|---:|---:|
| 55 | 2.25 | 17.4% | 4.37 | 10.3% |
| 65 | 2.46 | 16.4% | 5.16 | 8.8% |
| 75 | 2.67 | 15.5% | 5.93 | 7.6% |

Divisor 55 lowers both calculated ISFs and therefore makes correction stronger. The effect is modest for V1 and large for V2 near its logarithmic floor. V2 improves when divisor 55 is used, but remains worse than static ISF and V1 in the same-window test.

![Divisor sensitivity](charts/inv008/fig_divisor_sensitivity.png)
