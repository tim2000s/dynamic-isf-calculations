# Does the loop's ISF match the actual effect of insulin? (calibration audit)

92 users, 86,044 fasting overnight windows. error = actual end BG − the loop's IOB-predicted end BG (which uses the ISF it chose). >0 means insulin did less than the loop expected (ISF too aggressive); <0 means it did more (ISF too weak).

## Overall: median error **1.0 mg/dL**

## By glucose at decision

| BG band | n | median error | IQR |
|---|---|---|---|
| 80-100 | 7,026 | -3.2 | -14.0–8.6 |
| 100-120 | 14,792 | -2.0 | -14.0–10.2 |
| 120-145 | 21,860 | 0.4 | -13.0–16.6 |
| 145-175 | 19,922 | 5.0 | -13.8–25.2 |
| 175-230 | 17,478 | 5.0 | -21.8–32.0 |

## By ISF formula

| formula | n | median error | IQR |
|---|---|---|---|
| no_dynisf | 86,044 | 1.0 | -15.4–20.0 |
| dynisf_sigmoid | – | – | – |
| dynisf_log | – | – | – |

## By TDD band (per-user median error)

| TDD | median error |
|---|---|
| <25 | 10.4 |
| 25-45 | 0.5 |
| 45-70 | 7.4 |
| 70+ | 5.7 |

- Loop's own autosens ratio: median 1.0; correlation of prediction error with (autosens−1) = -0.05 (if autosens were fully correcting, residual error would be ~0 and uncorrelated).

![Calibration audit](charts/inv008/fig_calibration_audit.png)

*Fasting overnight, carb-screened, full-action (4h) horizon. The loop keeps dosing over the window, so the error is the net calibration of the loop's ISF-based forecast against the realised outcome, not a pure open-loop insulin response.*