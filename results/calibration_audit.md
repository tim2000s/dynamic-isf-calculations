# Does the loop's ISF match the actual effect of insulin? (calibration audit)

105 users, 93,989 fasting overnight windows. error = actual end BG − the loop's IOB-predicted end BG (which uses the ISF it chose). >0 means insulin did less than the loop expected (ISF too aggressive); <0 means it did more (ISF too weak).

## Overall: median error **2.0 mg/dL**

## By glucose at decision

| BG band | n | median error | IQR |
|---|---|---|---|
| 80-100 | 7,265 | -3.2 | -14.0–8.4 |
| 100-120 | 15,614 | -1.4 | -13.6–11.0 |
| 120-145 | 24,205 | 2.0 | -12.4–19.0 |
| 145-175 | 22,398 | 6.2 | -13.0–27.0 |
| 175-230 | 19,205 | 7.0 | -20.0–34.0 |

## By ISF formula

| formula | n | median error | IQR |
|---|---|---|---|
| no_dynisf | 90,030 | 1.4 | -15.0–21.0 |
| dynisf_sigmoid | 1,709 | 6.6 | -14.0–31.2 |
| dynisf_log | 2,250 | 21.7 | 3.0–38.0 |

## By TDD band (per-user median error)

| TDD | median error |
|---|---|
| <25 | 13.8 |
| 25-45 | 0.8 |
| 45-70 | 7.0 |
| 70+ | 5.7 |

- Loop's own autosens ratio: median 1.0; correlation of prediction error with (autosens−1) = -0.033 (if autosens were fully correcting, residual error would be ~0 and uncorrelated).

![Calibration audit](charts/inv008/fig_calibration_audit.png)

*Fasting overnight, carb-screened, full-action (4h) horizon. The loop keeps dosing over the window, so the error is the net calibration of the loop's ISF-based forecast against the realised outcome, not a pure open-loop insulin response.*