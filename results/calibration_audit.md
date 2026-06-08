# Does the loop's ISF match the actual effect of insulin? (calibration audit)

102 users, 88,993 fasting overnight windows. error = actual end BG − the loop's IOB-predicted end BG (which uses the ISF it chose). >0 means insulin did less than the loop expected (ISF too aggressive); <0 means it did more (ISF too weak).

## Overall: median error **1.7 mg/dL**

## By glucose at decision

| BG band | n | median error | IQR |
|---|---|---|---|
| 80-100 | 7,079 | -3.6 | -14.2–8.0 |
| 100-120 | 14,820 | -1.8 | -14.0–11.0 |
| 120-145 | 22,537 | 2.0 | -12.6–19.0 |
| 145-175 | 20,918 | 6.0 | -13.4–26.6 |
| 175-230 | 18,457 | 6.0 | -21.0–33.0 |

## By ISF formula

| formula | n | median error | IQR |
|---|---|---|---|
| no_dynisf | 84,982 | 1.0 | -15.6–21.0 |
| dynisf_sigmoid | 1,769 | 8.0 | -13.0–33.0 |
| dynisf_log | 2,242 | 21.6 | 3.0–38.0 |

## By TDD band (per-user median error)

| TDD | median error |
|---|---|
| <25 | 15.0 |
| 25-45 | 3.6 |
| 45-70 | 9.0 |
| 70+ | 0.0 |

- Loop's own autosens ratio: median 1.0; correlation of prediction error with (autosens−1) = -0.029 (if autosens were fully correcting, residual error would be ~0 and uncorrelated).

![Calibration audit](charts/inv008/fig_calibration_audit.png)

*Fasting overnight, carb-screened, full-action (4h) horizon. The loop keeps dosing over the window, so the error is the net calibration of the loop's ISF-based forecast against the realised outcome, not a pure open-loop insulin response.*