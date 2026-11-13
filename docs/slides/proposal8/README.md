# Slide deck for the dynamic ISF proposal

Source of the 20-slide deck that accompanies `docs/dynamic-isf-8-proposal.md`. Every number on
the slides is quoted from that document. The deck was written on 13 November 2026 and is
published as a private Slides artifact at https://claude.ai/artifact/PduSp6kFEsUZt42hCAtFTq.

`project/deck.json` holds the slide order, `project/deck.css` the shared styles, and
`project/slides/` one file per slide with its speaker notes in the closing `<aside>`.

The images are held by the artifact and are referenced by asset id. Each is a copy of a figure
in `charts/proposal8/`:

| Slide | Figure |
|---|---|
| published | fig1_published_curve.png |
| compare | fig2_scaled_to_profile.png |
| link | fig3_profile_or_daily_dose.png |
| multiply | fig4_multiply_or_shift.png |
| simulation | fig5_outcomes_in_simulation.png |
| lows | fig6_what_preceded_the_lows.png |

The colours and typeface follow diabettech.com: Open Sans, cyan #00aae1, slate #333f49 and pale
blue-grey #e9ecf6. All 20 slides were rendered locally in
headless Chrome at 1920 by 1080 on 13 November 2026 and read for overflow. The figures were not
loaded in those renders, so their appearance on the slide has not been checked.

The text follows the register of the paper: plain label headings, impersonal voice and no
contractions. `punctuation_pass.py` and `style_pass.py` record the two text passes made after
the first build, each as a list of exact replacements.
