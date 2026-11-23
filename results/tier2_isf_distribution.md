# Tier-2 v-next ISF — distribution over observed data

98 people, 5,605,658 readings. For each reading, ISF_Tier-2 = measured_ISF · √(median TDD / TDD) · g(BG) (raw, unclamped).

## Population

- Per-person median Tier-2 ISF: median 16.1 mg/dL per U (IQR 10.7–21.8).
- Pooled over every reading: median 16.5, p5–p95 6.2–43.5.
- Within a person the Tier-2 ISF spreads a median 3.28× from p10 to p90 (glucose range + TDD swings).

The per-person distribution chart (fig_tier2_per_person.png) stays in the private data folder.

![Population view](charts/inv008/fig_tier2_population.png)

## Per-person (sorted by median Tier-2 ISF)


*Raw Tier-2 calculation, no §8.2 clamp. The median sits near the person's measured ISF at typical TDD; readings above target pull ISF down (g(BG)<1, firmer corrections) and readings below target push it up. Measured ISF is well below profile ISF for most people, which is why Tier-2 doses more than the current profile — see the Tier-2 shadow evaluation.*