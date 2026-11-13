"""Punctuation pass over the proposal deck: exact replacements, each asserted to match once."""
import pathlib, sys
root = pathlib.Path(sys.argv[1]) / "project/slides"
E = {
"cover": [
 ("Diabettech · Where Diabetes and Technology meet", "Diabettech: Where Diabetes and Technology meet"),
 ("Tim Street, Diabettech · November 2026 · Draft for discussion", "Tim Street, Diabettech, November 2026. Draft for discussion."),
],
"origins": [
 ("Wilson C., talk to Loop and Learn, 6 February 2022. AndroidAPS source history.", "Chris Wilson’s talk to Loop and Learn, 6 February 2022, and the AndroidAPS source history."),
],
"published": [
 ("Romero-Ugalde et al., Diabetes 2025;74(Supplement 1):920-P, doi:10.2337/db25-920-P</p>", "Romero-Ugalde et al., Diabetes 2025;74(Supplement 1):920-P (doi:10.2337/db25-920-P).</p>"),
],
"cohorts": [
 (">people across the eight cohorts<", ">People across the eight cohorts<"),
 (">fall from 100 to 200 mg/dL within a person, by cohort<", ">Fall from 100 to 200 mg/dL within a person, by cohort<"),
 (">the published fall over the same range<", ">The published fall over the same range<"),
 (">REPLACE-BG (2.25 to 2.82), adults dosing with no algorithm<", ">REPLACE-BG, adults dosing with no algorithm (2.25 to 2.82)<"),
],
"link": [
 (">3,433 corrections from 150 adults in REPLACE-BG, each given by the person with no algorithm dosing.<", ">The test used 3,433 corrections from 150 adults in REPLACE-BG, each given by the person with no algorithm dosing.<"),
 (">on the entered ISF (0.55 to 1.32). The profile link predicts 1 and the TDD link predicts 0.<", ">Coefficient on the entered ISF (0.55 to 1.32), where the profile link predicts 1 and the TDD link predicts 0<"),
 (">on daily dose (−0.43 to 0.26). The profile link predicts 0 and the TDD link predicts −1.<", ">Coefficient on daily dose (−0.43 to 0.26), where the profile link predicts 0 and the TDD link predicts −1<"),
 ("REPLACE-BG, 150 adults giving their own corrections. Between-person comparison. Loop users and children not tested.", "A between-person comparison in 150 REPLACE-BG adults giving their own corrections. Loop users and children were not tested."),
],
"rule": [
 ("Curve: Romero-Ugalde et al., doi:10.2337/db25-920-P</p>", "The curve is from Romero-Ugalde et al. (doi:10.2337/db25-920-P).</p>"),
],
"worked": [
 ("One person: profile ISF 50 mg/dL per unit, 40 U a day, every setting at its default.", "This person has a profile ISF of 50 mg/dL per unit and uses 40 U a day, with every setting at its default."),
 ("Sensitivity in mg/dL per unit. AndroidAPS with a divisor of 55 at 100%. Trio at its default limits.", "Sensitivity in mg/dL per unit, with AndroidAPS on a divisor of 55 at 100% and Trio at its default limits."),
],
"compare": [
 ("Worked person: profile ISF 50, 40 U a day, default settings. The crossovers against AndroidAPS and Trio are this person’s only.", "The worked person has a profile ISF of 50 and uses 40 U a day on default settings. The crossovers against AndroidAPS and Trio apply to this person only."),
],
"correction": [
 ("Units to bring the same person down to 100 mg/dL, from the sensitivity alone.", "These are the units to bring the same person down to 100 mg/dL, worked from the sensitivity alone."),
 (">Profile ISF, unchanged<", ">Unchanged profile ISF<"),
 ("Units of insulin. Worked person: profile ISF 50, 40 U a day, default settings.", "Units of insulin for the worked person, with a profile ISF of 50 and 40 U a day on default settings."),
],
"limits": [
 ('<p class="stat-num">109 to 158</p><p class="stat-label">mg/dL where the curve sits inside the default limits, 49 of the 280 mg/dL it covers<', '<p class="stat-num">109 to 158 mg/dL</p><p class="stat-label">Where the curve sits inside the default limits, 49 of the 280 mg/dL it covers<'),
 ('<p class="stat-num">3.33×</p><p class="stat-label">the profile’s insulin that the curve asks for at 380 mg/dL, and 2.41× at 300<', '<p class="stat-num">3.33 times</p><p class="stat-label">The profile’s insulin, asked for by the curve at 380 mg/dL, and 2.41 times at 300<'),
 (">a correction from 300 to 100 mg/dL as the curve gives it<", ">A correction from 300 to 100 mg/dL as the curve gives it<"),
 (">the same correction inside the default limits<", ">The same correction inside the default limits<"),
],
"simulation": [
 (">40 simulated adults ran oref0 in TimSim for two months, and only the sensitivity rule differed between arms.<", ">In TimSim, 40 simulated adults ran oref0 for two months, and only the sensitivity rule differed between arms.<"),
 ("below 54 was up 0.17 (0.04 to 0.32), unresolved.</p>", "below 54 was up 0.17 (0.04 to 0.32), which is unresolved.</p>"),
 ("Margins 0.25 and 1.25 points, a quarter of the consensus limits (Battelino 2019). Two simulated months pooled.", "The margins are 0.25 and 1.25 points, a quarter of the consensus limits (Battelino et al. 2019), over two simulated months pooled."),
 ("above 250 was +0.09 (−0.12 to 0.29), cleared.", "above 250 was +0.09 (−0.12 to 0.29), which is cleared."),
],
"lows": [
 ("First simulated month. The extra time below 54 mg/dL against the unchanged profile was 0.43 points (0.24 to 0.61).", "In the first simulated month, the extra time below 54 mg/dL against the unchanged profile was 0.43 points (0.24 to 0.61)."),
 ("added +0.03.", "added 0.03."),
],
"stages": [
 (">Run, two simulated months<", ">Run for two simulated months<"),
 (">Code written, no records yet<", ">Code written and no records yet<"),
 (">5. Use, own device first, then released lines<", ">5. Use on my own device first and then on released lines<"),
 ("Shadow running is where the code now is: written, with no records yet.", "Shadow running is where things now stand, with the code written and no records yet."),
],
"open": [
 (">Tested only between adults giving their own corrections.", ">It was tested only between adults giving their own corrections."),
 (">Whether the steep part between 100 and 150 mg/dL should be dosed from, and how far to follow the fall above 200.<", ">It’s open whether the steep part between 100 and 150 mg/dL should be dosed from, and how far to follow the fall above 200.<"),
 ("From “What is not settled” in the paper.", "Taken from the section “What is not settled” in the paper."),
],
"asks": [
 ("keep the 0.7 to 1.2 limits, or have its own?", "keep the 0.7 to 1.2 limits or have its own?"),
],
"sources": [
 ("“An alternative proposal for dynamic ISF in open-source AID systems”, with the scripts behind each number, at github.com/tim2000s/dynamic-isf-calculations<", "The scripts behind each number are at github.com/tim2000s/dynamic-isf-calculations, with the paper, “An alternative proposal for dynamic ISF in open-source AID systems”.<"),
 ("Romero-Ugalde et al., Diabetes 2025;74(Supplement 1):920-P, doi:10.2337/db25-920-P<", "Romero-Ugalde et al., Diabetes 2025;74(Supplement 1):920-P (doi:10.2337/db25-920-P).<"),
 ("Wilson C., talk to Loop and Learn, 6 February 2022, youtube.com/watch?v=oL49FhOts3c<", "Chris Wilson’s talk to Loop and Learn on 6 February 2022 (youtube.com/watch?v=oL49FhOts3c).<"),
 ("Kovatchev et al. 2009 and Xie, Simglucose 2018, for the simulator. Battelino et al. 2019 for the consensus limits.", "Kovatchev et al. (2009) and Xie’s Simglucose (2018) for the simulator, and Battelino et al. (2019) for the consensus limits."),
],
}
for sid, pairs in E.items():
    f = root / f"{sid}.html"; t = f.read_text(encoding="utf-8")
    for old, new in pairs:
        assert t.count(old) == 1, (sid, old[:50], t.count(old))
        t = t.replace(old, new)
    f.write_text(t, encoding="utf-8")
print("edited", len(E), "slides,", sum(map(len, E.values())), "replacements")
