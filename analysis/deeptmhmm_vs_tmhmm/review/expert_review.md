# Expert review: TMHMM 2.0 vs DeepTMHMM 1.0 disagreements (80 fungal proteins)

Per-protein calls are in `expert_review.tsv`. This file gives the summary.

## Method

- One reviewer (Claude) judged each protein from:
  - the sequence: Kyte-Doolittle mean over each called segment, flank charges, N-terminal and C-terminal patterns;
  - the Pfam domains;
  - the best SwissProt hit;
  - general knowledge of protein families.
- No web lookups and no new predictors were used.
- The reason column marks calls that rest mainly on family knowledge with "FAMILY KNOWLEDGE".
- In this review, the C-terminal hydrophobic segment of a GPI-anchored protein is not counted as a TM helix. That segment is cleaved off when the GPI anchor is attached, so the mature protein has no TM there.
- Row 16 and row 52 are the same protein (F985D20D_004311-T1), which appears in groups 1 and 4. There are 79 unique proteins.

## Verdict tally

| Group | n | TMHMM | DeepTMHMM | both wrong | unclear |
|---|---|---|---|---|---|
| 1. C-terminal single helix, DeepTMHMM none | 16 | 5 | 7 | 1 | 3 |
| 2. N-terminal helix, SignalP OTHER, DeepTMHMM SP | 16 | 12 | 1 | 0 | 3 |
| 3. TMHMM >= 1 helix, DeepTMHMM globular | 16 | 8 | 1 | 2 | 5 |
| 4. PredGPI positive, helix counts differ | 16 | 4 | 11 | 0 | 1 |
| 5. TMHMM >= 7 helices, counts differ by >= 2 | 16 | 1 | 12 | 0 | 3 |
| **Total** | 80 | 30 | 32 | 3 | 15 |

Confidence:
- High: 31 calls.
- Medium: 21 calls.
- Low: 28 calls.
- Every "unclear" call has low confidence.

## Main patterns

1. **Type II signal anchors (group 2): TMHMM is usually right.**
   - 12 of 16 proteins favour TMHMM.
   - Most are Golgi or ER type II proteins:
     - glycosyltransferases: OCH1, MNN9/Anp1, LARGE-like GT49, GT25, GT90, Fringe/DUF604;
     - Grx6/7-type glutaredoxins;
     - MAPR;
     - V-ATPase subunit e.
   - These proteins have a positively charged cytosolic N-tail and an 18-25 aa hydrophobic core. DeepTMHMM labels them as signal peptides. SignalP 6 correctly says OTHER.
   - DeepTMHMM's SP call is right only when the gene model has a wrong upstream start (row 17; row 5 in group 1). There, an internal Met gives a classic SP.

2. **Globular calls on real anchors (group 3): DeepTMHMM misses single anchors.**
   - DeepTMHMM calls "globular" on several well-known single-anchor proteins:
     - the microsomal P450 N-terminal anchor;
     - the cytochrome b5 reductase N-terminal anchor;
     - the Tom20 signal anchor;
     - a clear coiled-coil tail-anchored protein.
   - TMHMM is right in 8 of 16 proteins.
   - TMHMM's typical false positive is a Rossmann GxGxxG or GxGxxP dinucleotide-binding motif called as a helix: D-amino-acid oxidase (row 41), FAD monooxygenase (row 51), and the second helix of row 47.

3. **GPI proteins (group 4): DeepTMHMM is usually right (11 of 16).** TMHMM makes two errors:
   - It calls the N-terminal SP h-region a helix (rows 53, 55, 58, 61, 62).
   - It calls the C-terminal GPI attachment signal a helix (rows 49, 50, 52, 62, 64).

   Exceptions:
   - In 3 proteins (rows 54, 59, 60), the count difference is a DeepTMHMM artifact. DeepTMHMM splits one helix in two with a 1-residue loop, and the helix positions agree. In these 3 proteins, PredGPI's call is also probably false, because there is no SP.
   - Row 63 is a soluble FMN reductase. DeepTMHMM calls a short helix in a probable gene-model extension.

4. **Large multi-pass proteins (group 5): DeepTMHMM is usually right (12 of 16).**
   - Eight are MFS transporters. DeepTMHMM gives 12 helices in each, which matches the family. TMHMM gives 8-10, because it misses amphipathic helices or merges helix pairs.
   - DeepTMHMM's raw helix count is often inflated: it splits single helices with 1-3 residue "loops" (rows 65, 70, 72, 73, 74). Count helices only after merging segments separated by 3 residues or fewer.
   - The one TMHMM win is ALG6 (row 79). DeepTMHMM calls only 6 helices for a protein family with more than 10.

5. **C-terminal single helix (group 1): the result is mixed.**
   - The two tools are about even: TMHMM 5, DeepTMHMM 7.
   - TMHMM is right for tail anchors: Tom7, metaxin, an FHA protein with a Gly/Ala tail anchor, and a GT4 C-terminal anchor.
   - DeepTMHMM is right when TMHMM's segment is either:
     - a GPI signal (yapsins), or
     - a weak segment (KD < 1.2) in a soluble protein or in a wrong 3' exon.
   - Five of the 16 "proteins" in this group are short ORFs (50-87 aa) with no Pfam and no homolog. Several have the hydrophobic segment at the N-terminus, not the C-terminus.

## Would "drop TMHMM's N-terminal helix when SignalP calls SP" help?

- **It helps in this sample.** SignalP 6 calls an SP in 6 proteins where TMHMM calls an N-terminal helix (rows 33, 53, 55, 58, 61, 62).
  - In 5 of them, the protein is a secreted or GPI-anchored protein, and the "helix" is the SP h-region.
  - In row 33, the SP call is supported by a K/R-rich n-region and a 15-aa h-region.
  - No case here shows SignalP calling an SP on a known signal anchor.
- **A rule based on DeepTMHMM's SP call would hurt.**
  - In group 2, DeepTMHMM calls an SP and SignalP says OTHER.
  - Dropping TMHMM's helix there would remove a real type II anchor in about 12 of 16 cases. Most of these are Golgi glycosyltransferases.
  - Do not use DeepTMHMM "SP" alone to drop an N-terminal helix.
- **A similar rule for the C-terminus would also help.**
  - Drop a C-terminal TMHMM helix when PredGPI is positive and an SP is present.
  - Supporting cases: rows 6, 16/52, 49, 50, 58, 62, 64.
- **Without an SP, the PredGPI call is unreliable.**
  - Rows 51, 54, 59, 60 and 63 are PredGPI-positive without an SP.
  - These are cytosolic or mitochondrial proteins, or multi-pass proteins.
  - Require an SP before trusting PredGPI.
- **A possible combined rule follows from this.** This sample suggests one, but has not tested it:
  - Use DeepTMHMM's helix set for proteins with more than 6 helices, after merging split helices.
  - Use TMHMM for single N-terminal anchors when SignalP says OTHER.
  - Remove SP and GPI-signal segments when SignalP or PredGPI support them.

## Gene-model issues seen

- **Wrong upstream start.** An internal Met gives a clean SP (rows 5, 17, 27, 53, 56).
- **Probable fusions.**
  - Row 46: Gta3 + Tom20 + SET.
  - Row 47: b5 reductase + an unrelated C-terminal region.
  - Row 70: a 923-aa dolichol kinase model.
  - Row 65: FAD monooxygenase + CtvD-like domain. This one may be real.
- **Spurious short ORFs** with no Pfam and no homolog (rows 1, 4, 10, 22, 42, 44).
- **Extensions outside the Pfam domain.** These may come from a wrong exon: the N-terminal extension in row 63, the C-terminal extensions in rows 11 and 32, and the low-complexity N-terminus in row 57.
- **Missing N-terminal SP** in apparent GPI proteins (rows 9, 15).

## Caveats

- **Small samples.** Each group has only 16 proteins. Tallies of 12/16 or 11/16 are not precise rates.
- **Biased sampling.** The groups were chosen for disagreement, so these numbers say nothing about overall accuracy.
- **No experimental truth.** There is no experimental topology for any of these proteins.
  - Many verdicts rest on what is known for the protein family in model organisms.
  - Remote homology (for example, identity below 35%) may not preserve topology.
- **Few homologs.** 26 of the 80 rows have neither a Pfam domain nor a SwissProt hit. For these, the call rests on hydrophobicity and charge alone. Most of them are "unclear" or "low".
- **Simple hydrophobicity measure.** Kyte-Doolittle averages are crude. Proline- or glycine-containing TMs, such as Tom7, score low but are real.
- **Unchecked gene models.** Gene-model notes are inferences from the sequence. They were not checked against RNA-seq or genome alignments.
- **Single reviewer.** No second reviewer has checked these calls.
