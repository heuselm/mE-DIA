* Rules *
| #  | Enzyme         | Sense | Cleaves At | Not If Followed By |
|----|----------------|-------|------------|---------------------|
| 0  | No_enzyme      | 0     | -          | -                   |
| 1  | Trypsin        | 1     | KR         | P                   |
| 2  | Trypsin_P      | 1     | KR         | -                   |
| 3  | Lys_C          | 1     | K          | P                   |
| 4  | Lys_N          | 0     | K          | -                   |
| 5  | Arg_C          | 1     | R          | P                   |
| 6  | Asp_N          | 0     | D          | -                   |
| 7  | CNBr           | 1     | M          | -                   |
| 8  | Glu_C          | 1     | DE         | P                   |
| 9  | PepsinA        | 1     | FL         | P                   |
| 10 | Chymotrypsin   | 1     | FWYL       | P                   |
| 11 | Elastase       | 1     | AV         | -                   |
| 11b | Elastase (extended), code `EX` | 1 | AVSGLI | -              |

*Description*
-  The first column of the parameter definition is the enzyme number. This number list must start from 0 and sequentially increase by 1. The second column is the enzyme name; no spaces are allowed in this name field. The third column is the digestion "sense" i.e. a value of "0" specifies cleavage N-teriminal to (before) the specified residues in column 4 and a value of "1" specifies cleavage C-terminal to (after) the specified residues in column 4. Column 4 contains the residue(s) that the enzyme cleaves at. Column 5 contains the flanking residue(s) that negate cleavage.

## Rule 11b: Elastase, extended (`EX`)

The elastase used in this project is **Promega Elastase, Cat. # V1891** (porcine
pancreas). Promega describe it as cleaving C-terminal to **Ala, Val, Ser, Gly, Leu
and Ile**:

- product page: https://www.promega.com/products/mass-spectrometry/proteases-and-surfactants/elastase/
- product information sheet (literature # 9PIV189): https://www.promega.com/resources/protocols/product-information-sheets/n/elastase-protocol

Rule 11 (`E`, A and V only) is unchanged, because changing it would alter the
meaning of every existing result folder. `EX` sits next to it: DIA-NN
`--cut A*,V*,S*,G*,L*,I*`, no proline rule, identical to `E` apart from the residues.

**`EX` has its own missed-cleavage default: 3.** DIA-NN counts every internal site
of the `--cut` rule as a missed cleavage, and A, V, S, G, L and I make up 44.8 % of
*S. pyogenes* residues (A and V alone: 14.5 %). So extending the residues also turns
the internal S/G/L/I of today's A/V peptides into missed cleavages. The table below
counts fully specific peptides of length 7–30 from an in-silico digest that
reproduces DIA-NN's own pass-1 elastase library (79,952 elution groups for A/V at
mc = 1 on *S. pyogenes*). "Retained" is the share of the elastase peptides the
current runs identified (A/V-ending, fully specific search) that stay inside the
six-residue search space. "Today" is A/V at mc = 1: 79,952 peptides (*S. pyogenes*)
and 221,933 (*E. coli*). Workstation measurement, issue #35.

| six-residue rule | *S. pyogenes*: peptides (× today) | retained | *E. coli*: peptides (× today) | retained |
|---|---|---|---|---|
| mc = 1 | 44,358 (0.6×) | 16 % | 105,182 (0.5×) | 16 % |
| mc = 2 | 142,852 (1.8×) | 43 % | 351,358 (1.6×) | 43 % |
| **mc = 3** | **305,095 (3.8×)** | **68 %** | **772,811 (3.5×)** | **70 %** |
| mc = 4 | 507,595 (6.3×) | 86 % | 1,311,189 (5.9×) | 87 % |
| mc = 5 | 725,039 (9.1×) | 94 % | 1,897,590 (8.6×) | 96 % |

At the global mc = 1 the six-residue rule would *shrink* the elastase search space
by about 45 % and exclude 84 % of today's elastase identifications. mc = 3 keeps
about 70 % at about 3.5–3.8× today's search space. Every other enzyme keeps the
global value. Override per run with `--enzyme-mc EX=N`.
