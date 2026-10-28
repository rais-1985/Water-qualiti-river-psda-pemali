# Data

The laboratory data are **not included** in this repository. They come from the river water quality
monitoring of the Water Resources Management Office (BPSDA / PSDA) Pemali Comal, Central Java, Indonesia,
analysed by the environmental laboratory of PT Cito Diagnostika Utama (CDULab). The test reports state that
results may not be reproduced without the written permission of the laboratory, so the data can only be
obtained from **BPSDA Pemali Comal** upon request and with its written permission.

## Expected format

One Excel file (`data/dataset.xlsx`, or any path given in the environment variable `WQ_DATA`),
sheet `Sheet1`, one row per sample, with these columns:

| Column | Content |
|---|---|
| `LOCATIONS` | sampling point, e.g. `SUNGAI GUNG HULU` |
| `YEAR` | 2023, 2024 or 2025 |
| `Suhu` | temperature (°C) |
| `pH` | pH |
| `DO` | dissolved oxygen (mg/L) |
| `Zat Padat Terlarut (TDS)` | TDS (mg/L) |
| `Zat Padat Tersuspensi (TSS)` | TSS (mg/L) |
| `Ammonia (NH₃-N)` | NH3-N (mg/L) |
| `BOD` | BOD (mg/L) |
| `COD` | COD (mg/L) |
| `Detergen (MBAS)` | detergent as MBAS (mg/L) |
| `Fosfat (PO₄)` | phosphate (mg/L) |
| `Minyak dan Lemak` | oil and grease (mg/L) |
| `Nitrat sebagai N (NO₃-N)` | NO3-N (mg/L) |
| `Sulfida` | sulfide (mg/L) |
| `Timbal (Pb) Terlarut` | dissolved lead (mg/L) |
| `Total Fecal Coliform` | fecal coliform (MPN/100 mL) |
| `TINGKAT PENCEMARAN` | reported status: `CEMAR RINGAN` or `CEMAR SEDANG` |

Censored values are written as text exactly as in the reports: `<0,04` (below the limit of quantitation)
or `>16000` (above the counting range). Decimal commas are accepted.
