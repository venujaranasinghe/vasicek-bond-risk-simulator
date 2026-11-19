# Data

## cbsl_yields_weekly.csv

Weekly yield rates on Sri Lankan government securities, January 2020 to July 2026.

- **Source:** Central Bank of Sri Lanka, Data Library, "Rates on Government Securities (Per cent per annum)"
- **URL:** https://www.cbsl.lk/eResearch/MoneyMarketRatesDefault.aspx?ReportId=6169
- **Accessed:** [date you downloaded it]
- **Rights:** Data © Central Bank of Sri Lanka, reproduced with attribution for non-commercial educational use. [Check the CBSL website terms and adjust this line.]

### Columns

| Column | Meaning | Unit |
|---|---|---|
| `date` | Week-end date as published (YYYY-MM-DD) | date |
| `tbill_91`, `tbill_182`, `tbill_364` | Treasury bill yield, 91/182/364 days | % per annum |
| `bond_2y` ... `bond_30y` | Treasury bond yield by maturity | % per annum |

### Processing notes

- Values were copied from the Data Library table (the export button was not usable) and transcribed to CSV with AI assistance, then spot-checked by hand against the website.
- Blank cell = no value published for that date and maturity.
- Cells published as "-" or "0" are stored as blank (2020-10-02 and 2020-10-23 for the 91-day bill; 2023-11-23 and 2025-12-18 for bonds).
- Bond yields exist only for weeks in which that maturity was auctioned.
- No values were altered, interpolated or filled.