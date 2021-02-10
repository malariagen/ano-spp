# Plasmodium rebalancing test set

24 samples submitted with increased volume to test 10 different primer concentrations:

Pool | P1/x | P2/x
-- | -- | --
1 | 1 | 1
2 | 1 | 3
3 | 1 | 9
4 | 3 | 1
5 | 3 | 3
6 | 3 | 9
7 | 9 | 1
8 | 9 | 3
9 | 9 | 9
10 | 0 | 0

Where ‘x’ is the final concentration of the P1 & P2 primers in the PCR relative to the average concentration of the primers in the balanced anopheles pool. Pool 10 is the control, having no plasmodium primers spiked in at all.

Replicate names given by Scott used, containing run ID (`P6`), source plate well, and pool number.

Sample names designed for easier visual inspection of the results: `{P1 conc}_{P2 conc}_{sample ID}`

Sample IDs:
- `PM` for plasmodium-mosquito mixes with additional data on DNA amount ratio, e.g. `e6` for 1 Plasmodium to 1,000,000 mosquito. Lab gambiae/coluzzii used as mosquito, Pf used as Plasmodium
- `C` control non-infected lab coluzzii
- `BS` dissected samples from Brandy, additional number indicates number of oocysts detected
- `NI` natural infected mosquitoes, additional data on species: `fun` - wild funestus suspected to be infected, `ste` - stephensi infected with Pb

Additional data on source samples can be found in `Notes`
