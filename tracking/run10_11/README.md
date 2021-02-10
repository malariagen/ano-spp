# Plasmodium rebalancing third experiment

Samples: An. stephensi were fed with a P. falciparum gametocyte culture, and two feed experiments were included - feed15, which had a prevalence of 60% and an average oocyst intensity of 5 per gut; and feed16, which had a prevalence of 95% and an average oocyst intensity of 18 per gut. As before, 8 control uninfected mosquitoes were included. Feed 15 was sampled at days 0, 8, and 13/14 and Feed 16 was sampled at day 0 and 9. Each time point included 16 individual mosquitoes for a total of 88 samples.

Primer concentrations were set as 80x P1/20x P2. Sample naming includes conditions, example `e2_f16_9d2_55s`:

- e[1|2] - extraction 1/2. Initial non-destructive and subsequent destructive. Both done with buffer C. 
- f1[5|6] - feed 15/16
- 9d2 - 9th day, sample 2
- 51/55 - PCR1 temperature. Default 55C, decreased to help with plasmodium amplification. 
- c/s - standard cycling/subcycling for PCR2. Subcycling differs from standard Sample Barcoding PCR cycling by the use of an oscillating temperature in the annealing/extension phase of the PCR, thus hoping to improve the primer balance (Liu & Sommer, 1998). This is proposed to help when amplifying multiple targets with widely varying %GC contents as the lower extension temperatures required for lower %GC targets are suboptimal for higher %GC targets and vice versa, so the oscillating temperature enables a compromise to be achieved between the requirements of the different targets. Subcycling PCR was carried out as follows: 95°C hold (PCR plate transferred directly from 4°C cooled Mosquito deck onto thermocycler block, then rest of protocol commenced); 31 cycles of 95°C for 20 seconds (denaturation) followed by annealing/extension using 4 cycles of 68°C for 15 seconds & 60°C for 15 seconds per cycle (i.e. a 4 cycle oscillation is nested within each of the 31 PCR cycles); 68°C for 3 minutes (final extension) followed by a 4°C hold.

For tagging, plates were iterated as follows: extractions - within plate; 1/2 - plate 10/11; 1-2/3-4 - 51/55; 1-4/5-8 - standard/subcycling.
