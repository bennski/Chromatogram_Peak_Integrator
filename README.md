# Chromatogram_Peak_Integrator
An example peak integration tool developed for chromatogram peak area analysis. In this example, I calculated the average and standard deviation of nucleobase concentrations in organic hazes.

start_finish_vals are based on the nucleobase standard peaks we ran prior to haze analysis, and are adjusted slightly to match the start and finish retention times for the organic haze sample analysis.

noise_floor values are chosen based on the minimum and maximum noise floor within the retention time gate for each nucleobase. This has to be done manually (by eye). Running the code with any values here will bring up the browser tool to allow you to do this easily.

This code was developed for Pearce, Hörst, Sebree, He (2023) Nature Geoscience (under review).
