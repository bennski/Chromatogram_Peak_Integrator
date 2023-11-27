# Chromatogram_Peak_Integrator
An example peak integration tool developed for chromatogram peak area analysis. In this example, I calculated the average and standard deviation of nucleobase peak areas from a standard.

"gates" are the programmed GC/MS/MS gates. These can be adjusted based on their retention times for your GC/MS and column.

"noise_regions" are chosen by running the program and finding the retention times in each gate that correspond to the noise floor. This has to be done manually, and is also a good time to make sure the automatically selected peak corresponds to the expected retention time from your standards. Running the code with any values here will bring up the browser tool to allow you to do this easily.

This code was developed for Pearce, Hörst, Sebree, He (2023) Planetary Science Journal (under review).
