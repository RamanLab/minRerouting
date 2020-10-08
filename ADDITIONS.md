## List of additions and changes done
- Added `pFBA` analysis codes and PSL, RSL analysis codes as separate folders inside the `core` folder.  
- Added `.mat` files inside the example folders. These files are read from the `.xml` files available on BiGG and saved as `.mat` files.  
- Added `cplex_direct` for L0 minimization with objective.
- Changed `disp()` functions in `getFastSL` to `fprintf` as the `disp()` function displayed weird characters as an effect of Character encoding.  
- Minor changes in code formating to follow the COBRA Toolbox guide. [#11](https://github.com/RamanLab/minRerouting/issues/11).

## Current Work
- Comparing results with [Serrano et al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003637#s5).
- Scrapping BiGG for Pathways of the double lethal reactions.
- Analysis of the reaction pathways.
- Analysis of the L0 minimization results.