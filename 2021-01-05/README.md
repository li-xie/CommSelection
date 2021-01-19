# Simulate community selection on HM communities:
- Spike with specified numbers of H and M clones.
- Calculate heritability with slope from least square linear regression of parent-offspring community function after removing outliners.
- Check heritability when selection stalls.
- Switch to a spiking strategy with highest heritability.
- Communities reproduce through pipetting.
- 4 phenotypes of M and 2 phenotypes of H can be modified by mutations.
## To run the codes:
- modify SixPhenoParaInit.m, which includes all the simulation parameters
- decision_func.m examines the past selected community function to decide whether to check heritability in 2 cycles. decision_func.m in 14, 15, 20-23 are different from those in other folders.  
