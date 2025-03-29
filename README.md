EOSsolver.py: A Python package for calculating physical properties through equation of state

Jinghuan Liu, Weijing Zeng and Xufeng Lin

1 College of Chemistry and Chemical Engineering, China University of Petroleum (East China), Qingdao, P. R. China, 266580
2 School of Mathematics and statistics, Wuhan University, Wuhan, P.R.China, 430070
3 State Key Laboratory of Heavy Oil Processing, China University of Petroleum (East China), Qingdao, P. R. China, 266580

Summary
	Equation of state, provided with a few empirical parameters, can give out plenty of properties of the fluid, making it wieldy accepted and frequently used in various fields, including, hydrodynamics, chemical engineering, etc.
Originating from the ideal gas equation, which is limited because of neglecting the atomic volume and interactions between molecules, the van der Waals equation of state introduces two empirical parameters to tackle these deficiencies, and successfully expand the utility of equation of state. After that, aiming at improving the availability and universality, tons of work have been put into the field and multiple versions of equations of state have been derived by many researchers by introducing more empirical parameters. Therefore, acceptable accuracy has been achieved for both pure and mixture, when describing alkane, carbon dioxide, etc.
However, it’s still inconvenient for researchers to utilize these equations manually, because it requires enormous work, particularly when dealing with mixtures, since the binary interaction factors are hard to determine, compressibility factors are hard to solve from cubic equations, etc. To address these difficulties, we developed a python package to help researchers get the physical properties that they want in a much easier and quicker way. The current version strictly follows the P-R equation of state and the corresponding mixing rule (Peng & Robinson, 1976). Instead of being obtained from the empirical measurements, the binary interaction parameters are calculated based on formulas considering both the solvent effect and temperature effect. (Moysan, Paradowski & Vidal, 1986). The current version only supports the computation of fugacity coefficients and compressibility factors. Developers tend to consummate their work in the near future.


References
Peng, D. Y., & Robinson, D. B. (1976). A new two-constant equation of state. Industrial & Engineering Chemistry Fundamentals, 15(1), 59-64.
Moysan, J. M., Paradowski, H., & Vidal, J. (1986). Prediction of phase behaviour of gas-containing systems with cubic equations of state. Chemical engineering science, 41(8), 2069-2074.
