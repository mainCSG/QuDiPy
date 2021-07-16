# QuDiPy Development Guidelines

The below guidelines are relatively obvious. The aim of this document is to standardize the way we write code for the qudipy module. Adopting a consistent framework will ensure that PRs are handled more easily and quickly.

1. Brandon (@theycallmesimba), Stephen (@bleutooth65), and Bohdan (@hromecB) are the current maintainers of the repository. All pull requests must be reviewed and approved by at least two of them before merging into master.

2. You should never work directly on the master branch. Always start your own development branch and submit a PR to merge with master.

3. Before merging you need to run all tutorials and unittests to ensure you have not broken anything.

4. External repositories are allowed as long as they are well maintained (for example numpy, scikit-learn, pandas, tqdm). This is subjective, so use your judgement.

5. Every function needs to have a docstring where there Parameters, Keyword Arguments, and Returns variables are clarified. An example of a good docstring is shown in qudipy/potential/potentialInterpolator.py. If you are using an external source as a reference (such as an academic paper or blog post), make sure to reference it in the docstring description. The general format for a docstring is

```
Parameters
----------
var_1 : type
    Description.
var_2 : type
    Description.
   
Keyword Arguments
-----------------
var_3 : type, optional
    Description.
        
Returns
-------
var_ret: type
     Description. 
```

6. If you have added code you believe accomplishes a mathematically non-trivial task, please add a short write up of the procedure to the QuDiPy math write up document: https://www.overleaf.com/3252553442tbqcmxntqvtk. You can then reference your write up in the appropriate function's docstring.

7. The qudipy module and submodules should always be loaded as:
	* `import qudipy as qd`
	* `import qudipy.potential as pot`
	* `import qudipy.chargestability as csd`
	* `import qudipy.qutils.math as qmath`
	* `import qudipy.qutils.matrices as matr`
	* `import qudipy.spinsimulator.spin_simulator as sps`

8. Common qudipy classes should be implemented as (although in somecases this may not always be appropriate):
	* `gparams = qd.potential.GridParameters()` or `gparams = pot.GridParameters()`
	* `consts = qd.Constants()`

9. Common external repositories should be loaded as:
	* `import numpy as np`
	* `import pandas as pd`
	* `import scipy as sp`
	* `from scipy import linalg as la`
	* `import matplotlib.pyplot as plt`
	* `import seaborn as sns`

8. The lowest python version supported is python 3.6.

9. Please report bugs as an issue on github.

10. Follow PEP 8 style guidelines.

11. Line length should not exceed 99 characters in the main code, and 72 characters in the docstrings.
	* **Long Inline Comments** - Example: 
	`elems_sorted = sorted(elem, reverse=True) # iterable sorted in the reversed order; necessary for the correct consecutive application of the project_up operations`
	can be changed to:
		```
		elems_sorted = sorted(elem, reverse=True)
		    # iterable sorted in the reversed order; necessary for the correct 
		    # consecutive application of the project_up operations
		```
	* **Long Function Calls and Definitions** - Example: 
	`self.__plot_heatmap(self.csd.csd_der, None, None, r'V$_1$', r'V$_2$', cbar=cbar_flag, cbar_kws=cbar_kws)`
	can be changed to:
		```
		self.__plot_heatmap(self.csd.csd_der, None, None, r'V$_1$', r'V$_2$', cbar=cbar_flag, 
				    cbar_kws=cbar_kws)
		```

		**Note**: The lines should end with a comma, and subsequent lines should ideally be indented to align with the opening bracket of the function call or definition. If indenting to align with the opening bracket makes the code less readable, indent four spaces further than the previous line.
	* **Long Arithmetic Operations** - Example: 
	`n_1 = 1/(1 - self.e_cm ** 2/(self.e_c1 * self.e_c2)) * 1/consts.e * (self.c_g1 * v_g1 * (1 - self.e_cm ** 2 / (self.e_c1 * self.e_c2)) + self.c_g2 * v_g2)`
	can be changed to:
		```
		n_1 = 1/(1 - self.e_cm ** 2/(self.e_c1 * self.e_c2)) * 1/consts.e * (self.c_g1 * v_g1 
			* (1 - self.e_cm ** 2 / (self.e_c1 * self.e_c2)) + self.c_g2 * v_g2)
		```
		**Note**: New lines should begin with the arithmetic operator (ie +, /, -, *) and be indented four spaces more than the first line of the statement
	* **Long Print Statements** - Example: 
	`n_1 = 1/(1 - self.e_cm ** 2/(self.e_c1 * self.e_c2)) * 1/consts.e * (self.c_g1 * v_g1 * (1 - self.e_cm ** 2 / (self.e_c1 * self.e_c2)) + self.c_g2 * v_g2)`
	can be changed to:
		```
		n_1 = 1/(1 - self.e_cm ** 2/(self.e_c1 * self.e_c2)) * 1/consts.e * (self.c_g1 * v_g1 
			* (1 - self.e_cm ** 2 / (self.e_c1 * self.e_c2)) + self.c_g2 * v_g2)
		```
		**Note**: The lines should end with the quotation mark (ie ", ') with no trailing space and subsequent lines should begin with the concatenation operator (ie. '+') followed by a quotation mark and space (ie " , or ' ) and ideally be indented to align with the opening bracket of the print statement. If indenting to align with the opening bracket makes the code less readable, indent four spaces further than the previous line.
	* **Long Conditional Statements** - Example: 
	`if (x_temp < self.csd.v_g1_min) or (x_temp > self.csd.v_g1_max) or (y_temp < self.csd.v_g2_min) or (y_temp > self.csd.v_g2_max)`
	can be changed to:
		```
		if (x_temp < self.csd.v_g1_min) or (x_temp > self.csd.v_g1_max) or \
		    (y_temp < self.csd.v_g2_min) or (y_temp > self.csd.v_g2_max)
		```

		**Note**: The line should end with the conditional operator (ie not/and/or) followed by a backslash character and the subsequent line should begin with the conditional statement, indented by four spaces

12. Naming should follow the following conventions (terminology consistent with: https://www.python.org/dev/peps/pep-0008/#id34):
	* **Module names** should be lowercase (with no spaces)
	* **Filenames** should be lower_case_with_underscores
	* **Class names** should be CapitalizedWords
	* **Function and method names** should be lower_case_with_underscores

13. Code should always prioritize readability first and optimality second.

14. Comment, comment, comment.

