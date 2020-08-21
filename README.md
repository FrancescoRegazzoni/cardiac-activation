# Mathematical models of active force generation in cardiomyocytes #

This repository contains codes implementing some microscale mathematical models of cardiomyocytes.
The first model only describes the dynamics of the so-called *regulatory units* (RU), namely Troponin and Tropomyosin, while the other models also incorporate the description of *crossbridges* (XB).

| Model    | # variables | RU dynamics        | XB dynamics        | Reference                                                                                                     |
|:---------|------------:|:-------------------|:-------------------|:--------------------------------------------------------------------------------------------------------------|
| RDQ18    | 2176        | :heavy_check_mark: | :x:                | [[1]](https://doi.org/10.1007/s10237-018-1049-0)                                                              |
| RDQ20-SE | 1408        | :heavy_check_mark: | :heavy_check_mark: | [[2]](http://hdl.handle.net/10589/152617), [[3]](https://arxiv.org/abs/2004.07910) (herein denoted as SE-ODE) |
| RDQ20-MF | 20          | :heavy_check_mark: | :heavy_check_mark: | [[2]](http://hdl.handle.net/10589/152617), [[3]](https://arxiv.org/abs/2004.07910) (herein denoted as MF-ODE) |

*Remark:* a reduced version of the RDQ18 model, based on Artificial Neural Networks and built with the Machine Learning library [model-learning](https://github.com/FrancescoRegazzoni/model-learning), is available at [cardiac-activation-ann](https://github.com/FrancescoRegazzoni/cardiac-activation-ann).

### Repository structure

- `models_matlab`: Matlab codes.
- `models_python`: Python codes.
- `params`: parameters of the models.

### References

- [1] F. Regazzoni, L. Dedè, A. Quarteroni ["Active contraction of cardiac cells: a reduced model for sarcomere dynamics with cooperative interactions"](https://doi.org/10.1007/s10237-018-1049-0), *Biomechanics and Modeling in Mechanobiology* (2018).
- [2] F. Regazzoni, ["Mathematical modeling and Machine Learning for the numerical simulation of cardiac electromechanics"](http://hdl.handle.net/10589/152617), *PhD Thesis - Politecnico di Milano* (2020).
- [3] F. Regazzoni, L. Dedè, A. Quarteroni ["Biophysically detailed mathematical models of multiscale cardiac active mechanics"](https://arxiv.org/abs/2004.07910), *submitted* (2020).

### Author

Francesco Regazzoni, MOX - Politecnico di Milano (<francesco.regazzoni@polimi.it>)
