# Special adjoint orbits through the combinatorics of Vogan diagrams

This repository contains codes to check if an adjoint orbit of a simple Lie algebra is special. Based on [Almost Kähler geometry of adjoint orbits of semisimple Lie groups](https://link.springer.com/article/10.1007/s00209-022-02995-9) and arXiv:...

## Jupyter notebook (SageMath)

This notebooks allows to play with the different types of Lie algebras (classical and exceptional) and check they speciality. It allows also to print Vogan diagrams.

## SageMath script

This script works for all types of rank. Just run ... with the following arguments:
- type: the Lie type of the simple Lie algebra (default 'A')
- rank: the rank of the Lie algebra (default 10)
- s: list of indices of painted vertices (default 1)
- cores: number of cores to be used in the computation (default 1)

The script outputs the underlying complex Lie algebra, its Dynkin diagram and its dimension and the indices of painted vertices. Then, if the diagram is special it returns the type (symplectic general type, symplectic Calabi-Yau or symplectic Fano), and the runtime.

## SageMath script to check speciality and info of a single adjoint orbit

It is better to use this code to check Lie algebras up to rank 15. Just run ... with the following arguments:
- type: the Lie type of the simple Lie algebra (default 'A')
- rank: the rank of the Lie algebra (default 10)
- s: list of indices of painted vertices (default 1)


## SageMath script to all special orbits for a given type and rank of simple Lie algebra

It is better to use this code to check Lie algebras up to rank 10. Just run ... with the following arguments:
- type: the Lie type of the simple Lie algebra (default 'A')
- rank: the rank of the Lie algebra (default 10)

The script returns all special adjoint orbits associated to the given Lie algebra, together with 
