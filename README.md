# Special adjoint orbits through the combinatorics of Vogan diagrams

This repository contains codes to check if an adjoint orbit of a simple Lie algebra is special. Based on [Almost Kähler geometry of adjoint orbits of semisimple Lie groups](https://link.springer.com/article/10.1007/s00209-022-02995-9) and arXiv:...

## Jupyter notebook (SageMath)

The notebook ``Special vector through the combinatorics of Vogan diagrams.ipynb`` allows to play with the different types of Lie algebras (classical and exceptional) and check they speciality and check the Lie algebra of the orbit and the stabilizer. It allows also to print Vogan diagrams.

## SageMath scripts

### Check speciality of a single adjoint orbit

This script works for all types of rank. Just run ``specialCombinatoricsVD.sage`` with the following arguments:
- ``type``: the Lie type of the simple Lie algebra (default 'A')
- ``rank``: the rank of the Lie algebra (default 10)
- ``s``: list of indices of painted vertices (default 1)
- ``cores``: number of cores to be used in the computation (default 1)

The script outputs the underlying complex Lie algebra, its Dynkin diagram and its dimension and the indices of painted vertices. Then, if the diagram is special it returns the type (symplectic general type, symplectic Calabi-Yau or symplectic Fano), the Lie algebra of the orbit and the stabilizer, the Hermitian scalar curvature and the special vector.

This code is the one used in arXiv:...


### Check speciality and info of a single adjoint orbit (small diagrams)

It is better to use this code to check Lie algebras up to rank 15. Just run ``specialVoganDiagrams.sage`` with the following arguments:
- ``type``: the Lie type of the simple Lie algebra (default 'A')
- ``rank``: the rank of the Lie algebra (default 10)
- ``s``: list of indices of painted vertices (default 1)

The script outputs the underlying complex Lie algebra, its Dynkin diagram and its dimension and the indices of painted vertices. Then, if the diagram is special it returns the type (symplectic general type, symplectic Calabi-Yau or symplectic Fano), the Lie algebra of the orbit and the stabilizer, the Hermitian scalar curvature and the special vector.

This code the one used in the paper [Almost Kähler geometry of adjoint orbits of semisimple Lie groups](https://link.springer.com/article/10.1007/s00209-022-02995-9)

## Check all special orbits for a given type and rank of simple Lie algebra (small diagrams)

It is better to use this code to check Lie algebras up to rank 10. Just run ``specialVoganDiagramSingle.sage`` with the following arguments:
- ``type``: the Lie type of the simple Lie algebra (default 'A')
- ``rank``: the rank of the Lie algebra (default 10)

The script returns all special adjoint orbits associated to the given Lie algebra, together with the underlying complex Lie algebra, its Dynkin diagram and its dimension and the indices of painted vertices. Then, if the diagram is special it returns the type (symplectic general type, symplectic Calabi-Yau or symplectic Fano), the Lie algebra of the orbit and the stabilizer, the Hermitian scalar curvature and the special vector.

This code the one used in the paper [Almost Kähler geometry of adjoint orbits of semisimple Lie groups](https://link.springer.com/article/10.1007/s00209-022-02995-9)
