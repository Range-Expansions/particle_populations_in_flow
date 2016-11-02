# Particle Populations in Flow

A cython instantiation of the algorithm discussed in 

Pigolotti S. et al. 
*Growth, competition and cooperation in spatial population genetics*. 
Theoretical Population Biology (2013), 84C, 72–86. https://doi.org/10.1016/j.tpb.2012.12.002

Currently, an arbitrary number of populations with differing growth rates
compete to eat nutrients while they are advected by fluid flow.

To install, use

```python
python setup.py develop
```

Installation requires Cython and CythonGSL (you can install CythonGSL 
with pip). 

See the examples folder for ipython notebooks showing how to use the code.