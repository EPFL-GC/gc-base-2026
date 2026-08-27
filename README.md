CS 457 Geometric Computing — Assignments
======================================

## Installation notes 

You can set the environment up by running the following command:

```
conda env create -f environment.yml
conda activate gc_course_env
```

Alternatively, the following commands should produce the same outcome:

```
conda create --name gc_course_env python=3.9
conda activate gc_course_env
conda install -c conda-forge gmsh jupyterlab pytest pytest-html pytest-timeout matplotlib svgwrite meshplot triangle
pip install pyvista libigl opencv-python
```
