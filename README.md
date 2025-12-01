# basic-py-fem

basic-py-fem aims to be a "Basic, simple Python-based Finite Element Method pre-processor and solver for 1D and 2D structural problems.". It can be used as a library for Finite Element Method enthusiasts or students. Nevertheless, it wasn't written to be used as a self-sustaining way to learn the Finite Element Method. Not to mention it shouldn't be used a legitimate professional Finite Element Method Software for any kind of structural validation whatsoever!

## Background

basic-py-fem is currently being developed as I am enrolled in the Computational Mechanics course as part of my BsC in Mechanical Engineering in Instituto Superior Técnico, Lisboa. As such, I'm trying to write the code according to what was given in my classes, so it can be understood by any student of Finite Element Method/Computational Mechanics.

---

## 🚀 Features
- Simulation in the x-y plane.
- Bar 1D elements.
- Beam 1D elements.

---

## 📦 Dependencies
- numpy

---

## 🧑‍💻 Usage

In order to use basic-py-fem, you should incorporate the "elements" folder and "finite_element_method.py" file into your project folder. From then, create a file/use own python file already in your project and import:
- Node class with "from elements.node import Nodde"
- Element1D class with "from elements.element_1d import Element1D"
- FiniteElementMethod class with "from finite_element_method import FiniteElementMethod"

Depending on which elements you will use, the import lines should be
- "from elements.beam2 import Beam2" for 2-node beam elements (linear)
- "from elements.bar3 import Bar2" for 2-node bar elements (bar)

From there, the usage is really simple!
- Create your nodes with Node class
- Create your elements with the corresponding class
- Create your finite element method global problem using FiniteElementMethod class
- Apply constraints and loads
- Solve the problem with the FiniteElementMethod.solve() method
