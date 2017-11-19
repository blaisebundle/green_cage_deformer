Crude implementation of [Green Coordinates](http://www.wisdom.weizmann.ac.il/~ylipman/GC/gc_techrep.pdf)

Using [Eigen](http://eigen.tuxfamily.org)

TODOs:
* extending logic to cage's exterior
* multi-index support
* deform membership handling
* exporting out green coordinates and original cage edges/normals
* spatial hashing to calculate only on modified verts/tris
* CUDA or OpenCL support
* code cleanup (separate green coords, etc) :P