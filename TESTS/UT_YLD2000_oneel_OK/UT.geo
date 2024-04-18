*Node, nset=NALL
      1,  0.1,  0.1,  0.1
      2,  0.1,  0.0,  0.1
      3,  0.1,  0.1,  0.0
      4,  0.1,  0.0,  0.0
      5,  0.0,  0.1,  0.1
      6,  0.0,  0.0,  0.1
      7,  0.0,  0.1,  0.0
      8,  0.0,  0.0,  0.0
*Element, type=C3D8R, elset=ELALL
1, 5, 6, 8, 7, 1, 2, 4, 3
*Nset, nset=YMAX, generate
 1,  5,  2
*Nset, nset=YMIN, generate
 2,  6,  2
*Nset, nset=XMIN, generate
 5,  8,  1
*Nset, nset=XMAX, generate
 1,  4,  1
*Nset, nset=ZMIN
 3, 4, 7, 8
*Nset, nset=LOAD
 7,
*Elset, elset=ONSET
 1,
*Nset, nset=BASE
 8,