# Finite Field gf.cpp

## gf.cpp
[x] additive inverse
[ ] multiplicative inverse
[x] additive table
[x] multiplicative table

The implementation doesn't use intrinsics, instead it uses dynamic array for the flexibility of creating field, which is not as optimized as the intrinsics approach. However, due to the complexity of multiplicative inverse in finite field (especially gf), it is not implemented.

multiplicative inverse reference:
1. [https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/](https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/).
2. [https://en.wikipedia.org/wiki/Finite_field_arithmetic#Multiplicative_inverse](https://en.wikipedia.org/wiki/Finite_field_arithmetic#Multiplicative_inverse).
3. [https://www.youtube.com/watch?v=AcJEyvzUT64&ab_channel=Ch-13ComputerScienceandEngineering](https://www.youtube.com/watch?v=AcJEyvzUT64&ab_channel=Ch-13ComputerScienceandEngineering).
