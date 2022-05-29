# Finite Field - n-hop part-time

## the intrinsics used in FiniteField project are ARM NEON intrinsics.

## gf.cpp
- [x] additive inverse
- [ ] multiplicative inverse
- [x] additive table
- [x] multiplicative table

The implementation doesn't use intrinsics, instead it uses dynamic array for the flexibility of creating field, which is not as optimized as the intrinsics approach. However, due to the complexity of multiplicative inverse in finite field (especially gf), it is not implemented.

multiplicative inverse reference:
1. [https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/](https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/)
2. [https://en.wikipedia.org/wiki/Finite_field_arithmetic#Multiplicative_inverse](https://en.wikipedia.org/wiki/Finite_field_arithmetic#Multiplicative_inverse)
3. [https://www.youtube.com/watch?v=AcJEyvzUT64&ab_channel=Ch-13ComputerScienceandEngineering](https://www.youtube.com/watch?v=AcJEyvzUT64&ab_channel=Ch-13ComputerScienceandEngineering)

## gf_intrinsics.cpp
The implementation uses a few intrinsics, experimented in hello_world.cpp.

## hello_world.cpp
This file experimented the intrinsics function, with some simple manipulation of intrinsics vector operations.
