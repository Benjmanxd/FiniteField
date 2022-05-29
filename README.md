# Finite Field - n-hop part-time

### the intrinsics used in FiniteField project are ARM NEON intrinsics.

## Basic finite field concepts (Galois field)
- base field (0,1), and extended field
- [primitive polynomial](https://en.wikipedia.org/wiki/Finite_field_arithmetic#Primitive_polynomials) -> for each extended fields, [checkout](https://www.partow.net/programming/polynomials/index.html)
- additive inverse -> mod calculation
- multiplicative inverse -> [Itoh–Tsujii inversion algorithm](https://en.wikipedia.org/wiki/Itoh%E2%80%93Tsujii_inversion_algorithm)
- [addition](https://en.wikipedia.org/wiki/Finite_field_arithmetic#C_programming_example) -> xor operation
- [multiplication](https://en.wikipedia.org/wiki/Finite_field_arithmetic#C_programming_example) -> Rijndael algorithm, [Russian peasant multiplication algorithm](https://en.wikipedia.org/wiki/Ancient_Egyptian_multiplication#Russian_peasant_multiplication)

## Intrinsics
- very optimized code that is architecture dependent
- [ARM NEON Intrinsics](https://developer.arm.com/architectures/instruction-sets/intrinsics/)
- [Intel Intrinsics](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#)
    - [tutorial](https://www.youtube.com/watch?v=x9Scb5Mku1g)
    - [code example 1](https://github.com/yrp604/arm64_intrinsics)
    - [code example 2](https://github.com/rogerou/Arm-neon-intrinsics)
    - [code example 3](https://github.com/intel/ARM_NEON_2_x86_SSE)
    - [code example 4](https://github.com/moepinet/libmoepgf)

## Finite Field further reading
- [Finite Field guide in cryptography and proofs](https://www.youtube.com/watch?v=ColSUxhpn6A)
- [MIT math notes](https://math.mit.edu/classes/18.783/2017/LectureNotes4.pdf)
- [Chinese Finite and Galois Field explanation](https://web.ntnu.edu.tw/~algo/FiniteField.html)
- [Finite Field Arithmetic notes](https://www.uotechnology.edu.iq/dep-eee/lectures/4th/Communication/Information%20theory/8.pdf)

## File explanations

### gf.cpp
- [x] additive inverse
- [ ] multiplicative inverse
- [x] additive table
- [x] multiplicative table

The implementation doesn't use intrinsics, instead it uses dynamic array for the flexibility of creating field, which is not as optimized as the intrinsics approach. However, due to the complexity of multiplicative inverse in finite field (especially gf), it is not implemented.

multiplicative inverse reference:
1. [reference1](https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/)
2. [reference2](https://en.wikipedia.org/wiki/Finite_field_arithmetic#Multiplicative_inverse)
3. [reference3](https://www.youtube.com/watch?v=AcJEyvzUT64&ab_channel=Ch-13ComputerScienceandEngineering)

### gf_intrinsics.cpp
The implementation uses a few intrinsics, experimented in hello_world.cpp.

### hello_world.cpp
This file experimented the intrinsics function, with some simple manipulation of intrinsics vector operations.

### gf8.txt
The gf(8) table generated using gf.cpp, only the multiplicative inverse field is incorrect.
