#include <iostream>
#include <arm_neon.h>

int main()
{
    int8_t vec[8];
    int8_t sth[8] = {1,2,3,4,5,6,7,8};
    int8_t sth2[8] = {2};
    int8x8_t a, b, c;
    a = vld1_s8(sth);
    b = vld1_s8(sth2);
    c = vadd_s8(a, b);
    vst1_s8(vec, c);
    for (int i = 0; i < 8; ++i) {
        std::cout << int(vec[i]) << std::endl;
    }
    std::cout << "Hello World\n";
}
