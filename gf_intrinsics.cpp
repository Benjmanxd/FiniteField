#include <iostream>
#include <arm_neon.h>
#include <iomanip>

int16_t PP [8] = { 1, 2, 7, 11, 19, 37, 67, 131 };

class GFLookup {
    private:
        int size;
        int8x8_t addInv;
        int8x8_t mulInv;
        uint16x8_t add[8];
        int8x8_t mul[8];

        void ConstructAddTable() {
            uint16x8_t xorVec, single, sequence = { 0, 1, 2, 3, 4 ,5, 6, 7 };
            for (unsigned short i = 0; i < size; ++i) {
                single = vld1q_dup_u16(&i);
                xorVec = veorq_u16(single, sequence);
                vst1q_u16((uint16_t*)&add[i], xorVec);
            }
        }
        void ConstructMulTable() {
            for (int i = 0; i < size; ++i) {

            }
        }

    public:
        GFLookup (int size) : size(size) {
            ConstructAddTable();
            ConstructMulTable();
        }
        void printTable() {
            std::cout << "=== Print Additive Table ===" << std::endl;
            std::cout << std::setw(4) << "x |";
            for (int i = 0; i < size; ++i)
                std::cout << std::setw(7) << i  << "|";
            std::cout << std::endl;
            for (unsigned short i = 0; i < size; ++i) {
                unsigned short *ptr = (unsigned short*)&add[i];
                std::cout << std::setw(2) << i << " |";
                for (unsigned short j = 0; j < size; ++j)
                    std::cout << std::setw(7) << std::bitset<4>(ptr[j]) << "|";
                std::cout << std::endl;
            }
        }
};

int main() {
    GFLookup gf(8);
    gf.printTable();
}
