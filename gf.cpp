#include <iostream>
#include <cmath>
#include <iomanip>

#define GF 2

const int IND = 4;
const int PP [] = { 1, 2, 7, 11, 19, 37, 67, 131, 285 };

enum GF_Size {
    BASE = 1,
    EXT2,
    EXT3,
    EXT4
};

enum table {
    ADD_MUL_INV,
    ADD,
    MUL,
    ALL
};

class GFLookup {
    private: 
        unsigned int order;
        unsigned int size;
        unsigned int *addInv;
        unsigned int *mulInv;
        unsigned int **add;
        unsigned int **mul;

        // === Construct add, addInv, mul, mulInv table ===
        void ConstructAddMulInvTable() {
            addInv = new unsigned int[size];
            mulInv = new unsigned int[size];
            for (int i = 0; i < size; ++i) {
                addInv[i] = AddInverse(i);
                mulInv[i] = MulInverse(i);
            }
        }
        void ConstructAddTable() {
            add = new unsigned int*[size];
            for (int i = 0; i < size; ++i) {
                add[i] = new unsigned int[size];
                for (int j = 0; j < size; ++j)
                    add[i][j] = i ^ j;
            }
        }
        void ConstructMulTable() {
            mul = new unsigned int*[size];
            int f = PP[order];
            for (int i = 0; i < size; ++i) {
                mul[i] = new unsigned int[size];
                for (int j = 0; j < size; ++j) {
                    if (i == 0 || j == 0) mul[i][j] = 0;
                    else mul[i][j] = FindMultiple(i, j);
                }
            }
        }
        // helper function to find additive inverse
        unsigned int AddInverse(unsigned int a) const {
            if (a == 0) return 0;
            return size - a;
        }
        // helper function to find multiplicative inverse
        // reference1: https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/
        // reference2: https://en.wikipedia.org/wiki/Finite_field_arithmetic#Multiplicative_inverse
        // reference3: https://www.youtube.com/watch?v=AcJEyvzUT64&ab_channel=Ch-13ComputerScienceandEngineering
        // not implemented - the MulInverse is wrong
        unsigned int MulInverse(unsigned int a) const {
            if (a == 0) return -1;
            unsigned int r = size - 1;
            unsigned int aRm1 = ((int)std::pow(a, r-1)) % size;
            unsigned int aR = aRm1 * a;
            unsigned int aRInv = MulInverse(aR, PP[order]);
            unsigned int aInv = aRInv * aRm1;
            return aInv;
        }
        // helper function to find GF multiple
        // reference: https://en.wikipedia.org/wiki/Finite_field_arithmetic#C_programming_example
        unsigned int FindMultiple(unsigned int a, unsigned int b) const {
            unsigned int p = 0;
            unsigned int m = PP[order];
            while (a != 0 && b != 0) {
                if (b & 1) p ^= a;
                if (a & 0x80) a = (a << 1) ^ m;
                else a <<= 1;
                b >>= 1;
            }
            return p;
        }

        // === print function ===
        void PrintAddMulInvTable() const {
            std::cout << "=== Print Additive / Multiplicative Inverse Table ===" << std::endl;
            std::cout << std::setw(5) << "x |" << std::setw(7) << "-x" << "|" << std::setw(7) << "x^-1" << "|" << std::endl;
            for (int i = 0; i < size; ++i) {
                std::cout << std::setw(3) << i << " |" << std::setw(7) << addInv[i] << "|" << std::setw(7);//std::bitset<IND>(addInv[i]) << "|" << std::setw(7);
                if (mulInv[i] == -1) std::cout << "UND" << "|" << std::endl;
                else std::cout << mulInv[i] << "|" << std::endl;//std::bitset<IND>(mulInv[i]) << std::endl;
            }
        }
        void PrintAddTable() const {
            std::cout << "=== Print Additive Table ===" << std::endl;
            std::cout << std::setw(5) << "x |";
            for (int i = 0; i < size; ++i)
                std::cout << std::setw(7) << i  << "|";
            std::cout << std::endl;
            for (int i = 0; i < size; ++i) {
                std::cout << std::setw(3) << i << " |";
                for (int j = 0; j < size; ++j)
                    std::cout << std::setw(7) << add[i][j] << "|";//std::bitset<IND>(add[i][j]) << "|";
                std::cout << std::endl;
            }
        }
        void PrintMulTable() const {
            std::cout << "=== Print Multiplicative Table ===" << std::endl;
            std::cout << std::setw(5) << "x |";
            for (int i = 0; i < size; ++i)
                std::cout << std::setw(7) << i  << "|";
            std::cout << std::endl;
            for (int i = 0; i < size; ++i) {
                std::cout << std::setw(3) << i << " |";
                for (int j = 0; j < size; ++j)
                    std::cout << std::setw(7) << mul[i][j] << "|";//std::bitset<IND>(mul[i][j]) << "|";
                std::cout << std::endl;
            }
        } 

    public:
        GFLookup(int s) : order(s), size(std::pow(GF, s)) {
            ConstructAddMulInvTable();
            ConstructAddTable();
            ConstructMulTable();
        }

        ~GFLookup() {
            delete [] addInv;
            delete [] mulInv;
            for (int i = 0; i < size; ++i) {
                delete add[i];
                delete mul[i];
            }
            delete [] add;
            delete [] mul;
        }

        void printTable(table x = ALL) const {
            if (x == ALL) {
                PrintAddMulInvTable();
                PrintAddTable();
                PrintMulTable();
            } else if (x == ADD_MUL_INV) 
                PrintAddMulInvTable();
            else if (x == ADD)
                PrintAddTable();
            else if (x == MUL)
                PrintMulTable();
        }  
        static unsigned int MulInverse(unsigned int a, unsigned int m) {
            if (a == 0)
                return -1;
            unsigned int m0 = m, y = 0, x = 1, n = 0;
            while (a > 1) {
                int q = a / m;
                int t = m;
                m = a % m, a = t;
                t = y;
                y = x - q * y;
                x = t;
                ++ n;
            }
            if (x < 0) x += m0;
            return x;
        }
};

int main() {
    //std::cout << "=== === GF Base Field === ===" << std::endl;
    //GFLookup GFBase(BASE);
    //GFBase.printTable();
    //std::cout << std::endl;
    //std::cout << "=== === GF Ext Field X2 === ===" << std::endl;
    //GFLookup GFX2(EXT2);
    //GFX2.printTable();
    //std::cout << std::endl;
    //std::cout << "=== === GF Ext Field X3 === ===" << std::endl;
    //GFLookup GFX3(EXT3);
    //GFX3.printTable();
    //std::cout << std::endl;
    //std::cout << "=== === GF Ext Field X4 === ===" << std::endl;
    //GFLookup GFX4(EXT4);
    //GFX4.printTable();
    std::cout << "=== === GF Ext Field X8 === ===" << std::endl;
    GFLookup GFX8(8);
    GFX8.printTable();
}
