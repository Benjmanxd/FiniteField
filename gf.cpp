#include <iostream>
#include <cmath>
#include <iomanip>

#define GF 2

const int IND = 4;
const int PP [] = { 1, 2, 7, 11, 19 };

enum size {
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
        int order;
        int size;
        int *addInv;
        int *mulInv;
        int **add;
        int **mul;

        // === Construct add, addInv, mul, mulInv table ===
        void ConstructAddMulInvTable() {
            addInv = new int[size];
            mulInv = new int[size];
            for (int i = 0; i < size; ++i) {
                addInv[i] = AddInverse(i);
                mulInv[i] = MulInverse(i);
            }
        }
        void ConstructAddTable() {
            add = new int*[size];
            for (int i = 0; i < size; ++i) {
                add[i] = new int[size];
                for (int j = 0; j < size; ++j)
                    add[i][j] = i ^ j;
            }
        }
        void ConstructMulTable() {
            mul = new int*[size];
            int f = PP[order];
            for (int i = 0; i < size; ++i) {
                mul[i] = new int[size];
                for (int j = 0; j < size; ++j) {
                    if (i == 0 || j == 0) mul[i][j] = 0;
                    else mul[i][j] = FindMultiple(i, j);
                }
            }
        }
        // helper function to find additive inverse
        int AddInverse(int a) const {
            if (a == 0)
                return 0;
            return size - a;
        }
        // helper function to find multiplicative inverse
        int MulInverse(int a) const {
            if (a == 0)
                return -1;
            int m = PP[order], m0 = PP[order], y = 0, x = 1, n = 0;
            while (a > 1) {
                int q = a / m;
                int t = m;
                m = a % m, a = t;
                t = y;
                y = x - q * y;
                x = t;
                ++ n;
            }
            if (x < 0)
                x += m0;
            return x;
        }
        // helper function to find GF multiple
        int FindMultiple(const int& i, const int& j) const {
            int f = PP[order];
            int res = i * j;
            int pow = order * 2;
            while (res >= std::pow(GF, order) && pow >= order) {
                if (res >= std::pow(GF, pow))
                    res = res ^ (f << (pow - order));
                -- pow;
            }
            return res;       
        }

        // === print function ===
        void PrintAddMulInvTable() const {
            std::cout << "=== Print Additive / Multiplicative Inverse Table ===" << std::endl;
            std::cout << std::setw(4) << "x |" << std::setw(7) << "-x" << "|" << std::setw(7) << "x^-1" << "|" << std::endl;
            for (int i = 0; i < size; ++i) {
                std::cout << std::setw(2) << i << " |" << std::setw(7) << std::bitset<IND>(addInv[i]) << "|" << std::setw(7);
                if (mulInv[i] == -1) std::cout << "UND" << std::endl;
                else std::cout << std::bitset<IND>(mulInv[i]) << std::endl;
            }
        }
        void PrintAddTable() const {
            std::cout << "=== Print Additive Table ===" << std::endl;
            std::cout << std::setw(4) << "x |";
            for (int i = 0; i < size; ++i)
                std::cout << std::setw(7) << i  << "|";
            std::cout << std::endl;
            for (int i = 0; i < size; ++i) {
                std::cout << std::setw(2) << i << " |";
                for (int j = 0; j < size; ++j)
                    std::cout << std::setw(7) << std::bitset<IND>(add[i][j]) << "|";
                std::cout << std::endl;
            }
        }
        void PrintMulTable() const {
            std::cout << "=== Print Multiplicative Table ===" << std::endl;
            std::cout << std::setw(4) << "x |";
            for (int i = 0; i < size; ++i)
                std::cout << std::setw(7) << i  << "|";
            std::cout << std::endl;
            for (int i = 0; i < size; ++i) {
                std::cout << std::setw(2) << i << " |";
                for (int j = 0; j < size; ++j)
                    std::cout << std::setw(7) << std::bitset<IND>(mul[i][j]) << "|";
                std::cout << std::endl;
            }
        } 

    public:
        GFLookup(int size) : order(size), size(std::pow(GF, size)) {
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
};

int main() {
    std::cout << "=== === GF Base Field === ===" << std::endl;
    GFLookup GFBase(BASE);
    GFBase.printTable();
    std::cout << std::endl;
    std::cout << "=== === GF Ext Field X2 === ===" << std::endl;
    GFLookup GFX2(EXT2);
    GFX2.printTable();
    std::cout << std::endl;
    std::cout << "=== === GF Ext Field X3 === ===" << std::endl;
    GFLookup GFX3(EXT3);
    GFX3.printTable();
    std::cout << std::endl;
    std::cout << "=== === GF Ext Field X4 === ===" << std::endl;
    GFLookup GFX4(EXT4);
    GFX4.printTable();
}
