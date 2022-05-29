#ifndef FINITEFIELD_HPP
#define FINITEFIELD_HPP

#include <string.h>
#include <iostream>
#include <iomanip>

#include "config.hpp"
#include "divtable_order8.hpp"
#include "multable_order8.hpp"


class FiniteField {
   public:
    // Note: All 2-d matrix are formed by column vectors unless otherwise
    // specified
    static SymbolType** CreateMat(int col, int row);

    static void DelMat(SymbolType** mat, int col);

    static inline SymbolType MulElem(SymbolType x, SymbolType y) {
        return mul_table_o8[(x << FIELD_ORDER) | y];
    }

    static inline SymbolType DivElem(SymbolType x, SymbolType y) {
        return div_table_o8[(x << FIELD_ORDER) | y];
    }

    static inline SymbolType InvElem(SymbolType x) {
        return inv_table_o8[x];
    }

    // return x dot y
    //NOTE: Here length is always batch_size, so we will hardcode the length to 4 or 16 in the SIMD implementation
    static SymbolType InnerProd(SymbolType* x, SymbolType* y, int length);

    // z[] = x[][] * y[]
    static void MulMat(SymbolType* z,
                       SymbolType** x,
                       SymbolType* y,
                       int col,
                       int row);

    // z[] += x[][] * y[]
    static void AddMulMat(SymbolType* z,
                          SymbolType** x,
                          SymbolType* y,
                          int col,
                          int row);

    // x[] = c * x[]
    static void MulVec(SymbolType* x, SymbolType c, int size);

    // z[] += c * x[]
    static void AddMulVec(SymbolType* z, SymbolType* x, SymbolType c, int size);

    // z[] += x[]
    static void AddVec(SymbolType* z, SymbolType* x, int size);

    // z[] -= c * x[]
    static void SubMulVec(SymbolType* z, SymbolType* x, SymbolType c, int size);

    // reduce the matrix to pseudo row-echelon form and return its rank
    static int GetRank(SymbolType** mat, int col, int row);

    // set inverse matrix [A, I] --> [I, A^-1], return rank
    static int GaussianElimination(SymbolType** mat,
                                   SymbolType** inv_mat,
                                   int col,
                                   int row);

    // Solve X for X * A = Y
    // X = Y * inv(A)
    static int GaussianSolve(SymbolType** A,
                             SymbolType** Y,
                             SymbolType** X,
                             int col_a,
                             int row_a,
                             int row_y);

    static void Print(SymbolType* mat, int size, char mode);

    static void Print(SymbolType** mat, int col, int row);
};

#endif
