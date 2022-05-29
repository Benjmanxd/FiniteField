#include "finitefield.hpp"
#include "config.hpp"

#include <algorithm> // min
#include <utility>   //swap (This function is no longer defined in header <algorithm>, but in <utility>)

#define USE_SPLIT_MUL  //Use this on systems with small cache

using namespace std;

//===================================================================================================
SymbolType** FiniteField::CreateMat(int col, int row) {
    row = PAD32(row);
    
    // SymbolType* pool =  new SymbolType[col * row];
    // memset(pool, 0, col * row * sizeof(SymbolType));
    SymbolType* pool = (SymbolType*)calloc (col * row,sizeof(SymbolType));
    
    SymbolType** mat = new SymbolType*[col];
    for (int i = 0; i < col; i++, pool += row) {
        mat[i] = pool;
    }
    return mat;
}

//===================================================================================================
void FiniteField::DelMat(SymbolType** mat, int col) {
    free(mat[0]);
    delete[] mat;
}

//===================================================================================================

int FiniteField::GetRank(SymbolType** mat, int col, int row) {
    int max_rank = min(row, col); // used from <algorithm>
    for (int i = 0, r = 0; i < max_rank; i++) {
        if (mat[i][r] == 0) {
            // look for non zero column in r-th row and swap
            int non_zero_col;
            for (; r < row; r++) {
                for (non_zero_col = i; non_zero_col < col; non_zero_col++) {
                    if (mat[non_zero_col][r]) {
                        goto FOUND_NON_ZERO_COLUMN;
                    }
                }
            }
            return i; // if non zero column not found 

        FOUND_NON_ZERO_COLUMN:
            // swap(mat[non_zero_col], mat[i]);
            if (i != non_zero_col) {
                // AddMulVec(mat[i], mat[non_zero_col], 1, row);
                AddVec(mat[i], mat[non_zero_col], row);
            }
        }

        {
            SymbolType c = InvElem(mat[i][r]);
            // SymbolType c = DivElem(1, mat[i][r]);

            // Simplified from MulVec(mat[i], c, row);
            MulVec(mat[i], c, row);
            // MulVec(mat[i] + r, c, row - r);
        }

        // reduce the i-th column to 0
        for (int k = i + 1; k < col; k++) {
            if (mat[k][r] != 0) {
                SymbolType c = mat[k][r];
                SubMulVec(mat[k],mat[i],c,row); // subtract row_k = row_k- c *row_i

                // for (int j = r; j < row; j++) {
                //     // subract element
                //     mat[k][j] ^= MulElem(mat[i][j], c);
                // }
            }
        }
    }
    return max_rank;
}

//===================================================================================================
int FiniteField::GaussianElimination(SymbolType** mat,
                                     SymbolType** inv_mat,
                                     int col,
                                     int row) {
    if (col < row) {
        printf("GaussianElimination: not supported\n");
        return 0;
    }

    for (int i = 0; i < row; i++) {
        if (mat[i][i] == 0) {
            // look for non zero column in i-th row and swap
            int non_zero_col;
            for (non_zero_col = i; non_zero_col < col; non_zero_col++) {
                if (mat[non_zero_col][i]) {
                    break;
                }
            }

            if (non_zero_col == col) {
                return i;
            }

            // swap(mat[non_zero_col], mat[i]);
            // swap(inv_mat[non_zero_col], inv_mat[i]);

            // AddMulVec(mat[i],mat[non_zero_col], 1, row);
            // AddMulVec(inv_mat[i],inv_mat[non_zero_col],  1, row);

            AddVec(mat[i],mat[non_zero_col],row);
            AddVec(inv_mat[i],inv_mat[non_zero_col],row);
        }

        {
            SymbolType c = InvElem(mat[i][i]);
            // SymbolType c = DivElem(1, mat[i][i]);

            
            MulVec(mat[i], c, row);
            // MulVec(mat[i] + i, c, row - i); // Simplified from MulVec(mat[i], c, row);

            MulVec(inv_mat[i], c, col);
        }

        // reduce the i-th column to 0
        for (int k = 0; k < col; k++) {
            if (k == i || mat[k][i] == 0) {
                continue;
            }
            SymbolType c = mat[k][i];
            SubMulVec(mat[k],mat[i],c,row);
            SubMulVec(inv_mat[k], inv_mat[i], c, col);

            // for (int j = i; j < row; j++) {
            //     // subract element
            //     mat[k][j] ^= MulElem(mat[i][j], c);
            // }
                    }
    }

    return row;
}

//===================================================================================================
int FiniteField::GaussianSolve(SymbolType** A,
                               SymbolType** Y,
                               SymbolType** X,
                               int col_a,
                               int row_a,
                               int row_y) {
    if (col_a < row_a) {
        printf("GaussianSolve: not supported\n");
        return 0;
    }

    SymbolType** inv = CreateMat(col_a, row_a);
    for (int i = 0; i < row_a; i++) {
        inv[i][i] = 1;
    }

    int rank = GaussianElimination(A, inv, col_a, row_a);
    if (rank == min(row_a, col_a)) {

        for (int i = 0; i < row_a; i++) {
            MulMat(X[i], Y, inv[i], row_a, row_y);
        }
    }

    DelMat(inv, col_a);

    return rank;
}


//===================================================================================================
void FiniteField::MulMat(SymbolType* z,
                         SymbolType** x,
                         SymbolType* y,
                         int col,
                         int row) {
    memset(z, 0, sizeof(SymbolType) * row);
    AddMulMat(z, x, y, col, row);
}

//===================================================================================================
void FiniteField::AddMulMat(SymbolType* z,
                            SymbolType** x,
                            SymbolType* y,
                            int col,
                            int row) {
    for (int i = 0; i < col; i++) {
        AddMulVec(z, x[i], y[i], row);
    }
}

//===================================================================================================
void FiniteField::SubMulVec(SymbolType* z,
                            SymbolType* x,
                            SymbolType c,
                            int size) {
    AddMulVec(z, x, c, size);
}

//===================================================================================================
//__NOT_USE_SIMD__NOT_USE_SIMD__NOT_USE_SIMD__NOT_USE_SIMD__NOT_USE_SIMD__NOT_USE_SIMD__NOT_USE_SIMD
#ifndef USE_SIMD
//===================================================================================================
SymbolType FiniteField::InnerProd(SymbolType* x, SymbolType* y, int size) {
    SymbolType result = 0;

    for (int i = 0; i < size; i++) {
        result ^= MulElem(x[i], y[i]);
    }
    return result;
}
//===================================================================================================
void FiniteField::AddMulVec(SymbolType* z,
                            SymbolType* x,
                            SymbolType c,
                            int size) {
    //size = PAD32(size);
    #ifdef USE_SPLIT_MUL
    const SymbolType* const mul_table_high = mul_table_o8_high_bits + (c << (FIELD_ORDER>>1));
    const SymbolType* const mul_table_low  = mul_table_o8_low_bits  + (c << (FIELD_ORDER>>1));
    #else
    const SymbolType* const mul_table = mul_table_o8 + (c << (FIELD_ORDER));
    #endif
    for (int i = 0; i < (size >> 4); ++i) {

        #ifdef USE_SPLIT_MUL
        *(z+0)  ^= mul_table_high[*(x+0) >> 4 ] ^ mul_table_low[*(x+0)  & 15];
        *(z+1)  ^= mul_table_high[*(x+1) >> 4 ] ^ mul_table_low[*(x+1)  & 15];
        *(z+2)  ^= mul_table_high[*(x+2) >> 4 ] ^ mul_table_low[*(x+2)  & 15];
        *(z+3)  ^= mul_table_high[*(x+3) >> 4 ] ^ mul_table_low[*(x+3)  & 15];
        *(z+4)  ^= mul_table_high[*(x+4) >> 4 ] ^ mul_table_low[*(x+4)  & 15];
        *(z+5)  ^= mul_table_high[*(x+5) >> 4 ] ^ mul_table_low[*(x+5)  & 15];
        *(z+6)  ^= mul_table_high[*(x+6) >> 4 ] ^ mul_table_low[*(x+6)  & 15];
        *(z+7)  ^= mul_table_high[*(x+7) >> 4 ] ^ mul_table_low[*(x+7)  & 15];

        *(z+8)  ^= mul_table_high[*(x+8) >> 4 ] ^ mul_table_low[*(x+8)  & 15];
        *(z+9)  ^= mul_table_high[*(x+9) >> 4 ] ^ mul_table_low[*(x+9)  & 15];
        *(z+10) ^= mul_table_high[*(x+10) >> 4] ^ mul_table_low[*(x+10) & 15];
        *(z+11) ^= mul_table_high[*(x+11) >> 4] ^ mul_table_low[*(x+11) & 15];
        *(z+12) ^= mul_table_high[*(x+12) >> 4] ^ mul_table_low[*(x+12) & 15];
        *(z+13) ^= mul_table_high[*(x+13) >> 4] ^ mul_table_low[*(x+13) & 15];
        *(z+14) ^= mul_table_high[*(x+14) >> 4] ^ mul_table_low[*(x+14) & 15];
        *(z+15) ^= mul_table_high[*(x+15) >> 4] ^ mul_table_low[*(x+15) & 15];
        
        #else
        *(z+0)  ^= mul_table[*(x+0)];
        *(z+1)  ^= mul_table[*(x+1)];
        *(z+2)  ^= mul_table[*(x+2)];
        *(z+3)  ^= mul_table[*(x+3)];
        *(z+4)  ^= mul_table[*(x+4)];
        *(z+5)  ^= mul_table[*(x+5)];
        *(z+6)  ^= mul_table[*(x+6)];
        *(z+7)  ^= mul_table[*(x+7)];

        *(z+8)  ^= mul_table[*(x+8)];
        *(z+9)  ^= mul_table[*(x+9)];
        *(z+10) ^= mul_table[*(x+10)];
        *(z+11) ^= mul_table[*(x+11)];
        *(z+12) ^= mul_table[*(x+12)];
        *(z+13) ^= mul_table[*(x+13)];
        *(z+14) ^= mul_table[*(x+14)];
        *(z+15) ^= mul_table[*(x+15)];
        #endif
        z += 16;
        x += 16;
    }

    for (int i = 0; i < (size & 15); ++i) {
        #ifdef USE_SPLIT_MUL
        *z ^= mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];
        z++;x++;
        #else        
        *z ^= mul_table[*x];
        z++;x++;
        #endif
    }
}

//===================================================================================================
void FiniteField::MulVec(SymbolType* x, SymbolType c, int size) {
    //size = PAD32(size);
    #ifdef USE_SPLIT_MUL
    const SymbolType* const mul_table_high = mul_table_o8_high_bits + (c << (FIELD_ORDER>>1));
    const SymbolType* const mul_table_low  = mul_table_o8_low_bits  + (c << (FIELD_ORDER>>1));
    #else
    const SymbolType* const mul_table = mul_table_o8 + (c << (FIELD_ORDER));
    #endif
    
    for (int i = 0; i < (size >> 4); ++i) {

        #ifdef USE_SPLIT_MUL
        *(x+0)  = mul_table_high[*(x+0) >> 4 ] ^ mul_table_low[*(x+0)  & 15];
        *(x+1)  = mul_table_high[*(x+1) >> 4 ] ^ mul_table_low[*(x+1)  & 15];
        *(x+2)  = mul_table_high[*(x+2) >> 4 ] ^ mul_table_low[*(x+2)  & 15];
        *(x+3)  = mul_table_high[*(x+3) >> 4 ] ^ mul_table_low[*(x+3)  & 15];
        *(x+4)  = mul_table_high[*(x+4) >> 4 ] ^ mul_table_low[*(x+4)  & 15];
        *(x+5)  = mul_table_high[*(x+5) >> 4 ] ^ mul_table_low[*(x+5)  & 15];
        *(x+6)  = mul_table_high[*(x+6) >> 4 ] ^ mul_table_low[*(x+6)  & 15];
        *(x+7)  = mul_table_high[*(x+7) >> 4 ] ^ mul_table_low[*(x+7)  & 15];

        *(x+8)  = mul_table_high[*(x+8) >> 4 ] ^ mul_table_low[*(x+8)  & 15];
        *(x+9)  = mul_table_high[*(x+9) >> 4 ] ^ mul_table_low[*(x+9)  & 15];
        *(x+10) = mul_table_high[*(x+10) >> 4] ^ mul_table_low[*(x+10) & 15];
        *(x+11) = mul_table_high[*(x+11) >> 4] ^ mul_table_low[*(x+11) & 15];
        *(x+12) = mul_table_high[*(x+12) >> 4] ^ mul_table_low[*(x+12) & 15];
        *(x+13) = mul_table_high[*(x+13) >> 4] ^ mul_table_low[*(x+13) & 15];
        *(x+14) = mul_table_high[*(x+14) >> 4] ^ mul_table_low[*(x+14) & 15];
        *(x+15) = mul_table_high[*(x+15) >> 4] ^ mul_table_low[*(x+15) & 15];
        #else
        *(x+0) = mul_table[*(x+0)];
        *(x+1) = mul_table[*(x+1)];
        *(x+2) = mul_table[*(x+2)];
        *(x+3) = mul_table[*(x+3)];
        *(x+4) = mul_table[*(x+4)];
        *(x+5) = mul_table[*(x+5)];
        *(x+6) = mul_table[*(x+6)];
        *(x+7) = mul_table[*(x+7)];

        *(x+8) = mul_table[*(x+8)];
        *(x+9) = mul_table[*(x+9)];
        *(x+10) = mul_table[*(x+10)];
        *(x+11) = mul_table[*(x+11)];
        *(x+12) = mul_table[*(x+12)];
        *(x+13) = mul_table[*(x+13)];
        *(x+14) = mul_table[*(x+14)];
        *(x+15) = mul_table[*(x+15)];

        #endif
        x += 16;
    }

    for (int i = 0; i < (size & 15); ++i) {
        #ifdef USE_SPLIT_MUL
        *x = mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];
        x++;
        #else
        *x = mul_table[*x];
        x++;
        #endif
    }
}


//===================================================================================================
void FiniteField::AddVec(SymbolType* z,
                            SymbolType* x,
                            int size) {

    //size = PAD32(size);
    for (int i = 0; i < (size >> 4); ++i) {

        *(z+0)  ^= *(x+0);
        *(z+1)  ^= *(x+1);
        *(z+2)  ^= *(x+2);
        *(z+3)  ^= *(x+3);
        *(z+4)  ^= *(x+4);
        *(z+5)  ^= *(x+5);
        *(z+6)  ^= *(x+6);
        *(z+7)  ^= *(x+7);

        *(z+8)  ^= *(x+8);
        *(z+9)  ^= *(x+9);
        *(z+10) ^= *(x+10);
        *(z+11) ^= *(x+11);
        *(z+12) ^= *(x+12);
        *(z+13) ^= *(x+13);
        *(z+14) ^= *(x+14);
        *(z+15) ^= *(x+15);
        z += 16;
        x += 16;
    }

    for (int i = 0; i < (size & 15); ++i) { 
        *z ^= *x;
        z++;x++;
    }
}


//===================================================================================================
//__USE_SIMD__USE_SIMD__USE_SIMD__USE_SIMD__USE_SIMD__USE_SIMD__USE_SIMD__USE_SIMD__USE_SIMD__USE_SIMD
#else  // if defined USE_SIMD
//===================================================================================================
// __AVX2____AVX2____AVX2____AVX2____AVX2____AVX2____AVX2____AVX2____AVX2____AVX2____AVX2____AVX2__
#if defined(__AVX2__)
#include <immintrin.h>

// The intrinsic _mm256_loadu2_m128i requires gcc version >= 10
#ifndef  _mm256_loadu2_m128i
#define _mm256_loadu2_m128i(x,y) (_mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128(x)),_mm_loadu_si128(y),1))
#endif
//===================================================================================================
void FiniteField::AddMulVec(SymbolType* z, SymbolType* x,
                            SymbolType c , int size)
{
    const SymbolType* const mul_table_high = mul_table_o8_high_bits + (c << (FIELD_ORDER>>1));
    const SymbolType* const mul_table_low  = mul_table_o8_low_bits  + (c << (FIELD_ORDER>>1));

    const __m256i v_mask_low = _mm256_set1_epi8(15);
    const __m256i v_map_high = _mm256_loadu2_m128i((const __m128i*)mul_table_high,(const __m128i*)mul_table_high);
    const __m256i v_map_low  = _mm256_loadu2_m128i((const __m128i*)mul_table_low ,(const __m128i*)mul_table_low );

    __m256i v_x;
    __m256i v_x_high;
    __m256i v_x_low;
    __m256i v_product;

    for (int i = 0; i < ( size >> 5 ); i++) { // Here we assume PAD32
        v_x       = _mm256_loadu_si256((const __m256i*)x);
        v_x_high  = _mm256_and_si256(_mm256_srli_epi16(v_x, 4),v_mask_low);  // right shift 4
        v_x_low   = _mm256_and_si256(v_x, v_mask_low);
        v_product = _mm256_xor_si256(_mm256_shuffle_epi8(v_map_high, v_x_high),
                                          _mm256_shuffle_epi8(v_map_low, v_x_low));

        _mm256_storeu_si256((__m256i*)z,
                             _mm256_xor_si256(v_product,
                                              _mm256_loadu_si256((const __m256i*)(z) ) 
                                             )
                            );

        z += 32;
        x += 32;
    }
    for (int i = 0; i < (size & 31); ++i) {
        *z ^= mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];
        z++;x++;
    }
}
//===================================================================================================
void FiniteField::MulVec(SymbolType* x, SymbolType c, int size) 
{
    const SymbolType* const mul_table_high = mul_table_o8_high_bits + (c << (FIELD_ORDER>>1));
    const SymbolType* const mul_table_low  = mul_table_o8_low_bits  + (c << (FIELD_ORDER>>1));

    const __m256i v_mask_low = _mm256_set1_epi8(15);
    const __m256i v_map_high = _mm256_loadu2_m128i((const __m128i*)mul_table_high,(const __m128i*)mul_table_high);
    const __m256i v_map_low  = _mm256_loadu2_m128i((const __m128i*)mul_table_low ,(const __m128i*)mul_table_low );
 
    __m256i v_x;
    __m256i v_x_high;
    __m256i v_x_low;
    __m256i v_product;
    
    for (int i = 0; i < ( size >> 5 ); i++) { // Here we assume PAD32
        v_x       = _mm256_loadu_si256((const __m256i*)x);
        v_x_high  = _mm256_and_si256(_mm256_srli_epi16(v_x, 4),v_mask_low);  // right shift 4
        v_x_low   = _mm256_and_si256(v_x, v_mask_low);
        v_product = _mm256_xor_si256(_mm256_shuffle_epi8(v_map_high, v_x_high),
                                          _mm256_shuffle_epi8(v_map_low, v_x_low));

        _mm256_storeu_si256((__m256i*)x,v_product);

        x += 32;
    }
    for (int i = 0; i < (size & 31); ++i) {
        *x = mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];
        x++;
    }
}
//===================================================================================================
void FiniteField::AddVec(SymbolType* z, SymbolType* x, int size)
{
    __m256i v_x;

    for (int i = 0; i < ( size >> 5 ); i++) { // Here we assume PAD32
        v_x       = _mm256_loadu_si256((const __m256i*)x);
        _mm256_storeu_si256((__m256i*)z,
                             _mm256_xor_si256(v_x,
                                              _mm256_loadu_si256((const __m256i*)(z) ) 
                                             )
                            );

        z += 32;
        x += 32;
    }
    for (int i = 0; i < (size & 31); ++i) { 
        *z ^= *x;
        z++;x++;
    }
}

//===================================================================================================
//__SSE4_2____SSE4_2____SSE4_2____SSE4_2____SSE4_2____SSE4_2____SSE4_2____SSE4_2____SSE4_2____SSE4_2__
#elif defined(__SSE4_2__)
#include <nmmintrin.h>
//===================================================================================================
void FiniteField::AddMulVec(SymbolType* z, SymbolType* x,
                            SymbolType c , int size)
{
    const SymbolType* const mul_table_high = mul_table_o8_high_bits + (c << (FIELD_ORDER>>1));
    const SymbolType* const mul_table_low  = mul_table_o8_low_bits  + (c << (FIELD_ORDER>>1));
    
    const __m128i v_mask_low = _mm_set1_epi8(15);
    const __m128i v_map_low  = _mm_load_si128((const __m128i*)mul_table_low );
    const __m128i v_map_high = _mm_load_si128((const __m128i*)mul_table_high);

    __m128i v_y_low;
    __m128i v_y_high;
    __m128i v_y_map_low;
    __m128i v_y_map_high;

    __m128i v_x_low;
    __m128i v_x_high;
    __m128i v_x_map_low;
    __m128i v_x_map_high;

    for (int i = 0; i < ( size >> 5 ); i++) { // Assume PAD32
        // Equivalent to *z ^= mul_table_high[*x >> 4] ^ mul_table[*x & 15];

        v_y_low  = _mm_loadu_si128( (const __m128i*)x );    // load 16 bytes of x
        v_y_high = _mm_srli_epi16(v_y_low, 4);                          // shift to put high bits
        v_y_map_low = v_map_low;                                        // multiplication table for the low part
        v_y_map_high = v_map_high;                                      // multiplication table for the high part

        v_x_low  = _mm_loadu_si128( (const __m128i*)(x+16) );
        v_x_high = _mm_srli_epi16(v_x_low, 4);  // shift to put high bits
        v_x_map_low  = v_map_low;
        v_x_map_high = v_map_high;

        //apply mask [.... 1111 0000 1111]
        v_y_low =  _mm_and_si128(v_y_low  ,v_mask_low);
        v_y_high = _mm_and_si128(v_y_high ,v_mask_low);
        v_x_low =  _mm_and_si128(v_x_low  ,v_mask_low);
        v_x_high = _mm_and_si128(v_x_high ,v_mask_low);

        //table lookup
        v_y_map_low =  _mm_shuffle_epi8(v_y_map_low , v_y_low );
        v_y_map_high = _mm_shuffle_epi8(v_y_map_high, v_y_high);
        v_x_map_low =  _mm_shuffle_epi8(v_x_map_low , v_x_low );
        v_x_map_high = _mm_shuffle_epi8(v_x_map_high, v_x_high);
        
        //add the result of the multiplication of low and high part
        v_y_map_low = _mm_xor_si128(v_y_map_high, v_y_map_low);
        v_x_map_low = _mm_xor_si128(v_x_map_high, v_x_map_low);

        // add the product to z
        v_y_high = _mm_loadu_si128( reinterpret_cast<__m128i*>(z) );
        v_y_low  = _mm_xor_si128(v_y_map_low, v_y_high);

        //store it into z
        _mm_storeu_si128( reinterpret_cast<__m128i*>(z), v_y_low);

        v_x_high = _mm_loadu_si128( reinterpret_cast<__m128i*>(z+16) );
        v_x_low  = _mm_xor_si128(v_x_map_low, v_x_high);
        
        _mm_storeu_si128( reinterpret_cast<__m128i*>(z+16), v_x_low);

        //move the pointer
        z += 32;
        x += 32;
    }
    for (int i = 0; i < (size & 31); ++i) {
        *z ^= mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];
        z++;x++;
    }
}
//===================================================================================================
void FiniteField::MulVec(SymbolType* x, SymbolType c, int size) 
{
    const SymbolType* const mul_table_high = mul_table_o8_high_bits + (c << (FIELD_ORDER>>1));
    const SymbolType* const mul_table_low  = mul_table_o8_low_bits  + (c << (FIELD_ORDER>>1));

    const __m128i v_mask_low = _mm_set1_epi8(15);
    const __m128i v_map_low  = _mm_load_si128((const __m128i*)mul_table_low );
    const __m128i v_map_high = _mm_load_si128((const __m128i*)mul_table_high);

    __m128i v_y_low;
    __m128i v_y_high;
    __m128i v_y_map_low;
    __m128i v_y_map_high;

    __m128i v_x_low;
    __m128i v_x_high;
    __m128i v_x_map_low;
    __m128i v_x_map_high;

    for (int i = 0; i < ( size >> 5 ); i++) { //Assume PAD32
        // Equivalent to *z ^= mul_table_high[*x >> 4] ^ mul_table[*x & 15];

        v_y_low  = _mm_loadu_si128( (const __m128i*)(x) );              // load 16 bytes of x
        v_y_high = _mm_srli_epi16(v_y_low, 4);                          // shift to put high bits
        v_y_map_low = v_map_low;                                        // multiplication table for the low part
        v_y_map_high = v_map_high;                                      // multiplication table for the high part

        v_x_low  = _mm_loadu_si128( (const __m128i*)(x+16) );
        v_x_high = _mm_srli_epi16(v_x_low, 4);  // shift to put high bits
        v_x_map_low  = v_map_low;
        v_x_map_high = v_map_high;

        //apply mask [.... 1111 0000 1111]
        v_y_low =  _mm_and_si128(v_y_low  ,v_mask_low);
        v_y_high = _mm_and_si128(v_y_high ,v_mask_low);
        v_x_low =  _mm_and_si128(v_x_low  ,v_mask_low);
        v_x_high = _mm_and_si128(v_x_high ,v_mask_low);

        //table lookup
        v_y_map_low =  _mm_shuffle_epi8(v_y_map_low , v_y_low );
        v_y_map_high = _mm_shuffle_epi8(v_y_map_high, v_y_high);
        v_x_map_low =  _mm_shuffle_epi8(v_x_map_low , v_x_low );
        v_x_map_high = _mm_shuffle_epi8(v_x_map_high, v_x_high);
        
        //add the result of the multiplication of low and high part
        v_y_map_low = _mm_xor_si128(v_y_map_high, v_y_map_low);
        v_x_map_low = _mm_xor_si128(v_x_map_high, v_x_map_low);

        _mm_storeu_si128(reinterpret_cast<__m128i*>(x), v_y_map_low);                 // store xmm1 back into z

        _mm_storeu_si128(reinterpret_cast<__m128i*>(x+16), v_x_map_low);

        x += 32;
    }
    for (int i = 0; i < (size & 31); ++i) {
        *x = mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];
        x++;
    }
}
//===================================================================================================
void FiniteField::AddVec(SymbolType* z, SymbolType* x, int size)
{
    __m128i v_x;

    for (int i = 0; i < ( size >> 5); i++) { // Assume PAD32
        // Equivalent to *z ^= mul_table_high[*x >> 4] ^ mul_table[*x & 15];

        v_x  = _mm_loadu_si128( (const __m128i*)x );
        _mm_storeu_si128( reinterpret_cast<__m128i*>(z),
                            _mm_xor_si128(v_x,
                                _mm_loadu_si128( reinterpret_cast<__m128i*>(z) )));

        v_x  = _mm_loadu_si128( (const __m128i*)(x+16) );
        _mm_storeu_si128( reinterpret_cast<__m128i*>(z+16),
                            _mm_xor_si128(v_x,
                                _mm_loadu_si128( reinterpret_cast<__m128i*>(z+16) )));

        z += 32;
        x += 32;
    }
    for (int i = 0; i < (size & 31); ++i) { 
        *z ^= *x;
        z++;x++;
    }
}

//===================================================================================================
//__aarch64____aarch64____aarch64____aarch64____aarch64____aarch64____aarch64____aarch64____aarch64__
#elif defined(__aarch64__)
#include <arm_neon.h>
//===================================================================================================
void FiniteField::AddMulVec(SymbolType* z, SymbolType* x,
                            SymbolType c , int size) 
{
    const SymbolType* const mul_table_high = mul_table_o8_high_bits + (c << (FIELD_ORDER>>1));
    const SymbolType* const mul_table_low  = mul_table_o8_low_bits  + (c << (FIELD_ORDER>>1));

    const uint8x16_t v_mask_low = vdupq_n_u8(15);
    const uint8x16_t v_map_low  = vld1q_u8(mul_table_low);
    const uint8x16_t v_map_high = vld1q_u8(mul_table_high);

    for (int i = 0; i < (size >> 4); i++) {
        // Equivalent to *z ^= mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];

        uint8x16_t v_x      = vld1q_u8(x);
        uint8x16_t v_x_high = vandq_u8(vshrq_n_u8(v_x, 4), v_mask_low);
        uint8x16_t v_x_low  = vandq_u8(v_x, v_mask_low);
        auto v_product      = veorq_u8(vqtbl1q_u8(v_map_low, v_x_low), vqtbl1q_u8(v_map_high, v_x_high));

        vst1q_u8(z, veorq_u8(v_product, vld1q_u8(z)));

        z += 16;
        x += 16;
    }
    for (int i = 0; i < (size & 15); ++i) {
        *z ^= mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];
        z++;x++;
    }
}
//===================================================================================================
void FiniteField::MulVec(SymbolType* x, SymbolType c, int size) 
{   
    const SymbolType* const mul_table_high = mul_table_o8_high_bits + (c << (FIELD_ORDER>>1));
    const SymbolType* const mul_table_low  = mul_table_o8_low_bits  + (c << (FIELD_ORDER>>1));

    const uint8x16_t v_mask_low = vdupq_n_u8(15);
    const uint8x16_t v_map_low  = vld1q_u8(mul_table_low);
    const uint8x16_t v_map_high = vld1q_u8(mul_table_high);

    for (int i = 0; i < (size >> 4); i++) {
        // Equivalent to *x = mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];

        uint8x16_t v_x      = vld1q_u8(x);
        uint8x16_t v_x_high = vandq_u8(vshrq_n_u8(v_x, 4), v_mask_low);
        uint8x16_t v_x_low  = vandq_u8(v_x, v_mask_low);
        auto v_product      = veorq_u8(vqtbl1q_u8(v_map_low, v_x_low), vqtbl1q_u8(v_map_high, v_x_high));

        vst1q_u8(x, v_product);

        x += 16;
    }
    for (int i = 0; i < (size & 15); ++i) {
        *x = mul_table_high[*x >> 4] ^ mul_table_low[*x & 15];
        x++;
    }
}

//===================================================================================================
void FiniteField::AddVec(SymbolType* z, SymbolType* x, int size) 
{   
    for (int i = 0; i < (size >> 4); i++) {

        uint8x16_t v_x = vld1q_u8(x);
        vst1q_u8(z, veorq_u8(v_x, vld1q_u8(z)));

        z += 16;
        x += 16;
    }
    for (int i = 0; i < (size & 15); ++i) { 
        *z ^= *x;
        z++;x++;
    }
}
//===================================================================================================
#endif


//===================================================================================================
#if defined(__PCLMUL__) && defined(__SSE2__)
#include <immintrin.h>
SymbolType FiniteField::InnerProd(SymbolType* x, SymbolType* y, int size) {
/** basic idea: to compute the inner product of [x_15, x_14, ..., x_0] and [y_15, y_14, ..., y_0],
    construct polynomials f(t) = x_0 + x_1*t + ... + x_15*(t^15) and g(t) = y_15 + y_14*t + ... + y_0*(t^15);
    then the desired inner product is the coefficient of t^15 in f(t)*g(t).
    the CLMUL instruction set computes products of polynomials over GF(2); adjustments are made to adapt this to GF(256).
**/
    SymbolType z[16];
    const __m128i ShuffleRev = _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7);
    switch (size)
    {
        case 16:{ // hard-coded for size = 16
            __m128i X = _mm_loadu_si128( (const __m128i*)(x) ); // load all 16 bytes of x
            __m128i Y = _mm_loadu_si128( (const __m128i*)(y) ); // load all 16 bytes of y
            Y = _mm_shuffle_epi8(Y, ShuffleRev);        // shuffle y = [y_15, y_14, ..., y_0] into [y_8, y_9, ..., y_15, y_0, y_1, ..., y_7]
                                                        // NOT shuffled into [y_0, y_1, ..., y_15], because CLMUL can only handle 64 bits at a time

            const __m128i mask = _mm_set1_epi16(255);   // masks every other byte ...00000000111111110000000011111111

            __m128i X_low  = _mm_and_si128(X, mask);                                // get even-indexed bytes of X
            __m128i Y_high = _mm_and_si128(_mm_srli_epi16(Y,8), mask);              // get odd-indexed bytes of Y
            __m128i result = _mm_xor_si128(_mm_clmulepi64_si128(X_low,Y_high,0),    // compute the carryless multiply of the lower 64 bits of X_low and Y_high, then the upper 64 bits, then XOR together
                                        _mm_clmulepi64_si128(X_low,Y_high,0x11));
            
            X_low  = _mm_and_si128(_mm_srli_epi16(X,8), mask);  // get odd-indexed bytes of X
            Y_high = _mm_and_si128(Y, mask);                    // get even-indexed bytes of Y
            result = _mm_xor_si128(result,_mm_xor_si128(_mm_clmulepi64_si128(X_low,Y_high,0),
                                                        _mm_clmulepi64_si128(X_low,Y_high,0x11)));
                            // compute the carryless multiply of the lower 64 bits of X_low and Y_high, then the upper 64 bits, then XOR together with result

            _mm_storeu_si128( reinterpret_cast<__m128i*>(z), result);
            return overflow_table[z[7]]^z[6]; // the desired result is in bits 48-62; bits 56-62 (z[7]) are in 'overflow' and a table is consulted
            break;
        }
        case 8: { // hard-coded for size = 8
            __m128i X = _mm_loadu_si128( (const __m128i*)(x) ); // load all 16 bytes of x
            __m128i Y = _mm_loadu_si128( (const __m128i*)(y) ); // load all 16 bytes of y
            Y = _mm_shuffle_epi8(Y, ShuffleRev);        // shuffle y = [y_15, y_14, ..., y_0] into [y_8, y_9, ..., y_15, y_0, y_1, ..., y_7]
                                                        // NOT shuffled into [y_0, y_1, ..., y_15], because CLMUL can only handle 64 bits at a time

            const __m128i mask = _mm_set1_epi16(255);   // masks every other byte ...00000000111111110000000011111111

            __m128i X_low  = _mm_and_si128(X, mask);                                // get even-indexed bytes of X
            __m128i Y_high = _mm_and_si128(_mm_srli_epi16(Y,8), mask);              // get odd-indexed bytes of Y
            __m128i result = _mm_clmulepi64_si128(X_low,Y_high,0);    // compute the carryless multiply of the lower 64 bits of X_low and Y_high, then the upper 64 bits, then XOR together
            
            X_low  = _mm_and_si128(_mm_srli_epi16(X,8), mask);  // get odd-indexed bytes of X
            Y_high = _mm_and_si128(Y, mask);                    // get even-indexed bytes of Y
            result = _mm_xor_si128(result,_mm_clmulepi64_si128(X_low,Y_high,0));
                            // compute the carryless multiply of the lower 64 bits of X_low and Y_high, then the upper 64 bits, then XOR together with result

            _mm_storeu_si128( reinterpret_cast<__m128i*>(z), result);
            return overflow_table[z[7]]^z[6]; // the desired result is in bits 48-62; bits 56-62 (z[7]) are in 'overflow' and a table is consulted
            break;
        }
        case 4:{ // hard-coded for size = 4
            const __m128i zero = _mm_set_epi8(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  0x00,  0x00,  0x00,  0x00,  0x00);
            __m128i X = _mm_loadu_si128( (const __m128i*)x );
            __m128i Y = _mm_loadu_si128( (const __m128i*)y );
            X = _mm_unpacklo_epi8(X, zero);         // if x is [x_0, x_1, x_2, x_3] then X is (0,0,0,0,0,0,0,0, x_0, 0, x_1, 0, x_2, 0, x_3, 0)
            Y = _mm_unpacklo_epi8(Y, zero);
            Y = _mm_shuffle_epi8(Y, ShuffleRev);    // if y is [y_0, y_1, y_2, y_3] then Y is (0,0,0,0,0,0,0,0, 0, y_3, 0, y_2, 0, y_1, 0, y_0)

            _mm_storeu_si128( reinterpret_cast<__m128i*>(z), _mm_clmulepi64_si128(X,Y,0));
            return overflow_table[z[8]]^z[7];
            break;
        }
        default:{ // in case size was not 4, 8, or 16
            const __m128i mask = _mm_set_epi8(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  0x00,  0x00,  0x00,  0x00,  0xFF);
            __m128i result     = _mm_set_epi8(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  0x00,  0x00,  0x00,  0x00,  0x00);
            __m128i X = _mm_loadu_si128( reinterpret_cast<__m128i*>(x) );
            __m128i Y = _mm_loadu_si128( reinterpret_cast<__m128i*>(y) );

            for (int i = 0; i < size; i++) { // Doing one by one
                result = _mm_xor_si128(result,_mm_clmulepi64_si128(_mm_and_si128(X, mask),_mm_and_si128(Y, mask),0));
                X = _mm_srli_si128(X,1);
                Y = _mm_srli_si128(Y,1);
            }
            _mm_storeu_si128( reinterpret_cast<__m128i*>(z), result);
            return overflow_table[z[1]]^z[0];
        }
    }

}
#elif defined(__aarch64__)
//===================================================================================================
//TODO
SymbolType FiniteField::InnerProd(SymbolType* x, SymbolType* y, int size)
{
    SymbolType result = 0;

    for (int i = 0; i < size; i++) {
        result ^= MulElem(x[i], y[i]);
    }
    return result;
}
//===================================================================================================
#else  // if we did not have PLCMUL, use the dumb way
SymbolType FiniteField::InnerProd(SymbolType* x, SymbolType* y, int size)
{
    SymbolType result = 0;
    for (int i = 0; i < size; i++) {
        result ^= MulElem(x[i], y[i]);
    }
    return result;
}
//===================================================================================================
#endif

#endif

void FiniteField::Print(SymbolType* mat, int size, char mode)
{
    if (mode != 'X' && mode != 'x' && mode != 'b') {
        std::cout << "invalid mode" << std::endl;
        return;
    }
    std::cout << std::setw(4) << "i |" << std::setw(7) << "mat[i]" << "|" << std::endl;
    for (int i = 0; i < size; ++i) {
        switch(mode) {
        case 'X':
        case 'x':
            printf("%2d |%7X|\n", i, mat[i]);
            break;
        
        case 'b':
            std::cout << std::setw(2) << i << " |" << std::setw(7) << std::bitset<4>(mat[i]) << "|" << std::endl;
            break;
        }
    }
}
