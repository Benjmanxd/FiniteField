#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <stdint.h>
#include <stdio.h>

typedef unsigned char SymbolType;
typedef uint16_t IdType;
#define MAX_DEGREE 256
#define FIELD_ORDER 8
#define FIELD_SIZE (1 << FIELD_ORDER)
#define MAX_OVERHEAD 4
#define PAD32(x) (((x) + 31) & -32)
#endif
