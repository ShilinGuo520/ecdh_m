#ifndef _ECDH_H_
#define _ECDH_H_

#include <inttypes.h>

//#define UINT64

typedef int8_t wordcount_t;
typedef int16_t bitcount_t;
typedef int8_t cmpresult_t;

typedef uint32_t uECC_word_t;
#ifdef UINT64
typedef uint64_t uECC_dword_t;
#endif
typedef const struct uECC_Curve_t * uECC_Curve;
#define uECC_VLI_API static

#define uECC_RNG_MAX_TRIES 64

#define num_words_secp192r1 6
#define num_bytes_secp192r1 24

#define HIGH_BIT_SET 0x80000000
#define uECC_WORD_BITS 32
#define uECC_WORD_BITS_SHIFT 5
#define uECC_WORD_BITS_MASK 0x01F

#define uECC_MAX_WORDS 24
#define uECC_WORD_SIZE 4

struct uECC_Curve_t {
    wordcount_t num_words;
    wordcount_t num_bytes;
    bitcount_t num_n_bits;
    uECC_word_t p[uECC_MAX_WORDS];
    uECC_word_t n[uECC_MAX_WORDS];
    uECC_word_t G[uECC_MAX_WORDS * 2];
    uECC_word_t b[uECC_MAX_WORDS];
    void (*double_jacobian)(uECC_word_t * X1,
                            uECC_word_t * Y1,
                            uECC_word_t * Z1,
                            uECC_Curve curve);
    void (*mod_sqrt)(uECC_word_t *a, uECC_Curve curve);
    void (*x_side)(uECC_word_t *result, const uECC_word_t *x, uECC_Curve curve);
    void (*mmod_fast)(uECC_word_t *result, uECC_word_t *product);
};

#define BYTES_TO_WORDS_8(a, b, c, d, e, f, g, h) 0x##d##c##b##a, 0x##h##g##f##e
#define BITS_TO_BYTES(num_bits) ((num_bits + 7) / 8)
#define BITS_TO_WORDS(num_bits) ((num_bits + ((uECC_WORD_SIZE * 8) - 1)) / (uECC_WORD_SIZE * 8))
#define EccPoint_isZero(point, curve) uECC_vli_isZero((point), (curve)->num_words * 2)

void vli_print(uint8_t *vli, unsigned int size) {
        unsigned int i;
    for(i=0; i<size; ++i) {
        printf("%02X ", (unsigned)vli[i]);
    }
}

#endif

