#pragma once

#ifndef i_am_noob

int dot_product_mod(int n, __m256i *a, __m256i *b, int mod){
    __m256i zero=_mm256_setzero_si256();
    __m256i res=zero;
    __m256i four_mod_mod=_mm256_set1_epi64x(4ll*mod*mod);
    for(int i=0; i<n; ++i){
        __m256i c=_mm256_mul_epu32(a[i],b[i]);
        res=_mm256_add_epi64(res,c);
        if(!(i&3)){
            __m256i msk=_mm256_cmpgt_epi64(res,four_mod_mod);
            __m256i val=_mm256_blendv_epi8(zero,four_mod_mod,msk);
            res=_mm256_sub_epi64(res,val);
        }
    }
    ll val[4]{};
    _mm256_storeu_si256((__m256i*)&val,res);
    ll tot=0;
    for(int i=0; i<4; ++i){
        tot+=val[i];
        tot%=mod;
    }
    return (int)tot;
}

#endif