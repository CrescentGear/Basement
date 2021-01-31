#ifndef RH_UTILITY_H
#define RH_UTILITY_H 

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
 extern "C" {
#endif


/*=====================================================================
 > Common
======================================================================*/

#define __map(val,i_min,i_max,o_min,o_max)   (double)( ( ((double)o_max)*(((double)val)-((double)i_min))+((double)o_min)*((double)(i_max)-(double)(val)) )/((double)(i_max)-(double)(i_min)) )

#define __abs(val)                           (((val)>0)?(val):(-(val)))

//#define __round(a)                           (int)((a)+0.5)>(int)(a)?((int)(a)+1):((int)(a))
//#define __round1000(a)                       (double)((__round((a)*1000.0))/1000.0)


#pragma anon_unions

/*=====================================================================
 > Algebra Reference 
======================================================================*/
long    __sign       (long x);
long    __sqrt       (long x);
double  __gussian    (long x,long __miu  ,double __sigma);
double  __gussian2D  (long x,long y      ,double __sigma);

struct __Kernel_t{
    uint16_t*   pBuffer; //  <- This buffer array should only be created by "__malloc".
    size_t      order;
    int32_t     sum;
};
typedef struct __Kernel_t       __Kernel_t;
__Kernel_t* __gussianKernel(double __sigma,size_t order,__Kernel_t* pKernel);
 
/*=====================================================================
 > Quantity Reference 
======================================================================*/
struct IntArray_t{
	size_t  index;
	int     value;
};
typedef struct IntArray_t IntArray_t;
IntArray_t __findMax_INT(const int* pValue,size_t num);
IntArray_t __findMin_INT(const int* pValue,size_t num);

struct UintArray_t{
	size_t        index;
	unsigned int  value;
};
typedef struct UintArray_t UintArray_t;

/*=====================================================================
 > Geometry Reference 
======================================================================*/
struct Vector2D_t{
	int x;
	int y;
};
typedef struct Vector2D_t Vector2D_t;
typedef struct Vector2D_t Point2D_t;

struct Vector3D_t{
	int x;
	int y;
	int z;
};
typedef struct Vector3D_t Vector3D_t;
typedef struct Vector3D_t Point3D_t;

Point3D_t __findPoint_LineCross        (const Point3D_t  line1[2] ,const Point3D_t  line2[2]);
Point3D_t __findPoint_VectorDistance   (const Point3D_t* A        ,const Point3D_t* B     ,int    dist_AP ); 
Point3D_t __findPoint_VectorProportion (const Point3D_t* A        ,const Point3D_t* B     ,double scale   );


int        __Vect2D_Dot    (const Vector2D_t* vect1,const Vector2D_t* vect2);
int        __Vect3D_Dot    (const Vector3D_t* vect1,const Vector3D_t* vect2);
Vector3D_t __Vect3D_Cross  (const Vector3D_t* vect1,const Vector3D_t* vect2);


int        __Dir_Line        (int xs,int ys,int xe,int ye);
int        __Point_toLine    (int xs,int ys,int xe,int ye,               int px,int py);
int        __Point_toTriangle(int x1,int y1,int x2,int y2,int x3,int y3, int px,int py);
int        __Point_toCircle  (int xc,int yc,int radius,                  int px,int py);

/*===========================================================================================================================
 > DSP Reference
============================================================================================================================*/
 struct __ComplexFLOAT_t{
     float_t  real;
     float_t  imag;
 };
 typedef struct __ComplexFLOAT_t __ComplexFLOAT_t;

 struct __ComplexINT_t{
     long  real;
     long  imag;
 };
 typedef struct __ComplexINT_t __ComplexINT_t;
  
void      __rDFT_Float    (const float_t*          src, float_t* dst_m, __ComplexFLOAT_t* dst_c, size_t dftLen);
void      __cDFT_Float    (const float complex*    src, float_t* dst_m, float complex*    dst_c, size_t dftLen);
void      __rIDFT_Float   (const float_t*          src, float_t* dst_m, __ComplexFLOAT_t* dst_c, size_t dftLen);
void      __cIDFT_Float   (const float complex*    src, float_t* dst_m, float complex*    dst_c, size_t dftLen);

/*=====================================================================
 > Image Processing Reference 
======================================================================*/

struct __PixelRGB565_t{
    uint16_t R : 5;
    uint16_t G : 6;
    uint16_t B : 5;
};
union __UNION_PixelRGB565_t{
    uint16_t R : 5;
    uint16_t G : 6;
    uint16_t B : 5;
    uint16_t data;
};
typedef struct  __PixelRGB565_t         __PixelRGB565_t;
typedef union   __UNION_PixelRGB565_t   __UNION_PixelRGB565_t;

 struct __PixelRGB888_t{
     uint8_t B ;
     uint8_t G ;
     uint8_t R ;
 };
union __UNION_PixelRGB888_t{
    struct{
        uint8_t B : 8;
        uint8_t G : 8;
        uint8_t R : 8;
    };
    uint32_t data;
};
typedef struct  __PixelRGB888_t         __PixelRGB888_t;
typedef union   __UNION_PixelRGB888_t   __UNION_PixelRGB888_t;

struct __ImageRGB565_t{
    __UNION_PixelRGB565_t* pBuffer;
    size_t      width;
    size_t      height;
};
typedef struct __ImageRGB565_t  __ImageRGB565_t;
 

struct __ImageRGB888_t{
    __UNION_PixelRGB888_t* pBuffer;
    size_t      width;
    size_t      height;
};
typedef struct __ImageRGB888_t  __ImageRGB888_t;

__ImageRGB888_t* __LoadBMP_ImgRGB888             (const char* __restrict__ path);
__ImageRGB888_t* __OutBMP_ImgRGB888              (const char* __restrict__ path,__ImageRGB888_t* p);
__ImageRGB888_t* __Create_ImgRGB888              (size_t width,size_t height);
__ImageRGB888_t* __CopyBMP_ImgRGB888             (const __ImageRGB888_t* src,__ImageRGB888_t* dst);

__ImageRGB888_t* __FreeBuffer_ImgRGB888          (__ImageRGB888_t* ptr);
void             __Free_ImgRGB888                (__ImageRGB888_t* ptr);

__ImageRGB888_t* __Filter_Gray_ImgRGB888         (const __ImageRGB888_t* src,__ImageRGB888_t* dst); //
__ImageRGB888_t* __Filter_Cold_ImgRGB888         (const __ImageRGB888_t* src,__ImageRGB888_t* dst); //
__ImageRGB888_t* __Filter_Warm_ImgRGB888         (const __ImageRGB888_t* src,__ImageRGB888_t* dst); //
__ImageRGB888_t* __Filter_OTUS_ImgRGB888         (const __ImageRGB888_t* src,__ImageRGB888_t* dst); //
 
__ImageRGB888_t* __Blur_Gussian_ImgRGB888        (const __ImageRGB888_t* src,__ImageRGB888_t* dst,unsigned int _0_65535_);
__ImageRGB888_t* __Blur_Average_ImgRGB888        (const __ImageRGB888_t* src,__ImageRGB888_t* dst,unsigned int _0_65535_);

__ImageRGB888_t* __Interpo_NstNeighbor_ImgRGB888 (const __ImageRGB888_t* src,__ImageRGB888_t* dst,size_t height,size_t width); //

__ImageRGB565_t* __Conv2D_ImgRGB565              (const __ImageRGB565_t* src,__ImageRGB565_t* dst,const __Kernel_t* k);
__ImageRGB888_t* __Conv2D_ImgRGB888              (const __ImageRGB888_t* src,__ImageRGB888_t* dst,const __Kernel_t* k);


/*=====================================================================
 > Memory Programming Reference 
======================================================================*/
#define __malloc(x)                    malloc(x)//__mallocHEAP(x)
#define __free(x)                      free(x)//__freeHEAP(x)
#define __VIRTUAL_HEAP_SIZE_BYTE       (1<<20)
 
void* __mallocHEAP(size_t size);
void  __freeHEAP(void* ptr);
void* __memsetWORD  (void* __b, int  value, size_t num);
void* __memsetDWORD (void* __b, long value, size_t num);
 
#ifdef __cplusplus
 }
#endif

#endif

