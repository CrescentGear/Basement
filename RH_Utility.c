#include "RH_Utility.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

/*=========================================
 > Common 
==========================================*/


/*=========================================
 > Algebra Reference 
==========================================*/

long __sign(long x){
    return (x>=0)?(1):(-1);
}

long __sqrt(long x){
    if(x <= 0) return 0;
    long l   = 1;
    long r   = x;
    long res = 0;
    while(l <= r){
        long mid=(l+r)>>1;
        if(mid <= x/mid){
          l   = mid+1;
          res = mid;
      }else{
          r = mid-1;
      }
    }
    if( ((res+1)*(res+1) - x) > (x - res*res) )
        return res;
    return (res+1);
}

#ifndef M_2_SQRTPI
#define M_2_SQRTPI  1.12837916709551257389615890312154517   /* 2/sqrt(pi)     */
#endif

#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880168872420969808   /* sqrt(2)        */
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif

double __gussian(long x,long __miu,long __sigma){
    return ((__sigma==0)?(0):(double)((M_2_SQRTPI/((__sigma<<1)*(M_SQRT2)))*exp(-(x-__miu)*(x-__miu)/(double)((__sigma*__sigma)<<1))));
}

double __gussian2D(long x,long y,long __sigma){
// Same Effect but slower,only suitable when __sigma is a value of double.
//    return ((__sigma==0)?(0):((double)((1/(2*M_PI*__sigma*__sigma))*exp(-((x*x)+(y*y))/(double)((__sigma*__sigma)*2)))));
    return ((__sigma==0)?(0):((double)((1/(M_PI*__sigma*(__sigma<<1)))*exp(-((x*x)+(y*y))/(double)((__sigma*__sigma<<1))))));
}

void __generate_GussianKernel(int order,int* buffer){

}

/*=========================================
 > Quantity Reference 
==========================================*/

struct IntArray_t __findMax_INT(const int* pValue,size_t num){
    int max = *pValue;
    int cnt = 0;
    while(num--){
        if(*pValue > max)
            max = *pValue;
        pValue++;
        cnt++;
    }
    struct IntArray_t result = {.index = cnt,.value = max};
    return result;
}

struct IntArray_t __findMin_INT(const int* pValue,size_t num){
    int min = *pValue;
    int cnt = 0;
    while(num--){
        if(*pValue < min)
            min = *pValue;
        pValue++;
        cnt++;
    }
    struct IntArray_t result = {.index = cnt,.value = min};
    return result;
}


/*=========================================
 > Geometry Reference 
==========================================*/

Point3D_t __findPoint_VectorDistance (const Point3D_t*  A,const Point3D_t*  B,int dist_AP){ 
    long dist_AB = lroundl(sqrt( (B->x - A->x)*(B->x - A->x) + \
                                (B->y - A->y)*(B->y - A->y) + \
                                (B->z - A->z)*(B->z - A->z)  ));
    Point3D_t result = {
        .x = (B->x - A->x)*dist_AP/dist_AB + A->x ,
        .y = (B->y - A->y)*dist_AP/dist_AB + A->y ,
        .z = (B->z - A->z)*dist_AP/dist_AB + A->z ,
    };

    return result;
}

Point3D_t __findPoint_VectorProportion (const Point3D_t*  A,const Point3D_t*  B,double scale){ 
    Point3D_t result = {
        .x = (B->x - A->x)*scale + A->x ,
        .y = (B->y - A->y)*scale + A->y ,
        .z = (B->z - A->z)*scale + A->z ,
    };
    return result;
}

int __Vect2D_Dot(const Vector2D_t* vect1,const Vector2D_t* vect2){
    return (int)((vect1->x*vect2->x)+(vect1->y*vect2->y));
}

int __Vect3D_Dot(const Vector3D_t* vect1,const Vector3D_t* vect2){
    return (int)((vect1->x*vect2->x)+(vect1->y*vect2->y)+(vect1->z*vect2->z));
}

Vector3D_t __Vect3D_Cross(const Vector3D_t* vect1,const Vector3D_t* vect2){
    Vector3D_t vecResult = {.x = ( vect1->y*vect2->z - vect1->z*vect2->y),\
                            .y = ( vect1->z*vect2->x - vect1->x*vect2->z),\
                            .z = ( vect1->x*vect2->y - vect1->y*vect2->x)};
    return vecResult;
}

 // -1    = Line is negative.
 //  0    = Line is horizontal.
 //  1    = Line is positive.
 // 65535 = Line is vertical.
int __Dir_Line(int xs,int ys,int xe,int ye){
    if(xs==xe)
        return 65535;
    if(ys==ye)
        return 0;

    return ((xe-xs)*(ye-ys)>0)?(1):(-1);
}

 // -1 = (px,py) is below the line.
 //  0 = (px,py) is at the line.
 //  1 = (px,py) is above the line.
int __Point_toLine(int xs,int ys,int xe,int ye,int px,int py){
    int param_1 = (xe>xs)?( (xe-xs)*py ):( (xs-xe)*py );
    int param_2 = (xe>xs)?( (ye-ys)*px+(ye*(xe-xs)-xe*(ye-ys)) ):( (ys-ye)*px+(ye*(xs-xe)-xe*(ys-ye)) );

    if(param_1 > param_2)
        return 1;
    else if(param_1 < param_2)
        return -1;
    else
        return 0;
}

 // -1 = (px,py) is outside the triangle
 //  0 = (px,py) is at the edge of triangle
 //  1 = (px,py) is inside the triangle
int __Point_toTriangle(int x1,int y1,int x2,int y2,int x3,int y3,int px,int py){
 // Condition:
 // P = A + u*(CA) + v*(BA)
 // u >= 0 && v >= 0 && u+v <= 1
 
 // Any point can be represented by: (PA) = u*(CA) + v*(BA)
 // 
 // When both multiply by (CA) and (BA):
 // (PA)·(CA) = u*[(CA)·(CA)] + v*[(BA)·(CA)]   
 // (PA)·(BA) = u*[(BA)·(CA)] + v*[(BA)·(BA)]  
 
 // Then:
 //         [(BA)·(BA)]*[(PA)·(CA)] - [(BA)·(CA)]*[(PA)·(BA)]
 // u = ---------------------------------------------------------
 //         [(CA)·(CA)]*[(BA)·(BA)] - [(CA)·(BA)]*[(BA)·(CA)]
 
 //         [(CA)·(CA)]*[(PA)·(BA)] - [(CA)·(CA)]*[(PA)·(CA)]
 // v = ---------------------------------------------------------
 //         [(CA)·(CA)]*[(BA)·(BA)] - [(CA)·(BA)]*[(BA)·(CA)]
 
 // Assume A = (x1,y1) | B = (x2,y2) | C = (x3,y3) :
    struct Vector2D_t v0 = {.x = x3-x1,.y = y3-y1};
    struct Vector2D_t v1 = {.x = x2-x1,.y = y2-y1};
    struct Vector2D_t v2 = {.x = px-x1,.y = py-y1};

    int v00 = __Vect2D_Dot(&v0,&v0);
    int v01 = __Vect2D_Dot(&v0,&v1);
    int v02 = __Vect2D_Dot(&v0,&v2);
    int v11 = __Vect2D_Dot(&v1,&v1);
    int v12 = __Vect2D_Dot(&v1,&v2);

    int u = v11*v02-v01*v12;
    int v = v00*v12-v01*v02;
    int d = v00*v11-v01*v01;
    if(u<0 || v<0)
        return -1;
    else if(u==0 || v==0)
        return 0;

    if(u+v > d)
        return -1;
    else if(u+v < d)
        return 1;
    else
        return 0;
}

 // -1 = (px,py) is outside the circle
 //  0 = (px,py) is at the edge of circle
 //  1 = (px,py) is inside the circle
int __Point_toCircle(int xc,int yc,int radius,int px,int py){
    int key = (xc-px)*(xc-px)+(yc-py)*(yc-py);
    int r_2 = radius*radius;
    if(key > r_2)
        return -1;
    else if(key < r_2)
        return 1;
    else
        return 0;
}

/*=========================================
 > Image Processing Reference 
==========================================*/
#if defined    (_WIN32)
#include <windows.h>
#elif defined  (__APPLE__)
typedef uint8_t   BYTE;
typedef uint16_t  WORD;
typedef uint32_t  DWORD;
typedef int32_t   LONG;
typedef float     FLOAT;
typedef FLOAT     *PFLOAT;
typedef BYTE      *PBYTE;
typedef BYTE      *LPBYTE;
typedef int       *PINT;
typedef int       *LPINT;
typedef WORD      *PWORD;
typedef WORD      *LPWORD;
typedef int32_t   *LPLONG;
typedef DWORD     *PDWORD;
typedef DWORD     *LPDWORD;
typedef void      *LPVOID;
#pragma pack(1)
typedef struct tagBITMAPFILEHEADER {
    WORD   bfType;
    DWORD  bfSize;
    WORD   bfReserved1;
    WORD   bfReserved2;
    DWORD  bfOffBits;
} BITMAPFILEHEADER;

#pragma pack(1)
typedef struct tagBITMAPINFOHEADER {
    DWORD  biSize;
    LONG   biWidth;
    LONG   biHeight;
    WORD   biPlanes;
    WORD   biBitCount;
    DWORD  biCompression;
    DWORD  biSizeImage;
    LONG   biXPelsPerMeter;
    LONG   biYPelsPerMeter;
    DWORD  biClrUsed;
    DWORD  biClrImportant;
} BITMAPINFOHEADER;

#endif
__ImageRGB888_t* __LoadBMP_ImgRGB888(const char* __restrict__ path){
    FILE* bmp;
    BITMAPFILEHEADER fileHead;
    BITMAPINFOHEADER infoHead;

    __ImageRGB888_t* pIMG = __malloc(sizeof(__ImageRGB888_t));
    pIMG->height  = 0;
    pIMG->width   = 0;
    pIMG->pBuffer = NULL;
    
    bmp = fopen(path, "r");
    if (bmp == NULL) {
        printf("open error\n");
        return pIMG;
    }
    fseek(bmp, 0L, SEEK_SET);
    fread(&fileHead, sizeof(BITMAPFILEHEADER), 1, bmp);
    fread(&infoHead, sizeof(BITMAPINFOHEADER), 1, bmp);

    if (fileHead.bfType != 0x4D42) {
        printf("This not a *.bmp file\n");
        return pIMG;
    }

    fseek(bmp, fileHead.bfOffBits, SEEK_SET);

    pIMG->pBuffer = (__UNION_PixelRGB888_t*)__malloc(infoHead.biWidth * infoHead.biHeight * sizeof(__UNION_PixelRGB888_t));
    
    for (int row = 0; row < infoHead.biHeight; row++) {
        for (int col = 0; col < infoHead.biWidth; col++) {
            fread(&(pIMG->pBuffer[(infoHead.biHeight - row - 1)*infoHead.biWidth + col].data), sizeof(__PixelRGB888_t), 1, bmp);

//            printf("[%d]: ",(infoHead.biHeight - row - 1)*infoHead.biWidth + col);
//            printf("R=%d G=%d B=%d\n",\
//                   pIMG->pBuffer[(infoHead.biHeight - row - 1)*infoHead.biWidth + col].R,\
//                   pIMG->pBuffer[(infoHead.biHeight - row - 1)*infoHead.biWidth + col].G,\
//                   pIMG->pBuffer[(infoHead.biHeight - row - 1)*infoHead.biWidth + col].B);
        }
        int eps = (4-(infoHead.biWidth*sizeof(__PixelRGB888_t))%4)%4;
        uint8_t dummyByte;
        while(eps--){
            fread(&dummyByte,sizeof(char) ,1 , bmp);
        }
    }
    fclose(bmp);

    pIMG->width   = infoHead.biWidth;
    pIMG->height  = infoHead.biHeight;

    return pIMG;
}

__ImageRGB888_t* __CopyBMP_ImgRGB888(const __ImageRGB888_t* src,__ImageRGB888_t* dst){
    memcpy(dst->pBuffer, src->pBuffer, (src->height)*(src->width)*sizeof(__UNION_PixelRGB888_t));
    dst->height = src->height;
    dst->width  = src->width;
    return dst;
}

__ImageRGB888_t* __Create_ImgRGB888(size_t width,size_t height){
    __ImageRGB888_t* pIMG = __malloc(sizeof(__ImageRGB888_t));
    if(pIMG == NULL)
        return NULL;
    pIMG->height          = height;
    pIMG->width           = width;
    pIMG->pBuffer         = __malloc((pIMG->height)*(pIMG->width)*sizeof(pIMG->pBuffer[0]));
    if(pIMG->pBuffer == NULL){
        __free(pIMG);
        return NULL;
    }
    memset(pIMG->pBuffer, 0, (pIMG->height)*(pIMG->width)*sizeof(pIMG->pBuffer[0]));
    return pIMG;
}

void __OutBMP_ImgRGB888(const char* __restrict__ path,__ImageRGB888_t* p){
    FILE* bmp;
    if(p == NULL && p->pBuffer == NULL)
        return;
    int eps = (4-(p->width*sizeof(__PixelRGB888_t))%4)%4;
    BITMAPFILEHEADER fileHead = {
        .bfOffBits      = 40 + 14   ,
        .bfReserved1    = 0         ,
        .bfReserved2    = 0         ,
        .bfSize         = (uint32_t)(p->height * p->width * sizeof(__PixelRGB888_t) + 54),
        .bfType         = 0x4D42    , 
    };
    BITMAPINFOHEADER infoHead = {
        .biSize          = 40        ,
        .biWidth         = (int)(p->width)  ,
        .biHeight        = (int)(p->height) ,
        .biPlanes        = 1         ,
        .biBitCount      = 8+8+8     ,
        .biCompression   = 0         , 
        .biSizeImage     = (uint32_t)(p->height*p->width*sizeof(__PixelRGB888_t) + eps*(p->height)) ,
        .biClrUsed       = 0         ,
        .biClrImportant  = 0         ,
        .biXPelsPerMeter = 0         ,
        .biYPelsPerMeter = 0         ,
    };
    printf("%d\n" ,infoHead.biSizeImage    );
    bmp = fopen(path,"wb");

    if(bmp == NULL) return;

    // RGB Sequence should be reversed.
    fseek(bmp,0L,SEEK_SET);
    fwrite(&fileHead ,1 ,sizeof(BITMAPFILEHEADER) , bmp);
    fwrite(&infoHead ,1 ,sizeof(BITMAPINFOHEADER) , bmp);
    fseek(bmp,54L,SEEK_SET);
    for (int row = 0; row < p->height; row++) {
        for (int col = 0; col < p->width; col++) {
            fwrite( &p->pBuffer[(infoHead.biHeight - row - 1) * infoHead.biWidth + col] ,sizeof(__PixelRGB888_t) ,1 ,bmp );
        }
        int eps = (4-(infoHead.biWidth*sizeof(__PixelRGB888_t))%4)%4;
        uint8_t dummyByte = 0x00;
        while(eps--){
            fwrite(&dummyByte,sizeof(char) ,1 , bmp);
        }
    }
    
    fclose(bmp);

}

__ImageRGB888_t* __FreeBuffer_ImgRGB888(__ImageRGB888_t* ptr){
    __free(ptr->pBuffer);
    ptr->height  = 0;
    ptr->width   = 0;
    ptr->pBuffer = NULL;
    return ptr;
}

void __Free_ImgRGB888(__ImageRGB888_t* ptr){
    __free(__FreeBuffer_ImgRGB888(ptr));
}

__ImageRGB888_t* __Filter_Gray_ImgRGB888(const __ImageRGB888_t* in,__ImageRGB888_t* out){
    
    if(in != NULL && out != NULL){
        if(in->pBuffer != NULL && out->pBuffer != NULL){
            for (int row = 0; row < in->height; row++) {
                for (int col = 0; col < in->width; col++) {
                    long temp = lroundl(0.299 * in->pBuffer[row * in->width + col].R + 0.587 * in->pBuffer[row * in->width + col].G + 0.114 * in->pBuffer[row * in->width + col].B);
                    out->pBuffer[row * in->width + col].data = (uint32_t)(((temp&0xff)<<16)|((temp&0xff)<<8)|((temp&0xff)));
                }
            }
        }
    }
    return out;
}

__ImageRGB888_t* __Filter_Warm_ImgRGB888(const __ImageRGB888_t* in,__ImageRGB888_t* out){
    return out;
}

__ImageRGB888_t* __Filter_Cold_ImgRGB888(const __ImageRGB888_t* in,__ImageRGB888_t* out){
    
    return out;
}

__ImageRGB888_t* __Blur_Gussian_ImgRGB888(const __ImageRGB888_t* src,__ImageRGB888_t* dst,unsigned int _0_100_){
    int      kernel[1+2+3+4+5] = {0};
    uint16_t gus_kernel[9][9];
    
    int sigma = 0;
        
    double temp = 12.8;
     
    kernel[0]  = (uint16_t)(100 / (_PI * temp));                       //[0][0]

    kernel[1]  = (uint16_t)(100 / (_PI * temp) * exp(-1  / (temp)));   //[1][0]
    kernel[2]  = (uint16_t)(100 / (_PI * temp) * exp(-2  / (temp)));   //[1][]
     
    kernel[3]  = (uint16_t)(100 / (_PI * temp) * exp(-4  / (temp)));   //[2][0]
    kernel[4]  = (uint16_t)(100 / (_PI * temp) * exp(-5  / (temp)));   //[2][1]
    kernel[5]  = (uint16_t)(100 / (_PI * temp) * exp(-8  / (temp)));   //[2][2]
    
    kernel[6]  = (uint16_t)(100 / (_PI * temp) * exp(-9  / (temp)));   //[3][0]
    kernel[7]  = (uint16_t)(100 / (_PI * temp) * exp(-10 / (temp)));   //[3][1]
    kernel[8]  = (uint16_t)(100 / (_PI * temp) * exp(-13 / (temp)));   //[3][2]
    kernel[9]  = (uint16_t)(100 / (_PI * temp) * exp(-18 / (temp)));   //[3][3]
    
    kernel[10] = (uint16_t)(100 / (_PI * temp) * exp(-16 / (temp)));   //[4][0]
    kernel[11] = (uint16_t)(100 / (_PI * temp) * exp(-17 / (temp)));   //[4][1]
    kernel[12] = (uint16_t)(100 / (_PI * temp) * exp(-20 / (temp)));   //[4][2]
    kernel[13] = (uint16_t)(100 / (_PI * temp) * exp(-25 / (temp)));   //[4][3]
    kernel[14] = (uint16_t)(100 / (_PI * temp) * exp(-32 / (temp)));   //[4][4]
    
    sigma += kernel[0];        // 1
    sigma += (kernel[1]) << 2; // 4
    sigma += (kernel[2]) << 2; // 4
    
    sigma += (kernel[3]) << 2; // 4
    sigma += (kernel[4]) << 3; // 8
    sigma += (kernel[5]) << 2; // 4
   
    sigma += (kernel[6]) << 2; // 4
    sigma += (kernel[7]) << 3; // 8
    sigma += (kernel[8]) << 3; // 8
    sigma += (kernel[9]) << 2; // 4
   
    sigma += (kernel[10]) << 2; // 4
    sigma += (kernel[11]) << 3; // 8
    sigma += (kernel[12]) << 3; // 8
    sigma += (kernel[13]) << 3; // 8
    sigma += (kernel[14]) << 2; // 4

    int center   = ((9-1)>>1);

    for(int i = 0;i <= center;i++){
        for(int j = 0;j <= i;j++){
            gus_kernel[center + i][center + j] = kernel[i+j];
            gus_kernel[center + i][center - j] = kernel[i+j];
            gus_kernel[center - i][center + j] = kernel[i+j];
            gus_kernel[center - i][center - j] = kernel[i+j];
            gus_kernel[center + j][center + i] = kernel[i+j];
            gus_kernel[center + j][center - i] = kernel[i+j];
            gus_kernel[center - j][center + i] = kernel[i+j];
            gus_kernel[center - j][center - i] = kernel[i+j];
        }
    }
    sigma = 0;
    for(int i=0;i<9;i++){
        for(int j=0;j<9;j++){
            printf("[%3d]",gus_kernel[i][j]);
            sigma += gus_kernel[i][j];
        }
        printf("\n");
    }

    printf("sigma = %d\n",sigma);
    
    __Kernel_t kernel_info = {
        .order   = 9,
        .pBuffer = gus_kernel[0]
    };
    
    
    
    return __Conv2D_ImgRGB888(src, dst, &kernel_info, sigma);
}

__ImageRGB565_t* __Conv2D_ImgRGB565(const __ImageRGB565_t* in,__ImageRGB565_t* out,const __Kernel_t* k,int coe){
    if( in == NULL || k == NULL){
        return out;
    }
    
    if(out == NULL){
        out = (__ImageRGB565_t*)__malloc(sizeof(__ImageRGB565_t));
        if(out == NULL) // Not enough space :-(
            return out;
        out->pBuffer = (__UNION_PixelRGB565_t*)__malloc(in->width * in->height * sizeof(__UNION_PixelRGB565_t));
        if(out->pBuffer == NULL) // Not enough space :-(
            return out;
    }
    
    if(out->pBuffer == NULL){
        out->pBuffer = (__UNION_PixelRGB565_t*)__malloc(in->width * in->height * sizeof(__UNION_PixelRGB565_t));
        if(out->pBuffer == NULL) // Not enough space :-(
            return out;
    }

    for(int j=0;j<in->height;j++){
        for(int i=0;i<in->width;i++){
            unsigned long tmp_R = 0,tmp_G = 0,tmp_B = 0;
            for(int n=0;n<k->order;n++){
                for(int m=0;m<k->order;m++){
                    int offset_y = (int)(j-(k->order>>1)+n);
                    int offset_x = (int)(i-(k->order>>1)+m);
                    if(offset_x<0||offset_y<0||offset_x>=in->width||offset_y>=in->height){
                         //...//
                        continue;
                    }
                    unsigned int select_R  = (in->pBuffer + offset_y*in->width + offset_x)->R;
                    unsigned int select_G  = (in->pBuffer + offset_y*in->width + offset_x)->G;
                    unsigned int select_B  = (in->pBuffer + offset_y*in->width + offset_x)->B;
                    int       selectKernel = *( k->pBuffer + n       * k->order + m       );
                    tmp_R += ( (select_R) * (selectKernel) );
                    tmp_G += ( (select_G) * (selectKernel) );
                    tmp_B += ( (select_B) * (selectKernel) );
                }
            }
            size_t offset = (j*in->width)+i;
            if(offset < out->width*out->height){
                (out->pBuffer+offset)->R = tmp_R/coe;//(coe==0)?((1<<5)-1):tmp_R/coe;
                (out->pBuffer+offset)->G = tmp_G/coe;//(coe==0)?((1<<6)-1):tmp_G/coe;
                (out->pBuffer+offset)->B = tmp_B/coe;//(coe==0)?((1<<5)-1):tmp_B/coe;
            }
        }
    }
    return out;
}

__ImageRGB888_t* __Conv2D_ImgRGB888(const __ImageRGB888_t* in,__ImageRGB888_t* out,const __Kernel_t* k,int coe){
    if( in == NULL || k == NULL){
        return out;
    }
    
    if(out == NULL){
        out = (__ImageRGB888_t*)__malloc(sizeof(__ImageRGB888_t));
        if(out == NULL) // Not enough space :-(
            return out;
        out->pBuffer = (__UNION_PixelRGB888_t*)__malloc(in->width * in->height * sizeof(__UNION_PixelRGB888_t));
        if(out->pBuffer == NULL) // Not enough space :-(
            return out;
    }
    
    if(out->pBuffer == NULL){
        out->pBuffer = (__UNION_PixelRGB888_t*)__malloc(in->width * in->height * sizeof(__UNION_PixelRGB888_t));
        if(out->pBuffer == NULL) // Not enough space :-(
            return out;
    }
    
    if(out == NULL){
        out = (__ImageRGB888_t*)__malloc(sizeof(__ImageRGB888_t));
        if(out == NULL) // Not enough space :-(
            return out;
        out->pBuffer = (__UNION_PixelRGB888_t*)__malloc(in->width * in->height * sizeof(__UNION_PixelRGB888_t));
    }

    for(int j=0;j<in->height;j++){
        for(int i=0;i<in->width;i++){
            
            unsigned long tmp_R = 0,tmp_G = 0,tmp_B = 0;
            for(int n=0;n<k->order;n++){
                for(int m=0;m<k->order;m++){
                    size_t offset_y = j-(k->order>>1)+n;
                    size_t offset_x = i-(k->order>>1)+m;
                    if(offset_x<0||offset_y<0||offset_x>=in->width||offset_y>=in->height){
                         //... //
                        continue;
                    }
                    unsigned int select_R  = (in->pBuffer + offset_y*in->width + offset_x)->R;
                    unsigned int select_G  = (in->pBuffer + offset_y*in->width + offset_x)->G;
                    unsigned int select_B  = (in->pBuffer + offset_y*in->width + offset_x)->B;
                    int       selectKernel = *( k->pBuffer + n       * k->order + m       );
                    tmp_R += ( (select_R) * (selectKernel) );
                    tmp_G += ( (select_G) * (selectKernel) );
                    tmp_B += ( (select_B) * (selectKernel) );
                }
            }
            size_t offset = (j*in->width)+i;
            if(offset < out->width*out->height){
                (out->pBuffer+offset)->R = tmp_R/coe;//(coe==0)?((1<<8)-1):tmp_R/coe;
                (out->pBuffer+offset)->G = tmp_G/coe;//(coe==0)?((1<<8)-1):tmp_G/coe;
                (out->pBuffer+offset)->B = tmp_B/coe;//(coe==0)?((1<<8)-1):tmp_B/coe;
            }
        }
    }
    
    return out;
}

/*=========================================
 > Memory Programming Reference 
==========================================*/

#define __ADR_DISTANCE( ptr1 , ptr2 )   (size_t)((__abs(ptr2 - ptr1)) - 1)

#pragma pack(1)
unsigned char __VERTUAL_HEAP[ __VIRTUAL_HEAP_SIZE_BYTE ];//__attribute__((at()));
static size_t __Allocated_Bytes__  = 0;

struct __MallocNode_t{
    unsigned long            index;
    size_t                   byte;
    void*                    ptr;
    struct __MallocNode_t*   pNextNode;
}*pHeapMemoryHeadNode = NULL;

/*--------------------------------------------------------------------------------------------------------
 * Memory Node Should be odered by the member of index
 *
 *
 *    Node1     ->>   Node2       ->>     Node3             ->>      Node4       Nodes
 *    Memory1   ->>   Memory2     ->>     Memory3           ->>      Memory4     Used Memory
 *       |               |                   |                          |
 *       |               |                   |                          |
 * [ xxxxxxxxx________xxxxxxx_______xxxxxxxxxxxxxxxxxxx____________xxxxxxxxxxx________ ] Virtual Heap
 * index=0                                                                   index=32768
 *
 --------------------------------------------------------------------------------------------------------*/

void* __mallocHEAP(size_t size){
    size_t size_need       = size;
    if( __Allocated_Bytes__ + size_need > __VIRTUAL_HEAP_SIZE_BYTE )
        return NULL;
    else{
        __Allocated_Bytes__ += size_need;
        // It doesn't mean there is enough space to allocate.
    }
    
    void* ptr = NULL;
    struct __MallocNode_t* pNode     = pHeapMemoryHeadNode;
    struct __MallocNode_t* pNewNode  = (struct __MallocNode_t*)malloc(sizeof(struct __MallocNode_t));
    struct __MallocNode_t* pForeward = NULL,*pBackward = NULL;
    size_t minDist                   = __VIRTUAL_HEAP_SIZE_BYTE;
    
    pNewNode->byte      = size_need;
    pNewNode->pNextNode = NULL;
    
    // Only for test.
    for(int i=0;i<__VIRTUAL_HEAP_SIZE_BYTE;i++)
        __VERTUAL_HEAP[i] = i;
    
    // Special Condition. There isn't any allocated memory.
    if(pNode == NULL){
        pHeapMemoryHeadNode = pNewNode;
        pNewNode->index     = 0;
        ptr                 = &__VERTUAL_HEAP[pNewNode->index];
        return ptr;
    }
    
    // Search the optimal memory block for users.
    while(pNode != NULL){
        size_t size_free = 0;
        // All nodes should be ordered by the member of "index". Which means...
        // "pNode->index" is always ahead of "pNextNode->index" or there will be a problem.
        if(pNode->pNextNode != NULL){
            size_free = (pNode->pNextNode->index) - (pNode->index + pNode->byte);
        }else{
            size_free = (__VIRTUAL_HEAP_SIZE_BYTE-1) - ((pNode->index) + (pNode->byte));
        }
        if( size_free - size_need < minDist && size_free >= size_need ){
            minDist             = size_free - size_need;
            ptr                 = &__VERTUAL_HEAP[ (pNode->index + pNode->byte) ];

            pForeward           = pNode;
            pBackward           = pNode->pNextNode;
            pNewNode->index     = (pForeward->index + pForeward->byte);
            pNewNode->pNextNode = pBackward;
            pNewNode->ptr       = ptr;
        }
        // Continue to search...
        pNode = pNode->pNextNode;
    }
    
    if(ptr != NULL && pForeward != NULL && pNewNode != NULL){
        // Found enough space to allocate
        pForeward->pNextNode = pNewNode;
        pNewNode->pNextNode  = pBackward;
    }else{
        // Fail to find enough space to allocate
        free(pNewNode);
        __Allocated_Bytes__ -= size_need;
    }
    return ptr;
}

void __freeHEAP(void* ptr){
    unsigned long index = (unsigned long)((unsigned char*)ptr - __VERTUAL_HEAP);
    struct __MallocNode_t* pNode     = pHeapMemoryHeadNode;
    struct __MallocNode_t* pForeward = NULL;
    while(pNode != NULL){
        if(pNode->index == index && pNode->ptr == ptr){
            if(pForeward != NULL){
                pForeward->pNextNode = pNode->pNextNode;
                __Allocated_Bytes__ -= pNode->byte;
                free(pNode);
            }
            break;
        }
        pForeward = pNode;
        pNode     = pNode->pNextNode;
    }
}

void* __memsetWORD(void* __b,int value,size_t num){
    uint16_t* src = (uint16_t*)__b;
    size_t remain = num;
    (*((uint16_t*)src)) = value;
    remain--;
    while(1){
        if(num<(remain<<1)){
            memmove((src+(num-remain)),src, (num-remain)*sizeof(uint16_t));
            remain-=(num-remain);
        }else{
            memmove((src+(num-remain)),src, remain*sizeof(uint16_t));
            break;
        }
    }
    return __b;
}

void* __memsetDWORD(void* __b,long value,size_t num){
    uint32_t* src = (uint32_t*)__b;
    size_t remain = num;
    (*((uint32_t*)src)) = (uint32_t)value;
    remain--;
    while(1){
        if(num<(remain<<1)){
            memmove((src+(num-remain)),src, (num-remain)*sizeof(uint32_t));
            remain-=(num-remain);
        }else{
            memmove((src+(num-remain)),src, remain*sizeof(uint32_t));
            break;
        }
    }
    return __b;
}



