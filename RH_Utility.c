#include "RH_Utility.h"
#include <string.h>
#include <stdint.h>
#include <math.h>

/*=========================================
 > Algebra Reference 
==========================================*/

int __sign(int x){
    return (x>=0)?(1):(-1);
}

int __sqrt(int x){
    if(x <= 0) return 0;
    int l   = 1;
    int r   = x;
    int res = 0;
    while(l <= r){
        int mid=(l+r)>>1;
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
    int dist_AB = lroundl(sqrt( (B->x - A->x)*(B->x - A->x) + \
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
#if 1//def _WINDOWS_

__ImageRGB888_t __LoadBMP_ImgRGB888(const char* __restrict__ path){
    FILE* bmp;
    BITMAPFILEHEADER fileHead;
    BITMAPINFOHEADER infoHead;
    __ImageRGB888_t  IMG = {
        .width   = 0,
        .height  = 0,
        .pBuffer = NULL
    };
    __PixelRGB888_t* pBuffer;

    bmp = fopen(path, "r");
    if (bmp == NULL) {
        printf("open error\n");
        return IMG;
    }
    fread(&fileHead, sizeof(BITMAPFILEHEADER), 1, bmp);
    fread(&infoHead, sizeof(BITMAPINFOHEADER), 1, bmp);

    if (fileHead.bfType != 0x4D42) {
        //This not a *.bmp file
        printf("not a *.bmp file\n");
        return IMG;
    }

    fseek(bmp, fileHead.bfOffBits, SEEK_SET);
    pBuffer = (__PixelRGB888_t*)malloc(infoHead.biWidth * infoHead.biHeight * sizeof(__PixelRGB888_t));
    fread(pBuffer, sizeof(__PixelRGB888_t), infoHead.biWidth * infoHead.biHeight, bmp);
    fclose(bmp);

    for (int row = 0; row < infoHead.biHeight; row++) {
        for (int col = 0; col < infoHead.biWidth; col++) {
            uint8_t R = pBuffer[row * infoHead.biWidth + col].R;
            uint8_t B = pBuffer[row * infoHead.biWidth + col].B;
            pBuffer[row * infoHead.biWidth + col].R = B;
            pBuffer[row * infoHead.biWidth + col].B = R;
        }
    }

    IMG.width   = infoHead.biWidth;
    IMG.height  = infoHead.biHeight;
    IMG.pBuffer = pBuffer;
    printf("w = %d,h = %d\n",IMG.width,IMG.height);
    printf("Success\n");
    return IMG;
}

void __OutBMP_ImgRGB888(const char* __restrict__ path,__ImageRGB888_t* p){
    FILE* bmp;
    BITMAPFILEHEADER fileHead = {
        .bfOffBits   = 14 + 40 + ...,
        .bfReserved1 = 0,
        .bfReserved2 = 0,
        .bfSize      = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+...+p->height*p->width,
        .bfType      = 0x4D42, 
    };
    BITMAPINFOHEADER infoHead = {
        .
    };
}

void __FreeBMP_ImgRGB888(__ImageRGB888_t* p){
    free(p->pBuffer);
    p->height  = 0;
    p->width   = 0;
    p->pBuffer = NULL;
}

#else

#endif

void __Conv2D_ImgRGB565(const __ImageRGB565_t* in,const __Kernel_t* k,__ImageRGB565_t* out,int coe){
    if( in == NULL || k == NULL || out == NULL)
        return;

    for(int j=0;j<in->height;j++){
        for(int i=0;i<in->width;i++){
            unsigned long tmp_R = 0,tmp_G = 0,tmp_B = 0;
            for(int n=0;n<k->order;n++){
                for(int m=0;m<k->order;m++){
                    int offset_y = j-(k->order>>1)+n;
                    int offset_x = i-(k->order>>1)+m;
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
            int offset = (j*in->width)+i;
            if(offset < out->width*out->height){
                (out->pBuffer+offset)->R = (coe==0)?(1<<5-1):tmp_R/coe;
                (out->pBuffer+offset)->G = (coe==0)?(1<<6-1):tmp_G/coe;
                (out->pBuffer+offset)->B = (coe==0)?(1<<5-1):tmp_B/coe;
            }
        }
    }
}

void __Conv2D_ImgRGB888(const __ImageRGB888_t* in,const __Kernel_t* k,__ImageRGB888_t* out,int coe){
    if( in == NULL || k == NULL || out == NULL)
        return;

    for(int j=0;j<in->height;j++){
        for(int i=0;i<in->width;i++){
            unsigned long tmp_R = 0,tmp_G = 0,tmp_B = 0;
            for(int n=0;n<k->order;n++){
                for(int m=0;m<k->order;m++){
                    int offset_y = j-(k->order>>1)+n;
                    int offset_x = i-(k->order>>1)+m;
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
            int offset = (j*in->width)+i;
            if(offset < out->width*out->height){
                (out->pBuffer+offset)->R = (coe==0)?(1<<8-1):tmp_R/coe;
                (out->pBuffer+offset)->G = (coe==0)?(1<<8-1):tmp_G/coe;
                (out->pBuffer+offset)->B = (coe==0)?(1<<8-1):tmp_B/coe;
            }
        }
    }
}

/*=========================================
 > Memory Programming Reference 
==========================================*/
 // unsigned char __VERTUAL_HEAP[ __VIRTUAL_HEAP_SIZE_BYTE ]; //__attribute__((at()));

 // void* __mallocHEAP(size_t size){
 //     static size_t byteCNT = 0;
 //     if( byteCNT+size > __VIRTUAL_HEAP_SIZE_BYTE )
 //         return NULL;
 //     else{
 //         byteCNT+=size;
 //          //... //
 //         return NULL;
 //     }
 // }

 // void __freeHEAP(void* ptr){

 // }

void* __memsetWORD(void* __b,WORD value,size_t num){
    WORD* src = (WORD*)__b;
    size_t remain = num;
    (*((WORD*)src)) = value;
    remain--;
    while(1){
        if(num<(remain<<1)){
            memmove((src+(num-remain)),src, (num-remain)*sizeof(WORD));
            remain-=(num-remain);
        }else{
            memmove((src+(num-remain)),src, remain*sizeof(WORD));
            break;
        }
    }
    return __b;
}

#if 1
#include <windows.h>
int main(int argc, char const *argv[]){
    const char* __restrict__ src  = "D:/Personal/Desktop/lenna.bmp";
    const char* __restrict__ des  = "D:/Personal/Desktop/output.bmp";
    __ImageRGB888_t IMG = __LoadBMP_ImgRGB888(src);
    __OutBMP_ImgRGB888(des,&IMG);
    __FreeBMP_ImgRGB888(&IMG);
}
#endif


