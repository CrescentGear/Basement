

#include "RH_Utility.h"

int main(int argc, char const *argv[]){
    // const char* __restrict__ src  = "C:\\Users\\asus\\Desktop\\lenna.bmp";
    // const char* __restrict__ des  = "C:\\Users\\asus\\Desktop\\output.bmp";

    const char* __restrict__ src  = "/Users/randle_h/desktop/lenna.bmp";
    const char* __restrict__ des  = "/Users/randle_h/desktop/output.bmp";

    __ImageRGB888_t IMG = __LoadBMP_ImgRGB888(src);
    __OutBMP_ImgRGB888(des,&IMG);
    __FreeBMP_ImgRGB888(&IMG);

}

