
/*=================================================================================================
 > This part of code will never be compiled.
===================================================================================================*/
#if (1==0)
// Dummy Segment for Gussian Kernel

    int      kernel[1+2+3+4+5] = {0};
    uint16_t gus_kernel[9][9];

    int sigma = 0;

    double temp = 12.8;

    kernel[0]  = (uint16_t)(100 / (M_PI * temp));                       //[0][0]

    kernel[1]  = (uint16_t)(100 / (M_PI * temp) * exp(-1  / (temp)));   //[1][0]
    kernel[2]  = (uint16_t)(100 / (M_PI * temp) * exp(-2  / (temp)));   //[1][1]

    kernel[3]  = (uint16_t)(100 / (M_PI * temp) * exp(-4  / (temp)));   //[2][0]
    kernel[4]  = (uint16_t)(100 / (M_PI * temp) * exp(-5  / (temp)));   //[2][1]
    kernel[5]  = (uint16_t)(100 / (M_PI * temp) * exp(-8  / (temp)));   //[2][2]

    kernel[6]  = (uint16_t)(100 / (M_PI * temp) * exp(-9  / (temp)));   //[3][0]
    kernel[7]  = (uint16_t)(100 / (M_PI * temp) * exp(-10 / (temp)));   //[3][1]
    kernel[8]  = (uint16_t)(100 / (M_PI * temp) * exp(-13 / (temp)));   //[3][2]
    kernel[9]  = (uint16_t)(100 / (M_PI * temp) * exp(-18 / (temp)));   //[3][3]

    kernel[10] = (uint16_t)(100 / (M_PI * temp) * exp(-16 / (temp)));   //[4][0]
    kernel[11] = (uint16_t)(100 / (M_PI * temp) * exp(-17 / (temp)));   //[4][1]
    kernel[12] = (uint16_t)(100 / (M_PI * temp) * exp(-20 / (temp)));   //[4][2]
    kernel[13] = (uint16_t)(100 / (M_PI * temp) * exp(-25 / (temp)));   //[4][3]
    kernel[14] = (uint16_t)(100 / (M_PI * temp) * exp(-32 / (temp)));   //[4][4]

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


#endif




