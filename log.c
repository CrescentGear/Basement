#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include<stdio.h>
#include "RH_Utility.h"
#ifdef __cplusplus
//extern "C" {
#endif

#define N 100
#define M 100



void Log(char message[N][M],char s[], ...)
{
	char *str = NULL;
	va_list sArgv;          // 申请参数列表变量
	va_start(sArgv, s); // 申明最后一个传递给函数的已知的固定参数

	int i = 0;
	for (i = 0; message[i][0] != 0; i++);
	strcpy(message[i], s);                 //锁定数组的可用位置

	str = va_arg(sArgv, char*);
	while (strcmp(str,"-1"))
	{
		for (i = 0; message[i][0] != 0; i++);
		strcpy(message[i], str);               //把数据存入数组
		str = va_arg(sArgv, char*);
	}
	va_end(sArgv);
}

void main()
{
	char s[M];
	scanf("%s",s);
	const char message[N][M] = { 0 };
	sprintf(message, "first");
	Log(message,s,"second","third","-1");      //"-1"代表结束输入
	for (int i = 0; message[i][0] != 0; i++)
	puts(message[i]);
}


