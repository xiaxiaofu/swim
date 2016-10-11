//  Specification：游泳暂停，转向，划水次数检测
//  Reference：by xxf
//  Author   ：xxf
//  Time     ：2016-09-10
//  Copyright：YF
//  需要改善：滤波算法，谷峰检测算法，次数与圈数的不同步，2016-09-12
//  缺   点：对泳池中途换泳姿计次不是很理想，圈数延迟问题（大概延迟半圈）
#include "Swimdetect.h"
#include <math.h>
#include <stdio.h>
//参数定义
#define pi 3.1415926
#define Fs 25       //sample frequence
//本文件全局变量
static int lg = Fs * 0.8;//陀螺仪窗长度0.4s
static int lm = Fs / 5 * 8;//磁力计窗长度8s
static float fifo4[4] = {0.0, -20.0, -40.0, 1.0};
static float q[4] = {1.0, 0.0, 0.0, 0.0};
static float fifot[3] = {.0};
static float fifoi[3] = {.0};
static float fifom[3] = {.0};
static float fifoacosq1[(int)(Fs*0.8)] = {.0};//20
static float fifoacosq2[(int)(Fs*0.8)] = {.0};//20
static float fifonm1[(int)(Fs/5*8)] = {.0};//40
static float fifonm2[(int)(Fs/5*8)] = {.0};//40 8s windows
static float KFg[2] = {.0};
static int flag = 1;
static int step = 2;
static int stroke = 0;
static int lap = 0;
static int lapt = 0;
static int c5 = 5;
static float KFmt[2] = {0, -5};
static float KFmn[2] = {.0};
static float maxKFmn = 150;
static float maxdet = 0;
static int cstatic = 0;
//辅助函数定义
static float swim_filter(float *fifo, float x, int l);
//函数引用
extern int report_change(int type, unsigned int value);
//检测函数
void swim_stroke_lap(float t, float gx, float gy, float gz, float mx, float my, float mz){
    //局部变量定义
    float det, s1, s2, s3, s4, s5, tq[4], nm, diffKFmn, dm;
    float gx1 = gx * (0.5 * (1.0 / Fs));
    float gy1 = gy * (0.5 * (1.0 / Fs));
    float gz1 = gz * (0.5 * (1.0 / Fs));
    float A[4][4] = {{0, -gx1, -gy1, -gz1}, {gx1, 0, gz1, -gy1}, {gy1, -gz1, 0, gx1}, {gz1, gy1, -gx1, 0}};
    tq[0] = q[0]; tq[1] = q[1]; tq[2] = q[2]; tq[3] = q[3];
    //----------------------count strokes-----------------------
    q[0] = q[0] + A[0][0]*tq[0] + A[0][1]*tq[1] + A[0][2]*tq[2] + A[0][3]*tq[3];
    q[1] = q[1] + A[1][0]*tq[0] + A[1][1]*tq[1] + A[1][2]*tq[2] + A[1][3]*tq[3];
    q[2] = q[2] + A[2][0]*tq[0] + A[2][1]*tq[1] + A[2][2]*tq[2] + A[2][3]*tq[3];
    q[3] = q[3] + A[3][0]*tq[0] + A[3][1]*tq[1] + A[3][2]*tq[2] + A[3][3]*tq[3];
    float tempsqrt = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    q[0] = q[0] / tempsqrt;
    q[1] = q[1] / tempsqrt;
    q[2] = q[2] / tempsqrt;
    q[3] = q[3] / tempsqrt;
    float acosq = 2 * acosf(q[0]) * 180 / pi;
    float acosq1 = swim_filter(fifoacosq1, acosq, lg);
    float acosq2 = swim_filter(fifoacosq2, acosq1, lg);
    for(int i = 0; i < 2; i++){
        fifoi[i] = fifoi[i+1];
    }
    fifoi[2] = acosq2;
    if((fifoi[2] >= fifoi[1] && fifoi[0] - fifoi[1] > 10e-7)||
       (fifoi[2] <= fifoi[1] && fifoi[1] - fifoi[0] > 10e-7)){
        KFg[0] = KFg[1]; KFg[1] = fifoi[1];
        det = fabsf(KFg[1] - KFg[0]);
        cstatic = cstatic + 1;
        if ((det > 15) || (det/fifo4[3] > 0.6 && det > 8 && det/maxdet > 0.1) ||
            (cstatic >= 4 && det > 8 && det/maxdet > 0.1)){
            cstatic  = 0;
            for(int i = 0; i < 3; i++){
                fifo4[i] = fifo4[i+1];
            }
            fifo4[3] = det;
            s1 = fabsf(fifo4[0] - fifo4[1]);
            s2 = fabsf(fifo4[2] - fifo4[3]);
            s3 = 2*s1/(fifo4[0] + fifo4[1]);
            s4 = 2*s2/(fifo4[2] + fifo4[3]);
            s5 = fabsf(fifo4[0] - fifo4[2]);
            if (step == 2){
                if (((s1 < 11 && s2 < 11) || (s3 < 0.18 && s4 < 0.18)) && (s5 < 48)){
                    if (flag == 1){
                        lap = lapt;
                        stroke = stroke + 3;
                        flag = 2;
                        report_change(STROKES,3);
                    }
                    else{
                        stroke = stroke + 1;
                        report_change(STROKES,1);
                    }
                    if (det > maxdet){
                        maxdet = det;
                    }
                    step = 1;
                }
                else{
                    flag = 1;
                }
            }
            else{
                step = step + 1;
            }
        }
        
    }
    //-----------------------count laps------------------------
    if(c5 == 5){//降频
        c5 = 1;
        nm = sqrt(mx*mx + my*my + mz*mz);
        float fnm1 = swim_filter(fifonm1, my, lm);
        float fnm2 = swim_filter(fifonm2, fnm1, lm);//second filter
        for(int i = 0; i < 2; i++){
            fifot[i] = fifot[i+1];
            fifom[i] = fifom[i+1];
        }
        fifot[2] = t;//no use
        fifom[2] = fnm2;//second filter        
        if(((fifom[2] >= fifom[1]) && (fifom[0] - fifom[1] > 10e-7)) ||
           ((fifom[2] <= fifom[1]) && (fifom[1] - fifom[0] > 10e-7))){//谷
            //KFmt:use to record 2 time
            //KFmn:use to record 2 norm of mag
            if(fifot[1] - KFmt[1] >= 5  && stroke > 0){
                KFmt[0] = KFmt[1];  KFmt[1] = fifot[1]; //record 3 Ku time
                KFmn[0] = KFmn[1];  KFmn[1] = fifom[1];//记录所有磁力计滤波的谷值
                diffKFmn = fabsf(KFmn[1] - KFmn[0]);
                if(diffKFmn > maxKFmn){
                    maxKFmn = diffKFmn;
                }
                dm = diffKFmn/maxKFmn;
                if((diffKFmn > 170 && dm > 0.55) || (dm > 0.9 && lapt == 0) ||  diffKFmn > 200){
                    lapt = lapt + 1;
                    report_change(LAPS,1);
                }
                else if(lapt == 0){
                    lapt = lapt + 1;
                    report_change(LAPS,1);
                }
                
//                printf("t:%f strokes：%d, laps:%d\n", KFmt[1], stroke, lapt);
            }
            else{
                //printf("DO NOTHING!\n");
            }
        }
    }
    else{
        c5 = c5 + 1;
    }
    
//    printf("t:%f strokes：%d, laps:%d\n", KFmt[1], stroke, lap);
    
    
}

float swim_filter(float *fifo, float x, int l){
    float mdw = .0;
    for(int i = 0; i < l-1; i++){
        fifo[i] = fifo[i+1];
        mdw = mdw + fifo[i];
    }
    fifo[l-1] = x;
    return (mdw + x)/l;
}

