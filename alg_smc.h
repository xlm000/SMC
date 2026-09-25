/**
 * @file alg_smc.h
 * @note 尚不完整，目前主要完善指数趋近律下的推导
 * @note 物理模型固定为Jθ¨=ku+d(t), J为转动惯量，θ¨为角加速度，ku为扭矩，若电机直接输出扭矩则k为1
 * @note 仅采用一阶滑模面，s=θ˙+cθ，c为滑模面系数
*/

#ifndef _ALG_SMC_H_
#define _ALG_SMC_H_

#define DT 0.001

#include "stdint.h"
#include <stdbool.h>

/* 切换函数(边界层)类型 */
typedef enum {
    SMC_SWITCH_SIGN    = 0,  /* 符号函数: 理想滑模, 抖振大 */
    SMC_SWITCH_SAT     = 1,  /* 饱和函数: 边界层法, 工程常用 */
    SMC_SWITCH_TANH    = 2,  /* 双曲正切: 连续平滑 */
    SMC_SWITCH_SIGMOID = 3   /* Sigmoid:  连续平滑 */
} smc_switch_t;

/* 趋近律类型 */
typedef enum {
    SMC_REACH_CONST = 0,  /* 等速:   s_dot = -k*sat(s)          */
    SMC_REACH_EXP   = 1,  /* 指数:   s_dot = -k*sat(s) - eps*s  */
    SMC_REACH_POWER = 2   /* 幂次:   s_dot = -k*|s|^a*sat(s)     */
} smc_reach_t;

typedef struct {
    float frequency;           // 控制频率(计算频率，一般freertos为1000Hz)(暂无用)

    float J;                   // 转动惯量
    float k_c;                 // 控制输入系数，k_c*u为扭矩，若电机直接输出扭矩则k为1
    float u_max;               // 控制输入上限

    float c;                   // 滑模面系数

    float k;                   // 切换增益，需大于扰动上界 
    float eps;                 // 指数趋近系数
    float alpha;               // 幂次趋近指数 (0<a<1)
    smc_reach_t   reach_type;  // 趋近律类型
    smc_switch_t  switch_type; // 切换函数类型

    int q;                   // 终端滑模指数参数上
    int p;                   // 终端滑模指数参数下，q、p均为奇数,q<p
    
    float error_eps;           // 死区，若误差小于该值则认为已收敛，输出0
    float phi;                 // sat边界层厚度
} smc_init_t;

typedef struct {
    float frequency;           // 控制频率(计算频率)

    float J;                   // 转动惯量
    float u;                   // 控制输入
    float k_c;                 // 控制输入系数，k_c*u为扭矩，若电机直接输出扭矩则k为1
    float u_max;               // 控制输入上限

    float s;                   // 滑模面
    float c;                   // 滑模面系数

    float ds;                  // 滑模面导数(趋近律)
    float k;                   // 切换增益，需大于扰动上界 
    float eps;                 // 指数趋近系数
    float alpha;               // 幂次趋近指数 (0<a<1)
    smc_reach_t   reach_type;  // 趋近律类型
    smc_switch_t  switch_type; // 切换函数类型

    float qp;                  // 终端滑模指数参数 qp = q/p
    
    float error_eps;           // 死区，若误差小于该值则认为已收敛，输出0
    float phi;                 // sat边界层厚度
    
    float target_angle;
    float target_angle_last;
    float target_angle_dot;
    float target_angle_ddot;

    float actual_angle_now;
    float actual_angle_last;
    float actual_angle_dot;
    float actual_angle_ddot;

    float error;
    float error_last;
    float error_dot;
    float error_qp;
} smc_instance_t;

int8_t sgn(float s);
float sat(float s, float phi);
smc_instance_t *smc_register(smc_init_t *smc_init);
float smc_switch_control(const smc_instance_t *smc_instance, float s);
float smc_tick_calculate(smc_instance_t *smc_instance, float actual_angle, float actual_angle_speed, float target_angle);

#endif