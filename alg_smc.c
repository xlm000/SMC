/**
 * @file alg_smc.c
 * @brief 
 * @note 尚不完整，目前主要完善指数趋近律下的推导
 * @note 物理模型固定为Jθ¨=ku+d(t), J为转动惯量，θ¨为角加速度，ku为扭矩，若电机直接输出扭矩则k为1
 * @note 仅采用一阶滑模面，s=θ˙+cθ，c为滑模面系数
*/


#include "alg_smc.h"
#include "robot_config.h"
#include <math.h>
#include <string.h>

/**
 * @brief 初始化注册
 */
smc_instance_t *smc_register(smc_init_t *smc_init)
{
    smc_instance_t *smc_instance = (smc_instance_t *)user_malloc(sizeof(smc_instance_t));
    if (smc_instance == NULL)
        return NULL;
    memset(smc_instance, 0, sizeof(smc_instance_t));
    smc_instance->J = smc_init->J;
    smc_instance->k_c = smc_init->k_c;
    smc_instance->u_max = smc_init->u_max;

    smc_instance->c = smc_init->c;

    smc_instance->k = smc_init->k;
    smc_instance->reach_type = smc_init->reach_type;
    smc_instance->switch_type = smc_init->switch_type;
    if (smc_instance->reach_type == SMC_REACH_EXP)
        smc_instance->eps = smc_init->eps;
    if (smc_instance->reach_type == SMC_REACH_POWER)
        smc_instance->alpha = smc_init->alpha;
    
    smc_instance->qp = smc_init->q / smc_init->p;

    smc_instance->phi = smc_init->phi;
    return smc_instance;
}

/**
 * @brief sign符号函数
 * @note 仅作保留，以备需要
 */
int8_t sgn(float s)
{
    if (s > 0)
		return 1;
	else if (s == 0)
		return 0;
	else
		return -1;
}

/**
 * @brief sat边界函数
 * @param phi 边界层（一般应设置为正数）
 * @note  
 */
float sat(float s, float phi)
{
    if (phi <= 0.0) 
        return (s > 0) - (s < 0);
    if (s >  phi) 
        return 1;
    if (s < -phi) 
        return -1;
    return s / phi;
}

/**
 * @brief 切换控制律计算
 * @note 仅作保留，以备需要
 */
float smc_switch_control(const smc_instance_t *smc_instance, float s)
{
    float k = smc_instance->k;
    float rho;
    switch (smc_instance->reach_type) {                  /* 1. 趋近律 -> 趋近速度 rho */
    case SMC_REACH_CONST: rho = k; break;
    case SMC_REACH_EXP:   rho = k + smc_instance->eps * fabs(s); break;
    case SMC_REACH_POWER: rho = k * pow(fabs(s), smc_instance->alpha); break;
    default: return 0;
    }
    switch (smc_instance->switch_type) {                     /* 2. 切换函数 -> u_sw */
    case SMC_SWITCH_SIGN:    return -rho * ((s > 0) - (s < 0));
    case SMC_SWITCH_SAT:     return -rho * sat(s, smc_instance->phi);
    case SMC_SWITCH_TANH:    return -rho * tanh(s / smc_instance->phi);
    case SMC_SWITCH_SIGMOID: return -rho * (s / (fabs(s) + smc_instance->phi));
    }
    return 0;
}

/**
 * @brief 滑模控制器计算
 * @param smc_instance 滑模控制器实例
 * @param angle_now 当前角度
 * @param target_angle 目标角度
 * @return 控制输入 u
 * @note 循环中调用此函数
 */
float smc_tick_calculate(smc_instance_t *smc_instance, float angle_now, float target_angle)
{
    smc_instance->angle_last = smc_instance->angle_now;
    smc_instance->angle_now = angle_now;
    smc_instance->angle_dot = smc_instance->angle_now - smc_instance->angle_last;

    smc_instance->target_angle_last = smc_instance->target_angle;
    smc_instance->target_angle = target_angle;
    smc_instance->target_angle_ddot = (smc_instance->target_angle - smc_instance->target_angle_last) - smc_instance->target_angle_dot;
    smc_instance->target_angle_dot = smc_instance->target_angle - smc_instance->target_angle_last;


    smc_instance->error_last = smc_instance->error;
    smc_instance->error = smc_instance->angle_now - smc_instance->target_angle;
    smc_instance->error_dot = smc_instance->angle_dot - smc_instance->target_angle_dot;

    if (fabsf(smc_instance->error) < smc_instance->error_eps) {
        smc_instance->u = 0;
        return smc_instance->u;
    }

    float error_qp;
    if (smc_instance->error < 0)
		error_qp = -pow(fabsf(smc_instance->error), smc_instance->qp);
	else
		error_qp = pow(fabsf(smc_instance->error), smc_instance->qp);
    
    smc_instance->s = smc_instance->error_dot + smc_instance->c * error_qp;
    //趋近律
    smc_instance->ds = -smc_instance->k * sat(smc_instance->s, smc_instance->phi) - smc_instance->eps * smc_instance->s;
    smc_instance->u = (smc_instance->target_angle_ddot - 
                       smc_instance->c * smc_instance->qp * smc_instance->error_dot * error_qp / smc_instance->error +
                       smc_instance->ds) * smc_instance->J * smc_instance->k_c;

    //奇异点附近使用普通滑模
    if (abs(smc_instance->error) < 1) {
        smc_instance->s = smc_instance->error_dot + smc_instance->c * smc_instance->error;
        smc_instance->ds = -smc_instance->k * sat(smc_instance->s, smc_instance->phi) - smc_instance->eps * smc_instance->s;
        smc_instance->u = (smc_instance->target_angle_ddot - 
                           smc_instance->c *  smc_instance->error_dot +
                           smc_instance->ds) * smc_instance->J * smc_instance->k_c;
    }

    if (smc_instance->u > smc_instance->u_max)
        smc_instance->u = smc_instance->u_max;
    else if (smc_instance->u < -smc_instance->u_max)
        smc_instance->u = -smc_instance->u_max; 

    return smc_instance->u;
}