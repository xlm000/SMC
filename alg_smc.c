/**
 * @file alg_smc.c
 * @brief 
 * @note 尚不完整，目前主要完善指数趋近律下的推导
 * @note 物理模型固定为Jθ¨=u/k+d(t), J为转动惯量，θ¨为角加速度，u/k为扭矩，若电机直接输出扭矩则k为1
 * @note 扭矩单位采用kg*m²/ms²（计算频率为1000Hz的正常情况下）！！！
 * @note 仅采用一阶滑模面，s=θ˙+cθ，c为滑模面系数
 * @note 使用弧度制！！！
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
    smc_instance->frequency = smc_init->frequency;

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

    smc_instance->error_eps = smc_init->error_eps;
    smc_instance->phi = smc_init->phi;
    return smc_instance;
}

/**
 * @brief sign符号函数
 * @note 仅作保留，以备需要
 */
static int8_t smc_sgn(float s)
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
static float smc_sat(float s, float phi)
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
    case SMC_REACH_EXP:   rho = k + smc_instance->eps * fabsf(s); break;
    case SMC_REACH_POWER: rho = k * pow(fabsf(s), smc_instance->alpha); break;
    default: return 0;
    }
    switch (smc_instance->switch_type) {                     /* 2. 切换函数 -> u_sw */
    case SMC_SWITCH_SIGN:    return -rho * ((s > 0) - (s < 0));
    case SMC_SWITCH_SAT:     return -rho * smc_sat(s, smc_instance->phi);
    case SMC_SWITCH_TANH:    return -rho * tanh(s / smc_instance->phi);
    case SMC_SWITCH_SIGMOID: return -rho * (s / (fabsf(s) + smc_instance->phi));
    }
    return 0;
}

/**
 * @brief 滑模控制器计算
 * @param smc_instance 滑模控制器实例
 * @param actual_angle 当前角度
 * @param actual_angle_speed 当前角速度(单位rpm)
 * @param target_angle 目标角度
 * @return 控制输入 u
 * @note 循环中调用此函数
 */
float smc_tick_calculate(smc_instance_t *smc_instance, float actual_angle, float actual_angle_speed, float target_angle)
{
    //角速度单位转换为rad/ms
    actual_angle_speed = actual_angle_speed * 2*3.1415926f/60.0f;
    //实际值更新
    smc_instance->actual_angle_last = smc_instance->actual_angle;
    smc_instance->actual_angle = actual_angle;
    //目标值更新
    smc_instance->target_angle_last = smc_instance->target_angle;
    smc_instance->target_angle = target_angle;
    smc_instance->target_angle_ddot = (smc_instance->target_angle - smc_instance->target_angle_last) - smc_instance->target_angle_dot;
    smc_instance->target_angle_dot = smc_instance->target_angle - smc_instance->target_angle_last;
    //最短路径规划
    if (smc_instance->actual_angle - smc_instance->target_angle > 3.1415926f)
        smc_instance->actual_angle -= 6.2831852f;
    else if (smc_instance->actual_angle - smc_instance->target_angle < -3.1415926f)
        smc_instance->actual_angle += 6.2831852f;
    //过零保护
    if (smc_instance->actual_angle - smc_instance->actual_angle_last > 3.1415926f)
        smc_instance->actual_angle_last += 6.2831852f;
    else if (smc_instance->actual_angle - smc_instance->actual_angle_last < -3.1415926f)
        smc_instance->actual_angle_last -= 6.2831852f;
    smc_instance->actual_angle_dot = smc_instance->actual_angle - smc_instance->actual_angle_last;

    //速度检查，以免差分值与实际值相差过大影响准确度(若偏差大于50%则使用实际值)
    if (fabsf(smc_instance->actual_angle_dot*smc_instance->frequency - actual_angle_speed) > 0.5*fabsf(actual_angle_speed))
        smc_instance->actual_angle_dot = actual_angle_speed/smc_instance->frequency;

    smc_instance->error_last = smc_instance->error;
    smc_instance->error = smc_instance->actual_angle - smc_instance->target_angle;
    smc_instance->error_dot = smc_instance->actual_angle_dot - smc_instance->target_angle_dot;

    if (smc_instance->error < 0)
		smc_instance->error_qp = -pow(fabsf(smc_instance->error), smc_instance->qp);
	else
		smc_instance->error_qp = pow(fabsf(smc_instance->error), smc_instance->qp);
    
    smc_instance->s = smc_instance->error_dot + smc_instance->c * smc_instance->error_qp;
    smc_instance->ds = -smc_instance->k * smc_sat(smc_instance->s, smc_instance->phi) - smc_instance->eps * smc_instance->s;
    smc_instance->u = (smc_instance->target_angle_ddot - 
                       smc_instance->c * smc_instance->qp * smc_instance->error_dot * smc_instance->error_qp / smc_instance->error +
                       smc_instance->ds) * smc_instance->J * smc_instance->k_c;

    //奇异点附近使用普通滑模
    if (fabsf(smc_instance->error) < 0.05) {
        smc_instance->s = smc_instance->error_dot + smc_instance->c * smc_instance->error;
        smc_instance->ds = - smc_instance->k * smc_sat(smc_instance->s, smc_instance->phi) - smc_instance->eps * smc_instance->s;
        smc_instance->u = (smc_instance->target_angle_ddot - 
                            smc_instance->c * smc_instance->error_dot +
                            smc_instance->ds) * smc_instance->J * smc_instance->k_c;
    }

    if (fabsf(smc_instance->error) < smc_instance->error_eps) {
        smc_instance->u = 0;
        return smc_instance->u;
    }

    if (smc_instance->u > smc_instance->u_max)
        smc_instance->u = smc_instance->u_max;
    else if (smc_instance->u < -smc_instance->u_max)
        smc_instance->u = -smc_instance->u_max; 

    return smc_instance->u;
}