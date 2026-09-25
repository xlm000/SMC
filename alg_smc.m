%% test_smc.m
%  alg_smc(2).c (最新版) 的 MATLAB 复刻版本 —— 纯差分/拍域设计, 弧度制
%  物理模型: J*theta_ddot = k_c*u + d(t)   (theta 单位为 rad)
%  误差定义: error = actual - target (本版约定), error_dot = actual_dot - target_dot
%  滑模面  : s = error_dot + c * error^(q/p)  (|error|>=0.05)
%            s = error_dot + c * error        (|error|< 0.05, 奇异点保护)
%  趋近律  : ds = -k*sat(s,phi) - eps*s       (指数趋近律, 拍域)
%
%  【差分代替微分 / 拍域说明】
%  控制器差分不除 dt, 被控对象按"拍"积分: 1 count = Ts = 1/frequency = 0.001 s。
%  此时 ds 是"滑模面每拍增量", 增益须满足离散稳定条件:
%    边界层内: s[k+1] = (1 - eps - k/phi) s[k]  →  eps + k/phi < 2 (建议 < 1.2)
%    滑模面上: e[k+1] = (1 - c) e[k]            →  0 < c < 2
%    抗扰:    k 大于扰动每拍对 s 的增量上界 ≈ (d_max/J)*(1 + c*qp*error_qp/error)
%
%  【符号约定验证】本版 error = actual - target, 控制律
%    u = (target_ddot - c*qp*error_dot*error_qp/error + ds) * J * k_c
%  与推导 u = (target_ddot - c*qp*error_dot*error_qp/error + ds) * J/k_c
%  符号一致; k_c = 1 时 J*k_c 与 J/k_c 数值相同, 已按 C 代码忠实翻译。
%  ⚠ 但若 k_c ≠ 1, 按模型 J*theta_ddot = k_c*u 应取 J/k_c, 请确认固件是否有意为之。

clear; clc;

%% ---------------- 配置参数(与 smc_init_t 对应, 拍域增益) ----------------
smc_config = struct( ...
    'frequency', 1000, ...  % 控制频率 (Hz), 1 count = 1/frequency s = 0.001 s
    'J',         0.0075,  ...  % 转动惯量
    'k_c',       1,    ...  % 控制输入系数, k_c*u 为扭矩
    'u_max',     8, ... % 控制输入上限
    'c',         7, ...  % 滑模面系数 (每拍误差衰减比例, 0<c<2)
    'k',         12, ...  % 切换增益 (需大于扰动每拍增量上界)
    'eps',       25, ...  % 指数趋近系数 (拍域中必须 < 2!)
    'alpha',     0,    ...  % 幂次趋近指数 (0<a<1), 本例未使用
    'reach_type',  1,  ...  % 趋近律类型: 1 = SMC_REACH_EXP
    'switch_type', 1,  ...  % 切换函数类型: 1 = SMC_SWITCH_SAT
    'q',         23,   ...  % 终端滑模指数参数(上)
    'p',         27,   ...  % 终端滑模指数参数(下), q、p 均为奇数
    'error_eps', 1e-3, ...  % 死区 (rad)
    'phi',       1);        % sat 边界层厚度 (拍域中 k/phi + eps < 2)

Ts = 1 / smc_config.frequency;   % 每拍物理时间 = 0.001 s
max_count = 200000;

%% ---------------- 枚举定义(与 alg_smc.h 对应) ----------------
% 注意: 脚本的局部函数无法访问脚本工作区, 故通过函数 smc_consts() 获取常量
C = smc_consts();
SMC_SWITCH_SIGN    = C.SMC_SWITCH_SIGN;
SMC_SWITCH_SAT     = C.SMC_SWITCH_SAT;
SMC_SWITCH_TANH    = C.SMC_SWITCH_TANH;
SMC_SWITCH_SIGMOID = C.SMC_SWITCH_SIGMOID;
SMC_REACH_CONST = C.SMC_REACH_CONST;
SMC_REACH_EXP   = C.SMC_REACH_EXP;
SMC_REACH_POWER = C.SMC_REACH_POWER;

%% ---------------- 仿真主循环 ----------------
TARGET = pi/3;            % 目标角度 60° (弧度制)
smc = smc_register(smc_config);
% 预填目标状态: 避免第一步目标从 0 突变到 TARGET 产生差分尖峰
smc.target_angle      = TARGET;
smc.target_angle_last = TARGET;

angle_now = 0;            % 实际角度 (rad)
angle_dot = 0;            % 被控对象"每拍角位移"状态(拍域速度)
count     = 0;

angle_log = zeros(max_count, 1);
u_log     = zeros(max_count, 1);

while abs(angle_now - TARGET) > smc_config.error_eps
    count = count + 1;
    if count > max_count
        warning('达到最大迭代次数 %d, 系统未收敛。', max_count);
        count = max_count;
        break;
    end

    [smc, u] = smc_tick_calculate(smc, angle_now, TARGET);
    % 注: C 版接口还有 actual_angle_speed 参数用于速度准确度判定, 按要求忽略

    angle_dot = angle_dot + (u / smc.J) * Ts;   % 拍域牛顿第二定律: 每拍速度增量 = u/J
    angle_dot = angle_dot * 0.99;          % 简单阻尼(每拍 1%, 时间常数约 100 拍 = 0.1 s)
    angle_now = angle_now + angle_dot * Ts;     % 拍域积分

    angle_log(count) = angle_now;
    u_log(count)     = u;

    if mod(count, 50) == 0
        fprintf('Angle: %f rad, U: %f, count: %d (t = %.3f s)\n', angle_now, u, count, count*Ts);
    end
end
fprintf('结束: count = %d (t = %.3f s), angle = %f rad\n', count, count*Ts, angle_now);

%% ---------------- 绘图 ----------------
t = (1:count) * Ts;
figure('Name', 'SMC Test (radian)');
subplot(2,1,1);
plot(t, angle_log(1:count)*180/pi, 'b'); hold on;
yline(TARGET*180/pi, 'r--', 'Target');
xlabel('时间 (s)'); ylabel('\theta (°)'); title('角度响应'); grid on;

subplot(2,1,2);
plot(t, u_log(1:count), 'b');
xlabel('时间 (s)'); ylabel('u'); title('控制输入'); grid on;

%% ========================================================================
%                            本文件内的局部函数
%  ========================================================================

function smc = smc_register(cfg)
% smc_register  与 alg_smc(1).c 中 smc_register 对应: 初始化控制器实例
    C = smc_consts();
    SMC_REACH_EXP = C.SMC_REACH_EXP; SMC_REACH_POWER = C.SMC_REACH_POWER;
    smc.frequency = cfg.frequency;
    smc.J         = cfg.J;
    smc.k_c       = cfg.k_c;
    smc.u_max     = cfg.u_max;
    smc.c         = cfg.c;
    smc.k         = cfg.k;
    smc.reach_type  = cfg.reach_type;
    smc.switch_type = cfg.switch_type;
    if smc.reach_type == SMC_REACH_EXP
        smc.eps = cfg.eps;
    end
    if smc.reach_type == SMC_REACH_POWER
        smc.alpha = cfg.alpha;
    end
    smc.qp        = cfg.q / cfg.p;
    smc.phi       = cfg.phi;
    smc.error_eps = cfg.error_eps;

    % ---- 状态变量清零 ----
    smc.u = 0;  smc.s = 0;  smc.ds = 0;
    smc.actual_angle      = 0;  smc.actual_angle_last = 0;
    smc.actual_angle_dot  = 0;
    smc.target_angle      = 0;  smc.target_angle_last = 0;
    smc.target_angle_dot  = 0;  smc.target_angle_ddot = 0;
    smc.error      = 0;  smc.error_last = 0;  smc.error_dot = 0;
    smc.error_qp   = 0;
end

function y = sat(s, phi)
% sat  饱和边界函数(与 alg_smc.c 中 sat 对应)
    if phi <= 0
        y = double(s > 0) - double(s < 0);
    elseif s >  phi
        y = 1;
    elseif s < -phi
        y = -1;
    else
        y = s / phi;
    end
end

function u_sw = smc_switch_control(smc, s)
% smc_switch_control  切换控制律(与 alg_smc.c 对应, 保留备用)
    C = smc_consts();
    SMC_REACH_CONST = C.SMC_REACH_CONST; SMC_REACH_EXP = C.SMC_REACH_EXP;
    SMC_REACH_POWER = C.SMC_REACH_POWER;
    SMC_SWITCH_SIGN = C.SMC_SWITCH_SIGN; SMC_SWITCH_SAT = C.SMC_SWITCH_SAT;
    SMC_SWITCH_TANH = C.SMC_SWITCH_TANH; SMC_SWITCH_SIGMOID = C.SMC_SWITCH_SIGMOID;
    k   = smc.k;
    rho = 0;
    switch smc.reach_type
        case SMC_REACH_CONST
            rho = k;
        case SMC_REACH_EXP
            rho = k + smc.eps * abs(s);
        case SMC_REACH_POWER
            rho = k * abs(s)^smc.alpha;
        otherwise
            u_sw = 0; return;
    end
    switch smc.switch_type
        case SMC_SWITCH_SIGN
            u_sw = -rho * (double(s > 0) - double(s < 0));
        case SMC_SWITCH_SAT
            u_sw = -rho * sat(s, smc.phi);
        case SMC_SWITCH_TANH
            u_sw = -rho * tanh(s / smc.phi);
        case SMC_SWITCH_SIGMOID
            u_sw = -rho * (s / (abs(s) + smc.phi));
        otherwise
            u_sw = 0;
    end
end

function [smc, u] = smc_tick_calculate(smc, actual_angle, target_angle)
% smc_tick_calculate  滑模控制器单步计算(与 alg_smc(2).c 对应, 纯差分/拍域)
% 每个控制周期调用一次。C 版签名中的 actual_angle_speed(速度准确度判定)按需求忽略。

    % ---- 实际角度状态(差分代替微分, 拍域) ----
    smc.actual_angle_last = smc.actual_angle;
    smc.actual_angle      = actual_angle;
    smc.actual_angle_dot  = smc.actual_angle - smc.actual_angle_last;

    % ---- 目标状态(注意顺序: 先算 ddot 再更新 dot, 与 C 代码一致) ----
    smc.target_angle_last = smc.target_angle;
    smc.target_angle      = target_angle;
    smc.target_angle_ddot = (smc.target_angle - smc.target_angle_last) - smc.target_angle_dot;
    smc.target_angle_dot  = smc.target_angle - smc.target_angle_last;

    % ---- 误差状态: error = actual - target ----
    smc.error_last = smc.error;
    smc.error      = smc.actual_angle - smc.target_angle;
    smc.error_dot  = smc.actual_angle_dot - smc.target_angle_dot;

    % ---- 死区 ----
    if abs(smc.error) < smc.error_eps
        smc.u = 0;
        u = smc.u;
        return;
    end

    % ---- 终端滑模面: s = error_dot + c * error^(q/p) ----
    if smc.error < 0
        smc.error_qp = -abs(smc.error)^smc.qp;
    else
        smc.error_qp =  abs(smc.error)^smc.qp;
    end

    smc.s  = smc.error_dot + smc.c * smc.error_qp;
    smc.ds = -smc.k * sat(smc.s, smc.phi) - smc.eps * smc.s;   % 指数趋近律
    smc.u  = (smc.target_angle_ddot ...
              - smc.c * smc.qp * smc.error_dot * smc.error_qp / smc.error ...
              + smc.ds) * smc.J * smc.k_c;    % 与 alg_smc(2).c 一致

    % ---- 奇异点(|error| < 0.05)退化为一阶线性滑模面 ----
    if abs(smc.error) < 0.05
        smc.s  = smc.error_dot + smc.c * smc.error;
        smc.ds = -smc.k * sat(smc.s, smc.phi) - smc.eps * smc.s;
        smc.u  = (smc.target_angle_ddot ...
                  - smc.c * smc.error_dot ...
                  + smc.ds) * smc.J * smc.k_c;    % 与 alg_smc(2).c 一致
    end

    % ---- 输出限幅 ----
    smc.u = max(min(smc.u, smc.u_max), -smc.u_max);
    u = smc.u;
end


function C = smc_consts()
% smc_consts  枚举常量定义(与 alg_smc.h 对应)
% 单独做成函数是因为: 脚本的局部函数不能访问脚本工作区中的变量
    C.SMC_SWITCH_SIGN    = 0;   % 符号函数
    C.SMC_SWITCH_SAT     = 1;   % 饱和函数
    C.SMC_SWITCH_TANH    = 2;   % 双曲正切
    C.SMC_SWITCH_SIGMOID = 3;   % Sigmoid
    C.SMC_REACH_CONST = 0;      % 等速趋近律
    C.SMC_REACH_EXP   = 1;      % 指数趋近律
    C.SMC_REACH_POWER = 2;      % 幂次趋近律
end
