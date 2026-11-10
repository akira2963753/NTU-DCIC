/******************************************************************************
* Copyright (C) 2026 Marco
*
* File Name:    define.vh
* Project:      NTU-DCIC HW4 Symbol Boundary Detector
* Module:       -
* Author:       Marco <harry2963753@gmail.com>
*
******************************************************************************/
`ifndef DEFINE_VH
    `define DEFINE_VH
    // Algorithm (keep in sync with 01_src/sbd_param.m)
    `define R_LEN 32
    `define L_LAG 128
    `define W_WIN 32
    // Word lengths
    `define W_IN 16
    `define W_PROD 33
    `define W_PWR 32
    `define W_ACC_C 38
    `define W_ACC_E 37
    `define W_MAG 23
    `define W_SH 4
    `define W_NUM 47
    `define W_DEN 46
    `define W_PHI 16
    `define F_PHI 14
    `define W_REM 61
    `define W_IDX 10
    `define W_CNT 6
`endif
