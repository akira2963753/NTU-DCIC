/******************************************************************************
* Copyright (C) 2026 Marco
*
* File Name:    SBD.sv
* Project:      NTU-DCIC HW4 Symbol Boundary Detector
* Module:       SBD
* Author:       Marco <harry2963753@gmail.com>
*
******************************************************************************/
`include "define.vh"

// Symbol boundary detector with the normalized metric
//   Phi_NM(m) = |sum_{r=0}^{R-1} z(m-r) z*(m-r-L)|^2 / (sum_{r=0}^{R-1} |z(m-r)|^2)^2 > eta
// DCORR -> NORM_SQ -> DIVIDER -> THRESH -> PEAK_SEARCH, one input sample per cycle
module SBD (
    input clk,
    input rst_n,
    input in_valid,
    input in_last,
    input signed [`W_IN-1:0] in_re,
    input signed [`W_IN-1:0] in_im,
    input [`W_PHI-1:0] eta,
    output logic phi_valid,
    output logic [`W_IDX-1:0] phi_idx,
    output logic [`W_PHI-1:0] phi,
    output logic phi_above,
    output logic det_valid,
    output logic [`W_IDX-1:0] det_idx,
    output logic [`W_PHI-1:0] det_phi
);

    //=============================================================
    //                      Internal Signals
    //=============================================================
    logic dc_valid, dc_last;
    logic [`W_IDX-1:0] dc_idx;
    logic signed [`W_ACC_C-1:0] dc_c_re, dc_c_im;
    logic [`W_ACC_E-1:0] dc_e;
    logic ns_valid, ns_last;
    logic [`W_IDX-1:0] ns_idx;
    logic [`W_NUM-1:0] ns_num;
    logic [`W_DEN-1:0] ns_den;
    logic dv_valid, dv_last;
    logic [`W_IDX-1:0] dv_idx;
    logic [`W_PHI-1:0] dv_phi;
    logic th_last;

    //=============================================================
    //                         Submodules
    //=============================================================
    DCORR u_DCORR (
        .clk(clk),
        .rst_n(rst_n),
        .in_valid(in_valid),
        .in_last(in_last),
        .in_re(in_re),
        .in_im(in_im),
        .out_valid(dc_valid),
        .out_last(dc_last),
        .out_idx(dc_idx),
        .out_c_re(dc_c_re),
        .out_c_im(dc_c_im),
        .out_e(dc_e)
    );

    NORM_SQ u_NORM_SQ (
        .clk(clk),
        .rst_n(rst_n),
        .in_valid(dc_valid),
        .in_last(dc_last),
        .in_idx(dc_idx),
        .in_c_re(dc_c_re),
        .in_c_im(dc_c_im),
        .in_e(dc_e),
        .out_valid(ns_valid),
        .out_last(ns_last),
        .out_idx(ns_idx),
        .out_num(ns_num),
        .out_den(ns_den)
    );

    DIVIDER u_DIVIDER (
        .clk(clk),
        .rst_n(rst_n),
        .in_valid(ns_valid),
        .in_last(ns_last),
        .in_idx(ns_idx),
        .in_num(ns_num),
        .in_den(ns_den),
        .out_valid(dv_valid),
        .out_last(dv_last),
        .out_idx(dv_idx),
        .out_phi(dv_phi)
    );

    THRESH u_THRESH (
        .clk(clk),
        .rst_n(rst_n),
        .in_valid(dv_valid),
        .in_last(dv_last),
        .in_idx(dv_idx),
        .in_phi(dv_phi),
        .eta(eta),
        .out_valid(phi_valid),
        .out_last(th_last),
        .out_idx(phi_idx),
        .out_phi(phi),
        .out_above(phi_above)
    );

    PEAK_SEARCH u_PEAK_SEARCH (
        .clk(clk),
        .rst_n(rst_n),
        .in_valid(phi_valid),
        .in_last(th_last),
        .in_idx(phi_idx),
        .in_phi(phi),
        .in_above(phi_above),
        .det_valid(det_valid),
        .det_idx(det_idx),
        .det_phi(det_phi)
    );

endmodule
