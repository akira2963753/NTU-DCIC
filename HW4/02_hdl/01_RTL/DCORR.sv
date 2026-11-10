/******************************************************************************
* Copyright (C) 2026 Marco
*
* File Name:    DCORR.sv
* Project:      NTU-DCIC HW4 Symbol Boundary Detector
* Module:       DCORR
* Author:       Marco <harry2963753@gmail.com>
*
******************************************************************************/
`include "define.vh"

// Delay correlator C(m) = sum_{r=0}^{R-1} z(m-r) z*(m-r-L) and window energy
// E(m) = sum_{r=0}^{R-1} |z(m-r)|^2, both in the recursive form
// acc(m) = acc(m-1) + x(m) - x(m-R). Three pipeline stages:
//   S1 input register, S2 z^-L delay line and products, S3 z^-R delay line and accumulators
module DCORR (
    input clk,
    input rst_n,
    input in_valid,
    input in_last,
    input signed [`W_IN-1:0] in_re,
    input signed [`W_IN-1:0] in_im,
    output logic out_valid,
    output logic out_last,
    output logic [`W_IDX-1:0] out_idx,
    output logic signed [`W_ACC_C-1:0] out_c_re,
    output logic signed [`W_ACC_C-1:0] out_c_im,
    output logic [`W_ACC_E-1:0] out_e
);

    //=============================================================
    //                 S1: Input Register & Index
    //=============================================================
    logic s1_valid, s1_last;
    logic [`W_IDX-1:0] s1_idx, smp_cnt;
    logic signed [`W_IN-1:0] s1_re, s1_im;

    always_ff @(posedge clk or negedge rst_n) begin : S1_REG
        if(!rst_n) begin
            s1_valid <= 1'b0;
            s1_last <= 1'b0;
            s1_idx <= '0;
            s1_re <= '0;
            s1_im <= '0;
            smp_cnt <= '0;
        end
        else begin
            s1_valid <= in_valid;
            s1_last <= in_valid & in_last;
            if(in_valid) begin
                s1_idx <= smp_cnt + 1'b1;
                s1_re <= in_re;
                s1_im <= in_im;
                smp_cnt <= (in_last)? '0 : smp_cnt + 1'b1;
            end
        end
    end

    //=============================================================
    //            S2: z^-L Delay Line & Complex Products
    //=============================================================
    logic signed [`W_IN-1:0] zl_re [0:`L_LAG-1];
    logic signed [`W_IN-1:0] zl_im [0:`L_LAG-1];
    logic signed [`W_PROD-1:0] p_re, p_im;
    logic [`W_PWR-1:0] p_pwr;
    logic s2_valid, s2_last;
    logic [`W_IDX-1:0] s2_idx;
    logic signed [`W_PROD-1:0] s2_p_re, s2_p_im;
    logic [`W_PWR-1:0] s2_pwr;

    // z(m) z*(m-L) = (a + jb)(c - jd) = (ac + bd) + j(bc - ad), |z(m)|^2 = a^2 + b^2
    always_comb begin : PRODUCT
        p_re = `W_PROD'(s1_re) * `W_PROD'(zl_re[`L_LAG-1]) + `W_PROD'(s1_im) * `W_PROD'(zl_im[`L_LAG-1]);
        p_im = `W_PROD'(s1_im) * `W_PROD'(zl_re[`L_LAG-1]) - `W_PROD'(s1_re) * `W_PROD'(zl_im[`L_LAG-1]);
        p_pwr = `W_PWR'(`W_PWR'(s1_re) * `W_PWR'(s1_re) + `W_PWR'(s1_im) * `W_PWR'(s1_im));
    end

    always_ff @(posedge clk or negedge rst_n) begin : S2_REG
        if(!rst_n) begin
            for(int i = 0; i < `L_LAG; i++) begin
                zl_re[i] <= '0;
                zl_im[i] <= '0;
            end
            s2_valid <= 1'b0;
            s2_last <= 1'b0;
            s2_idx <= '0;
            s2_p_re <= '0;
            s2_p_im <= '0;
            s2_pwr <= '0;
        end
        else begin
            s2_valid <= s1_valid;
            s2_last <= s1_last;
            if(s1_valid) begin
                zl_re[0] <= s1_re;
                zl_im[0] <= s1_im;
                for(int i = 1; i < `L_LAG; i++) begin
                    zl_re[i] <= zl_re[i-1];
                    zl_im[i] <= zl_im[i-1];
                end
                s2_idx <= s1_idx;
                s2_p_re <= p_re;
                s2_p_im <= p_im;
                s2_pwr <= p_pwr;
            end
        end
    end

    //=============================================================
    //            S3: z^-R Delay Line & Recursive Sums
    //=============================================================
    logic signed [`W_PROD-1:0] pr_re [0:`R_LEN-1];
    logic signed [`W_PROD-1:0] pr_im [0:`R_LEN-1];
    logic [`W_PWR-1:0] pr_pwr [0:`R_LEN-1];

    always_ff @(posedge clk or negedge rst_n) begin : S3_REG
        if(!rst_n) begin
            for(int i = 0; i < `R_LEN; i++) begin
                pr_re[i] <= '0;
                pr_im[i] <= '0;
                pr_pwr[i] <= '0;
            end
            out_valid <= 1'b0;
            out_last <= 1'b0;
            out_idx <= '0;
            out_c_re <= '0;
            out_c_im <= '0;
            out_e <= '0;
        end
        else begin
            out_valid <= s2_valid;
            out_last <= s2_last;
            if(s2_valid) begin
                pr_re[0] <= s2_p_re;
                pr_im[0] <= s2_p_im;
                pr_pwr[0] <= s2_pwr;
                for(int i = 1; i < `R_LEN; i++) begin
                    pr_re[i] <= pr_re[i-1];
                    pr_im[i] <= pr_im[i-1];
                    pr_pwr[i] <= pr_pwr[i-1];
                end
                out_idx <= s2_idx;
                out_c_re <= out_c_re + s2_p_re - pr_re[`R_LEN-1];
                out_c_im <= out_c_im + s2_p_im - pr_im[`R_LEN-1];
                out_e <= out_e + s2_pwr - pr_pwr[`R_LEN-1];
            end
        end
    end

endmodule
