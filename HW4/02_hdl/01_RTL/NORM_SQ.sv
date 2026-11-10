/******************************************************************************
* Copyright (C) 2026 Marco
*
* File Name:    NORM_SQ.sv
* Project:      NTU-DCIC HW4 Symbol Boundary Detector
* Module:       NORM_SQ
* Author:       Marco <harry2963753@gmail.com>
*
******************************************************************************/
`include "define.vh"

// Numerator |C|^2 and denominator E^2 of Phi_NM. The ratio is scale invariant, so
// |Re C|, |Im C| and E share one right shift that keeps them inside W_MAG bits
// before squaring (block normalization). Two pipeline stages:
//   S4 leading-one detection and shift, S5 squares
module NORM_SQ (
    input clk,
    input rst_n,
    input in_valid,
    input in_last,
    input [`W_IDX-1:0] in_idx,
    input signed [`W_ACC_C-1:0] in_c_re,
    input signed [`W_ACC_C-1:0] in_c_im,
    input [`W_ACC_E-1:0] in_e,
    output logic out_valid,
    output logic out_last,
    output logic [`W_IDX-1:0] out_idx,
    output logic [`W_NUM-1:0] out_num,
    output logic [`W_DEN-1:0] out_den
);

    // |C| <= sqrt(E(m) E(m-L)) <= R*2^31 = 2^36, so all magnitudes fit in W_ACC_E bits
    localparam int W_ABS = `W_ACC_E;

    //=============================================================
    //              S4: Leading-One Detection & Shift
    //=============================================================
    logic [W_ABS-1:0] a_re, a_im, a_or;
    logic [`W_SH-1:0] sh;
    logic [W_ABS-1:0] s_re, s_im, s_e;
    logic s4_valid, s4_last;
    logic [`W_IDX-1:0] s4_idx;
    logic [`W_MAG-1:0] s4_re, s4_im, s4_e;

    always_comb begin : MAGNITUDE
        a_re = (in_c_re < 0)? W_ABS'(-in_c_re) : W_ABS'(in_c_re);
        a_im = (in_c_im < 0)? W_ABS'(-in_c_im) : W_ABS'(in_c_im);
        a_or = a_re | a_im | in_e;
    end

    // sh = max(0, msb(a_or) - (W_MAG - 1)): the largest magnitude keeps its top W_MAG bits
    always_comb begin : LOD
        sh = '0;
        for(int i = `W_MAG; i < W_ABS; i++) begin
            if(a_or[i]) sh = `W_SH'(i - `W_MAG + 1);
        end
    end

    always_comb begin : SHIFT
        s_re = a_re >> sh;
        s_im = a_im >> sh;
        s_e = in_e >> sh;
    end

    always_ff @(posedge clk or negedge rst_n) begin : S4_REG
        if(!rst_n) begin
            s4_valid <= 1'b0;
            s4_last <= 1'b0;
            s4_idx <= '0;
            s4_re <= '0;
            s4_im <= '0;
            s4_e <= '0;
        end
        else begin
            s4_valid <= in_valid;
            s4_last <= in_last;
            if(in_valid) begin
                s4_idx <= in_idx;
                s4_re <= s_re[`W_MAG-1:0];
                s4_im <= s_im[`W_MAG-1:0];
                s4_e <= s_e[`W_MAG-1:0];
            end
        end
    end

    //=============================================================
    //                        S5: Squares
    //=============================================================
    logic [`W_NUM-1:0] num;
    logic [`W_DEN-1:0] den;

    always_comb begin : SQUARE
        num = `W_NUM'(s4_re) * `W_NUM'(s4_re) + `W_NUM'(s4_im) * `W_NUM'(s4_im);
        den = `W_DEN'(s4_e) * `W_DEN'(s4_e);
    end

    always_ff @(posedge clk or negedge rst_n) begin : S5_REG
        if(!rst_n) begin
            out_valid <= 1'b0;
            out_last <= 1'b0;
            out_idx <= '0;
            out_num <= '0;
            out_den <= '0;
        end
        else begin
            out_valid <= s4_valid;
            out_last <= s4_last;
            if(s4_valid) begin
                out_idx <= s4_idx;
                out_num <= num;
                out_den <= den;
            end
        end
    end

endmodule
