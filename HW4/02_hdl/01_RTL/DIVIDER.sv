/******************************************************************************
* Copyright (C) 2026 Marco
*
* File Name:    DIVIDER.sv
* Project:      NTU-DCIC HW4 Symbol Boundary Detector
* Module:       DIVIDER
* Author:       Marco <harry2963753@gmail.com>
*
******************************************************************************/
`include "define.vh"

// Pipelined restoring divider: phi = floor(num * 2^F_PHI / den) in U2.14
//   - one quotient bit per stage, MSB first, W_PHI stages (throughput 1 sample/cycle)
//   - saturates to 2^W_PHI - 1 when num >= 4*den, outputs 0 when den = 0
module DIVIDER (
    input clk,
    input rst_n,
    input in_valid,
    input in_last,
    input [`W_IDX-1:0] in_idx,
    input [`W_NUM-1:0] in_num,
    input [`W_DEN-1:0] in_den,
    output logic out_valid,
    output logic out_last,
    output logic [`W_IDX-1:0] out_idx,
    output logic [`W_PHI-1:0] out_phi
);

    localparam int NS = `W_PHI;

    //=============================================================
    //                     Pipeline Registers
    //=============================================================
    // r_*[k] is the register after stage k (k = 1..NS)
    logic [`W_REM-1:0] r_rem [1:NS];
    logic [`W_DEN-1:0] r_den [1:NS];
    logic [`W_PHI-1:0] r_quo [1:NS];
    logic r_sat [1:NS];
    logic r_zero [1:NS];
    logic r_vld [1:NS];
    logic r_lst [1:NS];
    logic [`W_IDX-1:0] r_idx [1:NS];

    //=============================================================
    //                        Stage Inputs
    //=============================================================
    // c_*[k] is the value entering stage k (k = 0 is the module input)
    logic [`W_REM-1:0] c_rem [0:NS-1];
    logic [`W_DEN-1:0] c_den [0:NS-1];
    logic [`W_PHI-1:0] c_quo [0:NS-1];
    logic c_sat [0:NS-1];
    logic c_zero [0:NS-1];
    logic c_vld [0:NS-1];
    logic c_lst [0:NS-1];
    logic [`W_IDX-1:0] c_idx [0:NS-1];
    logic [`W_REM-1:0] trial [0:NS-1];
    logic ge [0:NS-1];

    always_comb begin : STAGE_IN
        c_rem[0] = `W_REM'(in_num) << `F_PHI;
        c_den[0] = in_den;
        c_quo[0] = '0;
        c_sat[0] = ({1'b0, in_num} >= ({2'b0, in_den} << 2));
        c_zero[0] = (in_den == '0);
        c_vld[0] = in_valid;
        c_lst[0] = in_last;
        c_idx[0] = in_idx;
        for(int k = 1; k < NS; k++) begin
            c_rem[k] = r_rem[k];
            c_den[k] = r_den[k];
            c_quo[k] = r_quo[k];
            c_sat[k] = r_sat[k];
            c_zero[k] = r_zero[k];
            c_vld[k] = r_vld[k];
            c_lst[k] = r_lst[k];
            c_idx[k] = r_idx[k];
        end
    end

    //=============================================================
    //            Restoring Stages (stage k -> bit NS-1-k)
    //=============================================================
    always_comb begin : TRIAL
        for(int k = 0; k < NS; k++) begin
            trial[k] = `W_REM'(c_den[k]) << (NS - 1 - k);
            ge[k] = (c_rem[k] >= trial[k]);
        end
    end

    always_ff @(posedge clk or negedge rst_n) begin : STAGE_REG
        if(!rst_n) begin
            for(int k = 1; k <= NS; k++) begin
                r_rem[k] <= '0;
                r_den[k] <= '0;
                r_quo[k] <= '0;
                r_sat[k] <= 1'b0;
                r_zero[k] <= 1'b0;
                r_vld[k] <= 1'b0;
                r_lst[k] <= 1'b0;
                r_idx[k] <= '0;
            end
        end
        else begin
            for(int k = 0; k < NS; k++) begin
                r_vld[k+1] <= c_vld[k];
                r_lst[k+1] <= c_lst[k];
                if(c_vld[k]) begin
                    r_rem[k+1] <= (ge[k])? c_rem[k] - trial[k] : c_rem[k];
                    r_quo[k+1] <= c_quo[k] | (`W_PHI'(ge[k]) << (NS - 1 - k));
                    r_den[k+1] <= c_den[k];
                    r_sat[k+1] <= c_sat[k];
                    r_zero[k+1] <= c_zero[k];
                    r_idx[k+1] <= c_idx[k];
                end
            end
        end
    end

    //=============================================================
    //                    Saturation & Output
    //=============================================================
    always_comb begin : OUTPUT
        out_valid = r_vld[NS];
        out_last = r_lst[NS];
        out_idx = r_idx[NS];
        if(r_zero[NS]) out_phi = '0;
        else if(r_sat[NS]) out_phi = '1;
        else out_phi = r_quo[NS];
    end

endmodule
