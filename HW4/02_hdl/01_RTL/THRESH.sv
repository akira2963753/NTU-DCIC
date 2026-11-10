/******************************************************************************
* Copyright (C) 2026 Marco
*
* File Name:    THRESH.sv
* Project:      NTU-DCIC HW4 Symbol Boundary Detector
* Module:       THRESH
* Author:       Marco <harry2963753@gmail.com>
*
******************************************************************************/
`include "define.vh"

// Threshold judgement: above = (Phi_NM(m) > eta), one register stage
module THRESH (
    input clk,
    input rst_n,
    input in_valid,
    input in_last,
    input [`W_IDX-1:0] in_idx,
    input [`W_PHI-1:0] in_phi,
    input [`W_PHI-1:0] eta,
    output logic out_valid,
    output logic out_last,
    output logic [`W_IDX-1:0] out_idx,
    output logic [`W_PHI-1:0] out_phi,
    output logic out_above
);

    always_ff @(posedge clk or negedge rst_n) begin : THRESH_REG
        if(!rst_n) begin
            out_valid <= 1'b0;
            out_last <= 1'b0;
            out_idx <= '0;
            out_phi <= '0;
            out_above <= 1'b0;
        end
        else begin
            out_valid <= in_valid;
            out_last <= in_last;
            if(in_valid) begin
                out_idx <= in_idx;
                out_phi <= in_phi;
                out_above <= (in_phi > eta);
            end
        end
    end

endmodule
