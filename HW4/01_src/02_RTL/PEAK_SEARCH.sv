/******************************************************************************
* Copyright (C) 2026 Marco
*
* File Name:    PEAK_SEARCH.sv
* Project:      NTU-DCIC HW4 Symbol Boundary Detector
* Module:       PEAK_SEARCH
* Author:       Marco <harry2963753@gmail.com>
*
******************************************************************************/
`include "define.vh"

// Threshold-gated local peak searcher (same rule as 01_src/00_MATLAB/sbd_peak_search.m)
//   - a search starts at the first sample with above = 1 (once the previous region fell below eta)
//   - the candidate is replaced by a strictly larger above-threshold value, which restarts cnt
//   - the candidate is declared after W_WIN further samples without a larger value,
//     or at the last sample of the record (flush)
//   - after a declaration a new search needs one sample with above = 0
module PEAK_SEARCH (
    input clk,
    input rst_n,
    input in_valid,
    input in_last,
    input [`W_IDX-1:0] in_idx,
    input [`W_PHI-1:0] in_phi,
    input in_above,
    output logic det_valid,
    output logic [`W_IDX-1:0] det_idx,
    output logic [`W_PHI-1:0] det_phi
);

    //=============================================================
    //                         Registers
    //=============================================================
    logic run, hold;
    logic [`W_PHI-1:0] pk_phi;
    logic [`W_IDX-1:0] pk_idx;
    logic [`W_CNT-1:0] cnt;

    //=============================================================
    //                      Next-State Logic
    //=============================================================
    logic run_n, hold_n, emit;
    logic [`W_PHI-1:0] pk_phi_n;
    logic [`W_IDX-1:0] pk_idx_n;
    logic [`W_CNT-1:0] cnt_n;

    always_comb begin : NEXT
        run_n = run;
        hold_n = hold;
        pk_phi_n = pk_phi;
        pk_idx_n = pk_idx;
        cnt_n = cnt;
        emit = 1'b0;
        if(in_valid) begin
            if(!run) begin
                if(in_above && !hold) begin
                    run_n = 1'b1;
                    pk_phi_n = in_phi;
                    pk_idx_n = in_idx;
                    cnt_n = '0;
                end
                if(!in_above) hold_n = 1'b0;
            end
            else if(in_above && in_phi > pk_phi) begin
                pk_phi_n = in_phi;
                pk_idx_n = in_idx;
                cnt_n = '0;
            end
            else cnt_n = cnt + 1'b1;
            if(run_n && (cnt_n == `W_CNT'(`W_WIN) || in_last)) begin
                emit = 1'b1;
                run_n = 1'b0;
                hold_n = in_above;
            end
        end
    end

    always_ff @(posedge clk or negedge rst_n) begin : STATE_REG
        if(!rst_n) begin
            run <= 1'b0;
            hold <= 1'b0;
            pk_phi <= '0;
            pk_idx <= '0;
            cnt <= '0;
        end
        else begin
            run <= run_n;
            hold <= hold_n;
            pk_phi <= pk_phi_n;
            pk_idx <= pk_idx_n;
            cnt <= cnt_n;
        end
    end

    //=============================================================
    //                        Detection Output
    //=============================================================
    always_ff @(posedge clk or negedge rst_n) begin : DET_REG
        if(!rst_n) begin
            det_valid <= 1'b0;
            det_idx <= '0;
            det_phi <= '0;
        end
        else begin
            det_valid <= emit;
            if(emit) begin
                det_idx <= pk_idx_n;
                det_phi <= pk_phi_n;
            end
        end
    end

endmodule
